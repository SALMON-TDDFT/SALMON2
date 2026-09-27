#!/usr/bin/env python3
"""Compile-only isolation of frtpx crashes. Never link or execute the probes.

Run from any directory: python3 tools/diagnose_frtpx_lcfo.py --build build
All generated sources, module copies, objects and logs stay in a new temp folder.
The full source is compiled first, followed by controls and procedure-body probes.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

# Deliberately explicit: refuse changed source instead of guessing declaration ends.
STARTS = {
    "initialize_fragment": "nf=size(lcfo_counts);",
    "pack_coefficients": "if(.not.lcfo_rt_active)error stop",
    "lcfo_hse_direct_rotate": "if(.not.lcfo_rt_active.or..not.lcfo_direct_wf)return",
    "lcfo_hse_refresh": "started=wall_seconds();old_count=",
    "refresh_master": "started=wall_seconds()",
    "exchange_continuity": "ng=size(fragment_basis,1);",
    "lcfo_hse_stage": "if(stage==0)then",
    "lcfo_hse_add_action": "if(.not.allocated(hx))error stop",
    "wall_seconds": "call system_clock(tick,rate)",
    "report_timings": "local=0d0;local(1:4)=",
}


def variants(source):
    lines = source.splitlines(keepends=True)
    spans = {}
    current = None
    for i, line in enumerate(lines):
        match = re.match(r"  (?:subroutine|real\(8\) function) (\w+)\(", line)
        if match:
            if current is not None or match[1] not in STARTS:
                raise ValueError("Unexpected procedure structure; update diagnostic markers")
            current = match[1]
            begin = None
        if current and line.strip().startswith(STARTS[current]):
            if begin is not None:
                raise ValueError("Ambiguous body marker: " + current)
            begin = i
        if current and re.match(r"  end (subroutine|function)\s*$", line):
            if begin is None:
                raise ValueError("Missing body marker: " + current)
            spans[current] = (begin, i)
            current = None
    if set(spans) != set(STARTS) or current:
        raise ValueError("Source procedures changed; update diagnostic markers")

    def retain(names):
        result = list(lines)
        for name, (begin, end) in sorted(spans.items(), key=lambda item: -item[1][0]):
            if name not in names:
                result[begin:end] = ["    seconds=0d0\n" if name == "wall_seconds" else "    return\n"]
        return "".join(result)

    return [("all_stubs", retain(set()))] + [("only_" + name, retain({name})) for name in spans]


def refresh_blocks(source):
    """Split this known procedure into balanced top-level executable units.

    This is not a general Fortran parser. Reject unmatched constructs rather
    than emitting a misleading probe. Continued statements stay together.
    """
    lines = source.splitlines(keepends=True)
    start = next(i for i, line in enumerate(lines) if line.startswith("  subroutine refresh_master("))
    begin = next(i for i in range(start, len(lines)) if lines[i].strip() == "started=wall_seconds()")
    end = next(i for i in range(begin, len(lines)) if lines[i].strip() == "end subroutine")
    groups, unit, logical, depth = [], [], "", 0
    for line in lines[begin:end]:
        unit.append(line)
        code = line.strip()
        if not code or code.startswith("!"):
            continue
        logical += " " + code.lstrip("&").rstrip("&").strip()
        if code.endswith("&"):
            continue
        statement = logical.strip().lower()
        logical = ""
        if statement in ("endif", "enddo"):
            depth -= 1
        elif re.match(r"if\s*\(.*\)\s*then$", statement) or re.match(r"do\s+\w+\s*=", statement):
            depth += 1
        if depth < 0:
            raise ValueError("Unbalanced refresh_master block")
        if depth == 0:
            groups.append("".join(unit))
            unit = []
    if depth or logical:
        raise ValueError("Unbalanced/continued refresh_master block")
    if unit:
        groups[-1] += "".join(unit)
    return "".join(lines[:begin]), groups, "".join(lines[end:])


def reduce_indices(items, crashes):
    """Delta-debug balanced units; predicate must distinguish crashes from errors."""
    current, divisions = list(items), 2
    while len(current) >= 2:
        width = (len(current) + divisions - 1) // divisions
        reduced = False
        for first in range(0, len(current), width):
            candidate = current[:first] + current[first + width:]
            if crashes(candidate):
                current = candidate
                divisions = max(2, divisions - 1)
                reduced = True
                break
        if not reduced:
            if divisions >= len(current):
                break
            divisions = min(len(current), divisions * 2)
    return current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, required=True)
    parser.add_argument("--compiler", default="mpifrtpx")
    parser.add_argument("--gnu-check", action="store_true", help="Local probe syntax check with matching GNU build modules")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--reduce-refresh", action="store_true", help="Automatically isolate crashing refresh_master statement blocks")
    args = parser.parse_args()
    build = args.build.resolve()
    repo = Path(__file__).resolve().parents[1]
    source = repo / "src/xc/hse_lcfo_rt.f90"
    text = source.read_text()
    probes = variants(text)
    compiler = shutil.which(args.compiler)
    if not compiler:
        parser.error("Compiler not found: " + args.compiler)
    modules = sorted(p for p in build.iterdir() if p.suffix.lower() in (".mod", ".smod"))
    if not modules:
        parser.error("No module files at build root; use the existing configured build directory")
    output = Path(tempfile.mkdtemp(prefix="salmon-frtpx-diagnostic-"))
    includes = [build, build / "dependencies/libxc/include", build / "dependencies/fftw/include"]
    results = {"source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
               "compiler": compiler, "output": str(output), "cases": []}
    try:
        version = subprocess.run([compiler, "--version" if args.gnu_check else "-V"],
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 universal_newlines=True, timeout=30)
        results["version"] = version.stdout + version.stderr
    except subprocess.TimeoutExpired:
        results["version"] = "version command timed out"
    print("Output: " + str(output), flush=True)
    print(results["version"], flush=True)
    cases = [("full_O0", text, False), ("full_O0_no_openmp", text, True)]
    cases += [(label, body, False) for label, body in probes]
    def compile_case(label, body, no_openmp=False):
        folder = output / label
        folder.mkdir()
        # Each probe has its own module copies: no diagnostic .mod can replace
        # a production .mod, or contaminate a later probe.
        for module in modules:
            shutil.copy2(module, folder / module.name)
        probe = folder / "hse_lcfo_rt.f90"
        probe.write_text(body)
        if args.gnu_check:
            flags = ["-O0", "-cpp", "-ffree-line-length-none", "-I."]
            if not no_openmp:
                flags += ["-fopenmp"]
        else:
            flags = ["-Nfjomplib", "-Cpp", "-SCALAPACK", "-SSL2BLAMP",
                     "-O0", "-Kocl", "-Nlst=t", "-Koptmsg=2", "-Ncheck_std=03s", "-M", "."]
            if not no_openmp:
                flags.insert(0, "-Kopenmp")
        # Module search uses the private copies; only config/dependency headers
        # are read from the build. Do not pass -I<build> (which contains modules).
        shutil.copy2(build / "config.h", folder / "config.h")
        flags += ["-I" + str(path) for path in includes[1:] if path.is_dir()]
        command = [compiler] + flags + ["-c", str(probe), "-o", str(folder / "probe.o")]
        with (folder / "compile.log").open("w") as log:
            try:
                run = subprocess.run(command, cwd=folder, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout)
                status = run.returncode
            except subprocess.TimeoutExpired:
                status = "timeout"
        result = {"case": label, "status": status, "command": command}
        results["cases"].append(result)
        (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
        print("{}: status={}".format(label, status), flush=True)
        if status != 0:
            print("\n".join((folder / "compile.log").read_text(errors="replace").splitlines()[-8:]), flush=True)
        log_text = (folder / "compile.log").read_text(errors="replace")
        return status, "SIGSEGV" in log_text and status != 0

    if args.reduce_refresh:
        isolated = dict(probes)["only_refresh_master"]
        prefix, blocks, suffix = refresh_blocks(isolated)
        def make_source(indices):
            return prefix + "".join(blocks[i] for i in indices) + suffix
        full_status, full_crash = compile_case("refresh_full", isolated)
        empty_status, empty_crash = compile_case("refresh_empty", make_source([]))
        if args.gnu_check:
            for i in range(len(blocks)):
                compile_case("block_{:03d}".format(i), make_source([i]))
        elif not full_crash or empty_status != 0:
            print("Cannot reduce: full case must SIGSEGV and empty case must compile.", flush=True)
        else:
            cache = {}
            def crashes(indices):
                key = tuple(indices)
                if key not in cache:
                    if len(cache) >= 80:
                        raise RuntimeError("80-probe limit reached; inspect saved logs")
                    label = "reduce_{:03d}".format(len(cache))
                    status, crash = compile_case(label, make_source(indices))
                    cache[key] = crash
                return cache[key]
            remaining = reduce_indices(range(len(blocks)), crashes)
            final_source = make_source(remaining)
            # Recheck the result in a fresh folder to detect unstable reproduction.
            final_status, final_crash = compile_case("reduced_confirm", final_source)
            if not final_crash:
                raise RuntimeError("Reduced crash did not reproduce; inspect logs")
            (output / "reduced.f90").write_text(final_source)
            results["retained_blocks"] = remaining
            (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
            print("BEGIN RETAINED EXECUTABLE BLOCKS (declarations/imports unchanged)", flush=True)
            for i in remaining:
                print(blocks[i], flush=True)
            print("END RETAINED EXECUTABLE BLOCKS", flush=True)
            print("Reduced reproducer: " + str(output / "reduced.f90"), flush=True)
    else:
        for label, body, no_openmp in cases:
            compile_case(label, body, no_openmp)
    print("Compile-only diagnostics complete; no production objects changed.", flush=True)
    print("Please return the above statuses and compiler version; full logs: " + str(output), flush=True)
    if args.gnu_check and any(case["status"] != 0 for case in results["cases"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
