"""Run selected existing component checks; no builds, downloads, GS solves or CI.

Each case runs in a fresh MPI process because native attachment lifetimes differ.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="new result directory")
    parser.add_argument("--suite", nargs="+", choices=("handles", "memory-gs"), default=["handles"])
    parser.add_argument("--library", type=Path)
    parser.add_argument("--build", type=Path)
    parser.add_argument("--fixture", type=Path, action="append", default=[], help="existing GS bundle; repeat for other XC")
    parser.add_argument("--mpiexec", default="mpiexec")
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--dry-run", action="store_true", help="print case commands without creating files or launching MPI")
    args = parser.parse_args()
    if "handles" in args.suite and args.library is None:
        parser.error("handles requires --library")
    if "memory-gs" in args.suite and (args.build is None or not args.fixture):
        parser.error("memory-gs requires --build and at least one --fixture")
    root = Path(__file__).resolve().parent
    output = args.output.resolve()
    cases = []
    if "handles" in args.suite:
        for model in ("ld", "sbe"):
            for ranks in (2, 3):
                name = f"{model}-{ranks}"
                cases.append((name, ranks, [str(root/"test_material_handles_mpi.py"),
                    str(args.library.resolve()), str(output/(name+".json")), "--model", model]))
    if "memory-gs" in args.suite:
        for index, fixture in enumerate(args.fixture):
            for size, split in ((1, "orbital"), (2, "orbital"), (2, "spatial")):
                name = f"memory-gs-{index+1}-{split}-{size}"
                cases.append((name, 2*size, [str(root/"test_memory_gs_mpi.py"),
                    str(args.build.resolve()), str(fixture.resolve()), str(output/name),
                    "--size", str(size), "--split", split]))
    commands = [(name, [args.mpiexec, "-n", str(ranks), args.python, *tail]) for name, ranks, tail in cases]
    if args.dry_run:
        print(json.dumps(dict(commands), indent=2))
        return
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ)
    # Match documented single-thread regression controls without altering host settings.
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    repo = str(root.parent)
    env["PYTHONPATH"] = repo + (os.pathsep+env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    report = {"status": "running", "cases": [], "builds": 0, "SCF_calls": 0, "Actions_dispatched": 0}
    failed = False
    for name, command in commands:
        start = time.monotonic()
        item = {"name": name, "command": command, "status": "failed"}
        try:
            with (output/(name+".log")).open("x") as log:
                result = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT)
            item["returncode"] = result.returncode
            if result.returncode == 0:
                item["status"] = "pass"
        except OSError as exc:
            item["error"] = str(exc)
        item["seconds"] = time.monotonic()-start
        failed |= item["status"] != "pass"
        report["cases"].append(item)
        report["status"] = "failed" if failed else "pass"
        (output/"summary.json").write_text(json.dumps(report, indent=2)+"\n")
        print(name, item["status"], flush=True)
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
