#!/usr/bin/env python3
"""Static contract for the SALMON v2.3.0 and DG integration."""

from pathlib import Path
import re
import subprocess
import stat
import sys


ROOT = Path(__file__).resolve().parents[2]
MARKERS = ("<" * 7, "|" * 7, "=" * 7, ">" * 7)
CONFLICT_MARKER = re.compile(
    rb"^(" + b"|".join(re.escape(marker.encode("ascii")) for marker in MARKERS) + rb")",
    re.MULTILINE,
)


class ContractError(RuntimeError):
    """A required property of the integrated source tree is absent."""


def require(condition, message):
    if not condition:
        raise ContractError(message)


def salmon_project_version(text):
    project = re.search(r"\bproject\s*\(([^)]*)\)", text, re.S | re.I)
    if project is None or re.match(r"\s*SALMON(?:\s|$)", project.group(1), re.I) is None:
        return None
    version = re.search(r"\bVERSION\s+([^\s)]+)", project.group(1), re.I)
    return None if version is None else version.group(1)


def find_conflict_marker(data):
    match = CONFLICT_MARKER.search(data)
    return None if match is None else match.group(1).decode("ascii")


def first_conflict(surfaces):
    for relative, data in surfaces:
        marker = find_conflict_marker(data)
        if marker is not None:
            return relative, marker
    return None


def integration_text_paths():
    """Yield deterministic production/configuration surfaces relevant to integration."""
    for relative in ("CMakeLists.txt", "configure.py"):
        yield ROOT / relative
    for relative in ("cmakefiles", "gnu_makefiles", "src"):
        for path in sorted(ROOT.joinpath(relative).rglob("*")):
            try:
                mode = path.lstat().st_mode
            except OSError as error:
                raise ContractError(f"{path.relative_to(ROOT)}: cannot inspect: {error}") from error
            if stat.S_ISREG(mode):
                yield path


def integration_text_surfaces():
    for path in integration_text_paths():
        try:
            data = path.read_bytes()
        except OSError as error:
            raise ContractError(f"{path.relative_to(ROOT)}: cannot read: {error}") from error
        if b"\0" in data:
            continue
        yield path.relative_to(ROOT).as_posix(), data


def run_self_tests():
    require(
        salmon_project_version("project(SALMON VERSION 2.3.0 LANGUAGES Fortran)")
        == "2.3.0",
        "self-test: exact VERSION 2.3.0 was not accepted",
    )
    require(
        salmon_project_version("project(SALMON VERSION 2.3.01)") != "2.3.0",
        "self-test: VERSION 2.3.01 was incorrectly accepted",
    )
    require(
        salmon_project_version("project(SALMON VERSION 2.3.0.9)") != "2.3.0",
        "self-test: VERSION 2.3.0.9 was incorrectly accepted",
    )
    for marker in MARKERS:
        require(
            find_conflict_marker(f"{marker} integration fixture\n".encode("ascii"))
            == marker,
            f"self-test: conflict marker {marker!r} was not detected",
        )
    fixture = [("cmakefiles/create_test.cmake", ("|" * 7 + " ancestor\n").encode("ascii"))]
    require(
        first_conflict(fixture) == ("cmakefiles/create_test.cmake", "|" * 7),
        "self-test: sole cmakefiles/create_test.cmake ancestor marker was not detected",
    )


def require_tokens(relative, tokens):
    text = (ROOT / relative).read_text(encoding="utf-8", errors="replace")
    missing = [token for token in tokens if token.casefold() not in text.casefold()]
    require(not missing, f"{relative}: missing required token(s): {', '.join(missing)}")


def require_task7_mpi_guards():
    """Require production DG on MPI+ScaLAPACK and old arm on EigenExa."""
    patterns = (
        re.compile(r"^\s*use\s+(?:dg_|rt_dg_|lcfo_wannier_sawf|lcfo_flux\b)", re.I),
        re.compile(r"\b(?:run|initialize|update|propagate|build|write)_rt?_?dg_", re.I),
        re.compile(r"\brun_dg_", re.I),
    )
    for relative in ("src/gs/main_dft.f90", "src/rt/main_tddft.f90"):
        guards = [set()]
        violations = []
        eigenexa_violations = []
        for number, line in enumerate((ROOT / relative).read_text(errors="replace").splitlines(), 1):
            conditional = re.match(r"^\s*#\s*(?:if|ifdef|ifndef)\b(.*)$", line, re.I)
            if conditional:
                positive = {
                    name.upper()
                    for name in re.findall(r"(?:defined\s*\(\s*)?(USE_(?:MPI|SCALAPACK|EIGENEXA))", conditional.group(1), re.I)
                    if not re.search(rf"!\s*defined\s*\(\s*{name}\s*\)", conditional.group(1), re.I)
                }
                guards.append(guards[-1] | positive)
                continue
            if re.match(r"^\s*#\s*(?:else|elif)\b", line, re.I):
                guards[-1] = set() if len(guards) == 1 else set(guards[-2])
                continue
            if re.match(r"^\s*#\s*endif\b", line, re.I):
                if len(guards) > 1:
                    guards.pop()
                continue
            if any(pattern.search(line) for pattern in patterns) and not {"USE_MPI", "USE_SCALAPACK"}.issubset(guards[-1]):
                violations.append(number)
            if re.search(r"(?:use\s+eigenexa_module|\bcall\s+\w*_eigenexa\b)", line, re.I) and "USE_EIGENEXA" not in guards[-1]:
                eigenexa_violations.append(number)
        require(
            not violations,
            f"{relative}: production DG import/dispatch outside USE_MPI+USE_SCALAPACK guard at lines "
            + ", ".join(map(str, violations[:12])),
        )
        require(
            not eigenexa_violations,
            f"{relative}: old full-symmetry call outside USE_EIGENEXA guard at lines "
            + ", ".join(map(str, eigenexa_violations[:12])),
        )


def task7_preprocess(relative, definitions):
    command = ["cpp", "-P", "-traditional-cpp", f"-I{ROOT / 'gnu_makefiles'}"]
    command.extend(f"-D{name}" for name in definitions)
    command.append(str(ROOT / relative))
    result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
    require(result.returncode == 0, f"{relative}: cpp failed: {result.stderr.strip()}")
    return result.stdout


def require_task7_entrypoint_contracts():
    paths = {
        relative: (ROOT / relative).read_text(encoding="utf-8", errors="replace")
        for relative in (
            "src/gs/dc/dcdft.f90",
            "src/gs/main_dft.f90",
            "src/gs/scf_iteration_dft.f90",
            "src/io/inputoutput.f90",
            "src/rt/main_tddft.f90",
        )
    }
    for relative, text in paths.items():
        require(find_conflict_marker(text.encode()) is None, f"{relative}: unresolved conflict marker")

    input_text = paths["src/io/inputoutput.f90"].casefold()
    for token in (
        "lcfo_eigensolver", "lcfo_diag_chefsi_filter_degree",
        "dg_dc_seed_mode", "dg_fragment_wf_checkpoint_mode",
        "dg_fragment_w90_initial_projection", "dg_hybrid_symmetry_energy_window",
        "yn_dg_hybrid_continuation_scf", "yn_dg_hybrid_divided_scf",
    ):
        require(token in input_text, f"input schema lost {token}")
    require("lapack', 'eigenexa', or 'chefsi" in input_text, "LCFO solver validation lost v2.3 solver choices")
    require("requires a build with scalapack support" in input_text, "CheFSI capability diagnostic absent")
    require("production dg requires a build with mpi and scalapack support" in input_text,
            "production DG capability diagnostic absent")

    gs_text = paths["src/gs/main_dft.f90"].casefold()
    rt_text = paths["src/rt/main_tddft.f90"].casefold()
    require(gs_text.count("call run_dg_hybrid_divided_ground_state_for_main") == 1,
            "divided GS dispatch is not exactly once")
    require(gs_text.count("call run_dg_hybrid_continuation_ground_state_for_main") == 1,
            "continuation GS dispatch is not exactly once")
    require(rt_text.count("call run_dg_overlapping_wannier_coefficient_rt") == 1,
            "coefficient RT dispatch is not exactly once")
    require("384" not in gs_text and "384" not in input_text, "hard-coded 384-state policy is forbidden")

    for relative in ("src/gs/main_dft.f90", "src/rt/main_tddft.f90"):
        conventional = task7_preprocess(relative, ())
        require(re.search(r"^\s*use\s+(?:dg_|rt_dg_|lcfo_wannier_sawf)", conventional, re.I | re.M) is None,
                f"{relative}: conventional preprocessing still imports production DG modules")
        require(re.search(r"\bcall\s+run_dg_", conventional, re.I) is None,
                f"{relative}: conventional preprocessing still dispatches production DG")
        production = task7_preprocess(relative, ("USE_MPI", "USE_SCALAPACK"))
        require("run_dg_" in production.casefold(), f"{relative}: production preprocessing lost DG dispatch")


def require_task7_old_arm_eigenexa_contract():
    relative = "src/io/inputoutput.f90"
    message = "bare overlapping-wannier one-shot arm requires a build with eigenexa support"
    predicate = (
        "if(yn_dg_dc_overlapping_wannier=='y'.and."
        "yn_dg_hybrid_divided_scf/='y'.and."
        "yn_dg_hybrid_continuation_scf/='y'.and."
        "yn_dg_hybrid_scf/='y')"
    )
    production_off = task7_preprocess(relative, ("USE_MPI", "USE_SCALAPACK"))
    production_on = task7_preprocess(relative, ("USE_MPI", "USE_SCALAPACK", "USE_EIGENEXA"))
    conventional = task7_preprocess(relative, ())
    off_compact = re.sub(r"\s+|&", "", production_off.casefold())
    require(predicate in off_compact, "EigenExa-OFF input validation lacks the exact bare-arm predicate")
    require(message in production_off.casefold(), "EigenExa-OFF input validation lacks its early diagnostic")
    require(message not in production_on.casefold(), "EigenExa-ON must permit the bare one-shot arm")
    require(message not in conventional.casefold(), "conventional preprocessing must remain unaffected")


run_self_tests()
if "--self-test-only" in sys.argv[1:]:
    print("exact VERSION 2.3.0 accepted: PASS")
    print("near-miss VERSION 2.3.01 rejected: PASS")
    print("extended VERSION 2.3.0.9 rejected: PASS")
    print(f"conflict markers {' '.join(MARKERS)} detected: PASS")
    print("sole cmakefiles/create_test.cmake ancestor marker detected: PASS")
    raise SystemExit(0)

if "--task7-mpi-guard" in sys.argv[1:]:
    require_task7_mpi_guards()
    print("Task 7 DG MPI+ScaLAPACK and optional EigenExa compile guards: PASS")
    raise SystemExit(0)

if "--task7-entrypoints" in sys.argv[1:]:
    require_task7_entrypoint_contracts()
    require_task7_mpi_guards()
    print("Task 7 input, GS, SCF, and RT entrypoint contracts: PASS")
    raise SystemExit(0)

if "--task7-eigenexa-old-arm" in sys.argv[1:]:
    require_task7_old_arm_eigenexa_contract()
    print("Task 7 bare one-shot EigenExa input guard: PASS")
    raise SystemExit(0)

cmake = (ROOT / "CMakeLists.txt").read_text(encoding="utf-8", errors="replace")
actual_version = salmon_project_version(cmake)
require(
    actual_version == "2.3.0",
    f"CMakeLists.txt: SALMON project version is {actual_version!r}, expected '2.3.0'",
)

required = {
    "src/gs/dc/CMakeLists.txt": (
        "dg_fragment_scdm_gauge.f90",
        "dg_fragment_wf_checkpoint.f90",
    ),
    "src/gs/main_dft.f90": (
        "run_dg_hybrid_divided_ground_state_for_main",
        "dg_fragment_wf_checkpoint_mode",
    ),
    "src/rt/main_tddft.f90": ("run_dg_overlapping_wannier_coefficient_rt",),
}
for relative, tokens in required.items():
    require_tokens(relative, tokens)

conflict = first_conflict(integration_text_surfaces())
if conflict is not None:
    relative, marker = conflict
    raise ContractError(
        f"{relative}: unresolved Git conflict marker {marker!r}"
    )

require_task7_mpi_guards()

print("SALMON v2.3.0 + DG integration contract: PASS")
