#!/usr/bin/env python3
"""Focused executable/static checks for the Task 5 integration surfaces."""

from pathlib import Path
import argparse
import re
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
class ContractError(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise ContractError(message)


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, **kwargs)


def task5_compiler():
    compiler = shutil.which("mpif90") or shutil.which("mpifort")
    require(compiler is not None, "MPI Fortran compiler not found")
    return compiler


def write_config(directory, eigenexa):
    text = (ROOT / "src/config.h.in").read_text()
    substitutions = {
        "USE_MPI": True,
        "USE_SCALAPACK": False,
        "USE_EIGENEXA": eigenexa,
        "USE_LIBXC": False,
        "USE_FFTW": False,
        "USE_WANNIER90": False,
        "HAVE_SPGLIB": False,
        "USE_OPT_ARRAY_PADDING": True,
        "USE_OPT_EXPLICIT_VECTORIZATION": False,
    }
    for name, enabled in substitutions.items():
        replacement = f"#define {name}" if enabled else f"/* #undef {name} */"
        text = re.sub(rf"^#cmakedefine\s+{name}(?:\s+.*)?$", replacement, text, flags=re.M)
    text = re.sub(r"^#cmakedefine\s+(\w+)(?:\s+.*)?$", r"/* #undef \1 */", text, flags=re.M)
    (directory / "config.h").write_text(text)


def check_eigenexa_preprocessing():
    source = ROOT / "src/parallel/communication.f90"
    data = source.read_text()
    require(
        re.search(r'^\s*#include\s+["<]config\.h[">]', data, re.M) is not None,
        "communication.f90 does not consume generated config.h",
    )
    compiler = task5_compiler()
    for enabled, expected, rejected in (
        (True, "irequired = MPI_THREAD_MULTIPLE", "irequired = MPI_THREAD_FUNNELED"),
        (False, "irequired = MPI_THREAD_FUNNELED", "irequired = MPI_THREAD_MULTIPLE"),
    ):
        with tempfile.TemporaryDirectory(prefix="salmon-task5-eigenexa-") as tmp:
            work = Path(tmp)
            write_config(work, enabled)
            preprocessed = run(
                [compiler, "-cpp", "-E", "-I", str(work), str(source)],
                stdout=subprocess.PIPE,
            ).stdout
            require(expected in preprocessed, f"configured preprocessing omitted {expected}")
            require(rejected not in preprocessed, f"configured preprocessing retained {rejected}")
            compile_flags = [
                "-cpp", "-fallow-argument-mismatch", "-I", str(work),
                "-J", str(work), "-c",
            ]
            run([
                compiler, *compile_flags, str(ROOT / "src/misc/unusedvar.f90"),
                "-o", str(work / "unusedvar.o"),
            ])
            run([
                compiler, *compile_flags, str(ROOT / "src/misc/nvtx_wrapper.f90"),
                "-o", str(work / "nvtx_wrapper.o"),
            ])
            run([
                compiler, *compile_flags, str(source),
                "-o", str(work / "communication.o"),
            ])
    print("configured EigenExa ON/OFF preprocessing and compile probes: PASS")


def configured_cmake_sources(use_mpi, use_scalapack):
    with tempfile.TemporaryDirectory(prefix="salmon-task5-cmake-graph-") as tmp:
        build = Path(tmp)
        run(
            [
                "cmake", "-S", str(ROOT), "-B", str(build),
                f"-DUSE_MPI={'ON' if use_mpi else 'OFF'}",
                f"-DUSE_SCALAPACK={'ON' if use_scalapack else 'OFF'}",
                "-DUSE_EIGENEXA=OFF",
                "-DUSE_WANNIER90=OFF", "-DHAVE_SPGLIB=OFF",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
        dependency_file = build / "src/CMakeFiles/salmon.dir/DependInfo.cmake"
        require(dependency_file.is_file(), "CMake did not generate the SALMON dependency graph")
        matches = re.findall(r'"([^"\n]+/src/[^"\n]+\.f90)"', dependency_file.read_text())
    return ["./" + Path(item).relative_to(ROOT / "src").as_posix() for item in matches]


def affected_dg_lcfo(relative):
    return (
        relative.startswith("./common/dg_")
        or relative.startswith("./rt/dg/")
        or relative.startswith("./gs/dc/dg_")
        or relative.startswith("./gs/dc/lcfo_wannier_")
        or relative == "./gs/dc/lcfo_flux.f90"
        or relative in {
            "./gs/dc/dc_scf_convergence.f90",
            "./gs/dc/dc_fragment_occupation.f90",
        }
    )


def expanded_gnu_f90_sources(makefile):
    with tempfile.TemporaryDirectory(prefix="salmon-task5-make-expand-") as tmp:
        probe = Path(tmp) / "probe.mk"
        probe.write_text(
            f"include {ROOT / 'gnu_makefiles' / makefile}\n"
            "$(foreach item,$(F90_SRCS),$(info TASK5_F90=$(item)))\n"
            ".PHONY: task5_print\n"
            "task5_print: ;\n"
        )
        output = run(
            ["make", "-s", "-f", str(probe), "task5_print"],
            cwd=ROOT / "gnu_makefiles",
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    return [line.removeprefix("TASK5_F90=") for line in output.splitlines() if line.startswith("TASK5_F90=")]


def source_modules(relative):
    text = (ROOT / "src" / relative[2:]).read_text(errors="ignore")
    definitions = {
        match.group(1).casefold()
        for match in re.finditer(r"^\s*module\s+(?!procedure\b|subroutine\b|function\b)(\w+)", text, re.M | re.I)
    }
    unguarded = []
    cpp_depth = 0
    for line in text.splitlines():
        if re.match(r"^\s*#\s*(if|ifdef|ifndef)\b", line):
            cpp_depth += 1
            continue
        if re.match(r"^\s*#\s*endif\b", line):
            cpp_depth = max(0, cpp_depth - 1)
            continue
        if cpp_depth == 0:
            unguarded.append(line)
    imports = {
        match.group(1).casefold()
        for match in re.finditer(
            r"^\s*use(?:\s*,[^:]*)?(?:\s*::)?\s*(\w+)",
            "\n".join(unguarded), re.M | re.I,
        )
    }
    return definitions, imports


def authoritative_production_dg_sources():
    snapshot = "49fa4e2103c1cfcab2183ab60322df3c55ba5425"
    declared = []
    for relative, prefix in (
        ("src/common/CMakeLists.txt", "./common/"),
        ("src/gs/dc/CMakeLists.txt", "./gs/dc/"),
        ("src/rt/CMakeLists.txt", "./rt/"),
    ):
        text = run(
            ["git", "show", f"{snapshot}:{relative}"],
            cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ).stdout
        declared.extend(prefix + name for name in re.findall(r"[\w/]+\.f90", text))
    return sorted({item for item in declared if affected_dg_lcfo(item)})


def build_graph_inputs():
    cmake_production = configured_cmake_sources(True, True)
    cmake_mpi_only = configured_cmake_sources(True, False)
    cmake_nonmpi = configured_cmake_sources(False, False)
    authoritative = authoritative_production_dg_sources()
    require(authoritative, "DG snapshot declares no production DG sources")
    listed = expanded_gnu_f90_sources("Makefile.gnu-mpi")
    nonmpi_listed = expanded_gnu_f90_sources("Makefile.gnu-without-mpi")
    return cmake_production, cmake_mpi_only, cmake_nonmpi, authoritative, listed, nonmpi_listed


def check_gnu_policy():
    _, _, _, authoritative, listed, nonmpi_listed = build_graph_inputs()
    make_body = (ROOT / "gnu_makefiles/make.body").read_text()
    require("DG_F90_SRCS" not in make_body, "legacy GNU make still declares a production DG graph")
    leaked_mpi = sorted(set(authoritative).intersection(listed))
    leaked_nonmpi = sorted(set(authoritative).intersection(nonmpi_listed))
    require(not leaked_mpi, "legacy GNU MPI graph registers production DG: " + ", ".join(leaked_mpi))
    require(not leaked_nonmpi, "legacy GNU non-MPI graph registers production DG: " + ", ".join(leaked_nonmpi))
    require("eigen_subdiag_eigenexa.f90" not in make_body, "legacy GNU make registers EigenExa")
    print("legacy GNU graphs remain conventional-only and exclude DG/EigenExa: PASS")


def check_gnu_eigenexa_exclusion():
    _, _, _, _, listed, _ = build_graph_inputs()
    require(
        "./gs/eigen_subdiag_eigenexa.f90" not in listed,
        "legacy GNU make registers eigen_subdiag_eigenexa without EigenExa flags/includes/libs",
    )
    print("legacy GNU graph excludes unsupported EigenExa wrapper: PASS")


def check_gnu_module_order():
    _, _, _, _, listed, _ = build_graph_inputs()
    definitions = {}
    imports = {}
    for relative in listed:
        if "$(" in relative or not (ROOT / "src" / relative[2:]).is_file():
            continue
        source_definitions, source_imports = source_modules(relative)
        imports[relative] = source_imports
        for module in source_definitions:
            definitions[module] = relative
    position = {relative: index for index, relative in enumerate(listed)}
    ordering_errors = []
    listed_set = set(listed)
    for consumer, used_modules in imports.items():
        for module in used_modules:
            provider = definitions.get(module)
            if provider in listed_set and position[provider] >= position[consumer]:
                ordering_errors.append(f"{provider} !< {consumer} (module {module})")
    require(not ordering_errors, "GNU-make DG module ordering errors: " + "; ".join(ordering_errors))
    print("expanded sequential GNU Fortran module order: PASS")


def check_support_matrix():
    cmake_production, cmake_mpi_only, cmake_nonmpi, authoritative, _, _ = build_graph_inputs()
    missing = sorted(set(authoritative) - set(cmake_production))
    require(not missing, "MPI+ScaLAPACK CMake graph missing production DG: " + ", ".join(missing))
    leaked_mpi_only = sorted(set(authoritative).intersection(cmake_mpi_only))
    leaked_nonmpi = sorted(set(authoritative).intersection(cmake_nonmpi))
    require(not leaked_mpi_only, "MPI-only CMake graph registers production DG: " + ", ".join(leaked_mpi_only))
    require(not leaked_nonmpi, "non-MPI CMake graph registers production DG: " + ", ".join(leaked_nonmpi))
    chefsi = "./gs/dc/lcfo_diag_chefsi.f90"
    require(chefsi in cmake_production, "ScaLAPACK CMake graph lost conventional CheFSI")
    require(chefsi not in cmake_mpi_only, "CMake graph registers CheFSI without ScaLAPACK")
    conventional = {
        "./gs/dc/dc_fragment_geometry.f90", "./gs/dc/dcdft.f90",
        "./gs/dc/dcdft_soi.f90", "./gs/dc/lcfo.f90",
        "./gs/dc/lcfo_soi.f90", "./gs/dc/lcfo_soi_init.f90",
    }
    require(
        conventional.issubset(cmake_nonmpi) and conventional.issubset(cmake_mpi_only),
        "non-MPI CMake graph lost conventional DC/LCFO sources: "
        + ", ".join(sorted(conventional - set(cmake_nonmpi))),
    )
    print(
        f"CMake support matrix gates {len(authoritative)} production DG sources "
        "on MPI+ScaLAPACK and preserves conventional DC/LCFO: PASS"
    )


def check_gnu_dg_graph():
    check_gnu_policy()
    check_gnu_module_order()
    check_gnu_eigenexa_exclusion()
    check_support_matrix()
    make = shutil.which("make")
    require(make is not None, "make executable not found")
    run(
        [make, "-B", "-n", "-j8", "-f", "Makefile.gnu-mpi"],
        cwd=ROOT / "gnu_makefiles",
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
    )
    print(
        "legacy GNU policy and CMake production-DG capability matrix: PASS"
    )


def mpi_dependent_sources(sources):
    definitions = {}
    imports = {}
    for relative in sources:
        path = ROOT / "src" / relative[2:]
        if not path.is_file():
            continue
        source_definitions, source_imports = source_modules(relative)
        imports[relative] = source_imports
        for module in source_definitions:
            definitions[module] = relative
    dependent = {relative for relative, used in imports.items() if "mpi" in used}
    changed = True
    while changed:
        changed = False
        for consumer, used in imports.items():
            if consumer in dependent:
                continue
            if any(definitions.get(module) in dependent for module in used):
                dependent.add(consumer)
                changed = True
    return dependent


def check_dummy_alltoallv():
    compiler = shutil.which("gfortran")
    require(compiler is not None, "GNU Fortran compiler not found")
    with tempfile.TemporaryDirectory(prefix="salmon-task5-dummy-") as tmp:
        work = Path(tmp)
        driver = work / "driver.f90"
        driver.write_text(
            "program test_dummy_alltoallv\n"
            "  use communication, only: comm_alltoallv\n"
            "  implicit none\n"
            "  real(8) :: input(4), output(4)\n"
            "  integer :: counts(1), send_displs(1), recv_displs(1)\n"
            "  input = [10d0, 20d0, 30d0, 40d0]\n"
            "  output = -99d0\n"
            "  counts = 2; send_displs = 1; recv_displs = 1\n"
            "  call comm_alltoallv(input, counts, send_displs, output, counts, recv_displs, 0)\n"
            "  if (any(output(2:3) /= [20d0, 30d0])) error stop 1\n"
            "end program\n"
        )
        flags = ["-cpp", "-J", str(work), "-I", str(work)]
        run([compiler, *flags, "-c", str(ROOT / "src/misc/unusedvar.f90"), "-o", str(work / "unusedvar.o")])
        run([compiler, *flags, "-c", str(ROOT / "src/parallel/communication_dummy.f90"), "-o", str(work / "communication_dummy.o")])
        run([compiler, *flags, str(driver), str(work / "communication_dummy.o"), str(work / "unusedvar.o"), "-o", str(work / "driver")])
        run([str(work / "driver")])
    print("dummy one-rank alltoallv displacements/counts: PASS")


parser = argparse.ArgumentParser()
parser.add_argument(
    "check",
    choices=(
        "eigenexa", "gnu", "gnu-policy", "gnu-order", "gnu-eigenexa",
        "matrix", "dummy", "all",
    ),
    nargs="?", default="all",
)
options = parser.parse_args()
if options.check in ("eigenexa", "all"):
    check_eigenexa_preprocessing()
if options.check in ("gnu", "all"):
    check_gnu_dg_graph()
if options.check == "gnu-policy":
    check_gnu_policy()
if options.check == "gnu-order":
    check_gnu_module_order()
if options.check == "gnu-eigenexa":
    check_gnu_eigenexa_exclusion()
if options.check == "matrix":
    check_support_matrix()
if options.check in ("dummy", "all"):
    check_dummy_alltoallv()
