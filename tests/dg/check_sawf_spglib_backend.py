#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FIXTURES = ROOT / "tests/dg/fixtures/sawf_spglib"
SAWF_SOURCE = ROOT / "src/gs/dc/lcfo_wannier_sawf.f90"
WRAPPER_SOURCE = ROOT / "src/gs/dc/lcfo_wannier_spglib.c"
SPGLIB_BUILDER = ROOT / "cmakefiles/Builder/build_spglib.cmake"

parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", type=Path, help=argparse.SUPPRESS)
parser.parse_args()


def run(command, **kwargs):
    return subprocess.run(
        [str(item) for item in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    )


def write_mock_package(package_dir, version="2.1.0"):
    package_dir.mkdir()
    (package_dir / "SpglibConfig.cmake").write_text(
        "if (NOT TARGET Spglib::symspg)\n"
        f'  add_library(sawf_mock_spglib STATIC "{FIXTURES / "mock_spglib.c"}")\n'
        f'  target_include_directories(sawf_mock_spglib PUBLIC "{FIXTURES}")\n'
        "  add_library(Spglib::symspg ALIAS sawf_mock_spglib)\n"
        "endif ()\n"
    )
    (package_dir / "SpglibConfigVersion.cmake").write_text(
        f'set(PACKAGE_VERSION "{version}")\n'
        "if (PACKAGE_FIND_VERSION VERSION_GREATER PACKAGE_VERSION)\n"
        "  set(PACKAGE_VERSION_COMPATIBLE FALSE)\n"
        "else ()\n"
        "  set(PACKAGE_VERSION_COMPATIBLE TRUE)\n"
        "endif ()\n"
    )


def configure_build_test(source_dir, build_dir, options):
    configured = run(["cmake", "-S", source_dir, "-B", build_dir, *options])
    if configured.returncode != 0:
        raise SystemExit(f"CMake configure failed:\n{configured.stdout}")
    built = run(["cmake", "--build", build_dir, "--parallel", "2"])
    if built.returncode != 0:
        raise SystemExit(f"CMake build failed:\n{built.stdout}")
    tested = run(["ctest", "--test-dir", build_dir, "--output-on-failure"])
    if tested.returncode != 0:
        raise SystemExit(f"CTest failed:\n{tested.stdout}")


for required in [
    SAWF_SOURCE,
    WRAPPER_SOURCE,
    SPGLIB_BUILDER,
    FIXTURES / "CMakeLists.txt",
    FIXTURES / "driver.F90",
    FIXTURES / "mock_spglib.c",
    FIXTURES / "spglib.h",
]:
    if not required.is_file():
        raise SystemExit(f"missing SAWF Spglib backend input: {required}")

sawf_text = SAWF_SOURCE.read_text()
if not re.search(
    r"max_sawf_symmetry_operations\s*=\s*4096\b", sawf_text, re.I
):
    raise SystemExit("SAWF auto backend must fail fast above the 4096 operation limit")
closure = sawf_text.split("subroutine validate_closure", 1)[1].split(
    "end subroutine validate_closure", 1
)[0]
if "find_indexed_operation" not in closure or re.search(r"do\s+kop\s*=", closure, re.I):
    raise SystemExit("SAWF closure validation must use indexed O(N^2) lookup")
if not re.search(
    r"find_package\(Spglib\s+2\.1\s+CONFIG\s+REQUIRED\)",
    SPGLIB_BUILDER.read_text(),
    re.I,
):
    raise SystemExit("SAWF requires Spglib package version 2.1 or newer")

with tempfile.TemporaryDirectory(prefix="sawf-spglib-") as temp:
    temp_path = Path(temp)
    package_dir = temp_path / "package"
    write_mock_package(package_dir)

    configure_build_test(
        FIXTURES,
        temp_path / "focused-off",
        [f"-DSALMON_ROOT={ROOT}", "-DSAWF_TEST_HAVE_SPGLIB=OFF"],
    )
    configure_build_test(
        FIXTURES,
        temp_path / "focused-on",
        [
            f"-DSALMON_ROOT={ROOT}",
            "-DSAWF_TEST_HAVE_SPGLIB=ON",
            f"-DSpglib_DIR={package_dir}",
        ],
    )

    incompatible_package_dir = temp_path / "package-2.0"
    write_mock_package(incompatible_package_dir, version="2.0.0")
    incompatible_probe = temp_path / "incompatible-probe"
    incompatible_probe.mkdir()
    (incompatible_probe / "CMakeLists.txt").write_text(
        "cmake_minimum_required(VERSION 3.20)\n"
        "project(check_spglib_version NONE)\n"
        f'find_package(Spglib 2.1 CONFIG REQUIRED PATHS "{incompatible_package_dir}" NO_DEFAULT_PATH)\n'
    )
    incompatible = run(
        [
            "cmake", "-S", incompatible_probe, "-B", temp_path / "focused-on-2.0",
        ]
    )
    diagnostic = incompatible.stdout.lower()
    if (
        incompatible.returncode == 0
        or "spglib" not in diagnostic
        or "2.1" not in diagnostic
        or "2.0.0" not in diagnostic
    ):
        raise SystemExit(
            "focused configure must reject Spglib 2.0.0 when 2.1 is required:\n"
            + incompatible.stdout
        )

    full_build = temp_path / "salmon-mock-on"
    configured = run(
        [
            "cmake", "-S", ROOT, "-B", full_build,
            "-DUSE_SPGLIB=ON", "-DUSE_MPI=ON",
            f"-DSpglib_DIR={package_dir}",
        ]
    )
    if configured.returncode != 0:
        raise SystemExit(f"mock-ON SALMON configure failed:\n{configured.stdout}")
    built = run(["cmake", "--build", full_build, "--parallel", "2"])
    if built.returncode != 0:
        raise SystemExit(f"mock-ON SALMON build failed:\n{built.stdout}")

print("SAWF Spglib auto backend, portable focused build, and fallback behavior passed")
print(
    "LIMITATION: HAVE_SPGLIB=ON tests used an official-signature mock; "
    "a released real Spglib was not exercised."
)
