#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


root = Path(".").resolve()
cmake = (root / "CMakeLists.txt").read_text()
packages = (root / "cmakefiles/build_required_packages.cmake").read_text()
config = (root / "src/config.h.in").read_text()
salmon_global = (root / "src/io/salmon_global.f90").read_text()
inputoutput = (root / "src/io/inputoutput.f90").read_text()
lcfo_flux = (root / "src/gs/dc/lcfo_flux.f90").read_text()
spglib_builder_path = root / "cmakefiles/Builder/build_spglib.cmake"


def strip_cmake_comments(text):
    active_lines = []
    for line in text.splitlines():
        quote = None
        escaped = False
        for index, char in enumerate(line):
            if escaped:
                escaped = False
            elif char == "\\" and quote:
                escaped = True
            elif char in "\"'":
                quote = None if quote == char else char if quote is None else quote
            elif char == "#" and quote is None:
                line = line[:index]
                break
        active_lines.append(line.rstrip())
    return "\n".join(active_lines)


def normalized_statements(text):
    return [
        re.sub(r"\s+", " ", line.strip())
        for line in strip_cmake_comments(text).splitlines()
        if line.strip()
    ]


def require_line(pattern, statements, message):
    if not any(re.fullmatch(pattern, line, re.I) for line in statements):
        raise SystemExit(message)


def require_sequence(expected, statements, message):
    width = len(expected)
    for index in range(len(statements) - width + 1):
        if statements[index : index + width] == expected:
            return
    raise SystemExit(message)


cmake_statements = normalized_statements(cmake)
package_statements = normalized_statements(packages)

require_line(
    r'option_set\(USE_SPGLIB "[^"]+" OFF\)',
    cmake_statements,
    "USE_SPGLIB must be an active option that defaults to OFF",
)
require_sequence(
    [
        "if (USE_SPGLIB)",
        "include(${CMAKE_SOURCE_DIR}/cmakefiles/Builder/build_spglib.cmake)",
        "endif ()",
    ],
    package_statements,
    "USE_SPGLIB must actively and directly include build_spglib.cmake",
)

active_cmake_lines = strip_cmake_comments(cmake).splitlines()
package_include_index = next(
    (index for index, line in enumerate(active_cmake_lines) if "build_required_packages.cmake" in line),
    None,
)
config_template_index = next(
    (index for index, line in enumerate(active_cmake_lines) if "src/config.h.in" in line),
    None,
)
if package_include_index is None or config_template_index is None or package_include_index >= config_template_index:
    raise SystemExit("spglib package wiring must run before config.h generation")

if not spglib_builder_path.is_file():
    raise SystemExit("missing cmakefiles/Builder/build_spglib.cmake")
spglib_builder = spglib_builder_path.read_text()
builder_statements = normalized_statements(spglib_builder)
require_line(
    r"find_package\(Spglib 2\.1 CONFIG REQUIRED\)",
    builder_statements,
    "spglib wiring must require find_package(Spglib 2.1 CONFIG REQUIRED)",
)
require_sequence(
    [
        "if (NOT TARGET Spglib::symspg)",
        'message(FATAL_ERROR "Spglib package does not export required target Spglib::symspg")',
        "endif ()",
    ],
    builder_statements,
    "spglib wiring must require the official Spglib::symspg imported target",
)
require_line(
    r"set\(EXTERNAL_LIBS Spglib::symspg \$\{EXTERNAL_LIBS\}\)",
    builder_statements,
    "Spglib::symspg must be appended to EXTERNAL_LIBS",
)
require_line(
    r"set\(HAVE_SPGLIB TRUE\)",
    builder_statements,
    "HAVE_SPGLIB must be actively set to TRUE",
)
for forbidden in ["ExternalProject", "FetchContent", "GIT_REPOSITORY", "URL"]:
    if any(forbidden in line for line in builder_statements):
        raise SystemExit(f"spglib builder must not download dependencies: {forbidden}")

if "#cmakedefine HAVE_SPGLIB" not in config:
    raise SystemExit("config.h.in must expose HAVE_SPGLIB")

for declaration in [
    r"character\s*\([^)]*\)\s*::\s*wannier_site_symmetry",
    r"character\s*\(\s*256\s*\)\s*::\s*wannier_symmetry_file",
    r"real\s*\(\s*8\s*\)\s*::\s*wannier_symmetry_tolerance",
]:
    if not re.search(declaration, salmon_global, re.I):
        raise SystemExit(f"missing SAWF input declaration: {declaration}")

dc_namelist = inputoutput.split("namelist/dc/", 1)[1].split("namelist/dg_fragment/", 1)[0]
for token in ["wannier_site_symmetry", "wannier_symmetry_file", "wannier_symmetry_tolerance"]:
    if token not in dc_namelist:
        raise SystemExit(f"{token} must be in &dc")

for pattern, message in [
    (r"wannier_site_symmetry\s*=\s*['\"]off['\"]", "wannier_site_symmetry must default to off"),
    (r"wannier_symmetry_file\s*=\s*['\"]sym\.dat['\"]", "wannier_symmetry_file must default to sym.dat"),
    (r"wannier_symmetry_tolerance\s*=\s*1(?:\.0)?d-6", "wannier_symmetry_tolerance must default to 1d-6"),
    (r"call\s+comm_bcast\s*\(\s*wannier_site_symmetry", "wannier_site_symmetry must be broadcast"),
    (r"call\s+comm_bcast\s*\(\s*wannier_symmetry_file", "wannier_symmetry_file must be broadcast"),
    (r"call\s+comm_bcast\s*\(\s*wannier_symmetry_tolerance", "wannier_symmetry_tolerance must be broadcast"),
    (r"write\s*\([^\n]*\)[^\n]*['\"]wannier_site_symmetry['\"]", "variables.log must include wannier_site_symmetry"),
    (r"write\s*\([^\n]*\)[^\n]*['\"]wannier_symmetry_file['\"]", "variables.log must include wannier_symmetry_file"),
    (r"write\s*\([^\n]*\)[^\n]*['\"]wannier_symmetry_tolerance['\"]", "variables.log must include wannier_symmetry_tolerance"),
]:
    if not re.search(pattern, inputoutput, re.I):
        raise SystemExit(message)

validation = inputoutput.split("select case(trim(wannier_site_symmetry))", 1)[1].split("end select", 1)[0]
if not re.search(r"case\s*\(\s*['\"]off['\"]\s*,\s*['\"]auto['\"]\s*,\s*['\"]file['\"]\s*\)", validation, re.I):
    raise SystemExit("wannier_site_symmetry must accept off, auto, and file")
for message in [
    "wannier_site_symmetry must be off, auto, or file",
    "wannier_site_symmetry='file' requires nonblank wannier_symmetry_file",
    "wannier_symmetry_tolerance must be positive",
    "wannier_site_symmetry='auto' requires SALMON built with USE_SPGLIB=ON",
]:
    if not re.search(r"call\s+sawf_input_fatal\s*\([^\n]*" + re.escape(message), inputoutput, re.I):
        raise SystemExit(f"SAWF validation must use sawf_input_fatal: {message}")

sawf_loader = lcfo_flux.split("subroutine generate_sawf_dmn", 1)[1].split(
    "end subroutine generate_sawf_dmn", 1
)[0]
if "wannier_symmetry_file" not in sawf_loader:
    raise SystemExit("file-mode SAWF loader must use wannier_symmetry_file")
if re.search(r"import_run_root_dir\s*\(\s*\)\s*//\s*['\"]sym\.dat['\"]", sawf_loader, re.I):
    raise SystemExit("file-mode SAWF loader must not hardcode run-root sym.dat")
if not re.search(
    r"wannier_symmetry_file\s*\(\s*1\s*:\s*1\s*\)\s*==\s*['\"]/['\"]",
    sawf_loader,
    re.I,
):
    raise SystemExit("file-mode SAWF loader must preserve absolute symmetry paths")
if not re.search(
    r"trim\s*\(\s*import_run_root_dir\s*\(\s*\)\s*\)\s*//\s*"
    r"trim\s*\(\s*wannier_symmetry_file\s*\)",
    sawf_loader,
    re.I,
):
    raise SystemExit("file-mode SAWF loader must resolve relative paths against the import run root")

fatal_helper = inputoutput.split("subroutine sawf_input_fatal", 1)[1].split(
    "end subroutine sawf_input_fatal", 1
)[0]
for pattern, message in [
    (r"comm_is_root\s*\(\s*nproc_id_global\s*\)", "SAWF fatal messages must be root-only"),
    (r"call\s+end_parallel", "SAWF fatal handling must finalize parallel execution"),
    (r"stop\s+1", "SAWF fatal handling must exit nonzero"),
]:
    if not re.search(pattern, fatal_helper, re.I):
        raise SystemExit(message)


def check_mock_spglib_configure():
    with tempfile.TemporaryDirectory(prefix="salmon-spglib-mock-") as temp:
        temp_path = Path(temp)
        package_dir = temp_path / "package"
        build_dir = temp_path / "build"
        package_dir.mkdir()
        mock_library = package_dir / "libmock_spglib.a"
        mock_library.touch()
        (package_dir / "SpglibConfig.cmake").write_text(
            "add_library(Spglib::symspg UNKNOWN IMPORTED)\n"
            f'set_target_properties(Spglib::symspg PROPERTIES IMPORTED_LOCATION "{mock_library}")\n'
        )
        (package_dir / "SpglibConfigVersion.cmake").write_text(
            'set(PACKAGE_VERSION "2.1.0")\n'
            'set(PACKAGE_VERSION_COMPATIBLE TRUE)\n'
        )
        result = subprocess.run(
            [
                "cmake",
                "-S",
                str(root),
                "-B",
                str(build_dir),
                "-DUSE_SPGLIB=ON",
                f"-DSpglib_DIR={package_dir}",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode != 0:
            raise SystemExit(f"mock spglib configure failed:\n{result.stdout}")
        generated_config = (build_dir / "config.h").read_text()
        if "#define HAVE_SPGLIB" not in generated_config:
            raise SystemExit("mock spglib configure did not define HAVE_SPGLIB")
        link_command = (build_dir / "src/CMakeFiles/salmon.dir/link.txt").read_text()
        if str(mock_library) not in link_command:
            raise SystemExit("mock Spglib::symspg target is absent from the SALMON link command")

        incompatible_dir = temp_path / "incompatible"
        incompatible_build = temp_path / "incompatible-build"
        incompatible_dir.mkdir()
        (incompatible_dir / "SpglibConfig.cmake").write_text(
            "add_library(Spglib::symspg UNKNOWN IMPORTED)\n"
        )
        (incompatible_dir / "SpglibConfigVersion.cmake").write_text(
            'set(PACKAGE_VERSION "2.0.0")\n'
            'set(PACKAGE_VERSION_COMPATIBLE FALSE)\n'
        )
        incompatible = subprocess.run(
            [
                "cmake", "-S", str(root), "-B", str(incompatible_build),
                "-DUSE_SPGLIB=ON", f"-DSpglib_DIR={incompatible_dir}",
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if incompatible.returncode == 0 or "2.1" not in incompatible.stdout:
            raise SystemExit(
                "Spglib 2.0 package must be rejected with the 2.1 requirement:\n"
                + incompatible.stdout
            )


def check_invalid_input_exit(executable):
    cases = [
        ("wannier_site_symmetry='invalid'", "wannier_site_symmetry must be off, auto, or file"),
        (
            "wannier_site_symmetry='file', wannier_symmetry_file=' '",
            "wannier_site_symmetry='file' requires nonblank wannier_symmetry_file",
        ),
        ("wannier_symmetry_tolerance=0d0", "wannier_symmetry_tolerance must be positive"),
        ("wannier_site_symmetry='auto'", "requires SALMON built with USE_SPGLIB=ON"),
    ]
    for setting, expected_message in cases:
        result = subprocess.run(
            [str(executable)],
            input=f"&dc\n {setting}\n/\n",
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if result.returncode == 0:
            raise SystemExit(f"invalid SAWF input exited zero: {setting}\n{result.stdout}")
        if expected_message not in result.stdout:
            raise SystemExit(f"invalid SAWF input missed diagnostic: {setting}\n{result.stdout}")


parser = argparse.ArgumentParser()
parser.add_argument("--salmon-executable", type=Path)
args = parser.parse_args()

check_mock_spglib_configure()
if args.salmon_executable:
    check_invalid_input_exit(args.salmon_executable.resolve())

print("SAWF Task 1 build, input, configure, and fatal-exit wiring is present")
