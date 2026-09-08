#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument(
    "--build-dir",
    type=Path,
    default=ROOT / "build-mpi-eigenexa-wannier-lib",
    help="configured SALMON build directory providing compiler metadata and objects",
)
args = parser.parse_args()
BUILD = args.build_dir.resolve()
FIXTURES = ROOT / "tests/dg/fixtures/sawf_symmetry"
SAWF_SOURCE = ROOT / "src/gs/dc/lcfo_wannier_sawf.f90"
SYMMETRY_SOURCE = ROOT / "src/symmetry/symmetry.f90"
COMMUNICATION_OBJECT = (
    BUILD / "src/CMakeFiles/salmon.dir/parallel/communication.f90.o"
)
NVTX_OBJECT = BUILD / "src/CMakeFiles/salmon.dir/misc/nvtx.f90.o"

DRIVER = r"""
program check_sawf_symmetry_file_driver
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use sym_sub, only: read_symmetry_file, use_symmetry, DISPLAY, SymMatA, SymMatB
  use lcfo_wannier_sawf, only: t_sawf_symop, load_sawf_symmetry_file
  implicit none
  character(1024) :: filename, case_name, message
  real(8), allocatable :: raw_ops(:,:,:), frac_pos(:,:)
  integer, allocatable :: species(:)
  type(t_sawf_symop), allocatable :: ops(:)
  real(8) :: lattice(3,3), tolerance, expected_R(3,3)
  logical :: ok

  call get_command_argument(1, filename)
  call get_command_argument(2, case_name)
  lattice = 0.0d0
  lattice(1,1) = 1.0d0
  lattice(2,2) = 1.0d0
  lattice(3,3) = 1.0d0
  expected_R = lattice
  tolerance = 1.0d-8

  if (trim(case_name) == 'reader') then
    use_symmetry = .true.
    DISPLAY = .true.
    allocate(SymMatA(3,4,1), SymMatB(3,4,1))
    SymMatA = 7.0d0
    SymMatB = 9.0d0
    call read_symmetry_file(trim(filename), raw_ops, ok, message)
    if (.not. ok) then
      write(*,'(a)') trim(message)
      stop 2
    end if
    if (size(raw_ops,3) /= 1) error stop 'reader operation count'
    if (.not. use_symmetry .or. .not. DISPLAY) error stop 'reader changed flags'
    if (any(SymMatA /= 7.0d0) .or. any(SymMatB /= 9.0d0)) &
      error stop 'reader changed legacy matrices'
    write(*,'(a)') 'PASS reader'
    stop
  end if

  select case (trim(case_name))
  case ('periodic', 'large_translation')
    allocate(frac_pos(3,2), species(2))
    frac_pos = 0.0d0
    frac_pos(1,2) = 0.5d0
    species = 1
  case ('large_frac_pos')
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    frac_pos(1,1) = 1500000000.125d0
    species = 1
    tolerance = 0.3d0
  case ('species')
    allocate(frac_pos(3,2), species(2))
    frac_pos = 0.0d0
    frac_pos(1,2) = 0.5d0
    species = [1, 2]
  case ('matching')
    allocate(frac_pos(3,3), species(3))
    frac_pos = 0.0d0
    frac_pos(1,:) = [0.01d0, 0.08d0, 0.49d0]
    species = 1
    tolerance = 0.11d0
  case ('skew_pbc', 'skew_rotation')
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    if (trim(case_name) == 'skew_pbc') frac_pos(1:2,1) = -0.245d0
    species = 1
    lattice = 0.0d0
    lattice(1,1) = 1.0d0
    lattice(1,2) = 0.5d0
    lattice(2,2) = sqrt(3.0d0)/2.0d0
    lattice(3,3) = 1.0d0
    if (trim(case_name) == 'skew_pbc') tolerance = 0.6d0
  case ('large_unimodular')
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    species = 1
    lattice = 0.0d0
    lattice(1,1) = 1.0d0
    lattice(1,2) = -2000.0d0
    lattice(2,2) = 1.0d0
    lattice(3,3) = 1.0d0
  case ('zero_atoms')
    allocate(frac_pos(3,0), species(0))
  case ('nan_tolerance', 'inf_tolerance', 'nan_lattice', 'inf_lattice', &
        'nan_frac_pos', 'inf_frac_pos')
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    species = 1
    select case (trim(case_name))
    case ('nan_tolerance')
      tolerance = ieee_value(0.0d0, ieee_quiet_nan)
    case ('inf_tolerance')
      tolerance = ieee_value(0.0d0, ieee_positive_inf)
    case ('nan_lattice')
      lattice(1,1) = ieee_value(0.0d0, ieee_quiet_nan)
    case ('inf_lattice')
      lattice(1,1) = ieee_value(0.0d0, ieee_positive_inf)
    case ('nan_frac_pos')
      frac_pos(1,1) = ieee_value(0.0d0, ieee_quiet_nan)
    case ('inf_frac_pos')
      frac_pos(1,1) = ieee_value(0.0d0, ieee_positive_inf)
    end select
  case ('unrepresentable_pbc_bounds')
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    frac_pos(2,1) = -0.25d0
    species = 1
    tolerance = 1.0d0
    lattice = 0.0d0
    lattice(1,1) = 1.0d0
    lattice(1,2) = 1.0d0
    lattice(2,2) = 1.0d-3
    lattice(3,3) = 1.0d0
  case default
    allocate(frac_pos(3,1), species(1))
    frac_pos = 0.0d0
    species = 1
  end select

  call load_sawf_symmetry_file(trim(filename), lattice, frac_pos, species, &
    tolerance, ops, ok, message)
  if (.not. ok) then
    write(*,'(a)') trim(message)
    stop 2
  end if
  if (size(ops) < 1) error stop 'no normalized operations'
  if (.not. allocated(ops(1)%atom_map)) error stop 'atom map not allocated'
  if (any(ops(1)%W /= reshape([1,0,0,0,1,0,0,0,1],[3,3]))) &
    error stop 'identity is not normalized'
  if (maxval(abs(ops(1)%R - expected_R)) > 1.0d-12) &
    error stop 'identity Cartesian rotation'
  if (ops(1)%integer_residual > 1.0d-12 .or. &
      ops(1)%metric_residual > 1.0d-12 .or. &
      ops(1)%orthogonality_residual > 1.0d-12 .or. &
      ops(1)%determinant_residual > 1.0d-12 .or. &
      ops(1)%atom_map_residual > 1.0d-12) error stop 'identity residual fields'
  if (trim(case_name) == 'periodic') then
    if (any(ops(2)%atom_map /= [2,1])) error stop 'periodic translation atom map'
  end if
  if (trim(case_name) == 'large_translation') then
    if (abs(ops(2)%tau(1) - 0.5d0) > 1.0d-12) &
      error stop 'large translation normalization'
    if (any(ops(2)%atom_map /= [2,1])) error stop 'large translation atom map'
  end if
  if (trim(case_name) == 'large_frac_pos') then
    if (abs(ops(2)%atom_map_residual - 0.25d0) > 1.0d-12) &
      error stop 'large fractional closest-image residual'
  end if
  if (trim(case_name) == 'matching') then
    if (any(ops(2)%atom_map /= [2,1,3])) error stop 'augmenting-path atom map'
    if (abs(ops(2)%atom_map_residual - 0.09d0) > 1.0d-12) &
      error stop 'augmenting-path atom residual'
  end if
  if (trim(case_name) == 'skew_pbc') then
    if (abs(ops(2)%atom_map_residual - sqrt(0.2503d0)) > 1.0d-12) then
      write(*,'(a,es24.16)') 'actual skew residual=', ops(2)%atom_map_residual
      error stop 'skew closest-image residual'
    end if
  end if
  if (trim(case_name) == 'skew_rotation') then
    expected_R = 0.0d0
    expected_R(1,1) = 0.5d0
    expected_R(1,2) = sqrt(3.0d0)/2.0d0
    expected_R(2,1) = sqrt(3.0d0)/2.0d0
    expected_R(2,2) = -0.5d0
    expected_R(3,3) = 1.0d0
    if (maxval(abs(ops(2)%R - expected_R)) > 1.0d-12) &
      error stop 'skew Cartesian A W A^-1 convention'
  end if
  if (trim(case_name) == 'large_unimodular') then
    expected_R = 0.0d0
    expected_R(1,2) = 1.0d0
    expected_R(2,1) = 1.0d0
    expected_R(3,3) = 1.0d0
    if (maxval(abs(ops(2)%R - expected_R)) > 1.0d-10) &
      error stop 'large unimodular Cartesian rotation'
  end if
  write(*,'(a,1x,i0)') 'PASS', size(ops)
end program check_sawf_symmetry_file_driver
"""


def run(command, **kwargs):
    return subprocess.run(
        [str(item) for item in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    )


if not SAWF_SOURCE.is_file():
    raise SystemExit(f"missing SAWF operation module: {SAWF_SOURCE}")
if not COMMUNICATION_OBJECT.is_file():
    raise SystemExit(f"missing build artifact: {COMMUNICATION_OBJECT}")

compiler = (BUILD / "CMakeCache.txt").read_text()
compiler_line = next(
    (line for line in compiler.splitlines() if line.startswith("CMAKE_Fortran_COMPILER:")),
    "",
)
if "=" not in compiler_line:
    raise SystemExit("unable to determine CMAKE_Fortran_COMPILER")
fortran = compiler_line.split("=", 1)[1]
compiler_configs = sorted((BUILD / "CMakeFiles").glob("*/CMakeFortranCompiler.cmake"))
compiler_config = compiler_configs[0] if compiler_configs else None
compiler_id = ""
if compiler_config:
    match = re.search(
        r'set\(CMAKE_Fortran_COMPILER_ID\s+"([^"]+)"\)',
        compiler_config.read_text(),
    )
    if match:
        compiler_id = match.group(1)
compile_flags = []
if compiler_id == "GNU":
    compile_flags.extend(["-std=f2008", "-fcheck=bounds", "-ftrapv"])

with tempfile.TemporaryDirectory(prefix="sawf-symmetry-") as temp:
    temp_path = Path(temp)
    driver_path = temp_path / "driver.f90"
    driver_path.write_text(DRIVER)
    compile_result = run(
        [
            fortran,
            "-I",
            BUILD,
            "-J",
            temp_path,
            *compile_flags,
            "-c",
            SYMMETRY_SOURCE,
            SAWF_SOURCE,
            driver_path,
        ],
        cwd=temp_path,
    )
    if compile_result.returncode != 0:
        raise SystemExit(f"focused SAWF driver compile failed:\n{compile_result.stdout}")
    executable = temp_path / "check_sawf_symmetry_file"
    link_result = run(
        [
            fortran,
            temp_path / "symmetry.o",
            temp_path / "lcfo_wannier_sawf.o",
            temp_path / "driver.o",
            COMMUNICATION_OBJECT,
            NVTX_OBJECT,
            "-o",
            executable,
        ]
    )
    if link_result.returncode != 0:
        raise SystemExit(f"focused SAWF driver link failed:\n{link_result.stdout}")

    passing = [
        ("identity.sym", "reader", "PASS reader"),
        ("commented_identity.sym", "reader", "PASS reader"),
        ("commented_identity.sym", "identity", "PASS 1"),
        ("identity.sym", "identity", "PASS 1"),
        ("cubic_inversion.sym", "inversion", "PASS 2"),
        ("rotation_90.sym", "rotation", "PASS 4"),
        ("periodic_translation.sym", "periodic", "PASS 2"),
        ("ambiguous_inversion.sym", "matching", "PASS 2"),
        ("skew_inversion.sym", "skew_pbc", "PASS 2"),
        ("skew_reflection.sym", "skew_rotation", "PASS 2"),
        ("large_unimodular.sym", "large_unimodular", "PASS 2"),
        ("large_translation.sym", "large_translation", "PASS 2"),
        ("cubic_inversion.sym", "large_frac_pos", "PASS 2"),
    ]
    for fixture, case_name, expected in passing:
        result = run([executable, FIXTURES / fixture, case_name])
        if result.returncode != 0 or expected not in result.stdout:
            raise SystemExit(
                f"valid fixture failed: {fixture}\nexit={result.returncode}\n{result.stdout}"
            )

    failing = [
        ("malformed_incomplete.sym", "identity", "incomplete"),
        ("noninteger_rotation.sym", "identity", "near integer"),
        ("missing_identity.sym", "identity", "identity"),
        ("periodic_translation.sym", "species", "same-species"),
        ("nonclosed_group.sym", "rotation", "closure"),
        ("nan_rotation.sym", "identity", "finite"),
        ("inf_rotation.sym", "identity", "finite"),
        ("out_of_range_rotation.sym", "identity", "representable"),
        ("nonunimodular.sym", "identity", "determinant"),
        ("unsupported_int64_determinant.sym", "identity", "determinant intermediate exceeds int64 range"),
        ("empty.sym", "identity", "no operations"),
        ("does_not_exist.sym", "identity", "cannot open"),
        ("extra_field.sym", "identity", "extra field"),
        ("malformed_token.sym", "identity", "malformed"),
        ("identity.sym", "zero_atoms", "at least one atom"),
        ("identity.sym", "nan_tolerance", "tolerance must be finite"),
        ("identity.sym", "inf_tolerance", "tolerance must be finite"),
        ("identity.sym", "nan_lattice", "lattice entries must be finite"),
        ("identity.sym", "inf_lattice", "lattice entries must be finite"),
        ("identity.sym", "nan_frac_pos", "fractional atom positions must be finite"),
        ("identity.sym", "inf_frac_pos", "fractional atom positions must be finite"),
        (
            "cubic_inversion.sym",
            "unrepresentable_pbc_bounds",
            "lattice too ill-conditioned for exact periodic search",
        ),
    ]
    for fixture, case_name, expected in failing:
        result = run([executable, FIXTURES / fixture, case_name])
        if result.returncode == 0 or expected.lower() not in result.stdout.lower():
            raise SystemExit(
                f"invalid fixture diagnostic failed: {fixture}\n"
                f"exit={result.returncode}\n{result.stdout}"
            )

print("SAWF symmetry file reader and normalized operation validation passed")
