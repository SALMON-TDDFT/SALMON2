#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", type=Path,
                    default=ROOT / "build-mpi-eigenexa-wannier-lib")
args = parser.parse_args()
BUILD = args.build_dir.resolve()
SAWF = ROOT / "src/gs/dc/lcfo_wannier_sawf.f90"
LCFO = ROOT / "src/gs/dc/lcfo_flux.f90"
W90_ORACLE = (BUILD / "wannier90/src/wannier90-project/pwscf/v6.5/"
              "pw2wannier90.f90")

DRIVER = r"""
program check_sawf_real_harmonics_driver
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  use lcfo_wannier_sawf, only: t_sawf_symop, t_sawf_projection_channel, &
    build_sawf_spd_projection_map, validate_sawf_spd_projection_map, &
    build_sawf_wannier_representation, sawf_real_harmonic_value, &
    sawf_spd_projection_count, write_sawf_projection_block
  implicit none
  type(t_sawf_projection_channel), allocatable :: channels(:), incomplete(:)
  type(t_sawf_symop), allocatable :: ops(:), product_op(:)
  real(8), allocatable :: dmat(:,:,:), dproduct(:,:,:)
  real(8) :: eye(3,3), inversion(3,3), c4(3,3), mirror(3,3), generic(3,3)
  real(8) :: axis(3), point(3), rotated(3), lhs, rhs
  real(8) :: pi, angle, c, s, one_c, maxerr
  integer :: i, j, l, m, op, source_atom, target_atom, source, target, count, unit
  logical :: ok
  character(512) :: message, off_file, on_file
  real(8) :: frac_pos(3,2)

  call get_command_argument(1, off_file)
  call get_command_argument(2, on_file)

  pi = acos(-1.0d0)
  eye = 0.0d0
  eye(1,1)=1.0d0; eye(2,2)=1.0d0; eye(3,3)=1.0d0
  inversion = -eye
  c4 = reshape([0.0d0,1.0d0,0.0d0, -1.0d0,0.0d0,0.0d0, &
                0.0d0,0.0d0,1.0d0], [3,3])
  mirror = eye
  mirror(1,1) = -1.0d0
  axis = [1.0d0, 2.0d0, 3.0d0]
  axis = axis / sqrt(sum(axis**2))
  angle = 0.371d0
  c = cos(angle); s = sin(angle); one_c = 1.0d0-c
  generic = reshape([ &
    c+axis(1)*axis(1)*one_c, axis(2)*axis(1)*one_c+axis(3)*s, axis(3)*axis(1)*one_c-axis(2)*s, &
    axis(1)*axis(2)*one_c-axis(3)*s, c+axis(2)*axis(2)*one_c, axis(3)*axis(2)*one_c+axis(1)*s, &
    axis(1)*axis(3)*one_c+axis(2)*s, axis(2)*axis(3)*one_c-axis(1)*s, c+axis(3)*axis(3)*one_c], [3,3])

  call sawf_spd_projection_count(1, count, ok, message)
  if (.not. ok .or. count /= 9) error stop 'H-containing SAWF count must be complete spd'
  call sawf_spd_projection_count(huge(0)/9+1, count, ok, message)
  if (ok .or. index(message,'overflow') == 0) error stop 'projection count overflow rejection'
  call build_sawf_spd_projection_map(-1, channels, ok, message)
  if (ok .or. allocated(channels)) error stop 'negative atom count rejection'
  call build_sawf_spd_projection_map(2, channels, ok, message)
  if (.not. ok) then; write(*,'(a)') trim(message); error stop 'map construction'; endif
  if (size(channels) /= 18) error stop 'map size'
  do source_atom=1,2
    do i=1,9
      source = (source_atom-1)*9+i
      if (channels(source)%atom /= source_atom) error stop 'atom-major ordering'
      select case(i)
      case(1)
        if (channels(source)%l /= 0 .or. channels(source)%m /= 1) error stop 's ordering'
      case(2:4)
        if (channels(source)%l /= 1 .or. channels(source)%m /= i-1) error stop 'p ordering'
      case(5:9)
        if (channels(source)%l /= 2 .or. channels(source)%m /= i-4) error stop 'd ordering'
      end select
    end do
  end do
  call validate_sawf_spd_projection_map(channels, 2, ok, message)
  if (.not. ok) then; write(*,'(a)') trim(message); error stop 'map validation'; endif
  allocate(incomplete(17)); incomplete=channels(1:17)
  call validate_sawf_spd_projection_map(incomplete, 2, ok, message)
  if (ok .or. index(message,'complete s+p+d') == 0) error stop 'incomplete-shell rejection'

  allocate(ops(5))
  do op=1,size(ops)
    allocate(ops(op)%atom_map(2))
    ops(op)%atom_map = [1,2]
  end do
  ops(1)%R=eye; ops(2)%R=inversion; ops(3)%R=c4
  ops(4)%R=mirror; ops(5)%R=generic
  ops(5)%atom_map=[2,1]
  call build_sawf_wannier_representation(ops, channels, dmat, ok, message)
  if (.not. ok) then; write(*,'(a)') trim(message); error stop 'D_wann build'; endif

  if (maxval(abs(dmat(:,:,1)-identity18())) > 2.0d-13) error stop 'identity'
  do source_atom=1,2
    source=(source_atom-1)*9
    if (abs(dmat(source+1,source+1,2)-1.0d0)>1d-13) error stop 's inversion'
    if (maxval(abs(dmat(source+2:source+4,source+2:source+4,2)+eye))>1d-13) error stop 'p inversion'
    if (maxval(abs(dmat(source+5:source+9,source+5:source+9,2)-identity5()))>1d-13) error stop 'd inversion'
  end do

  point = [0.314d0, -0.271d0, 0.619d0]
  if (abs(sawf_real_harmonic_value(1,1,point)-point(3))>1d-15 .or. &
      abs(sawf_real_harmonic_value(1,2,point)-point(1))>1d-15 .or. &
      abs(sawf_real_harmonic_value(1,3,point)-point(2))>1d-15) &
    error stop 'Wannier90 p mr order is not z,x,y'
  if (abs(sawf_real_harmonic_value(2,1,point)- &
      (2d0*point(3)**2-point(1)**2-point(2)**2)/sqrt(6d0))>1d-15 .or. &
      abs(sawf_real_harmonic_value(2,2,point)-sqrt(2d0)*point(3)*point(1))>1d-15 .or. &
      abs(sawf_real_harmonic_value(2,3,point)-sqrt(2d0)*point(2)*point(3))>1d-15 .or. &
      abs(sawf_real_harmonic_value(2,4,point)-(point(1)**2-point(2)**2)/sqrt(2d0))>1d-15 .or. &
      abs(sawf_real_harmonic_value(2,5,point)-sqrt(2d0)*point(1)*point(2))>1d-15) &
    error stop 'Wannier90 d mr order or normalized scaling'
  maxerr=0.0d0
  do op=1,size(ops)
    rotated=matmul(ops(op)%R,point)
    do source_atom=1,2
      target_atom=ops(op)%atom_map(source_atom)
      do l=0,2
        do m=1,2*l+1
          target=(target_atom-1)*9+shell_offset(l)+m
          lhs=sawf_real_harmonic_value(l,m,rotated)
          rhs=0.0d0
          do i=1,2*l+1
            source=(source_atom-1)*9+shell_offset(l)+i
            rhs=rhs+dmat(target,source,op)*sawf_real_harmonic_value(l,i,point)
          end do
          maxerr=max(maxerr,abs(lhs-rhs))
        end do
      end do
    end do
    if (maxval(abs(matmul(transpose(dmat(:,:,op)),dmat(:,:,op))-identity18()))>5d-13) &
      error stop 'non-unitary D_wann'
  end do
  if (maxerr>5d-13) error stop 'Wannier90 wws(jw,iw) orientation or polynomial transform'

  ops(1)%R(1,1)=ieee_value(0.0d0,ieee_quiet_nan)
  call build_sawf_wannier_representation(ops(1:1),channels,dproduct,ok,message)
  if (ok .or. index(message,'non-finite') == 0) error stop 'NaN rotation rejection'
  ops(1)%R=eye; ops(1)%R(2,2)=ieee_value(0.0d0,ieee_positive_inf)
  call build_sawf_wannier_representation(ops(1:1),channels,dproduct,ok,message)
  if (ok .or. index(message,'non-finite') == 0) error stop 'Inf rotation rejection'
  ops(1)%R=eye

  allocate(product_op(1)); allocate(product_op(1)%atom_map(2))
  product_op(1)%R=matmul(generic,c4)
  product_op(1)%atom_map=[2,1]
  call build_sawf_wannier_representation(product_op,channels,dproduct,ok,message)
  if (.not.ok) then; write(*,'(a)') trim(message); error stop 'product D_wann build'; endif
  if (maxval(abs(dproduct(:,:,1)-matmul(dmat(:,:,5),dmat(:,:,3))))>8d-13) &
    error stop 'representation group multiplication'

  frac_pos(:,1)=[0.0d0,0.25d0,0.5d0]
  frac_pos(:,2)=[0.75d0,0.5d0,0.25d0]
  open(newunit=unit,file=trim(off_file),status='replace',action='write')
  call write_sawf_projection_block(unit,'off',frac_pos,ok,message)
  close(unit)
  if (.not.ok) then; write(*,'(a)') trim(message); error stop 'off formatter'; endif
  open(newunit=unit,file=trim(on_file),status='replace',action='write')
  call write_sawf_projection_block(unit,'file',frac_pos,ok,message)
  close(unit)
  if (.not.ok) then; write(*,'(a)') trim(message); error stop 'on formatter'; endif
  write(*,'(a,es12.4)') 'PASS SAWF real harmonics maxerr=',maxerr

contains
  integer function shell_offset(lval)
    integer,intent(in)::lval
    select case(lval); case(0); shell_offset=0; case(1); shell_offset=1; case(2); shell_offset=4; end select
  end function
  function identity5() result(a)
    real(8)::a(5,5); integer::k
    a=0d0; do k=1,5; a(k,k)=1d0; enddo
  end function
  function identity18() result(a)
    real(8)::a(18,18); integer::k
    a=0d0; do k=1,18; a(k,k)=1d0; enddo
  end function
end program
"""


def run(command, **kwargs):
    return subprocess.run([str(x) for x in command], text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, **kwargs)


cache = (BUILD / "CMakeCache.txt").read_text()
line = next((x for x in cache.splitlines() if x.startswith("CMAKE_Fortran_COMPILER:")), "")
if "=" not in line:
    raise SystemExit("unable to determine Fortran compiler")
fc = line.split("=", 1)[1]
compiler_configs = sorted((BUILD / "CMakeFiles").glob("*/CMakeFortranCompiler.cmake"))
compiler_id = ""
if compiler_configs:
    match = re.search(r'set\(CMAKE_Fortran_COMPILER_ID\s+"([^"]+)"\)',
                      compiler_configs[0].read_text())
    if match:
        compiler_id = match.group(1)
strict_flags = ["-std=f2008", "-fcheck=all"] if compiler_id == "GNU" else []

with tempfile.TemporaryDirectory(prefix="sawf-harmonics-") as td:
    tmp = Path(td)
    (tmp / "config.h").write_text("")
    (tmp / "sym_stub.f90").write_text("""
module sym_sub
contains
  subroutine read_symmetry_file(filename, matrices, ok, message)
    character(*),intent(in)::filename
    real(8),allocatable,intent(out)::matrices(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    allocate(matrices(3,4,0)); ok=.false.; message='stub'
  end subroutine
end module
""")
    (tmp / "driver.f90").write_text(DRIVER)
    result = run([fc, "-cpp", *strict_flags, "-I", tmp, "-J", tmp,
                  tmp / "sym_stub.f90", SAWF, tmp / "driver.f90", "-o", tmp / "check"])
    if result.returncode:
        raise SystemExit(result.stdout)
    off_file = tmp / "off_projection.txt"
    on_file = tmp / "on_projection.txt"
    result = run([tmp / "check", off_file, on_file])
    if result.returncode:
        raise SystemExit(result.stdout)
    print(result.stdout.strip())
    off_lines = off_file.read_text().splitlines()
    on_lines = on_file.read_text().splitlines()
    if off_lines != ["random"]:
        raise SystemExit(f"off projection output changed: {off_lines!r}")
    if len(on_lines) != 2 or not all(line.endswith(":s;p;d") for line in on_lines):
        raise SystemExit(f"SAWF projection output is not per-atom s;p;d: {on_lines!r}")
    coords = []
    for line in on_lines:
        if not line.startswith("f="):
            raise SystemExit(f"SAWF projection line is not fractional: {line!r}")
        coords.append(tuple(float(x) for x in line[2:-6].split(",")))
    if coords != [(0.0, 0.25, 0.5), (0.75, 0.5, 0.25)]:
        raise SystemExit(f"SAWF projection atom order changed: {coords!r}")
    print("PASS actual off/on .win projection block output")

source = LCFO.read_text()
required = [
    "build_sawf_spd_projection_map",
    "sawf_real_harmonic_value",
    "sawf_spd_projection_count",
    "write_sawf_projection_block",
    "SAWF pseudo_channels ordering: complete atom-major shell lmax=",
]
missing = [token for token in required if token not in source]
if missing:
    raise SystemExit("lcfo_flux integration missing: " + ", ".join(missing))
if not re.search(r"if\s*\(\s*complete_shell\s*\)\s*then\s*.*"
                 r"if\s*\(\s*proj_l_raw\(ip_raw\)\s*>\s*projection_lmax\s*\)\s*cycle.*"
                 r"selected_raw\(selected_count\)\s*=\s*ip_raw",
                 source, re.I | re.S):
    raise SystemExit("SAWF .amn must retain the shared atom-major channel order")
seed_start = source.find("subroutine write_wannier_seed_files")
seed_end = source.find("end subroutine write_wannier_seed_files", seed_start)
seed_body = source[seed_start:seed_end]
collective_pattern = re.compile(
    r"call\s+comm_bcast\(projection_failure\s*,\s*dc%icomm_tot\s*,\s*0\s*\).*"
    r"call\s+comm_bcast\(projection_failure_message\s*,\s*dc%icomm_tot\s*,\s*0\s*\).*"
    r"if\s*\(projection_failure\s*/=\s*0\s*\)\s*call\s+lcfo_sawf_fatal",
    re.I | re.S,
)
if not collective_pattern.search(seed_body):
    raise SystemExit("root projection failure must be broadcast before collective fatal")
writer_start = source.find("subroutine write_pseudo_channel_projection_block")
writer_end = source.find("end subroutine write_pseudo_channel_projection_block", writer_start)
writer_body = source[writer_start:writer_end]
if "error stop" in writer_body.lower() or re.search(r"\bstop\b", writer_body, re.I):
    raise SystemExit("root-only projection writer must return failure, not stop")
if "stat=allocation_status" not in writer_body:
    raise SystemExit("root-only projection allocation must use stat=")
validation_start = source.find("subroutine validate_sawf_projection_configuration")
validation_end = source.find("end subroutine validate_sawf_projection_configuration", validation_start)
validation_body = source[validation_start:validation_end]
if "build_sawf_spd_projection_map" in validation_body:
    raise SystemExit("SAWF configuration validation must use checked count without map allocation")
if "sawf_spd_projection_count" not in validation_body:
    raise SystemExit("SAWF configuration validation must use checked complete-shell count")
map_start = source.find("subroutine build_pseudo_channel_ao_candidate_map")
map_end = source.find("end subroutine build_pseudo_channel_ao_candidate_map", map_start)
map_body = source[map_start:map_end]
if not re.search(r"call\s+comm_get_max\(failure_flag\s*,\s*dc%icomm_tot\s*\).*"
                 r"if\s*\(failure_flag\s*/=\s*0\s*\)\s*call\s+lcfo_sawf_fatal",
                 map_body, re.I | re.S):
    raise SystemExit("rank-local SAWF map failure must be reduced before collective fatal")
amn_start = source.find("subroutine write_wannier_amn_file_pseudo_channels")
amn_end = source.find("end subroutine write_wannier_amn_file_pseudo_channels", amn_start)
amn_body = source[amn_start:amn_end]
if re.search(r"\berror\s+stop\b|\bstop\s+['\"]", amn_body, re.I):
    raise SystemExit("pseudo-channel AMN path must not contain rank-local stop/error stop")
if amn_body.count("stat=allocation_status") < 8:
    raise SystemExit("all pseudo-channel AMN allocations must use stat=/errmsg=")
if amn_body.count("call comm_get_max(allocation_failure, dc%icomm_tot)") < 2:
    raise SystemExit("each pseudo-channel AMN allocation phase must reduce rank-local failure")
if amn_body.count("call comm_get_max(zero_norm_failure, dc%icomm_tot)") < 2:
    raise SystemExit("both pseudo-channel zero-norm phases must validate collectively")
selected_collective = re.compile(
    r"call\s+compute_pseudo_channel_projection_chunk\(.*?"
    r"zero_norm_failure\s*=.*?call\s+comm_get_max\(zero_norm_failure\s*,\s*dc%icomm_tot\s*\).*?"
    r"if\s*\(zero_norm_failure\s*/=\s*0\s*\)\s*call\s+lcfo_sawf_fatal.*?"
    r"if\s*\(dc%id_tot\s*==\s*0\s*\)\s*then",
    re.I | re.S,
)
if not selected_collective.search(amn_body):
    raise SystemExit("selected zero norm must be rejected collectively before root-only AMN write")
print("PASS lcfo_flux SAWF ordering integration")

oracle = W90_ORACLE.read_text()
oracle_tokens = [
    "wws(jw,iw,isym)=", "<Rotated Y(jw)|Not rotated Y(iw)>",
    "IF (mr==1) ylm(ir) = p_z", "IF (mr==2) ylm(ir) = px",
    "IF (mr==3) ylm(ir) = py", "IF (mr==1) ylm(ir) = dz2",
    "IF (mr==2) ylm(ir) = dxz", "IF (mr==3) ylm(ir) = dyz",
    "IF (mr==4) ylm(ir) = dx2my2", "IF (mr==5) ylm(ir) = dxy",
]
missing = [token for token in oracle_tokens if token not in oracle]
if missing:
    raise SystemExit("local Wannier90 convention oracle changed: " + ", ".join(missing))
print("PASS local Wannier90 wws orientation and mr ordering oracle")
