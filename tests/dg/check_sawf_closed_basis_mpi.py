#!/usr/bin/env python3
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build-mpi-eigenexa-wannier-lib"

DRIVER = r"""
program check_sawf_closed_basis_mpi
  use mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use lcfo_wannier_sawf_basis, only: append_sawf_mapped_basis
  implicit none
  real(8) :: source(4,2),candidate(4,4),before(4,4),received(4,2)
  integer :: point_map(4),rank,nproc,peer,ierr,local_failure,global_failure
  logical :: ok
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc /= 2) call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  peer=1-rank
  point_map=[4,3,2,1]
  source=0d0
  source(1,1)=10d0*dble(rank+1)+1d0
  source(2,2)=10d0*dble(rank+1)+2d0
  candidate=-7d0
  call append_sawf_mapped_basis(source,point_map,candidate,2,ok,message)
  call require(ok,'mapped append: '//trim(message),rank)
  call require(candidate(4,2)==source(1,1) .and. candidate(3,3)==source(2,2), &
    'forward target(point_map)=source orientation',rank)
  call MPI_Sendrecv(candidate(:,2:3),8,MPI_DOUBLE_PRECISION,peer,91, &
    received,8,MPI_DOUBLE_PRECISION,peer,91,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
  call require(received(4,1)==10d0*dble(peer+1)+1d0 .and. &
    received(3,2)==10d0*dble(peer+1)+2d0,'rank-pair mapped block exchange',rank)

  before=candidate
  point_map=[4,3,3,1]
  call append_sawf_mapped_basis(source,point_map,candidate,1,ok,message)
  call require(.not.ok .and. index(message,'permutation')>0,'duplicate map rejection',rank)
  call require(all(candidate==before),'invalid map must not modify candidate',rank)

  point_map=[4,3,2,1]
  if(rank==1) source(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call append_sawf_mapped_basis(source,point_map,candidate,1,ok,message)
  local_failure=merge(0,1,ok)
  if(local_failure/=0) write(*,'(a,i0,2a)') '[SAWF-CLOSED-BASIS-LOCAL] rank=',rank, &
    ' details=',trim(message)
  call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  call require(global_failure==1,'collective failure reduction',rank)
  if(rank==0) write(*,'(a)') 'PASS SAWF mapped fragment blocks and rank-local failure'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text,owner)
    logical,intent(in) :: condition
    character(*),intent(in) :: text
    integer,intent(in) :: owner
    if(.not.condition) then
      write(*,'(a,i0,2a)') 'rank=',owner,' ',trim(text)
      call MPI_Abort(MPI_COMM_WORLD,3,ierr)
    end if
  end subroutine require
end program check_sawf_closed_basis_mpi
"""


def run(command, **kwargs):
    return subprocess.run(
        [str(item) for item in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    )


cache = (BUILD / "CMakeCache.txt").read_text()
compiler = next(
    line.split("=", 1)[1]
    for line in cache.splitlines()
    if line.startswith("CMAKE_Fortran_COMPILER:")
)

with tempfile.TemporaryDirectory(prefix="sawf-closed-basis-mpi-") as directory:
    directory = Path(directory)
    source = directory / "src"
    binary = directory / "build"
    source.mkdir()
    (source / "driver.f90").write_text(DRIVER)
    (source / "CMakeLists.txt").write_text(
        f"""cmake_minimum_required(VERSION 3.18)
project(sawf_closed_basis_mpi LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)
add_executable(check_closed_basis_mpi
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf_basis.f90'}
  ${{CMAKE_CURRENT_SOURCE_DIR}}/driver.f90)
target_link_libraries(check_closed_basis_mpi PRIVATE MPI::MPI_Fortran)
if(TARGET LAPACK::LAPACK)
  target_link_libraries(check_closed_basis_mpi PRIVATE LAPACK::LAPACK)
else()
  target_link_libraries(check_closed_basis_mpi PRIVATE ${{LAPACK_LIBRARIES}})
endif()
"""
    )
    result = run(
        [
            "cmake",
            "-S",
            source,
            "-B",
            binary,
            f"-DCMAKE_Fortran_COMPILER={compiler}",
        ]
    )
    if result.returncode:
        raise SystemExit(result.stdout)
    result = run(["cmake", "--build", binary, "-j", "2"])
    if result.returncode:
        raise SystemExit(result.stdout)
    cache_text = (binary / "CMakeCache.txt").read_text()
    match = re.search(r"^MPIEXEC_EXECUTABLE:FILEPATH=(.+)$", cache_text, re.M)
    if not match:
        raise SystemExit("MPIEXEC_EXECUTABLE was not found")
    result = run([match.group(1), "-n", "2", binary / "check_closed_basis_mpi"])
    if result.returncode:
        raise SystemExit(result.stdout)
    output = result.stdout.strip()
    if "[SAWF-CLOSED-BASIS-LOCAL] rank=1" not in output:
        raise SystemExit("rank-local failure diagnostic was not emitted before reduction")
    print(output)
