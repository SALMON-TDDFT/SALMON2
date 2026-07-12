#!/usr/bin/env python3
"""Run the production SAWF alignment-failure collective on two MPI ranks."""

from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build-mpi-eigenexa-wannier-lib"

COMMUNICATION_STUB = r"""
module communication
  use mpi
  implicit none
  interface comm_get_max
    module procedure comm_get_max_integer
  end interface
contains
  subroutine comm_get_max_integer(value,communicator)
    integer, intent(inout) :: value
    integer, intent(in) :: communicator
    integer :: reduced,ierr
    call MPI_Allreduce(value,reduced,1,MPI_INTEGER,MPI_MAX,communicator,ierr)
    if(ierr/=MPI_SUCCESS) call MPI_Abort(communicator,91,ierr)
    value=reduced
  end subroutine comm_get_max_integer
end module communication
"""

DRIVER = r"""
program check_sawf_fragment_symmetry_map_mpi
  use mpi
  use lcfo_wannier_sawf, only: t_sawf_symop
  use lcfo_wannier_sawf_band, only: validate_sawf_fragment_symmetry_map
  use lcfo_wannier_sawf_collective, only: reduce_sawf_fragment_alignment_failure
  implicit none
  type(t_sawf_symop) :: op
  integer :: mesh(3),origin(3,8),shape(3,8),buffer(3),rank,nrank,ierr
  integer :: ix,iy,iz,ifrag,max_targets,local_failure,global_failure
  integer, allocatable :: source_to_target(:)
  logical :: grid_ok,fragment_ok,center_available
  real(8) :: residual,center(3)
  character(256) :: message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank/=2) call MPI_Abort(MPI_COMM_WORLD,90,ierr)
  mesh=[32,32,32]; shape=16; ifrag=0
  do iz=0,1
    do iy=0,1
      do ix=0,1
        ifrag=ifrag+1
        origin(:,ifrag)=[16*ix,16*iy,16*iz]
      end do
    end do
  end do
  op%W=0; op%W(1,2)=1; op%W(2,1)=1; op%W(3,3)=1; op%tau=0d0
  buffer=[6,6,6]
  if(rank==1) buffer=[6,7,6]
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  local_failure=merge(0,1,grid_ok .and. fragment_ok)
  call reduce_sawf_fragment_alignment_failure(local_failure,MPI_COMM_WORLD,rank,2, &
    grid_ok,fragment_ok,max_targets,residual,message,global_failure)
  if(rank==0 .and. local_failure/=0) call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  if(rank==1 .and. local_failure/=1) call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  if(global_failure/=1) call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  if(rank==0) write(*,'(a)') 'PASS SAWF MPI shared alignment failure reduction'
  call MPI_Finalize(ierr)
end program check_sawf_fragment_symmetry_map_mpi
"""


def run(command, **kwargs):
    return subprocess.run(
        [str(item) for item in command], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, **kwargs
    )


def main():
    cache = (BUILD / "CMakeCache.txt").read_text()
    compiler = next(
        line.split("=", 1)[1] for line in cache.splitlines()
        if line.startswith("CMAKE_Fortran_COMPILER:")
    )
    with tempfile.TemporaryDirectory(prefix="sawf-fragment-map-mpi-") as temporary:
        temporary = Path(temporary)
        source = temporary / "src"
        build = temporary / "build"
        source.mkdir()
        (source / "config.h").write_text("")
        (source / "sym_stub.f90").write_text(
            """module sym_sub
contains
subroutine read_symmetry_file(f,m,o,s)
character(*),intent(in)::f
real(8),allocatable,intent(out)::m(:,:,:)
logical,intent(out)::o
character(*),intent(out)::s
allocate(m(3,4,0)); o=.false.; s='stub'
end subroutine read_symmetry_file
end module sym_sub
"""
        )
        (source / "communication.f90").write_text(COMMUNICATION_STUB)
        (source / "driver.f90").write_text(DRIVER)
        (source / "CMakeLists.txt").write_text(
            f"""cmake_minimum_required(VERSION 3.18)
project(sawf_fragment_map_mpi LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
find_package(MPI REQUIRED COMPONENTS Fortran)
add_library(sawf
  {source / 'sym_stub.f90'}
  {source / 'communication.f90'}
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf.f90'}
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf_band.f90'}
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf_collective.f90'})
target_include_directories(sawf PRIVATE {source})
target_link_libraries(sawf PUBLIC MPI::MPI_Fortran)
if(TARGET LAPACK::LAPACK)
  target_link_libraries(sawf PUBLIC LAPACK::LAPACK)
else()
  target_link_libraries(sawf PUBLIC ${{LAPACK_LIBRARIES}})
endif()
add_executable(check_fragment_map_mpi {source / 'driver.f90'})
target_link_libraries(check_fragment_map_mpi PRIVATE sawf MPI::MPI_Fortran)
"""
        )
        result = run(["cmake", "-S", source, "-B", build, f"-DCMAKE_Fortran_COMPILER={compiler}"])
        if result.returncode:
            raise SystemExit(result.stdout)
        result = run(["cmake", "--build", build, "-j", "2"])
        if result.returncode:
            raise SystemExit(result.stdout)
        mpi_match = re.search(r"^MPIEXEC_EXECUTABLE:FILEPATH=(.+)$", (build / "CMakeCache.txt").read_text(), re.M)
        if not mpi_match:
            raise SystemExit("MPIEXEC_EXECUTABLE was not found")
        result = run([mpi_match.group(1), "-n", "2", build / "check_fragment_map_mpi"])
        if result.returncode:
            raise SystemExit(result.stdout)
        assert result.stdout.count("[DC-LCFO-SAWF-ALIGN] rank=1 operation=2") == 1
        assert "PASS SAWF MPI shared alignment failure reduction" in result.stdout
        print(result.stdout.strip())

    shared = (ROOT / "src/gs/dc/lcfo_wannier_sawf_collective.f90").read_text()
    diagnostic = shared.index("[DC-LCFO-SAWF-ALIGN]")
    collective = shared.index("call comm_get_max")
    assert diagnostic < collective
    flux = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()
    assert "call reduce_sawf_fragment_alignment_failure" in flux
    print("PASS production and MPI test share alignment failure orchestration")


if __name__ == "__main__":
    main()
