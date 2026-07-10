#!/usr/bin/env python3
from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build-mpi-eigenexa-wannier-lib"

DRIVER = r"""
program check_sawf_closed_basis
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use lcfo_wannier_sawf_basis, only: close_sawf_candidate_basis
  implicit none
  real(8), allocatable :: basis(:,:), singular_values(:),candidate_transform(:,:)
  real(8) :: candidate(4,4), duplicate(4,4), identity_candidate(4,2)
  real(8) :: gram(4,4), mapped(4,4), overlap(4,4), nan_value
  integer :: nbasis, i
  logical :: ok
  character(512) :: message

  candidate=0d0
  candidate(:,1)=[1d0,1d0,0d0,0d0]/sqrt(2d0)
  candidate(:,2)=[0d0,0d0,1d0,0d0]
  candidate(:,3)=[0d0,0d0,1d0,1d0]/sqrt(2d0)
  candidate(:,4)=[0d0,1d0,0d0,0d0]
  call close_sawf_candidate_basis(candidate,4,4,1d0,1d-12,4, &
    basis,nbasis,singular_values,ok,message,candidate_transform=candidate_transform)
  call require(ok,'full inversion orbit closure failed: '//trim(message))
  call require(nbasis==4,'full inversion orbit rank')
  call require(maxval(abs(basis-matmul(candidate,candidate_transform)))<1d-13, &
    'published candidate-to-closed transform')
  gram=matmul(transpose(basis),basis)
  call require(maxval(abs(gram-identity4()))<1d-12,'closed basis orthonormality')
  mapped=basis([4,3,2,1],:)
  overlap=matmul(transpose(basis),mapped)
  call require(maxval(abs(matmul(transpose(overlap),overlap)-identity4()))<1d-12, &
    'closed basis inversion representation')

  duplicate=0d0
  duplicate(1,1)=1d0; duplicate(2,2)=1d0
  duplicate(:,3)=duplicate(:,1); duplicate(:,4)=duplicate(:,2)
  call close_sawf_candidate_basis(duplicate,4,4,1d0,1d-12,4, &
    basis,nbasis,singular_values,ok,message)
  call require(ok .and. nbasis==2,'dependent orbit pruning')
  call require(all(singular_values(1:2)>1d0),'descending retained singular values')

  identity_candidate=0d0
  identity_candidate(1,1)=1d0; identity_candidate(2,2)=1d0
  call close_sawf_candidate_basis(identity_candidate,4,2,1d0,1d-12,2, &
    basis,nbasis,singular_values,ok,message)
  call require(ok .and. nbasis==2,'identity basis rank')
  call require(maxval(abs(matmul(transpose(basis),basis)-identity2()))<1d-12, &
    'identity basis unchanged span')

  call close_sawf_candidate_basis(candidate,4,4,1d0,1d-12,3, &
    basis,nbasis,singular_values,ok,message)
  call require(.not.ok .and. index(message,'capacity')>0,'capacity rejection')
  call require(.not.allocated(basis) .and. .not.allocated(singular_values) .and. nbasis==0, &
    'capacity failure must not publish partial output')

  nan_value=ieee_value(0d0,ieee_quiet_nan)
  candidate(1,1)=nan_value
  call close_sawf_candidate_basis(candidate,4,4,1d0,1d-12,4, &
    basis,nbasis,singular_values,ok,message)
  call require(.not.ok .and. index(message,'non-finite')>0,'non-finite rejection')
  call require(.not.allocated(basis) .and. nbasis==0,'non-finite transactional failure')

  write(*,'(a)') 'PASS SAWF rank-revealing symmetry-orbit closure'
contains
  function identity4() result(a)
    real(8) :: a(4,4)
    a=0d0
    do i=1,4; a(i,i)=1d0; end do
  end function identity4
  function identity2() result(a)
    real(8) :: a(2,2)
    a=0d0; a(1,1)=1d0; a(2,2)=1d0
  end function identity2
  subroutine require(condition,text)
    logical,intent(in) :: condition
    character(*),intent(in) :: text
    if(.not.condition) then
      write(*,'(a)') trim(text)
      error stop 1
    end if
  end subroutine require
end program check_sawf_closed_basis
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

with tempfile.TemporaryDirectory(prefix="sawf-closed-basis-") as directory:
    directory = Path(directory)
    source = directory / "src"
    binary = directory / "build"
    source.mkdir()
    (source / "driver.f90").write_text(DRIVER)
    (source / "CMakeLists.txt").write_text(
        f"""cmake_minimum_required(VERSION 3.18)
project(sawf_closed_basis LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check_closed_basis
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf_basis.f90'}
  ${{CMAKE_CURRENT_SOURCE_DIR}}/driver.f90)
if(TARGET LAPACK::LAPACK)
  target_link_libraries(check_closed_basis PRIVATE LAPACK::LAPACK)
else()
  target_link_libraries(check_closed_basis PRIVATE ${{LAPACK_LIBRARIES}})
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
    result = run([binary / "check_closed_basis"])
    if result.returncode:
        raise SystemExit(result.stdout)
    print(result.stdout.strip())
