#!/usr/bin/env python3
"""Compile the CheFSI module against minimal typed dependency interfaces."""

from pathlib import Path
import shutil
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
compiler = shutil.which("gfortran")
if compiler is None:
    raise SystemExit("gfortran is required")

stubs = r"""
module structures
  implicit none
  type :: s_dcdft
    integer :: i_frag=1,icomm_frag=0,icomm_tot=0,id_frag=0,id_tot=0
    integer :: isize_frag=1,isize_tot=1,n_frag=1,nstate_tot=1
  end type
end module

module communication
  implicit none
  interface comm_get_max
    module procedure comm_get_max_integer
    module procedure comm_get_max_array1d_double
  end interface
  interface comm_summation
    module procedure comm_summation_double
    module procedure comm_summation_array1d_double
    module procedure comm_summation_array2d_double
    module procedure comm_summation_array1d_integer
    module procedure comm_summation_array2d_integer
  end interface
contains
  subroutine comm_bcast(value,comm,root)
    class(*), intent(inout) :: value(..)
    integer, intent(in) :: comm,root
  end subroutine
  integer function comm_create_group(comm,color,key)
    integer, intent(in) :: comm,color,key
    comm_create_group=comm
  end function
  subroutine comm_free_group(comm)
    integer, intent(inout) :: comm
  end subroutine
  subroutine comm_get_max_integer(value,comm)
    integer, intent(inout) :: value
    integer, intent(in) :: comm
  end subroutine
  subroutine comm_get_max_array1d_double(input,output,n,comm)
    real(8), intent(in) :: input(:)
    real(8), intent(out) :: output(:)
    integer, intent(in) :: n,comm
  end subroutine
  integer function comm_isend(value,dest,tag,comm)
    class(*), intent(in) :: value(..)
    integer, intent(in) :: dest,tag,comm
    comm_isend=0
  end function
  integer function comm_irecv(value,source,tag,comm)
    class(*) :: value(..)
    integer, intent(in) :: source,tag,comm
    comm_irecv=0
  end function
  subroutine comm_summation_double(input,output,comm)
    real(8), intent(in) :: input
    real(8), intent(out) :: output
    integer, intent(in) :: comm
  end subroutine
  subroutine comm_summation_array1d_double(input,output,n,comm)
    real(8), intent(in) :: input(:)
    real(8), intent(out) :: output(:)
    integer, intent(in) :: n,comm
  end subroutine
  subroutine comm_summation_array2d_double(input,output,n,comm)
    real(8), intent(in) :: input(:,:)
    real(8), intent(out) :: output(:,:)
    integer, intent(in) :: n,comm
  end subroutine
  subroutine comm_summation_array1d_integer(input,output,n,comm)
    integer, intent(in) :: input(:)
    integer, intent(out) :: output(:)
    integer, intent(in) :: n,comm
  end subroutine
  subroutine comm_summation_array2d_integer(input,output,n,comm)
    integer, intent(in) :: input(:,:)
    integer, intent(out) :: output(:,:)
    integer, intent(in) :: n,comm
  end subroutine
  subroutine comm_wait_all(request)
    integer, intent(inout) :: request(:)
  end subroutine
end module

module eigen_subdiag_sub
  implicit none
contains
  subroutine eigen_dsyev(matrix,eigenvalue,eigenvector)
    real(8), intent(inout) :: matrix(:,:)
    real(8), intent(out) :: eigenvalue(:),eigenvector(:,:)
  end subroutine
end module

module timer
  implicit none
  integer, parameter :: LOG_CHEFSI_SETUP=1,LOG_CHEFSI_TOTAL=2
  integer, parameter :: LOG_CHEFSI_H_APPLY=3,LOG_CHEFSI_H_COMM_POST=4
  integer, parameter :: LOG_CHEFSI_H_DIAG=5,LOG_CHEFSI_H_HALO=6
  integer, parameter :: LOG_CHEFSI_H_RECV_WAIT=7,LOG_CHEFSI_H_SEND_WAIT=8
  integer, parameter :: LOG_CHEFSI_FILTER=9,LOG_CHEFSI_ORTHO=10
  integer, parameter :: LOG_CHEFSI_PROJECT=11,LOG_CHEFSI_PROJECT_EIGEN=12
  integer, parameter :: LOG_CHEFSI_RAYLEIGH_RITZ=13,LOG_CHEFSI_REDISTRIBUTE=14
  integer, parameter :: LOG_CHEFSI_RESIDUAL=15,LOG_CHEFSI_ROTATE=16
  integer, parameter :: LOG_CHEFSI_LANCZOS=17,LOG_CHEFSI_EXPORT=18
contains
  subroutine timer_begin(id)
    integer, intent(in) :: id
  end subroutine
  subroutine timer_end(id)
    integer, intent(in) :: id
  end subroutine
end module
"""

with tempfile.TemporaryDirectory(prefix="lcfo-chefsi-interface-") as name:
    build = Path(name)
    stub_path = build / "stubs.f90"
    stub_path.write_text(stubs)
    subprocess.run(
        [
            compiler,
            "-std=f2018",
            "-ffree-line-length-none",
            "-J",
            str(build),
            "-I",
            str(build),
            "-c",
            str(stub_path),
            str(ROOT / "src/gs/dc/lcfo_diag_chefsi.f90"),
        ],
        cwd=build,
        check=True,
    )

print("PASS CheFSI requested-state interface compiles")
