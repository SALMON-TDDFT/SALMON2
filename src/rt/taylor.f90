!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module taylor_sub
  implicit none
#ifdef USE_OPENACC
  private :: flat_set, flat_axpy
#endif

contains

subroutine taylor(mg,system,info,stencil,srg,tspsi_in,tspsi_out,sshtpsi,   &
                  ppg,V_local,rt)
  use inputoutput, only: n_hamil
  use structures
  use hamiltonian, only: hpsi
  use sendrecv_grid, only: s_sendrecv_grid
  use nvtx_wrapper
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_stencil),intent(in) :: stencil
  type(s_sendrecv_grid),intent(inout) :: srg
  type(s_orbital),intent(inout) :: tspsi_in
  type(s_orbital),intent(inout) :: tspsi_out
  type(s_orbital),intent(inout) :: sshtpsi
  type(s_pp_grid),intent(in) :: ppg
  type(s_scalar) ,intent(in) :: V_local(system%nspin)
  type(s_rt),     intent(in) :: rt
  integer :: nn,ix,iy,iz
  integer :: ik,io,is,nspin
#ifdef USE_OPENACC
  logical :: flat
  integer(8) :: n
  integer :: x0,y0,z0,o0,k0
#endif
  call nvtxStartRange('taylor', __LINE__)

  nspin = system%nspin
#ifdef USE_OPENACC
  ! Undivided grid: the arrays hold no halo, so update the whole im=1 slab in one flat loop.
  flat = .not. info%if_divide_rspace
  n = size(tspsi_out%zwf, kind=8) / size(tspsi_out%zwf, 7, kind=8)
  x0 = mg%is_array(1); y0 = mg%is_array(2); z0 = mg%is_array(3); o0 = info%io_s; k0 = info%ik_s
#endif

  do nn=1,n_hamil
    if(mod(nn,2)==1)then
      call hpsi(tspsi_in,sshtpsi,info,mg,V_local,system,stencil,srg,ppg)
#ifdef USE_OPENACC
      if (flat) then
        if (nn==1) then
          call flat_set(n, rt%zc(nn), tspsi_in%zwf(x0,y0,z0,1,o0,k0,1), sshtpsi%zwf(x0,y0,z0,1,o0,k0,1), &
                        tspsi_out%zwf(x0,y0,z0,1,o0,k0,1))
        else
          call flat_axpy(n, rt%zc(nn), sshtpsi%zwf(x0,y0,z0,1,o0,k0,1), tspsi_out%zwf(x0,y0,z0,1,o0,k0,1))
        end if
      else if (nn==1) then
#else
      if (nn==1) then
#endif
#ifdef USE_OPENACC
!$acc kernels
!$acc loop collapse(6) private(ik,io,is,iz,iy,ix)
#else
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
#endif
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          do is=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_in%zwf(ix,iy,iz,is,io,ik,1)+ &
                                                   rt%zc(nn)*sshtpsi%zwf(ix,iy,iz,is,io,ik,1)
            end do
            end do
            end do
          end do
        end do
        end do
#ifdef USE_OPENACC
!$acc end kernels
#endif
      else
#ifdef USE_OPENACC
!$acc parallel loop collapse(6) private(ik,io,is,iz,iy,ix)
#else
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
#endif
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,info%io_e
          do is=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                                   rt%zc(nn)*sshtpsi%zwf(ix,iy,iz,is,io,ik,1)
            end do
            end do
            end do
          end do
        end do
        end do
#ifdef USE_OPENACC
!$acc end parallel
#endif
      end if
    else
      call hpsi(sshtpsi,tspsi_in,info,mg,V_local,system,stencil,srg,ppg)
#ifdef USE_OPENACC
      if (flat) then
        call flat_axpy(n, rt%zc(nn), tspsi_in%zwf(x0,y0,z0,1,o0,k0,1), tspsi_out%zwf(x0,y0,z0,1,o0,k0,1))
      else
!$acc parallel loop collapse(6) private(ik,io,is,iz,iy,ix)
#else
!$OMP parallel do collapse(5) private(ik,io,is,iz,iy,ix)
#endif
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
        do is=1,nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            tspsi_out%zwf(ix,iy,iz,is,io,ik,1)=tspsi_out%zwf(ix,iy,iz,is,io,ik,1)+ &
                                               rt%zc(nn)*tspsi_in%zwf(ix,iy,iz,is,io,ik,1)
          end do
          end do
          end do
        end do
      end do
      end do
#ifdef USE_OPENACC
!$acc end parallel
      end if
#endif
    end if
  end do

  call nvtxEndRange
end subroutine taylor

#ifdef USE_OPENACC
! y = a + z*b over the n contiguous elements of an im=1 slab, passed by its first element;
! the padding is zero in every operand and stays zero.
subroutine flat_set(n,z,a,b,y)
  implicit none
  integer(8),intent(in) :: n
  complex(8),intent(in) :: z,a(n),b(n)
  complex(8),intent(out) :: y(n)
  integer(8) :: i
!$acc parallel loop copyin(a,b) copyout(y)
  do i=1,n
    y(i) = a(i) + z*b(i)
  end do
end subroutine flat_set

! y = y + z*b over n contiguous elements.
subroutine flat_axpy(n,z,b,y)
  implicit none
  integer(8),intent(in) :: n
  complex(8),intent(in) :: z,b(n)
  complex(8),intent(inout) :: y(n)
  integer(8) :: i
!$acc parallel loop copyin(b) copy(y)
  do i=1,n
    y(i) = y(i) + z*b(i)
  end do
end subroutine flat_axpy
#endif

end module taylor_sub

