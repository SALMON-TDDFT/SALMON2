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
module mixing_sub
  implicit none

contains

!===================================================================================================================================
subroutine simple_mixing(mg,system,c1,c2,rho_s,mixing)
  use structures, only: s_rgrid, s_dft_system, s_scalar, s_mixing  
  implicit none
  type(s_rgrid),intent(in) :: mg
  type(s_dft_system),intent(in) :: system
  real(8),intent(in) :: c1,c2
  type(s_scalar),intent(inout) :: rho_s(system%nspin)
  type(s_mixing),intent(inout) :: mixing
  
  integer :: ix,iy,iz
  
  if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  elseif(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
      mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  !rho = c1*rho + c2*matmul( psi**2, occ )
  if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(1)%f(ix,iy,iz) = c1*mixing%rho_in(mixing%num_rho_stock)%f(ix,iy,iz) &
                              + c2*mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)
      mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
  else if(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      rho_s(1)%f(ix,iy,iz) = c1*mixing%rho_s_in(mixing%num_rho_stock,1)%f(ix,iy,iz) &
                              + c2*mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)
      rho_s(2)%f(ix,iy,iz) = c1*mixing%rho_s_in(mixing%num_rho_stock,2)%f(ix,iy,iz) &
                              + c2*mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)
      mixing%rho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz) = rho_s(1)%f(ix,iy,iz)
      mixing%rho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz) = rho_s(2)%f(ix,iy,iz)
    end do
    end do
    end do
  end if
  
  
  return
  
end subroutine simple_mixing

!===================================================================================================================================

subroutine wrapper_broyden(comm,mg,system,rho_s,iter,mixing)
  use structures, only: s_rgrid,s_dft_system,s_scalar,s_mixing
  use broyden_sub
  implicit none
  integer,intent(in) :: comm
  type(s_rgrid) :: mg
  type(s_dft_system),intent(in) :: system
  type(s_scalar),intent(inout) :: rho_s(system%nspin)
  integer,intent(in) :: iter
  type(s_mixing),intent(inout) :: mixing
  integer :: ix,iy,iz,is
  integer :: i
  real(8) :: vecr(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
  real(8) :: vecr_in(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),mixing%num_rho_stock+1)
  real(8) :: vecr_out(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),mixing%num_rho_stock+1)

  if(system%nspin==1)then

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$omp parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       vecr(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do

#ifdef USE_OPENACC
!$acc parallel loop private(i,iz,iy,ix) collapse(3)
#else
!$omp parallel do private(i,iz,iy,ix) collapse(3)
#endif
    do i=1,mixing%num_rho_stock+1
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          vecr_in(ix,iy,iz,i) =mixing%rho_in(i)%f(ix,iy,iz)
          vecr_out(ix,iy,iz,i)=mixing%rho_out(i)%f(ix,iy,iz)
       end do
       end do
       end do
    end do

    call broyden(mixing%alpha_mb,vecr,vecr_in,vecr_out,mg%num(1)*mg%num(2)*mg%num(3),iter,    &
                 mixing%num_rho_stock,mixing%num_rho_stock,comm,&
                 mixing%flag_mix_zero)

#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix) collapse(2)
#else
!$omp parallel do private(iz,iy,ix) collapse(2)
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
       rho_s(1)%f(ix,iy,iz)= vecr(ix,iy,iz)
    end do
    end do
    end do

#ifdef USE_OPENACC
!$acc parallel loop private(i,iz,iy,ix) collapse(3)
#else
!$omp parallel do private(i,iz,iy,ix) collapse(3)
#endif
    do i=1,mixing%num_rho_stock+1
       do iz=mg%is(3),mg%ie(3)
       do iy=mg%is(2),mg%ie(2)
       do ix=mg%is(1),mg%ie(1)
          mixing%rho_in(i)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
          mixing%rho_out(i)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
       end do
       end do
       end do
    end do

  else if(system%nspin==2)then
    
    do is=1,2
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         vecr(ix,iy,iz)=rho_s(is)%f(ix,iy,iz)
      end do
      end do
      end do
  
!$omp parallel do private(i,iz,iy,ix) collapse(3)
      do i=1,mixing%num_rho_stock+1
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           vecr_in(ix,iy,iz,i)=mixing%rho_s_in(i,is)%f(ix,iy,iz)
           vecr_out(ix,iy,iz,i)=mixing%rho_s_out(i,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do

      call broyden(mixing%alpha_mb,vecr,vecr_in, vecr_out, mg%num(1)*mg%num(2)*mg%num(3),iter,  &
                   mixing%num_rho_stock,mixing%num_rho_stock,comm,&
                   mixing%flag_mix_zero )

!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
         rho_s(is)%f(ix,iy,iz)= vecr(ix,iy,iz)
      end do
      end do
      end do

!$omp parallel do private(i,iz,iy,ix) collapse(3)
      do i=1,mixing%num_rho_stock+1
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           mixing%rho_s_in(i,is)%f(ix,iy,iz)=vecr_in(ix,iy,iz,i)
           mixing%rho_s_out(i,is)%f(ix,iy,iz)=vecr_out(ix,iy,iz,i)
        end do
        end do
        end do
      end do
    end do

  end if

end subroutine wrapper_broyden

!===================================================================================================================================

subroutine pulay(mg,info,system,rho_s,iter,mixing)
  use salmon_global, only: nmemory_p
  use structures, only: s_rgrid,s_parallel_info,s_dft_system,s_scalar,s_mixing,allocate_scalar,deallocate_scalar
  use communication, only: comm_summation
  implicit none
  type(s_rgrid),intent(in)            :: mg
  type(s_parallel_info),intent(in) :: info
  type(s_dft_system),intent(in)       :: system
  type(s_scalar),intent(inout)        :: rho_s(system%nspin)
  integer,intent(in)                  :: iter
  type(s_mixing),intent(inout)        :: mixing
  integer :: nsize
  integer, allocatable :: ipiv(:)
  type(s_scalar) :: x,y
  real(8), allocatable :: b1(:)
  real(8), allocatable :: a1(:,:)
  real(8), allocatable :: a0(:,:)
  integer :: i,j,i0,j0
  integer :: is,ix,iy,iz
  integer :: ierr
  real(8) :: ss
  real(8) :: rc

  if(iter==1.or.nmemory_p==1)then

    call simple_mixing(mg,system,1.d0-mixing%beta_p,mixing%beta_p,rho_s,mixing)

  else
!pulay mixing

    if(system%nspin == 1)then
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_out(mixing%num_rho_stock)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
      end do
      end do
      end do
    elseif(system%nspin == 2)then
!$omp parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_s_out(mixing%num_rho_stock,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
        mixing%rho_s_out(mixing%num_rho_stock,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
      end do
      end do
      end do
    end if
  
  !rho = c1*rho + c2*matmul( psi**2, occ )

    if(iter>=nmemory_p)then
      nsize=nmemory_p
    else
      nsize=iter
    end if
  
    allocate( ipiv(nsize) )
    allocate( b1(nsize) )
    allocate( a1(nsize,nsize) )
    allocate( a0(nsize,nsize) )
    call allocate_scalar(mg,x)
    call allocate_scalar(mg,y)
  
  
    b1(:)   = 0.d0
    a1(:,:) = 0.d0
    a0(:,:) = 0.d0
  
    do j0=1 ,nsize
    do i0=j0,nsize
      i=mixing%num_rho_stock-nsize+i0
      j=mixing%num_rho_stock-nsize+j0
      ss=0.d0
      if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2) reduction(+:ss)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          ss=ss+(mixing%rho_out(i)%f(ix,iy,iz)-mixing%rho_in(i)%f(ix,iy,iz))* &
                (mixing%rho_out(j)%f(ix,iy,iz)-mixing%rho_in(j)%f(ix,iy,iz))
        end do
        end do
        end do
      else
!$omp parallel do private(ix,iy,iz) collapse(3) reduction(+:ss)
        do is=1,system%nspin
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            ss=ss+(mixing%rho_s_out(i,is)%f(ix,iy,iz)-mixing%rho_s_in(i,is)%f(ix,iy,iz))* &
                  (mixing%rho_s_out(j,is)%f(ix,iy,iz)-mixing%rho_s_in(j,is)%f(ix,iy,iz))
          end do
          end do
          end do
        end do
      end if
      a0(i0,j0)=ss
      a0(j0,i0)=ss
    end do
    end do

    call comm_summation(a0,a1,nsize*nsize,info%icomm_r)

    b1(1:nsize) = 1.d0
  
    call dgesv(nsize,1,a1,nsize,ipiv,b1,nsize,ierr)

    rc=1.d0/sum( b1(1:nsize) )
    b1(1:nsize)=rc*b1(1:nsize)
 
    do is=1,system%nspin
!$omp parallel do private(ix,iy,iz) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        x%f(ix,iy,iz)=0.d0
        y%f(ix,iy,iz)=0.d0
      end do
      end do
      end do
  
      do i0=1,nsize
        i=mixing%num_rho_stock-nsize+i0
        if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            x%f(ix,iy,iz)=x%f(ix,iy,iz)+b1(i0)*mixing%rho_in(i)%f(ix,iy,iz)
            y%f(ix,iy,iz)=y%f(ix,iy,iz)+b1(i0)*mixing%rho_out(i)%f(ix,iy,iz)
          end do
          end do
          end do
        else
!$omp parallel do private(ix,iy,iz) collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            x%f(ix,iy,iz)=x%f(ix,iy,iz)+b1(i0)*mixing%rho_s_in(i,is)%f(ix,iy,iz)
            y%f(ix,iy,iz)=y%f(ix,iy,iz)+b1(i0)*mixing%rho_s_out(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
        end if
      end do
  
      if(system%nspin==1)then
!$omp parallel do private(ix,iy,iz) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz) = max(1.d-20,  &
                                                                   x%f(ix,iy,iz) + mixing%beta_p*( y%f(ix,iy,iz)-x%f(ix,iy,iz) ))
          rho_s(is)%f(ix,iy,iz) = mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)
        end do
        end do
        end do
      else
!$omp parallel do private(ix,iy,iz) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_in(mixing%num_rho_stock+1,is)%f(ix,iy,iz) = max(1.d-20,  &
                                                                     x%f(ix,iy,iz) + mixing%beta_p*( y%f(ix,iy,iz)-x%f(ix,iy,iz) ))
          rho_s(is)%f(ix,iy,iz) = mixing%rho_s_in(mixing%num_rho_stock+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end if
  
    end do

    deallocate( a0,a1,b1,ipiv )
    call deallocate_scalar(x)
    call deallocate_scalar(y)

  end if

end subroutine

!===================================================================================================================================

subroutine init_mixing(nspin,mg,mixing)
  use salmon_global, only: mixrate,alpha_mb,beta_p
  use structures
  implicit none
  integer      ,intent(in) :: nspin
  type(s_rgrid),intent(in) :: mg
  type(s_mixing)           :: mixing
  !
  integer :: i,j

  mixing%mixrate=mixrate
  mixing%alpha_mb=alpha_mb
  mixing%beta_p=beta_p
  mixing%convergence_value_prev=1.d10

  allocate(mixing%rho_in(1:mixing%num_rho_stock+1))
  allocate(mixing%rho_out(1:mixing%num_rho_stock+1))
  do i=1,mixing%num_rho_stock+1
    allocate(mixing%rho_in(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(mixing%rho_out(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    mixing%rho_in(i)%f(:,:,:) =0.d0
    mixing%rho_out(i)%f(:,:,:)=0.d0
  end do

  if(nspin==2)then
    allocate(mixing%rho_s_in(1:mixing%num_rho_stock+1,2))
    allocate(mixing%rho_s_out(1:mixing%num_rho_stock+1,2))
    do j=1,2
      do i=1,mixing%num_rho_stock+1
        allocate(mixing%rho_s_in(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
        allocate(mixing%rho_s_out(i,j)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
        mixing%rho_s_in(i,j)%f(:,:,:) =0.d0
        mixing%rho_s_out(i,j)%f(:,:,:)=0.d0
      end do
    end do
  end if

  mixing%flag_mix_zero=.false.

end subroutine init_mixing

!===================================================================================================================================

subroutine copy_density(Miter,nspin,mg,rho_s,mixing)
  use structures, only: s_rgrid, s_scalar, s_mixing
  implicit none
  integer       ,intent(in) :: Miter,nspin
  type(s_rgrid), intent(in) :: mg
  type(s_scalar),intent(in) :: rho_s(nspin)
  type(s_mixing),intent(inout) :: mixing
  !
  integer :: iiter
  integer :: is
  integer :: ix,iy,iz

  if(Miter==1)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
    end do
    end do
    end do
    if(nspin==2)then
!$OMP parallel do private(iz,iy,ix) collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_s_in(mixing%num_rho_stock+1,1)%f(ix,iy,iz)=rho_s(1)%f(ix,iy,iz)
        mixing%rho_s_in(mixing%num_rho_stock+1,2)%f(ix,iy,iz)=rho_s(2)%f(ix,iy,iz)
      end do
      end do
      end do
    end if
  end if

  do iiter=1,mixing%num_rho_stock
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_in(iiter)%f(ix,iy,iz)=mixing%rho_in(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do
  do iiter=1,mixing%num_rho_stock-1
!$OMP parallel do private(iz,iy,ix) collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      mixing%rho_out(iiter)%f(ix,iy,iz)=mixing%rho_out(iiter+1)%f(ix,iy,iz)
    end do
    end do
    end do
  end do

  if(nspin==2)then
    do iiter=1,mixing%num_rho_stock
      do is=1,2
!$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_in(iiter,is)%f(ix,iy,iz)=mixing%rho_s_in(iiter+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
    do iiter=1,mixing%num_rho_stock-1
      do is=1,2
!$OMP parallel do private(iz,iy,ix) collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          mixing%rho_s_out(iiter,is)%f(ix,iy,iz)=mixing%rho_s_out(iiter+1,is)%f(ix,iy,iz)
        end do
        end do
        end do
      end do
    end do
  end if

end subroutine copy_density

!===================================================================================================================================
subroutine check_mixing_half(Miter,convergence_value,mixing)
  use salmon_global, only: method_mixing,update_mixing_ratio
  use structures, only: s_mixing
  use parallelization, only: nproc_id_global, nproc_group_global
  use communication, only: comm_is_root, comm_bcast
  implicit none
  integer, intent(in) :: Miter
  real(8), intent(in) :: convergence_value
  type(s_mixing), intent(inout) :: mixing
  integer :: icheck
  
  if(comm_is_root(nproc_id_global)) then
    if(convergence_value > update_mixing_ratio * mixing%convergence_value_prev)then
      icheck=1
    else
      icheck=0
    end if
  end if

  call comm_bcast(icheck,nproc_group_global)

  if(icheck==1)then
    select case(method_mixing)
    case('simple')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" mixrate decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%mixrate, mixing%mixrate*0.5d0, Miter
      end if
      mixing%mixrate=mixing%mixrate*0.5d0
    case('broyden')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" alpha_mb decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%alpha_mb, mixing%alpha_mb*0.5d0, Miter 
      end if
      mixing%alpha_mb=mixing%alpha_mb*0.5d0
    case('pulay')
      if(comm_is_root(nproc_id_global)) then
        write(*,'(" beta_p decreased from",e16.8," to",e16.8," at iter = ", i6,"." )')  &
          mixing%beta_p, mixing%beta_p*0.5d0, Miter 
      end if
      mixing%beta_p=mixing%beta_p*0.5d0
    end select
  end if

  mixing%convergence_value_prev=convergence_value

end subroutine check_mixing_half

subroutine reset_mixing_rate(mixing)
  use salmon_global, only: mixrate, alpha_mb, beta_p
  use structures, only: s_mixing
  implicit none
  type(s_mixing), intent(inout) :: mixing

  mixing%mixrate  = mixrate
  mixing%alpha_mb = alpha_mb
  mixing%beta_p   = beta_p
  mixing%convergence_value_prev = 1.d10

end subroutine reset_mixing_rate


!===================================================================================================================================
end module mixing_sub
