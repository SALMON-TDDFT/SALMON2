!
!  Copyright 2017-2020 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
module broyden_sub

contains

subroutine broyden(alpha_mb,vecr,vecr_in,vecr_out,nl,iter,iter_mod,nstock,icomm,flag_mix_zero)
  use salmon_global, only: nmemory_mb
  use communication, only: comm_summation
  use salmon_math
  implicit none
  real(8), intent(in) :: alpha_mb
  integer, intent(in) :: nl
  integer, intent(in) :: iter
  integer, intent(in) :: iter_mod
  integer, intent(in) :: nstock
  integer, intent(in), optional  :: icomm
  real(8), intent(inout) :: vecr(1:nl)
  real(8), intent(inout) :: vecr_in(1:nl,1:nstock+1)
  real(8), intent(inout) :: vecr_out(1:nl,1:nstock+1)
  logical, intent(in),optional :: flag_mix_zero
  integer,parameter :: iter_mb=0
  real(8),parameter :: omega0=0.01d0
  integer :: iter_s,iter_e
  integer :: i,j,k
  real(8),allocatable :: vecf(:,:)
  real(8) :: vecr_tmp(1:nl)
  real(8),allocatable :: del_vecf(:,:),del_vecx(:,:)
  real(8),allocatable :: omega_mb(:)
  real(8),allocatable :: aa(:,:),aa_tmp1(:,:)
  real(8),allocatable :: beta(:,:)
  real(8),allocatable :: ss_tmp1(:), ss_tmp2(:)
  real(8) :: ss, amix
  integer :: nnegative,nnegative_tmp

  amix = alpha_mb
  if(present(flag_mix_zero)) then
     if(flag_mix_zero) amix=0d0
  endif

#ifdef USE_OPENACC
!$acc parallel loop private(i)
#else
!$omp parallel do private(i)
#endif
  do i=1,nl
    vecr_out(i,iter_mod)=vecr(i)
  end do

  if (iter <= iter_mb+1) then
    allocate(vecf(1:nl,iter:iter))
#ifdef USE_OPENACC
!$acc parallel loop private(i)
#else
!$omp parallel do private(i)
#endif
    do i=1,nl
      vecf(i,iter) = vecr_out(i,iter_mod) - vecr_in(i,iter_mod)
    end do
#ifdef USE_OPENACC
!$acc parallel loop private(i)
#else
!$omp parallel do private(i)
#endif
    do i=1,nl
      vecr_in(i,iter_mod+1) = vecr_in(i,iter_mod) + amix*vecf(i,iter)
    end do
    deallocate(vecf)
  else
    iter_s=max(iter_mb+1+(iter_mod-iter),iter_mod-nmemory_mb)
    iter_e=iter_mod-1
    allocate(vecf(1:nl,iter_s:iter_e+1))
    allocate(del_vecf(1:nl,iter_s:iter_e))
    allocate(del_vecx(1:nl,iter_s:iter_e))
    allocate(omega_mb(iter_s:iter_e))
    allocate(beta(iter_s:iter_e,iter_s:iter_e))
    allocate(aa(1:iter_e-iter_s+1,1:iter_e-iter_s+1))
    allocate(aa_tmp1(1:iter_e-iter_s+1,1:iter_e-iter_s+1))
    if(present(icomm)) then
      allocate(ss_tmp1(iter_s:iter_e))
      allocate(ss_tmp2(iter_s:iter_e))
    end if

    omega_mb(:) = 1.d0 !???

!$omp parallel do private(i,j) collapse(2)
    do i=iter_s,iter_mod
    do j=1,nl
      vecf(j,i) = vecr_out(j,i) - vecr_in(j,i)
    end do
    end do

!$omp parallel do private(i) collapse(2)
    do i=iter_s,iter_e
    do j=1,nl
      del_vecx(j,i)=vecr_in(j,i+1)-vecr_in(j,i)
      del_vecf(j,i)=vecf   (j,i+1)-vecf   (j,i)
    end do
    end do

    if(present(icomm)) then
!$omp parallel do private(i,j,ss)
      do i=iter_s,iter_e
        ss = 0d0
        do j=1,nl
          ss = ss + del_vecf(j,i)**2
        end do
        ss_tmp1(i) = ss
      end do
      call comm_summation(ss_tmp1,ss_tmp2,iter_e-iter_s+1,icomm)

!$omp parallel do private(i,j,ss) collapse(2)
      do i=iter_s,iter_e
      do j=1,nl
        ss = sqrt(ss_tmp2(i))
        del_vecx(j,i)=del_vecx(j,i)/ss
        del_vecf(j,i)=del_vecf(j,i)/ss
      end do
      end do
    else 
      do i=iter_s,iter_e
        ss=sum(del_vecf(1:nl,i)**2)
        del_vecx(1:nl,i)=del_vecx(1:nl,i)/ss
        del_vecf(1:nl,i)=del_vecf(1:nl,i)/ss
      end do
    end if

    if(present(icomm)) then
!$omp parallel do collapse(2) private(i,j,k,ss)
      do i=1,iter_e-iter_s+1
      do j=1,iter_e-iter_s+1
        ss = 0d0
        do k=1,nl
          ss = ss + del_vecf(k,iter_s-1+i)*del_vecf(k,iter_s-1+j)
        end do
        aa_tmp1(i,j) = ss
      end do
      end do
      call comm_summation(aa_tmp1,aa,(iter_e-iter_s+1)**2,icomm)

!$omp parallel do collapse(2) private(i,j)
      do i=1,iter_e-iter_s+1
      do j=1,iter_e-iter_s+1
        aa(i,j)=omega_mb(iter_s-1+i)*omega_mb(iter_s-1+j)*aa(i,j)
        if(i==j)then
          aa(i,j)=aa(i,j)+omega0**2
        end if
      end do
      end do
    else
      do i=1,iter_e-iter_s+1
        do j=1,iter_e-iter_s+1
          aa(i,j)=omega_mb(iter_s-1+i)*omega_mb(iter_s-1+j)*sum(del_vecf(1:nl,iter_s-1+i)*del_vecf(1:nl,iter_s-1+j))
          if(i==j)then
            aa(i,j)=aa(i,j)+omega0**2
          end if
        end do
      end do
    end if

    call matrix_inverse(aa)

!$omp parallel do private(i,j) collapse(2)
    do i=iter_s,iter_e
    do j=iter_s,iter_e
      beta(i,j) = aa(i-iter_s+1, j-iter_s+1)
    end do
    end do

    if(present(icomm)) then
!$omp parallel do private(i,j,ss)
      do i=iter_s,iter_e
        ss = 0d0
        do j=1,nl
          ss = ss + del_vecf(j,i)*vecf(j,iter_mod)
        end do
        ss_tmp1(i) = ss
      end do
      call comm_summation(ss_tmp1,ss_tmp2,iter_e-iter_s+1,icomm)

      vecr_tmp(1:nl)=0.d0
!$omp parallel do private(k,i,j,ss)
      do k=1,nl
        ss = 0d0
        do j=iter_s,iter_e
        do i=iter_s,iter_e
          ss = ss + omega_mb(i) * omega_mb(j) * beta(i,j) * ss_tmp2(i) * (amix * del_vecf(k,j) + del_vecx(k,j))
        end do
        end do
        vecr_tmp(k) = ss
      end do
    else 
      vecr_tmp(1:nl)=0.d0
      do i=iter_s,iter_e
      do j=iter_s,iter_e
         vecr_tmp(1:nl)=vecr_tmp(1:nl)&
              &+omega_mb(i)*omega_mb(j)*beta(i,j)*sum(del_vecf(:,i)*vecf(:,iter_mod))*(amix*del_vecf(1:nl,j)+del_vecx(1:nl,j))
      end do
      end do
    end if

!$omp parallel do private(i)
    do i=1,nl
      vecr_in(i,iter_mod+1)=vecr_in(i,iter_mod)+amix*vecf(i,iter_mod)-vecr_tmp(i)
    end do

    if(present(icomm)) then
      nnegative_tmp=0
!$omp parallel do private(i) reduction(+:nnegative_tmp)
      do i=1,nl
        if(vecr_in(i,iter_mod+1) < 0.d0) then
          nnegative_tmp=nnegative_tmp+1
        end if
      end do
      call comm_summation(nnegative_tmp,nnegative,icomm)
    else
      nnegative=0
      do i=1,nl
        if(vecr_in(i,iter_mod+1) < 0.d0) then
          nnegative=nnegative+1
        end if
      end do
    end if
    if(nnegative > 0) then
!$omp parallel do private(i)
      do i=1,nl
        vecr_in(i,iter_mod+1)=vecr_in(i,iter_mod)+amix*vecf(i,iter_mod)
      end do
    end if

    deallocate(vecf,del_vecf,del_vecx,omega_mb,beta,aa,aa_tmp1)
    if(present(icomm)) then
      deallocate(ss_tmp1,ss_tmp2)
    end if
  end if

#ifdef USE_OPENACC
!$acc parallel loop private(i)
#else
!$omp parallel do private(i)
#endif
  do i=1,nl
    vecr(i)=vecr_in(i,iter_mod+1)
  end do

  return
end subroutine broyden

end module broyden_sub
