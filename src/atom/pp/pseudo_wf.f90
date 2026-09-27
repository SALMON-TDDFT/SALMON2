!
!  Copyright 2018-2026 SALMON developers
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
module pseudo_wf
  implicit none

  private
  public :: calc_pseudo_wf

  real(8),parameter :: dr_wf=0.03d0   ! spacing of the uniform radial grid for the solver (Bohr)
  real(8),parameter :: rmax_wf=30.0d0 ! outer boundary of the solver grid (Bohr)
  integer,parameter :: nr_wf=nint(rmax_wf/dr_wf)-1 ! # of inner grid points r(j)=j*dr_wf, j=1..nr_wf

contains

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
subroutine calc_pseudo_wf(pp,ik)
  use structures,only : s_pp_info
  use math_constants,only : pi
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  integer :: n,nloc,nnl,i,j,ll,p,p0,nb,info
  real(8) :: const,e
  real(8),allocatable :: r(:),vloc(:),beta(:,:),hmat(:,:),u(:)
  real(8),allocatable :: rwf(:,:),eigval(:),occup(:)
  integer,allocatable :: sgn(:)

  n=nr_wf
  allocate(r(n),vloc(n),hmat(n,n),u(n))
  allocate(rwf(n,0:pp%mlps(ik)),eigval(0:pp%mlps(ik)),occup(0:pp%mlps(ik)))
  do j=1,n
    r(j)=dr_wf*j
  end do

! range of valid vloctbl and udvtbl (see making_ps_with(out)_masking)
  if (pp%has_proj_pp(ik)) then
    nloc=pp%mr(ik)
    nnl=min(pp%mr(ik)+1,pp%nrmax)
  else
    nloc=pp%nrps(ik)
    nnl=pp%nrps(ik)
  end if

! local potential; -Z/r outside the table
  do j=1,n
    if (r(j) <= pp%rad(nloc,ik)) then
      vloc(j)=interp_linear(nloc,pp%rad(1:nloc,ik),pp%vloctbl(1:nloc,ik),r(j))
    else
      vloc(j)=-pp%zps(ik)/r(j)
    end if
  end do

  p0=0
  do ll=0,pp%mlps(ik)
    const=sqrt((2*ll+1)/(4*pi))

! radial projectors b_p(r) of this angular momentum
    nb=pp%nproj(ll,ik)
    allocate(beta(n,max(nb,1)),sgn(max(nb,1)))
    beta=0d0; sgn=0
    do p=1,nb
      sgn(p)=pp%inorm(p0+p-1,ik)
      if (sgn(p) == 0) cycle
      do j=1,n
        if (r(j) <= pp%rad(nnl,ik)) then
          beta(j,p)=interp_linear(nnl,pp%rad(1:nnl,ik),pp%udvtbl(1:nnl,p0+p-1,ik),r(j))*r(j)**(ll+1)/const
        end if
      end do
    end do
    p0=p0+nb

! Hamiltonian matrix for u(r) in Hartree
    hmat=0d0
    do j=1,n
      hmat(j,j)=1d0/dr_wf**2 + 0.5d0*ll*(ll+1)/r(j)**2 + vloc(j)
      if (j < n) then
        hmat(j,j+1)=-0.5d0/dr_wf**2
        hmat(j+1,j)=-0.5d0/dr_wf**2
      end if
    end do
    do p=1,nb
      if (sgn(p) == 0) cycle
      do i=1,n
        do j=1,n
          hmat(j,i)=hmat(j,i)+sgn(p)*beta(j,p)*beta(i,p)*dr_wf
        end do
      end do
    end do
    deallocate(beta,sgn)

    call calc_single_eigenpair(n,hmat,e,u,info)
    if (info /= 0) then
      write(*,*) "calc_pseudo_wf: diagonalization failed: ik,l,info=",ik,ll,info
      stop "calc_pseudo_wf"
    end if

! normalize int u^2 dr = 1 and fix the sign so that the largest component is positive
    u=u/sqrt(sum(u**2)*dr_wf)
    if (u(maxloc(abs(u),1)) < 0d0) u=-u

    eigval(ll)=e
    rwf(:,ll)=u(:)
  end do

  call set_occupation(pp%mlps(ik),pp%zps(ik),eigval,occup)

  write(*,*) "calc_pseudo_wf: ik=",ik
  do ll=0,pp%mlps(ik)
    write(*,'(a,i3,a,f14.8,a,f8.4)') "   l=",ll,"  eigval(Ha)=",eigval(ll),"  occup=",occup(ll)
  end do

! Preserve file-provided pseudo wavefunctions; interpolate only for elements without them.
  if (.not. pp%has_wf_pp(ik)) call set_upp_from_rwf(pp,ik,n,r,rwf)

  deallocate(r,vloc,hmat,u,rwf,eigval,occup)
  return


contains 


!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! Fill zps electrons into bound channels (eigval<0) in ascending order of eigval, at most 2(2l+1) per channel.
subroutine set_occupation(mlps1,zps1,eigval,occup)
  implicit none
  integer,intent(in) :: mlps1,zps1
  real(8),intent(in) :: eigval(0:mlps1)
  real(8),intent(out) :: occup(0:mlps1)
  logical :: done(0:mlps1)
  integer :: l,lmin
  real(8) :: nelec_tmp

  occup=0d0
  done=eigval >= 0d0
  nelec_tmp=dble(zps1)
  do while (nelec_tmp >= 1d0 .and. .not. all(done))
    lmin=-1
    do l=0,mlps1
      if (done(l)) cycle
      if (lmin < 0) then
        lmin=l
      else if (eigval(l) < eigval(lmin)) then
        lmin=l
      end if
    end do
    occup(lmin)=min(nelec_tmp,dble(2*(2*lmin+1)))
    nelec_tmp=nelec_tmp-occup(lmin)
    done(lmin)=.true.
  end do
  if (nelec_tmp >= 1d0) write(*,*) "Warning: electrons not assigned to any bound channel in calc_pseudo_wf:",nelec_tmp

  return
end subroutine set_occupation



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! Lowest eigenvalue and eigenvector of the real symmetric matrix a (a is destroyed).
subroutine calc_single_eigenpair(n,a,e,v,info)
  implicit none

  integer, intent(in) :: n
  real(8), intent(inout) :: a(n,n)
  real(8), intent(out) :: e, v(n)
  integer, intent(out) :: info

  integer :: m

  real(8) :: w(n), z(n,1)
  real(8) :: work(10*n)
  integer :: iwork(5*n)
  integer :: ifail(n)

  call dsyevx('V','I','U',n,a,n, &
              0d0,0d0,1,1,0d0, &
              m,w,z,n,work,size(work),iwork,ifail,info)

  if (info == 0 .and. m >= 1) then
     e = w(1)
     v(:) = z(:,1)
  end if

  return
end subroutine calc_single_eigenpair



!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! Linear interpolation of f given on the increasing grid x(1:n); x0 must be in [x(1),x(n)].
function interp_linear(n,x,f,x0) result(y)
  implicit none
  integer,intent(in) :: n
  real(8),intent(in) :: x(n),f(n),x0
  real(8) :: y
  integer :: ilo,ihi,imid

  ilo=1; ihi=n
  do while (ihi-ilo > 1)
    imid=(ilo+ihi)/2
    if (x(imid) <= x0) then
      ilo=imid
    else
      ihi=imid
    end if
  end do
  if (x(ihi) == x(ilo)) then
    y=f(ilo)
  else
    y=f(ilo)+(f(ihi)-f(ilo))*(x0-x(ilo))/(x(ihi)-x(ilo))
  end if

  return
end function interp_linear




end subroutine calc_pseudo_wf

!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
! Overwrite the work array pp%upp with the local solver result on the radial mesh of the element.
! As for the input pseudo wavefunctions, upp(i,l) is u_l at rad(i+1,ik) and upp(0,l)=0.
subroutine set_upp_from_rwf(pp,ik,n,r,rwf)
  use structures,only : s_pp_info
  implicit none
  type(s_pp_info),intent(inout) :: pp
  integer,intent(in) :: ik
  integer,intent(in) :: n
  real(8),intent(in) :: r(n),rwf(n,0:pp%mlps(ik))
  integer :: i,ll

  if (pp%has_wf_pp(ik)) return

  pp%upp(:,:)=0d0
  do ll=0,pp%mlps(ik)
    do i=1,min(ubound(pp%upp,1),pp%nrmax-1)
      if (pp%rad(i+1,ik) > r(n)) exit
      pp%upp(i,ll)=interp_linear_origin(n,r,rwf(:,ll),pp%rad(i+1,ik))
    end do
  end do

  return
contains 

function interp_linear_origin(n,r,u,x0) result(y)
  implicit none
  integer,intent(in) :: n
  real(8),intent(in) :: r(n),u(n),x0
  real(8) :: y
  integer :: j

  j=int(x0/dr_wf)
  if (j <= 0) then
    y=u(1)*x0/r(1)
  else if (j >= n) then
    y=u(n)
  else
    y=u(j)+(u(j+1)-u(j))*(x0-r(j))/(r(j+1)-r(j))
  end if

  return
end function interp_linear_origin


end subroutine set_upp_from_rwf





end module pseudo_wf
