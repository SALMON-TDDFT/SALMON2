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
module atomic_solver
  implicit none
  private
  public :: calc_pseudo_wavefunction
  public :: assign_orbital_occupations

  real(8), parameter :: r_max_atom=30d0
  integer, parameter :: nr_atom=1000
  ! Interior points; Dirichlet boundaries are at zero and r_max_atom.
  real(8), parameter :: dr_atom=r_max_atom/(nr_atom+1)
contains


subroutine calc_pseudo_wavefunction(pp,ik,with_masking)
  use structures, only: s_pp_info
  use math_constants, only: pi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  ! Fixed-density, scalar norm-conserving atomic solver (no PAW overlap or SO).
  ! Solve one state per angular momentum, including all its diagonal KB projectors.
  ! Call after making_ps_*; rho_pp_tbl must already contain a valence density.
  ! For generated orbitals, upp follows the post-making_ps_* convention.
  ! upp_f always stores untransformed u=r*R for the PDOS consumer.
  ! Do not subsequently overwrite upp_f with upp.
  type(s_pp_info), intent(inout) :: pp
  integer, intent(in) :: ik
  logical, intent(in), optional :: with_masking
  logical :: masked
  real(8) :: c_l,radius
  real(8) :: r(nr_atom),rho(nr_atom),beta(nr_atom)
  real(8) :: v_eff(nr_atom),v_h(nr_atom),v_xc(nr_atom),v_loc(nr_atom)
  real(8) :: u(nr_atom,0:pp%mlps(ik)),eig(0:pp%mlps(ik))
  integer :: occ(0:pp%mlps(ik))
  real(8) :: norm
  real(8), allocatable :: h(:,:)
  integer :: nloc,nnl,l,p,p0,iproj,i,j,info

  masked=.false.
  if (present(with_masking)) masked=with_masking

  write(*, *) "Recalculating pseudo wavefunction (atomic_solver)"
  call check_input_tables()

  ! Build the uniform mesh and fixed-density effective potential.
  do i=1,nr_atom
    r(i)=i*dr_atom
  end do
  rho=interp_rho_in(pp,ik,nr_atom,r)
  v_loc=interp_local_potential(pp,ik,nloc,nr_atom,r) 
  v_h=calc_hartree_potential(nr_atom,r,rho)
  v_xc= calc_exchange_correlation_potential(nr_atom,rho)
  v_eff= v_loc + v_h + v_xc

  ! Solve each angular channel with all of its KB projectors.
  allocate(h(nr_atom,nr_atom))
  p0=0
  do l=0,pp%mlps(ik)
    h=0d0
    do i=1,nr_atom
      h(i,i)=1d0/dr_atom**2+0.5d0*l*(l+1d0)/r(i)**2+v_eff(i)
      if (i < nr_atom) then
        h(i,i+1)=-0.5d0/dr_atom**2
        h(i+1,i)=-0.5d0/dr_atom**2
      end if
    end do
    do iproj=1,pp%nproj(l,ik)
      p=p0+iproj-1
      if (pp%inorm(p,ik) == 0) cycle
      beta=interp_projector(pp,ik,l,p,nnl,nr_atom,r)
      do j=1,nr_atom
        do i=1,nr_atom
          h(i,j)=h(i,j)+pp%inorm(p,ik)*beta(i)*beta(j)*dr_atom
        end do
      end do
    end do
    p0=p0+pp%nproj(l,ik)
    call calc_single_eigenpair(nr_atom,h,eig(l),u(:,l),info)
    if (info /= 0) then
      write(*,*) "atomic_solver: diagonalization failed: ik,l,info=",ik,l,info
      error stop "atomic_solver: diagonalization failed"
    end if
    norm=sqrt(sum(u(:,l)**2)*dr_atom)
    if (.not. ieee_is_finite(norm) .or. norm <= 0d0) error stop "atomic_solver: invalid wavefunction norm"
    u(:,l)=u(:,l)/norm
    if (u(maxloc(abs(u(:,l)),1),l) < 0d0) u(:,l)=-u(:,l)
  end do

  occ = assign_orbital_occupations(pp%mlps(ik), eig, int(pp%zps(ik)))
  do l = 0, pp%mlps(ik)
    write(*,'(2x,a,i4,a,i3,a,f16.9,a,i)') "ik=",ik," l=",l," eigval(Ha)=",eig(l), " occup=",occ(l)
  end do

  if (.not. pp%has_wf_pp(ik)) then
    ! First store raw u=r*R on the PP mesh, then apply the upp convention.
    pp%upp=0d0
    pp%upp_f(:,:,ik)=0d0
    do l=0,pp%mlps(ik)
      pp%upp_f(:,l,ik)=interp_radial_wavefunction(pp,ik,nr_atom,r,u(:,l))
      pp%upp(:,l)=pp%upp_f(:,l,ik)
    end do

    if (.not. masked) then
      do l=0,pp%mlps(ik)
        c_l=sqrt((2d0*l+1d0)/(4d0*pi))
        do i=1,min(ubound(pp%upp,1),pp%nrmax-1)
          radius=pp%rad(i+1,ik)
          if (radius > r_max_atom) exit
          pp%upp(i,l)=pp%upp(i,l)*c_l/radius**(l+1)
        end do
      end do
    end if
    pp%upp(0,:)=pp%upp(1,:)

    pp%has_wf_pp(ik) = .true.
  end if

  ! Existing density, derivative and DFT+U orbital tables are not rebuilt.
contains

  subroutine check_input_tables()
    implicit none
    if (.not. pp%has_rho_pp(ik)) error stop "atomic_solver: valence density is required"
    nloc=pp%nrps(ik)
    if (pp%has_proj_pp(ik)) nloc=pp%mr(ik)
    ! Use the same nonlocal domain and coordinates as prep_pp/calc_uv.
    nnl=pp%nrps(ik)
    if (nloc < 2 .or. nloc > pp%nrmax .or. nnl < 2 .or. nnl > pp%nrmax) &
      error stop "atomic_solver: invalid potential table range"
    if (pp%mr(ik) < 2 .or. pp%mr(ik) > pp%nrmax) &
      error stop "atomic_solver: invalid density table range"
    if (any(pp%nproj(0:pp%mlps(ik),ik) < 0)) error stop "atomic_solver: negative projector count"
    if (sum(pp%nproj(0:pp%mlps(ik),ik)) > size(pp%udvtbl,2)) &
      error stop "atomic_solver: projector count exceeds table capacity"
  end subroutine check_input_tables

end subroutine calc_pseudo_wavefunction

function interp_radial_wavefunction(pp,ik,nr,r,u) result(upp_raw)
  use structures, only: s_pp_info
  use salmon_math, only: interp_linear
  implicit none
  ! Interpolate raw u=r*R, including its zero Dirichlet boundaries.
  ! Work-array index i corresponds to rad(i+1,ik); unused entries stay zero.
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik,nr
  real(8), intent(in) :: r(nr),u(nr)
  real(8) :: upp_raw(0:ubound(pp%upp,1))
  real(8) :: x(0:nr+1),wf(0:nr+1),radius
  integer :: i

  x(0)=0d0
  x(1:nr)=r
  x(nr+1)=r_max_atom
  wf(0)=0d0
  wf(1:nr)=u
  wf(nr+1)=0d0

  upp_raw=0d0
  do i=1,min(ubound(pp%upp,1),pp%nrmax-1)
    radius=pp%rad(i+1,ik)
    if (radius > r_max_atom) exit
    upp_raw(i)=interp_linear(nr+2,x,wf,radius)
  end do
end function interp_radial_wavefunction

function interp_local_potential(pp,ik,nt,nr,r) result(v)
  use structures, only: s_pp_info
  use salmon_math, only: interp_linear
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik,nt,nr
  real(8), intent(in) :: r(nr)
  real(8) :: v(nr)
  integer :: i
  do i=1,nr
    if (r(i) <= pp%rad(1,ik)) then
      v(i)=pp%vloctbl(1,ik)
    else if (r(i) <= pp%rad(nt,ik)) then
      v(i)=interp_linear(nt,pp%rad(1:nt,ik),pp%vloctbl(1:nt,ik),r(i))
    else
      v(i)=-pp%zps(ik)/r(i)
    end if
  end do
end function interp_local_potential

function interp_rho_in(pp,ik,nr,r) result(rho)
  use structures, only: s_pp_info
  use math_constants, only: pi
  use salmon_math, only: interp_linear
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik,nr
  real(8), intent(in) :: r(nr)
  real(8) :: rho(nr),tmp
  integer :: i,nt
  nt=pp%mr(ik)
  if (any(.not. ieee_is_finite(pp%rho_pp_tbl(1:nt,ik))) .or. &
      any(pp%rho_pp_tbl(1:nt,ik) < 0d0)) error stop "atomic_solver: invalid valence density"
  do i=1,nr
    tmp=0d0
    if (r(i) <= pp%rad(1,ik)) then
      tmp=pp%rho_pp_tbl(1,ik)
    else if (r(i) <= pp%rad(nt,ik)) then
      tmp=interp_linear(nt,pp%rad(1:nt,ik),pp%rho_pp_tbl(1:nt,ik),r(i))
    end if
    rho(i)=tmp/(4d0*pi*r(i)**2)
  end do
end function interp_rho_in

function calc_hartree_potential(nr,r,rho) result(v)
  use math_constants, only: pi
  implicit none
  integer, intent(in) :: nr
  real(8), intent(in) :: r(nr),rho(nr)
  real(8) :: v(nr),q(nr),s(nr),dr
  integer :: i
  ! Cumulative trapezoidal integrals; rho*r**2 is zero at the origin.
  q(1)=0.5d0*r(1)**3*rho(1)
  do i=2,nr
    dr=r(i)-r(i-1)
    q(i)=q(i-1)+0.5d0*dr*(rho(i-1)*r(i-1)**2+rho(i)*r(i)**2)
  end do
  ! The density is truncated at the last solver point.
  s(nr)=0d0
  do i=nr-1,1,-1
    dr=r(i+1)-r(i)
    s(i)=s(i+1)+0.5d0*dr*(rho(i)*r(i)+rho(i+1)*r(i+1))
  end do
  v=4d0*pi*(q/r+s)
end function calc_hartree_potential

function calc_exchange_correlation_potential(nr,rho) result(v)
  use builtin_pz, only: exc_cor_pz
  implicit none
  integer, intent(in) :: nr
  real(8), intent(in) :: rho(nr)
  real(8) :: v(nr),exc(nr),eexc(nr)

  ! rho is the total density; exc_cor_pz expects one spin component
  ! of an unpolarized system and returns the XC potential in vexc.
  call exc_cor_pz(nr,0.5d0*rho,exc,eexc,v)
end function calc_exchange_correlation_potential

function interp_projector(pp,ik,l,p,nt,nr,r) result(beta)
  use structures, only: s_pp_info
  use math_constants, only: pi
  use salmon_math, only: interp_linear
  implicit none
  type(s_pp_info), intent(in) :: pp
  integer, intent(in) :: ik,l,p,nt,nr
  real(8), intent(in) :: r(nr)
  real(8) :: beta(nr),tmp,c_l
  integer :: i
  c_l=sqrt((2d0*l+1d0)/(4d0*pi))
  do i=1,nr
    tmp=0d0
    if (r(i) <= pp%radnl(1,ik)) then
      tmp=pp%udvtbl(1,p,ik)
    else if (r(i) <= pp%radnl(nt,ik)) then
      tmp=interp_linear(nt,pp%radnl(1:nt,ik),pp%udvtbl(1:nt,p,ik),r(i))
    end if
    ! udvtbl already includes sqrt(abs(KB coefficient)); only its sign remains.
    beta(i)=tmp*r(i)**(l+1)/c_l
  end do
end function interp_projector

function assign_orbital_occupations(lmax,eps,nelec) result(occ)
  implicit none
  integer, intent(in) :: lmax,nelec
  real(8), intent(in) :: eps(0:lmax)
  integer :: occ(0:lmax)
  integer :: remaining,l,lmin

  occ=0
  remaining=nelec
  do while (remaining > 0)
    lmin=-1
    do l=0,lmax
      if (occ(l) > 0) cycle
      if (lmin < 0) then
        lmin=l
      else if (eps(l) < eps(lmin)) then
        lmin=l
      end if
    end do
    occ(lmin)=min(remaining,2*(lmin+1))
    remaining=remaining-occ(lmin)
  end do
end function assign_orbital_occupations

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
  else if (info == 0) then
     info = -1
  end if

  return
end subroutine calc_single_eigenpair


end module atomic_solver
