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

module init_gs
  implicit none

contains

!===================================================================================================================================

SUBROUTINE init_wf(lg,mg,system,info,pp,spsi)
  use structures
  use inputoutput, only: au_length_aa, method_init_wf
  use salmon_global, only: yn_periodic,natom,Rion,yn_jm
  use gram_schmidt_orth
  implicit none

  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_pp_info)         ,intent(in) :: pp
  type(s_orbital)                     :: spsi
  !
  integer :: ik,io,is,a,ix,iy,iz,ip, ig,ngauss
  real(8) :: xx,yy,zz,x1,y1,z1,rr,Xmax,Ymax,Zmax,q(3)
  real(8) :: Xzero,Yzero,Zzero

  call init_wf_rand

  select case(method_init_wf)
  case ('gauss'  ) ; ngauss=1
  case ('gauss2' ) ; ngauss=2
  case ('gauss3' ) ; ngauss=3
  case ('gauss4' ) ; ngauss=4
  case ('gauss5' ) ; ngauss=5
  case ('gauss10') ; ngauss=10
  end select

  ! get offset (0-th element) : Xzero means center
  if (yn_periodic == 'y') then
    Xzero = lg%coordinate(lg%num(1)/2+mod(lg%num(1),2),1)
    Yzero = lg%coordinate(lg%num(2)/2+mod(lg%num(2),2),2)
    Zzero = lg%coordinate(lg%num(3)/2+mod(lg%num(3),2),3)
  else
    Xzero = 0d0
    Yzero = 0d0
    Zzero = 0d0
  end if

  if(system%if_real_orbital) then

    select case(method_init_wf)

    case ('random')
      call gen_random_rwf

    case ('pp')
      call gen_pp_orbitals(lg,mg,system,pp,spsi)

    case ('gauss','gauss2','gauss3','gauss4','gauss5','gauss10')
      if(yn_jm=='n') then
         Xmax=0.d0 ; Ymax=0.d0 ; Zmax=0.d0
         do a=1,natom
           if ( abs(system%Rion(1,a)) > Xmax ) Xmax=abs(system%Rion(1,a))
           if ( abs(system%Rion(2,a)) > Ymax ) Ymax=abs(system%Rion(2,a))
           if ( abs(system%Rion(3,a)) > Zmax ) Zmax=abs(system%Rion(3,a))
         end do

         Xmax=Xmax-Xzero+1.d0/au_length_aa
         Ymax=Ymax-Yzero+1.d0/au_length_aa
         Zmax=Zmax-Zzero+1.d0/au_length_aa
      else
         Xmax=( lg%coordinate(lg%ie(1),1)-Xzero )*0.8d0
         Ymax=( lg%coordinate(lg%ie(2),2)-Yzero )*0.8d0
         Zmax=( lg%coordinate(lg%ie(3),3)-Zzero )*0.8d0
      end if

      do is=1,system%nspin
      do io=1,system%no
      do ig=1,ngauss
        call random_number(q)
        x1=Xmax*(2.d0*q(1)-1.d0)
        y1=Ymax*(2.d0*q(2)-1.d0)
        z1=Zmax*(2.d0*q(3)-1.d0)
        if(info%io_s <= io .and. io <= info%io_e) then
#ifdef USE_OPENACC
!$acc parallel loop collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
#else
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
#endif
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            xx=lg%coordinate(ix,1)-Xzero
            yy=lg%coordinate(iy,2)-Yzero
            zz=lg%coordinate(iz,3)-Zzero
            rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
            if(ig==1) then
              spsi%rwf(ix,iy,iz,is,io,1,1) = exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
            else
              spsi%rwf(ix,iy,iz,is,io,1,1) = spsi%rwf(ix,iy,iz,is,io,1,1) &
                   + exp(-0.5d0*(rr*au_length_aa)**2)*(au_length_aa)**(3/2)
            endif
          end do
          end do
          end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif
        end if
      end do !ig
      end do
      end do

    end select

  else

    select case(method_init_wf)

    case ('random')
      call gen_random_zwf

    case ('pp')
      call gen_pp_orbitals(lg,mg,system,pp,spsi)

    case ('gauss','gauss2','gauss3','gauss4','gauss5','gauss10')
      if(yn_jm=='n') then
        Xmax = sqrt(sum(system%primitive_a(1:3,1)**2))-Xzero
        Ymax = sqrt(sum(system%primitive_a(1:3,2)**2))-Yzero
        Zmax = sqrt(sum(system%primitive_a(1:3,3)**2))-Zzero
      else
         Xmax=( lg%coordinate(lg%ie(1),1)-Xzero )*0.8d0
         Ymax=( lg%coordinate(lg%ie(2),2)-Yzero )*0.8d0
         Zmax=( lg%coordinate(lg%ie(3),3)-Zzero )*0.8d0
      endif

      do is=1,system%nspin
      do ik=1,system%nk
      do io=1,system%no
      do ig=1,ngauss
        call random_number(q)
        x1=Xmax*q(1)
        y1=Ymax*q(2)
        z1=Zmax*q(3)
        if(info%ik_s <= ik .and. ik <= info%ik_e .and.   &
           info%io_s <= io .and. io <= info%io_e) then
!$OMP parallel do collapse(2) private(iz,iy,ix,xx,yy,zz,rr)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            xx=lg%coordinate(ix,1)-Xzero
            yy=lg%coordinate(iy,2)-Yzero
            zz=lg%coordinate(iz,3)-Zzero
            rr=sqrt((xx-x1)**2+(yy-y1)**2+(zz-z1)**2)
            if(ig==1) then
               spsi%zwf(ix,iy,iz,is,io,ik,1) = exp(-0.5d0*rr**2)
            else
               spsi%zwf(ix,iy,iz,is,io,ik,1) = spsi%zwf(ix,iy,iz,is,io,ik,1) + exp(-0.5d0*rr**2)
            endif
          end do
          end do
          end do
!$omp end parallel do
        end if
      end do !ig
      end do
      end do
      end do

    end select

  end if

  call gram_schmidt(system, mg, info, spsi)

  return

CONTAINS

  ! cf. RSDFT
  subroutine init_wf_rand
    use salmon_global, only: iseed_number_change
    implicit none
    integer :: s,k,n,i,llen
    integer,allocatable :: iseed(:)

    call random_seed(size = n)
    allocate(iseed(n))
    llen = product(lg%num)
    iseed(:) = (info%ik_s * system%no + info%io_s - 1) * llen &
             + (mg%is(3) - lg%is(3) + 1) * lg%num(2) * lg%num(1) &
             + (mg%is(2) - lg%is(2) + 1) * lg%num(1) &
             + (mg%is(1) - lg%is(1) + 1) + iseed_number_change
    call random_seed(put = iseed)
    deallocate(iseed)
  end subroutine

  subroutine gen_random_rwf
    implicit none
    real(8) :: u(3),v(3)

    do ip=lbound(spsi%rwf,7),ubound(spsi%rwf,7)
    do ik=lbound(spsi%rwf,6),ubound(spsi%rwf,6)
    do io=lbound(spsi%rwf,5),ubound(spsi%rwf,5)
    do is=lbound(spsi%rwf,4),ubound(spsi%rwf,4)
    do iz=lbound(spsi%rwf,3),ubound(spsi%rwf,3)
    do iy=lbound(spsi%rwf,2),ubound(spsi%rwf,2)
    do ix=lbound(spsi%rwf,1),ubound(spsi%rwf,1)
      v = dble([ix, iy, iz])
      call random_number(u)
      spsi%rwf(ix,iy,iz,is,io,ik,ip) = product(sign(u(1:3),v(1:3)))
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine

  subroutine gen_random_zwf
    implicit none
    real(8) :: u(2)

    do ip=lbound(spsi%zwf,7),ubound(spsi%zwf,7)
    do ik=lbound(spsi%zwf,6),ubound(spsi%zwf,6)
    do io=lbound(spsi%zwf,5),ubound(spsi%zwf,5)
    do is=lbound(spsi%zwf,4),ubound(spsi%zwf,4)
    do iz=lbound(spsi%zwf,3),ubound(spsi%zwf,3)
    do iy=lbound(spsi%zwf,2),ubound(spsi%zwf,2)
    do ix=lbound(spsi%zwf,1),ubound(spsi%zwf,1)
      call random_number(u)
      spsi%zwf(ix,iy,iz,is,io,ik,ip) = dcmplx(u(1),u(2))
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end subroutine

END SUBROUTINE init_wf

!===================================================================================================================================

! Build the table of pseudo-atomic orbitals (index -> atom, l, m, energy), sorted by energy in ascending order.
! The sorted table is also written to dbg_sort_test.txt for debugging.
SUBROUTINE make_pp_orbital_table(pp,n_ao,iatom_ao,l_ao,m_ao,eig_ao)
  use structures, only: s_pp_info
  use salmon_global, only: natom,kion,nstate,nelem
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root,comm_sync_all,comm_bcast
  implicit none
  type(s_pp_info),intent(in) :: pp
  integer,intent(out) :: n_ao
  integer,allocatable,intent(out) :: iatom_ao(:),l_ao(:),m_ao(:)
  real(8),allocatable,intent(out) :: eig_ao(:)
  integer :: ia,ik,l,m,i,j,iao,fh
  integer :: it,lt,mt
  real(8) :: et

  ! Orbitals can be generated only from the eigenvalues of the atomic solver.
  do ik=1,nelem
    if (.not. (pp%has_wf_pp(ik) .and. pp%has_pseudo_wf(ik))) then
      if (comm_is_root(nproc_id_global)) then
        write(*,*) "Error: pseudo atomic orbitals (eigenvalue and occupation) are not available: ik=",ik
        write(*,*) "       Set yn_pseudo_atomic_orbital='y' to use method_init_wf='pp'."
      end if
      call comm_sync_all(nproc_group_global)
      stop "method_init_wf='pp' requires yn_pseudo_atomic_orbital='y'"
    end if
  end do

  ! Each l channel has (2*l+1) spatial orbitals; spin is handled by the spin index of spsi.
  n_ao=0
  do ia=1,natom
    n_ao=n_ao+sum([(2*l+1,l=0,pp%mlps(kion(ia)))])
  end do
  if (n_ao < nstate) then
    if (comm_is_root(nproc_id_global)) then
      write(*,*) "Warning: number of pseudo atomic orbitals is smaller than nstate: n_ao,nstate=",n_ao,nstate
    end if
  end if

  ! Only the root process builds, sorts and writes the table; it is broadcast afterward.
  allocate(iatom_ao(n_ao),l_ao(n_ao),m_ao(n_ao),eig_ao(n_ao))
  if (comm_is_root(nproc_id_global)) then
    iao=0
    do ia=1,natom
      ik=kion(ia)
      do l=0,pp%mlps(ik)
      do m=-l,l
        iao=iao+1
        iatom_ao(iao)=ia
        l_ao(iao)=l
        m_ao(iao)=m
        eig_ao(iao)=pp%pseudo_wf_eigval(l,ik)
      end do
      end do
    end do

    ! Stable insertion sort by energy (ascending).
    do i=2,n_ao
      it=iatom_ao(i); lt=l_ao(i); mt=m_ao(i); et=eig_ao(i)
      j=i-1
      do while (j >= 1)
        if (eig_ao(j) <= et) exit
        iatom_ao(j+1)=iatom_ao(j); l_ao(j+1)=l_ao(j); m_ao(j+1)=m_ao(j); eig_ao(j+1)=eig_ao(j)
        j=j-1
      end do
      iatom_ao(j+1)=it; l_ao(j+1)=lt; m_ao(j+1)=mt; eig_ao(j+1)=et
    end do

    open(newunit=fh,file='method_init_wf_pp.txt',status='replace')
    write(fh,'(a,2i8)') '# n_ao, nstate =',n_ao,nstate
    write(fh,'(a)') '# index  iatom  l   m   energy(Ha)'
    do i=1,n_ao
      write(fh,'(i7,3i5,f20.10)') i,iatom_ao(i),l_ao(i),m_ao(i),eig_ao(i)
    end do
    close(fh)
  end if

  call comm_bcast(iatom_ao,nproc_group_global)
  call comm_bcast(l_ao,nproc_group_global)
  call comm_bcast(m_ao,nproc_group_global)
  call comm_bcast(eig_ao,nproc_group_global)
END SUBROUTINE make_pp_orbital_table

!===================================================================================================================================

! Map the sorted pseudo-atomic orbitals onto the real-space grid of each process.
! Orbital io (io <= n_ao) is the atomic orbital (iatom_ao,l_ao,m_ao)(io), summed over periodic replicas.
! Orbitals with io > n_ao are random. A tiny random noise is added to every atomic orbital only as a trigger
! to break the symmetry of the system; it is kept small because high-frequency noise is costly to damp in CG.
SUBROUTINE gen_pp_orbitals(lg,mg,system,pp,spsi)
  use structures, only: s_rgrid,s_dft_system,s_pp_info,s_orbital
  use salmon_global, only: kion,iperiodic
  use salmon_math, only: Ylm
  use math_constants, only: pi
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_rgrid)     ,intent(in) :: lg,mg
  type(s_dft_system),intent(in) :: system
  type(s_pp_info)   ,intent(in) :: pp
  type(s_orbital)               :: spsi
  integer :: n_ao
  integer,allocatable :: iatom_ao(:),l_ao(:),m_ao(:)
  real(8),allocatable :: eig_ao(:)
  real(8),parameter :: noise_amp=0.0d0
  integer :: io,ik,is,ia,ikion,l,m,nc,irep1,irep2,irep3,i1,i2,i3
  integer :: io_s,io_e,ik_s,ik_e,is_s,is_e
  integer :: j,fh
  real(8) :: rshift(3),hgs(3),u,v,w,x,y,z,r,rinv,c_l,rc(3),q
  real(8) :: phi(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))

  ! orbital table: index -> (atom, l, m), sorted by the energy
  call make_pp_orbital_table(pp,n_ao,iatom_ao,l_ao,m_ao,eig_ao)

  hgs=system%Hgs
  if (iperiodic == 3) then
    rshift=-hgs
    nc=1
  else
    do i1=1,3
      if (mod(lg%num(i1),2) == 1) then
        rshift(i1)=0d0
      else
        rshift(i1)=-0.5d0*hgs(i1)
      end if
    end do
    nc=0
  end if

  if (allocated(spsi%rwf)) then
    is_s=lbound(spsi%rwf,4); is_e=ubound(spsi%rwf,4)
    io_s=lbound(spsi%rwf,5); io_e=ubound(spsi%rwf,5)
    ik_s=lbound(spsi%rwf,6); ik_e=ubound(spsi%rwf,6)
  else
    is_s=lbound(spsi%zwf,4); is_e=ubound(spsi%zwf,4)
    io_s=lbound(spsi%zwf,5); io_e=ubound(spsi%zwf,5)
    ik_s=lbound(spsi%zwf,6); ik_e=ubound(spsi%zwf,6)
  end if

  ! debug output (root process only)
  if (comm_is_root(nproc_id_global)) then
    ! open(newunit=fh,file='method_init_wf_pp.txt',status='old',position='append')
    ! write(fh,'(a)') '# spsi bounds (is_s,is_e,io_s,io_e,ik_s,ik_e)'
    ! write(fh,'(6i8)') is_s,is_e,io_s,io_e,ik_s,ik_e
    ! write(fh,'(a)') '# system%Rion(1:3,ia), ia=1..nion'
    ! do ia=1,size(system%Rion,2)
    !   write(fh,'(i7,3es24.15)') ia,system%Rion(1:3,ia)
    ! end do
    ! write(fh,'(a)') '# system%primitive_a(1:3,j), j=1..3 (columns)'
    ! do j=1,3
    !   write(fh,'(i7,3es24.15)') j,system%primitive_a(1:3,j)
    ! end do
    ! close(fh)
  end if

  do io=io_s,io_e
    phi=0d0
    if (io <= n_ao) then
      ia=iatom_ao(io); ikion=kion(ia); l=l_ao(io); m=m_ao(io)
      ! The radial table is used up to rad(nrmax); rad must be extended monotonically beyond mr.
      if (pp%rad(pp%nrmax,ikion) <= pp%rad(pp%mr(ikion),ikion)) &
        stop "gen_pp_orbitals: radial grid is not extended beyond mr"
      c_l=sqrt((2d0*l+1d0)/(4d0*pi))
      do irep1=-nc,nc
      do irep2=-nc,nc
      do irep3=-nc,nc
        rc(1:3)=system%Rion(1:3,ia)+irep1*system%primitive_a(1:3,1) &
                                   +irep2*system%primitive_a(1:3,2) &
                                   +irep3*system%primitive_a(1:3,3)
!$OMP parallel do collapse(2) default(shared) private(i1,i2,i3,u,v,w,x,y,z,r,rinv)
        do i3=mg%is(3),mg%ie(3)
        do i2=mg%is(2),mg%ie(2)
        do i1=mg%is(1),mg%ie(1)
          u=i1*hgs(1)+rshift(1)
          v=i2*hgs(2)+rshift(2)
          w=i3*hgs(3)+rshift(3)
          x=u*system%rmatrix_A(1,1)+v*system%rmatrix_A(1,2)+w*system%rmatrix_A(1,3)-rc(1)
          y=u*system%rmatrix_A(2,1)+v*system%rmatrix_A(2,2)+w*system%rmatrix_A(2,3)-rc(2)
          z=u*system%rmatrix_A(3,1)+v*system%rmatrix_A(3,2)+w*system%rmatrix_A(3,3)-rc(3)
          r=sqrt(x*x+y*y+z*z)+1d-50
          if (r >= pp%rad(pp%nrmax,ikion)) cycle
          rinv=1d0/r
          phi(i1,i2,i3)=phi(i1,i2,i3)+c_l*radial_orbital(r)*Ylm(x*rinv,y*rinv,z*rinv,l,m)
        end do
        end do
        end do
!$omp end parallel do
      end do
      end do
      end do
      do i3=mg%is(3),mg%ie(3)
      do i2=mg%is(2),mg%ie(2)
      do i1=mg%is(1),mg%ie(1)
        call random_number(q)
        phi(i1,i2,i3)=phi(i1,i2,i3)+noise_amp*(2d0*q-1d0)
      end do
      end do
      end do
    else
      do i3=mg%is(3),mg%ie(3)
      do i2=mg%is(2),mg%ie(2)
      do i1=mg%is(1),mg%ie(1)
        call random_number(q)
        phi(i1,i2,i3)=2d0*q-1d0
      end do
      end do
      end do
    end if

    do ik=ik_s,ik_e
    do is=is_s,is_e
      if (allocated(spsi%rwf)) then
        spsi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,io,ik,1)=phi
      else
        spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),is,io,ik,1)=cmplx(phi,0d0,kind=8)
      end if
    end do
    end do
  end do

CONTAINS

  ! u_l(r)/r from the raw u=r*R table pp%upp_f; upp_f(i) refers to rad(i+1).
  ! The caller guarantees r < rad(nrmax); the first interval (r < rad(2)) is extended by a constant u/r.
  function radial_orbital(r) result(f)
    use salmon_math, only: interp_linear
    implicit none
    real(8),intent(in) :: r
    real(8) :: f

    if (r < pp%rad(2,ikion)) then
      f=pp%upp_f(1,l,ikion)/pp%rad(2,ikion)
    else
      f=interp_linear(pp%nrmax-1,pp%rad(2:pp%nrmax,ikion),pp%upp_f(1:pp%nrmax-1,l,ikion),r)/r
    end if
  end function radial_orbital

END SUBROUTINE gen_pp_orbitals

end module init_gs
