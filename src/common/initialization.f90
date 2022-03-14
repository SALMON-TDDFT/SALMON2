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

#include "config.h"

module initialization_sub
  implicit none
  integer,parameter,private :: Nd=4

contains

!===================================================================================================================================

subroutine init_dft(comm,info,lg,mg,system,stencil,fg,poisson,srg,srg_scalar,ofile)
  use structures
  use salmon_global, only: iperiodic,layout_multipole, &
                           nproc_k,nproc_ob,nproc_rgrid,method_poisson
  use sendrecv_grid
  use init_communicator
  use init_poisson_sub
  use checkpoint_restart_sub, only: init_dir_out_restart
  use sym_rho_sub, only: init_sym_rho
  implicit none
  integer      ,intent(in) :: comm
  type(s_parallel_info)    :: info
  type(s_rgrid)            :: lg,mg
  type(s_dft_system)       :: system
  type(s_stencil)          :: stencil
  type(s_reciprocal_grid)  :: fg
  type(s_poisson)          :: poisson
  type(s_sendrecv_grid)    :: srg,srg_scalar
  type(s_ofile)            :: ofile
  !
  integer,dimension(2,3) :: neig

! electron system
  call init_dft_system(lg,system,stencil)

! process distribution
  info%npk       = nproc_k
  info%nporbital = nproc_ob
  info%nprgrid   = nproc_rgrid
  call init_process_distribution(system,comm,info)
  call init_communicator_dft(comm,info)

! parallelization
  call check_ffte_condition(info,lg)
  call init_grid_parallel(info,lg,mg) ! lg --> mg
  call init_parallel_dft(system,info)
  call create_sendrecv_neig(neig, info) ! neighboring node array
  ! sendrecv_grid object for wavefunction updates
  call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)
  ! sendrecv_grid object for scalar field updates
  call init_sendrecv_grid(srg_scalar, mg, 1, info%icomm_rko, neig)

! symmetry

  call init_sym_rho( lg%num, mg%is, mg%ie, info%icomm_r )

! for Poisson equation

  poisson%iterVh = 0 ! Iteration counter
  select case(iperiodic)
  case(0)
    if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(lg,mg,system,poisson)
    if(method_poisson=='ft') then
      call init_reciprocal_grid_isolated_ft(lg,mg,fg,system,info,poisson)
    end if
  case(3)
    call init_reciprocal_grid(lg,mg,fg,system,info,poisson)
  end select
  call set_ig_bound(lg,mg,poisson)

! output files

  call init_dir_out_restart(ofile)

end subroutine init_dft

!===================================================================================================================================

subroutine init_dft_system(lg,system,stencil)
  use structures
  use lattice
  use salmon_global, only: al_vec1,al_vec2,al_vec3,al,spin,natom,nelem,nstate,iperiodic,num_kgrid,num_rgrid,dl, &
  & nproc_rgrid,Rion,Rion_red,nelec,calc_mode,temperature,nelec_spin,yn_spinorbit, &
  & iflag_atom_coor,ntype_atom_coor_reduced,quiet
  use sym_sub, only: init_sym_sub
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  implicit none
  type(s_rgrid)      :: lg
  type(s_dft_system) :: system
  type(s_stencil)    :: stencil
  !
  integer :: ii,jj
  real(8) :: rsize(3),hgs(3),cnmat(0:12,12),bnmat(4,4)

  if(al_vec1(2)==0d0 .and. al_vec1(3)==0d0 .and. al_vec2(1)==0d0 .and. &
     al_vec2(3)==0d0 .and. al_vec3(1)==0d0 .and. al_vec3(2)==0d0) then
    stencil%if_orthogonal = .true.
    if(al(1)*al(2)*al(3)==0d0) then
      if(num_rgrid(1)*num_rgrid(2)*num_rgrid(3)==0 .or. dl(1)*dl(2)*dl(3)==0d0) then
        stop "error: invalid cell"
      end if
      al = dl * dble(num_rgrid)
    end if
    system%primitive_a = 0d0
    system%primitive_a(1,1) = al(1)
    system%primitive_a(2,2) = al(2)
    system%primitive_a(3,3) = al(3)
    rsize = al
  else
    stencil%if_orthogonal = .false.
    system%primitive_a(1:3,1) = al_vec1
    system%primitive_a(1:3,2) = al_vec2
    system%primitive_a(1:3,3) = al_vec3
    rsize(1) = sqrt(sum(al_vec1**2))
    rsize(2) = sqrt(sum(al_vec2**2))
    rsize(3) = sqrt(sum(al_vec3**2))
  end if

  if(sum(abs(dl)) == 0d0) then
    hgs(1:3) = rsize(1:3) / dble(num_rgrid(1:3))
  else
    hgs(1:3) = dl(1:3)
  end if
  call init_grid_whole(rsize,hgs,lg)
  system%hgs = hgs
  system%ngrid = lg%num(1) * lg%num(2) * lg%num(3)

  call init_lattice(system,stencil)
  call init_sym_sub( system%primitive_a, system%primitive_b )
  call init_kvector(num_kgrid,system)

  if(calc_mode=='RT') then
    system%if_real_orbital = .false.
  else if(calc_mode=='GS') then
    select case(iperiodic)
    case(0)
      system%if_real_orbital = .true.
    case(3)
      if(num_kgrid(1)*num_kgrid(2)*num_kgrid(3)==1 .and. stencil%if_orthogonal) then
        system%if_real_orbital = .true.
      else
        system%if_real_orbital = .false.
      end if
    end select
    if ( yn_spinorbit=='y' ) system%if_real_orbital=.false.
  end if
  if ((.not. quiet) .and. comm_is_root(nproc_id_global)) then
     write(*,*) "  use of real value orbitals = ", system%if_real_orbital
  endif

  system%nion = natom

  if(spin=='unpolarized') then
    system%nspin=1
  else !if(spin=='polarized') then
    system%nspin=2
  end if

  if(calc_mode=='RT'.and. temperature<-1.d-12)then
    if( yn_spinorbit=='y' )then
       system%no = nelec
    else if(system%nspin==2.and.sum(nelec_spin(:))>0)then
       system%no = maxval(nelec_spin(:))
    else
       if(mod(nelec,2)==0)then
          system%no = nelec/2
       else
          system%no = (nelec+1)/2
       end if
    end if
  else
    system%no = nstate
  end if

  allocate(system%mass(1:nelem))
  if ( allocated(system%Rion) ) deallocate(system%Rion)
  if ( allocated(system%rocc) ) deallocate(system%rocc)
  if ( allocated(system%Velocity) ) deallocate(system%Velocity)
  if ( allocated(system%Force) ) deallocate(system%Force)
  allocate(system%Rion(3,system%nion),system%rocc(system%no,system%nk,system%nspin))
  allocate(system%Velocity(3,system%nion),system%Force(3,system%nion))
  system%Velocity(:,:) =0d0

  if(iflag_atom_coor==ntype_atom_coor_reduced) then
    Rion = matmul(system%primitive_a,Rion_red) ! [ a1, a2, a3 ] * R_ion
  end if
  system%Rion = Rion

! initial value of occupation
  system%rocc = 0d0
  if ( yn_spinorbit=='y' ) then
    system%rocc(1:nelec,:,:) = 1.0d0
  else
    select case(system%nspin)
    case(1)
      system%rocc(1:nelec/2,:,1) = 2d0
    case(2)
      if ( nelec > 0 ) then
        if ( mod(nelec,2) == 0 ) then
          system%rocc(1:nelec/2,:,1:2) = 1d0
        else
          system%rocc(1:(nelec-1)/2,:,1:2) = 1d0
          system%rocc((nelec-1)/2+1,:,1  ) = 1d0
        end if
      else if ( any(nelec_spin>0) ) then
        system%rocc(1:nelec_spin(1),:,1) = 1d0
        system%rocc(1:nelec_spin(2),:,2) = 1d0
      else
        write(*,*) "nelect or nelec_spin should be specified in input"
      end if
    end select
  end if

  call set_bn(bnmat)
  call set_cn(cnmat)
  if(stencil%if_orthogonal) then
    stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  else
    if(nproc_rgrid(1)*nproc_rgrid(2)*nproc_rgrid(3)/=1) &
      stop "error: nonorthogonal lattice and r-space parallelization"
    stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*  &
                      & ( stencil%coef_F(1)/Hgs(1)**2 + stencil%coef_F(2)/Hgs(2)**2 + stencil%coef_F(3)/Hgs(3)**2 )
  end if
  do jj=1,3
    do ii=1,4
      stencil%coef_lap(ii,jj) = cnmat(ii,4)/hgs(jj)**2
      stencil%coef_nab(ii,jj) = bnmat(ii,4)/hgs(jj)
    end do
  end do

  call set_gridcoordinate(lg,system)

  system%vec_Ac = 0d0 ! initial value

  system%vec_Ac_ext = 0d0 ! initial value for output
  system%vec_E = 0d0     ! initial value for output
  system%vec_E_ext = 0d0 ! initial value for output

  return
end subroutine init_dft_system

!===================================================================================================================================

subroutine init_process_distribution(system,icomm1,info)
  use structures, only: s_parallel_info,s_dft_system
  use parallelization, only: nproc_id_global, nproc_group_global
  use salmon_global, only: theory
  use communication, only: comm_is_root,comm_bcast
  use set_numcpu
  implicit none
  type(s_dft_system),intent(in)       :: system
  integer, intent(in)                 :: icomm1 ! Communicator for single DFT system.
  type(s_parallel_info),intent(inout) :: info
  logical :: if_stop

  if((info%nporbital + sum(info%nprgrid)) == 0) then
    ! Process distribution is automatically decided by SALMON.
    if (system%ngrid > 16**3) then
      call set_numcpu_general(iprefer_domain_distribution,system%nk,system%no,icomm1,info)
    else
      select case(theory)
      case('dft','dft_band','dft_md','dft2tddft')
        call set_numcpu_general(iprefer_k_distribution,system%nk,system%no,icomm1,info)
      case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multi_scale_maxwell_tddft')
        call set_numcpu_general(iprefer_orbital_distribution,system%nk,system%no,icomm1,info)
      case default
        stop 'invalid theory @ initialization'
      end select
    end if
  else
    ! Process distribution is explicitly specified by user.
  end if

  if (comm_is_root(nproc_id_global)) then
    if_stop = .not. check_numcpu(icomm1, info)
  end if
  call comm_bcast(if_stop, nproc_group_global)
  if (if_stop) stop 'fail: check_numcpu'

#ifdef USE_SCALAPACK
  info%flag_blacs_gridinit = .false.
#endif
#ifdef USE_EIGENEXA
  info%flag_eigenexa_init = .false.
#endif
end subroutine init_process_distribution

!===================================================================================================================================

subroutine init_parallel_dft(system,info)
  use structures
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info)         :: info
  !
  integer :: io,nproc_k,nproc_ob,nproc_domain_orbital(3),m,na

  nproc_k              = info%npk
  nproc_ob             = info%nporbital
  nproc_domain_orbital = info%nprgrid

! for single-cell calculations
  info%im_s = 1
  info%im_e = 1
  info%numm = 1

! # of k points
  info%ik_s = (info%id_k * system%nk) / nproc_k + 1
  info%ik_e = ((info%id_k+1) * system%nk) / nproc_k
  info%numk = info%ik_e - info%ik_s + 1

! # of orbitals
  info%io_s = (info%id_o * system%no) / nproc_ob + 1
  info%io_e = ((info%id_o+1) * system%no) / nproc_ob
  info%numo = info%io_e - info%io_s + 1

  allocate(info%io_s_all(0:nproc_ob-1))
  allocate(info%io_e_all(0:nproc_ob-1))
  allocate(info%numo_all(0:nproc_ob-1))
  do m=0,nproc_ob-1
    info%io_s_all(m) = m * system%no / nproc_ob + 1
    info%io_e_all(m) = (m+1) * system%no / nproc_ob
    info%numo_all(m) = info%io_e_all(m) - info%io_s_all(m) + 1
  end do
  info%numo_max=maxval(info%numo_all(0:nproc_ob-1))

! flags
  info%if_divide_rspace = nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3).ne.1
  info%if_divide_orbit  = nproc_ob.ne.1

! process ID corresponding to the orbital index io
  if ( allocated(info%irank_io) ) deallocate(info%irank_io)
  allocate(info%irank_io(1:system%no))
  do io=1, system%no
    if(mod(io*nproc_ob,system%no)==0)then
      info%irank_io(io) = io*nproc_ob/system%no - 1
    else
      info%irank_io(io) = io*nproc_ob/system%no
    end if
  end do

! #ia: atom index (communicator=info%icomm_ko)
  info%ia_s = int((system%nion * info%id_ko) / info%isize_ko) + 1
  info%ia_e = int((system%nion * (info%id_ko + 1)) / info%isize_ko)
  na = info%ia_e - info%ia_s + 1
  if (info%id_ko == info%isize_ko-1) info%ia_e = system%nion

end subroutine init_parallel_dft

!===================================================================================================================================

subroutine init_grid_whole(rsize,hgs,lg)
  use structures, only: s_rgrid
  use salmon_global, only: iperiodic,dl,num_rgrid,theory,al_em,dl_em,num_rgrid_em
  implicit none
  real(8),intent(in) :: rsize(3),hgs(3)
  type(s_rgrid)      :: lg
  !
  real(8),parameter :: epsilon=1.d-10
  integer :: j

  lg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  lg%nd = Nd

  select case(iperiodic)
    case(0)
      lg%ie(:)=int((rsize(:)+epsilon)/2.d0/Hgs(:))
      do j=1,3
        if(mod(int(rsize(j)/Hgs(j)+1.d-12),2)==1)then
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))
        else
          lg%is(j)=-(int((rsize(j)+epsilon)/2.d0/Hgs(j)))+1
        end if
      end do
    case(3)
      lg%is(:)=1
      lg%ie(:)=int((rsize(:)+epsilon)/Hgs(:))
  end select
  lg%num(:)=lg%ie(:)-lg%is(:)+1

  lg%is_overlap(1:3) = lg%is(1:3) - Nd
  lg%ie_overlap(1:3) = lg%ie(1:3) + Nd
  lg%is_array(1:3)   = lg%is_overlap(1:3)
  lg%ie_array(1:3)   = lg%ie_overlap(1:3)

  if ( allocated(lg%idx) ) deallocate(lg%idx)
  if ( allocated(lg%idy) ) deallocate(lg%idy)
  if ( allocated(lg%idz) ) deallocate(lg%idz)
  allocate(lg%idx(lg%is_overlap(1):lg%ie_overlap(1)) &
          ,lg%idy(lg%is_overlap(2):lg%ie_overlap(2)) &
          ,lg%idz(lg%is_overlap(3):lg%ie_overlap(3)))

  do j=lg%is_overlap(1),lg%ie_overlap(1)
    lg%idx(j) = j
  end do
  do j=lg%is_overlap(2),lg%ie_overlap(2)
    lg%idy(j) = j
  end do
  do j=lg%is_overlap(3),lg%ie_overlap(3)
    lg%idz(j) = j
  end do

  select case(theory)
  case('maxwell')
    if(sum(abs(dl_em)) <= 1d-12) then
      if( maxval(abs(num_rgrid_em-lg%num)) > 0) stop "error: num_rgrid_em /= lg%num"
    else
      if( maxval(abs((al_em/dl_em)-dble(lg%num))) > 1d-4 ) stop "error: abs((al_em/dl_em)-dble(lg%num)) is too large"
    end if
  case default
    if(sum(abs(dl)) <= 1d-12) then
      if( maxval(abs(num_rgrid-lg%num)) > 0) stop "error: num_rgrid /= lg%num"
    else
      if( maxval(abs((rsize/dl)-dble(lg%num))) > 1d-4 ) stop "error: abs((rsize/dl)-dble(lg%num)) is too large"
    end if
  end select

  return
end subroutine init_grid_whole

!===================================================================================================================================

subroutine check_ffte_condition(info,lg)
  use structures
  use salmon_global, only: yn_ffte
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_rgrid),        intent(in) :: lg
  integer :: mx,my,mz
  integer :: j,lg_num_tmp,ii

  if (yn_ffte == 'y') then
    mx = mod(lg%num(1), info%nprgrid(2))
    my = mod(lg%num(2), info%nprgrid(2))
    if (mx /= 0 .or. my /= 0) stop 'Both lg%num(1) and lg%num(2) must be divisible by nproc_domain_orbital(2)'

    my = mod(lg%num(2), info%nprgrid(3))
    mz = mod(lg%num(3), info%nprgrid(3))
    if (my /= 0 .or. mz /= 0) stop 'Both lg%num(2) and lg%num(3) must be divisible by nproc_domain_orbital(3)'

    ! this code treats the situation that lg%num(1:3) is less than or equal to 48,828,125
    do j=1,3
      lg_num_tmp=lg%num(j)
      do ii=1,26
        if(mod(lg_num_tmp,2)==0)then
          lg_num_tmp=lg_num_tmp/2
        end if
      end do

      do ii=1,17
        if(mod(lg_num_tmp,3)==0)then
          lg_num_tmp=lg_num_tmp/3
        end if
      end do

      do ii=1,11
        if(mod(lg_num_tmp,5)==0)then
          lg_num_tmp=lg_num_tmp/5
        end if
      end do

      if(lg_num_tmp/=1) stop "When using FFTE, prime factors for number of grids must be combination of 2, 3 or 5."
    end do
  end if
end subroutine check_ffte_condition

!===================================================================================================================================

subroutine init_grid_parallel(info,lg,mg)
  use communication, only: comm_is_root
  use salmon_global, only: yn_periodic, quiet
  use structures, only: s_rgrid,s_parallel_info
  implicit none
  type(s_parallel_info),intent(in)    :: info
  type(s_rgrid),        intent(inout) :: lg
  type(s_rgrid),        intent(inout) :: mg
  !
  integer :: myrank,nproc,nproc_domain_orbital(3),nproc_k,nproc_ob
  integer :: i1,i2,i3,i4,i5,ibox,j,nsize,npo(3)
  
  myrank = info%id_rko
  nproc  = info%isize_rko

  nproc_k              = info%npk
  nproc_ob             = info%nporbital
  nproc_domain_orbital = info%nprgrid

  if ( allocated(mg%is_all) ) deallocate(mg%is_all)
  if ( allocated(mg%ie_all) ) deallocate(mg%ie_all)
  allocate(mg%is_all(3,0:nproc-1),mg%ie_all(3,0:nproc-1))

! +-------------------------------+
! | mg: r-space grid for orbitals |
! +-------------------------------+

  mg%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  mg%nd = Nd

  do i5=0,nproc_k-1
  do i4=0,nproc_ob-1
  do i3=0,nproc_domain_orbital(3)-1
  do i2=0,nproc_domain_orbital(2)-1
  do i1=0,nproc_domain_orbital(1)-1
    ibox = info%imap(i1,i2,i3,i4,i5)
    npo = [i1,i2,i3]
    do j=1,3
      nsize = (lg%num(j) + nproc_domain_orbital(j) - 1) / nproc_domain_orbital(j)
      mg%is_all(j,ibox) = lg%is(j) + nsize * npo(j)
      mg%ie_all(j,ibox) = mg%is_all(j,ibox) + nsize - 1
      if (mg%ie_all(j,ibox) > lg%ie(j)) then
        mg%ie_all(j,ibox) = lg%ie(j)
      end if
    end do
  end do
  end do
  end do
  end do
  end do

  mg%is(:) = mg%is_all(:,myrank)
  mg%ie(:) = mg%ie_all(:,myrank)

  mg%num(:) = mg%ie(:)-mg%is(:)+1

  mg%is_overlap(1:3) = mg%is(1:3)-nd
  mg%ie_overlap(1:3) = mg%ie(1:3)+nd

  if ( allocated(mg%idx) ) deallocate(mg%idx)
  if ( allocated(mg%idy) ) deallocate(mg%idy)
  if ( allocated(mg%idz) ) deallocate(mg%idz)
  allocate(mg%idx(mg%is_overlap(1):mg%ie_overlap(1)) &
          ,mg%idy(mg%is_overlap(2):mg%ie_overlap(2)) &
          ,mg%idz(mg%is_overlap(3):mg%ie_overlap(3)))

  if(yn_periodic=='y' .and. product(nproc_domain_orbital)==1) then
    if((.not. quiet) .and. comm_is_root(myrank)) &
      & write(*,*) "r-space parallelization: off"
    mg%is_array(1:3) = mg%is(1:3)
    mg%ie_array(1:3) = mg%ie(1:3)
    do j=mg%is_overlap(1),mg%ie_overlap(1)
      mg%idx(j) = mod(j+mg%num(1)-1,mg%num(1))+1
    end do
    do j=mg%is_overlap(2),mg%ie_overlap(2)
      mg%idy(j) = mod(j+mg%num(2)-1,mg%num(2))+1
    end do
    do j=mg%is_overlap(3),mg%ie_overlap(3)
      mg%idz(j) = mod(j+mg%num(3)-1,mg%num(3))+1
    end do
  else
    mg%is_array(1:3) = mg%is(1:3)-nd
    mg%ie_array(1:3) = mg%ie(1:3)+nd
    do j=mg%is_overlap(1),mg%ie_overlap(1)
      mg%idx(j) = j
    end do
    do j=mg%is_overlap(2),mg%ie_overlap(2)
      mg%idy(j) = j
    end do
    do j=mg%is_overlap(3),mg%ie_overlap(3)
      mg%idz(j) = j
    end do
  end if

  if(mg%num(1)<nd .or.mg%num(2)<nd .or.mg%num(3)<nd)then
    stop "The system is small. Please use less number of processors."
  end if

#ifdef USE_OPT_ARRAY_PADDING
  lg%ie_array(2)=lg%ie_array(2) + 1
  mg%ie_array(2)=mg%ie_array(2) + 1
#endif

end subroutine init_grid_parallel

!===================================================================================================================================

subroutine init_reciprocal_grid(lg,mg,fg,system,info,poisson)
  use structures
  use math_constants,  only : pi,zi
  use phys_constants, only: cspeed_au
  use salmon_global, only: dt,yn_ffte,aEwald,theory,cutoff_G2_emfield
#ifdef USE_FFTW
  use salmon_global, only: yn_fftw
  use, intrinsic :: iso_c_binding
  use mpi
#endif
  implicit none
#ifdef USE_FFTW
  include 'fftw3-mpi.f03'
#endif
  type(s_rgrid)          ,intent(in)    :: lg
  type(s_rgrid)          ,intent(in)    :: mg
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_dft_system)     ,intent(in)    :: system
  type(s_parallel_info)  ,intent(in)    :: info
  type(s_poisson)        ,intent(inout) :: poisson
  !
  real(8) :: brl(3,3)
  integer :: ix,iy,iz,kx,ky,kz
  integer :: iix,iiy,iiz
  real(8) :: G2,g(3)
  complex(8) :: tmp

  brl(:,:)=system%primitive_b(:,:)

#ifdef USE_FFTW
  if(yn_fftw=='y') then
    allocate(fg%if_Gzero (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%vec_G(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%coef     (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%exp_ewald(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
  else
#endif
    allocate(fg%if_Gzero (lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%vec_G(1:3,lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%coef     (lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%exp_ewald(lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
#ifdef USE_FFTW
  end if
#endif
  fg%if_Gzero = .false.
  fg%vec_G = 0d0
  fg%coef = 0d0
  fg%exp_ewald = 0d0

  if(theory=='single_scale_maxwell_tddft' .and. yn_ffte=='y') then
  ! for single-scale Maxwell-TDDFT
    allocate(fg%coef_nabla(lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),3))
    allocate(fg%coef_gxgy0(lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%cos_cGdt  (lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(fg%sin_cGdt  (lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    fg%coef_nabla = 0d0
    fg%coef_gxgy0 = 1d0
    fg%cos_cGdt   = 0d0
    fg%sin_cGdt   = 0d0
  end if

#ifdef USE_FFTW
  if(yn_fftw=='y') then
    do iz=mg%is(3),mg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
  
      if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) fg%if_Gzero(ix,iy,iz) = .true.
      iix=ix-1-lg%num(1)*(1+sign(1,(ix-1-(lg%num(1)+1)/2)))/2
      iiy=iy-1-lg%num(2)*(1+sign(1,(iy-1-(lg%num(2)+1)/2)))/2
      iiz=iz-1-lg%num(3)*(1+sign(1,(iz-1-(lg%num(3)+1)/2)))/2
      g(1) = dble(iix)*brl(1,1) + dble(iiy)*brl(1,2) + dble(iiz)*brl(1,3)
      g(2) = dble(iix)*brl(2,1) + dble(iiy)*brl(2,2) + dble(iiz)*brl(2,3)
      g(3) = dble(iix)*brl(3,1) + dble(iiy)*brl(3,2) + dble(iiz)*brl(3,3)
      fg%vec_G(1,ix,iy,iz) = g(1)
      fg%vec_G(2,ix,iy,iz) = g(2)
      fg%vec_G(3,ix,iy,iz) = g(3)
      G2 = g(1)**2+g(2)**2+g(3)**2
      if(fg%if_Gzero(ix,iy,iz)) then
        fg%coef(ix,iy,iz) = 0.d0
      else
        fg%coef(ix,iy,iz) = 4.d0*pi/G2
      end if
      fg%exp_ewald(ix,iy,iz) = exp(-G2/(4d0*aEwald))
  
    end do
    end do
    end do
  else
#endif
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=lg%is(1),lg%ie(1)
  
      if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) fg%if_Gzero(ix,iy,iz) = .true.
      iix=ix-1-lg%num(1)*(1+sign(1,(ix-1-(lg%num(1)+1)/2)))/2
      iiy=iy-1-lg%num(2)*(1+sign(1,(iy-1-(lg%num(2)+1)/2)))/2
      iiz=iz-1-lg%num(3)*(1+sign(1,(iz-1-(lg%num(3)+1)/2)))/2
      g(1) = dble(iix)*brl(1,1) + dble(iiy)*brl(1,2) + dble(iiz)*brl(1,3)
      g(2) = dble(iix)*brl(2,1) + dble(iiy)*brl(2,2) + dble(iiz)*brl(2,3)
      g(3) = dble(iix)*brl(3,1) + dble(iiy)*brl(3,2) + dble(iiz)*brl(3,3)
      fg%vec_G(1,ix,iy,iz) = g(1)
      fg%vec_G(2,ix,iy,iz) = g(2)
      fg%vec_G(3,ix,iy,iz) = g(3)
      G2 = g(1)**2+g(2)**2+g(3)**2
      if(fg%if_Gzero(ix,iy,iz)) then
        fg%coef(ix,iy,iz) = 0.d0
      else
        fg%coef(ix,iy,iz) = 4.d0*pi/G2
      end if
      fg%exp_ewald(ix,iy,iz) = exp(-G2/(4d0*aEwald))
      if(theory=='single_scale_maxwell_tddft' .and. yn_ffte=='y') then
      ! for single-scale Maxwell-TDDFT
        fg%coef_nabla(ix,iy,iz,1) = -zi*g(1)
        fg%coef_nabla(ix,iy,iz,2) = -zi*g(2)
        fg%coef_nabla(ix,iy,iz,3) = -zi*g(3)
        if(ix==1.and.iy==1) fg%coef_gxgy0(ix,iy,iz) = 0d0
        if(cutoff_G2_emfield > 0d0 .and. G2 > cutoff_G2_emfield) fg%coef_gxgy0(ix,iy,iz) = 0d0
        fg%cos_cGdt(ix,iy,iz) = cos(cspeed_au*sqrt(G2)*dt)
        fg%sin_cGdt(ix,iy,iz) = sin(cspeed_au*sqrt(G2)*dt)
      end if
  
    enddo
    enddo
    enddo
#ifdef USE_FFTW
  end if
#endif

  if(yn_ffte=='n') then
  ! discrete Fourier transform (general)

    allocate(fg%egx(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
    allocate(fg%egxc(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
    allocate(fg%egy(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
    allocate(fg%egyc(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
    allocate(fg%egz(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))
    allocate(fg%egzc(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))

    allocate(poisson%ff1x(lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(poisson%ff1y(mg%is(1):mg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(poisson%ff1z(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),lg%is(3):lg%ie(3)))
    allocate(poisson%ff2x(lg%is(1):lg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(poisson%ff2y(mg%is(1):mg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(poisson%ff2z(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),lg%is(3):lg%ie(3)))

  !$OMP parallel do private(ix,kx,tmp)
    do ix=lg%is(1),lg%ie(1)
      do kx=lg%is(1),lg%ie(1)
        tmp = exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(lg%num(1))))
        fg%egx(kx,ix)  = tmp
        fg%egxc(kx,ix) = conjg(tmp)
      end do
    end do
  !$OMP parallel do private(iy,ky,tmp)
    do iy=lg%is(2),lg%ie(2)
      do ky=lg%is(2),lg%ie(2)
        tmp = exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(lg%num(2))))
        fg%egy(ky,iy)  = tmp
        fg%egyc(ky,iy) = conjg(tmp)
      end do
    end do
  !$OMP parallel do private(iz,kz,tmp)
    do iz=lg%is(3),lg%ie(3)
      do kz=lg%is(3),lg%ie(3)
        tmp = exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(lg%num(3))))
        fg%egz(kz,iz)  = tmp
        fg%egzc(kz,iz) = conjg(tmp)
      end do
    end do

  else
  ! FFTE

    allocate(poisson%a_ffte(lg%num(1),mg%num(2),mg%num(3)))
    allocate(poisson%b_ffte(lg%num(1),mg%num(2),mg%num(3)))

  ! FFTE initialization step
    call PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,lg%num(1),lg%num(2),lg%num(3), &
                      info%isize_y,info%isize_z,0, &
                      info%icomm_y,info%icomm_z)

  end if

  allocate(poisson%zrhoG_ele(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

#ifdef USE_FFTW
  if(yn_fftw=='y') then
    call fftw_mpi_init()
    allocate(poisson%fftw1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
    allocate(poisson%fftw2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),mg%is(3):mg%ie(3)))
!$OMP parallel do private(ix,iy,iz)
    do iz=mg%is(3),mg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      poisson%fftw1(ix,iy,iz)=0.d0
      poisson%fftw2(ix,iy,iz)=0.d0
    end do
    end do
    end do
  end if
#endif

  return
end subroutine init_reciprocal_grid

!===================================================================================================================================

subroutine init_reciprocal_grid_isolated_ft(lg,mg,fg,system,info,poisson)
  use structures, only: s_rgrid,s_reciprocal_grid,s_dft_system,s_parallel_info,s_poisson
  use math_constants,  only : pi,zi
  use salmon_global, only: yn_ffte,aEwald
#ifdef USE_FFTW
  use salmon_global, only: yn_fftw
  use, intrinsic :: iso_c_binding
  use mpi
#endif
  implicit none
#ifdef USE_FFTW
  include 'fftw3-mpi.f03'
#endif
  type(s_rgrid)          ,intent(in)    :: lg
  type(s_rgrid)          ,intent(in)    :: mg
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_dft_system)     ,intent(in)    :: system
  type(s_parallel_info)  ,intent(in)    :: info
  type(s_poisson)        ,intent(inout) :: poisson
  !
  real(8) :: brl(3,3)
  integer :: ix,iy,iz
  integer :: kx,ky,kz
  integer :: iix,iiy,iiz
  real(8) :: G2,g(3)
  integer :: imgx_s,imgx_e
  integer :: imgy_s,imgy_e
  integer :: imgz_s,imgz_e
  integer :: ifgx_s,ifgx_e
  integer :: ifgy_s,ifgy_e
  integer :: ifgz_s,ifgz_e
  complex(8) :: tmp

#ifdef USE_FFTW
  if(yn_fftw=='n')then
#endif
    if(yn_ffte=='n') then
      ifgx_s = (mg%is(1)-lg%is(1))*2+1
      ifgx_e = (mg%is(1)-lg%is(1))*2+mg%num(1)*2
      ifgy_s = (mg%is(2)-lg%is(2))*2+1
      ifgy_e = (mg%is(2)-lg%is(2))*2+mg%num(2)*2
      ifgz_s = (mg%is(3)-lg%is(3))*2+1
      ifgz_e = (mg%is(3)-lg%is(3))*2+mg%num(3)*2
    else
      if(mod(info%nporbital,4)==0)then
        ! start and end point of reciprocal grids for x, y, z
        ifgx_s = 1
        ifgx_e = 2*lg%num(1)
        if(info%id_y_isolated_ffte >= info%isize_y_isolated_ffte/2) then
          ifgy_s = mg%is(2)-lg%is(2)+1+lg%num(2)
        else
          ifgy_s = mg%is(2)-lg%is(2)+1
        end if
        ifgy_e = ifgy_s+mg%num(2)-1
        if(info%id_z_isolated_ffte >= info%isize_z_isolated_ffte/2) then
          ifgz_s = mg%is(3)-lg%is(3)+1+lg%num(3)
        else
          ifgz_s = mg%is(3)-lg%is(3)+1
        end if
        ifgz_e = ifgz_s+mg%num(3)-1
      else
        ! start and end point of reciprocal grids for x, y, z
        ifgx_s = 1
        ifgx_e = 2*lg%num(1)
        ifgy_s = 1
        ifgy_e = 2*lg%num(2)
        ifgz_s = 1
        ifgz_e = 2*lg%num(3)
      end if
    end if
  
    if(yn_ffte=='n') then
      imgx_s = mg%is(1)-lg%is(1)+1
      imgx_e = mg%is(1)-lg%is(1)+mg%num(1)
      imgy_s = mg%is(2)-lg%is(2)+1
      imgy_e = mg%is(2)-lg%is(2)+mg%num(2)
      imgz_s = mg%is(3)-lg%is(3)+1
      imgz_e = mg%is(3)-lg%is(3)+mg%num(3)
      
      allocate(fg%egx(1:2*lg%num(1),1:2*lg%num(1)))  ! (kx, ix)
      allocate(fg%egxc(1:2*lg%num(1),1:2*lg%num(1))) ! (kx, ix)
      allocate(fg%egy(1:2*lg%num(2),1:2*lg%num(2)))  ! (ky, iy)
      allocate(fg%egyc(1:2*lg%num(2),1:2*lg%num(2))) ! (ky, iy)
      allocate(fg%egz(1:2*lg%num(3),1:2*lg%num(3)))  ! (kz, iz)
      allocate(fg%egzc(1:2*lg%num(3),1:2*lg%num(3))) ! (kz, iz)
  
      allocate(poisson%ff1x(1:2*lg%num(1),ifgy_s:ifgy_e,ifgz_s:ifgz_e)) ! (ix, ky, kz)
      allocate(poisson%ff1y(imgx_s:imgx_e,1:2*lg%num(2),ifgz_s:ifgz_e)) ! (ix, iy, kz)
      allocate(poisson%ff1z(imgx_s:imgx_e,imgy_s:imgy_e,1:2*lg%num(3))) ! (ix, iy, iz)
      allocate(poisson%ff2x(1:2*lg%num(1),ifgy_s:ifgy_e,ifgz_s:ifgz_e)) ! (ix, ky, kz)
      allocate(poisson%ff2y(imgx_s:imgx_e,1:2*lg%num(2),ifgz_s:ifgz_e)) ! (ix, iy, kz)
      allocate(poisson%ff2z(imgx_s:imgx_e,imgy_s:imgy_e,1:2*lg%num(3))) ! (ix, iy, iz)
      
      allocate(poisson%ff1(1:2*lg%num(1),ifgy_s:ifgy_e,ifgz_s:ifgz_e))  ! (kx, ky, kz)
      allocate(poisson%ff2(1:2*lg%num(1),ifgy_s:ifgy_e,ifgz_s:ifgz_e))  ! (kx, ky, kz)
  
      allocate(poisson%ff3x(1:2*lg%num(1),imgy_s:imgy_e,imgz_s:imgz_e)) ! (kx, iy, iz)
      allocate(poisson%ff3y(ifgx_s:ifgx_e,1:2*lg%num(2),imgz_s:imgz_e)) ! (kx, ky, iz)
      allocate(poisson%ff3z(ifgx_s:ifgx_e,ifgy_s:ifgy_e,1:2*lg%num(3))) ! (kx, ky, kz)
      allocate(poisson%ff4x(1:2*lg%num(1),imgy_s:imgy_e,imgz_s:imgz_e)) ! (kx, iy, iz)
      allocate(poisson%ff4y(ifgx_s:ifgx_e,1:2*lg%num(2),imgz_s:imgz_e)) ! (kx, ky, iz)
      allocate(poisson%ff4z(ifgx_s:ifgx_e,ifgy_s:ifgy_e,1:2*lg%num(3))) ! (kx, ky, kz)
  
    !$OMP parallel do private(ix,kx,tmp)
      do ix=1,lg%num(1)
        do kx=1,2*lg%num(1)
          !tmp = exp(zI*(2.d0*Pi*dble(lg%coordinate(ix,1)*system%hgs(1)*(kx-1))/dble(2*lg%num(1))))
          tmp = exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(2*lg%num(1))))
          fg%egx(kx,ix)  = tmp
          fg%egxc(kx,ix) = conjg(tmp)
          !tmp = exp(zI*(2.d0*Pi*dble((lg%coordinate(ix,1)*system%hgs(1)+lg%num(1))*(kx-1))/dble(2*lg%num(1))))
          tmp = exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(2*lg%num(1))))
          fg%egx(kx,ix+lg%num(1))  = tmp
          fg%egxc(kx,ix+lg%num(1)) = conjg(tmp)
        end do
      end do
    !$OMP parallel do private(iy,ky,tmp)
      do iy=1,lg%num(2)
        do ky=1,2*lg%num(2)
          !tmp = exp(zI*(2.d0*Pi*dble(lg%coordinate(iy,2)*system%hgs(2)*(ky-1))/dble(2*lg%num(2))))
          tmp = exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(2*lg%num(2))))
          fg%egy(ky,iy)  = tmp
          fg%egyc(ky,iy) = conjg(tmp)
          !tmp = exp(zI*(2.d0*Pi*dble((lg%coordinate(iy,2)*system%hgs(2)+lg%num(2))*(ky-1))/dble(2*lg%num(2))))
          tmp = exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(2*lg%num(2))))
          fg%egy(ky,iy+lg%num(2))  = tmp
          fg%egyc(ky,iy+lg%num(2)) = conjg(tmp)
        end do
      end do
    !$OMP parallel do private(iz,kz,tmp)
      do iz=1,lg%num(3)
        do kz=1,2*lg%num(3)
          !tmp = exp(zI*(2.d0*Pi*dble(lg%coordinate(iz,3)*system%hgs(3)*(kz-1))/dble(2*lg%num(3))))
          tmp = exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(2*lg%num(3))))
          fg%egz(kz,iz)  = tmp
          fg%egzc(kz,iz) = conjg(tmp)
          !tmp = exp(zI*(2.d0*Pi*dble((lg%coordinate(iz,3)*system%hgs(3)+lg%num(3))*(kz-1))/dble(2*lg%num(3))))
          tmp = exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(2*lg%num(3))))
          fg%egz(kz,iz+lg%num(3))  = tmp
          fg%egzc(kz,iz+lg%num(3)) = conjg(tmp)
        end do
      end do
    else
    ! FFTE
      if(mod(info%nporbital,4)==0)then
        allocate(poisson%a_ffte(2*lg%num(1),mg%num(2),mg%num(3)))
        allocate(poisson%b_ffte(2*lg%num(1),mg%num(2),mg%num(3)))
      else
        allocate(poisson%a_ffte(2*lg%num(1),2*lg%num(2),2*lg%num(3)))
        allocate(poisson%b_ffte(2*lg%num(1),2*lg%num(2),2*lg%num(3)))
      end if
      ! FFTE initialization step
      call PZFFT3DV_MOD(poisson%a_ffte,poisson%b_ffte,2*lg%num(1),2*lg%num(2),2*lg%num(3), &
                        info%isize_y_isolated_ffte,info%isize_z_isolated_ffte,0, &
                        info%icomm_y_isolated_ffte,info%icomm_z_isolated_ffte)
    end if
#ifdef USE_FFTW
  else if(yn_fftw=='y')then

    if(mod(info%nporbital,2)==0)then
      ! start and end point of reciprocal grids for x, y, z
      ifgx_s = 1
      ifgx_e = 2*lg%num(1)
      ifgy_s = 1
      ifgy_e = 2*lg%num(2)
      if(info%iaddress_isolated_fftw(4)==1) then
        ifgz_s = mg%is(3)-lg%is(3)+1+lg%num(3)
      else
        ifgz_s = mg%is(3)-lg%is(3)+1
      end if
      ifgz_e = ifgz_s+mg%num(3)-1
    else
      ! start and end point of reciprocal grids for x, y, z
      ifgx_s = 1
      ifgx_e = 2*lg%num(1)
      ifgy_s = 1
      ifgy_e = 2*lg%num(2)
      ifgz_s = 1
      ifgz_e = 2*lg%num(3)
    end if

    if(mod(info%nporbital,2)==0)then
      call fftw_mpi_init()
      allocate(poisson%fftw1(2*lg%num(1),2*lg%num(2),mg%num(3)))
      allocate(poisson%fftw2(2*lg%num(1),2*lg%num(2),mg%num(3)))
!$OMP parallel do private(ix,iy,iz)
      do iz=1,mg%num(3)
      do iy=1,2*lg%num(2)
      do ix=1,2*lg%num(1)
        poisson%fftw1(ix,iy,iz)=0.d0
        poisson%fftw2(ix,iy,iz)=0.d0
      end do
      end do
      end do
    else
      allocate(poisson%fftw1(2*lg%num(1),2*lg%num(2),2*lg%num(3)))
      allocate(poisson%fftw2(2*lg%num(1),2*lg%num(2),2*lg%num(3)))
!$OMP parallel do private(ix,iy,iz)
      do iz=1,2*lg%num(3)
      do iy=1,2*lg%num(2)
      do ix=1,2*lg%num(1)
        poisson%fftw1(ix,iy,iz)=0.d0
        poisson%fftw2(ix,iy,iz)=0.d0
      end do
      end do
      end do
    end if
  end if
#endif

  allocate(poisson%zrhoG_ele(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e))
  poisson%zrhoG_ele = 0.d0
  
  brl(:,:)=system%primitive_b(:,:)

  allocate(fg%if_Gzero (ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e))
  allocate(fg%vec_G(1:3,ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e))
  allocate(fg%coef     (ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e))
  allocate(fg%exp_ewald(ifgx_s:ifgx_e,ifgy_s:ifgy_e,ifgz_s:ifgz_e))
  fg%if_Gzero = .false.
  fg%vec_G = 0d0
  fg%coef = 0d0
  fg%exp_ewald = 0d0

  do iz=ifgz_s,ifgz_e
  do iy=ifgy_s,ifgy_e
  do ix=ifgx_s,ifgx_e

    iix=ix-1-2*lg%num(1)*(1+sign(1,(ix-1-(2*lg%num(1)+1)/2)))/2
    iiy=iy-1-2*lg%num(2)*(1+sign(1,(iy-1-(2*lg%num(2)+1)/2)))/2
    iiz=iz-1-2*lg%num(3)*(1+sign(1,(iz-1-(2*lg%num(3)+1)/2)))/2
    
    if(iix**2+iiy**2+iiz**2 == 0) fg%if_Gzero(ix,iy,iz) = .true.

    g(1) = dble(iix)*brl(1,1)/2.d0
    g(2) = dble(iiy)*brl(2,2)/2.d0
    g(3) = dble(iiz)*brl(3,3)/2.d0
    fg%vec_G(1,ix,iy,iz) = g(1)
    fg%vec_G(2,ix,iy,iz) = g(2)
    fg%vec_G(3,ix,iy,iz) = g(3)
    G2 = g(1)**2+g(2)**2+g(3)**2
    
    if(fg%if_Gzero(ix,iy,iz)) then
      fg%coef(ix,iy,iz) = 0.d0
    else
      fg%coef(ix,iy,iz) = 4.d0*pi/G2
    end if
    fg%exp_ewald(ix,iy,iz) = exp(-G2/(4d0*aEwald))

  enddo
  enddo
  enddo

  return
end subroutine init_reciprocal_grid_isolated_ft

!===================================================================================================================================

subroutine set_gridcoordinate(lg,system)
  use structures, only: s_rgrid,s_dft_system
  use salmon_global, only: iperiodic
  implicit none
  type(s_rgrid),     intent(inout) :: lg
  type(s_dft_system),intent(in)    :: system
  integer :: ix,iy,iz

  allocate(lg%coordinate(minval(lg%is_overlap(1:3)):maxval(lg%ie_overlap(1:3)),3))

  select case(iperiodic)
  case(0)
    select case(mod(lg%num(1),2))
      case(1)
!$OMP parallel do
        do ix=lg%is_overlap(1),lg%ie_overlap(1)
          lg%coordinate(ix,1)=dble(ix)*system%hgs(1)
        end do
      case(0)
!$OMP parallel do
        do ix=lg%is_overlap(1),lg%ie_overlap(1)
          lg%coordinate(ix,1)=(dble(ix)-0.5d0)*system%hgs(1)
        end do
    end select

    select case(mod(lg%num(2),2))
      case(1)
!$OMP parallel do
        do iy=lg%is_overlap(2),lg%ie_overlap(2)
          lg%coordinate(iy,2)=dble(iy)*system%hgs(2)
        end do
      case(0)
!$OMP parallel do
      do iy=lg%is_overlap(2),lg%ie_overlap(2)
        lg%coordinate(iy,2)=(dble(iy)-0.5d0)*system%hgs(2)
      end do
    end select

    select case(mod(lg%num(3),2))
      case(1)
!$OMP parallel do
        do iz=lg%is_overlap(3),lg%ie_overlap(3)
          lg%coordinate(iz,3)=dble(iz)*system%hgs(3)
        end do
      case(0)
!$OMP parallel do
        do iz=lg%is_overlap(3),lg%ie_overlap(3)
          lg%coordinate(iz,3)=(dble(iz)-0.5d0)*system%hgs(3)
        end do
    end select
  case(3)
!$OMP parallel do
    do ix=lg%is_overlap(1),lg%ie_overlap(1)
      lg%coordinate(ix,1)=dble(ix-1)*system%hgs(1)
    end do
!$OMP parallel do
    do iy=lg%is_overlap(2),lg%ie_overlap(2)
      lg%coordinate(iy,2)=dble(iy-1)*system%hgs(2)
    end do
!$OMP parallel do
    do iz=lg%is_overlap(3),lg%ie_overlap(3)
      lg%coordinate(iz,3)=dble(iz-1)*system%hgs(3)
    end do
  end select

end subroutine set_gridcoordinate

!===================================================================================================================================

subroutine init_nion_div(system,lg,mg,info)
  use structures, only: s_dft_system, s_reciprocal_grid, s_rgrid, s_parallel_info
  use communication, only: comm_summation, comm_get_groupinfo, comm_is_root
 !use parallelization, only: nproc_id_global
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_rgrid)     ,intent(in) :: lg
  type(s_rgrid)     ,intent(in) :: mg
  type(s_parallel_info)         :: info
  !
  logical :: flag_cuboid
 !integer :: k,irank,nproc
  integer :: ia,j,nc,ix,iy,iz,iia,nion_total
  real(8) :: hgs(3), Rion_tmp(3,system%nion), al0(3,3)
  real(8) :: r_mg_min(3), r_mg_max(3), al_min(3), al_max(3), al_len(3)

  al0(:,:) = system%primitive_a(:,:)
  if( abs(al0(1,2)).ge.1d-10 .or. &
      abs(al0(1,3)).ge.1d-10 .or. &
      abs(al0(2,3)).ge.1d-10 )  then
     flag_cuboid=.false.
  else
     flag_cuboid = .true.
  endif


  if(flag_cuboid) then

  nc = 2
  hgs(:)    = system%hgs(:)

  r_mg_min(:) = (mg%is(:)-1) * hgs(:)
  r_mg_max(:) =  mg%ie(:)    * hgs(:)

  al_min(:) = 0d0
  al_max(:) = system%hgs(:)*dble(lg%num(:))
  al_len(:) = al_max(:) - al_min(:)

  info%nion_mg = 0

  do ia = 1,system%nion
     j=1
     do ix = -nc, nc
        Rion_tmp(j,ia) = system%Rion(j,ia) + ix*al_len(j)
        if( Rion_tmp(j,ia) .ge. al_min(j) .and. Rion_tmp(j,ia) .lt. al_max(j) ) exit
     enddo
     j=2
     do iy = -nc, nc
        Rion_tmp(j,ia) = system%Rion(j,ia) + iy*al_len(j)
        if( Rion_tmp(j,ia) .ge. al_min(j) .and. Rion_tmp(j,ia) .lt. al_max(j) ) exit
     enddo
     j=3
     do iz = -nc, nc
        Rion_tmp(j,ia) = system%Rion(j,ia) + iz*al_len(j)
        if( Rion_tmp(j,ia) .ge. al_min(j) .and. Rion_tmp(j,ia) .lt. al_max(j) ) exit
     enddo

     if( (Rion_tmp(1,ia).ge.r_mg_min(1) .and. Rion_tmp(1,ia).lt.r_mg_max(1)) .and. &
         (Rion_tmp(2,ia).ge.r_mg_min(2) .and. Rion_tmp(2,ia).lt.r_mg_max(2)) .and. &
         (Rion_tmp(3,ia).ge.r_mg_min(3) .and. Rion_tmp(3,ia).lt.r_mg_max(3)) ) then
        info%nion_mg = info%nion_mg + 1
     endif
  enddo

  allocate( info%ia_mg(info%nion_mg) )

  iia = 0
  do ia = 1,system%nion
     if( (Rion_tmp(1,ia).ge.r_mg_min(1) .and. Rion_tmp(1,ia).lt.r_mg_max(1)) .and. &
         (Rion_tmp(2,ia).ge.r_mg_min(2) .and. Rion_tmp(2,ia).lt.r_mg_max(2)) .and. &
         (Rion_tmp(3,ia).ge.r_mg_min(3) .and. Rion_tmp(3,ia).lt.r_mg_max(3)) ) then
        iia = iia + 1
        info%ia_mg(iia) = ia
     endif
  enddo

  if( info%nion_mg .ne. iia ) stop "Error1 in dividing atom in mg domain"

  else !(flag_cuboid=.false.)

     ! assuming r-space parallelization is not available in nonorthogonal lattice cell
     info%nion_mg = system%nion
     allocate( info%ia_mg(info%nion_mg) )
     do ia=1,system%nion
        info%ia_mg(ia) = ia
     enddo

  endif

  !check
  call comm_summation(info%nion_mg, nion_total, info%icomm_r)
  if( nion_total .ne. system%nion ) stop "Error2 in dividing atom in mg domain"

  !write(*,*) "  #nion_mg=", system%nion_mg
  !write(*,*) "  #check nion_total=", nion_total


  !!(divide nion with all processes: not used now)
  !call comm_get_groupinfo(fg%icomm_G,irank,nproc)
  !
  !if(nproc .le. system%nion) then
  !   k = mod(system%nion,nproc)
  !   if(k==0) then
  !      system%nion_r = system%nion / nproc
  !   else
  !      system%nion_r = system%nion / nproc + 1
  !   endif
  !   system%nion_s = system%nion_r * irank + 1
  !   system%nion_e = system%nion_s + system%nion_r - 1
  !   if (irank == nproc-1) system%nion_e = system%nion
  !   if (system%nion_e .gt. system%nion) system%nion_e = -1
  !   if (system%nion_s .gt. system%nion) then
  !      system%nion_s =  0
  !      system%nion_e = -1
  !   endif
  !
  !else
  !   if(irank+1.le.system%nion) then
  !      system%nion_s = irank + 1
  !      system%nion_e = system%nion_s
  !   else
  !      system%nion_s = 0
  !      system%nion_e = -1
  !   endif
  !endif

end subroutine init_nion_div

!===================================================================================================================================

subroutine set_bN(bnmat)
  implicit none
  real(8) :: bnmat(4,4)

  bNmat(1,1)=1.d0/2.d0

  bNmat(1,2)=2.d0/3.d0
  bNmat(2,2)=-1.d0/12.d0

  bNmat(1,3)=3.d0/4.d0
  bNmat(2,3)=-3.d0/20.d0
  bNmat(3,3)=1.d0/60.d0

  bNmat(1,4)=4.d0/5.d0
  bNmat(2,4)=-1.d0/5.d0
  bNmat(3,4)=4.d0/105.d0
  bNmat(4,4)=-1.d0/280.d0

end subroutine set_bN

subroutine set_cN(cnmat)
  implicit none
  real(8) :: cnmat(0:12,12)

  cNmat(0,1)=-2.d0
  cNmat(1,1)=1.d0

  cNmat(0,2)=-5.d0/2.d0
  cNmat(1,2)=4.d0/3.d0
  cNmat(2,2)=-1.d0/12.d0

  cNmat(0,3)=-49.d0/18.d0
  cNmat(1,3)=3.d0/2.d0
  cNmat(2,3)=-3.d0/20.d0
  cNmat(3,3)=1.d0/90.d0

  cNmat(0,5)=-5269.d0/1800.d0
  cNmat(1,5)=5.d0/3.d0
  cNmat(2,5)=-5.d0/21.d0
  cNmat(3,5)=5.d0/126.d0
  cNmat(4,5)=-5.d0/1008.d0
  cNmat(5,5)=1.d0/3150.d0

  cNmat(0,4)=-205.d0/72.d0
  cNmat(1,4)=8.d0/5.d0
  cNmat(2,4)=-1.d0/5.d0
  cNmat(3,4)=8.d0/315.d0
  cNmat(4,4)=-1.d0/560.d0

  cNmat(0,6)=-5369.d0/1800.d0
  cNmat(1,6)=12.d0/7.d0
  cNmat(2,6)=-15.d0/56.d0
  cNmat(3,6)=10.d0/189.d0
  cNmat(4,6)=-1.d0/112.d0
  cNmat(5,6)=2.d0/1925.d0
  cNmat(6,6)=-1.d0/16632.d0

  cNmat(0,7)=-266681.d0/88200.d0
  cNmat(1,7)=7.d0/4.d0
  cNmat(2,7)=-7.d0/24.d0
  cNmat(3,7)=7.d0/108.d0
  cNmat(4,7)=-7.d0/528.d0
  cNmat(5,7)=7.d0/3300.d0
  cNmat(6,7)=-7.d0/30888.d0
  cNmat(7,7)=1.d0/84084.d0

  cNmat(0,8)=-1077749.d0/352800.d0
  cNmat(1,8)=16.d0/9.d0
  cNmat(2,8)=-14.d0/45.d0
  cNmat(3,8)=112.d0/1485.d0
  cNmat(4,8)=-7.d0/396.d0
  cNmat(5,8)=112.d0/32175.d0
  cNmat(6,8)=-2.d0/3861.d0
  cNmat(7,8)=16.d0/315315.d0
  cNmat(8,8)=-1.d0/411840.d0

  cNmat(0,9)=-9778141.d0/3175200.d0
  cNmat(1,9)=9.d0/5.d0
  cNmat(2,9)=-18.d0/55.d0
  cNmat(3,9)=14.d0/165.d0
  cNmat(4,9)=-63.d0/2860.d0
  cNmat(5,9)=18.d0/3575.d0
  cNmat(6,9)=-2.d0/2145.d0
  cNmat(7,9)=9.d0/70070.d0
  cNmat(8,9)=-9.d0/777920.d0
  cNmat(9,9)=1.d0/1969110.d0

  cNmat(0,10)=-1968329.d0/635040.d0
  cNmat(1,10)=20.d0/11.d0
  cNmat(2,10)=-15.d0/44.d0
  cNmat(3,10)=40.d0/429.d0
  cNmat(4,10)=-15.d0/572.d0
  cNmat(5,10)=24.d0/3575.d0
  cNmat(6,10)=-5.d0/3432.d0
  cNmat(7,10)=30.d0/119119.d0
  cNmat(8,10)=-5.d0/155584.d0
  cNmat(9,10)=10.d0/3741309.d0
  cNmat(10,10)=-1.d0/9237800.d0

  cNmat(0,11)=-239437889.d0/76839840.d0
  cNmat(1,11)=11.d0/6.d0
  cNmat(2,11)=-55.d0/156.d0
  cNmat(3,11)=55.d0/546.d0
  cNmat(4,11)=-11.d0/364.d0
  cNmat(5,11)=11.d0/1300.d0
  cNmat(6,11)=-11.d0/5304.d0
  cNmat(7,11)=55.d0/129948.d0
  cNmat(8,11)=-55.d0/806208.d0
  cNmat(9,11)=11.d0/1360476.d0
  cNmat(10,11)=-11.d0/17635800.d0
  cNmat(11,11)=1.d0/42678636.d0

  cNmat(0,12)=-240505109.d0/76839840.d0
  cNmat(1,12)=24.d0/13.d0
  cNmat(2,12)=-33.d0/91.d0
  cNmat(3,12)=88.d0/819.d0
  cNmat(4,12)=-99.d0/2912.d0
  cNmat(5,12)=396.d0/38675.d0
  cNmat(6,12)=-11.d0/3978.d0
  cNmat(7,12)=132.d0/205751.d0
  cNmat(8,12)=-33.d0/268736.d0
  cNmat(9,12)=44.d0/2380833.d0
  cNmat(10,12)=-3.d0/1469650.d0
  cNmat(11,12)=12.d0/81800719.d0
  cNmat(12,12)=-1.d0/194699232.d0

end subroutine set_cN

end module initialization_sub
