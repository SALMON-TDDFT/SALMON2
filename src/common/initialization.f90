!
!  Copyright 2019 SALMON developers
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

subroutine init_dft(comm,pinfo,info,info_field,lg,mg,ng,system,stencil,fg,poisson,srg,srg_ng,ofile)
  use structures
  use salmon_global, only: iperiodic,layout_multipole, &
                           nproc_k,nproc_ob,nproc_domain_orbital,nproc_domain_general
  use sendrecv_grid
  use init_communicator
  use init_poisson_sub
  use checkpoint_restart_sub, only: init_dir_out_restart
  use sym_rho_sub, only: init_sym_rho
  implicit none
  integer      ,intent(in) :: comm
  type(s_process_info)     :: pinfo
  type(s_orbital_parallel) :: info
  type(s_field_parallel)   :: info_field
  type(s_rgrid)            :: lg,mg,ng
  type(s_dft_system)       :: system
  type(s_stencil)          :: stencil
  type(s_reciprocal_grid)  :: fg
  type(s_poisson)          :: poisson
  type(s_sendrecv_grid)    :: srg,srg_ng
  type(s_ofile)            :: ofile
  !
  integer,dimension(2,3) :: neig,neig_ng

! electron system
  call init_dft_system(lg,system,stencil)

! process distribution
  pinfo%npk              = nproc_k
  pinfo%nporbital        = nproc_ob
  pinfo%npdomain_orbital = nproc_domain_orbital
  pinfo%npdomain_general = nproc_domain_general
  call init_process_distribution(system,comm,pinfo)
  call init_communicator_dft(comm,pinfo,info,info_field)

! parallelization
  call check_ffte_condition(pinfo,lg)
  call init_grid_parallel(info%id_rko,info%isize_rko,pinfo,info,info_field,lg,mg,ng) ! lg --> mg & ng
  call init_orbital_parallel_singlecell(system,info,pinfo)
  ! sendrecv_grid object for wavefunction updates
  call create_sendrecv_neig_mg(neig, info, pinfo, iperiodic) ! neighboring node array
  call init_sendrecv_grid(srg, mg, info%numo*info%numk*system%nspin, info%icomm_rko, neig)
  ! sendrecv_grid object for scalar potential updates
  call create_sendrecv_neig_ng(neig_ng, pinfo, info_field, iperiodic) ! neighboring node array
  call init_sendrecv_grid(srg_ng, ng, 1, info%icomm_rko, neig_ng)

! symmetry

  call init_sym_rho( lg%num, mg%is, mg%ie, info%icomm_r )

! for Poisson equation

  poisson%iterVh = 0 ! Iteration counter
  select case(iperiodic)
  case(0)
    if(layout_multipole==2.or.layout_multipole==3) call make_corr_pole(lg,ng,poisson)
  case(3)
    call init_reciprocal_grid(lg,ng,fg,system,info_field,poisson)
  end select
  call set_ig_bound(lg,ng,poisson)

! output files

  call init_dir_out_restart(ofile)

end subroutine init_dft

!===================================================================================================================================

subroutine init_dft_system(lg,system,stencil)
  use structures
  use lattice
  use salmon_global, only: al_vec1,al_vec2,al_vec3,al,ispin,natom,nelem,nstate,iperiodic,num_kgrid,num_rgrid,dl, &
  & nproc_domain_orbital,rion,rion_red,nelec,calc_mode,temperature,nelec_spin, &
  & iflag_atom_coor,ntype_atom_coor_reduced,epdir_re1,nstate_spin
  use sym_sub, only: init_sym_sub
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
  call init_sym_sub( system%primitive_a, system%primitive_b, epdir_re1 )
  call init_kvector(num_kgrid,system)

  system%iperiodic = iperiodic
  system%nion = natom

  if(ispin==0)then
    system%nspin=1
  else
    system%nspin=2
  end if

  if(calc_mode=='RT'.and.temperature<-1.d-12)then
    if(system%nspin==2.and.sum(nelec_spin(:))>0)then
      system%no = maxval(nelec_spin(:))
    else
      if(mod(nelec,2)==0)then
        system%no = nelec/2
      else
        system%no = (nelec+1)/2
      end if
    end if
  else
    system%no = max( nstate, maxval(nstate_spin) )
  end if

  allocate(system%mass(1:nelem))
  if ( allocated(system%Rion) ) deallocate(system%Rion)
  if ( allocated(system%rocc) ) deallocate(system%rocc)
  if ( allocated(system%Velocity) ) deallocate(system%Velocity)
  if ( allocated(system%Force) ) deallocate(system%Force)
  allocate(system%Rion(3,system%nion),system%rocc(system%no,system%nk,system%nspin))
  allocate(system%Velocity(3,system%nion),system%Force(3,system%nion))

  if(iflag_atom_coor==ntype_atom_coor_reduced) then
    Rion = matmul(system%primitive_a,Rion_red) ! [ a1, a2, a3 ] * R_ion
  end if
  system%Rion = Rion

! initial value of occupation
  system%rocc = 0d0
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

  call set_bn(bnmat)
  call set_cn(cnmat)
  if(stencil%if_orthogonal) then
    stencil%coef_lap0 = -0.5d0*cNmat(0,Nd)*(1.d0/Hgs(1)**2+1.d0/Hgs(2)**2+1.d0/Hgs(3)**2)
  else
    if(nproc_domain_orbital(1)*nproc_domain_orbital(2)*nproc_domain_orbital(3)/=1) &
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

subroutine init_process_distribution(system,icomm1,pinfo)
  use structures, only: s_process_info,s_dft_system
  use parallelization, only: nproc_id_global, nproc_group_global
  use salmon_global, only: theory,ispin
  use communication, only: comm_is_root,comm_bcast
  use set_numcpu
  implicit none
  type(s_dft_system),intent(in)    :: system
  integer, intent(in)              :: icomm1 ! Communicator for single DFT system.
  type(s_process_info),intent(out) :: pinfo
  logical :: if_stop

  if((pinfo%nporbital + sum(pinfo%npdomain_orbital) + sum(pinfo%npdomain_general)) == 0) then
    ! Process distribution is automatically decided by SALMON.
    if (system%ngrid > 16**3) then
      call set_numcpu_general(iprefer_domain_distribution,system%nk,system%no,icomm1,pinfo)
    else
      select case(theory)
      case('dft','dft_band','dft_md','dft2tddft')
        call set_numcpu_general(iprefer_k_distribution,system%nk,system%no,icomm1,pinfo)
      case('tddft_response','tddft_pulse','single_scale_maxwell_tddft','multiscale_experiment')
        call set_numcpu_general(iprefer_orbital_distribution,system%nk,system%no,icomm1,pinfo)
      case default
        stop 'invalid theory'
      end select
    end if
  else
    ! Process distribution is explicitly specified by user.
  end if

  if (comm_is_root(nproc_id_global)) then
    if_stop = .not. check_numcpu(icomm1, pinfo)
  end if
  call comm_bcast(if_stop, nproc_group_global)
  if (if_stop) stop 'fail: check_numcpu'

  if (ispin==1) then
    pinfo%nporbital_spin(1)=(pinfo%nporbital+1)/2
    pinfo%nporbital_spin(2)= pinfo%nporbital   /2
  else
    pinfo%nporbital_spin = 0
  end if
  pinfo%npdomain_general_dm(1:3)=pinfo%npdomain_general(1:3)/pinfo%npdomain_orbital(1:3)
end subroutine init_process_distribution

!===================================================================================================================================

subroutine init_orbital_parallel_singlecell(system,info,pinfo)
  use structures
  implicit none
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel)      :: info
  type(s_process_info)          :: pinfo
  !
  integer :: io,nproc_k,nproc_ob,nproc_domain_orbital(3),m

  nproc_k              = pinfo%npk
  nproc_ob             = pinfo%nporbital
  nproc_domain_orbital = pinfo%npdomain_orbital

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

end subroutine init_orbital_parallel_singlecell

!===================================================================================================================================

subroutine init_grid_whole(rsize,hgs,lg)
  use structures, only: s_rgrid
  use salmon_global, only: iperiodic,dl,num_rgrid,theory,al_em,dl_em
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
    if( maxval(abs((al_em/dl_em)-dble(lg%num))) > 1d-4 ) stop "error: abs((al_em/dl_em)-dble(lg%num)) is too large"
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

subroutine check_ffte_condition(pinfo,lg)
  use structures
  use salmon_global, only: yn_ffte
  implicit none
  type(s_process_info),intent(in) :: pinfo
  type(s_rgrid),       intent(in) :: lg
  integer :: mx,my,mz
  integer :: j,lg_num_tmp,ii

  if (yn_ffte == 'y') then
    mx = mod(lg%num(1), pinfo%npdomain_general(2))
    my = mod(lg%num(2), pinfo%npdomain_general(2))
    if (mx /= 0 .or. my /= 0) stop 'Both lg%num(1) and lg%num(2) must be divisible by nproc_domain_general(2)'

    my = mod(lg%num(2), pinfo%npdomain_general(3))
    mz = mod(lg%num(3), pinfo%npdomain_general(3))
    if (my /= 0 .or. mz /= 0) stop 'Both lg%num(2) and lg%num(3) must be divisible by nproc_domain_general(3)'

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

subroutine init_grid_parallel(myrank,nproc,pinfo,info,info_field,lg,mg,ng)
  use communication, only: comm_is_root
  use salmon_global, only: yn_periodic
  use structures, only: s_process_info,s_rgrid,s_orbital_parallel,s_field_parallel
  implicit none
  integer,             intent(in)    :: myrank,nproc
  type(s_process_info),intent(in)    :: pinfo
  type(s_orbital_parallel),intent(in):: info
  type(s_field_parallel),intent(in)  :: info_field
  type(s_rgrid),       intent(inout) :: lg
  type(s_rgrid),       intent(inout) :: mg,ng
  !
  integer :: nproc_domain_orbital(3),nproc_domain_general(3),nproc_k,nproc_ob
  integer :: i1,i2,i3,i4,i5,ibox,j,nsize,npo(3)

  nproc_k              = pinfo%npk
  nproc_ob             = pinfo%nporbital
  nproc_domain_orbital = pinfo%npdomain_orbital
  nproc_domain_general = pinfo%npdomain_general

  if ( allocated(mg%is_all) ) deallocate(mg%is_all)
  if ( allocated(mg%ie_all) ) deallocate(mg%ie_all)
  if ( allocated(ng%is_all) ) deallocate(ng%is_all)
  if ( allocated(ng%ie_all) ) deallocate(ng%ie_all)
  allocate(mg%is_all(3,0:nproc-1),mg%ie_all(3,0:nproc-1))
  allocate(ng%is_all(3,0:nproc-1),ng%ie_all(3,0:nproc-1))

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
      if (mg%ie_all(j,ibox) > lg%num(j)) then
        mg%ie_all(j,ibox) = lg%num(j)
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
    if(comm_is_root(myrank)) write(*,*) "r-space parallelization: off"
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

! +-----------------------------+
! | ng: r-space grid for fields |
! +-----------------------------+

  ng%ndir = 3 ! high symmetry nonorthogonal lattice is not implemented
  ng%nd = Nd

  ng%is_all = mg%is_all
  ng%ie_all = mg%ie_all

  ng%is(1:3)  = ng%is_all(:,myrank)
  ng%ie(1:3)  = ng%ie_all(:,myrank)
  ng%num(1:3) = ng%ie(1:3)-ng%is(1:3)+1

  ng%is_overlap(1:3) = ng%is(1:3)-nd
  ng%ie_overlap(1:3) = ng%ie(1:3)+nd
  ng%is_array(1:3) = ng%is(1:3)-nd
  ng%ie_array(1:3) = ng%ie(1:3)+nd

  if ( allocated(ng%idx) ) deallocate(ng%idx)
  if ( allocated(ng%idy) ) deallocate(ng%idy)
  if ( allocated(ng%idz) ) deallocate(ng%idz)
  allocate(ng%idx(ng%is_overlap(1):ng%ie_overlap(1)) &
          ,ng%idy(ng%is_overlap(2):ng%ie_overlap(2)) &
          ,ng%idz(ng%is_overlap(3):ng%ie_overlap(3)))
  do j=ng%is_overlap(1),ng%ie_overlap(1)
    ng%idx(j) = j
  end do
  do j=ng%is_overlap(2),ng%ie_overlap(2)
    ng%idy(j) = j
  end do
  do j=ng%is_overlap(3),ng%ie_overlap(3)
    ng%idz(j) = j
  end do

  if(yn_periodic=='n'.and.(ng%num(1)<nd.or.ng%num(2)<nd.or.ng%num(3)<nd))then
    stop "The system is small. Please use less number of processors."
  end if

#ifdef USE_OPT_ARRAY_PADDING
  lg%ie_array(2)=lg%ie_array(2) + 1
  mg%ie_array(2)=mg%ie_array(2) + 1
  ng%ie_array(2)=ng%ie_array(2) + 1
#endif

! final check...
  do j=1,3
    if (mg%is(j) /= ng%is(j) .or. mg%ie(j) /= ng%ie(j)) then
      stop 'invalid range: mg and ng, please check init_grid_parallel subroutine'
    end if
  end do

end subroutine init_grid_parallel

!===================================================================================================================================

subroutine init_reciprocal_grid(lg,ng,fg,system,info_field,poisson)
  use inputoutput,     only : nelem,yn_ffte
  use math_constants,  only : pi,zi
  use structures,      only : s_rgrid,s_reciprocal_grid,s_dft_system,s_field_parallel,s_poisson
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: ng
  type(s_reciprocal_grid),intent(inout) :: fg
  type(s_dft_system),intent(in) :: system
  type(s_field_parallel),intent(in) :: info_field
  type(s_poisson),intent(inout) :: poisson
  real(8) :: brl(3,3)
  integer :: n,myrank,nproc
  integer :: ix,iy,iz
  integer :: iix,iiy,iiz
  real(8) :: G2
  integer :: kx,ky,kz
  integer :: ky2,kz2
  complex(8) :: tmp

  brl(:,:)=system%primitive_b(:,:)

  fg%icomm_G = info_field%icomm_all
  fg%ng = system%ngrid

  nproc = info_field%isize_all
  myrank = info_field%id_all

  fg%ig_s = myrank*(fg%ng/nproc)+1
  fg%ig_e = (myrank+1)*(fg%ng/nproc)
  if(myrank==nproc-1) fg%ig_e = fg%ng

  allocate(fg%Gx(fg%ng),fg%Gy(fg%ng),fg%Gz(fg%ng))

  fg%Gx = 0.d0
  fg%Gy = 0.d0
  fg%Gz = 0.d0

  n=0
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
    n=n+1
    if((ix-1)**2+(iy-1)**2+(iz-1)**2 == 0) fg%iGzero=n
    iix=ix-1-lg%num(1)*(1+sign(1,(ix-1-(lg%num(1)+1)/2)))/2
    iiy=iy-1-lg%num(2)*(1+sign(1,(iy-1-(lg%num(2)+1)/2)))/2
    iiz=iz-1-lg%num(3)*(1+sign(1,(iz-1-(lg%num(3)+1)/2)))/2
    fg%Gx(n) = dble(iix)*brl(1,1) + dble(iiy)*brl(1,2) + dble(iiz)*brl(1,3)
    fg%Gy(n) = dble(iix)*brl(2,1) + dble(iiy)*brl(2,2) + dble(iiz)*brl(2,3)
    fg%Gz(n) = dble(iix)*brl(3,1) + dble(iiy)*brl(3,2) + dble(iiz)*brl(3,3)
  enddo
  enddo
  enddo
  allocate(poisson%ff1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(poisson%ff1x(lg%is(1):lg%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
  allocate(poisson%ff1y(ng%is(1):ng%ie(1),lg%is(2):lg%ie(2),ng%is(3):ng%ie(3)))
  allocate(poisson%ff1z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
  allocate(poisson%ff2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(poisson%ff2x(lg%is(1):lg%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3)))
  allocate(poisson%ff2y(ng%is(1):ng%ie(1),lg%is(2):lg%ie(2),ng%is(3):ng%ie(3)))
  allocate(poisson%ff2z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
  allocate(poisson%egx(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
  allocate(poisson%egxc(lg%is(1):lg%ie(1),lg%is(1):lg%ie(1)))
  allocate(poisson%egy(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
  allocate(poisson%egyc(lg%is(2):lg%ie(2),lg%is(2):lg%ie(2)))
  allocate(poisson%egz(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))
  allocate(poisson%egzc(lg%is(3):lg%ie(3),lg%is(3):lg%ie(3)))
  allocate(poisson%trho2z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
  allocate(poisson%trho3z(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),lg%is(3):lg%ie(3)))
!$OMP parallel do private(ix,kx,tmp)
  do ix=lg%is(1),lg%ie(1)
    do kx=lg%is(1),lg%ie(1)
      tmp = exp(zI*(2.d0*Pi*dble((ix-1)*(kx-1))/dble(lg%num(1))))
      poisson%egx(kx,ix)  = tmp
      poisson%egxc(kx,ix) = conjg(tmp)
    end do
  end do
!$OMP parallel do private(iy,ky,tmp)
  do iy=lg%is(2),lg%ie(2)
    do ky=lg%is(2),lg%ie(2)
      tmp = exp(zI*(2.d0*Pi*dble((iy-1)*(ky-1))/dble(lg%num(2))))
      poisson%egy(ky,iy)  = tmp
      poisson%egyc(ky,iy) = conjg(tmp)
    end do
  end do
!$OMP parallel do private(iz,kz,tmp)
  do iz=lg%is(3),lg%ie(3)
    do kz=lg%is(3),lg%ie(3)
      tmp = exp(zI*(2.d0*Pi*dble((iz-1)*(kz-1))/dble(lg%num(3))))
      poisson%egz(kz,iz)  = tmp
      poisson%egzc(kz,iz) = conjg(tmp)
    end do
  end do

  select case(yn_ffte)
  case('n')
    ! no operation
  case('y')
    if(.not.allocated(poisson%coef))then
      allocate(poisson%coef      (lg%num(1),ng%num(2),ng%num(3)))
      allocate(poisson%a_ffte    (lg%num(1),ng%num(2),ng%num(3)))
      allocate(poisson%a_ffte_tmp(lg%num(1),ng%num(2),ng%num(3)))
      allocate(poisson%b_ffte    (lg%num(1),ng%num(2),ng%num(3)))
    end if
    poisson%coef=0.d0

    do kz = 1,ng%num(3)
    do ky = 1,ng%num(2)
    do kx = lg%is(1),lg%ie(1)
      kz2=kz+ng%is(3)-1
      ky2=ky+ng%is(2)-1
      n = (kz2-lg%is(3))*lg%num(1)*lg%num(2) &
        + (ky2-lg%is(2))*lg%num(1) &
        + (kx -lg%is(1)) &
        + 1
      G2=fg%Gx(n)**2+fg%Gy(n)**2+fg%Gz(n)**2
      if(kx==1.and.ky2==1.and.kz2==1)then
        poisson%coef(kx,ky,kz)=0.d0
      else
        poisson%coef(kx,ky,kz)=4.d0*pi/G2
      end if
    end do
    end do
    end do
  end select

  allocate(fg%zrhoG_ion(fg%ng),fg%zrhoG_ele(fg%ng),fg%zrhoG_ele_tmp(fg%ng),fg%zVG_ion(fg%ng,nelem))

  return
end subroutine init_reciprocal_grid

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
