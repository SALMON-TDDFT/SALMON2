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
!-----------------------------------------------------------------------------------------
module fdtd_eh
  implicit none
  
  type ls_fdtd_eh
    integer             :: Nd              !number of additional grid in mpi
    integer             :: iter_sta        !start of time-iteration
    integer             :: iter_end        !end of time-iteration
    integer             :: iter_now        !present iteration Number
    integer             :: ifn             !file number for inputing or saving data
    integer             :: ioddeven(3)     !odd or even grid paterns
    integer             :: ipml_l          !pml parameter
    real(8)             :: pml_m           !pml parameter
    real(8)             :: pml_r           !pml parameter
    real(8),allocatable :: c_ebloch_x(:)   !coeff. of bloch boundary condition used for de/dx
    real(8),allocatable :: c_ebloch_y(:)   !coeff. of bloch boundary condition used for de/dy 
    real(8),allocatable :: c_ebloch_z(:)   !coeff. of bloch boundary condition used for de/dz
    real(8),allocatable :: c_hbloch_x(:)   !coeff. of bloch boundary condition used for dh/dx
    real(8),allocatable :: c_hbloch_y(:)   !coeff. of bloch boundary condition used for dh/dy
    real(8),allocatable :: c_hbloch_z(:)   !coeff. of bloch boundary condition used for dh/dz
    real(8)             :: e_max           !maximum e for observation
    real(8)             :: h_max           !maximum e for observation
    real(8)             :: uVperm_from_au  !convert parameter for V/m
    real(8)             :: uAperm_from_au  !convert parameter for A/m
    real(8),allocatable :: rep(:)          !relative permittivity(non-dispersion)
    real(8),allocatable :: rmu(:)          !relative permeability(non-dispersion)
    real(8),allocatable :: sig(:)          !conductivity(non-dispersion)
    real(8),allocatable :: coo(:,:)        !grid coordinate
    integer,allocatable :: iobs_po_pe(:)   !processor element at observation point
    integer,allocatable :: iobs_li_pe(:,:) !processor element at observation line
    integer,allocatable :: iobs_pl_pe(:,:) !processor element at observation plane
    integer,allocatable :: iobs_po_id(:,:) !id at observation point
    integer,allocatable :: iobs_num_ene(:) !number of energy point used for obs_plane_ene_em option
    !array used for obs_plane_ene_em option
    !(r1,r2,obs_num_em,iobs_num_ene,with or without window function)
    complex(8),allocatable :: obs_ex_xy_ene(:,:,:,:,:),obs_ex_yz_ene(:,:,:,:,:),obs_ex_xz_ene(:,:,:,:,:), &
                              obs_ey_xy_ene(:,:,:,:,:),obs_ey_yz_ene(:,:,:,:,:),obs_ey_xz_ene(:,:,:,:,:), &
                              obs_ez_xy_ene(:,:,:,:,:),obs_ez_yz_ene(:,:,:,:,:),obs_ez_xz_ene(:,:,:,:,:), &
                              obs_hx_xy_ene(:,:,:,:,:),obs_hx_yz_ene(:,:,:,:,:),obs_hx_xz_ene(:,:,:,:,:), &
                              obs_hy_xy_ene(:,:,:,:,:),obs_hy_yz_ene(:,:,:,:,:),obs_hy_xz_ene(:,:,:,:,:), &
                              obs_hz_xy_ene(:,:,:,:,:),obs_hz_yz_ene(:,:,:,:,:),obs_hz_xz_ene(:,:,:,:,:)
    integer             :: inc_num         !number of incident current source
    integer,allocatable :: inc_po_pe(:)    !processor element at incident current source point
    integer,allocatable :: inc_li_pe(:,:)  !processor element at incident current source line
    integer,allocatable :: inc_pl_pe(:,:)  !processor element at incident current source plane
    integer,allocatable :: inc_po_id(:,:)  !id at incident current source point
    character(16)       :: inc_dist1       !spatial distribution type of inc.
    character(16)       :: inc_dist2       !spatial distribution type of inc.
    real(8)             :: c2_inc1_xyz(3)  !coeff. for inc.1, xyz(1:3) means propa. direc. of the inc.
    real(8)             :: c2_inc2_xyz(3)  !coeff. for inc.2, xyz(1:3) means propa. direc. of the inc.
    real(8),allocatable :: ex_y(:,:,:),c1_ex_y(:,:,:),c2_ex_y(:,:,:),ex_z(:,:,:),c1_ex_z(:,:,:),c2_ex_z(:,:,:) !e
    real(8),allocatable :: ey_z(:,:,:),c1_ey_z(:,:,:),c2_ey_z(:,:,:),ey_x(:,:,:),c1_ey_x(:,:,:),c2_ey_x(:,:,:) !e
    real(8),allocatable :: ez_x(:,:,:),c1_ez_x(:,:,:),c2_ez_x(:,:,:),ez_y(:,:,:),c1_ez_y(:,:,:),c2_ez_y(:,:,:) !e
    integer             :: iex_y_is(3),iex_y_ie(3),iex_z_is(3),iex_z_ie(3)                                     !e
    integer             :: iey_z_is(3),iey_z_ie(3),iey_x_is(3),iey_x_ie(3)                                     !e
    integer             :: iez_x_is(3),iez_x_ie(3),iez_y_is(3),iez_y_ie(3)                                     !e
    real(8),allocatable :: hx_y(:,:,:),c1_hx_y(:,:,:),c2_hx_y(:,:,:),hx_z(:,:,:),c1_hx_z(:,:,:),c2_hx_z(:,:,:) !h
    real(8),allocatable :: hy_z(:,:,:),c1_hy_z(:,:,:),c2_hy_z(:,:,:),hy_x(:,:,:),c1_hy_x(:,:,:),c2_hy_x(:,:,:) !h
    real(8),allocatable :: hz_x(:,:,:),c1_hz_x(:,:,:),c2_hz_x(:,:,:),hz_y(:,:,:),c1_hz_y(:,:,:),c2_hz_y(:,:,:) !h
    integer             :: ihx_y_is(3),ihx_y_ie(3),ihx_z_is(3),ihx_z_ie(3)                                     !h
    integer             :: ihy_z_is(3),ihy_z_ie(3),ihy_x_is(3),ihy_x_ie(3)                                     !h
    integer             :: ihz_x_is(3),ihz_x_ie(3),ihz_y_is(3),ihz_y_ie(3)                                     !h
    real(8),allocatable :: ex_s(:,:,:),ey_s(:,:,:),ez_s(:,:,:)                   !e for save
    real(8),allocatable :: hx_s(:,:,:),hy_s(:,:,:),hz_s(:,:,:)                   !h for save
    real(8),allocatable :: c2_jx(:,:,:),c2_jy(:,:,:),c2_jz(:,:,:)                !coeff. for general curr. dens.
    integer             :: num_ld                                                !LD: number of LD media
    integer             :: max_pole_num_ld                                       !LD: maximum of pole_num_ld
    integer,allocatable :: media_ld(:)                                           !LD: imedia number for drude model
                                                                                 !    (media id)
    integer,allocatable :: idx_ld(:,:,:,:),idy_ld(:,:,:,:),idz_ld(:,:,:,:)       !LD: id for each component
                                                                                 !    (x,y,z,num_ld)
    real(8),allocatable :: rjx_ld(:,:,:,:,:),rjy_ld(:,:,:,:,:),rjz_ld(:,:,:,:,:) !LD: poparization current density J
                                                                                 !    (x,y,z,max_pole_num_ld,num_ld)
    real(8),allocatable :: rjx_sum_ld(:,:,:),rjy_sum_ld(:,:,:),rjz_sum_ld(:,:,:) !LD: sum of J
                                                                                 !    (x,y,z)
    real(8),allocatable :: px_ld(:,:,:,:,:),py_ld(:,:,:,:,:),pz_ld(:,:,:,:,:)    !LD: poparization vector P
                                                                                 !    (x,y,z,max_pole_num_ld,num_ld)
    real(8),allocatable :: px_sum_ld(:,:,:),py_sum_ld(:,:,:),pz_sum_ld(:,:,:)    !LD: sum of P
                                                                                 !    (x,y,z)
    real(8),allocatable :: c1_j_ld(:,:),c2_j_ld(:,:),c3_j_ld(:,:)                !LD: coefficient for J
                                                                                 !    (max_pole_num_ld,num_ld)
    real(8),allocatable :: rmedia(:,:,:)                                !Material information for tmp.
    real(8),allocatable :: time_lr(:)                                   !LR: time
    integer             :: iter_lr                                      !LR: time iteration for save
    real(8),allocatable :: fr_lr(:,:)                                   !LR: Re[f]
    real(8),allocatable :: fi_lr(:,:)                                   !LR: Im[f]
    real(8),allocatable :: px_lr(:,:,:), py_lr(:,:,:), pz_lr(:,:,:)     !LR: poparization vector
    real(8),allocatable :: dip_lr(:,:)                                  !LR: dipolemoment
    real(8),allocatable :: rjx_lr(:,:,:),rjy_lr(:,:,:),rjz_lr(:,:,:)    !LR: poparization current density
    real(8),allocatable :: curr_lr(:,:)                                 !LR: average current density
    real(8),allocatable :: e_lr(:,:)                                    !LR: average electric field
    real(8),allocatable :: er_lr(:,:)                                   !LR: Re[e_lr]
    real(8),allocatable :: ei_lr(:,:)                                   !LR: Im[e_lr]
  end type ls_fdtd_eh
  
  private :: calc_es_and_hs
  private :: allocate_poynting
  private :: calc_poynting_vector
  private :: calc_poynting_vector_div
  integer,private,parameter :: unit2=3000
  character(11),private,parameter :: file_unit2='ttm_rt.data'

contains
  
  !===========================================================================================
  != initialize eh-FDTD ======================================================================
  subroutine eh_init(fs,fe)
    use salmon_global,   only: nt_em,al_em,dl_em,num_rgrid_em,dt_em,boundary_em,yn_periodic,base_directory,&
                               media_num,shape_file,epsilon_em,mu_em,sigma_em,media_type,&
                               pole_num_ld,omega_p_ld,f_ld,gamma_ld,omega_ld,&
                               obs_num_em,obs_loc_em,obs_plane_ene_em,yn_obs_plane_integral_em,&
                               media_id_pml,media_id_source1,media_id_source2,&
                               wave_input,trans_longi,e_impulse,nenergy,&
                               source_loc1,ek_dir1,epdir_re1,epdir_im1,ae_shape1,&
                               phi_cep1,I_wcm2_1,E_amplitude1,&
                               source_loc2,ek_dir2,epdir_re2,epdir_im2,ae_shape2,&
                               phi_cep2,I_wcm2_2,E_amplitude2,&
                               bloch_k_em,bloch_real_imag_em,yn_make_shape
    use inputoutput,     only: utime_from_au,ulength_from_au,uenergy_from_au,unit_system,&
                               uenergy_to_au,ulength_to_au,ucharge_to_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_bcast,comm_sync_all
    use structures,      only: s_fdtd_system
    use phys_constants,  only: cspeed_au
    use math_constants,  only: pi
    use common_maxwell,  only: set_coo_em,find_point_em,input_shape_em
    use ttm           ,  only: use_ttm, init_ttm_parameters, init_ttm_grid, init_ttm_alloc
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    integer                           :: ii,ij,ik,ix,iy,iz,icount,icount_ld,iflag_lr,iflag_pml
    real(8)                           :: dt_cfl
    character(1)                      :: dir
    character(2)                      :: plane_name
    character(16)                     :: tmp_name1,tmp_name2
    character(128)                    :: save_name
    
    !*** set initial parameter and value **********************************************************************!
    fe%Nd        = 1
    fe%iter_sta  = 1
    fe%iter_end  = nt_em
    fe%ifn       = 600
    fe%ipml_l    = 8
    fe%pml_m     = 4.0d0
    fe%pml_r     = 1.0d-7
    if( ( (al_em(1)*al_em(2)*al_em(3) >0.0d0).and.               &
          (dl_em(1)*dl_em(2)*dl_em(3) >0.0d0).and.               &
          (num_rgrid_em(1)*num_rgrid_em(2)*num_rgrid_em(3) >0) ) &
        .or.                                                     &
        ( (al_em(1)*al_em(2)*al_em(3)==0.0d0).and.               &
          (dl_em(1)*dl_em(2)*dl_em(3)==0.0d0)                  ) &
        .or.                                                     &
        ( (al_em(1)*al_em(2)*al_em(3)==0.0d0).and.               &
          (num_rgrid_em(1)*num_rgrid_em(2)*num_rgrid_em(3)==0) ) &
        .or.                                                     &
        ( (dl_em(1)*dl_em(2)*dl_em(3)==0.0d0).and.               &
          (num_rgrid_em(1)*num_rgrid_em(2)*num_rgrid_em(3)==0) ) ) then
      if(comm_is_root(nproc_id_global)) &
      write(*,*) "Only two of al_em, dl_em, and num_rgrid_em must be set."
      stop
    end if
    if(al_em(1)*al_em(2)*al_em(3)==0.0d0) then
      fs%rlsize(:) = dl_em(:) * dble(num_rgrid_em(:))
      fs%hgs(:)    = dl_em(:)
    else
      if(num_rgrid_em(1)*num_rgrid_em(2)*num_rgrid_em(3)>0) then
        fs%rlsize(:) = al_em(:)
        fs%hgs(:)    = fs%rlsize(:) / dble(num_rgrid_em(:))
      else
        fs%rlsize(:) = al_em(:)
        fs%hgs(:)    = dl_em(:)
      end if
    end if
    do ii=1,3
    do ij=1,2
      select case(boundary_em(ii,ij))
      case('default')
        if(yn_periodic=='n') then
          fs%a_bc(ii,ij) = 'pml'
          iflag_pml      = 1
        elseif(yn_periodic=='y') then
          fs%a_bc(ii,ij) = 'periodic'
        end if
      case('abc')
        fs%a_bc(ii,ij) = 'pml'
        iflag_pml      = 1
      case('pec')
        if(comm_is_root(nproc_id_global).and.(yn_periodic=='y')) &
        write(*,*) "For yn_periodic = y, boundary_em must be default or abc."
        stop
        fs%a_bc(ii,ij) = 'pec'
      case('periodic')
        if(comm_is_root(nproc_id_global).and.(yn_periodic=='n')) &
        write(*,*) "For yn_periodic = n, boundary_em must be default, abc, or pec."
        stop
      end select
    end do
    end do
    select case(unit_system)
    case('au','a.u.')
      fe%uVperm_from_au=1.0d0
      fe%uAperm_from_au=1.0d0
    case('A_eV_fs')
      !see E_amplitude1 or E_amplitude2 in src/io/iunputoutput.f90
      fe%uVperm_from_au=1/(uenergy_to_au/ulength_to_au/ucharge_to_au)
      fe%uAperm_from_au=fe%uVperm_from_au
    end select
    
    !*** prepare mpi, gird, and sendrecv environments *********************************************************!
    call eh_mpi_grid_sr(fs,fe)
    
    !*** set coordinate ***************************************************************************************!
    do ii=1,3
      if(mod(int(fs%rlsize(ii)/fs%hgs(ii)+1.d-12),2)==1)then
        fe%ioddeven(ii)=1
      else
        fe%ioddeven(ii)=2
      end if
    end do 
    allocate(fe%coo(minval(fs%lg%is(:))-fe%Nd:maxval(fs%lg%ie(:))+fe%Nd,3))
    call set_coo_em(fe%Nd,fe%ioddeven(:),fs%lg%is(:),fs%lg%ie(:),fs%hgs(:),fe%coo(:,:),yn_periodic)
    
    !*** make or input fdtd shape *****************************************************************************!
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'i3d',i3d=fs%imedia)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rmedia)
    if(media_num>0) then
      if(comm_is_root(nproc_id_global)) write(*,*)
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
      
      !make or input
      if(yn_make_shape=='y') then
        !make input file
        call eh_make_shape(fs,fe)
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "shape file is made by input keywords."
        end if
      elseif(yn_make_shape=='n') then
        !check file format and input shape file
        if(index(shape_file,".cube", back=.true.)/=0) then
          if(comm_is_root(nproc_id_global)) then
            write(*,*) "shape file is inputed by .cube format."
          end if
          call input_shape_em(shape_file,fe%ifn,fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fs%imedia,'cu')
        elseif(index(shape_file,".mp", back=.true.)/=0) then
          if(comm_is_root(nproc_id_global)) then
            write(*,*) "shape file is inputed by .mp format."
            write(*,*) "This version works for only .cube format.."
          end if
          stop
        else
          if(comm_is_root(nproc_id_global)) then
            write(*,*) "shape file must be .cube or .mp formats."
          end if
          stop
        end if
      end if
      
      !send and receive for shape file
      fe%rmedia(:,:,:)=dble(fs%imedia(:,:,:))
      call eh_sendrecv(fs,fe,'r')
      fs%imedia(:,:,:)=int(fe%rmedia(:,:,:)+1d-3)
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    end if
    
    !*** set and check dt *************************************************************************************!
    dt_cfl=1.0d0/( &
           cspeed_au*sqrt( (1.0d0/fs%hgs(1))**2.0d0+(1.0d0/fs%hgs(2))**2.0d0+(1.0d0/fs%hgs(3))**2.0d0 ) &
           )
    if(dt_em==0.0d0) then
      dt_em=dt_cfl*0.99d0
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "**************************"
        write(*,*) "From CFL condition, dt_em is determined by", dt_em*utime_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
        write(*,*) "**************************"
      end if
    elseif(dt_em>=dt_cfl) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "**************************"
        write(*,*) "To sufficient CFL condition, dt_em must be set smaller than", dt_cfl*utime_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
        write(*,*) "**************************"
      end if
      stop
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "**************************"
        write(*,*) "dt_em =", dt_em*utime_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
        write(*,*) "**************************"
      end if
    end if
    call comm_bcast(dt_em,nproc_group_global)
    
    !*** basic allocation in eh-FDTD **************************************************************************!
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ex_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ex_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ex_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ex_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ey_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ey_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ey_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ey_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ez_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ez_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ez_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ez_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hx_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hx_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hx_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hx_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hy_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hy_z)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hy_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hy_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hz_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hz_x)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hz_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hz_y)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jx)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jy)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jz)
    
    !*** prepare Lorentz-Drude ********************************************************************************!
    fe%num_ld=0
    do ii=0,media_num
      select case(media_type(ii))
      case('lorentz-drude')
        fe%num_ld=fe%num_ld+1
      end select
    end do
    if(fe%num_ld>0) then
      !set counter, make media_ld, make max_pole_num_ld, and check pole_num_ld condition
      icount_ld=1; fe%max_pole_num_ld=0;
      allocate(fe%media_ld(fe%num_ld))
      fe%media_ld(:)=0;
      do ii=0,media_num
        select case(media_type(ii))
        case('lorentz-drude')
          fe%media_ld(icount_ld)=ii
          icount_ld=icount_ld+1
          if(fe%max_pole_num_ld<pole_num_ld(ii)) fe%max_pole_num_ld=pole_num_ld(ii)
          if(pole_num_ld(ii)<=0) then
            if(comm_is_root(nproc_id_global)) &
              write(*,*) "For media_type = lorentz-drude, pole_num_ld must be equal to or larger than 1."
            stop
          end if
        end select
      end do
      
      !reset counter
      icount_ld=1
      
      !allocate Lorentz-Drude variables
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'i4d',num4d=fe%num_ld,i4d=fe%idx_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'i4d',num4d=fe%num_ld,i4d=fe%idy_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'i4d',num4d=fe%num_ld,i4d=fe%idz_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%rjx_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%rjy_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%rjz_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjx_sum_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjy_sum_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjz_sum_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%px_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%py_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r5d',num4d=fe%max_pole_num_ld,num5d=fe%num_ld,r5d=fe%pz_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%px_sum_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%py_sum_ld)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%pz_sum_ld)
      allocate(fe%c1_j_ld(fe%max_pole_num_ld,fe%num_ld),&
               fe%c2_j_ld(fe%max_pole_num_ld,fe%num_ld),&
               fe%c3_j_ld(fe%max_pole_num_ld,fe%num_ld))
      fe%c1_j_ld(:,:)=0.0d0; fe%c2_j_ld(:,:)=0.0d0; fe%c3_j_ld(:,:)=0.0d0;
    end if
    
    !*** set fdtd coeffient ***********************************************************************************!
    allocate(fe%rep(0:media_num),fe%rmu(0:media_num),fe%sig(0:media_num))
    fe%rep(:)=1.0d0; fe%rmu(:)=1.0d0; fe%sig(:)=0.0d0;
    do ii=0,media_num
      fe%rep(ii)=epsilon_em(ii); fe%rmu(ii)=mu_em(ii); fe%sig(ii)=sigma_em(ii);
    end do
    do ii=0,media_num
      call eh_coeff
    end do

    !*** check TTM On/Off, and initialization *****************************************************************!
    call init_ttm_parameters( dt_em )
    if( use_ttm )then
       call init_ttm_grid( fs%hgs, fs%mg%is_array, fs%mg%is, fs%mg%ie, fs%imedia )
       call init_ttm_alloc( fs%srg_ng, fs%mg )
       if( comm_is_root(nproc_id_global) )then
          open(unit2,file=file_unit2)
          write(unit2,'("#",99(1X,I0,":",A))') &
               1, "Time[fs]",                &
               2, "EM_energy[eV]",           &
               3, "EM_energy_in_medium[eV]", &
               4, "T_ele[K]",                &
               5, "T_lat[K]"
          close(unit2)
       end if
    end if !use_ttm

    !*** write media information ******************************************************************************!
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      write(*,'(A,I6)')          ' media_num   = ',media_num
      do ii=0,media_num
        write(*,*) "=========================="
        write(*,'(A,I3,A)')      ' id ',ii, ':'
        call eh_check_media_type(fe%rep(ii),fe%rmu(ii),fe%sig(ij),media_type(ii),tmp_name1)
        write(*,'(A,A)')         ' media_type  =  ', trim(tmp_name1)
        select case(media_type(ii))
        case('lorentz-drude')
          write(*,'(A,I6)')      ' pole_num_ld = ', pole_num_ld(ii)
          write(*,'(A,ES12.5)')  ' omega_p_ld  = ', omega_p_ld(ii)*uenergy_from_au
          do ij=1,pole_num_ld(ii)
            if(ij==1) then
              write(*,'(A,ES12.5)')  ' f_ld        = ', f_ld(ii,ij)
            else
              write(*,'(A,ES12.5)')  '               ', f_ld(ii,ij)
            end if
          end do
          do ij=1,pole_num_ld(ii)
            if(ij==1) then
              write(*,'(A,ES12.5)')  ' gamma_ld    = ', gamma_ld(ii,ij)*uenergy_from_au
            else
              write(*,'(A,ES12.5)')  '               ', gamma_ld(ii,ij)*uenergy_from_au
            end if
          end do
          do ij=1,pole_num_ld(ii)
            if(ij==1) then
              write(*,'(A,ES12.5)')  ' omega_ld    = ', omega_ld(ii,ij)*uenergy_from_au
            else
              write(*,'(A,ES12.5)')  '               ', omega_ld(ii,ij)*uenergy_from_au
            end if
          end do
        end select
        write(*,'(A,ES12.5)')  ' epsilon_em  = ', fe%rep(ii)
        write(*,'(A,ES12.5)')  ' mu_em       = ', fe%rmu(ii)
        write(*,'(A,ES12.5)')  ' sigma_em    = ', fe%sig(ii)
      end do
      write(*,*) "**************************"
    end if
    
    !*** set calculation area *********************************************************************************!
    fe%iex_y_is(:)=fs%mg%is(:); fe%iex_y_ie(:)=fs%mg%ie(:);
    fe%iex_z_is(:)=fs%mg%is(:); fe%iex_z_ie(:)=fs%mg%ie(:);
    fe%iey_z_is(:)=fs%mg%is(:); fe%iey_z_ie(:)=fs%mg%ie(:);
    fe%iey_x_is(:)=fs%mg%is(:); fe%iey_x_ie(:)=fs%mg%ie(:);
    fe%iez_x_is(:)=fs%mg%is(:); fe%iez_x_ie(:)=fs%mg%ie(:);
    fe%iez_y_is(:)=fs%mg%is(:); fe%iez_y_ie(:)=fs%mg%ie(:);
    fe%ihx_y_is(:)=fs%mg%is(:); fe%ihx_y_ie(:)=fs%mg%ie(:);
    fe%ihx_z_is(:)=fs%mg%is(:); fe%ihx_z_ie(:)=fs%mg%ie(:);
    fe%ihy_z_is(:)=fs%mg%is(:); fe%ihy_z_ie(:)=fs%mg%ie(:);
    fe%ihy_x_is(:)=fs%mg%is(:); fe%ihy_x_ie(:)=fs%mg%ie(:);
    fe%ihz_x_is(:)=fs%mg%is(:); fe%ihz_x_ie(:)=fs%mg%ie(:);
    fe%ihz_y_is(:)=fs%mg%is(:); fe%ihz_y_ie(:)=fs%mg%ie(:);
    if((fs%a_bc(1,1)=='pml').and.(fs%mg%is(1)==fs%lg%is(1))) then !x, bottom
      fe%iey_x_is(1)=fs%mg%is(1)+1; fe%iez_x_is(1)=fs%mg%is(1)+1;
    end if
    if((fs%a_bc(1,2)=='pml').and.(fs%mg%ie(1)==fs%lg%ie(1))) then !x, top
      fe%iex_y_ie(1)=fs%mg%ie(1)-1; fe%iex_z_ie(1)=fs%mg%ie(1)-1;
      fe%iey_x_ie(1)=fs%mg%ie(1)-1; fe%iez_x_ie(1)=fs%mg%ie(1)-1;
      fe%ihy_z_ie(1)=fs%mg%ie(1)-1; fe%ihy_x_ie(1)=fs%mg%ie(1)-1;
      fe%ihz_x_ie(1)=fs%mg%ie(1)-1; fe%ihz_y_ie(1)=fs%mg%ie(1)-1;
    end if
    if((fs%a_bc(2,1)=='pml').and.(fs%mg%is(2)==fs%lg%is(2))) then !y, bottom
      fe%iex_y_is(2)=fs%mg%is(2)+1; fe%iez_y_is(2)=fs%mg%is(2)+1;
    end if
    if((fs%a_bc(2,2)=='pml').and.(fs%mg%ie(2)==fs%lg%ie(2))) then !y, top
      fe%iex_y_ie(2)=fs%mg%ie(2)-1; fe%iey_z_ie(2)=fs%mg%ie(2)-1;
      fe%iey_x_ie(2)=fs%mg%ie(2)-1; fe%iez_y_ie(2)=fs%mg%ie(2)-1;
      fe%ihx_y_ie(2)=fs%mg%ie(2)-1; fe%ihx_z_ie(2)=fs%mg%ie(2)-1;
      fe%ihz_x_ie(2)=fs%mg%ie(2)-1; fe%ihz_y_ie(2)=fs%mg%ie(2)-1;
    end if
    if((fs%a_bc(3,1)=='pml').and.(fs%mg%is(3)==fs%lg%is(3))) then !z, bottom
      fe%iex_z_is(3)=fs%mg%is(3)+1; fe%iey_z_is(3)=fs%mg%is(3)+1;
    end if
    if((fs%a_bc(3,2)=='pml').and.(fs%mg%ie(3)==fs%lg%ie(3))) then !z, top
      fe%iex_z_ie(3)=fs%mg%ie(3)-1; fe%iey_z_ie(3)=fs%mg%ie(3)-1;
      fe%iez_x_ie(3)=fs%mg%ie(3)-1; fe%iez_y_ie(3)=fs%mg%ie(3)-1;
      fe%ihx_y_ie(3)=fs%mg%ie(3)-1; fe%ihx_z_ie(3)=fs%mg%ie(3)-1;
      fe%ihy_z_ie(3)=fs%mg%ie(3)-1; fe%ihy_x_ie(3)=fs%mg%ie(3)-1;
    end if
    
    !*** set pml **********************************************************************************************!
    call eh_set_pml(1,fe%c1_ey_x,fe%c2_ey_x,fe%c1_ez_x,fe%c2_ez_x,&
                      fe%c1_hy_x,fe%c2_hy_x,fe%c1_hz_x,fe%c2_hz_x) !x direction
    call eh_set_pml(2,fe%c1_ez_y,fe%c2_ez_y,fe%c1_ex_y,fe%c2_ex_y,&
                      fe%c1_hz_y,fe%c2_hz_y,fe%c1_hx_y,fe%c2_hx_y) !y direction
    call eh_set_pml(3,fe%c1_ex_z,fe%c2_ex_z,fe%c1_ey_z,fe%c2_ey_z,&
                      fe%c1_hx_z,fe%c2_hx_z,fe%c1_hy_z,fe%c2_hy_z) !z direction
    if(iflag_pml==1) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        do ii=1,3
          if(ii==1) then
            dir='x'
          elseif(ii==2) then
            dir='y'
          elseif(ii==3) then
            dir='z'
          end if
          if(fs%a_bc(ii,1)=='pml') write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                                   ' PML has been set for ',dir,'-direction: ',&
                                   fe%coo(fs%lg%is(ii),ii)*ulength_from_au,' to ',&
                                   fe%coo(fs%lg%is(ii)+fe%ipml_l,ii)*ulength_from_au,'.'
          if(fs%a_bc(ii,2)=='pml') write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                                   ' PML has been set for ',dir,'-direction: ',&
                                   fe%coo(fs%lg%ie(ii)-fe%ipml_l,ii)*ulength_from_au,' to ',&
                                   fe%coo(fs%lg%ie(ii),ii)*ulength_from_au,'.'
        end do
        write(*,*) "**************************"
      end if
    end if
    
    !*** set bloch boundary condition *************************************************************************!
    allocate( fe%c_ebloch_x(fs%mg%is_array(1):fs%mg%ie_array(1)), &
              fe%c_ebloch_y(fs%mg%is_array(2):fs%mg%ie_array(2)), &
              fe%c_ebloch_z(fs%mg%is_array(3):fs%mg%ie_array(3)), &
              fe%c_hbloch_x(fs%mg%is_array(1):fs%mg%ie_array(1)), &
              fe%c_hbloch_y(fs%mg%is_array(2):fs%mg%ie_array(2)), &
              fe%c_hbloch_z(fs%mg%is_array(3):fs%mg%ie_array(3)) )
    fe%c_ebloch_x(:) = 1.0d0; fe%c_ebloch_y(:) = 1.0d0; fe%c_ebloch_z(:) = 1.0d0;
    fe%c_hbloch_x(:) = 1.0d0; fe%c_hbloch_y(:) = 1.0d0; fe%c_hbloch_z(:) = 1.0d0;
    if( (sum(abs(bloch_k_em(:)))/=0.0d0)    &
        .or.(bloch_real_imag_em(1)=='imag') &
        .or.(bloch_real_imag_em(2)=='imag') &
        .or.(bloch_real_imag_em(3)=='imag') ) then !--- use bloch boundary condition --------------------------!
      !check yn_periodic
      if(yn_periodic/='y') then
        if(comm_is_root(nproc_id_global)) &
        write(*,*) "When bloch boundary conditions for electromagnetic problem are used by |bloch_k_em| > 0,"
        write(*,*) "yn_periodic must be y."
        stop
      end if
      
      !set x direction
      if(fs%mg%is(1)==fs%lg%is(1)) then !bottom
        select case(bloch_real_imag_em(1))
        case('real')
          fe%c_hbloch_x( fs%mg%is(1) - 1 ) = cos( -bloch_k_em(1)*fs%rlsize(1) );
        case('imag')
          fe%c_hbloch_x( fs%mg%is(1) - 1 ) = sin( -bloch_k_em(1)*fs%rlsize(1) );
        end select
      end if
      if(fs%mg%ie(1)==fs%lg%ie(1)) then !top
        select case(bloch_real_imag_em(1))
        case('real')
          fe%c_ebloch_x( fs%mg%ie(1) + 1 ) = cos(  bloch_k_em(1)*fs%rlsize(1) );
        case('imag')
          fe%c_ebloch_x( fs%mg%ie(1) + 1 ) = sin(  bloch_k_em(1)*fs%rlsize(1) );
        end select
      end if
      
      !set y direction
      if(fs%mg%is(2)==fs%lg%is(2)) then !bottom
        select case(bloch_real_imag_em(2))
        case('real')
          fe%c_hbloch_y( fs%mg%is(2) - 1 ) = cos( -bloch_k_em(2)*fs%rlsize(2) );
        case('imag')
          fe%c_hbloch_y( fs%mg%is(2) - 1 ) = sin( -bloch_k_em(2)*fs%rlsize(2) );
        end select
      end if
      if(fs%mg%ie(2)==fs%lg%ie(2)) then !top
        select case(bloch_real_imag_em(2))
        case('real')
          fe%c_ebloch_y( fs%mg%ie(2) + 1 ) = cos(  bloch_k_em(2)*fs%rlsize(2) );
        case('imag')
          fe%c_ebloch_y( fs%mg%ie(2) + 1 ) = sin(  bloch_k_em(2)*fs%rlsize(2) );
        end select
      end if
      
      !set z direction
      if(fs%mg%is(3)==fs%lg%is(3)) then !bottom
        select case(bloch_real_imag_em(3))
        case('real')
          fe%c_hbloch_z( fs%mg%is(3) - 1 ) = cos( -bloch_k_em(3)*fs%rlsize(3) );
        case('imag')
          fe%c_hbloch_z( fs%mg%is(3) - 1 ) = sin( -bloch_k_em(3)*fs%rlsize(3) );
        end select
      end if
      if(fs%mg%ie(3)==fs%lg%ie(3)) then !top
        select case(bloch_real_imag_em(3))
        case('real')
          fe%c_ebloch_z( fs%mg%ie(3) + 1 ) = cos(  bloch_k_em(3)*fs%rlsize(3) );
        case('imag')
          fe%c_ebloch_z( fs%mg%ie(3) + 1 ) = sin(  bloch_k_em(3)*fs%rlsize(3) );
        end select
      end if
      
      !output message
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        write(*,'(A)'        ) ' Bloch boundary conditions for electromagnetic problem have been set:'
        write(*,'(A,3ES12.5)') ' bloch_k_em =',bloch_k_em(:)/ulength_from_au
        write(*,*) "**************************"
      end if
    end if
    
    !*** prepare observation **********************************************************************************!
    if(obs_num_em>0) then
      !set initial
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_s)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_s)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_s)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_s)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_s)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_s)
      allocate(fe%iobs_po_id(obs_num_em,3)) !1:x,        2:y,        3:z
      allocate(fe%iobs_po_pe(obs_num_em))
      allocate(fe%iobs_li_pe(obs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
      allocate(fe%iobs_pl_pe(obs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
      fe%iobs_po_id(:,:)=0; fe%iobs_po_pe(:)=0; fe%iobs_li_pe(:,:)=0; fe%iobs_pl_pe(:,:)=0; 
      fe%e_max=0.0d0; fe%h_max=0.0d0;
      
      !search observation point
      do ii=1,obs_num_em
        call find_point_em(obs_loc_em(ii,:),fe%iobs_po_id(ii,:),&
                           fe%iobs_po_pe(ii),fe%iobs_li_pe(ii,:),fe%iobs_pl_pe(ii,:),fs%mg%is(:),fs%mg%ie(:),&
                           minval(fs%lg%is)-fe%Nd,maxval(fs%lg%ie)+fe%Nd,fe%coo(:,:))
      end do
      
      !prepare for obs_plane_ene_em option(spatial distribution at each energy point)
      allocate(fe%iobs_num_ene(obs_num_em))
      fe%iobs_num_ene(:)=0
      do ii=1,obs_num_em
        icount = 0
        do ij=1,size(obs_plane_ene_em,2)
          if(obs_plane_ene_em(ii,ij)>=0.0d0) then
            icount = icount + 1
          end if
        end do
        fe%iobs_num_ene(ii) = icount
      end do
      if(sum(fe%iobs_num_ene(:))>0) then
        allocate( fe%obs_ex_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ex_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ex_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ey_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ey_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ey_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ez_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ez_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_ez_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hx_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hx_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hx_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hy_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hy_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hy_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hz_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hz_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_hz_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2) )
        fe%obs_ex_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ex_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ex_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ey_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ey_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ey_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ez_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ez_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_ez_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hx_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hx_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hx_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hy_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hy_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hy_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hz_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hz_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_hz_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
      end if
      
      !write information
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        if(obs_num_em==1) then
          write(*,*) "Observation point is placed at"
        else
          write(*,*) "Observation points are placed at"
        end if
      end if
      do ii=1,obs_num_em
        if(fe%iobs_po_pe(ii)==1) then
          !coordinate
          write(*,'(I3,A,3ES14.5)') ii,":",(fe%coo(fe%iobs_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
          
          !ex
          write(*,'(6X,A)') "Ex is averaged from two spatial grids whose material parameters are:"
          ij = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)  )
          call eh_check_media_type(fe%rep(ij),fe%rmu(ij),fe%sig(ij),media_type(ij),tmp_name1)
          ik = fs%imedia(fe%iobs_po_id(ii,1)-1,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)  )
          call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
          
          !ey
          write(*,'(6X,A)') "Ey is averaged from two spatial grids whose material parameters are:"
          ik = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)-1,fe%iobs_po_id(ii,3)  )
          call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
          
          !ez
          write(*,'(6X,A)') "Ez is averaged from two spatial grids whose material parameters are:"
          ik = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)-1)
          call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
        end if
        call comm_sync_all
      end do
      call comm_sync_all
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "**************************"
        do ii=1,obs_num_em
          !set header for point data
          write(save_name,*) ii
          save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//'_at_point_rt.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Real time calculation:" 
          write(fe%ifn,'(A)') "# E: Electric field" 
          write(fe%ifn,'(A)') "# H: Magnetic field" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[a.u.]",                &
                  2, "E_x[a.u.]",                 &
                  3, "E_y[a.u.]",                 &
                  4, "E_z[a.u.]",                 &
                  5, "H_x[a.u.]",                 &
                  6, "H_y[a.u.]",                 &
                  7, "H_z[a.u.]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[fs]",                  &
                  2, "E_x[V/Angstrom]",           &
                  3, "E_y[V/Angstrom]",           &
                  4, "E_z[V/Angstrom]",           &
                  5, "H_x[A/Angstrom]",           &
                  6, "H_y[A/Angstrom]",           &
                  7, "H_z[A/Angstrom]"
          end select
          close(fe%ifn)
          
          !set header for plane integral data
          if(yn_obs_plane_integral_em(ii)=='y') then
            do ij=1,3
              !set plane name
              if(ij==1)     then !xy plane
                plane_name='xy';
              elseif(ij==2) then !yz plane
                plane_name='yz';
              elseif(ij==3) then !xz plane
                plane_name='xz';
              end if
              
              !set header
              write(save_name,*) ii
              save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//&
                        '_'//plane_name//'_integral_rt.data'
              open(fe%ifn,file=save_name)
              write(fe%ifn,'(A)') "# Real time calculation:" 
              write(fe%ifn,'(A)') "# IE: Plane integration of electric field" 
              write(fe%ifn,'(A)') "# IH: Plane integration of magnetic field" 
              write(fe%ifn,'(A)') "# IP: Plane integration of Poynting vector" 
              select case(unit_system)
              case('au','a.u.')
                write(fe%ifn,'("#",99(1X,I0,":",A))') &
                      1, "Time[a.u.]",                &
                      2, "IE_x[a.u.]",                &
                      3, "IE_y[a.u.]",                &
                      4, "IE_z[a.u.]",                &
                      5, "IH_x[a.u.]",                &
                      6, "IH_y[a.u.]",                &
                      7, "IH_z[a.u.]",                &
                      8, "IP_x[a.u.]",                &
                      9, "IP_y[a.u.]",                &
                      10,"IP_z[a.u.]"
              case('A_eV_fs')
                write(fe%ifn,'("#",99(1X,I0,":",A))') &
                      1, "Time[fs]",                  &
                      2, "IE_x[V*Angstrom]",          &
                      3, "IE_y[V*Angstrom]",          &
                      4, "IE_z[V*Angstrom]",          &
                      5, "IH_x[A*Angstrom]",          &
                      6, "IH_y[A*Angstrom]",          &
                      7, "IH_z[A*Angstrom]",          &
                      8, "IP_x[eV/fs]",               &
                      9, "IP_y[eV/fs]",               &
                      10,"IP_z[eV/fs]"
              end select
              close(fe%ifn)
            end do
          end if
        end do
      end if
    end if
    
    !*** check incident current source condition **************************************************************!
    select case(wave_input)
    case('source')
      !check linear response
      if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid ae_shape1/2:"
          write(*,*) "For ae_shape1/2 = impulse, wave_input must be default(do not set source)."
        end if
        stop
      end if
      
      !check source1 and source2
      call eh_check_inc(1,ek_dir1,epdir_re1,epdir_im1,fe%inc_dist1)
      call eh_check_inc(2,ek_dir2,epdir_re2,epdir_im2,fe%inc_dist2)
    case('point_hs','point_ss','x-line_hs','x-line_ss','y-line_hs','y-line_ss','z-line_hs','z-line_ss',&
         'point_ip','x-line_ip','y-line_ip','z-line_ip','xy-plane_ip','yz-plane_ip','xz-plane_ip')
      !these selection are for debug
      fe%inc_dist1=wave_input; fe%inc_dist2='none';
      if(comm_is_root(nproc_id_global)) write(*,*) trim(wave_input), " source is used."
    case default
      fe%inc_dist1='none'; fe%inc_dist2='none';
      if(ae_shape1/='impulse'.and.ae_shape2/='impulse') then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid wave_input:"
          write(*,*) "For theory = maxwell, wave_input must be source"
          write(*,*) "or ae_shape1 and/or ae_shape2 must be impulse."
        end if
        stop
      end if
    end select
    
    !*** prepare incident current source **********************************************************************!
    if((fe%inc_dist1=='none').and.(fe%inc_dist2=='none')) then
      fe%inc_num=0
    else
      fe%inc_num=2
    end if
    if(fe%inc_num>0) then
      !set initial
      allocate(fe%inc_po_id(fe%inc_num,3)) !1:x,        2:y,        3:z
      allocate(fe%inc_po_pe(fe%inc_num))
      allocate(fe%inc_li_pe(fe%inc_num,3)) !1:x-line,   2:y-line,   3:z-line
      allocate(fe%inc_pl_pe(fe%inc_num,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
      fe%inc_po_id(:,:)=0; fe%inc_po_pe(:)=0; fe%inc_li_pe(:,:)=0; fe%inc_pl_pe(:,:)=0; 
      do ii=1,3
        fe%c2_inc1_xyz(ii)=(cspeed_au/fe%rep(media_id_source1)*dt_em) &
                            /(1.0d0+2.0d0*pi*fe%sig(media_id_source1)/fe%rep(media_id_source1)*dt_em) &
                            *2.0d0/( fs%hgs(ii)*sqrt(fe%rmu(media_id_source1)/fe%rep(media_id_source1)) )
        fe%c2_inc2_xyz(ii)=(cspeed_au/fe%rep(media_id_source2)*dt_em) &
                            /(1.0d0+2.0d0*pi*fe%sig(media_id_source2)/fe%rep(media_id_source2)*dt_em) &
                            *2.0d0/( fs%hgs(ii)*sqrt(fe%rmu(media_id_source2)/fe%rep(media_id_source2)) )
      end do
      
      !search incident current source point and check others
      if(fe%inc_dist1/='none') then
        ii=1
        call find_point_em(source_loc1(:),fe%inc_po_id(ii,:),&
                           fe%inc_po_pe(ii),fe%inc_li_pe(ii,:),fe%inc_pl_pe(ii,:),fs%mg%is(:),fs%mg%ie(:),&
                           minval(fs%lg%is(:))-fe%Nd,maxval(fs%lg%ie(:))+fe%Nd,fe%coo(:,:))
        select case(wave_input)
        case('point_ip','x-line_ip','y-line_ip','z-line_ip','xy-plane_ip','yz-plane_ip','xz-plane_ip')
        case default
          call eh_check_iw_parameter(ii,phi_cep1,I_wcm2_1,E_amplitude1,ae_shape1)
        end select
      end if
      if(fe%inc_dist2/='none') then
        ii=2
        call find_point_em(source_loc2(:),fe%inc_po_id(ii,:),&
                           fe%inc_po_pe(ii),fe%inc_li_pe(ii,:),fe%inc_pl_pe(ii,:),fs%mg%is(:),fs%mg%ie(:),&
                           minval(fs%lg%is(:))-fe%Nd,maxval(fs%lg%ie(:))+fe%Nd,fe%coo(:,:))
        select case(wave_input)
        case('point_ip','x-line_ip','y-line_ip','z-line_ip','xy-plane_ip','yz-plane_ip','xz-plane_ip')
        case default
          call eh_check_iw_parameter(ii,phi_cep2,I_wcm2_2,E_amplitude2,ae_shape2)
        end select
      end if
      
      !write information
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        if((fe%inc_dist1=='none').or.(fe%inc_dist2=='none')) then
          write(*,*) "Incident current source is placed at"
        else
          write(*,*) "Incident current sources are placed at"
        end if
        if(fe%inc_dist1/='none') then
          ii=1
          write(*,'(I8,A,3ES14.5,A)') ii,":",(fe%coo(fe%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
          write(*,'(A,3ES14.5)') " ek_dir1:",ek_dir1
        end if
        if(fe%inc_dist2/='none') then
          ii=2
          write(*,'(I8,A,3ES14.5,A)') ii,":",(fe%coo(fe%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
          write(*,'(A,3ES14.5)') " ek_dir2:",ek_dir2
        end if
        write(*,*) "**************************"
      end if
    end if
    
    !*** impulsive source option used for debug ***************************************************************!
    select case(wave_input)
    case('point_ip')
      if(fe%inc_po_pe(1)==1) then
        ix=fe%inc_po_id(1,1); iy=fe%inc_po_id(1,2); iz=fe%inc_po_id(1,3);
        fe%ex_y(ix,iy,iz)=epdir_re1(1)/2.0d0
        fe%ex_z(ix,iy,iz)=epdir_re1(1)/2.0d0
        fe%ey_z(ix,iy,iz)=epdir_re1(2)/2.0d0
        fe%ey_x(ix,iy,iz)=epdir_re1(2)/2.0d0
        fe%ez_x(ix,iy,iz)=epdir_re1(3)/2.0d0
        fe%ez_y(ix,iy,iz)=epdir_re1(3)/2.0d0
      end if
    case('x-line_ip')
      if(fe%inc_li_pe(1,1)==1) then
        iy=fe%inc_po_id(1,2); iz=fe%inc_po_id(1,3);
        fe%ex_y(fe%iex_y_is(1):fe%iex_y_ie(1),iy,iz)=epdir_re1(1)/2.0d0
        fe%ex_z(fe%iex_z_is(1):fe%iex_z_ie(1),iy,iz)=epdir_re1(1)/2.0d0
        fe%ey_z(fe%iey_z_is(1):fe%iey_z_ie(1),iy,iz)=epdir_re1(2)/2.0d0
        fe%ey_x(fe%iey_x_is(1):fe%iey_x_ie(1),iy,iz)=epdir_re1(2)/2.0d0
        fe%ez_x(fe%iez_x_is(1):fe%iez_x_ie(1),iy,iz)=epdir_re1(3)/2.0d0
        fe%ez_y(fe%iez_y_is(1):fe%iez_y_ie(1),iy,iz)=epdir_re1(3)/2.0d0
      end if
    case('y-line_ip')
      if(fe%inc_li_pe(1,2)==1) then
        ix=fe%inc_po_id(1,1); iz=fe%inc_po_id(1,3);
        fe%ex_y(ix,fe%iex_y_is(2):fe%iex_y_ie(2),iz)=epdir_re1(1)/2.0d0
        fe%ex_z(ix,fe%iex_z_is(2):fe%iex_z_ie(2),iz)=epdir_re1(1)/2.0d0
        fe%ey_z(ix,fe%iey_z_is(2):fe%iey_z_ie(2),iz)=epdir_re1(2)/2.0d0
        fe%ey_x(ix,fe%iey_x_is(2):fe%iey_x_ie(2),iz)=epdir_re1(2)/2.0d0
        fe%ez_x(ix,fe%iez_x_is(2):fe%iez_x_ie(2),iz)=epdir_re1(3)/2.0d0
        fe%ez_y(ix,fe%iez_y_is(2):fe%iez_y_ie(2),iz)=epdir_re1(3)/2.0d0
      end if
    case('z-line_ip')
      if(fe%inc_li_pe(1,3)==1) then
        ix=fe%inc_po_id(1,1); iy=fe%inc_po_id(1,2);
        fe%ex_y(ix,iy,fe%iex_y_is(3):fe%iex_y_ie(3))=epdir_re1(1)/2.0d0
        fe%ex_z(ix,iy,fe%iex_z_is(3):fe%iex_z_ie(3))=epdir_re1(1)/2.0d0
        fe%ey_z(ix,iy,fe%iey_z_is(3):fe%iey_z_ie(3))=epdir_re1(2)/2.0d0
        fe%ey_x(ix,iy,fe%iey_x_is(3):fe%iey_x_ie(3))=epdir_re1(2)/2.0d0
        fe%ez_x(ix,iy,fe%iez_x_is(3):fe%iez_x_ie(3))=epdir_re1(3)/2.0d0
        fe%ez_y(ix,iy,fe%iez_y_is(3):fe%iez_y_ie(3))=epdir_re1(3)/2.0d0
      end if
    case('xy-plane_ip')
      if(fe%inc_pl_pe(1,1)==1) then
        iz=fe%inc_po_id(1,3)
        do iy=fe%iex_z_is(2),fe%iex_z_ie(2)
        do ix=fe%iex_z_is(1),fe%iex_z_ie(1)
          fe%ex_z(ix,iy,iz)=epdir_re1(1)
        end do
        end do
        do iy=fe%iey_z_is(2),fe%iey_z_ie(2)
        do ix=fe%iey_z_is(1),fe%iey_z_ie(1)
          fe%ey_z(ix,iy,iz)=epdir_re1(2)
        end do
        end do
      end if
    case('yz-plane_ip')
      if(fe%inc_pl_pe(1,2)==1) then
        ix=fe%inc_po_id(1,1)
        do iz=fe%iey_x_is(3),fe%iey_x_ie(3)
        do iy=fe%iey_x_is(2),fe%iey_x_ie(2)
          fe%ey_x(ix,iy,iz)=epdir_re1(2)
        end do
        end do
        do iz=fe%iez_x_is(3),fe%iez_x_ie(3)
        do iy=fe%iez_x_is(2),fe%iez_x_ie(2)
          fe%ez_x(ix,iy,iz)=epdir_re1(3)
        end do
        end do
      end if
    case('xz-plane_ip')
      if(fe%inc_pl_pe(1,3)==1) then
        iy=fe%inc_po_id(1,2)
        do iz=fe%iex_y_is(3),fe%iex_y_ie(3)
        do ix=fe%iex_y_is(1),fe%iex_y_ie(1)
          fe%ex_y(ix,iy,iz)=epdir_re1(1)
        end do
        end do
        do iz=fe%iez_y_is(3),fe%iez_y_ie(3)
        do ix=fe%iez_y_is(1),fe%iez_y_ie(1)
          fe%ez_y(ix,iy,iz)=epdir_re1(3)
        end do
        end do
      end if
    end select
    
    !*** prepare linear response ******************************************************************************!
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      !check condition
      iflag_lr=0
      if(yn_periodic=='y'.and.trans_longi/='tr') iflag_lr=1
      do ii=0,media_num
        if(fe%rep(ii)/=1.0d0.or.fe%rmu(ii)/=1.0d0.or.fe%sig(ii)/=0.0d0) iflag_lr=1
        if(ii==0) then
          select case(media_type(ii))
          case('vacuum')
            continue
          case default
            iflag_lr=1
          end select
        else
          select case(media_type(ii))
          case('lorentz-drude')
            continue
          case default
            iflag_lr=1
          end select
        end if
      end do
      if(iflag_lr==1) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "Invalid input keywords:"
          write(*,*) "When you execute linear response calculation by ae_shape1=impulse and/or ae_shape2=impulse,"
          write(*,*) "epsilon_em and mu_em must be 1.0d0."
          write(*,*) "sigma_em must be 0.0d0."
          write(*,*) "media_type(i) must be lorentz-drude, where i > 0."
          if(yn_periodic=='y') write(*,*) "trans_longi must be tr."
        end if
        stop
      end if
      
      !set initial current density
      if(fe%num_ld>0) then
        do ii=1,fe%num_ld
        do ij=1,pole_num_ld(fe%media_ld(ii))
          do iz=fs%mg%is(3),fs%mg%ie(3)
          do iy=fs%mg%is(2),fs%mg%ie(2)
          do ix=fs%mg%is(1),fs%mg%ie(1)
            if(fe%idx_ld(ix,iy,iz,ii)==1) then
              if(ae_shape1=='impulse') then
                fe%rjx_ld(ix,iy,iz,ij,ii)=fe%rjx_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re1(1)+epdir_im1(1))
              end if
              if(ae_shape2=='impulse') then
                fe%rjx_ld(ix,iy,iz,ij,ii)=fe%rjx_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re2(1)+epdir_im2(1))
              end if
            end if
            if(fe%idy_ld(ix,iy,iz,ii)==1) then
              if(ae_shape1=='impulse') then
                fe%rjy_ld(ix,iy,iz,ij,ii)=fe%rjy_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re1(2)+epdir_im1(2))
              end if
              if(ae_shape2=='impulse') then
                fe%rjy_ld(ix,iy,iz,ij,ii)=fe%rjy_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re2(2)+epdir_im2(2))
              end if
            end if
            if(fe%idz_ld(ix,iy,iz,ii)==1) then
              if(ae_shape1=='impulse') then
                fe%rjz_ld(ix,iy,iz,ij,ii)=fe%rjz_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re1(3)+epdir_im1(3))
              end if
              if(ae_shape2=='impulse') then
                fe%rjz_ld(ix,iy,iz,ij,ii)=fe%rjz_ld(ix,iy,iz,ij,ii) &
                                         -(f_ld(fe%media_ld(ii),ij)*omega_p_ld(fe%media_ld(ii))**2.0d0) &
                                         /(4.0d0*pi)*e_impulse*(epdir_re2(3)+epdir_im2(3))
              end if
            end if
          end do
          end do
          end do
        end do
        end do
      end if
      
      !initialize and allocate
      allocate(fe%time_lr(nt_em))
      fe%time_lr(:)=0.0d0
      fe%iter_lr=1
      allocate(fe%fr_lr(nenergy,3),fe%fi_lr(nenergy,3))
      fe%fr_lr(:,:)=0.0d0; fe%fi_lr(:,:)=0.0d0;
      if(yn_periodic=='n') then
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%px_lr)
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%py_lr)
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%pz_lr)
        allocate(fe%dip_lr(nt_em,3))
        fe%dip_lr(:,:)=0.0d0
      elseif(yn_periodic=='y') then
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjx_lr)
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjy_lr)
        call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rjz_lr)
        allocate(fe%curr_lr(nt_em,3),fe%e_lr(nt_em,3))
        fe%curr_lr(:,:)=0.0d0; fe%e_lr(:,:)=0.0d0;
        allocate(fe%er_lr(0:nenergy,3),fe%ei_lr(0:nenergy,3))
        fe%er_lr(:,:)=0.0d0; fe%ei_lr(:,:)=0.0d0;
      end if
    end if
    
    !*** deallocate unused variables **************************************************************************!
    deallocate(fs%imedia,fe%rmedia);
    
    !*** write start ******************************************************************************************!
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      write(*,*) "FDTD start"
      write(*,*) "**************************"
      write(*,*) "timestep"
      write(*,*) "-------------------------------------------------------"
    end if
    
    return
  contains
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ allocation in eh-FDTD +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_allocate(ista,iend,allocate_mode,num4d,num5d,i3d,r3d,i4d,r5d)
      implicit none
      integer,            intent(in)           :: ista(3),iend(3)
      character(3),       intent(in)           :: allocate_mode
      integer,            intent(in), optional :: num4d,num5d
      integer,allocatable,intent(out),optional :: i3d(:,:,:)
      real(8),allocatable,intent(out),optional :: r3d(:,:,:)
      integer,allocatable,intent(out),optional :: i4d(:,:,:,:)
      real(8),allocatable,intent(out),optional :: r5d(:,:,:,:,:)
      
      select case(allocate_mode)
      case('i3d')
        allocate(i3d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3)))
        i3d(:,:,:)=0
      case('r3d')
        allocate(r3d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3)))
        r3d(:,:,:)=0.0d0
      case('i4d')
        allocate(i4d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3),num4d))
        i4d(:,:,:,:)=0
      case('r5d')
        allocate(r5d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3),num4d,num5d))
        r5d(:,:,:,:,:)=0.0d0
      end select
      
      return
    end subroutine eh_allocate
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ set fdtd coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_coeff
      implicit none
      real(8)  :: c1_e,c2_e_x,c2_e_y,c2_e_z,c1_h,c2_h_x,c2_h_y,c2_h_z,c2_j,&
                  c1_e_mid,c2_e_x_mid,c2_e_y_mid,c2_e_z_mid,c2_j_mid
      
      !prepare coefficient
      c1_e  =(1.0d0-2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em) &
             /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)
      c2_e_x=(cspeed_au/fe%rep(ii)*dt_em) &
             /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(1)
      c2_e_y=(cspeed_au/fe%rep(ii)*dt_em) &
             /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(2)
      c2_e_z=(cspeed_au/fe%rep(ii)*dt_em) &
             /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(3)
      call comm_bcast(c1_e,  nproc_group_global)
      call comm_bcast(c2_e_x,nproc_group_global)
      call comm_bcast(c2_e_y,nproc_group_global)
      call comm_bcast(c2_e_z,nproc_group_global)
      c1_h=1.0d0
      c2_h_x=cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(1)
      c2_h_y=cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(2)
      c2_h_z=cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(3)
      call comm_bcast(c1_h,  nproc_group_global)
      call comm_bcast(c2_h_x,nproc_group_global)
      call comm_bcast(c2_h_y,nproc_group_global)
      call comm_bcast(c2_h_z,nproc_group_global)
      c2_j=(4.0d0*pi/fe%rep(ii)*dt_em) &
           /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)
      call comm_bcast(c2_j,nproc_group_global)
      
      !check media_type
      select case(media_type(ii))
      case('pec')
        c1_e=0.0d0; c2_e_x=0.0d0; c2_e_y=0.0d0; c2_e_z=0.0d0;
      case('lorentz-drude')
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          if(fs%imedia(ix,iy,iz)==ii) then
            if(fs%imedia(ix+1,iy,iz)==ii) then !x
              fe%idx_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
              fe%idx_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
              do ij=1,fe%num_ld
                if(fe%media_ld(ij)==fs%imedia(ix+1,iy,iz)) then
                  fe%idx_ld(ix,iy,iz,ij)=1;
                end if
              end do
            end if
            if(fs%imedia(ix,iy+1,iz)==ii) then !y
              fe%idy_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
              fe%idy_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
              do ij=1,fe%num_ld
                if(fe%media_ld(ij)==fs%imedia(ix,iy+1,iz)) then
                  fe%idy_ld(ix,iy,iz,ij)=1;
                end if
              end do
            end if
            if(fs%imedia(ix,iy,iz+1)==ii) then !z
              fe%idz_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
              fe%idz_ld(ix,iy,iz,icount_ld)=1;
            elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
              do ij=1,fe%num_ld
                if(fe%media_ld(ij)==fs%imedia(ix,iy,iz+1)) then
                  fe%idz_ld(ix,iy,iz,ij)=1;
                end if
              end do
            end if
          end if
        end do
        end do
        end do
        do ij=1,pole_num_ld(ii)
          fe%c1_j_ld(ij,icount_ld)=(1.0d0-gamma_ld(ii,ij)*dt_em/2.0d0) &
                                   / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
          fe%c2_j_ld(ij,icount_ld)=(f_ld(ii,ij)*(omega_p_ld(ii)**2.0d0)*dt_em/(4.0d0*pi)) &
                                   / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
          fe%c3_j_ld(ij,icount_ld)=((omega_ld(ii,ij)**2.0d0)*dt_em) &
                                   / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
        end do
        icount_ld=icount_ld+1
      end select
      
      !set coefficient
      if(ii==0) then
        fe%c1_ex_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ex_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_e_y
        fe%c1_ex_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ex_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_e_z
            
        fe%c1_ey_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ey_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_e_z
        fe%c1_ey_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ey_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_e_x
          
        fe%c1_ez_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ez_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_e_x
        fe%c1_ez_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_e
        fe%c2_ez_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_e_y
          
        fe%c1_hx_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hx_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_h_y
        fe%c1_hx_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hx_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_h_z
          
        fe%c1_hy_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hy_z(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_h_z
        fe%c1_hy_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hy_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_h_x
        
        fe%c1_hz_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hz_x(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_h_x
        fe%c1_hz_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c1_h
        fe%c2_hz_y(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=c2_h_y
                    
        fe%c2_jx(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_j
        fe%c2_jy(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_j
        fe%c2_jz(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))=-c2_j
      else
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          if(fs%imedia(ix,iy,iz)==ii) then
            !ex and jx
            if(fs%imedia(ix+1,iy,iz)==ii) then
              fe%c1_ex_y(ix,iy,iz)=c1_e; fe%c2_ex_y(ix,iy,iz)= c2_e_y;
              fe%c1_ex_z(ix,iy,iz)=c1_e; fe%c2_ex_z(ix,iy,iz)=-c2_e_z;
              fe%c2_jx(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
              fe%c1_ex_y(ix,iy,iz)=c1_e; fe%c2_ex_y(ix,iy,iz)= c2_e_y;
              fe%c1_ex_z(ix,iy,iz)=c1_e; fe%c2_ex_z(ix,iy,iz)=-c2_e_z;
              fe%c2_jx(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
              c1_e_mid  = (1.0d0-2.0d0*pi*fe%sig(fs%imedia(ix+1,iy,iz))/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix+1,iy,iz))/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em)
              c2_e_y_mid= (cspeed_au/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix+1,iy,iz))/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /fs%hgs(2)
              c2_e_z_mid= (cspeed_au/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix+1,iy,iz))/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /fs%hgs(3)
              c2_j_mid  = (4.0d0*pi/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix+1,iy,iz))/fe%rep(fs%imedia(ix+1,iy,iz))*dt_em)
              fe%c1_ex_y(ix,iy,iz)=c1_e_mid; fe%c2_ex_y(ix,iy,iz)= c2_e_y_mid;
              fe%c1_ex_z(ix,iy,iz)=c1_e_mid; fe%c2_ex_z(ix,iy,iz)=-c2_e_z_mid;
              fe%c2_jx(ix,iy,iz)=-c2_j_mid;
              if(fe%num_ld>0) then !LD update
                do ij=1,fe%num_ld
                  if(fe%media_ld(ij)==fs%imedia(ix+1,iy,iz)) then
                    fe%idx_ld(ix,iy,iz,ij)=1;
                  end if
                end do
              end if
            end if
            
            !ey and jy
            if(fs%imedia(ix,iy+1,iz)==ii) then
              fe%c1_ey_z(ix,iy,iz)=c1_e; fe%c2_ey_z(ix,iy,iz)= c2_e_z;
              fe%c1_ey_x(ix,iy,iz)=c1_e; fe%c2_ey_x(ix,iy,iz)=-c2_e_x;
              fe%c2_jy(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
              fe%c1_ey_z(ix,iy,iz)=c1_e; fe%c2_ey_z(ix,iy,iz)= c2_e_z;
              fe%c1_ey_x(ix,iy,iz)=c1_e; fe%c2_ey_x(ix,iy,iz)=-c2_e_x;
              fe%c2_jy(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
              c1_e_mid  = (1.0d0-2.0d0*pi*fe%sig(fs%imedia(ix,iy+1,iz))/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy+1,iz))/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em)
              c2_e_z_mid= (cspeed_au/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy+1,iz))/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /fs%hgs(3)
              c2_e_x_mid= (cspeed_au/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy+1,iz))/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /fs%hgs(1)
              c2_j_mid  = (4.0d0*pi/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy+1,iz))/fe%rep(fs%imedia(ix,iy+1,iz))*dt_em)
              fe%c1_ey_z(ix,iy,iz)=c1_e_mid; fe%c2_ey_z(ix,iy,iz)= c2_e_z_mid;
              fe%c1_ey_x(ix,iy,iz)=c1_e_mid; fe%c2_ey_x(ix,iy,iz)=-c2_e_x_mid;
              fe%c2_jy(ix,iy,iz)=-c2_j_mid;
              if(fe%num_ld>0) then !LD update
                do ij=1,fe%num_ld
                  if(fe%media_ld(ij)==fs%imedia(ix,iy+1,iz)) then
                    fe%idy_ld(ix,iy,iz,ij)=1;
                  end if
                end do
              end if
            end if
            
            !ez and jz
            if(fs%imedia(ix,iy,iz+1)==ii) then
              fe%c1_ez_x(ix,iy,iz)=c1_e; fe%c2_ez_x(ix,iy,iz)= c2_e_x;
              fe%c1_ez_y(ix,iy,iz)=c1_e; fe%c2_ez_y(ix,iy,iz)=-c2_e_y;
              fe%c2_jz(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
              fe%c1_ez_x(ix,iy,iz)=c1_e; fe%c2_ez_x(ix,iy,iz)= c2_e_x;
              fe%c1_ez_y(ix,iy,iz)=c1_e; fe%c2_ez_y(ix,iy,iz)=-c2_e_y;
              fe%c2_jz(ix,iy,iz)=-c2_j;
            elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
              c1_e_mid  = (1.0d0-2.0d0*pi*fe%sig(fs%imedia(ix,iy,iz+1))/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy,iz+1))/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em)
              c2_e_x_mid= (cspeed_au/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy,iz+1))/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /fs%hgs(1)
              c2_e_y_mid= (cspeed_au/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy,iz+1))/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /fs%hgs(2)
              c2_j_mid  = (4.0d0*pi/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em) &
                         /(1.0d0+2.0d0*pi*fe%sig(fs%imedia(ix,iy,iz+1))/fe%rep(fs%imedia(ix,iy,iz+1))*dt_em)
              fe%c1_ez_x(ix,iy,iz)=c1_e_mid; fe%c2_ez_x(ix,iy,iz)= c2_e_x_mid;
              fe%c1_ez_y(ix,iy,iz)=c1_e_mid; fe%c2_ez_y(ix,iy,iz)=-c2_e_y_mid;
              fe%c2_jz(ix,iy,iz)=-c2_j_mid;
              if(fe%num_ld>0) then !LD update
                do ij=1,fe%num_ld
                  if(fe%media_ld(ij)==fs%imedia(ix,iy,iz+1)) then
                    fe%idz_ld(ix,iy,iz,ij)=1;
                  end if
                end do
              end if
            end if
            
            !hx
            fe%c1_hx_y(ix,iy-1:iy,iz-1:iz)=c1_h; fe%c2_hx_y(ix,iy-1:iy,iz-1:iz)=-c2_h_y;
            fe%c1_hx_z(ix,iy-1:iy,iz-1:iz)=c1_h; fe%c2_hx_z(ix,iy-1:iy,iz-1:iz)= c2_h_z;
            
            !hy
            fe%c1_hy_z(ix-1:ix,iy,iz-1:iz)=c1_h; fe%c2_hy_z(ix-1:ix,iy,iz-1:iz)=-c2_h_z;
            fe%c1_hy_x(ix-1:ix,iy,iz-1:iz)=c1_h; fe%c2_hy_x(ix-1:ix,iy,iz-1:iz)= c2_h_x;
            
            !hz
            fe%c1_hz_x(ix-1:ix,iy-1:iy,iz)=c1_h; fe%c2_hz_x(ix-1:ix,iy-1:iy,iz)=-c2_h_x;
            fe%c1_hz_y(ix-1:ix,iy-1:iy,iz)=c1_h; fe%c2_hz_y(ix-1:ix,iy-1:iy,iz)= c2_h_y;
          end if
        end do
        end do
        end do
      end if
      
      return
    end subroutine eh_coeff
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ allocation in eh-FDTD +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_check_media_type(rep,rmu,sig,cm_type,cm_name)
      implicit none
      real(8),       intent(in) :: rep,rmu,sig
      character(16),intent(in)  :: cm_type
      character(16),intent(out) :: cm_name
      
      select case(cm_type)
      case ('vacuum')
        if(rep/=1d0 .or. rmu/=1d0 .or. sig/=0d0) then
          cm_name = 'constant media'
        else
          cm_name = cm_type
        end if
      case default
        cm_name = cm_type
      end select
      
      return
    end subroutine eh_check_media_type
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ set pml +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_set_pml(idir,c1_e1,c2_e1,c1_e2,c2_e2,c1_h1,c2_h1,c1_h2,c2_h2)
      implicit none
      integer,intent(in)  :: idir
      real(8),intent(out) :: c1_e1(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c2_e1(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c1_e2(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c2_e2(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c1_h1(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c2_h1(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c1_h2(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3)),&
                             c2_h2(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                   fs%mg%is_array(2):fs%mg%ie_array(2),&
                                   fs%mg%is_array(3):fs%mg%ie_array(3))
      integer :: ista,iend
      real(8) :: pml_del,s_max_bot,s_max_top
      real(8) :: s_l_bot(fe%ipml_l+1),sh_l_bot(fe%ipml_l), &
                 s_l_top(fe%ipml_l+1),sh_l_top(fe%ipml_l), &
                 c1_pml_bot(fe%ipml_l+1),c2_pml_bot(fe%ipml_l+1),c1_pml_bot_h(fe%ipml_l),c2_pml_bot_h(fe%ipml_l), &
                 c1_pml_top(fe%ipml_l+1),c2_pml_top(fe%ipml_l+1),c1_pml_top_h(fe%ipml_l),c2_pml_top_h(fe%ipml_l)
      
      !set pml conductivity
      pml_del=fs%hgs(idir)
      s_max_bot=-(fe%pml_m+1.0d0)*log(fe%pml_r)/(2.0d0*dble(fe%ipml_l)*pml_del) &
                *cspeed_au/(4.0d0*pi)*sqrt(fe%rep(media_id_pml(idir,1))/fe%rmu(media_id_pml(idir,1)));
      s_max_top=-(fe%pml_m+1.0d0)*log(fe%pml_r)/(2.0d0*dble(fe%ipml_l)*pml_del) &
                *cspeed_au/(4.0d0*pi)*sqrt(fe%rep(media_id_pml(idir,2))/fe%rmu(media_id_pml(idir,2)));
      do ii=1,(fe%ipml_l+1)
        s_l_bot(ii)=s_max_bot*(&
                              (dble(fe%ipml_l)*pml_del-(dble(ii)-1.0d0)*pml_del)/(dble(fe%ipml_l)*pml_del) &
                               )**fe%pml_m;
        s_l_top(ii)=s_max_top*(&
                              (dble(fe%ipml_l)*pml_del-(dble(ii)-1.0d0)*pml_del)/(dble(fe%ipml_l)*pml_del) &
                               )**fe%pml_m;
      end do
      do ii=1,fe%ipml_l
        sh_l_bot(ii)=(fe%rmu(media_id_pml(idir,1))/fe%rep(media_id_pml(idir,1))) &
                     *s_max_bot*(&
                                (dble(fe%ipml_l)*pml_del-(dble(ii)-0.5d0)*pml_del)/(dble(fe%ipml_l)*pml_del) &
                                 )**fe%pml_m;
        sh_l_top(ii)=(fe%rmu(media_id_pml(idir,2))/fe%rep(media_id_pml(idir,2))) &
                     *s_max_top*(&
                                (dble(fe%ipml_l)*pml_del-(dble(ii)-0.5d0)*pml_del)/(dble(fe%ipml_l)*pml_del) &
                                 )**fe%pml_m;
      end do
      
      !set pml coefficient
      do ii=1,(fe%ipml_l+1)
        c1_pml_bot(ii)=(1.0d0-2.0d0*pi*s_l_bot(ii)/fe%rep(media_id_pml(idir,1))*dt_em) &
                       /(1.0d0+2.0d0*pi*s_l_bot(ii)/fe%rep(media_id_pml(idir,1))*dt_em)
        c2_pml_bot(ii)=(cspeed_au/fe%rep(media_id_pml(idir,1))*dt_em) &
                       /(1.0d0+2.0d0*pi*s_l_bot(ii)/fe%rep(media_id_pml(idir,1))*dt_em)/pml_del
        c1_pml_top(ii)=(1.0d0-2.0d0*pi*s_l_top(ii)/fe%rep(media_id_pml(idir,2))*dt_em) &
                       /(1.0d0+2.0d0*pi*s_l_top(ii)/fe%rep(media_id_pml(idir,2))*dt_em)
        c2_pml_top(ii)=(cspeed_au/fe%rep(media_id_pml(idir,2))*dt_em) &
                       /(1.0d0+2.0d0*pi*s_l_top(ii)/fe%rep(media_id_pml(idir,2))*dt_em)/pml_del
      end do
      call comm_bcast(c1_pml_bot,nproc_group_global)
      call comm_bcast(c2_pml_bot,nproc_group_global)
      call comm_bcast(c1_pml_top,nproc_group_global)
      call comm_bcast(c2_pml_top,nproc_group_global)
      do ii=1,fe%ipml_l
        c1_pml_bot_h(ii)=(1.0d0-2.0d0*pi*sh_l_bot(ii)/fe%rmu(media_id_pml(idir,1))*dt_em) &
                         /(1.0d0+2.0d0*pi*sh_l_bot(ii)/fe%rmu(media_id_pml(idir,1))*dt_em)
        c2_pml_bot_h(ii)=(cspeed_au/fe%rmu(media_id_pml(idir,1))*dt_em) &
                         /(1.0d0+2.0d0*pi*sh_l_bot(ii)/fe%rmu(media_id_pml(idir,1))*dt_em)/pml_del
        c1_pml_top_h(ii)=(1.0d0-2.0d0*pi*sh_l_top(ii)/fe%rmu(media_id_pml(idir,2))*dt_em) &
                         /(1.0d0+2.0d0*pi*sh_l_top(ii)/fe%rmu(media_id_pml(idir,2))*dt_em)
        c2_pml_top_h(ii)=(cspeed_au/fe%rmu(media_id_pml(idir,2))*dt_em) &
                         /(1.0d0+2.0d0*pi*sh_l_top(ii)/fe%rmu(media_id_pml(idir,2))*dt_em)/pml_del
      end do
      call comm_bcast(c1_pml_bot_h,nproc_group_global)
      call comm_bcast(c2_pml_bot_h,nproc_group_global)
      call comm_bcast(c1_pml_top_h,nproc_group_global)
      call comm_bcast(c2_pml_top_h,nproc_group_global)
      
      !set pml(bottom)
      if((fs%a_bc(idir,1)=='pml').and.(fs%mg%is(idir)<=(fs%lg%is(idir)+fe%ipml_l))) then
        !e
        iend=fs%lg%is(idir)+fe%ipml_l
        if(fs%mg%ie(idir)<iend) then
          iend=fs%mg%ie(idir)
        end if
        icount=1
        do ii=fs%mg%is(idir),iend
          if(idir==1) then
            c1_e1(ii,:,:)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e1(ii,:,:)=-c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_e2(ii,:,:)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e2(ii,:,:)= c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          elseif(idir==2) then
            c1_e1(:,ii,:)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e1(:,ii,:)=-c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_e2(:,ii,:)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e2(:,ii,:)= c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          elseif(idir==3) then
            c1_e1(:,:,ii)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e1(:,:,ii)=-c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_e2(:,:,ii)= c1_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_e2(:,:,ii)= c2_pml_bot(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          end if
          icount=icount+1
        end do
        
        !h
        if(iend==(fs%lg%is(idir)+fe%ipml_l)) then
          iend=iend-1
        end if
        icount=1
        do ii=fs%mg%is(idir),iend
          if(idir==1) then
            c1_h1(ii,:,:)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h1(ii,:,:)= c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_h2(ii,:,:)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h2(ii,:,:)=-c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          elseif(idir==2) then
            c1_h1(:,ii,:)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h1(:,ii,:)= c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_h2(:,ii,:)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h2(:,ii,:)=-c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          elseif(idir==3) then
            c1_h1(:,:,ii)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h1(:,:,ii)= c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c1_h2(:,:,ii)= c1_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
            c2_h2(:,:,ii)=-c2_pml_bot_h(fs%mg%is(idir)-fs%lg%is(idir)+icount)
          end if
          icount=icount+1
        end do
      end if
      
      !set pml(top)
      if((fs%a_bc(idir,2)=='pml').and.(fs%mg%ie(idir)>=(fs%lg%ie(idir)-fe%ipml_l))) then
        !e
        ista=fs%lg%ie(idir)-fe%ipml_l
        if(fs%mg%is(idir)>ista) then
          ista=fs%mg%is(idir)
        end if
        icount=1
        do ii=ista,fs%mg%ie(idir)
          if(idir==1) then
            c1_e1(ii,:,:)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e1(ii,:,:)=-c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_e2(ii,:,:)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e2(ii,:,:)= c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          elseif(idir==2) then
            c1_e1(:,ii,:)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e1(:,ii,:)=-c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_e2(:,ii,:)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e2(:,ii,:)= c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          elseif(idir==3) then
            c1_e1(:,:,ii)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e1(:,:,ii)=-c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_e2(:,:,ii)= c1_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_e2(:,:,ii)= c2_pml_top((fe%ipml_l+1)-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          end if
          icount=icount+1
        end do
        
        !h
        if(fs%mg%ie(idir)==fs%lg%ie(idir)) then
          iend=fs%mg%ie(idir)-1
        else
          iend=fs%mg%ie(idir)
        end if
        icount=1
        do ii=ista,iend
          if(idir==1) then
            c1_h1(ii,:,:)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h1(ii,:,:)= c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_h2(ii,:,:)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h2(ii,:,:)=-c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          elseif(idir==2) then
            c1_h1(:,ii,:)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h1(:,ii,:)= c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_h2(:,ii,:)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h2(:,ii,:)=-c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          elseif(idir==3) then
            c1_h1(:,:,ii)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h1(:,:,ii)= c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c1_h2(:,:,ii)= c1_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
            c2_h2(:,:,ii)=-c2_pml_top_h(fe%ipml_l-(ista-(fs%lg%ie(idir)-fe%ipml_l)+(icount-1)))
          end if
          icount=icount+1
        end do
      end if
      
      return
    end subroutine eh_set_pml
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check incident current source condition +++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_check_inc(ipulse,ek_dir,epdir_re,epdir_im,inc_dist)
      implicit none
      integer,      intent(in)  :: ipulse
      real(8),      intent(in)  :: ek_dir(3),epdir_re(3),epdir_im(3)
      character(16),intent(out) :: inc_dist
      
      if    (ek_dir(1)==0.0d0.and.ek_dir(2)==0.0d0.and.ek_dir(3)==0.0d0) then 
        inc_dist='none'
      elseif(ek_dir(1)==1.0d0.and.ek_dir(2)==0.0d0.and.ek_dir(3)==0.0d0) then !x-direction propagation
        if(epdir_re(1)/=0.0d0.or.epdir_im(1)/=0.0d0) then
          if(comm_is_root(nproc_id_global)) then
            if     (ipulse==1) then
              write(*,*) "invalid epdir_re1(1) and epdir_im1(1):"
              write(*,*) "For theory = maxwell and ek_dir1(1) = 1.0d0, epdir_re1(1) and epdir_im1(1) must be 0.0d0."
            elseif (ipulse==2) then
              write(*,*) "invalid epdir_re2(1) and epdir_im2(1):"
              write(*,*) "For theory = maxwell and ek_dir2(1) = 1.0d0, epdir_re2(1) and epdir_im2(1) must be 0.0d0."
            end if
          end if
          stop
        else
          inc_dist='yz-plane'
        end if
      elseif(ek_dir(1)==0.0d0.and.ek_dir(2)==1.0d0.and.ek_dir(3)==0.0d0) then !y-direction propagation
        if(epdir_re(2)/=0.0d0.or.epdir_im(2)/=0.0d0) then
          if(comm_is_root(nproc_id_global)) then
            if     (ipulse==1) then
              write(*,*) "invalid epdir_re1(2) and epdir_im1(2):"
              write(*,*) "For theory = maxwell and ek_dir1(2) = 1.0d0, epdir_re1(2) and epdir_im1(2) must be 0.0d0."
            elseif (ipulse==2) then
              write(*,*) "invalid epdir_re2(2) and epdir_im2(2):"
              write(*,*) "For theory = maxwell and ek_dir2(2) = 1.0d0, epdir_re2(2) and epdir_im2(2) must be 0.0d0."
            end if
          end if
          stop
        else
          inc_dist='xz-plane'
        end if
      elseif(ek_dir(1)==0.0d0.and.ek_dir(2)==0.0d0.and.ek_dir(3)==1.0d0) then !z-direction propagation
        if(epdir_re(3)/=0.0d0.or.epdir_im(3)/=0.0d0) then
          if(comm_is_root(nproc_id_global)) then
            if     (ipulse==1) then
              write(*,*) "invalid epdir_re1(3) and epdir_im1(3):"
              write(*,*) "For theory = maxwell and ek_dir1(3) = 1.0d0, epdir_re1(3) and epdir_im1(3) must be 0.0d0."
            elseif (ipulse==2) then
              write(*,*) "invalid epdir_re2(3) and epdir_im2(3):"
              write(*,*) "For theory = maxwell and ek_dir2(3) = 1.0d0, epdir_re2(3) and epdir_im2(3) must be 0.0d0."
            end if
          end if
          stop
        else
          inc_dist='xy-plane'
        end if
      else
        if(comm_is_root(nproc_id_global)) then
          if     (ipulse==1) then
            write(*,*) "invalid ek_dir1:"
            write(*,*) "For theory = maxwell, ek_dir1 is only allowed by"
            write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
          elseif (ipulse==2) then
            write(*,*) "invalid ek_dir2:"
            write(*,*) "For theory = maxwell, ek_dir2 is only allowed by"
            write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
          end if
        end if
        stop
      end if
      
      return
    end subroutine eh_check_inc
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check incident wave parameter +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_check_iw_parameter(ipulse,phi_cep,I_wcm2,E_amplitude,ae_shape)
      implicit none
      integer,      intent(in)  :: ipulse
      real(8),      intent(in)  :: phi_cep
      real(8),      intent(in)  :: I_wcm2
      real(8),      intent(out) :: E_amplitude
      character(16),intent(in)  :: ae_shape
      real(8)                   :: diff_cep
      
      !check ae_shape
      select case(ae_shape)
      case('Ecos2','Acos2')
        continue
      case default
        if(comm_is_root(nproc_id_global)) then
          if     (ipulse==1) then
            write(*,*) 'set ae_shape1 to "Ecos2" or "Acos2".'
          elseif (ipulse==2) then
            write(*,*) 'set ae_shape2 to "Ecos2" or "Acos2".'
          end if
        end if
        stop
      end select
      
      !check phi_cep
      diff_cep=(phi_cep-0.25d0)*2.d0-int((phi_cep-0.25d0)*2.d0)
      if(ae_shape=='Ecos2'.and.abs(diff_cep)>=1.d-12)then
        if(comm_is_root(nproc_id_global)) then
          if     (ipulse==1) then
            write(*,*) 'phi_cep1 must be equal to 0.25+0.5*i when "Ecos2" is specified for ae_shape1.'
          elseif (ipulse==2) then
            write(*,*) 'phi_cep2 must be equal to 0.25+0.5*i when "Ecos2" is specified for ae_shape2.'
          end if
        end if 
        stop
      end if
      
      !set E_amplitude
      if(I_wcm2/=-1d0) &
        E_amplitude=sqrt(I_wcm2)*1.0d2*2.74492d1/(5.14223d11) !I[W/cm^2]->E[a.u.]
        
      return
    end subroutine eh_check_iw_parameter
    
  end subroutine eh_init
  
  !===========================================================================================
  != calculate eh-FDTD =======================================================================
  subroutine eh_calc(fs,fe)
    use salmon_global,   only: dt_em,pole_num_ld,obs_num_em,obs_samp_em,yn_obs_plane_em,yn_obs_plane_integral_em,&
                               base_directory,t1_t2,t1_start,&
                               E_amplitude1,tw1,omega1,phi_cep1,epdir_re1,epdir_im1,ae_shape1,&
                               E_amplitude2,tw2,omega2,phi_cep2,epdir_re2,epdir_im2,ae_shape2
    use inputoutput,     only: utime_from_au,uenergy_from_au
    use parallelization, only: nproc_id_global,nproc_size_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_fdtd_system
    use math_constants,  only: pi
    use ttm,             only: use_ttm,ttm_penetration,ttm_main,ttm_get_temperatures
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    integer                           :: iter,ii,ij,ix,iy,iz
    character(128)                    :: save_name
    !for ttm
    integer :: jx,jy,jz,unit1=4000
    real(8),allocatable :: Spoynting(:,:,:,:), divS(:,:,:)
    real(8),allocatable :: u_energy(:,:,:), u_energy_p(:,:,:)
    real(8),allocatable :: work(:,:,:), work1(:,:,:), work2(:,:,:)
    real(8)             :: dV, tmp(4)
    real(8), parameter  :: hartree_kelvin_relationship = 3.1577502480407d5 ! [K] (+/- 6.1e-07)

    !time-iteration
    do iter=fe%iter_sta,fe%iter_end
      !update iter_now
      fe%iter_now=iter
      if(comm_is_root(nproc_id_global))then
        write(*,*) fe%iter_now
      end if
      
      !update lorentz-drude
      if(fe%num_ld>0) then
        call eh_update_ld
      end if
      
      !update e
      call eh_fd(fe%iex_y_is,fe%iex_y_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ex_y,fe%c2_ex_y,fe%ex_y,fe%hz_x,fe%hz_y,'e','y',cb_y = fe%c_hbloch_y) !ex_y
      call eh_fd(fe%iex_z_is,fe%iex_z_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ex_z,fe%c2_ex_z,fe%ex_z,fe%hy_z,fe%hy_x,'e','z',cb_z = fe%c_hbloch_z) !ex_z
      call eh_fd(fe%iey_z_is,fe%iey_z_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ey_z,fe%c2_ey_z,fe%ey_z,fe%hx_y,fe%hx_z,'e','z',cb_z = fe%c_hbloch_z) !ey_z
      call eh_fd(fe%iey_x_is,fe%iey_x_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ey_x,fe%c2_ey_x,fe%ey_x,fe%hz_x,fe%hz_y,'e','x',cb_x = fe%c_hbloch_x) !ey_x
      call eh_fd(fe%iez_x_is,fe%iez_x_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ez_x,fe%c2_ez_x,fe%ez_x,fe%hy_z,fe%hy_x,'e','x',cb_x = fe%c_hbloch_x) !ez_x
      call eh_fd(fe%iez_y_is,fe%iez_y_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_ez_y,fe%c2_ez_y,fe%ez_y,fe%hx_y,fe%hx_z,'e','y',cb_y = fe%c_hbloch_y) !ez_y
      if(fe%inc_num>0) then !add incident current source
        if(fe%inc_dist1/='none') call eh_add_inc(1,E_amplitude1,tw1,omega1,phi_cep1,&
                                                    epdir_re1,epdir_im1,fe%c2_inc1_xyz,ae_shape1,fe%inc_dist1)
        if(fe%inc_dist2/='none') call eh_add_inc(2,E_amplitude2,tw2,omega2,phi_cep2,&
                                                    epdir_re2,epdir_im2,fe%c2_inc2_xyz,ae_shape2,fe%inc_dist2)
      end if
      if(fe%num_ld>0) then
        call eh_add_curr(fe%rjx_sum_ld(:,:,:),fe%rjy_sum_ld(:,:,:),fe%rjz_sum_ld(:,:,:))
      end if
      call eh_sendrecv(fs,fe,'e')
      
      !calculate linear response
      if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
        call eh_calc_lr
      end if
      
      !store old h
      if( (obs_num_em>0).and.(mod(iter,obs_samp_em)==0) )then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=(fs%mg%is_array(3)),(fs%mg%ie_array(3))
        do iy=(fs%mg%is_array(2)),(fs%mg%ie_array(2))
        do ix=(fs%mg%is_array(1)),(fs%mg%ie_array(1))
          fe%hx_s(ix,iy,iz)=fe%hx_y(ix,iy,iz)+fe%hx_z(ix,iy,iz)
          fe%hy_s(ix,iy,iz)=fe%hy_z(ix,iy,iz)+fe%hy_x(ix,iy,iz)
          fe%hz_s(ix,iy,iz)=fe%hz_x(ix,iy,iz)+fe%hz_y(ix,iy,iz)
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
      
      !update h
      call eh_fd(fe%ihx_y_is,fe%ihx_y_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hx_y,fe%c2_hx_y,fe%hx_y,fe%ez_x,fe%ez_y,'h','y',cb_y = fe%c_ebloch_y) !hx_y
      call eh_fd(fe%ihx_z_is,fe%ihx_z_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hx_z,fe%c2_hx_z,fe%hx_z,fe%ey_z,fe%ey_x,'h','z',cb_z = fe%c_ebloch_z) !hx_z
      call eh_fd(fe%ihy_z_is,fe%ihy_z_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hy_z,fe%c2_hy_z,fe%hy_z,fe%ex_y,fe%ex_z,'h','z',cb_z = fe%c_ebloch_z) !hy_z
      call eh_fd(fe%ihy_x_is,fe%ihy_x_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hy_x,fe%c2_hy_x,fe%hy_x,fe%ez_x,fe%ez_y,'h','x',cb_x = fe%c_ebloch_x) !hy_x
      call eh_fd(fe%ihz_x_is,fe%ihz_x_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hz_x,fe%c2_hz_x,fe%hz_x,fe%ey_z,fe%ey_x,'h','x',cb_x = fe%c_ebloch_x) !hz_x
      call eh_fd(fe%ihz_y_is,fe%ihz_y_ie,      fs%mg%is,fs%mg%ie,fe%Nd,&
                 fe%c1_hz_y,fe%c2_hz_y,fe%hz_y,fe%ex_y,fe%ex_z,'h','y',cb_y = fe%c_ebloch_y) !hz_y
      call eh_sendrecv(fs,fe,'h')
      
      !ttm
      if( use_ttm )then
         !Poynting vector
         if ( .not.allocated(Spoynting) ) then
            call allocate_poynting(fs, Spoynting, divS, u_energy)
            call allocate_poynting(fs, u=u_energy_p)
         end if
         call calc_es_and_hs(fs, fe)
         call calc_poynting_vector(fs, fe, Spoynting)
         call calc_poynting_vector_div(fs, Spoynting, divS)
         u_energy(:,:,:) = u_energy(:,:,:) - divS(:,:,:)*dt_em
         
         u_energy_p = u_energy
         call ttm_penetration( fs%mg%is, u_energy_p )
         
         call ttm_main( fs%srg_ng, fs%mg, u_energy_p )
         
         if( mod(iter,obs_samp_em) == 0 )then
            ix=fs%lg%is(1)-fe%Nd; jx=fs%lg%ie(1)+fe%Nd
            iy=fs%lg%is(2)-fe%Nd; jy=fs%lg%ie(2)+fe%Nd
            iz=fs%lg%is(3)-fe%Nd; jz=fs%lg%ie(3)+fe%Nd
            allocate( work(ix:jx,iy:jy,iz:jz) ); work=0.0d0
            allocate( work1(ix:jx,iy:jy,iz:jz) ); work1=0.0d0
            allocate( work2(ix:jx,iy:jy,iz:jz) ); work2=0.0d0
            call ttm_get_temperatures( (/ix,iy,iz/), work, work2 )
            call comm_summation( work, work1, size(work), nproc_group_global )
            work=work2
            call comm_summation( work, work2, size(work), nproc_group_global )
            if( comm_is_root(nproc_id_global) )then
               unit1=unit1+1
               do iz=lbound(work1,3),ubound(work1,3)
               do ix=lbound(work1,1),ubound(work1,1)
                  write(unit1,'(1x,2i8,2g20.10)') ix,iz,work1(ix,0,iz),work2(ix,0,iz)
               end do
               write(unit1,*)
               end do
            end if
            
            dV=fs%hgs(1)*fs%hgs(2)*fs%hgs(3)
            tmp(3)=sum(u_energy)*dV*uenergy_from_au
            tmp(4)=sum(u_energy_p)*dV*uenergy_from_au
            call comm_summation( tmp(3:4), tmp(1:2), 2, nproc_group_global )
            tmp(3)=sum(work1)/count(work1/=0.0d0)*hartree_kelvin_relationship
            tmp(4)=sum(work2)/count(work2/=0.0d0)*hartree_kelvin_relationship
            if( comm_is_root(nproc_id_global) )then
               open(unit2,file=file_unit2,status='old',position='append')
               write(unit2,"(F16.8,99(1X,E23.15E3))",advance='no') &
                    iter*dt_em*utime_from_au, tmp(1:2), tmp(3:4)
               close(unit2)
            end if
            deallocate( work2 )
            deallocate( work1 )
            deallocate( work )
         end if
      end if !use_ttm

      !observation
      if( (obs_num_em>0).and.(mod(iter,obs_samp_em)==0) )then
        !prepare e and h for save
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=(fs%mg%is_array(3)),(fs%mg%ie_array(3))
        do iy=(fs%mg%is_array(2)),(fs%mg%ie_array(2))
        do ix=(fs%mg%is_array(1)),(fs%mg%ie_array(1))
          fe%ex_s(ix,iy,iz)=fe%ex_y(ix,iy,iz)+fe%ex_z(ix,iy,iz)
          fe%ey_s(ix,iy,iz)=fe%ey_z(ix,iy,iz)+fe%ey_x(ix,iy,iz)
          fe%ez_s(ix,iy,iz)=fe%ez_x(ix,iy,iz)+fe%ez_y(ix,iy,iz)
          fe%hx_s(ix,iy,iz)=( fe%hx_s(ix,iy,iz)+(fe%hx_y(ix,iy,iz)+fe%hx_z(ix,iy,iz)) )/2.0d0
          fe%hy_s(ix,iy,iz)=( fe%hy_s(ix,iy,iz)+(fe%hy_z(ix,iy,iz)+fe%hy_x(ix,iy,iz)) )/2.0d0
          fe%hz_s(ix,iy,iz)=( fe%hz_s(ix,iy,iz)+(fe%hz_x(ix,iy,iz)+fe%hz_y(ix,iy,iz)) )/2.0d0
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        call eh_sendrecv(fs,fe,'s')
        
        !save data
        do ii=1,obs_num_em
          !point
          if(fe%iobs_po_pe(ii)==1) then
            write(save_name,*) ii
            save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//'_at_point_rt.data'
            open(fe%ifn,file=save_name,status='old',position='append')
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no')                                          &
                  dble(iter)*dt_em*utime_from_au,                                                         &
                  fe%ex_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                  fe%ey_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                  fe%ez_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                  fe%hx_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au, &
                  fe%hy_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au, &
                  fe%hz_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au
            close(fe%ifn)
          end if
          
          !plane
          if(yn_obs_plane_em(ii)=='y') then
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ex_s,'ex')
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ey_s,'ey')
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ez_s,'ez')
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hx_s,'hx')
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hy_s,'hy')
            call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                               fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hz_s,'hz')
          end if
          
          !plane integral
          if(yn_obs_plane_integral_em(ii)=='y') then
            call eh_save_plane_integral(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,fe%uAperm_from_au, &
                                        dble(iter)*dt_em*utime_from_au,fs%mg%is,fs%mg%ie,fs%hgs,fe%Nd,fe%ifn,ii,     &
                                        fe%ex_s,fe%ey_s,fe%ez_s,fe%hx_s,fe%hy_s,fe%hz_s)
          end if
          
          !obs_plane_ene_em option(spatial distribution at each energy point)
          if(fe%iobs_num_ene(ii)>0) then
            call eh_calc_plane_ene(fs,fe,ii,iter)
          end if
        end do
        
        !check maximum
        call eh_update_max
      end if
    end do
    
    return
  contains
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ update lorentz-drude ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_update_ld
      implicit none
      
      !initialize
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        fe%rjx_sum_ld(ix,iy,iz)=0.0d0; fe%rjy_sum_ld(ix,iy,iz)=0.0d0; fe%rjz_sum_ld(ix,iy,iz)=0.0d0;
        fe%px_sum_ld(ix,iy,iz) =0.0d0; fe%py_sum_ld(ix,iy,iz) =0.0d0; fe%pz_sum_ld(ix,iy,iz) =0.0d0;
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !update ld polarization vector
      do ii=1,fe%num_ld
      do ij=1,pole_num_ld(fe%media_ld(ii))
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          fe%px_ld(ix,iy,iz,ij,ii)=fe%px_ld(ix,iy,iz,ij,ii)+dt_em*fe%rjx_ld(ix,iy,iz,ij,ii)
          fe%py_ld(ix,iy,iz,ij,ii)=fe%py_ld(ix,iy,iz,ij,ii)+dt_em*fe%rjy_ld(ix,iy,iz,ij,ii)
          fe%pz_ld(ix,iy,iz,ij,ii)=fe%pz_ld(ix,iy,iz,ij,ii)+dt_em*fe%rjz_ld(ix,iy,iz,ij,ii)
          fe%px_sum_ld(ix,iy,iz)=fe%px_sum_ld(ix,iy,iz)+fe%px_ld(ix,iy,iz,ij,ii)
          fe%py_sum_ld(ix,iy,iz)=fe%py_sum_ld(ix,iy,iz)+fe%py_ld(ix,iy,iz,ij,ii)
          fe%pz_sum_ld(ix,iy,iz)=fe%pz_sum_ld(ix,iy,iz)+fe%pz_ld(ix,iy,iz,ij,ii)
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end do
      end do
      
      !update ld polarization  current
      do ii=1,fe%num_ld
      do ij=1,pole_num_ld(fe%media_ld(ii))
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          fe%rjx_ld(ix,iy,iz,ij,ii)= fe%c1_j_ld(ij,ii)*fe%rjx_ld(ix,iy,iz,ij,ii) &
                                    +fe%c2_j_ld(ij,ii)*( fe%ex_y(ix,iy,iz)+fe%ex_z(ix,iy,iz) ) &
                                    *dble(fe%idx_ld(ix,iy,iz,ii)) &
                                    -fe%c3_j_ld(ij,ii)*fe%px_ld(ix,iy,iz,ij,ii)
          fe%rjy_ld(ix,iy,iz,ij,ii)= fe%c1_j_ld(ij,ii)*fe%rjy_ld(ix,iy,iz,ij,ii) &
                                    +fe%c2_j_ld(ij,ii)*( fe%ey_z(ix,iy,iz)+fe%ey_x(ix,iy,iz) ) &
                                    *dble(fe%idy_ld(ix,iy,iz,ii)) &
                                    -fe%c3_j_ld(ij,ii)*fe%py_ld(ix,iy,iz,ij,ii)
          fe%rjz_ld(ix,iy,iz,ij,ii)= fe%c1_j_ld(ij,ii)*fe%rjz_ld(ix,iy,iz,ij,ii) &
                                    +fe%c2_j_ld(ij,ii)*( fe%ez_x(ix,iy,iz)+fe%ez_y(ix,iy,iz) ) &
                                    *dble(fe%idz_ld(ix,iy,iz,ii)) &
                                    -fe%c3_j_ld(ij,ii)*fe%pz_ld(ix,iy,iz,ij,ii)
          fe%rjx_sum_ld(ix,iy,iz)=fe%rjx_sum_ld(ix,iy,iz)+fe%rjx_ld(ix,iy,iz,ij,ii)
          fe%rjy_sum_ld(ix,iy,iz)=fe%rjy_sum_ld(ix,iy,iz)+fe%rjy_ld(ix,iy,iz,ij,ii)
          fe%rjz_sum_ld(ix,iy,iz)=fe%rjz_sum_ld(ix,iy,iz)+fe%rjz_ld(ix,iy,iz,ij,ii)
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end do
      end do
      
      return
    end subroutine eh_update_ld
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ calculate linear response +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_calc_lr
      use salmon_global, only: yn_periodic
      implicit none
      real(8) :: sum_lr_x,sum_lr_y,sum_lr_z
      real(8) :: sum_lr(3),sum_lr2(3)
      
      !update time
      fe%time_lr(fe%iter_lr)=dble(fe%iter_lr)*dt_em
      
      if(yn_periodic=='n') then
        !initialize polarization vector
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          fe%px_lr(ix,iy,iz)=0.0d0; fe%py_lr(ix,iy,iz)=0.0d0; fe%pz_lr(ix,iy,iz)=0.0d0;
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        
        !add all polarization vector
        if(fe%num_ld>0) then
!$omp parallel
!$omp do private(ix,iy,iz)
          do iz=fs%mg%is(3),fs%mg%ie(3)
          do iy=fs%mg%is(2),fs%mg%ie(2)
          do ix=fs%mg%is(1),fs%mg%ie(1)
            fe%px_lr(ix,iy,iz)=fe%px_lr(ix,iy,iz)+fe%px_sum_ld(ix,iy,iz);
            fe%py_lr(ix,iy,iz)=fe%py_lr(ix,iy,iz)+fe%py_sum_ld(ix,iy,iz);
            fe%pz_lr(ix,iy,iz)=fe%pz_lr(ix,iy,iz)+fe%pz_sum_ld(ix,iy,iz);
          end do
          end do
          end do
!$omp end do
!$omp end parallel
        end if
        
        !calculate dipolemoment
        sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          sum_lr_x=sum_lr_x+fe%px_lr(ix,iy,iz)
          sum_lr_y=sum_lr_y+fe%py_lr(ix,iy,iz)
          sum_lr_z=sum_lr_z+fe%pz_lr(ix,iy,iz)
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%dip_lr(fe%iter_lr,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3)
      elseif(yn_periodic=='y') then
        !initialize current density
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          fe%rjx_lr(ix,iy,iz)=0.0d0; fe%rjy_lr(ix,iy,iz)=0.0d0; fe%rjz_lr(ix,iy,iz)=0.0d0;
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        
        !add all current density
        if(fe%num_ld>0) then
!$omp parallel
!$omp do private(ix,iy,iz)
          do iz=fs%mg%is(3),fs%mg%ie(3)
          do iy=fs%mg%is(2),fs%mg%ie(2)
          do ix=fs%mg%is(1),fs%mg%ie(1)
            fe%rjx_lr(ix,iy,iz)=fe%rjx_lr(ix,iy,iz)+fe%rjx_sum_ld(ix,iy,iz)
            fe%rjy_lr(ix,iy,iz)=fe%rjy_lr(ix,iy,iz)+fe%rjy_sum_ld(ix,iy,iz)
            fe%rjz_lr(ix,iy,iz)=fe%rjz_lr(ix,iy,iz)+fe%rjz_sum_ld(ix,iy,iz)
          end do
          end do
          end do
!$omp end do
!$omp end parallel
        end if
        
        !calculate average current density
        sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          sum_lr_x=sum_lr_x+fe%rjx_lr(ix,iy,iz)
          sum_lr_y=sum_lr_y+fe%rjy_lr(ix,iy,iz)
          sum_lr_z=sum_lr_z+fe%rjz_lr(ix,iy,iz)
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%curr_lr(fe%iter_lr,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3) &
                                 /(fs%rlsize(1)*fs%rlsize(2)*fs%rlsize(3))
        
        !calculate average electric field
        sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          sum_lr_x=sum_lr_x+( fe%ex_y(ix,iy,iz)+fe%ex_z(ix,iy,iz) )
          sum_lr_y=sum_lr_y+( fe%ey_z(ix,iy,iz)+fe%ey_x(ix,iy,iz) )
          sum_lr_z=sum_lr_z+( fe%ez_x(ix,iy,iz)+fe%ez_y(ix,iy,iz) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%e_lr(fe%iter_lr,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3) &
                              /(fs%rlsize(1)*fs%rlsize(2)*fs%rlsize(3))
      end if
      
      !update time iteration
      fe%iter_lr=fe%iter_lr+1
      
      return
    end subroutine eh_calc_lr
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ add incident current source +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_add_inc(iord,amp,tw,omega,cep,ep_r,ep_i,c2_inc_xyz,aes,typ)
      implicit none
      integer,intent(in)       :: iord
      real(8),intent(in)       :: amp,tw,omega,cep
      real(8),intent(in)       :: ep_r(3),ep_i(3),c2_inc_xyz(3)
      character(16),intent(in) :: aes,typ
      real(8)                  :: t_sta,t,theta1,theta2_r,theta2_i,alpha,beta,gamma,tf_r,tf_i
      real(8)                  :: add_inc(3)
      
      !calculate time factor and adding current
      if(iord==1) then
        t_sta=t1_start
      elseif(iord==2) then
        t_sta=t1_start+t1_t2
      end if
      t=(dble(iter)-0.5d0)*dt_em-t_sta
      theta1=pi/tw*(t-0.5d0*tw)                         !for cos(theta1)**2
      alpha =pi/tw                                      !for cos(theta1)**2
      theta2_r=omega*(t-0.5d0*tw)+cep*2d0*pi            !for cos(theta2)
      theta2_i=omega*(t-0.5d0*tw)+cep*2d0*pi+3d0/2d0*pi !for cos(theta2), where this is translated to sin.
      beta=omega                                        !for cos(theta2)
      if(t>=0.0d0.and.t<=tw) then
        gamma=1.0d0
      else
        gamma=0.0d0
      end if
      if(aes=='Ecos2')then
        tf_r=cos(theta1)**2*cos(theta2_r)*gamma
        tf_i=cos(theta1)**2*cos(theta2_i)*gamma
      else if(aes=='Acos2')then
        tf_r=-(-alpha*sin(2.d0*theta1)*cos(theta2_r)   &
               -beta*cos(theta1)**2*sin(theta2_r))/beta*gamma
        tf_i=-(-alpha*sin(2.d0*theta1)*cos(theta2_i)   &
               -beta*cos(theta1)**2*sin(theta2_i))/beta*gamma
      else
        tf_r=0.0d0; tf_i=0.0d0;
      end if
!      tf_r=exp(-0.5d0*(( ((dble(iter)-0.5d0)*dt_em-10.0d0*tw1)/tw1 )**2.0d0) ) !test time factor
      add_inc(:)=amp*(tf_r*ep_r(:)+tf_i*ep_i(:))
      
      if(typ=='point_hs') then
        if(fe%inc_po_pe(iord)==1) then
          ix=fe%inc_po_id(iord,1); iy=fe%inc_po_id(iord,2); iz=fe%inc_po_id(iord,3);
          fe%ex_y(ix,iy,iz)=add_inc(1)/2.0d0
          fe%ex_z(ix,iy,iz)=add_inc(1)/2.0d0
          fe%ey_z(ix,iy,iz)=add_inc(2)/2.0d0
          fe%ey_x(ix,iy,iz)=add_inc(2)/2.0d0
          fe%ez_x(ix,iy,iz)=add_inc(3)/2.0d0
          fe%ez_y(ix,iy,iz)=add_inc(3)/2.0d0
        end if
      elseif(typ=='point_ss') then
        if(fe%inc_po_pe(iord)==1) then
          ix=fe%inc_po_id(iord,1); iy=fe%inc_po_id(iord,2); iz=fe%inc_po_id(iord,3);
          fe%ex_y(ix,iy,iz)=fe%ex_y(ix,iy,iz)+add_inc(1)/2.0d0
          fe%ex_z(ix,iy,iz)=fe%ex_z(ix,iy,iz)+add_inc(1)/2.0d0
          fe%ey_z(ix,iy,iz)=fe%ey_z(ix,iy,iz)+add_inc(2)/2.0d0
          fe%ey_x(ix,iy,iz)=fe%ey_x(ix,iy,iz)+add_inc(2)/2.0d0
          fe%ez_x(ix,iy,iz)=fe%ez_x(ix,iy,iz)+add_inc(3)/2.0d0
          fe%ez_y(ix,iy,iz)=fe%ez_y(ix,iy,iz)+add_inc(3)/2.0d0
        end if
      elseif(typ=='x-line_hs') then
        if(fe%inc_li_pe(iord,1)==1) then
          iy=fe%inc_po_id(iord,2); iz=fe%inc_po_id(iord,3);
          fe%ex_y(fe%iex_y_is(1):fe%iex_y_ie(1),iy,iz)=add_inc(1)/2.0d0
          fe%ex_z(fe%iex_z_is(1):fe%iex_z_ie(1),iy,iz)=add_inc(1)/2.0d0
          fe%ey_z(fe%iey_z_is(1):fe%iey_z_ie(1),iy,iz)=add_inc(2)/2.0d0
          fe%ey_x(fe%iey_x_is(1):fe%iey_x_ie(1),iy,iz)=add_inc(2)/2.0d0
          fe%ez_x(fe%iez_x_is(1):fe%iez_x_ie(1),iy,iz)=add_inc(3)/2.0d0
          fe%ez_y(fe%iez_y_is(1):fe%iez_y_ie(1),iy,iz)=add_inc(3)/2.0d0
        end if
      elseif(typ=='x-line_ss') then
        if(fe%inc_li_pe(iord,1)==1) then
          iy=fe%inc_po_id(iord,2); iz=fe%inc_po_id(iord,3);
          fe%ex_y(fe%iex_y_is(1):fe%iex_y_ie(1),iy,iz)=fe%ex_y(fe%iex_y_is(1):fe%iex_y_ie(1),iy,iz)+add_inc(1)/2.0d0
          fe%ex_z(fe%iex_z_is(1):fe%iex_z_ie(1),iy,iz)=fe%ex_z(fe%iex_z_is(1):fe%iex_z_ie(1),iy,iz)+add_inc(1)/2.0d0
          fe%ey_z(fe%iey_z_is(1):fe%iey_z_ie(1),iy,iz)=fe%ey_z(fe%iey_z_is(1):fe%iey_z_ie(1),iy,iz)+add_inc(2)/2.0d0
          fe%ey_x(fe%iey_x_is(1):fe%iey_x_ie(1),iy,iz)=fe%ey_x(fe%iey_x_is(1):fe%iey_x_ie(1),iy,iz)+add_inc(2)/2.0d0
          fe%ez_x(fe%iez_x_is(1):fe%iez_x_ie(1),iy,iz)=fe%ez_x(fe%iez_x_is(1):fe%iez_x_ie(1),iy,iz)+add_inc(3)/2.0d0
          fe%ez_y(fe%iez_y_is(1):fe%iez_y_ie(1),iy,iz)=fe%ez_y(fe%iez_y_is(1):fe%iez_y_ie(1),iy,iz)+add_inc(3)/2.0d0
        end if
      elseif(typ=='y-line_hs') then
        if(fe%inc_li_pe(iord,2)==1) then
          ix=fe%inc_po_id(iord,1); iz=fe%inc_po_id(iord,3);
          fe%ex_y(ix,fe%iex_y_is(2):fe%iex_y_ie(2),iz)=add_inc(1)/2.0d0
          fe%ex_z(ix,fe%iex_z_is(2):fe%iex_z_ie(2),iz)=add_inc(1)/2.0d0
          fe%ey_z(ix,fe%iey_z_is(2):fe%iey_z_ie(2),iz)=add_inc(2)/2.0d0
          fe%ey_x(ix,fe%iey_x_is(2):fe%iey_x_ie(2),iz)=add_inc(2)/2.0d0
          fe%ez_x(ix,fe%iez_x_is(2):fe%iez_x_ie(2),iz)=add_inc(3)/2.0d0
          fe%ez_y(ix,fe%iez_y_is(2):fe%iez_y_ie(2),iz)=add_inc(3)/2.0d0
        end if
      elseif(typ=='y-line_ss') then
        if(fe%inc_li_pe(iord,2)==1) then
          ix=fe%inc_po_id(iord,1); iz=fe%inc_po_id(iord,3);
          fe%ex_y(ix,fe%iex_y_is(2):fe%iex_y_ie(2),iz)=fe%ex_y(ix,fe%iex_y_is(2):fe%iex_y_ie(2),iz)+add_inc(1)/2.0d0
          fe%ex_z(ix,fe%iex_z_is(2):fe%iex_z_ie(2),iz)=fe%ex_z(ix,fe%iex_z_is(2):fe%iex_z_ie(2),iz)+add_inc(1)/2.0d0
          fe%ey_z(ix,fe%iey_z_is(2):fe%iey_z_ie(2),iz)=fe%ey_z(ix,fe%iey_z_is(2):fe%iey_z_ie(2),iz)+add_inc(2)/2.0d0
          fe%ey_x(ix,fe%iey_x_is(2):fe%iey_x_ie(2),iz)=fe%ey_x(ix,fe%iey_x_is(2):fe%iey_x_ie(2),iz)+add_inc(2)/2.0d0
          fe%ez_x(ix,fe%iez_x_is(2):fe%iez_x_ie(2),iz)=fe%ez_x(ix,fe%iez_x_is(2):fe%iez_x_ie(2),iz)+add_inc(3)/2.0d0
          fe%ez_y(ix,fe%iez_y_is(2):fe%iez_y_ie(2),iz)=fe%ez_y(ix,fe%iez_y_is(2):fe%iez_y_ie(2),iz)+add_inc(3)/2.0d0
        end if
      elseif(typ=='z-line_hs') then
        if(fe%inc_li_pe(iord,3)==1) then
          ix=fe%inc_po_id(iord,1); iy=fe%inc_po_id(iord,2);
          fe%ex_y(ix,iy,fe%iex_y_is(3):fe%iex_y_ie(3))=add_inc(1)/2.0d0
          fe%ex_z(ix,iy,fe%iex_z_is(3):fe%iex_z_ie(3))=add_inc(1)/2.0d0
          fe%ey_z(ix,iy,fe%iey_z_is(3):fe%iey_z_ie(3))=add_inc(2)/2.0d0
          fe%ey_x(ix,iy,fe%iey_x_is(3):fe%iey_x_ie(3))=add_inc(2)/2.0d0
          fe%ez_x(ix,iy,fe%iez_x_is(3):fe%iez_x_ie(3))=add_inc(3)/2.0d0
          fe%ez_y(ix,iy,fe%iez_y_is(3):fe%iez_y_ie(3))=add_inc(3)/2.0d0
        end if
      elseif(typ=='z-line_ss') then
        if(fe%inc_li_pe(iord,3)==1) then
          ix=fe%inc_po_id(iord,1); iy=fe%inc_po_id(iord,2);
          fe%ex_y(ix,iy,fe%iex_y_is(3):fe%iex_y_ie(3))=fe%ex_y(ix,iy,fe%iex_y_is(3):fe%iex_y_ie(3))+add_inc(1)/2.0d0
          fe%ex_z(ix,iy,fe%iex_z_is(3):fe%iex_z_ie(3))=fe%ex_z(ix,iy,fe%iex_z_is(3):fe%iex_z_ie(3))+add_inc(1)/2.0d0
          fe%ey_z(ix,iy,fe%iey_z_is(3):fe%iey_z_ie(3))=fe%ey_z(ix,iy,fe%iey_z_is(3):fe%iey_z_ie(3))+add_inc(2)/2.0d0
          fe%ey_x(ix,iy,fe%iey_x_is(3):fe%iey_x_ie(3))=fe%ey_x(ix,iy,fe%iey_x_is(3):fe%iey_x_ie(3))+add_inc(2)/2.0d0
          fe%ez_x(ix,iy,fe%iez_x_is(3):fe%iez_x_ie(3))=fe%ez_x(ix,iy,fe%iez_x_is(3):fe%iez_x_ie(3))+add_inc(3)/2.0d0
          fe%ez_y(ix,iy,fe%iez_y_is(3):fe%iez_y_ie(3))=fe%ez_y(ix,iy,fe%iez_y_is(3):fe%iez_y_ie(3))+add_inc(3)/2.0d0
        end if
      elseif(typ=='xy-plane') then !z propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(fe%inc_pl_pe(iord,1)==1) then
          iz=fe%inc_po_id(iord,3)
!$omp parallel
!$omp do private(ix,iy)
          do iy=fe%iex_z_is(2),fe%iex_z_ie(2)
          do ix=fe%iex_z_is(1),fe%iex_z_ie(1)
            fe%ex_z(ix,iy,iz)=fe%ex_z(ix,iy,iz)+c2_inc_xyz(3)*add_inc(1)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy)
          do iy=fe%iey_z_is(2),fe%iey_z_ie(2)
          do ix=fe%iey_z_is(1),fe%iey_z_ie(1)
            fe%ey_z(ix,iy,iz)=fe%ey_z(ix,iy,iz)+c2_inc_xyz(3)*add_inc(2)
          end do
          end do
!$omp end do
!$omp end parallel
        end if
      elseif(typ=='yz-plane') then !x propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(fe%inc_pl_pe(iord,2)==1) then
          ix=fe%inc_po_id(iord,1)
!$omp parallel
!$omp do private(iy,iz)
          do iz=fe%iey_x_is(3),fe%iey_x_ie(3)
          do iy=fe%iey_x_is(2),fe%iey_x_ie(2)
            fe%ey_x(ix,iy,iz)=fe%ey_x(ix,iy,iz)+c2_inc_xyz(1)*add_inc(2)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz)
          do iz=fe%iez_x_is(3),fe%iez_x_ie(3)
          do iy=fe%iez_x_is(2),fe%iez_x_ie(2)
            fe%ez_x(ix,iy,iz)=fe%ez_x(ix,iy,iz)+c2_inc_xyz(1)*add_inc(3)
          end do
          end do
!$omp end do
!$omp end parallel
        end if
      elseif(typ=='xz-plane') then !y propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if(fe%inc_pl_pe(iord,3)==1) then
          iy=fe%inc_po_id(iord,2)
!$omp parallel
!$omp do private(ix,iz)
          do iz=fe%iex_y_is(3),fe%iex_y_ie(3)
          do ix=fe%iex_y_is(1),fe%iex_y_ie(1)
            fe%ex_y(ix,iy,iz)=fe%ex_y(ix,iy,iz)+c2_inc_xyz(2)*add_inc(1)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz)
          do iz=fe%iez_y_is(3),fe%iez_y_ie(3)
          do ix=fe%iez_y_is(1),fe%iez_y_ie(1)
            fe%ez_y(ix,iy,iz)=fe%ez_y(ix,iy,iz)+c2_inc_xyz(2)*add_inc(3)
          end do
          end do
!$omp end do
!$omp end parallel
        end if
      end if
      
      return
    end subroutine eh_add_inc
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ add current +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_add_curr(rjx,rjy,rjz)
      implicit none
      real(8),intent(in) :: rjx(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                fs%mg%is_array(2):fs%mg%ie_array(2),&
                                fs%mg%is_array(3):fs%mg%ie_array(3)),&
                            rjy(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                fs%mg%is_array(2):fs%mg%ie_array(2),&
                                fs%mg%is_array(3):fs%mg%ie_array(3)),&
                            rjz(fs%mg%is_array(1):fs%mg%ie_array(1),&
                                fs%mg%is_array(2):fs%mg%ie_array(2),&
                                fs%mg%is_array(3):fs%mg%ie_array(3))
      
      !ex
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iex_y_is(3),fe%iex_y_ie(3)
      do iy=fe%iex_y_is(2),fe%iex_y_ie(2)
      do ix=fe%iex_y_is(1),fe%iex_y_ie(1)
        fe%ex_y(ix,iy,iz)=fe%ex_y(ix,iy,iz)+fe%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iex_z_is(3),fe%iex_z_ie(3)
      do iy=fe%iex_z_is(2),fe%iex_z_ie(2)
      do ix=fe%iex_z_is(1),fe%iex_z_ie(1)
        fe%ex_z(ix,iy,iz)=fe%ex_z(ix,iy,iz)+fe%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !ey
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iey_z_is(3),fe%iey_z_ie(3)
      do iy=fe%iey_z_is(2),fe%iey_z_ie(2)
      do ix=fe%iey_z_is(1),fe%iey_z_ie(1)
        fe%ey_z(ix,iy,iz)=fe%ey_z(ix,iy,iz)+fe%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iey_x_is(3),fe%iey_x_ie(3)
      do iy=fe%iey_x_is(2),fe%iey_x_ie(2)
      do ix=fe%iey_x_is(1),fe%iey_x_ie(1)
        fe%ey_x(ix,iy,iz)=fe%ey_x(ix,iy,iz)+fe%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !ez
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iez_x_is(3),fe%iez_x_ie(3)
      do iy=fe%iez_x_is(2),fe%iez_x_ie(2)
      do ix=fe%iez_x_is(1),fe%iez_x_ie(1)
        fe%ez_x(ix,iy,iz)=fe%ez_x(ix,iy,iz)+fe%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fe%iez_y_is(3),fe%iez_y_ie(3)
      do iy=fe%iez_y_is(2),fe%iez_y_ie(2)
      do ix=fe%iez_y_is(1),fe%iez_y_ie(1)
        fe%ez_y(ix,iy,iz)=fe%ez_y(ix,iy,iz)+fe%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      return
    end subroutine eh_add_curr
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check and update maximum of e and h +++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_update_max
      implicit none
      real(8) :: e_abs(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3)),&
                 h_abs(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3))
      real(8) :: e_max_tmp(0:nproc_size_global-1), h_max_tmp(0:nproc_size_global-1),&
                 e_max_tmp2(0:nproc_size_global-1),h_max_tmp2(0:nproc_size_global-1)
      
      e_max_tmp(:)=0.0d0; h_max_tmp(:)=0.0d0;
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        e_abs(ix,iy,iz)=sqrt( fe%ex_s(ix,iy,iz)**2.0d0 + fe%ey_s(ix,iy,iz)**2.0d0 + fe%ez_s(ix,iy,iz)**2.0d0 )
        h_abs(ix,iy,iz)=sqrt( fe%hx_s(ix,iy,iz)**2.0d0 + fe%hy_s(ix,iy,iz)**2.0d0 + fe%hz_s(ix,iy,iz)**2.0d0 )
        if(e_max_tmp(nproc_id_global)<e_abs(ix,iy,iz)) e_max_tmp(nproc_id_global)=e_abs(ix,iy,iz)
        if(h_max_tmp(nproc_id_global)<h_abs(ix,iy,iz)) h_max_tmp(nproc_id_global)=h_abs(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(e_max_tmp,e_max_tmp2,nproc_size_global,nproc_group_global)
      call comm_summation(h_max_tmp,h_max_tmp2,nproc_size_global,nproc_group_global)
      e_max_tmp2(:)=e_max_tmp2(:)*fe%uVperm_from_au
      h_max_tmp2(:)=h_max_tmp2(:)*fe%uAperm_from_au
      if(fe%e_max<maxval(e_max_tmp2(:))) fe%e_max=maxval(e_max_tmp2(:))
      if(fe%h_max<maxval(h_max_tmp2(:))) fe%h_max=maxval(h_max_tmp2(:))
      
      return
    end subroutine eh_update_max
    
  end subroutine eh_calc
  
  !===========================================================================================
  != finalize eh-FDTD ========================================================================
  subroutine eh_finalize(fs,fe)
    use salmon_global,   only: dt_em,unit_system,yn_periodic,ae_shape1,ae_shape2,e_impulse,sysname, &
                               nt_em,nenergy,de,base_directory,obs_num_em,obs_samp_em,obs_plane_ene_em,yn_obs_plane_em
    use inputoutput,     only: utime_from_au,ulength_from_au,uenergy_from_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_fdtd_system
    use math_constants,  only: pi,zi
    implicit none
    type(s_fdtd_system),intent(in)    :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    complex(8),allocatable            :: z_ex1(:,:),z_ey1(:,:),z_ez1(:,:),z_hx1(:,:),z_hy1(:,:),z_hz1(:,:),&
                                         z_ex2(:,:),z_ey2(:,:),z_ez2(:,:),z_hx2(:,:),z_hy2(:,:),z_hz2(:,:)
    integer                           :: ii,ij,ik,il,i1,i1s,i2,i2s
    complex(8)                        :: z_tmp(3)
    character(2)                      :: plane_name
    character(128)                    :: iobs_name,iene_name,wf_name,save_name
    
    !output linear response(matter dipole pm and current jm are outputted: pm = -dip and jm = -curr)
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      if(yn_periodic=='n') then
        !output time-dependent dipole data
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Real time calculation:" 
          write(fe%ifn,'(A)') "# ddm_e: Change of dipole moment (electrons/plus definition)" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[a.u.]",                &
                  2, "ddm_e_x[a.u.]",             &
                  3, "ddm_e_y[a.u.]",             &
                  4, "ddm_e_z[a.u.]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[fs]",                  &
                  2, "ddm_e_x[Angstrom]",         &
                  3, "ddm_e_y[Angstrom]",         &
                  4, "ddm_e_z[Angstrom]"
          end select
          do ii=1,nt_em
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
                  fe%time_lr(ii)*utime_from_au,                  &
                  fe%dip_lr(ii,1:3)*ulength_from_au
            write(fe%ifn,*)
          end do
          close(fe%ifn)
        end if
        
        !output response data
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,1),fe%fr_lr(:,1),fe%fi_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,2),fe%fr_lr(:,2),fe%fi_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,3),fe%fr_lr(:,3),fe%fi_lr(:,3))
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Fourier-transform spectra:" 
          write(fe%ifn,'(A)') "# alpha: Polarizability" 
          write(fe%ifn,'(A)') "# df/dE: Strength function" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Energy[a.u.]",              &
                  2, "Re(alpha_x)[a.u.]",         &
                  3, "Im(alpha_x)[a.u.]",         &
                  4, "Re(alpha_y)[a.u.]",         &
                  5, "Im(alpha_y)[a.u.]",         &
                  6, "Re(alpha_z)[a.u.]",         &
                  7, "Im(alpha_z)[a.u.]",         &
                  8, "df_x/dE[none]",             &
                  9, "df_y/dE[none]",             &
                  10,"df_z/dE[none]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Energy[eV]",                &
                  2, "Re(alpha_x)[Augstrom^2/V]", &
                  3, "Im(alpha_x)[Augstrom^2/V]", &
                  4, "Re(alpha_y)[Augstrom^2/V]", &
                  5, "Im(alpha_y)[Augstrom^2/V]", &
                  6, "Re(alpha_z)[Augstrom^2/V]", &
                  7, "Im(alpha_z)[Augstrom^2/V]", &
                  8, "df_x/dE[none]",             &
                  9, "df_y/dE[none]",             &
                  10,"df_z/dE[none]"
          end select
          do ii=1,nenergy
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no')                      &
                  dble(ii)*de*uenergy_from_au,                                        &
                  fe%fr_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fr_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fr_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  2.0d0*dble(ii)*de/pi*fe%fi_lr(ii,1:3)/(-e_impulse)
            write(fe%ifn,*)
          end do
          close(fe%ifn)
        end if
      elseif(yn_periodic=='y') then
        !output time-dependent average matter current density data and average electric field data
        if(comm_is_root(nproc_id_global)) then
          !average matter current density data
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Real time calculation:" 
          write(fe%ifn,'(A)') "# Jm: Matter current density (electrons)" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[a.u.]",                &
                  2, "Jm_x[a.u.]",                &
                  3, "Jm_y[a.u.]",                &
                  4, "Jm_z[a.u.]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[fs]",                  &
                  2, "Jm_x[1/fs*Angstrom^2]",     &
                  3, "Jm_y[1/fs*Angstrom^2]",     &
                  4, "Jm_z[1/fs*Angstrom^2]"
          end select
          do ii=1,nt_em
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
                  fe%time_lr(ii)*utime_from_au,                  &
                 -fe%curr_lr(ii,:)*((ulength_from_au/utime_from_au)/ulength_from_au**3)
            write(fe%ifn,*)
          end do
          close(fe%ifn)
          
          !average electric field data
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt_e.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Real time calculation:" 
          write(fe%ifn,'(A)') "# E: Averaged electric field in the unit cell" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[a.u.]",                &
                  2, "E_x[a.u.]",                 &
                  3, "E_y[a.u.]",                 &
                  4, "E_z[a.u.]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Time[fs]",                  &
                  2, "E_x[V/Angstrom]",           &
                  3, "E_y[V/Angstrom]",           &
                  4, "E_z[V/Angstrom]"
          end select
          do ii=1,nt_em
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
                 (fe%time_lr(ii)+0.5d0*dt_em)*utime_from_au,     &
                  fe%e_lr(ii,:)*fe%uVperm_from_au
            write(fe%ifn,*)
          end do
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! 0.5d0*dt_em is introduced to adjust actual time of electric field !
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          close(fe%ifn)
        end if
        
        !output response data
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,1),fe%fr_lr(:,1),fe%fi_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,2),fe%fr_lr(:,2),fe%fi_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,3),fe%fr_lr(:,3),fe%fi_lr(:,3))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr+0.5d0*dt_em,fe%e_lr(:,1),fe%er_lr(:,1),fe%ei_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr+0.5d0*dt_em,fe%e_lr(:,2),fe%er_lr(:,2),fe%ei_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr+0.5d0*dt_em,fe%e_lr(:,3),fe%er_lr(:,3),fe%ei_lr(:,3))
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
          open(fe%ifn,file=save_name)
          write(fe%ifn,'(A)') "# Fourier-transform spectra:" 
          write(fe%ifn,'(A)') "# sigma: Conductivity" 
          write(fe%ifn,'(A)') "# eps: Dielectric constant" 
          select case(unit_system)
          case('au','a.u.')
            write(fe%ifn,'("#",99(1X,I0,":",A))') &
                  1, "Energy[a.u.]",              &
                  2, "Re(sigma_x)[a.u.]",         &
                  3, "Im(sigma_x)[a.u.]",         &
                  4, "Re(sigma_y)[a.u.]",         &
                  5, "Im(sigma_y)[a.u.]",         &
                  6, "Re(sigma_z)[a.u.]",         &
                  7, "Im(sigma_z)[a.u.]",         &
                  8, "Re(eps_x)[none]",           &
                  9, "Im(eps_x)[none]",           &
                  10,"Re(eps_y)[none]",           &
                  11,"Im(eps_y)[none]",           &
                  12,"Re(eps_z)[none]",           &
                  13,"Im(eps_z)[none]"
          case('A_eV_fs')
            write(fe%ifn,'("#",99(1X,I0,":",A))')    &
                  1, "Energy[eV]",                   &
                  2, "Re(sigma_x)[1/fs*V*Angstrom]", &
                  3, "Im(sigma_x)[1/fs*V*Angstrom]", &
                  4, "Re(sigma_y)[1/fs*V*Angstrom]", &
                  5, "Im(sigma_y)[1/fs*V*Angstrom]", &
                  6, "Re(sigma_z)[1/fs*V*Angstrom]", &
                  7, "Im(sigma_z)[1/fs*V*Angstrom]", &
                  8, "Re(eps_x)[none]",              &
                  9, "Im(eps_x)[none]",              &
                  10,"Re(eps_y)[none]",              &
                  11,"Im(eps_y)[none]",              &
                  12,"Re(eps_z)[none]",              &
                  13,"Im(eps_z)[none]"
          end select
          do ii=1,nenergy
            z_tmp(:)=(fe%fr_lr(ii,:)+zi*fe%fi_lr(ii,:)) / (-e_impulse+fe%er_lr(ii,:)+zi*fe%ei_lr(ii,:))
            write(fe%ifn,"(F16.8,99(1X,E23.15E3))",advance='no')                           &
                  dble(ii)*de*uenergy_from_au,                                             &
                  real( z_tmp(1))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  aimag(z_tmp(1))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  real( z_tmp(2))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  aimag(z_tmp(2))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  real( z_tmp(3))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  aimag(z_tmp(3))*(1.0d0/utime_from_au/fe%uVperm_from_au/ulength_from_au), &
                  1.0d0-4.0d0*pi*aimag(z_tmp(1))/(dble(ii)*de),                            &
                  4.0d0*pi*real(z_tmp(1))/(dble(ii)*de),                                   &
                  1.0d0-4.0d0*pi*aimag(z_tmp(2))/(dble(ii)*de),                            &
                  4.0d0*pi*real(z_tmp(2))/(dble(ii)*de),                                   &
                  1.0d0-4.0d0*pi*aimag(z_tmp(3))/(dble(ii)*de),                            &
                  4.0d0*pi*real(z_tmp(3))/(dble(ii)*de)
            write(fe%ifn,*)
          end do
        end if
      end if
    end if
    
    !output for obs_plane_ene_em option
    if(obs_num_em>0) then
      do ii=1,obs_num_em !observation point---------------------------------------------------------------------------!
        if(fe%iobs_num_ene(ii)>0) then
          do ij=1,fe%iobs_num_ene(ii) !energy-------------------------------------------------------------------------!
            do ik=1,3 !plane------------------------------------------------------------------------------------------!
              !set plane
              if(ik==1)     then !xy plane
                i1s=1; i2s=2; plane_name='xy';
              elseif(ik==2) then !yz plane
                i1s=2; i2s=3; plane_name='yz';
              elseif(ik==3) then !xz plane
                i1s=1; i2s=3; plane_name='xz';
              end if
              do il=1,2 !with or without window function--------------------------------------------------------------!
                !allocation
                allocate(z_ex1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_ex2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_ey1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_ey2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_ez1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_ez2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hx1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hx2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hy1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hy2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hz1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_hz2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)) )
                z_ex1(:,:)=(0.0d0,0.0d0); z_ey1(:,:)=(0.0d0,0.0d0); z_ez1(:,:)=(0.0d0,0.0d0);
                z_ex2(:,:)=(0.0d0,0.0d0); z_ey2(:,:)=(0.0d0,0.0d0); z_ez2(:,:)=(0.0d0,0.0d0);
                z_hx1(:,:)=(0.0d0,0.0d0); z_hy1(:,:)=(0.0d0,0.0d0); z_hz1(:,:)=(0.0d0,0.0d0);
                z_hx2(:,:)=(0.0d0,0.0d0); z_hy2(:,:)=(0.0d0,0.0d0); z_hz2(:,:)=(0.0d0,0.0d0);
                         
                !collect lg-data
                do i2=fs%mg%is(i2s),fs%mg%ie(i2s)
                do i1=fs%mg%is(i1s),fs%mg%ie(i1s)
                  if(ik==1)     then !xy plane
                    z_ex1(i1,i2)=fe%obs_ex_xy_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2)=fe%obs_ey_xy_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2)=fe%obs_ez_xy_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2)=fe%obs_hx_xy_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2)=fe%obs_hy_xy_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2)=fe%obs_hz_xy_ene(i1,i2,ii,ij,il)
                  elseif(ik==2) then !yz plane
                    z_ex1(i1,i2)=fe%obs_ex_yz_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2)=fe%obs_ey_yz_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2)=fe%obs_ez_yz_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2)=fe%obs_hx_yz_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2)=fe%obs_hy_yz_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2)=fe%obs_hz_yz_ene(i1,i2,ii,ij,il)
                  elseif(ik==3) then !xz plane
                    z_ex1(i1,i2)=fe%obs_ex_xz_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2)=fe%obs_ey_xz_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2)=fe%obs_ez_xz_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2)=fe%obs_hx_xz_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2)=fe%obs_hy_xz_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2)=fe%obs_hz_xz_ene(i1,i2,ii,ij,il)
                  end if
                end do
                end do
                call comm_summation(z_ex1,z_ex2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_ey1,z_ey2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_ez1,z_ez2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hx1,z_hx2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hy1,z_hy2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hz1,z_hz2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                
                !output lg-data
                if(comm_is_root(nproc_id_global)) then
                  write(iobs_name,*) ii
                  write(iene_name,*) ij
                  if(il==1)     then
                    wf_name='with_wf'
                  elseif(il==2) then
                    wf_name='without_wf'
                  end if
                  save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(iobs_name))//&
                            '_ene'//trim(adjustl(iene_name))//'_'//plane_name//'_'//trim(wf_name)//'.data'
                  open(fe%ifn,file=save_name)
                  if(il==1)     then
                    write(fe%ifn,'(A)') "# Fourier-transformed spatial distribution with window function:"
                  elseif(il==2) then
                    write(fe%ifn,'(A)') "# Fourier-transformed spatial distribution without window function:"
                  end if
                  if(ik==1)     then !xy plane
                    write(fe%ifn,'(A)') "# ID_1: Grid ID along x axis"
                    write(fe%ifn,'(A)') "# ID_2: Grid ID along y axis"
                  elseif(ik==2) then !yz plane
                    write(fe%ifn,'(A)') "# ID_1: Grid ID along y axis"
                    write(fe%ifn,'(A)') "# ID_2: Grid ID along z axis"
                  elseif(ik==3) then !xz plane
                    write(fe%ifn,'(A)') "# ID_1: Grid ID along x axis"
                    write(fe%ifn,'(A)') "# ID_2: Grid ID along z axis"
                  end if
                  write(fe%ifn,'(A)') "# E: Electric field"
                  write(fe%ifn,'(A)') "# H: Magnetic field"
                  select case(unit_system)
                  case('au','a.u.')
                    write(fe%ifn,'(A,E23.15E3,A)') "# Sampling energy: ",obs_plane_ene_em(ii,ij),' a.u.'
                    write(fe%ifn,'("#",99(1X,I0,":",A))') &
                          1, "ID_1[none]",                &
                          2, "ID_2[none]",                &
                          3, "Re(E_x)[a.u.]",             &
                          4, "Im(E_x)[a.u.]",             &
                          5, "Re(E_y)[a.u.]",             &
                          6, "Im(E_y)[a.u.]",             &
                          7, "Re(E_z)[a.u.]",             &
                          8, "Im(E_z)[a.u.]",             &
                          9, "Re(H_x)[a.u.]",             &
                          10,"Im(H_x)[a.u.]",             &
                          11,"Re(H_y)[a.u.]",             &
                          12,"Im(H_y)[a.u.]",             &
                          13,"Re(H_z)[a.u.]",             &
                          14,"Im(H_z)[a.u.]"
                  case('A_eV_fs')
                    write(fe%ifn,'(A,E23.15E3,A)') "# Sampling energy: ",obs_plane_ene_em(ii,ij)*uenergy_from_au,' eV'
                    write(fe%ifn,'("#",99(1X,I0,":",A))') &
                          1, "ID_1[none]",                &
                          2, "ID_2[none]",                &
                          3, "Re(E_x)[V/Angstrom*fs]",    &
                          4, "Im(E_x)[V/Angstrom*fs]",    &
                          5, "Re(E_y)[V/Angstrom*fs]",    &
                          6, "Im(E_y)[V/Angstrom*fs]",    &
                          7, "Re(E_z)[V/Angstrom*fs]",    &
                          8, "Im(E_z)[V/Angstrom*fs]",    &
                          9, "Re(H_x)[A/Angstrom*fs]",    &
                          10,"Im(H_x)[A/Angstrom*fs]",    &
                          11,"Re(H_y)[A/Angstrom*fs]",    &
                          12,"Im(H_y)[A/Angstrom*fs]",    &
                          13,"Re(H_z)[A/Angstrom*fs]",    &
                          14,"Im(H_z)[A/Angstrom*fs]"
                  end select
                  do i2=fs%lg%is(i2s),fs%lg%ie(i2s)
                  do i1=fs%lg%is(i1s),fs%lg%ie(i1s)
                    write(fe%ifn,'(I8,I8,99(1X,E23.15E3))')                   &
                          i1,i2,                                              &
                          real( z_ex2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          aimag(z_ex2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          real( z_ey2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          aimag(z_ey2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          real( z_ez2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          aimag(z_ez2(i1,i2))*fe%uVperm_from_au*utime_from_au,&
                          real( z_hx2(i1,i2))*fe%uAperm_from_au*utime_from_au,&
                          aimag(z_hx2(i1,i2))*fe%uAperm_from_au*utime_from_au,&
                          real( z_hy2(i1,i2))*fe%uAperm_from_au*utime_from_au,&
                          aimag(z_hy2(i1,i2))*fe%uAperm_from_au*utime_from_au,&
                          real( z_hz2(i1,i2))*fe%uAperm_from_au*utime_from_au,&
                          aimag(z_hz2(i1,i2))*fe%uAperm_from_au*utime_from_au
                  end do
                  end do
                  close(fe%ifn)
                end if
                
                !deallocation
                deallocate(z_ex1,z_ey1,z_ez1,z_hx1,z_hy1,z_hz1,&
                           z_ex2,z_ey2,z_ez2,z_hx2,z_hy2,z_hz2)
              end do !with or without window function-----------------------------------------------------------------!
            end do !plane---------------------------------------------------------------------------------------------!
          end do !energy----------------------------------------------------------------------------------------------!
        end if
      end do !observation point---------------------------------------------------------------------------------------!
    end if
    
    !make information file for observation
    if(obs_num_em>0) then
      if(comm_is_root(nproc_id_global)) then
        open(fe%ifn,file=trim(base_directory)//"/obs0_info.data")
        write(fe%ifn,'(A,A23)')         'unit_system          =',trim(unit_system)
        write(fe%ifn,'(A,A23)')         'yn_periodic          =',yn_periodic
        write(fe%ifn,'(A,E23.15E3)')    'dt_em                =',dt_em*utime_from_au
        write(fe%ifn,'(A,I23)')         'nt_em                =',(fe%iter_end-fe%iter_sta+1)
        write(fe%ifn,'(3(A,E23.15E3))') 'al_em                =',fs%rlsize(1)*ulength_from_au,', ',&
                                                                 fs%rlsize(2)*ulength_from_au,', ',&
                                                                 fs%rlsize(3)*ulength_from_au
        write(fe%ifn,'(3(A,E23.15E3))') 'dl_em                =',fs%hgs(1)*ulength_from_au,', ',&
                                                                 fs%hgs(2)*ulength_from_au,', ',&
                                                                 fs%hgs(3)*ulength_from_au
        write(fe%ifn,'(3(A,I23))')      'lg_sta               =',fs%lg%is(1),', ',fs%lg%is(2),', ',fs%lg%is(3)
        write(fe%ifn,'(3(A,I23))')      'lg_end               =',fs%lg%ie(1),', ',fs%lg%ie(2),', ',fs%lg%ie(3)
        write(fe%ifn,'(A,I23)')         'obs_num_em           =',obs_num_em
        write(fe%ifn,'(A,I23)')         'obs_samp_em          =',obs_samp_em
        do ii=1,obs_num_em
          write(fe%ifn,'(A,I3,A,A23)')  'yn_obs_plane_em(',&
                                                        ii,') =',yn_obs_plane_em(ii)
        end do
        write(fe%ifn,'(A,E23.15E3)')    'e_max                =',fe%e_max
        write(fe%ifn,'(A,E23.15E3)')    'h_max                =',fe%h_max
        close(fe%ifn)
      end if
    end if
    
    !deallocate
    deallocate(fe%ex_y,fe%c1_ex_y,fe%c2_ex_y,fe%ex_z,fe%c1_ex_z,fe%c2_ex_z,&
               fe%ey_z,fe%c1_ey_z,fe%c2_ey_z,fe%ey_x,fe%c1_ey_x,fe%c2_ey_x,&
               fe%ez_x,fe%c1_ez_x,fe%c2_ez_x,fe%ez_y,fe%c1_ez_y,fe%c2_ez_y,&
               fe%hx_y,fe%c1_hx_y,fe%c2_hx_y,fe%hx_z,fe%c1_hx_z,fe%c2_hx_z,&
               fe%hy_z,fe%c1_hy_z,fe%c2_hy_z,fe%hy_x,fe%c1_hy_x,fe%c2_hy_x,&
               fe%hz_x,fe%c1_hz_x,fe%c2_hz_x,fe%hz_y,fe%c1_hz_y,fe%c2_hz_y)
    
    !write end
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "-------------------------------------------------------"
      write(*,*) "**************************"
      write(*,*) "FDTD end"
      write(*,*) "**************************"
    end if
    
    return
  end subroutine eh_finalize
  
  !===========================================================================================
  != prepare mpi, grid, and sendrecv enviroments==============================================
  subroutine eh_mpi_grid_sr(fs,fe)
    use salmon_global,     only: nproc_rgrid,nproc_k,nproc_ob
    use parallelization,   only: nproc_group_global
    use set_numcpu,        only: set_numcpu_general,iprefer_domain_distribution
    use init_communicator, only: init_communicator_dft
    use sendrecv_grid,     only: create_sendrecv_neig,init_sendrecv_grid
    use structures,        only: s_fdtd_system, s_parallel_info
    use initialization_sub
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    type(s_parallel_info)             :: info
    integer                           :: neig_ng_eh(1:2,1:3)
    integer                           :: ii
    
    !set mpi condition
    if((nproc_k==0).and.(nproc_ob==0).and.(sum(nproc_rgrid(:))==0)) then
      call set_numcpu_general(iprefer_domain_distribution,1,1,nproc_group_global,info)
    else
      info%npk       = nproc_k
      info%nporbital = nproc_ob
      info%nprgrid   = nproc_rgrid
    end if
    call init_communicator_dft(nproc_group_global,info)
    
    !initialize r-grid
    call init_grid_whole(fs%rlsize,fs%hgs,fs%lg)
    call init_grid_parallel(info,fs%lg,fs%mg) ! lg --> mg
    !### This process about ng is temporal. #####################!
    !### With modifying set_ng to be applied to arbitrary Nd, ###!
    !### this process will be removed.###########################!
    fs%mg%is_overlap(1:3)=fs%mg%is(1:3)-fe%Nd
    fs%mg%ie_overlap(1:3)=fs%mg%ie(1:3)+fe%Nd
    fs%mg%is_array(1:3)  =fs%mg%is(1:3)-fe%Nd
    fs%mg%ie_array(1:3)  =fs%mg%ie(1:3)+fe%Nd
    !############################################################!
    
    !prepare for setting sendrecv environment
    if(allocated(fs%mg%idx)) deallocate(fs%mg%idx)
    if(allocated(fs%mg%idy)) deallocate(fs%mg%idy)
    if(allocated(fs%mg%idz)) deallocate(fs%mg%idz)
    allocate(fs%mg%idx(fs%mg%is_overlap(1):fs%mg%ie_overlap(1)), &
             fs%mg%idy(fs%mg%is_overlap(2):fs%mg%ie_overlap(2)), &
             fs%mg%idz(fs%mg%is_overlap(3):fs%mg%ie_overlap(3)))
    do ii=fs%mg%is_overlap(1),fs%mg%ie_overlap(1)
      fs%mg%idx(ii)=ii
    end do
    do ii=fs%mg%is_overlap(2),fs%mg%ie_overlap(2)
      fs%mg%idy(ii)=ii
    end do
    do ii=fs%mg%is_overlap(3),fs%mg%ie_overlap(3)
      fs%mg%idz(ii)=ii
    end do
    fs%mg%Nd=fe%Nd
    
    !set sendrecv environment
    call create_sendrecv_neig(neig_ng_eh,info) ! neighboring node array
    call init_sendrecv_grid(fs%srg_ng,fs%mg,1,info%icomm_r,neig_ng_eh)
    
    return
  end subroutine eh_mpi_grid_sr
  
  !===========================================================================================
  != send and receive eh =====================================================================
  subroutine eh_sendrecv(fs,fe,var)
    use sendrecv_grid,  only: update_overlap_real8
    use structures,     only: s_fdtd_system
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh), intent(inout)   :: fe
    character(1),intent(in)           :: var
    integer                           :: ix,iy,iz
    real(8),allocatable               :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
    
    if(var=='e') then
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ex_y)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ex_z)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ey_z)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ey_x)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ez_x)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%ez_y)
    elseif(var=='h') then
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hx_y)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hx_z)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hy_z)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hy_x)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hz_x)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hz_y)
    elseif(var=='r') then
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%rmedia)
    elseif(var=='s') then
      !allocate temporary variable
      allocate(f1(fs%mg%is_array(1):fs%mg%ie_array(1),&
                  fs%mg%is_array(2):fs%mg%ie_array(2),&
                  fs%mg%is_array(3):fs%mg%ie_array(3)),&
               f2(fs%mg%is_array(1):fs%mg%ie_array(1),&
                  fs%mg%is_array(2):fs%mg%ie_array(2),&
                  fs%mg%is_array(3):fs%mg%ie_array(3)),&
               f3(fs%mg%is_array(1):fs%mg%ie_array(1),&
                  fs%mg%is_array(2):fs%mg%ie_array(2),&
                  fs%mg%is_array(3):fs%mg%ie_array(3)))
      f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
      
      !spatially adjust e for save
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        f1(ix,iy,iz)=( fe%ex_s(ix,iy,iz)+fe%ex_s(ix-1,iy,iz) )/2.0d0
        f2(ix,iy,iz)=( fe%ey_s(ix,iy,iz)+fe%ey_s(ix,iy-1,iz) )/2.0d0
        f3(ix,iy,iz)=( fe%ez_s(ix,iy,iz)+fe%ez_s(ix,iy,iz-1) )/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      fe%ex_s(:,:,:)=f1(:,:,:); fe%ey_s(:,:,:)=f2(:,:,:); fe%ez_s(:,:,:)=f3(:,:,:);
      f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
      
      !spatially adjust h for save
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        f1(ix,iy,iz)=( fe%hx_s(ix,iy,iz)+fe%hx_s(ix,iy-1,iz) )/2.0d0
        f2(ix,iy,iz)=( fe%hy_s(ix,iy,iz)+fe%hy_s(ix,iy,iz-1) )/2.0d0
        f3(ix,iy,iz)=( fe%hz_s(ix,iy,iz)+fe%hz_s(ix-1,iy,iz) )/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      fe%hx_s(:,:,:)=f1(:,:,:); fe%hy_s(:,:,:)=f2(:,:,:); fe%hz_s(:,:,:)=f3(:,:,:);
      f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hx_s)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hy_s)
      call update_overlap_real8(fs%srg_ng,fs%mg,fe%hz_s)
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        f1(ix,iy,iz)=( fe%hx_s(ix,iy,iz)+fe%hx_s(ix,iy,iz-1) )/2.0d0
        f2(ix,iy,iz)=( fe%hy_s(ix,iy,iz)+fe%hy_s(ix-1,iy,iz) )/2.0d0
        f3(ix,iy,iz)=( fe%hz_s(ix,iy,iz)+fe%hz_s(ix,iy-1,iz) )/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      fe%hx_s(:,:,:)=f1(:,:,:); fe%hy_s(:,:,:)=f2(:,:,:); fe%hz_s(:,:,:)=f3(:,:,:);
      
      !deallocate temporary variable
      deallocate(f1,f2,f3)
    end if
    
    return
  end subroutine eh_sendrecv
  
  !===========================================================================================
  != calculate finite difference in eh =======================================================
  subroutine eh_fd(ista,iend,ng_is,ng_ie,Nd,c1,c2,f1,f2,f3,var,dir,&
                   cb_x,cb_y,cb_z)
    implicit none
    integer,     intent(in)          :: ista(3),iend(3),ng_is(3),ng_ie(3)
    integer,     intent(in)          :: Nd
    real(8),     intent(in)          :: c1(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd),&
                                        c2(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    real(8),     intent(inout)       :: f1(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    real(8),     intent(in)          :: f2(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd),&
                                        f3(ng_is(1)-Nd:ng_ie(1)+Nd, &
                                           ng_is(2)-Nd:ng_ie(2)+Nd, &
                                           ng_is(3)-Nd:ng_ie(3)+Nd)
    character(1),intent(in)          :: var,dir
    real(8),     intent(in),optional :: cb_x(ng_is(1)-Nd:ng_ie(1)+Nd), &
                                        cb_y(ng_is(2)-Nd:ng_ie(2)+Nd), &
                                        cb_z(ng_is(3)-Nd:ng_ie(3)+Nd)
    integer :: ix,iy,iz
    
    if(var=='e') then
      if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_x(ix-1)*(f2(ix-1,iy,iz)+f3(ix-1,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_y(iy-1)*(f2(ix,iy-1,iz)+f3(ix,iy-1,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       (f2(ix,iy,iz)+f3(ix,iy,iz)) - cb_z(iz-1)*(f2(ix,iy,iz-1)+f3(ix,iy,iz-1)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(var=='h') then
      if(dir=='x') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_x(ix+1)*(f2(ix+1,iy,iz)+f3(ix+1,iy,iz)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='y') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_y(iy+1)*(f2(ix,iy+1,iz)+f3(ix,iy+1,iz)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(dir=='z') then
!$omp parallel
!$omp do private(ix,iy,iz)
        do iz=ista(3),iend(3)
        do iy=ista(2),iend(2)
        do ix=ista(1),iend(1)
          f1(ix,iy,iz)= c1(ix,iy,iz)*f1(ix,iy,iz) &
                       +c2(ix,iy,iz)*( &
                       cb_z(iz+1)*(f2(ix,iy,iz+1)+f3(ix,iy,iz+1)) - (f2(ix,iy,iz)+f3(ix,iy,iz)) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    end if
    
    return
  end subroutine eh_fd
  
  !===========================================================================================
  != save plane data =========================================================================
  subroutine eh_save_plane(id,ipl,conv,ng_is,ng_ie,lg_is,lg_ie,Nd,ifn,iobs,iter,f,var)
    use salmon_global,   only: base_directory
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    implicit none
    integer,intent(in)      :: id(3),ipl(3)
    real(8),intent(in)      :: conv
    integer,intent(in)      :: ng_is(3),ng_ie(3),lg_is(3),lg_ie(3)
    integer,intent(in)      :: Nd,ifn,iobs,iter
    real(8),intent(in)      :: f(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                 ng_is(2)-Nd:ng_ie(2)+Nd,&
                                 ng_is(3)-Nd:ng_ie(3)+Nd)
    character(2),intent(in) :: var
    real(8),allocatable     :: save_pl(:,:),save_pl2(:,:)
    integer          :: ii,inum,i1,i1s,i2,i2s
    character(2)     :: plane_name
    character(128)   :: iobs_name,iter_name,save_name
    
    do ii=1,3
      !allocate
      if(ii==1)     then !xy plane
        i1s=1; i2s=2; plane_name='xy';
      elseif(ii==2) then !yz plane
        i1s=2; i2s=3; plane_name='yz';
      elseif(ii==3) then !xz plane
        i1s=1; i2s=3; plane_name='xz';
      end if
      allocate(save_pl(lg_is(i1s):lg_ie(i1s),lg_is(i2s):lg_ie(i2s)),&
               save_pl2(lg_is(i1s):lg_ie(i1s),lg_is(i2s):lg_ie(i2s)))
      save_pl(:,:)=0.0d0; save_pl2(:,:)=0.0d0
      inum=(lg_ie(i1s)-lg_is(i1s)+1)*(lg_ie(i2s)-lg_is(i2s)+1)
      
      !prepare save data
      if(ipl(ii)==1) then
        if(ii==1) then     !xy plane
!$omp parallel
!$omp do private(i1,i2)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            save_pl(i1,i2)=f(i1,i2,id(3))
          end do
          end do
!$omp end do
!$omp end parallel
        elseif(ii==2) then !yz plane
!$omp parallel
!$omp do private(i1,i2)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            save_pl(i1,i2)=f(id(1),i1,i2)
          end do
          end do
!$omp end do
!$omp end parallel
        elseif(ii==3) then !xz plane
!$omp parallel
!$omp do private(i1,i2)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            save_pl(i1,i2)=f(i1,id(2),i2)
          end do
          end do
!$omp end do
!$omp end parallel
        end if
      end if
      call comm_summation(save_pl,save_pl2,inum,nproc_group_global)
      
      !make save data
      if(comm_is_root(nproc_id_global)) then
        write(iobs_name,*) iobs
        write(iter_name,*) iter
        save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(iobs_name))//'_'//var//&
                  '_'//plane_name//'_'//trim(adjustl(iter_name))//'.data'
        open(ifn,file=save_name)
        do i2=lg_is(i2s),lg_ie(i2s)
        do i1=lg_is(i1s),lg_ie(i1s)
          write(ifn,'(I8,I8,1X,E23.15E3)') i1,i2,save_pl2(i1,i2)*conv
        end do
        end do
        close(ifn)
      end if
      
      !deallocate
      deallocate(save_pl,save_pl2)
    end do
    
    return
  end subroutine eh_save_plane
  
  !===========================================================================================
  != save plane data =========================================================================
  subroutine eh_save_plane_integral(id,ipl,conv_e,conv_h,ti,ng_is,ng_ie,dl,Nd,ifn,iobs,ex,ey,ez,hx,hy,hz)
    use salmon_global,   only: base_directory
    use inputoutput,     only: ulength_from_au,uenergy_from_au,utime_from_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    implicit none
    integer,intent(in) :: id(3),ipl(3)
    real(8),intent(in) :: conv_e,conv_h,ti
    integer,intent(in) :: ng_is(3),ng_ie(3)
    real(8),intent(in) :: dl(3)
    integer,intent(in) :: Nd,ifn,iobs
    real(8),intent(in) :: ex(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd),&
                          ey(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd),&
                          ez(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd),&
                          hx(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd),&
                          hy(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd),&
                          hz(ng_is(1)-Nd:ng_ie(1)+Nd, &
                             ng_is(2)-Nd:ng_ie(2)+Nd, &
                             ng_is(3)-Nd:ng_ie(3)+Nd)
    real(8)          :: ds,&
                        ex_sum1,ex_sum2,ey_sum1,ey_sum2,ez_sum1,ez_sum2,&
                        hx_sum1,hx_sum2,hy_sum1,hy_sum2,hz_sum1,hz_sum2,&
                        px_sum1,px_sum2,py_sum1,py_sum2,pz_sum1,pz_sum2
    integer          :: ii,i1,i1s,i2,i2s
    character(2)     :: plane_name
    character(128)   :: save_name
    
    do ii=1,3
      !set plane
      if(ii==1)     then !xy plane
        i1s=1; i2s=2; plane_name='xy';
      elseif(ii==2) then !yz plane
        i1s=2; i2s=3; plane_name='yz';
      elseif(ii==3) then !xz plane
        i1s=1; i2s=3; plane_name='xz';
      end if
      ds=dl(i1s)*dl(i2s)*(ulength_from_au**2.0d0);

      !prepare integral data
      ex_sum1=0.0d0; ex_sum2=0.0d0; ey_sum1=0.0d0; ey_sum2=0.0d0; ez_sum1=0.0d0; ez_sum2=0.0d0;
      hx_sum1=0.0d0; hx_sum2=0.0d0; hy_sum1=0.0d0; hy_sum2=0.0d0; hz_sum1=0.0d0; hz_sum2=0.0d0;
      px_sum1=0.0d0; px_sum2=0.0d0; py_sum1=0.0d0; py_sum2=0.0d0; pz_sum1=0.0d0; pz_sum2=0.0d0;
      if(ipl(ii)==1) then
        if(ii==1) then     !xy plane
!$omp parallel
!$omp do private(i1,i2) reduction(+:ex_sum1,ey_sum1,ez_sum1,hx_sum1,hy_sum1,hz_sum1,px_sum1,py_sum1,pz_sum1)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            ex_sum1=ex_sum1+ex(i1,i2,id(3)); ey_sum1=ey_sum1+ey(i1,i2,id(3)); ez_sum1=ez_sum1+ez(i1,i2,id(3));
            hx_sum1=hx_sum1+hx(i1,i2,id(3)); hy_sum1=hy_sum1+hy(i1,i2,id(3)); hz_sum1=hz_sum1+hz(i1,i2,id(3));
            px_sum1=px_sum1+ey(i1,i2,id(3))*hz(i1,i2,id(3));
            py_sum1=py_sum1+ez(i1,i2,id(3))*hx(i1,i2,id(3));
            pz_sum1=pz_sum1+ex(i1,i2,id(3))*hy(i1,i2,id(3));
          end do
          end do
!$omp end do
!$omp end parallel
        elseif(ii==2) then !yz plane
!$omp parallel
!$omp do private(i1,i2) reduction(+:ex_sum1,ey_sum1,ez_sum1,hx_sum1,hy_sum1,hz_sum1,px_sum1,py_sum1,pz_sum1)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            ex_sum1=ex_sum1+ex(id(1),i1,i2); ey_sum1=ey_sum1+ey(id(1),i1,i2); ez_sum1=ez_sum1+ez(id(1),i1,i2);
            hx_sum1=hx_sum1+hx(id(1),i1,i2); hy_sum1=hy_sum1+hy(id(1),i1,i2); hz_sum1=hz_sum1+hz(id(1),i1,i2);
            px_sum1=px_sum1+ey(id(1),i1,i2)*hz(id(1),i1,i2);
            py_sum1=py_sum1+ez(id(1),i1,i2)*hx(id(1),i1,i2);
            pz_sum1=pz_sum1+ex(id(1),i1,i2)*hy(id(1),i1,i2);
          end do
          end do
!$omp end do
!$omp end parallel
        elseif(ii==3) then !xz plane
!$omp parallel
!$omp do private(i1,i2) reduction(+:ex_sum1,ey_sum1,ez_sum1,hx_sum1,hy_sum1,hz_sum1,px_sum1,py_sum1,pz_sum1)
          do i2=ng_is(i2s),ng_ie(i2s)
          do i1=ng_is(i1s),ng_ie(i1s)
            ex_sum1=ex_sum1+ex(i1,id(2),i2); ey_sum1=ey_sum1+ey(i1,id(2),i2); ez_sum1=ez_sum1+ez(i1,id(2),i2);
            hx_sum1=hx_sum1+hx(i1,id(2),i2); hy_sum1=hy_sum1+hy(i1,id(2),i2); hz_sum1=hz_sum1+hz(i1,id(2),i2);
            px_sum1=px_sum1+ey(i1,id(2),i2)*hz(i1,id(2),i2);
            py_sum1=py_sum1+ez(i1,id(2),i2)*hx(i1,id(2),i2);
            pz_sum1=pz_sum1+ex(i1,id(2),i2)*hy(i1,id(2),i2);
          end do
          end do
!$omp end do
!$omp end parallel
        end if
      end if
      call comm_summation(ex_sum1,ex_sum2,nproc_group_global)
      call comm_summation(ey_sum1,ey_sum2,nproc_group_global)
      call comm_summation(ez_sum1,ez_sum2,nproc_group_global)
      call comm_summation(hx_sum1,hx_sum2,nproc_group_global)
      call comm_summation(hy_sum1,hy_sum2,nproc_group_global)
      call comm_summation(hz_sum1,hz_sum2,nproc_group_global)
      call comm_summation(px_sum1,px_sum2,nproc_group_global)
      call comm_summation(py_sum1,py_sum2,nproc_group_global)
      call comm_summation(pz_sum1,pz_sum2,nproc_group_global)
      ex_sum2=ex_sum2*ds; ey_sum2=ey_sum2*ds; ez_sum2=ez_sum2*ds;
      hx_sum2=hx_sum2*ds; hy_sum2=hy_sum2*ds; hz_sum2=hz_sum2*ds;
      px_sum2=px_sum2*ds; py_sum2=py_sum2*ds; pz_sum2=pz_sum2*ds;
      
      !save plane integral data
      if(comm_is_root(nproc_id_global)) then
        write(save_name,*) iobs
        save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//&
                  '_'//plane_name//'_integral_rt.data'
        open(ifn,file=save_name,status='old',position='append')
        write(ifn,"(F16.8,99(1X,E23.15E3))",advance='no') &
              ti,ex_sum2*conv_e,ey_sum2*conv_e,ez_sum2*conv_e,hx_sum2*conv_h,hy_sum2*conv_h,hz_sum2*conv_h, &
              px_sum2*(uenergy_from_au/utime_from_au),&
              py_sum2*(uenergy_from_au/utime_from_au),&
              pz_sum2*(uenergy_from_au/utime_from_au)
        close(ifn)
      end if
    end do
    
    return
  end subroutine eh_save_plane_integral
  
  !===========================================================================================
  != calc plane ene data =====================================================================
  subroutine eh_calc_plane_ene(fs,fe,iobs,iter)
    use salmon_global,   only: nt_em,dt_em,obs_samp_em,obs_plane_ene_em
    use structures,      only: s_fdtd_system
    use math_constants,  only: zi
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    integer,            intent(in)    :: iobs,iter
    integer    :: ii,ij,i1,i2
    real(8)    :: t,t_max,wf
    complex(8) :: f_factor
    
    !update time-information, window function
    t        = dble(iter) *dt_em
    t_max    = dble(nt_em)*dt_em
    wf       = 1.0d0 -3.0d0*(t/t_max)**2.0d0 +2.0d0*(t/t_max)**3.0d0
    
    !Fourier transformation
    do ii=1,3
      if(fe%iobs_pl_pe(iobs,ii)==1) then
        do ij=1,fe%iobs_num_ene(iobs)
          !update f_factor
          f_factor = dt_em*dble(obs_samp_em)*exp(zi*obs_plane_ene_em(iobs,ij)*t)
          
          !update time-integration
          if(ii==1)     then !xy plane
!$omp parallel
!$omp do private(i1,i2)
            do i2=fs%mg%is(2),fs%mg%ie(2)
            do i1=fs%mg%is(1),fs%mg%ie(1)
              fe%obs_ex_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ex_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ex_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_ex_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ex_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ex_s(i1,i2,fe%iobs_po_id(iobs,3))
              fe%obs_ey_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ey_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ey_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_ey_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ey_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ey_s(i1,i2,fe%iobs_po_id(iobs,3))
              fe%obs_ez_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ez_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ez_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_ez_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ez_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ez_s(i1,i2,fe%iobs_po_id(iobs,3))
              fe%obs_hx_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hx_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hx_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_hx_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hx_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hx_s(i1,i2,fe%iobs_po_id(iobs,3))
              fe%obs_hy_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hy_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hy_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_hy_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hy_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hy_s(i1,i2,fe%iobs_po_id(iobs,3))
              fe%obs_hz_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hz_xy_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hz_s(i1,i2,fe%iobs_po_id(iobs,3))*wf
              fe%obs_hz_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hz_xy_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hz_s(i1,i2,fe%iobs_po_id(iobs,3))
            end do
            end do
!$omp end do
!$omp end parallel
          elseif(ii==2) then !yz plane
!$omp parallel
!$omp do private(i1,i2)
            do i2=fs%mg%is(3),fs%mg%ie(3)
            do i1=fs%mg%is(2),fs%mg%ie(2)
              fe%obs_ex_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ex_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ex_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_ex_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ex_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ex_s(fe%iobs_po_id(iobs,1),i1,i2)
              fe%obs_ey_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ey_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ey_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_ey_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ey_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ey_s(fe%iobs_po_id(iobs,1),i1,i2)
              fe%obs_ez_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ez_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ez_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_ez_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ez_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ez_s(fe%iobs_po_id(iobs,1),i1,i2)
              fe%obs_hx_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hx_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hx_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_hx_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hx_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hx_s(fe%iobs_po_id(iobs,1),i1,i2)
              fe%obs_hy_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hy_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hy_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_hy_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hy_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hy_s(fe%iobs_po_id(iobs,1),i1,i2)
              fe%obs_hz_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hz_yz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hz_s(fe%iobs_po_id(iobs,1),i1,i2)*wf
              fe%obs_hz_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hz_yz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hz_s(fe%iobs_po_id(iobs,1),i1,i2)
            end do
            end do
!$omp end do
!$omp end parallel
          elseif(ii==3) then !xz plane
!$omp parallel
!$omp do private(i1,i2)
            do i2=fs%mg%is(3),fs%mg%ie(3)
            do i1=fs%mg%is(1),fs%mg%ie(1)
              fe%obs_ex_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ex_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ex_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_ex_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ex_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ex_s(i1,fe%iobs_po_id(iobs,2),i2)
              fe%obs_ey_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ey_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ey_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_ey_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ey_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ey_s(i1,fe%iobs_po_id(iobs,2),i2)
              fe%obs_ez_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ez_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%ez_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_ez_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ez_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%ez_s(i1,fe%iobs_po_id(iobs,2),i2)
              fe%obs_hx_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hx_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hx_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_hx_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hx_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hx_s(i1,fe%iobs_po_id(iobs,2),i2)
              fe%obs_hy_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hy_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hy_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_hy_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hy_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hy_s(i1,fe%iobs_po_id(iobs,2),i2)
              fe%obs_hz_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hz_xz_ene(i1,i2,iobs,ij,1) &
                                                 +f_factor*fe%hz_s(i1,fe%iobs_po_id(iobs,2),i2)*wf
              fe%obs_hz_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hz_xz_ene(i1,i2,iobs,ij,2) &
                                                 +f_factor*fe%hz_s(i1,fe%iobs_po_id(iobs,2),i2)
            end do
            end do
!$omp end do
!$omp end parallel
          end if
        end do
      end if
    end do
    
    return
  end subroutine eh_calc_plane_ene
  
  !===========================================================================================
  != Fourier transformation in eh ============================================================
  subroutine eh_fourier(nt,ne,dt,de,ti,ft,fr,fi)
    use salmon_global,  only: yn_wf_em
    use math_constants, only: zi
    implicit none
    integer,intent(in)   :: nt,ne
    real(8),intent(in)   :: dt,de
    real(8),intent(in)   :: ti(nt),ft(nt)
    real(8),intent(out)  :: fr(ne),fi(ne)
    integer              :: ie,it
    real(8)              :: ft_wf(nt)
    real(8)              :: hw
    complex(8)           :: zf
    
    !apply window function
    if(yn_wf_em=='y') then
      do it=1,nt
        ft_wf(it)=ft(it)*( 1.0d0 -3.0d0*(ti(it)/maxval(ti(:)))**2.0d0 +2.0d0*(ti(it)/maxval(ti(:)))**3.0d0 )
      end do
    else
      ft_wf(:)=ft(:)
    end if
    
    !Fourier transformation
    do ie=1,ne
      hw=dble(ie)*de; zf=(0.0d0,0.0d0);
!$omp parallel
!$omp do private(it) reduction( + : zf )
      do it=1,nt
        zf=zf+exp(zi*hw*ti(it))*ft_wf(it)
      end do
!$omp end do
!$omp end parallel
      zf=zf*dt; fr(ie)=real(zf,8); fi(ie)=aimag(zf)
    end do
    
    return
  end subroutine eh_fourier
  
  !===========================================================================================
  != For ttm =================================================================================
  subroutine calc_es_and_hs(fs, fe)
    use structures, only: s_fdtd_system
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh), intent(inout)   :: fe
    integer :: ix,iy,iz
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=(fs%mg%is_array(3)),(fs%mg%ie_array(3))
    do iy=(fs%mg%is_array(2)),(fs%mg%ie_array(2))
    do ix=(fs%mg%is_array(1)),(fs%mg%ie_array(1))
       fe%ex_s(ix,iy,iz)=fe%ex_y(ix,iy,iz)+fe%ex_z(ix,iy,iz)
       fe%ey_s(ix,iy,iz)=fe%ey_z(ix,iy,iz)+fe%ey_x(ix,iy,iz)
       fe%ez_s(ix,iy,iz)=fe%ez_x(ix,iy,iz)+fe%ez_y(ix,iy,iz)
       fe%hx_s(ix,iy,iz)=( fe%hx_s(ix,iy,iz)+(fe%hx_y(ix,iy,iz)+fe%hx_z(ix,iy,iz)) )*0.5d0
       fe%hy_s(ix,iy,iz)=( fe%hy_s(ix,iy,iz)+(fe%hy_z(ix,iy,iz)+fe%hy_x(ix,iy,iz)) )*0.5d0
       fe%hz_s(ix,iy,iz)=( fe%hz_s(ix,iy,iz)+(fe%hz_x(ix,iy,iz)+fe%hz_y(ix,iy,iz)) )*0.5d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    call eh_sendrecv(fs,fe,'s')
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%ex_s )
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%ey_s )
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%ez_s )
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%hx_s )
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%hy_s )
    call update_overlap_real8( fs%srg_ng, fs%mg, fe%hz_s )
  end subroutine calc_es_and_hs

  !===========================================================================================
  != For ttm =================================================================================
  subroutine allocate_poynting(fs, S, divS, u)
    use structures, only: s_fdtd_system
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    real(8),allocatable,optional,intent(inout) :: S(:,:,:,:)
    real(8),allocatable,optional,intent(inout) :: divS(:,:,:)
    real(8),allocatable,optional,intent(inout) :: u(:,:,:)
    integer :: is1,is2,is3,ie1,ie2,ie3
    is1=fs%mg%is_array(1); ie1=fs%mg%ie_array(1)
    is2=fs%mg%is_array(2); ie2=fs%mg%ie_array(2)
    is3=fs%mg%is_array(3); ie3=fs%mg%ie_array(3)
    if ( present(S) ) then
       if ( allocated(S) ) deallocate(S)
       allocate( S(is1:ie1,is2:ie2,is3:ie3,3) ); S=0.0d0
    end if
    is1=fs%mg%is(1); ie1=fs%mg%ie(1)
    is2=fs%mg%is(2); ie2=fs%mg%ie(2)
    is3=fs%mg%is(3); ie3=fs%mg%ie(3)
    if ( present(divS) ) then
       if ( allocated(divS) ) deallocate(divS)
       allocate( divS(is1:ie1,is2:ie2,is3:ie3) ); divS=0.0d0
    end if
    if ( present(u) ) then
       if ( allocated(u) ) deallocate(u)
       allocate( u(is1:ie1,is2:ie2,is3:ie3) ); u=0.0d0
    end if
  end subroutine allocate_poynting

  !===========================================================================================
  != For ttm =================================================================================
  subroutine calc_poynting_vector(fs, fe, Spoynting)
    use structures, only: s_fdtd_system
    use sendrecv_grid, only: update_overlap_real8
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh), intent(inout)   :: fe
    real(8), intent(inout) :: Spoynting(:,:,:,:)
    integer :: ix,iy,iz,jx,jy,jz,ix0,iy0,iz0,ix1,iy1,iz1
    ix0=fs%mg%is_array(1)
    iy0=fs%mg%is_array(2)
    iz0=fs%mg%is_array(3)
    ix1=fs%mg%ie_array(1)
    iy1=fs%mg%ie_array(2)
    iz1=fs%mg%ie_array(3)
!$omp parallel
!$omp do private(ix,iy,iz,jx,jy,jz)
    do iz = iz0, iz1
       jz=iz-iz0+1
    do iy = iy0, iy1
       jy=iy-iy0+1
    do ix = ix0, ix1
       jx=ix-ix0+1
       Spoynting(jx,jy,jz,1) = fe%ey_s(ix,iy,iz)*fe%hz_s(ix,iy,iz) - fe%ez_s(ix,iy,iz)*fe%hy_s(ix,iy,iz)
       Spoynting(jx,jy,jz,2) = fe%ez_s(ix,iy,iz)*fe%hx_s(ix,iy,iz) - fe%ex_s(ix,iy,iz)*fe%hz_s(ix,iy,iz)
       Spoynting(jx,jy,jz,3) = fe%ex_s(ix,iy,iz)*fe%hy_s(ix,iy,iz) - fe%ey_s(ix,iy,iz)*fe%hx_s(ix,iy,iz)
    end do
    end do
    end do
!$omp end do
!$omp end parallel
  end subroutine calc_poynting_vector

  !===========================================================================================
  != For ttm =================================================================================
  subroutine calc_poynting_vector_div(fs, Spoynting, divS)
    use structures, only: s_fdtd_system
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    real(8), intent(in) :: Spoynting(:,:,:,:)
    real(8), intent(out) :: divS(:,:,:)
    integer :: ix,iy,iz,jx,jy,jz,ix0,iy0,iz0,jx0,jy0,jz0
    real(8) :: cx,cy,cz
    cx = 0.5d0/fs%hgs(1)
    cy = 0.5d0/fs%hgs(2)
    cz = 0.5d0/fs%hgs(3)
    ix0=fs%mg%is(1)-1
    iy0=fs%mg%is(2)-1
    iz0=fs%mg%is(3)-1
    jx0=fs%mg%is_array(1)-1
    jy0=fs%mg%is_array(2)-1
    jz0=fs%mg%is_array(3)-1
!$omp parallel
!$omp do private(jx,jy,jz)
    do iz=fs%mg%is(3),fs%mg%ie(3)
       jz=iz-jz0
    do iy=fs%mg%is(2),fs%mg%ie(2)
       jy=iy-jy0
    do ix=fs%mg%is(1),fs%mg%ie(1)
       jx=ix-jx0
       divS(ix-ix0,iy-iy0,iz-iz0) &
            = ( Spoynting(jx+1,jy,jz,1) - Spoynting(jx-1,jy,jz,1) )*cx &
            + ( Spoynting(jx,jy+1,jz,2) - Spoynting(jx,jy-1,jz,2) )*cy &
            + ( Spoynting(jx,jy,jz+1,3) - Spoynting(jx,jy,jz-1,3) )*cz
    end do
    end do
    end do
!$omp end do
!$omp end parallel
  end subroutine calc_poynting_vector_div
  
  !===========================================================================================
  != make shape ==============================================================================
  subroutine eh_make_shape(fs,fe)
    use salmon_global,   only: yn_output_shape,n_s,typ_s,id_s,inf_s,ori_s,rot_s,&
                               yn_copy_x,yn_copy_y,yn_copy_z,rot_type
    use inputoutput,     only: ulength_from_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_fdtd_system
    use math_constants,  only: pi
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    integer             :: icopy_num(3),i1d_tmp1(6),i1d_tmp2(6)
    real(8)             :: rot_s_d(1000,3)
    real(8),allocatable :: rmove_x(:,:),rmove_y(:,:),rmove_z(:,:)
    integer(8) :: ii
    integer    :: ij,ik,ix,iy,iz,ip_x,ip_y,ip_z,il,l_max
    real(8)    :: x,y,z,x_o,y_o,z_o,x_tmp,y_tmp,z_tmp,adj_err,cal_tmp
    
    !convert from degree to radian or keep radian
    if(trim(rot_type)=='degree') then
      rot_s_d(:,:) = (rot_s(:,:)/360.0d0) * (2.0d0*pi)
    else
      rot_s_d(:,:) = rot_s(:,:)
    end if
    
    !make icopy_num and l_max
    if    (yn_copy_x=='y') then
      icopy_num(1) = 3
    elseif(yn_copy_x=='n') then
      icopy_num(1) = 1     
    end if
    if    (yn_copy_y=='y') then
      icopy_num(2) = 3
    elseif(yn_copy_y=='n') then
      icopy_num(2) = 1     
    end if
    if    (yn_copy_z=='y') then
      icopy_num(3) = 3
    elseif(yn_copy_z=='n') then
      icopy_num(3) = 1     
    end if
    l_max = icopy_num(1)*icopy_num(2)*icopy_num(3)
    
    !make move matrix
    allocate(rmove_x(n_s,l_max),rmove_y(n_s,l_max),rmove_z(n_s,l_max))
    do ii=1,n_s
      il = 1
      do iz=1,icopy_num(3)
      do iy=1,icopy_num(2)
      do ix=1,icopy_num(1)
        ip_x = ix-1; ip_y = iy-1; ip_z = iz-1;
        if(yn_copy_x=='y') ip_x = ip_x-1;
        if(yn_copy_y=='y') ip_y = ip_y-1;
        if(yn_copy_z=='y') ip_z = ip_z-1;
        x_tmp = fs%rlsize(1) * dble(ip_x)
        y_tmp = fs%rlsize(2) * dble(ip_y)
        z_tmp = fs%rlsize(3) * dble(ip_z)
        call rotate_x(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_x(ii,il))
        call rotate_y(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_y(ii,il))
        call rotate_z(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),rmove_z(ii,il))
        il = il+1
      end do
      end do
      end do
    end do
    
    !set adjust parameter
    adj_err = 1.0d-6
    
    !make shape
    cal_tmp = 0.0d0
    do ii=1,n_s
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        !move origin
        x_tmp = fe%coo(ix,1) - ori_s(ii,1)
        y_tmp = fe%coo(iy,2) - ori_s(ii,2)
        z_tmp = fe%coo(iz,3) - ori_s(ii,3)
        
        !rotation
        call rotate_x(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),x_o)
        call rotate_y(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),y_o)
        call rotate_z(x_tmp,y_tmp,z_tmp,rot_s_d(ii,:),z_o)
        
        !copy loop
        do il=1,l_max
          !determine point
          x = x_o + rmove_x(ii,il)
          y = y_o + rmove_y(ii,il)
          z = z_o + rmove_z(ii,il)
          
          !determine shape
          if    (trim(typ_s(ii))=='ellipsoid')            then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0 + (z/(inf_s(ii,3)/2.0d0))**2.0d0
            if(cal_tmp<=1.0d0) fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='half-ellipsoid')       then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0 + (z/(inf_s(ii,3)      ))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-adj_err)) fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-cylinder')    then
            cal_tmp = (x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0)) fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='triangular-cylinder')  then
            if( (x>= -inf_s(ii,1)/2.0d0).and.(x<=inf_s(ii,1)/2.0d0).and.                  &
                (y>= -inf_s(ii,2)/3.0d0).and.                                             &
                (y<=( inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x + inf_s(ii,2)*2.0d0/3.0d0 )).and. &
                (y<=(-inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x + inf_s(ii,2)*2.0d0/3.0d0 )).and. &
                (z>= -inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0) )                     &
              fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='rectangular-cylinder') then
            if( (x>=-inf_s(ii,1)/2.0d0).and.(x<=inf_s(ii,1)/2.0d0).and. &
                (y>=-inf_s(ii,2)/2.0d0).and.(y<=inf_s(ii,2)/2.0d0).and. &
                (z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0) )    &
              fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-cone')        then
            if(inf_s(ii,3)-z/=0.0d0) then 
              cal_tmp= (x/( inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) ))**2.0d0 &
                     + (y/( inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) ))**2.0d0
            else
              cal_tmp=10
            end if
            if( (cal_tmp<=1.0d0).and.(z>=-adj_err).and.(z<=inf_s(ii,3)) ) fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='triangular-cone')      then
            if( (x>= -inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (x<=  inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (y>= -inf_s(ii,2)/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3))        .and. &
                (y<=( inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x                           &
                     +inf_s(ii,2)*2.0d0/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) )).and. &
                (y<=(-inf_s(ii,2)/(inf_s(ii,1)/2.0d0)*x                           &
                     +inf_s(ii,2)*2.0d0/3.0d0*(inf_s(ii,3)-z)/inf_s(ii,3) )).and. &
                (z>=-adj_err).and.(z<=inf_s(ii,3)) )                               &
              fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='rectangular-cone')     then
            if( (x>= -inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (x<=  inf_s(ii,1)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (y>= -inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (y<=  inf_s(ii,2)/2.0d0*(inf_s(ii,3)-z)/inf_s(ii,3)).and. &
                (z>=-adj_err).and.(z<=inf_s(ii,3)) )                      &
              fs%imedia(ix,iy,iz) = id_s(ii)
          elseif(trim(typ_s(ii))=='elliptic-ring')        then
            cal_tmp=(x/(inf_s(ii,1)/2.0d0))**2.0d0 + (y/(inf_s(ii,2)/2.0d0))**2.0d0
            if((cal_tmp<=1.0d0).and.(z>=-inf_s(ii,3)/2.0d0).and.(z<=inf_s(ii,3)/2.0d0)) then
              cal_tmp=(x/(inf_s(ii,4)/2.0d0))**2.0d0 + (y/(inf_s(ii,5)/2.0d0))**2.0d0
              if(cal_tmp>=1.0d0) fs%imedia(ix,iy,iz) = id_s(ii)
            end if
          end if
        end do
      end do
      end do
      end do
    end do
    
    !output cube file
    if(yn_output_shape=='y') then
      !open cube file and write basic information
      if(comm_is_root(nproc_id_global)) then
        open(fe%ifn,file='./shape.cube')
        write(fe%ifn,'(1X,A)') "An input shape described in a cube file."
        write(fe%ifn,'(1X,A)') "A hydrogen atom is used to set the origin of the model."
        write(fe%ifn,'(i5,3f23.15)') 1,fe%coo(fs%lg%is(1),1)*ulength_from_au, &
                                       fe%coo(fs%lg%is(2),2)*ulength_from_au, &
                                       fe%coo(fs%lg%is(3),3)*ulength_from_au
        write(fe%ifn,'(i5,3f23.15)') fs%lg%num(1),fs%hgs(1)*ulength_from_au,0.0d0,0.0d0
        write(fe%ifn,'(i5,3f23.15)') fs%lg%num(2),0.0d0,fs%hgs(2)*ulength_from_au,0.0d0
        write(fe%ifn,'(i5,3f23.15)') fs%lg%num(3),0.0d0,0.0d0,fs%hgs(3)*ulength_from_au
        write(fe%ifn,'(i5,4f23.15)') 1,1.0d0,0.0d0,0.0d0,0.0d0
      end if
      
      !write 3d data
      ix=fs%lg%is(1); iy=fs%lg%is(2); iz=fs%lg%is(3);
      ij=1; i1d_tmp1(:)=0; i1d_tmp2(:)=0;
      do ii=1,fs%lg%num(1)*fs%lg%num(2)*fs%lg%num(3)
        !collect data
        if( fs%mg%is(1)<=ix .and. ix<=fs%mg%ie(1) .and. &
            fs%mg%is(2)<=iy .and. iy<=fs%mg%ie(2) .and. &
            fs%mg%is(3)<=iz .and. iz<=fs%mg%ie(3) ) then
          i1d_tmp1(ij)=fs%imedia(ix,iy,iz)
        end if
        
        !(write data and reset i1d & ij) or (update ij)
        if(mod(ij,6)==0) then
          call comm_summation(i1d_tmp1(:),i1d_tmp2(:),6,nproc_group_global)
          if(comm_is_root(nproc_id_global)) write(fe%ifn,'(6(1X,I2))', advance="yes") i1d_tmp2(:)
          i1d_tmp1(:)=0; i1d_tmp2(:)=0;
          ij=1;
        else
          ij=ij+1;
        end if
        
        !update iz
        iz=iz+1;
        if(iz>fs%lg%ie(3)) iz=fs%lg%is(3);
        
        !update iy
        if(iz==fs%lg%is(3)) iy=iy+1;
        if(iy>fs%lg%ie(2)) iy=fs%lg%is(2);
        
        !update ix
        if(iz==fs%lg%is(3) .and. iy==fs%lg%is(2)) ix=ix+1;
        
        !final output for special case
        if(ii==fs%lg%num(1)*fs%lg%num(2)*fs%lg%num(3) .and. ij>1) then
          call comm_summation(i1d_tmp1(:),i1d_tmp2(:),6,nproc_group_global)
          do ik=1,(ij-1)
            if(comm_is_root(nproc_id_global)) write(fe%ifn,'(1X,I2)', advance="no") i1d_tmp2(ik)
          end do
        end if
      end do
      
      !close file
      if(comm_is_root(nproc_id_global)) close(fe%ifn)
    end if
    
    return
  contains
    
    !+ CONTAINED IN eh_make_shape ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation around x-axis ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_x(x_in,y_in,z_in,rot_in,x_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: x_out
      
      x_out = cos(rot_in(3))*( cos(rot_in(2))*x_in - sin(rot_in(2))*(cos(rot_in(1))*z_in - sin(rot_in(1))*y_in) ) &
             +sin(rot_in(3))*( sin(rot_in(1))*z_in + cos(rot_in(1))*y_in )
      
    end subroutine rotate_x
    
    !+ CONTAINED IN eh_make_shape ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation around y-axis ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_y(x_in,y_in,z_in,rot_in,y_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: y_out
      
      y_out = -sin(rot_in(3))*( cos(rot_in(2))*x_in - sin(rot_in(2))*( cos(rot_in(1))*z_in - sin(rot_in(1))*y_in) ) &
              +cos(rot_in(3))*( sin(rot_in(1))*z_in + cos(rot_in(1))*y_in )
      
    end subroutine rotate_y
    
    !+ CONTAINED IN eh_make_shape ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ rotation around z-axis ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rotate_z(x_in,y_in,z_in,rot_in,z_out)
      implicit none
      real(8),intent(in)  :: x_in,y_in,z_in
      real(8),intent(in)  :: rot_in(3)
      real(8),intent(out) :: z_out
      
      z_out = sin(rot_in(2))*x_in + cos(rot_in(2))*( cos(rot_in(1))*z_in - sin(rot_in(1))*y_in )
      
    end subroutine rotate_z
    
  end subroutine eh_make_shape
    
end module fdtd_eh
