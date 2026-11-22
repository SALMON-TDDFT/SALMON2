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
  use filesystem,only:resolve_path
  use lorentz_drude, only: s_ld_state, ld_allocate, ld_coefficients, ld_store_old, ld_update, ld_finalize
  use fdtd_eh_kernel, only: eh_fd
  use subpixel_smoothing,only:s_subpixel_profile,read_subpixel_profile,allocate_subpixel_profile, &
    subpixel_at,subpixel_bulk,subpixel_spd_bound
  implicit none
  
  ! Batch-only diagnostic resources; embedded contexts leave these unallocated.
  type s_eh_plane_output
    integer :: axis=0,index=0,unit=0
    real(8),allocatable :: buffer(:,:,:)
  end type s_eh_plane_output

  type ls_fdtd_eh
    type(s_eh_plane_output),allocatable :: plane_output(:)
    logical :: tensor_enabled=.false.,tensor_failed=.false.,subpixel_enabled=.false.
    real(8)::subpixel_reference=0d0
    real(8),allocatable::subpixel_corr(:,:,:,:,:,:)
    logical,allocatable::subpixel_cross(:)
    real(8) :: tensor_rtol=1d-10
    integer :: tensor_maxiter=200
    real(8),allocatable :: tensor_raw(:,:,:,:),tensor_corr(:,:,:,:),tensor_rates(:,:,:)
    real(8),allocatable :: tensor_weight(:,:,:,:,:),tensor_diag(:,:,:,:)
    integer,allocatable :: tensor_owner(:,:,:,:)
    logical             :: snapshot_io_warned=.false. ! owned diagnostic warning only
    logical             :: plane_integral_io_warned=.false. ! owned diagnostic warning only
    logical             :: point_io_warned=.false. ! owned diagnostic warning; not physical checkpoint state
    logical             :: flag_pml        !pml  option: work for a_bc = 'pml'
    logical             :: flag_lr         !LR   option: work for ae_shape1 = 'impulse'or ae_shape2 = 'impulse'
    logical             :: flag_obs        !obs. option: work for obs_num_em > 0
    logical             :: flag_ase        !ase  option: work for ase_num_em > 0
    logical             :: flag_art        !art  option: work for art_num_em > 0
    logical             :: flag_save       !save option: work for flag_obs or flag_ase or flag_art
    logical             :: flag_o_restart  !output data_for_restart: work for checkpoint_interval > 0
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
    !(r1,r2,obs_num_em,max(iobs_num_ene),w/ or w/o window function)
    complex(8),allocatable :: obs_ex_xy_ene(:,:,:,:,:),obs_ex_yz_ene(:,:,:,:,:),obs_ex_xz_ene(:,:,:,:,:), &
                              obs_ey_xy_ene(:,:,:,:,:),obs_ey_yz_ene(:,:,:,:,:),obs_ey_xz_ene(:,:,:,:,:), &
                              obs_ez_xy_ene(:,:,:,:,:),obs_ez_yz_ene(:,:,:,:,:),obs_ez_xz_ene(:,:,:,:,:), &
                              obs_hx_xy_ene(:,:,:,:,:),obs_hx_yz_ene(:,:,:,:,:),obs_hx_xz_ene(:,:,:,:,:), &
                              obs_hy_xy_ene(:,:,:,:,:),obs_hy_yz_ene(:,:,:,:,:),obs_hy_xz_ene(:,:,:,:,:), &
                              obs_hz_xy_ene(:,:,:,:,:),obs_hz_yz_ene(:,:,:,:,:),obs_hz_xz_ene(:,:,:,:,:), &
                              obs_jx_xy_ene(:,:,:,:,:),obs_jx_yz_ene(:,:,:,:,:),obs_jx_xz_ene(:,:,:,:,:), &
                              obs_jy_xy_ene(:,:,:,:,:),obs_jy_yz_ene(:,:,:,:,:),obs_jy_xz_ene(:,:,:,:,:), &
                              obs_jz_xy_ene(:,:,:,:,:),obs_jz_yz_ene(:,:,:,:,:),obs_jz_xz_ene(:,:,:,:,:)
    integer             :: inc_num        !number of incident current source
    integer,allocatable :: inc_po_pe(:)   !processor element at incident current source point
    integer,allocatable :: inc_li_pe(:,:) !processor element at incident current source line
    integer,allocatable :: inc_pl_pe(:,:) !processor element at incident current source plane
    integer,allocatable :: inc_po_id(:,:) !id at incident current source point
    character(16)       :: inc_dist1      !spatial distribution type of inc.
    character(16)       :: inc_dist2      !spatial distribution type of inc.
    real(8),allocatable :: gbeam_width_xy1(:,:),gbeam_width_yz1(:,:),gbeam_width_xz1(:,:),& !2D beam spot of inc.
                           gbeam_width_xy2(:,:),gbeam_width_yz2(:,:),gbeam_width_xz2(:,:),& !2D beam spot of inc.
                           gbeam_width_x1(:),   gbeam_width_y1(:),   gbeam_width_z1(:)   ,& !1D beam spot of inc.
                           gbeam_width_x2(:),   gbeam_width_y2(:),   gbeam_width_z2(:)      !1D beam spot of inc.
    real(8)             :: c2_inc1_xyz(3) !coeff. for inc.1, xyz(1:3) means propa. direc. of the inc.
    real(8)             :: c2_inc2_xyz(3) !coeff. for inc.2, xyz(1:3) means propa. direc. of the inc.
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
    real(8),allocatable :: ex_s(:,:,:),ey_s(:,:,:),ez_s(:,:,:)     !e for save
    real(8),allocatable :: hx_s(:,:,:),hy_s(:,:,:),hz_s(:,:,:)     !h for save
    real(8),allocatable :: c2_jx(:,:,:),c2_jy(:,:,:),c2_jz(:,:,:)  !coeff. for general curr. dens.
    real(8),allocatable :: rmedia(:,:,:)                           !Material information for tmp.
    real(8),allocatable :: time_lr(:)                              !LR: time
    real(8),allocatable :: fr_lr(:,:)                              !LR: Re[f]
    real(8),allocatable :: fi_lr(:,:)                              !LR: Im[f]
    real(8),allocatable :: dip_lr(:,:)                             !LR: dipolemoment
    real(8),allocatable :: curr_lr(:,:)                            !LR: average current density
    real(8),allocatable :: e_lr(:,:)                               !LR: average electric field
    real(8),allocatable :: er_lr(:,:)                              !LR: Re[e_lr]
    real(8),allocatable :: ei_lr(:,:)                              !LR: Im[e_lr]
    real(8),allocatable :: ase_ene(:)                              !ase: sampling energy axis(ase_num_em)
    real(8)             :: cspeed_smedia_ase                       !ase: light speed for surrounding media
    integer,allocatable :: iase_mask_xy(:,:)                       !ase: mask on xy plane(x,y)
    integer,allocatable :: iase_mask_yz(:,:)                       !ase: mask on yz plane(y,z)
    integer,allocatable :: iase_mask_xz(:,:)                       !ase: mask on xz plane(x,z)
    integer             :: id_on_box_surf(3,2)                     !ase: id on closed surface[box shape]
                                                                   !     === 1st index ===
                                                                   !     1:ix on yz plane
                                                                   !     2:iy on xz plane
                                                                   !     3:iz on xy plane
                                                                   !     === 2nd index ===
                                                                   !     1:bottom and 2:top
    integer             :: ipe_on_box_surf(3,2)                    !ase: processor element on closed surface[box shape]
                                                                   !     === 1st index ===
                                                                   !     1:ix on yz plane
                                                                   !     2:iy on xz plane
                                                                   !     3:iz on xy plane
                                                                   !     === 2nd index ===
                                                                   !     1:bottom and 2:top
    !ase: inc. ex-z and hx-z along propagation axis(p.-axis,ase_num_em,w/ or w/o window function)
    complex(8),allocatable :: ase_ex_inc_pa_ene(:,:,:),ase_ey_inc_pa_ene(:,:,:),ase_ez_inc_pa_ene(:,:,:), &
                              ase_hx_inc_pa_ene(:,:,:),ase_hy_inc_pa_ene(:,:,:),ase_hz_inc_pa_ene(:,:,:)
    !ase: inc. ex-z and hx-z on bottom or top(ase_num_em,1:bottom and 2:top,w/ or w/o window function)
    complex(8),allocatable :: ase_ex_inc_bt_ene(:,:,:),ase_ey_inc_bt_ene(:,:,:),ase_ez_inc_bt_ene(:,:,:), &
                              ase_hx_inc_bt_ene(:,:,:),ase_hy_inc_bt_ene(:,:,:),ase_hz_inc_bt_ene(:,:,:)
    !ase: sca. ex-z and hx-z on xy, yz,and xz plane(i1,i2,ase_num_em,1:bottom and 2:top,w/ or w/o window function)
    complex(8),allocatable :: ase_ex_xy_sca_ene(:,:,:,:,:),ase_ex_xz_sca_ene(:,:,:,:,:), &
                              ase_ey_xy_sca_ene(:,:,:,:,:),ase_ey_yz_sca_ene(:,:,:,:,:), &
                              ase_ez_yz_sca_ene(:,:,:,:,:),ase_ez_xz_sca_ene(:,:,:,:,:), &
                              ase_hx_xy_sca_ene(:,:,:,:,:),ase_hx_xz_sca_ene(:,:,:,:,:), &
                              ase_hy_xy_sca_ene(:,:,:,:,:),ase_hy_yz_sca_ene(:,:,:,:,:), &
                              ase_hz_yz_sca_ene(:,:,:,:,:),ase_hz_xz_sca_ene(:,:,:,:,:)
    real(8),allocatable :: art_ene(:)        !art: sampling energy axis(art_num_em)
    real(8)             :: cspeed_smedia_art !art: light speed for surrounding media
    integer             :: id_on_plane(2)    !art: id on plane(1:bottom and 2:top)
    integer             :: ipe_on_plane(2)   !art: processor element on plane(1:bottom and 2:top)
    !art: inc. ex-z and hx-z on bottom or top(art_num_em,1:bottom and 2:top,w/ or w/o window function)
    complex(8),allocatable :: art_ex_inc_bt_ene(:,:,:),art_ey_inc_bt_ene(:,:,:),art_ez_inc_bt_ene(:,:,:), &
                              art_hx_inc_bt_ene(:,:,:),art_hy_inc_bt_ene(:,:,:),art_hz_inc_bt_ene(:,:,:)
    !art: tot. ex-z and hx-z on bottom or top plane(i1,i2,art_num_em,1:bottom and 2:top,w/ or w/o window function)
    complex(8),allocatable :: art_ex_tot_bt_ene(:,:,:,:,:),art_ey_tot_bt_ene(:,:,:,:,:),art_ez_tot_bt_ene(:,:,:,:,:), &
                              art_hx_tot_bt_ene(:,:,:,:,:),art_hy_tot_bt_ene(:,:,:,:,:),art_hz_tot_bt_ene(:,:,:,:,:)
  end type ls_fdtd_eh
  
  public  :: eh_init
  public  :: eh_calc
  public  :: eh_finalize
  private :: eh_mpi_grid_sr
  private :: eh_sendrecv
  public  :: eh_fd
  private :: eh_calc_e_inc
  private :: eh_calc_polarization_inc
  private :: eh_calc_wf
  public  :: eh_save_plane
  public  :: eh_calc_plane_ene
  private :: eh_fourier
  private :: eh_point_output,eh_integral_output

contains

  subroutine eh_integral_output(name,header,time,values,warned)
    use salmon_global,only:unit_system
    use parallelization,only:nproc_id_global
    character(*),intent(in) :: name
    logical,intent(in) :: header
    real(8),optional,intent(in) :: time,values(9)
    logical,optional,intent(inout) :: warned
    integer :: unit,ios,close_status
    logical :: report
    if(header)then
      open(newunit=unit,file=resolve_path(name),iostat=ios)
    else
      open(newunit=unit,file=resolve_path(name),status='old',position='append',iostat=ios)
    endif
    if(ios==0)then
      if(header)then
              if(ios==0)write(unit,'(A)',iostat=ios) "# Real time calculation:"
              if(ios==0)write(unit,'(A)',iostat=ios) "# IE: Plane integration of electric field"
              if(ios==0)write(unit,'(A)',iostat=ios) "# IH: Plane integration of magnetic field"
              if(ios==0)write(unit,'(A)',iostat=ios) "# IP: Plane integration of Poynting vector"
              select case(unit_system)
              case('au','a.u.')
                if(ios==0)write(unit,'("#",99(1X,I0,":",A))',iostat=ios) &
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
                if(ios==0)write(unit,'("#",99(1X,I0,":",A))',iostat=ios) &
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
      else
        write(unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=ios)time,values
      endif
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
    endif
    if(ios/=0)then
      report=.true.
      if(present(warned))report=.not.warned
      if(report)write(*,'(A,I0,A,I0,2A)')'SALMONPY_PLANE_INTEGRAL_WARNING rank=',nproc_id_global, &
        ' iostat=',ios,' file=',trim(name)//'; calculation continues; output may be incomplete'
      if(present(warned))warned=.true.
    endif
  end subroutine eh_integral_output

  subroutine eh_point_output(fe,name,header,time,values)
    use salmon_global,only:unit_system
    use parallelization,only:nproc_id_global
    type(ls_fdtd_eh),intent(inout) :: fe
    character(*),intent(in) :: name
    logical,intent(in) :: header
    real(8),optional,intent(in) :: time,values(6)
    integer :: unit,ios,close_status
    if(header)then
      open(newunit=unit,file=resolve_path(name),iostat=ios)
    else
      open(newunit=unit,file=resolve_path(name),status='old',position='append',iostat=ios)
    endif
    if(ios==0)then
      if(header)then
          if(ios==0)write(unit,'(A)',iostat=ios) "# Real time calculation:"
          if(ios==0)write(unit,'(A)',iostat=ios) "# E: Electric field"
          if(ios==0)write(unit,'(A)',iostat=ios) "# H: Magnetic field"
          select case(unit_system)
          case('au','a.u.')
            if(ios==0)write(unit,'("#",99(1X,I0,":",A))',iostat=ios) &
                  1, "Time[a.u.]",                &
                  2, "E_x[a.u.]",                 &
                  3, "E_y[a.u.]",                 &
                  4, "E_z[a.u.]",                 &
                  5, "H_x[a.u.]",                 &
                  6, "H_y[a.u.]",                 &
                  7, "H_z[a.u.]"
          case('A_eV_fs')
            if(ios==0)write(unit,'("#",99(1X,I0,":",A))',iostat=ios) &
                  1, "Time[fs]",                  &
                  2, "E_x[V/Angstrom]",           &
                  3, "E_y[V/Angstrom]",           &
                  4, "E_z[V/Angstrom]",           &
                  5, "H_x[A/Angstrom]",           &
                  6, "H_y[A/Angstrom]",           &
                  7, "H_z[A/Angstrom]"
          end select
      else
        write(unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=ios)time,values
      endif
      ! Always close our own unit, retaining an earlier write error.
      close(unit,iostat=close_status)
      if(ios==0)ios=close_status
    endif
    if(ios/=0.and..not.fe%point_io_warned)then
      write(*,'(A,I0,A,I0,2A)')'SALMONPY_POINT_OUTPUT_WARNING rank=',nproc_id_global, &
        ' iostat=',ios,' file=',trim(name)//'; calculation continues; output may be incomplete'
      fe%point_io_warned=.true.
    endif
  end subroutine eh_point_output

  !===========================================================================================
  != initialize eh-FDTD ======================================================================
  subroutine eh_planes_open(fs,fe)
    use iso_c_binding,only:c_int32_t,c_double
    use structures,only:s_fdtd_system
    use parallelization,only:nproc_id_global
    use salmon_global,only:obs_samp_em,nt_em,dt_em
    type(s_fdtd_system),intent(in)::fs
    type(ls_fdtd_eh),intent(inout)::fe
    character(64)::magic
    character(32)::file_magic
    character(128)::name
    integer::input,ios,n,m,d,offset,tangent(2),j,a,idx,u
    integer(c_int32_t)::header(21)
    logical::exists,owns
    inquire(file=resolve_path('salmonpy-planes.dat'),exist=exists)
    if(.not.exists)return
    if(.not.fe%flag_save)then
      write(*,*)'SALMONPY_PLANE_WARNING: saved E/H unavailable';return
    endif
    open(newunit=input,file=resolve_path('salmonpy-planes.dat'),status='old',iostat=ios)
    if(ios/=0)return
    read(input,'(a)',iostat=ios)magic
    if(ios/=0.or.trim(magic)/='SALMONPY_NATIVE_PLANES_V1')then
      close(input);write(*,*)'SALMONPY_PLANE_WARNING: invalid plane asset';return
    endif
    n=0
    read(input,*,iostat=ios)n
    if(ios/=0.or.n<1)then
      close(input);return
    endif
    allocate(fe%plane_output(n),stat=ios)
    if(ios/=0)then
      close(input);write(*,*)'SALMONPY_PLANE_WARNING: allocation failed';return
    endif
    do m=1,n
      d=0;offset=-1
      read(input,*,iostat=ios)d,offset
      if(ios/=0)exit
      if(d<1.or.d>3)exit
      idx=fs%lg%is(d)+offset
      if(idx<=fs%lg%is(d).or.idx>fs%lg%ie(d))exit
      fe%plane_output(m)%axis=d;fe%plane_output(m)%index=idx
      owns=(idx>=fs%mg%is(d).and.idx<=fs%mg%ie(d)).or. &
           (idx-1>=fs%mg%is(d).and.idx-1<=fs%mg%ie(d))
      if(.not.owns)cycle
      j=0
      do a=1,3
        if(a==d)cycle
        j=j+1;tangent(j)=a
      enddo
      allocate(fe%plane_output(m)%buffer(6, &
        fs%mg%ie(tangent(1))-fs%mg%is(tangent(1))+1, &
        fs%mg%ie(tangent(2))-fs%mg%is(tangent(2))+1),stat=ios)
      if(ios/=0)cycle
      write(name,'("salmonpy_plane_",i4.4,"_rank_",i8.8,".bin")')m,nproc_id_global
      open(newunit=u,file=resolve_path(trim(name)),status='new',access='stream',form='unformatted',iostat=ios)
      if(ios/=0)then
        write(*,*)'SALMONPY_PLANE_WARNING: cannot create ',trim(name);cycle
      endif
      fe%plane_output(m)%unit=u
      header=int([16909060,1,d,idx,idx-1,fs%lg%is,fs%lg%ie,fs%mg%is,fs%mg%ie, &
                  obs_samp_em,nt_em,nproc_id_global,fe%iter_sta-1],c_int32_t)
      file_magic='SALMONPY_NATIVE_PLANES_V1'
      write(u,iostat=ios)file_magic,header,real(dt_em,c_double)
      if(ios/=0)then
        close(u,iostat=ios);fe%plane_output(m)%unit=0
      endif
    enddo
    close(input)
    if(fe%iter_sta==1)call eh_planes_write(fs,fe,0,.true.)
  end subroutine eh_planes_open

  subroutine eh_planes_write(fs,fe,iter,initial)
    use iso_c_binding,only:c_int32_t
    use structures,only:s_fdtd_system
    type(s_fdtd_system),intent(in)::fs
    type(ls_fdtd_eh),intent(inout)::fe
    integer,intent(in)::iter
    logical,optional,intent(in)::initial
    integer::m,d,u,ios,q,idx,tangent(2),j,a,b,ijk(3),i1,i2,i3
    logical::first
    if(.not.allocated(fe%plane_output))return
    first=.false.
    if(present(initial))first=initial
    do m=1,size(fe%plane_output)
      u=fe%plane_output(m)%unit
      if(u==0)cycle
      d=fe%plane_output(m)%axis;idx=fe%plane_output(m)%index
      j=0
      do a=1,3
        if(a==d)cycle
        j=j+1;tangent(j)=a
      enddo
      write(u,iostat=ios)int(iter,c_int32_t)
      do q=0,1
        if(ios/=0)exit
        if(idx-q<fs%mg%is(d).or.idx-q>fs%mg%ie(d))cycle
        ijk(d)=idx-q
        do b=1,size(fe%plane_output(m)%buffer,3)
          ijk(tangent(2))=fs%mg%is(tangent(2))+b-1
          do a=1,size(fe%plane_output(m)%buffer,2)
            ijk(tangent(1))=fs%mg%is(tangent(1))+a-1
            i1=ijk(1);i2=ijk(2);i3=ijk(3)
            if(first)then
              fe%plane_output(m)%buffer(:,a,b)=[ &
                fe%ex_y(i1,i2,i3)+fe%ex_z(i1,i2,i3),fe%ey_z(i1,i2,i3)+fe%ey_x(i1,i2,i3), &
                fe%ez_x(i1,i2,i3)+fe%ez_y(i1,i2,i3),fe%hx_y(i1,i2,i3)+fe%hx_z(i1,i2,i3), &
                fe%hy_z(i1,i2,i3)+fe%hy_x(i1,i2,i3),fe%hz_x(i1,i2,i3)+fe%hz_y(i1,i2,i3)]
            else
              fe%plane_output(m)%buffer(:,a,b)=[fe%ex_s(i1,i2,i3),fe%ey_s(i1,i2,i3),fe%ez_s(i1,i2,i3), &
                fe%hx_s(i1,i2,i3),fe%hy_s(i1,i2,i3),fe%hz_s(i1,i2,i3)]
            endif
          enddo
        enddo
        write(u,iostat=ios)fe%plane_output(m)%buffer
      enddo
      if(ios/=0)then
        write(*,*)'SALMONPY_PLANE_WARNING: diagnostic write failed at step ',iter
        close(u,iostat=ios);fe%plane_output(m)%unit=0
      endif
    enddo
  end subroutine eh_planes_write

  subroutine eh_planes_close(fe)
    type(ls_fdtd_eh),intent(inout)::fe
    integer::m,ios
    if(.not.allocated(fe%plane_output))return
    do m=1,size(fe%plane_output)
      if(fe%plane_output(m)%unit/=0)close(fe%plane_output(m)%unit,iostat=ios)
    enddo
    deallocate(fe%plane_output)
  end subroutine eh_planes_close

  subroutine eh_init(fs,fe,ld,status,grid_ready,excluded_response_points &
#ifdef SALMONPY_RECOVERY_TESTING
                     ,allocation_failure_slot &
#endif
                     )
    use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
    use salmon_global,   only: al_em,dl_em,num_rgrid_em,at_em,dt_em,nt_em,boundary_em,yn_periodic,         &
                               base_directory,media_num,shape_file,epsilon_em,mu_em,sigma_em,media_type,   &
                               pole_num_ld,omega_p_ld,f_ld,gamma_ld,omega_ld,                              &
                               obs_num_em,obs_loc_em,obs_plane_ene_em,yn_obs_plane_integral_em,obs_samp_em,&
                               media_id_pml,media_id_source1,media_id_source2,                             &
                               wave_input,trans_longi,e_impulse,nenergy,                                   &
                               source_loc1,ek_dir1,epdir_re1,epdir_im1,ae_shape1,omega1,                   &
                               phi_cep1,I_wcm2_1,E_amplitude1,gbeam_sigma_plane1,gbeam_sigma_line1,        &
                               source_loc2,ek_dir2,epdir_re2,epdir_im2,ae_shape2,omega2,                   &
                               phi_cep2,I_wcm2_2,E_amplitude2,gbeam_sigma_plane2,gbeam_sigma_line2,        &
                               bloch_k_em,bloch_real_imag_em,yn_make_shape,yn_output_shape,                &
                               ase_num_em,ase_ene_min_em,ase_ene_max_em,ase_wav_min_em,ase_wav_max_em,     &
                               ase_smedia_id_em,ase_box_cent_em,ase_box_size_em,                           &
                               art_num_em,art_ene_min_em,art_ene_max_em,art_wav_min_em,art_wav_max_em,     &
                               art_smedia_id_em,art_plane_bot_em,art_plane_top_em,                         &
                               yn_restart,directory_read_data,checkpoint_interval
    use inputoutput,     only: utime_from_au,ulength_from_au,uenergy_from_au,unit_system,&
                               uenergy_to_au,ulength_to_au,ucharge_to_au
    use parallelization, only: nproc_id_global,nproc_group_global,nproc_size_global
    use communication,   only: comm_is_root,comm_bcast,comm_summation,comm_sync_all,comm_get_max,comm_get_min
    use structures,      only: s_fdtd_system
    use filesystem,      only: create_directory,directory_exists
    use misc_routines,   only: get_wtime
    use phys_constants,  only: cspeed_au
    use math_constants,  only: pi
    use common_maxwell,  only: set_coo_em,find_point_line_plane_em,find_point_id_only_em,make_shape_em,          &
                               input_i3d_em,output_i3d_em,stop_em,input_r_txt_em,output_r_txt_em,input_r_bin_em, &
                               txtfile_copy_em
    implicit none
    integer :: restart_io,restart_bad,restart_total
    type(s_ld_state), intent(inout) :: ld
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    ! Optional native early-return contract. Later I/O/allocation stages retain
    ! their existing behavior; callers without status keep standalone diagnostics.
    integer,optional,intent(in)::excluded_response_points(:,:)
    integer,optional,intent(out) :: status
    logical,optional,intent(out) :: grid_ready
    integer             :: nsg_p          !yn_restart='y': nproc_size_global used in the previous calc.
    integer             :: nt_em_p        !yn_restart='y': nt_em             used in the previous calc.
    integer             :: iobs_samp_em_p !yn_restart='y': obs_samp_em       used in the previous calc.
    integer,allocatable :: is_p(:,:)      !yn_restart='y': mg%is             used in the previous calc.
    integer,allocatable :: ie_p(:,:)      !yn_restart='y': mg%ie             used in the previous calc.
    integer,allocatable :: is_a_p(:,:)    !yn_restart='y': mg%is_array       used in the previous calc.
    integer,allocatable :: ie_a_p(:,:)    !yn_restart='y': mg%ie_array       used in the previous calc.
    logical             :: flag_same_p    !yn_restart='y': same parallel condition compared with the previous calc.
    integer             :: ii,ij,ik,ix,iy,iz,icount,icount_ld,iroot1,iroot2,ibt
    integer             :: itmp1(3),itmp2(3)
    integer             :: is5(5),ie5(5)
    real(8)             :: dt_cfl,elapsed_time,tmp_min,tmp_max,tmp_x1,tmp_x2,tmp_y1,tmp_y2,tmp_z1,tmp_z2
    real(8),allocatable :: tmp_c0_ld(:,:),tmp_c1_ld(:,:),tmp_c2_ld(:,:),tmp_c3_ld(:,:)
    real(8),allocatable :: tmp_c4_ld(  :),tmp_c5_ld(  :)
    real(8),allocatable :: diagonal_eps(:,:),diagonal_f(:,:,:)
    real(8),allocatable :: diagonal_c2(:,:,:),diagonal_c4(:,:),diagonal_c5(:,:)
    type(s_subpixel_profile)::subpixel_profile
    integer :: observation_unit
    integer :: allocation_status,ld_allocation_status
#ifdef SALMONPY_RECOVERY_TESTING
    integer,optional,intent(in) :: allocation_failure_slot
    integer :: allocation_slot,ld_failure_slot,auxiliary_slot
#endif
    integer :: diagonal_status,diagonal_axis,diagonal_media,shape_status
    character(1)        :: dir
    character(2)        :: plane_name
    character(16)       :: tmp_name1,tmp_name2
    character(256)      :: save_name
    character(256)      :: tmp_c(2)
    logical             :: flag_stop
    
    if(present(status))status=0
    allocation_status=0
#ifdef SALMONPY_RECOVERY_TESTING
    allocation_slot=0
    auxiliary_slot=0
#endif
    if(present(grid_ready))grid_ready=.false.

    !*** set flag *********************************************************************************************!
    !pml option: use absorbing baundary(PML) condition(this flag would be updated later)
    fe%flag_pml = .false.
    
    !LD option: calculate Lorentz-Drude model
    ld%num_ld=0
    do ii=0,media_num
      select case(media_type(ii))
      case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
        ld%num_ld=ld%num_ld+1
      end select
    end do
    if(ld%num_ld>0) then
      ld%enabled = .true.
    else
      ld%enabled = .false.
    end if
    
    !LR option: execute linear response calculation under the dipole approximation
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      fe%flag_lr = .true.
    else
      fe%flag_lr = .false.
    end if
    
    !observation option: set observation point
    if(obs_num_em>0) then
      fe%flag_obs = .true.
    else
      fe%flag_obs = .false.
    end if
    
    !ase option: calculate absorption-, scattering-, and extinction-cross-sections
    if(ase_num_em>0) then
      fe%flag_ase = .true.
    else
      fe%flag_ase = .false.
    end if
    
    !art option: calculate absorption-, reflection-, and transmission-rates
    if(art_num_em>0) then
      fe%flag_art = .true.
    else
      fe%flag_art = .false.
    end if
    
    !save option: calculate variables used for save(typically, ex-z_s and hx-z_s)
    if(fe%flag_obs .or. fe%flag_ase .or. fe%flag_art) then
      fe%flag_save = .true.
    else
      fe%flag_save = .false.
    end if
    
    !output data_for_restart
    if(checkpoint_interval>0) then
      fe%flag_o_restart = .true.
    else
      fe%flag_o_restart = .false.
    end if
    
    !same parallel condition compared with the previous calc.(this flag will be updated later for yn_restart='y')
    flag_same_p = .true.
    
    !temporarily used flag for stop(this flag would be updated later)
    flag_stop = .false.
    
    !*** set initial parameter and value **********************************************************************!
    fe%Nd     = 1
    fe%ifn    = 600
    fe%ipml_l = 8
    fe%pml_m  = 4.0d0
    fe%pml_r  = 1.0d-7
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
      if(present(status))then
        status=1101;return
      endif
      call stop_em('Only two of al_em, dl_em, and num_rgrid_em must be set.')
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
          fe%flag_pml    = .true.
        elseif(yn_periodic=='y') then
          fs%a_bc(ii,ij) = 'periodic'
        end if
      case('abc')
        fs%a_bc(ii,ij) = 'pml'
        fe%flag_pml    = .true.
      case('pec')
        if(yn_periodic=='y') then
          if(present(status))then
            status=1102;return
          endif
          call stop_em('For yn_periodic = y, boundary_em must be default, periodic or abc.')
        end if
        fs%a_bc(ii,ij) = 'pec'
      case('periodic')
        if(yn_periodic=='n'.and.boundary_em(ii,3-ij)/='periodic') then
          if(present(status))then
            status=1102;return
          endif
          call stop_em('Periodic boundary_em faces must be paired when yn_periodic=n.')
        end if
        fs%a_bc(ii,ij) = 'periodic'
      end select
    end do
    end do
    select case(unit_system)
    case('au','a.u.')
      fe%uVperm_from_au=1.0d0
      fe%uAperm_from_au=1.0d0
    case('A_eV_fs')
      !see E_amplitude1 or E_amplitude2 in src/io/iunputoutput.f90
      fe%uVperm_from_au=1.0d0/(uenergy_to_au/ulength_to_au/ucharge_to_au)
      fe%uAperm_from_au=fe%uVperm_from_au
    end select
    
    !*** prepare mpi, gird, and sendrecv environments *********************************************************!
    if(present(status))then
      call eh_mpi_grid_sr(fs,fe,status)
      if(status/=0)return
    else
      call eh_mpi_grid_sr(fs,fe)
    endif
    if(present(grid_ready))grid_ready=.true.
    
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
    
    !*** check and set time information ***********************************************************************!
    !set maximum dt in CFL condition
    dt_cfl=1.0d0/( &
           cspeed_au*sqrt( (1.0d0/fs%hgs(1))**2.0d0+(1.0d0/fs%hgs(2))**2.0d0+(1.0d0/fs%hgs(3))**2.0d0 ) )
    
    ! Explicit diagonal tensor records are read once into init-owned buffers.
    if(any(media_type(0:media_num)=='ld-diagonal'))then
      call eh_read_diagonal(diagonal_status)
      if(diagonal_status/=0)then
        if(present(status))then
          status=merge(1120,1110,diagonal_status==2);return
        endif
        if(diagonal_status==2)error stop 'Required LD material coefficient allocation failed'
        call stop_em('Invalid or missing salmonpy-ld-diagonal.dat tensor input.')
      endif
      if(comm_is_root(nproc_id_global))write(*,'(A)')'SALMONPY_FDTD_LD_DIAGONAL_V1'
    endif
    if(any(media_type(0:media_num)=='ld-tensor'))then
      call eh_read_tensor(diagonal_status)
      if(diagonal_status/=0)then
        if(present(status))then
          status=merge(1120,1110,diagonal_status==2);return
        endif
        if(diagonal_status==2)error stop 'Required LD material coefficient allocation failed'
        call stop_em('Invalid or missing salmonpy-ld-tensor.dat input.')
      endif
      if(comm_is_root(nproc_id_global))write(*,'(A)')'SALMONPY_FDTD_LD_TENSOR_V1'
    endif
    if(any(media_type(0:media_num)=='ld-subpixel'))then
      call eh_read_subpixel(diagonal_status)
      if(diagonal_status/=0)then
        if(present(status))then
          status=1113;return
        endif
        call stop_em('Invalid or missing salmonpy-subpixel.dat profile.')
      endif
    endif
    ! Keep vacuum/PML in the bound; faster positive backgrounds tighten it.
    do ii=0,media_num
      if(epsilon_em(ii)>0d0.and.mu_em(ii)>0d0)then
        if(epsilon_em(ii)*mu_em(ii)<1d0) &
          dt_cfl=min(dt_cfl,sqrt(epsilon_em(ii)*mu_em(ii))/( &
          cspeed_au*sqrt(sum((1d0/fs%hgs)**2))))
      endif
    enddo

    if(allocated(diagonal_eps))then
      do ii=0,media_num
        if(media_type(ii)/='ld-diagonal'.and.media_type(ii)/='ld-tensor')cycle
        do ij=1,3
          if(diagonal_eps(ij,ii)<1d0)dt_cfl=min(dt_cfl,sqrt(diagonal_eps(ij,ii))/( &
            cspeed_au*sqrt(sum((1d0/fs%hgs)**2))))
        enddo
      enddo
    endif
    !check and set at_em, dt_em, and nt_em
    if(comm_is_root(nproc_id_global)) write(*,*)
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    if( ( (at_em >0.0d0) .and. (dt_em >0.0d0) .and. (nt_em>0) ) &
        .or.                                                    &
        ( (at_em==0.0d0) .and. (nt_em==0)                     ) ) then
      if(present(status))then
        status=1103;return
      endif
      call stop_em(   'Two of at_em, dt_em, and nt_em must be set.',                    &
                   c2='Otherwise, both at_em and nt_em or either of those must be set.',&
                   c3='(then, dt_em is automatically determined)')
    elseif(dt_em==0.0d0) then
      if( (at_em>0.0d0) .and. (nt_em>0) ) then
        dt_em = at_em/dble(nt_em)
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "From at_em and nt_em, dt_em is determined as", dt_em*utime_from_au
          write(*,*) "in the unit system, ",trim(unit_system),"."
        end if
      else
        dt_em=dt_cfl*0.999d0
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "From CFL condition, dt_em is determined as", dt_em*utime_from_au
          write(*,*) "in the unit system, ",trim(unit_system),"."
        end if
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "dt_em =", dt_em*utime_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
      end if
    end if
    call comm_bcast(dt_em,nproc_group_global)
    if(dt_em>=dt_cfl) then
      write(tmp_c(1),*) dt_cfl*utime_from_au
      write(tmp_c(2),*) 'in the unit system, ',trim(unit_system),'.'
      if(present(status))then
        status=1104;return
      endif
      call stop_em('To sufficient CFL condition, dt_em must be set smaller than',&
                   c2=trim(adjustl(tmp_c(1))), c3=trim(adjustl(tmp_c(2))) )
    end if
    if    (at_em==0.0d0) then
      at_em = dt_em * dble(nt_em)
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "From dt_em and nt_em, at_em is determined as", at_em*utime_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
      end if
    elseif(nt_em==0    ) then
      nt_em = nint(at_em/dt_em)
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "From at_em and dt_em, nt_em is determined as", nt_em
      end if
    end if
    
    !check and set obs_samp_em
    if(obs_samp_em==0) then
      if    (fe%flag_ase) then
        if     (ase_ene_max_em>0.0d0) then
          obs_samp_em = nint( ((2*pi/ase_ene_max_em)/20.0d0) / dt_em )
          if(obs_samp_em<1) obs_samp_em=1
          if(comm_is_root(nproc_id_global)) write(*,*) "From ase_ene_max_em, obs_samp_em is determined as", obs_samp_em
        else if(ase_wav_min_em>0.0d0) then
          obs_samp_em = nint( ((ase_wav_min_em/cspeed_au)/20.0d0) / dt_em )
          if(obs_samp_em<1) obs_samp_em=1
          if(comm_is_root(nproc_id_global)) write(*,*) "From ase_wav_min_em, obs_samp_em is determined as", obs_samp_em
        end if
      elseif(fe%flag_art) then
        if     (art_ene_max_em>0.0d0) then
          obs_samp_em = nint( ((2*pi/art_ene_max_em)/20.0d0) / dt_em )
          if(obs_samp_em<1) obs_samp_em=1
          if(comm_is_root(nproc_id_global)) write(*,*) "From art_ene_max_em, obs_samp_em is determined as", obs_samp_em
        else if(art_wav_min_em>0.0d0) then
          obs_samp_em = nint( ((art_wav_min_em/cspeed_au)/20.0d0) / dt_em )
          if(obs_samp_em<1) obs_samp_em=1
          if(comm_is_root(nproc_id_global)) write(*,*) "From art_wav_min_em, obs_samp_em is determined as", obs_samp_em
        end if
      elseif(omega1>0.0d0) then
        obs_samp_em = nint( ((2*pi/omega1)/20.0d0) / dt_em )
        if(obs_samp_em<1) obs_samp_em=1
        if(comm_is_root(nproc_id_global)) write(*,*) "From omega1, obs_samp_em is determined as", obs_samp_em
      elseif(omega2>0.0d0) then
        obs_samp_em = nint( ((2*pi/omega2)/20.0d0) / dt_em )
        if(obs_samp_em<1) obs_samp_em=1
        if(comm_is_root(nproc_id_global)) write(*,*) "From omega2, obs_samp_em is determined as", obs_samp_em
      else
        obs_samp_em=1
        if(comm_is_root(nproc_id_global)) write(*,*) "obs_samp_em is determined as", obs_samp_em
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "obs_samp_em =", obs_samp_em
      end if
    end if
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    
    !set time-step(start and end)
    fe%iter_sta = 1
    fe%iter_end = nt_em
    
    !*** check and partially input restart date ***************************************************************!
    if(yn_restart=='y') then
      restart_io=0
      !check restart date
        if(comm_is_root(nproc_id_global)) then
          write(*,*)
          write(*,*) "**************************"
          write(*,*) "yn_restart='y':"
          write(*,*) "restart data is loaded from ",trim(adjustl(directory_read_data))
          write(*,'(A)') " --- CAUTION ---"
          write(*,'(A)') " All input keywords except yn_restart, directory_read_data, and checkpoint_interval"
          write(*,'(A)') " are recomended to be same with those used in the previous calculation."
          write(*,'(A)') " Otherwise, unexpected errors would appear."
          write(*,*) "**************************"
        end if
      !input nproc_size_global, nt_em, and obs_samp_em used in the previous calculation
      call input_r_txt_em(fe%ifn,1,1,'nproc_size_global',i0d=nsg_p  , &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_txt_em(fe%ifn,1,1,'nt_em'            ,i0d=nt_em_p, &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_txt_em(fe%ifn,1,1,'obs_samp_em'      ,i0d=iobs_samp_em_p, &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      if(nt_em/=nt_em_p) then
        call stop_em('When yn_restart="y",',c2='nt_em must be same with that used in the previous calculation.')
      end if
      if(obs_samp_em/=iobs_samp_em_p) then
        call stop_em('When yn_restart="y",',c2='obs_samp_em must be same with that used in the previous calculation.')
      end if
      
      !input array size
      allocate(  is_p(1:3,0:(nsg_p-1)),  ie_p(1:3,0:(nsg_p-1)),&
               is_a_p(1:3,0:(nsg_p-1)),ie_a_p(1:3,0:(nsg_p-1)) )
      do ii=0,(nsg_p-1)
        write(tmp_name1,*) ii
        tmp_name2 = 'mg_is_'//trim(adjustl(tmp_name1))
        call input_r_txt_em(fe%ifn,1,3,trim(adjustl(tmp_name2)),i1d=is_p(:,ii), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        tmp_name2 = 'mg_ie_'//trim(adjustl(tmp_name1));
        call input_r_txt_em(fe%ifn,1,3,trim(adjustl(tmp_name2)),i1d=ie_p(:,ii), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        tmp_name2 = 'mg_is_array_'//trim(adjustl(tmp_name1))
        call input_r_txt_em(fe%ifn,1,3,trim(adjustl(tmp_name2)),i1d=is_a_p(:,ii), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        tmp_name2 = 'mg_ie_array_'//trim(adjustl(tmp_name1))
        call input_r_txt_em(fe%ifn,1,3,trim(adjustl(tmp_name2)),i1d=ie_a_p(:,ii), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

      end do
      
      !set flag
      if(nproc_size_global==nsg_p) then
        ii=0; ij=0;
        if( fs%mg%is(1)/=is_p(1,nproc_id_global).or.fs%mg%ie(1)/=ie_p(1,nproc_id_global).or.&
            fs%mg%is(2)/=is_p(2,nproc_id_global).or.fs%mg%ie(2)/=ie_p(2,nproc_id_global).or.&
            fs%mg%is(3)/=is_p(3,nproc_id_global).or.fs%mg%ie(3)/=ie_p(3,nproc_id_global) ) then
          ii=1
        end if
        call comm_summation(ii,ij,nproc_group_global)
        if(ij==0) then
          flag_same_p = .true.
        else
          flag_same_p = .false.
        end if
      else
        flag_same_p = .false.
      end if
    end if
    
    !*** prepare to store restart date ************************************************************************!
    if(fe%flag_o_restart) then
      !check, modify, and write checkpoint_interval
      if(comm_is_root(nproc_id_global)) write(*,*)
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
      if(checkpoint_interval>0) then
        if(checkpoint_interval<obs_samp_em) then
          call stop_em('checkpoint_interval must be larger than obs_samp_em.')
        elseif(mod(checkpoint_interval,obs_samp_em)/=0) then
          checkpoint_interval = nint(dble(checkpoint_interval)/dble(obs_samp_em))*obs_samp_em
          if(comm_is_root(nproc_id_global)) then
            write(*,*) "checkpoint_interval must be a multiple of obs_samp_em."
            write(*,*) "Therefore, checkpoint_interval is modified as",checkpoint_interval
          end if
        else
          if(comm_is_root(nproc_id_global)) write(*,*) "checkpoint_interval =",checkpoint_interval
        end if
      end if
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "data_for_restart will be outputed."
        write(*,*) "**************************"
      end if
      
      ! Replace the creation barrier with one status agreement. A root-only
      ! filesystem failure must not leave other ranks waiting or report success.
      restart_io=0
      if(comm_is_root(nproc_id_global))call create_directory('data_for_restart/',iret=restart_io)
      restart_bad=0
      if(restart_io/=0)restart_bad=1
      call comm_summation(restart_bad,restart_total,nproc_group_global)
      if(restart_total/=0)then
        if(comm_is_root(nproc_id_global))write(*,'(A,I0)') &
          'SALMONPY_RESTART_DIRECTORY_ERROR status=',restart_io
        error stop 'Required EM restart directory creation failed'
      endif

    end if
    
    !*** make or input shape data *****************************************************************************!
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'i3d',i3d=fs%imedia,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%rmedia,allocation_status=allocation_status)
    ! Agree once per allocation phase, before shape/halo communication.
    if(present(status))then
      call comm_get_max(allocation_status,nproc_group_global)
      if(allocation_status/=0)then
        status=1118;return
      endif
    endif
    if(media_num>0) then
      if(comm_is_root(nproc_id_global)) write(*,*)
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
      elapsed_time=get_wtime()
      
      !make or input
      if(yn_make_shape=='y') then
        !make shape data
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "shape data is made by input keywords."
        end if
        call make_shape_em(fs%mg,fs%lg,fs%rlsize(:),fs%imedia,fe%Nd,fe%coo(:,:))
        if(yn_output_shape=='y') then
          save_name = trim(adjustl(base_directory))//'/shape'
          call output_i3d_em(fs%mg,fs%lg,600,1,fs%imedia,'cu',save_name,fe%Nd,fe%coo(:,:))
        end if
      elseif(yn_make_shape=='n') then
        ! Readers finish locally and close their owned units before agreement;
        ! no rank enters halo exchange if any peer failed to obtain its grid.
        if(present(status))then
          shape_status=7
          if(index(shape_file,'.cube',back=.true.)/=0)then
            call input_i3d_em(shape_file,fe%ifn,fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie, &
                              fe%Nd,fs%imedia,'cu',shape_status)
          endif
          call comm_get_max(shape_status,nproc_group_global)
          if(shape_status/=0)then
            status=merge(1116,1115,shape_status==5);return
          endif
        else
        !check file format and input shape data
        if(index(shape_file,".cube", back=.true.)/=0) then
          if(comm_is_root(nproc_id_global)) then
            write(*,*) "shape_file is inputed by .cube format."
          end if
          call input_i3d_em(shape_file,fe%ifn,fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fs%imedia,'cu')
        elseif(index(shape_file,".mp", back=.true.)/=0) then
          call stop_em(   'shape_file is inputed by .mp format.', &
                       c2='This version works for only .cube format.')
        else
          call stop_em('shape_file must be .cube or .mp formats.')
        end if
        endif
      end if
      
      !send and receive for shape data
      fe%rmedia(:,:,:)=dble(fs%imedia(:,:,:))
      call eh_sendrecv(fs,fe,'r')
      fs%imedia(:,:,:)=int(fe%rmedia(:,:,:)+1d-3)
      elapsed_time=get_wtime()-elapsed_time
      if(comm_is_root(nproc_id_global)) then
        write(*,'(A,f16.8)') " elapsed time for making shape data [s] = ", elapsed_time
        write(*,*) "**************************"
      end if
    end if
    
    !*** basic allocation in eh-FDTD **************************************************************************!
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_y   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ex_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ex_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_z   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ex_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ex_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_z   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ey_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ey_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_x   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ey_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ey_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_x   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ez_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ez_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_y   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_ez_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_ez_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_y   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hx_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hx_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_z   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hx_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hx_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_z   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hy_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hy_z,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_x   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hy_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hy_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_x   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hz_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hz_x,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_y   ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c1_hz_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_hz_y,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jx  ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jy  ,allocation_status=allocation_status)
    call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%c2_jz  ,allocation_status=allocation_status)
    if(fe%flag_save) then
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ex_s,allocation_status=allocation_status)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ey_s,allocation_status=allocation_status)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%ez_s,allocation_status=allocation_status)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hx_s,allocation_status=allocation_status)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hy_s,allocation_status=allocation_status)
      call eh_allocate(fs%mg%is_array,fs%mg%ie_array,'r3d',r3d=fe%hz_s,allocation_status=allocation_status)
    end if
    
    if(present(status))then
      call comm_get_max(allocation_status,nproc_group_global)
      if(allocation_status/=0)then
        status=1118;return
      endif
    endif

    !*** prepare Lorentz-Drude ********************************************************************************!
    if(ld%enabled) then
      !set counter, make media_ld, make max_pole_num_ld, and check pole_num_ld & other conditions
      icount_ld=1; ld%max_pole_num_ld=0;
      allocate(ld%media_ld(ld%num_ld))
      ld%media_ld(:)=0;
      do ii=0,media_num
        select case(media_type(ii))
        case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
          ld%media_ld(icount_ld) = ii
          icount_ld              = icount_ld+1
          if(ld%max_pole_num_ld<pole_num_ld(ii)) ld%max_pole_num_ld=pole_num_ld(ii)
          if(pole_num_ld(ii)<=0) then
            call stop_em('For media_type = lorentz-drude, pole_num_ld must be equal to or larger than 1.')
          end if
          if((epsilon_em(ii)<=0.0d0).or.(mu_em(ii)/=1.0d0).or.(sigma_em(ii)/=0.0d0)) then
            call stop_em(   'For media_type = lorentz-drude,', &
                         c2='epsilon_em(ii) > 0, mu_em(ii) = 1.0d0, and sigma_em(ii) = 0.0d0 must be set.')
          end if
        end select
      end do
      
      !reset LD counter
      icount_ld = 1
      
      !set temporal LD coefficients (will be used later)
      allocate(tmp_c0_ld(ld%max_pole_num_ld,ld%num_ld),tmp_c1_ld(ld%max_pole_num_ld,ld%num_ld),&
               tmp_c2_ld(ld%max_pole_num_ld,ld%num_ld),tmp_c3_ld(ld%max_pole_num_ld,ld%num_ld),&
               tmp_c4_ld(                   ld%num_ld),tmp_c5_ld(                   ld%num_ld))
      tmp_c0_ld(:,:) = 0.0d0; tmp_c1_ld(:,:) = 0.0d0; tmp_c2_ld(:,:) = 0.0d0; tmp_c3_ld(:,:) = 0.0d0;
      tmp_c4_ld(  :) = 0.0d0; tmp_c5_ld(  :) = 0.0d0;
      call ld_coefficients(ld%media_ld,pole_num_ld,gamma_ld,omega_ld,f_ld,omega_p_ld,dt_em,pi, &
                           tmp_c0_ld,tmp_c1_ld,tmp_c2_ld,tmp_c3_ld,tmp_c4_ld,tmp_c5_ld,epsilon_em)
      if(allocated(diagonal_eps))then
        allocate(diagonal_c2(ld%max_pole_num_ld,ld%num_ld,3), &
          diagonal_c4(ld%num_ld,3),diagonal_c5(ld%num_ld,3))
        diagonal_c2=0d0;diagonal_c4=0d0;diagonal_c5=0d0
        diagonal_status=0
        do diagonal_axis=1,3
          do diagonal_media=1,ld%num_ld
            ii=ld%media_ld(diagonal_media)
            do ij=1,pole_num_ld(ii)
              diagonal_c2(ij,diagonal_media,diagonal_axis)=diagonal_f(diagonal_axis,ii,ij) &
                *omega_p_ld(ii)**2*dt_em/(8d0*pi*tmp_c0_ld(ij,diagonal_media))
              if(.not.ieee_is_finite(diagonal_c2(ij,diagonal_media,diagonal_axis)))diagonal_status=1
              diagonal_c4(diagonal_media,diagonal_axis)=diagonal_c4(diagonal_media,diagonal_axis) &
                +diagonal_c2(ij,diagonal_media,diagonal_axis)
              diagonal_c5(diagonal_media,diagonal_axis)=diagonal_c5(diagonal_media,diagonal_axis) &
                +diagonal_c2(ij,diagonal_media,diagonal_axis)
            enddo
            diagonal_c4(diagonal_media,diagonal_axis)=diagonal_eps(diagonal_axis,ii) &
              +2d0*pi*dt_em*diagonal_c4(diagonal_media,diagonal_axis)
            diagonal_c5(diagonal_media,diagonal_axis)=diagonal_eps(diagonal_axis,ii) &
              -2d0*pi*dt_em*diagonal_c5(diagonal_media,diagonal_axis)
            if(.not.ieee_is_finite(diagonal_c4(diagonal_media,diagonal_axis)))diagonal_status=1
            if(.not.ieee_is_finite(diagonal_c5(diagonal_media,diagonal_axis)))diagonal_status=1
            if(diagonal_c4(diagonal_media,diagonal_axis)<=0d0)diagonal_status=1
          enddo
        enddo
        call comm_bcast(diagonal_status,nproc_group_global)
        if(diagonal_status/=0)then
          if(present(status))then
            status=1111;return
          endif
          call stop_em('Nonfinite diagonal LD update coefficients.')
        endif
      endif
      ! LD owns its allocations and zero initialization, including halo values.
      if(present(status))then
#ifdef SALMONPY_RECOVERY_TESTING
        ld_failure_slot=0
        if(present(allocation_failure_slot))ld_failure_slot=max(0,-allocation_failure_slot)
        call ld_allocate(ld,fs%mg%is_array,fs%mg%ie_array,ld_allocation_status,failure_slot=ld_failure_slot)
#else
        call ld_allocate(ld,fs%mg%is_array,fs%mg%ie_array,ld_allocation_status)
#endif
        ! No communication/use of LD arrays before the material-vector gate.
        ! Carry a failure there instead of adding another startup agreement.
        if(ld_allocation_status/=0)allocation_status=1119
      else
      call ld_allocate(ld,fs%mg%is_array,fs%mg%ie_array)
      endif
    end if
    
    !*** set fdtd coeffient and light speed for media ID = 0 **************************************************!
    call eh_allocate_vector(fe%rep,0,media_num)
    call eh_allocate_vector(fe%rmu,0,media_num)
    call eh_allocate_vector(fe%sig,0,media_num)
    if(present(status))then
      call comm_get_max(allocation_status,nproc_group_global)
      if(allocation_status/=0)then
        status=allocation_status;return
      endif
    endif
    fe%rep(:)=1.0d0; fe%rmu(:)=1.0d0; fe%sig(:)=0.0d0;
    do ii=0,media_num
      fe%rep(ii)=epsilon_em(ii); fe%rmu(ii)=mu_em(ii); fe%sig(ii)=sigma_em(ii);
    end do
    do ii=0,media_num
      call eh_coeff
    end do
    call eh_coeff_h_geometry

    if(comm_is_root(nproc_id_global)) write(*,'(A)') 'SALMONPY_FDTD_INTERFACE_V2'
    


    !*** write media information ******************************************************************************!
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      write(*,'(A,I6)')          ' media_num   = ',media_num
      do ii=0,media_num
        write(*,*) "=========================="
        write(*,'(A,I3,A)')      ' id ',ii, ':'
        call eh_check_media_type(fe%rep(ii),fe%rmu(ii),fe%sig(ii),media_type(ii),tmp_name1)
        write(*,'(A,A)')         ' media_type  =  ', trim(tmp_name1)
        select case(media_type(ii))
        case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
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
    if(((fs%a_bc(1,1)=='pml').or.(fs%a_bc(1,1)=='pec')).and.(fs%mg%is(1)==fs%lg%is(1))) then !x, bottom
      fe%iey_x_is(1)=fs%mg%is(1)+1; fe%iez_x_is(1)=fs%mg%is(1)+1;
      if(fs%a_bc(1,1)=='pec') then
        fe%iey_z_is(1)=fs%mg%is(1)+1; fe%iez_y_is(1)=fs%mg%is(1)+1;
        fe%ihx_y_is(1)=fs%mg%is(1)+1; fe%ihx_z_is(1)=fs%mg%is(1)+1;
      end if
    end if
    if(((fs%a_bc(1,2)=='pml').or.(fs%a_bc(1,2)=='pec')).and.(fs%mg%ie(1)==fs%lg%ie(1))) then !x, top
      fe%iey_x_ie(1)=fs%mg%ie(1)-1; fe%iez_x_ie(1)=fs%mg%ie(1)-1;
      fe%iex_y_ie(1)=fs%mg%ie(1)-1; fe%iex_z_ie(1)=fs%mg%ie(1)-1;
      fe%ihy_z_ie(1)=fs%mg%ie(1)-1; fe%ihy_x_ie(1)=fs%mg%ie(1)-1;
      fe%ihz_x_ie(1)=fs%mg%ie(1)-1; fe%ihz_y_ie(1)=fs%mg%ie(1)-1;
      if(fs%a_bc(1,2)=='pec') then
        fe%iey_z_ie(1)=fs%mg%ie(1)-1; fe%iez_y_ie(1)=fs%mg%ie(1)-1;
        fe%ihx_y_ie(1)=fs%mg%ie(1)-1; fe%ihx_z_ie(1)=fs%mg%ie(1)-1;
      end if
    end if
    if(((fs%a_bc(2,1)=='pml').or.(fs%a_bc(2,1)=='pec')).and.(fs%mg%is(2)==fs%lg%is(2))) then !y, bottom
      fe%iex_y_is(2)=fs%mg%is(2)+1; fe%iez_y_is(2)=fs%mg%is(2)+1;
      if(fs%a_bc(2,1)=='pec') then
        fe%iex_z_is(2)=fs%mg%is(2)+1; fe%iez_x_is(2)=fs%mg%is(2)+1;
        fe%ihy_x_is(2)=fs%mg%is(2)+1; fe%ihy_z_is(2)=fs%mg%is(2)+1;
      end if
    end if
    if(((fs%a_bc(2,2)=='pml').or.(fs%a_bc(2,2)=='pec')).and.(fs%mg%ie(2)==fs%lg%ie(2))) then !y, top
      fe%iex_y_ie(2)=fs%mg%ie(2)-1; fe%iez_y_ie(2)=fs%mg%ie(2)-1;
      fe%iey_z_ie(2)=fs%mg%ie(2)-1; fe%iey_x_ie(2)=fs%mg%ie(2)-1;
      fe%ihx_y_ie(2)=fs%mg%ie(2)-1; fe%ihx_z_ie(2)=fs%mg%ie(2)-1;
      fe%ihz_x_ie(2)=fs%mg%ie(2)-1; fe%ihz_y_ie(2)=fs%mg%ie(2)-1;
      if(fs%a_bc(2,2)=='pec') then
        fe%iex_z_ie(2)=fs%mg%ie(2)-1; fe%iez_x_ie(2)=fs%mg%ie(2)-1;
        fe%ihy_x_ie(2)=fs%mg%ie(2)-1; fe%ihy_z_ie(2)=fs%mg%ie(2)-1;
      end if
    end if
    if(((fs%a_bc(3,1)=='pml').or.(fs%a_bc(3,1)=='pec')).and.(fs%mg%is(3)==fs%lg%is(3))) then !z, bottom
      fe%iex_z_is(3)=fs%mg%is(3)+1; fe%iey_z_is(3)=fs%mg%is(3)+1;
      if(fs%a_bc(3,1)=='pec') then
        fe%iex_y_is(3)=fs%mg%is(3)+1; fe%iey_x_is(3)=fs%mg%is(3)+1;
        fe%ihz_x_is(3)=fs%mg%is(3)+1; fe%ihz_y_is(3)=fs%mg%is(3)+1;
      end if
    end if
    if(((fs%a_bc(3,2)=='pml').or.(fs%a_bc(3,2)=='pec')).and.(fs%mg%ie(3)==fs%lg%ie(3))) then !z, top
      fe%iex_z_ie(3)=fs%mg%ie(3)-1; fe%iey_z_ie(3)=fs%mg%ie(3)-1;
      fe%iez_x_ie(3)=fs%mg%ie(3)-1; fe%iez_y_ie(3)=fs%mg%ie(3)-1;
      fe%ihx_y_ie(3)=fs%mg%ie(3)-1; fe%ihx_z_ie(3)=fs%mg%ie(3)-1;
      fe%ihy_z_ie(3)=fs%mg%ie(3)-1; fe%ihy_x_ie(3)=fs%mg%ie(3)-1;
      if(fs%a_bc(3,2)=='pec') then
        fe%iex_y_ie(3)=fs%mg%ie(3)-1; fe%iey_x_ie(3)=fs%mg%ie(3)-1;
        fe%ihz_x_ie(3)=fs%mg%ie(3)-1; fe%ihz_y_ie(3)=fs%mg%ie(3)-1;
      end if
    end if
    
    !*** set pml **********************************************************************************************!
    call eh_set_pml(1,fe%c1_ey_x,fe%c2_ey_x,fe%c1_ez_x,fe%c2_ez_x,&
                      fe%c1_hy_x,fe%c2_hy_x,fe%c1_hz_x,fe%c2_hz_x) !x direction
    call eh_set_pml(2,fe%c1_ez_y,fe%c2_ez_y,fe%c1_ex_y,fe%c2_ex_y,&
                      fe%c1_hz_y,fe%c2_hz_y,fe%c1_hx_y,fe%c2_hx_y) !y direction
    call eh_set_pml(3,fe%c1_ex_z,fe%c2_ex_z,fe%c1_ey_z,fe%c2_ey_z,&
                      fe%c1_hx_z,fe%c2_hx_z,fe%c1_hy_z,fe%c2_hy_z) !z direction
    if(fe%flag_pml) then
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
    call eh_allocate_vector(fe%c_ebloch_x,fs%mg%is_array(1),fs%mg%ie_array(1))
    call eh_allocate_vector(fe%c_ebloch_y,fs%mg%is_array(2),fs%mg%ie_array(2))
    call eh_allocate_vector(fe%c_ebloch_z,fs%mg%is_array(3),fs%mg%ie_array(3))
    call eh_allocate_vector(fe%c_hbloch_x,fs%mg%is_array(1),fs%mg%ie_array(1))
    call eh_allocate_vector(fe%c_hbloch_y,fs%mg%is_array(2),fs%mg%ie_array(2))
    call eh_allocate_vector(fe%c_hbloch_z,fs%mg%is_array(3),fs%mg%ie_array(3))
    if(present(status))then
      call comm_get_max(allocation_status,nproc_group_global)
      if(allocation_status/=0)then
        status=1120;return
      endif
    endif
    fe%c_ebloch_x(:) = 1.0d0; fe%c_ebloch_y(:) = 1.0d0; fe%c_ebloch_z(:) = 1.0d0;
    fe%c_hbloch_x(:) = 1.0d0; fe%c_hbloch_y(:) = 1.0d0; fe%c_hbloch_z(:) = 1.0d0;
    if( (sum(abs(bloch_k_em(:)))/=0.0d0)    &
        .or.(bloch_real_imag_em(1)=='imag') &
        .or.(bloch_real_imag_em(2)=='imag') &
        .or.(bloch_real_imag_em(3)=='imag') ) then !--- use bloch boundary condition --------------------------!
      !check yn_periodic
      if(yn_periodic/='y') then
        call stop_em(   'When bloch boundary conditions for electromagnetic problem are used by |bloch_k_em| > 0,',&
                     c2='yn_periodic must be y.')
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
    if(fe%flag_obs) then
      !set initial
      allocate(fe%iobs_po_id(obs_num_em,3)) !1:x,        2:y,        3:z
      allocate(fe%iobs_po_pe(obs_num_em))
      allocate(fe%iobs_li_pe(obs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
      allocate(fe%iobs_pl_pe(obs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
      fe%iobs_po_id(:,:)=0; fe%iobs_po_pe(:)=0; fe%iobs_li_pe(:,:)=0; fe%iobs_pl_pe(:,:)=0; 
      fe%e_max=0.0d0; fe%h_max=0.0d0;
      
      !search observation point
      do ii=1,obs_num_em
        call find_point_line_plane_em(obs_loc_em(ii,:),fe%iobs_po_id(ii,:),                     &
                                      fe%iobs_po_pe(ii),fe%iobs_li_pe(ii,:),fe%iobs_pl_pe(ii,:),&
                                      fs%mg%is(:),fs%mg%ie(:),                                  &
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
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jx_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jx_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jx_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jy_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jy_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jy_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jz_xy_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jz_yz_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),&
                                          obs_num_em,maxval(fe%iobs_num_ene(:)),2),&
                  fe%obs_jz_xz_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),&
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
        fe%obs_jx_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jx_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jx_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jy_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jy_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jy_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jz_xy_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jz_yz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
        fe%obs_jz_xz_ene(:,:,:,:,:)=(0.0d0,0.0d0);
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
        !initialize
        iroot1 = 0; iroot2 = 0;
        
        !find pe
        if(fe%iobs_po_pe(ii)==1) then
          iroot1 = nproc_id_global
        end if
        
        !determine pe
        call comm_summation(iroot1,iroot2,nproc_group_global)
        
        !coordinate
        if(comm_is_root(nproc_id_global)) then
          write(*,'(I3,A,3ES14.5)') ii,":",(fe%coo(fe%iobs_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        end if
        
        !ex
        if(nproc_id_global==iroot2) then
          ij = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)  )
          ik = fs%imedia(fe%iobs_po_id(ii,1)-1,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)  )
        end if
        call comm_bcast(ij,nproc_group_global,iroot2)
        call comm_bcast(ik,nproc_group_global,iroot2)
        call eh_check_media_type(fe%rep(ij),fe%rmu(ij),fe%sig(ij),media_type(ij),tmp_name1)
        call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
        if(comm_is_root(nproc_id_global)) then
          write(*,'(6X,A)') "Ex is averaged from two spatial grids whose parameters are:"
          write(*,'(9X,A,I6,A,I6)')         'media ID    = ', ij              , ', ',ik
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
        end if
        
        !ey
        if(nproc_id_global==iroot2) then
          ik = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)-1,fe%iobs_po_id(ii,3)  )
        end if
        call comm_bcast(ik,nproc_group_global,iroot2)
        call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
        if(comm_is_root(nproc_id_global)) then
          write(*,'(6X,A)') "Ey is averaged from two spatial grids whose parameters are:"
          write(*,'(9X,A,I6,A,I6)')         'media ID    = ', ij              , ', ',ik
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
        end if
        
        !ez
        if(nproc_id_global==iroot2) then
          ik = fs%imedia(fe%iobs_po_id(ii,1)  ,fe%iobs_po_id(ii,2)  ,fe%iobs_po_id(ii,3)-1)
        end if
        call comm_bcast(ik,nproc_group_global,iroot2)
        call eh_check_media_type(fe%rep(ik),fe%rmu(ik),fe%sig(ik),media_type(ik),tmp_name2)
        if(comm_is_root(nproc_id_global)) then
          write(*,'(6X,A)') "Ez is averaged from two spatial grids whose parameters are:"
          write(*,'(9X,A,I6,A,I6)')         'media ID    = ', ij              , ', ',ik
          write(*,'(9X,A,A,A,A)')           'media_type  =  ', trim(tmp_name1), ', ',trim(tmp_name2)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'epsilon_em  = ', fe%rep(ij)      , ', ',fe%rep(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'mu_em       = ', fe%rmu(ij)      , ', ',fe%rmu(ik)
          write(*,'(9X,A,ES12.5,A,ES12.5)') 'sigma_em    = ', fe%sig(ij)      , ', ',fe%sig(ik)
        end if
      end do
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "**************************"
        do ii=1,obs_num_em
          !set header for point data
          write(save_name,*) ii
          save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//'_at_point_rt.data'
          call eh_point_output(fe,save_name,.true.)
          
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
              call eh_integral_output(save_name,.true.,warned=fe%plane_integral_io_warned)
            end do
          end if
        end do
      end if
    end if
    
    !*** check and prepare incident current source ************************************************************!
    !check condition
    select case(wave_input)
    case('source')
      !check linear response
      if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
        call stop_em(   'Invalid ae_shape1/2:',&
                     c2='For ae_shape1/2 = impulse, wave_input must be default(do not set source).')
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
        call stop_em(   'Invalid wave_input:',      &
                     c2='wave_input must be source',&
                     c3='or ae_shape1 and/or ae_shape2 must be impulse.')
      end if
    end select
    
    !prepare incident current source
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
        call find_point_line_plane_em(source_loc1(:),fe%inc_po_id(ii,:),                     &
                                      fe%inc_po_pe(ii),fe%inc_li_pe(ii,:),fe%inc_pl_pe(ii,:),&
                                      fs%mg%is(:),fs%mg%ie(:),                               &
                                      minval(fs%lg%is(:))-fe%Nd,maxval(fs%lg%ie(:))+fe%Nd,fe%coo(:,:))
        select case(wave_input)
        case('point_ip','x-line_ip','y-line_ip','z-line_ip','xy-plane_ip','yz-plane_ip','xz-plane_ip')
        case default
          call eh_check_iw_parameter(ii,phi_cep1,I_wcm2_1,E_amplitude1,ae_shape1)
        end select
      end if
      if(fe%inc_dist2/='none') then
        ii=2
        call find_point_line_plane_em(source_loc2(:),fe%inc_po_id(ii,:),                     &
                                      fe%inc_po_pe(ii),fe%inc_li_pe(ii,:),fe%inc_pl_pe(ii,:),&
                                      fs%mg%is(:),fs%mg%ie(:),                               &
                                      minval(fs%lg%is(:))-fe%Nd,maxval(fs%lg%ie(:))+fe%Nd,fe%coo(:,:))
        select case(wave_input)
        case('point_ip','x-line_ip','y-line_ip','z-line_ip','xy-plane_ip','yz-plane_ip','xz-plane_ip')
        case default
          call eh_check_iw_parameter(ii,phi_cep2,I_wcm2_2,E_amplitude2,ae_shape2)
        end select
      end if
      
      !set gauss function used for beam spot
      allocate( fe%gbeam_width_xy1(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(2):fs%mg%ie_array(2)), &
                fe%gbeam_width_xy2(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(2):fs%mg%ie_array(2)), &
                fe%gbeam_width_yz1(fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3)), &
                fe%gbeam_width_yz2(fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3)), &
                fe%gbeam_width_xz1(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(3):fs%mg%ie_array(3)), &
                fe%gbeam_width_xz2(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(3):fs%mg%ie_array(3)) )
      allocate( fe%gbeam_width_x1( fs%mg%is_array(1):fs%mg%ie_array(1) ), &
                fe%gbeam_width_x2( fs%mg%is_array(1):fs%mg%ie_array(1) ), &
                fe%gbeam_width_y1( fs%mg%is_array(2):fs%mg%ie_array(2) ), &
                fe%gbeam_width_y2( fs%mg%is_array(2):fs%mg%ie_array(2) ), &
                fe%gbeam_width_z1( fs%mg%is_array(3):fs%mg%ie_array(3) ), &
                fe%gbeam_width_z2( fs%mg%is_array(3):fs%mg%ie_array(3) ) )
      fe%gbeam_width_xy1(:,:) = 1.0d0; fe%gbeam_width_xy2(:,:) = 1.0d0;
      fe%gbeam_width_yz1(:,:) = 1.0d0; fe%gbeam_width_yz2(:,:) = 1.0d0;
      fe%gbeam_width_xz1(:,:) = 1.0d0; fe%gbeam_width_xz2(:,:) = 1.0d0;
      fe%gbeam_width_x1(:)    = 1.0d0; fe%gbeam_width_x2(:)    = 1.0d0; 
      fe%gbeam_width_y1(:)    = 1.0d0; fe%gbeam_width_y2(:)    = 1.0d0; 
      fe%gbeam_width_z1(:)    = 1.0d0; fe%gbeam_width_z2(:)    = 1.0d0; 
      call eh_make_gauss(source_loc1(:),gbeam_sigma_plane1,gbeam_sigma_line1,      &
                         fe%gbeam_width_xy1,fe%gbeam_width_yz1,fe%gbeam_width_xz1, &
                         fe%gbeam_width_x1 ,fe%gbeam_width_y1 ,fe%gbeam_width_z1  )
      call eh_make_gauss(source_loc2(:),gbeam_sigma_plane2,gbeam_sigma_line2,      &
                         fe%gbeam_width_xy2,fe%gbeam_width_yz2,fe%gbeam_width_xz2, &
                         fe%gbeam_width_x2 ,fe%gbeam_width_y2 ,fe%gbeam_width_z2  )
      
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
    
    if(subpixel_profile%nm>0)then
      call eh_subpixel_coefficients(diagonal_status)
      if(diagonal_status/=0)then
        if(present(status))then
          status=1114;return
        endif
        call stop_em('Invalid subpixel material update coefficients.')
      endif
    else
      if(fe%tensor_enabled)call eh_tensor_weights
    endif
    if(fe%tensor_enabled)then
      if(any(bloch_k_em/=0d0))then
        if(present(status))then
          status=1112;return
        endif
        call stop_em('Nonzero Bloch phase with off-diagonal tensor is not implemented.')
      endif
    endif
    !*** prepare linear response ******************************************************************************!
    if(fe%flag_lr) then
      !check condition
      if(yn_periodic=='y'.and.trans_longi/='tr') flag_stop = .true.
      do ii=0,media_num
        if(fe%rmu(ii)/=1.0d0.or.fe%sig(ii)/=0.0d0) flag_stop = .true.
        if(ii==0.and.fe%rep(ii)/=1.0d0)flag_stop=.true.
        if(ii==0) then
          select case(media_type(ii))
          case('vacuum')
            continue
          case default
            flag_stop = .true.
          end select
        else
          select case(media_type(ii))
          case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
            continue
          case default
            flag_stop = .true.
          end select
        end if
      end do
      if(flag_stop) then
        call stop_em(   'Invalid input keywords:',&
                     c2='When ae_shape1=impulse and/or ae_shape2=impulse,',&
                     c3='Vacuum epsilon_em and all mu_em must be 1.0d0.',&
                     c4='sigma_em must be 0.0d0.',&
                     c5='media_type(i) must be lorentz-drude, where i > 0.',&
                     c6='trans_longi must be tr for yn_periodic=y.')
      end if
      
      !set impulse
      if(ae_shape1=='impulse') then
        tmp_x1 = e_impulse*(epdir_re1(1)+epdir_im1(1))
        tmp_y1 = e_impulse*(epdir_re1(2)+epdir_im1(2))
        tmp_z1 = e_impulse*(epdir_re1(3)+epdir_im1(3))
      else
        tmp_x1 = 0.0d0
        tmp_y1 = 0.0d0
        tmp_z1 = 0.0d0
      end if
      if(ae_shape2=='impulse') then
        tmp_x2 = e_impulse*(epdir_re2(1)+epdir_im2(1))
        tmp_y2 = e_impulse*(epdir_re2(2)+epdir_im2(2))
        tmp_z2 = e_impulse*(epdir_re2(3)+epdir_im2(3))
      else
        tmp_x2 = 0.0d0
        tmp_y2 = 0.0d0
        tmp_z2 = 0.0d0
      end if
      
      if(fe%subpixel_enabled)then
        call eh_subpixel_impulse
      else
      !set initial current density
      if(ld%enabled) then
        do ii = 1,ld%num_ld
        do ij = 1,pole_num_ld(ld%media_ld(ii))
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          ! Only the edge's owning LD material contributes to this pole.
          if(ld%c1_jx_ld(ix,iy,iz,ij)>0.0d0.and.fs%imedia(ix,iy,iz)==ld%media_ld(ii)) then
            ld%rjx_ld(ix,iy,iz,ij)   = ld%rjx_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(1,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_x1
            ld%rjx_ld(ix,iy,iz,ij)   = ld%rjx_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(1,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_x2
            ld%rjx_fdtd_ld(ix,iy,iz) = ld%rjx_fdtd_ld(ix,iy,iz)                                          &
                                      +0.5d0*( (1.0d0+ld%c1_jx_ld(ix,iy,iz,ij)) * ld%rjx_ld(ix,iy,iz,ij) &
                                                     -ld%c3_jx_ld(ix,iy,iz,ij)  * ld%px_ld(ix,iy,iz,ij) )
          end if
          ! Only the edge's owning LD material contributes to this pole.
          if(ld%c1_jy_ld(ix,iy,iz,ij)>0.0d0.and.fs%imedia(ix,iy,iz)==ld%media_ld(ii)) then
            ld%rjy_ld(ix,iy,iz,ij)   = ld%rjy_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(2,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_y1
            ld%rjy_ld(ix,iy,iz,ij)   = ld%rjy_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(2,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_y2
            ld%rjy_fdtd_ld(ix,iy,iz) = ld%rjy_fdtd_ld(ix,iy,iz)                                          &
                                      +0.5d0*( (1.0d0+ld%c1_jy_ld(ix,iy,iz,ij)) * ld%rjy_ld(ix,iy,iz,ij) &
                                                     -ld%c3_jy_ld(ix,iy,iz,ij)  * ld%py_ld(ix,iy,iz,ij) )
          end if
          ! Only the edge's owning LD material contributes to this pole.
          if(ld%c1_jz_ld(ix,iy,iz,ij)>0.0d0.and.fs%imedia(ix,iy,iz)==ld%media_ld(ii)) then
            ld%rjz_ld(ix,iy,iz,ij)   = ld%rjz_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(3,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_z1
            ld%rjz_ld(ix,iy,iz,ij)   = ld%rjz_ld(ix,iy,iz,ij)                                        &
                                      -(eh_axis_strength(3,ld%media_ld(ii),ij)*omega_p_ld(ld%media_ld(ii))**2.0d0) &
                                      /(4.0d0*pi)*tmp_z2
            ld%rjz_fdtd_ld(ix,iy,iz) = ld%rjz_fdtd_ld(ix,iy,iz)                                          &
                                      +0.5d0*( (1.0d0+ld%c1_jz_ld(ix,iy,iz,ij)) * ld%rjz_ld(ix,iy,iz,ij) &
                                                     -ld%c3_jz_ld(ix,iy,iz,ij)  * ld%pz_ld(ix,iy,iz,ij) )
          end if
        end do
        end do
        end do
        end do
        end do
      end if
      
      if(fe%tensor_enabled)call eh_tensor_impulse

      endif
      !initialize and allocate
      allocate(fe%time_lr(nt_em))
      fe%time_lr(:)=0.0d0
      allocate(fe%fr_lr(nenergy,3),fe%fi_lr(nenergy,3))
      fe%fr_lr(:,:)=0.0d0; fe%fi_lr(:,:)=0.0d0;
      if(yn_periodic=='n') then
        allocate(fe%dip_lr(nt_em,3))
        fe%dip_lr(:,:)=0.0d0
      elseif(yn_periodic=='y') then
        allocate(fe%curr_lr(nt_em,3),fe%e_lr(nt_em,3))
        fe%curr_lr(:,:)=0.0d0; fe%e_lr(:,:)=0.0d0;
        allocate(fe%er_lr(0:nenergy,3),fe%ei_lr(0:nenergy,3))
        fe%er_lr(:,:)=0.0d0; fe%ei_lr(:,:)=0.0d0;
      end if
    end if
    
    !*** prepare ase option ***********************************************************************************!
    if(fe%flag_ase) then
      !check condition
      if(sum(ase_box_size_em(:))<=0.0d0 ) call stop_em('When ase_num_em > 0, ase_box_size_em must be > 0.')
      if(sum(ek_dir2(:))         >0.0d0 ) call stop_em('When ase_num_em > 0, only pulse1 can be used.'    )
      if(trim(ae_shape2)        /='none') call stop_em('When ase_num_em > 0, only pulse1 can be used.'    )
      if(yn_periodic            =='y'   ) call stop_em('When ase_num_em > 0, yn_periodic must be n.'      )
      if(ase_smedia_id_em/=media_id_source1) then
        call stop_em('When ase_num_em > 0, ase_smedia_id_em and media_id_source1 must be same.')
      end if
      if(fe%sig(ase_smedia_id_em)>0.0d0) then
        call stop_em('When ase_num_em > 0, sigma_em(ase_smedia_id_em) must be 0.0d0.'  )
      end if
      if(ae_shape1=='impulse') then
        call stop_em('When ase_num_em > 0, impulse can not be used.')
      end if
      if( (gbeam_sigma_plane1(1)>0.0d0).or.(gbeam_sigma_plane1(2)>0.0d0).or.(gbeam_sigma_plane1(3)>0.0d0).or. &
          (gbeam_sigma_plane2(1)>0.0d0).or.(gbeam_sigma_plane2(2)>0.0d0).or.(gbeam_sigma_plane2(3)>0.0d0).or. &
          (gbeam_sigma_line1(1) >0.0d0).or.(gbeam_sigma_line1(2) >0.0d0).or.(gbeam_sigma_line1(3) >0.0d0).or. &
          (gbeam_sigma_line2(1) >0.0d0).or.(gbeam_sigma_line2(2) >0.0d0).or.(gbeam_sigma_line2(3) >0.0d0) ) then
        call stop_em(   'When ase_num_em > 0,',      &
                     c2='gauss beam cannot be used.',&
                     c3='gbeam_sigma_plane1/2 and gbeam_sigma_line1/2 must be <= 0.0d0')
      end if
      if    ((ase_ene_min_em <0.0d0).and.(ase_ene_max_em <0.0d0).and.&
             (ase_wav_min_em <0.0d0).and.(ase_wav_max_em <0.0d0))         then
        flag_stop = .true.
      elseif((ase_ene_min_em>=0.0d0).and.(ase_ene_max_em>=0.0d0).and.&
             (ase_wav_min_em>=0.0d0).and.(ase_wav_max_em>=0.0d0))         then
        flag_stop = .true.
      elseif((ase_ene_min_em>=0.0d0).and.(ase_wav_min_em>=0.0d0))         then
        flag_stop = .true.
      elseif((ase_ene_min_em>=0.0d0).and.(ase_wav_max_em>=0.0d0)) then
        flag_stop = .true.
      elseif((ase_wav_min_em>=0.0d0).and.(ase_ene_max_em>=0.0d0)) then
        flag_stop = .true.
      elseif((ase_ene_min_em>=0.0d0).and.(ase_ene_max_em<ase_ene_min_em)) then
        flag_stop = .true.
      elseif((ase_wav_min_em>=0.0d0).and.(ase_wav_max_em<ase_wav_min_em)) then
        flag_stop = .true.
      end if
      if(flag_stop) then
        call stop_em(   'When ase_num_em > 0, you must set either',     &
                     c2='[ase_ene_min_em and ase_ene_min_em >= 0.0d0]', &
                     c3='or',                                           &
                     c4='[ase_wav_min_em and ase_wav_min_em >= 0.0d0].',&
                     c5='ase_ene/wav_max_em must be larger than ase_ene/wav_min_em.')
      end if
      if    (ek_dir1(1)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,1),1) >= (-fs%rlsize(1)/4.0d0) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='source_loc1(1) must be set < -al_em(1)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,1),1) >= (ase_box_cent_em(1)-(ase_box_size_em(1)/2.0d0)) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='source_loc1(1) must be set < ase_box_cent_em(1)-(ase_box_size_em(1)/2.0d0).')
        end if
      elseif(ek_dir1(2)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,2),2) >= (-fs%rlsize(2)/4.0d0) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to y-axis,', &
                       c2='source_loc1(2) must be set < -al_em(2)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,2),2) >= (ase_box_cent_em(2)-(ase_box_size_em(2)/2.0d0)) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to y-axis,', &
                       c2='source_loc1(2) must be set < ase_box_cent_em(2)-(ase_box_size_em(2)/2.0d0).')
        end if
      elseif(ek_dir1(3)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,3),3) >= (-fs%rlsize(3)/4.0d0) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to z-axis,', &
                       c2='source_loc1(3) must be set < -al_em(3)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,3),3) >= (ase_box_cent_em(3)-(ase_box_size_em(3)/2.0d0)) ) then
          call stop_em(   'When ase_num_em > 0 with incident pulse propagating to z-axis,', &
                       c2='source_loc1(3) must be set < ase_box_cent_em(3)-(ase_box_size_em(3)/2.0d0).')
        end if
      end if
      
      !set light speed for surrounding media
      fe%cspeed_smedia_ase = cspeed_au / sqrt(fe%rep(ase_smedia_id_em)*fe%rmu(ase_smedia_id_em))
      
      !make sampling energy axis
      if    (ase_ene_min_em>=0.0d0) then
        tmp_min = ase_ene_min_em; tmp_max = ase_ene_max_em;
      elseif(ase_wav_min_em>=0.0d0) then
        tmp_min = ase_wav_min_em; tmp_max = ase_wav_max_em;
      end if
      allocate(fe%ase_ene(ase_num_em)); fe%ase_ene(:)=0.0d0;
      do ii=1,ase_num_em
        if (ase_num_em==1) then
          fe%ase_ene(ii) = tmp_min
        else
          fe%ase_ene(ii) = tmp_min + ( (tmp_max-tmp_min)/(ase_num_em-1) )*(dble(ii)-1)
        end if
        if(ase_wav_min_em>=0.0d0) then
          !convert from wavelength to energy
          fe%ase_ene(ii) = 2.0d0*pi*cspeed_au/fe%ase_ene(ii)
        end if
      end do
      
      !specify propagation axis and allocate incident ex-z and hx-z
      if    (ek_dir1(1)==1.0d0) then
        allocate( fe%ase_ex_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2),&
                  fe%ase_ey_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2),&
                  fe%ase_ez_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2),&
                  fe%ase_hx_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2),&
                  fe%ase_hy_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2),&
                  fe%ase_hz_inc_pa_ene(fs%mg%is(1):fs%mg%ie(1),ase_num_em,2) )
      elseif(ek_dir1(2)==1.0d0) then
        allocate( fe%ase_ex_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2),&
                  fe%ase_ey_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2),&
                  fe%ase_ez_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2),&
                  fe%ase_hx_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2),&
                  fe%ase_hy_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2),&
                  fe%ase_hz_inc_pa_ene(fs%mg%is(2):fs%mg%ie(2),ase_num_em,2) )
      elseif(ek_dir1(3)==1.0d0) then
        allocate( fe%ase_ex_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2),&
                  fe%ase_ey_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2),&
                  fe%ase_ez_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2),&
                  fe%ase_hx_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2),&
                  fe%ase_hy_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2),&
                  fe%ase_hz_inc_pa_ene(fs%mg%is(3):fs%mg%ie(3),ase_num_em,2) )
      end if
      fe%ase_ex_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      fe%ase_ey_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      fe%ase_ez_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      fe%ase_hx_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      fe%ase_hy_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      fe%ase_hz_inc_pa_ene(:,:,:)=(0.0d0,0.0d0);
      
      !allocation
      allocate( fe%ase_ex_inc_bt_ene(ase_num_em,2,2),&
                fe%ase_ey_inc_bt_ene(ase_num_em,2,2),&
                fe%ase_ez_inc_bt_ene(ase_num_em,2,2),&
                fe%ase_hx_inc_bt_ene(ase_num_em,2,2),&
                fe%ase_hy_inc_bt_ene(ase_num_em,2,2),&
                fe%ase_hz_inc_bt_ene(ase_num_em,2,2) )
      fe%ase_ex_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%ase_ey_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%ase_ez_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%ase_hx_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%ase_hy_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%ase_hz_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      allocate( fe%ase_ex_xy_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),ase_num_em,2,2),&
                fe%ase_ey_xy_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),ase_num_em,2,2),&
                fe%ase_hx_xy_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),ase_num_em,2,2),&
                fe%ase_hy_xy_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),ase_num_em,2,2) )
      fe%ase_ex_xy_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_ey_xy_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hx_xy_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hy_xy_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      allocate( fe%ase_ey_yz_sca_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_ez_yz_sca_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_hy_yz_sca_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_hz_yz_sca_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2) )
      fe%ase_ey_yz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_ez_yz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hy_yz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hz_yz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      allocate( fe%ase_ex_xz_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_ez_xz_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_hx_xz_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2),&
                fe%ase_hz_xz_sca_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),ase_num_em,2,2) )
      fe%ase_ex_xz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_ez_xz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hx_xz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      fe%ase_hz_xz_sca_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      allocate( fe%iase_mask_xy(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2)),&
                fe%iase_mask_yz(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3)),&
                fe%iase_mask_xz(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3)) )
      fe%iase_mask_xy(:,:) = 0; fe%iase_mask_yz(:,:) = 0; fe%iase_mask_xz(:,:) = 0;
      
      !set mask function
!$omp parallel
!$omp do private(ix,iy)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        if( (abs(fe%coo(ix,1)-ase_box_cent_em(1))<=(ase_box_size_em(1)/2.0d0)).and.&
            (abs(fe%coo(iy,2)-ase_box_cent_em(2))<=(ase_box_size_em(2)/2.0d0)) ) then
          fe%iase_mask_xy(ix,iy) = 1;
        end if
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
        if( (abs(fe%coo(iy,2)-ase_box_cent_em(2))<=(ase_box_size_em(2)/2.0d0)).and.&
            (abs(fe%coo(iz,3)-ase_box_cent_em(3))<=(ase_box_size_em(3)/2.0d0)) ) then
          fe%iase_mask_yz(iy,iz) = 1;
        end if
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        if( (abs(fe%coo(ix,1)-ase_box_cent_em(1))<=(ase_box_size_em(1)/2.0d0)).and.&
            (abs(fe%coo(iz,3)-ase_box_cent_em(3))<=(ase_box_size_em(3)/2.0d0)) ) then
          fe%iase_mask_xz(ix,iz) = 1;
        end if
      end do
      end do
!$omp end do
!$omp end parallel
      
      !set id on closed surface[box shape]
      !--- ix on yz plane ----!
      call find_point_id_only_em(fs%lg,itmp1,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1)-(ase_box_size_em(1)/2.0d0),   &
                                    ase_box_cent_em(2),                              &
                                    ase_box_cent_em(3) /) )
      call find_point_id_only_em(fs%lg,itmp2,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1)+(ase_box_size_em(1)/2.0d0),   &
                                    ase_box_cent_em(2),                              &
                                    ase_box_cent_em(3) /) )
      fe%id_on_box_surf(1,1) = itmp1(1) !bottom
      fe%id_on_box_surf(1,2) = itmp2(1) !top
      !--- iy on xz plane ----!
      call find_point_id_only_em(fs%lg,itmp1,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1),                              &
                                    ase_box_cent_em(2)-(ase_box_size_em(2)/2.0d0),   &
                                    ase_box_cent_em(3) /) )
      call find_point_id_only_em(fs%lg,itmp2,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1),                              &
                                    ase_box_cent_em(2)+(ase_box_size_em(2)/2.0d0),   &
                                    ase_box_cent_em(3) /) )
      fe%id_on_box_surf(2,1) = itmp1(2) !bottom
      fe%id_on_box_surf(2,2) = itmp2(2) !top
      !--- iz on xy plane ----!
      call find_point_id_only_em(fs%lg,itmp1,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1),                              &
                                    ase_box_cent_em(2),                              &
                                    ase_box_cent_em(3)-(ase_box_size_em(3)/2.0d0) /) )
      call find_point_id_only_em(fs%lg,itmp2,fe%Nd,fe%coo(:,:),                      &
                                 (/ ase_box_cent_em(1),                              &
                                    ase_box_cent_em(2),                              &
                                    ase_box_cent_em(3)+(ase_box_size_em(3)/2.0d0) /) )
      fe%id_on_box_surf(3,1) = itmp1(3) !bottom
      fe%id_on_box_surf(3,2) = itmp2(3) !top
      
      !set ipe on closed surface[box shape]
      do ii=1,3
        do ibt = 1,2
          if((fs%mg%is(ii)<=fe%id_on_box_surf(ii,ibt)).and.(fe%id_on_box_surf(ii,ibt)<=fs%mg%ie(ii))) then
            fe%ipe_on_box_surf(ii,ibt)=1
          else
            fe%ipe_on_box_surf(ii,ibt)=0
          end if
        end do
      end do
      
      !check media ID on closed surface[box shape](that must be 0)
      ii=0; ij=0;
!$omp parallel
!$omp do private(ix,iy,ik)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        if(fe%iase_mask_xy(ix,iy)==1) then
          do ik=1,2
            if((fs%mg%is(3)<=fe%id_on_box_surf(3,ik)).and.fe%id_on_box_surf(3,ik)<=fs%mg%ie(3)) then
              if(fs%imedia(ix,iy,fe%id_on_box_surf(3,ik))/=ase_smedia_id_em) ii=1
            end if
          end do
        end if
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz,ik)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
        if(fe%iase_mask_yz(iy,iz)==1) then
          do ik=1,2
            if((fs%mg%is(1)<=fe%id_on_box_surf(1,ik)).and.fe%id_on_box_surf(1,ik)<=fs%mg%ie(1)) then
              if(fs%imedia(fe%id_on_box_surf(1,ik),iy,iz)/=ase_smedia_id_em) ii=1
            end if
          end do
        end if
      end do
      end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz,ik)
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        if(fe%iase_mask_xz(ix,iz)==1) then
          do ik=1,2
            if((fs%mg%is(2)<=fe%id_on_box_surf(2,ik)).and.fe%id_on_box_surf(2,ik)<=fs%mg%ie(2)) then
              if(fs%imedia(ix,fe%id_on_box_surf(2,ik),iz)/=ase_smedia_id_em) ii=1
            end if
          end do
        end if
      end do
      end do
!$omp end do
!$omp end parallel
      call comm_summation(ii,ij,nproc_group_global)
      if(ij>0) then
        call stop_em(   'When ase_num_em > 0, media ID on closed surface[box shape] must be ase_smedia_id_em.', &
                     c2='ase_box_cent_em and/or ase_box_size_em must be reconsidered.')
      end if
      
      !write information
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        write(*,*) "ase_num_em > 0:"
        write(*,*) "Absorption-, Scattering-, and Extinction-cross-sections"
        write(*,*) "will be outputed at the end of calculation."
        write(*,*) "NOTE: Those are normalized by the spectral distribution of the incident pulse."
        write(*,*)
        write(*,*) "Closed surface(box shape) ranges from"
        write(*,'(A,ES12.5,A,ES12.5)') " x = ",fe%coo(fe%id_on_box_surf(1,1),1)*ulength_from_au,&
                                       " to " ,fe%coo(fe%id_on_box_surf(1,2),1)*ulength_from_au
        write(*,'(A,ES12.5,A,ES12.5)') " y = ",fe%coo(fe%id_on_box_surf(2,1),2)*ulength_from_au,&
                                       " to " ,fe%coo(fe%id_on_box_surf(2,2),2)*ulength_from_au
        write(*,'(A,ES12.5,A,ES12.5)') " z = ",fe%coo(fe%id_on_box_surf(3,1),3)*ulength_from_au,&
                                       " to " ,fe%coo(fe%id_on_box_surf(3,2),3)*ulength_from_au
        write(*,*) "in the unit system, ",trim(unit_system),"."
        if(obs_samp_em==1) then
          write(*,*)
          write(*, '(A)') " --- CAUTION ---"
          write(*, '(A)') " In this calculation, obs_samp_em is set to 1."
          write(*,'(2A)') " It is strongly recommended to increase obs_samp_em",&
                          " as far as your intending sampling period allows"
          write(*, '(A)') " because the calculation cost of this ase-option can be reduced with larger obs_samp_em."
        end if
        write(*,*) "**************************"
      end if
    end if
    
    !*** prepare art option ***********************************************************************************!
    if(fe%flag_art) then
      !check condition
      if(sum(ek_dir2(:))  >0.0d0 ) call stop_em('When art_num_em > 0, only pulse1 can be used.'  )
      if(trim(ae_shape2) /='none') call stop_em('When art_num_em > 0, only pulse1 can be used.'  )
      if(yn_periodic     =='n'   ) call stop_em('When art_num_em > 0, yn_periodic must be y.'    )
      if(art_smedia_id_em/=media_id_source1) then
        call stop_em('When art_num_em > 0, art_smedia_id_em and media_id_source1 must be same.')
      end if
      if(fe%sig(art_smedia_id_em)>0.0d0 ) then
        call stop_em('When art_num_em > 0, sigma_em(art_smedia_id_em) must be 0.0d0.')
      end if
      if(ae_shape1=='impulse') then
        call stop_em('When art_num_em > 0, impulse can not be used.')
      end if
      if( (gbeam_sigma_plane1(1)>0.0d0).or.(gbeam_sigma_plane1(2)>0.0d0).or.(gbeam_sigma_plane1(3)>0.0d0).or. &
          (gbeam_sigma_plane2(1)>0.0d0).or.(gbeam_sigma_plane2(2)>0.0d0).or.(gbeam_sigma_plane2(3)>0.0d0).or. &
          (gbeam_sigma_line1(1) >0.0d0).or.(gbeam_sigma_line1(2) >0.0d0).or.(gbeam_sigma_line1(3) >0.0d0).or. &
          (gbeam_sigma_line2(1) >0.0d0).or.(gbeam_sigma_line2(2) >0.0d0).or.(gbeam_sigma_line2(3) >0.0d0) ) then
        call stop_em(   'When art_num_em > 0,',      &
                     c2='gauss beam cannot be used.',&
                     c3='gbeam_sigma_plane1/2 and gbeam_sigma_line1/2 must be <= 0.0d0.')
      end if
      if    ((art_ene_min_em <0.0d0).and.(art_ene_max_em <0.0d0).and.&
             (art_wav_min_em <0.0d0).and.(art_wav_max_em <0.0d0))         then
        flag_stop = .true.
      elseif((art_ene_min_em>=0.0d0).and.(art_ene_max_em>=0.0d0).and.&
             (art_wav_min_em>=0.0d0).and.(art_wav_max_em>=0.0d0))         then
        flag_stop = .true.
      elseif((art_ene_min_em>=0.0d0).and.(art_wav_min_em>=0.0d0))         then
        flag_stop = .true.
      elseif((art_ene_min_em>=0.0d0).and.(art_wav_max_em>=0.0d0)) then
        flag_stop = .true.
      elseif((art_wav_min_em>=0.0d0).and.(art_ene_max_em>=0.0d0)) then
        flag_stop = .true.
      elseif((art_ene_min_em>=0.0d0).and.(art_ene_max_em<art_ene_min_em)) then
        flag_stop = .true.
      elseif((art_wav_min_em>=0.0d0).and.(art_wav_max_em<art_wav_min_em)) then
        flag_stop = .true.
      end if
      if(flag_stop) then
        call stop_em(   'When art_num_em > 0, you must set either',     &
                     c2='[art_ene_min_em and art_ene_min_em >= 0.0d0]', &
                     c3='or',                                           &
                     c4='[art_wav_min_em and art_wav_min_em >= 0.0d0].',&
                     c5='art_ene/wav_max_em must be larger than art_ene/wav_min_em.')
      end if
      if    (ek_dir1(1)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,1),1) >= (fs%rlsize(1)/4.0d0) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='source_loc1(1) must be set < al_em(1)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,1),1) >= art_plane_bot_em(1) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='source_loc1(1) must be set < art_plane_bot_em(1).')
        end if
        if( art_plane_top_em(1) <= art_plane_bot_em(1) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='art_plane_top_em(1) must be set > art_plane_bot_em(1).')
        end if
      elseif(ek_dir1(2)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,2),2) >= (fs%rlsize(2)/4.0d0) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to y-axis,', &
                       c2='source_loc1(2) must be set < al_em(2)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,2),2) >= art_plane_bot_em(2) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to y-axis,', &
                       c2='source_loc1(2) must be set < art_plane_bot_em(2).')
        end if
        if( art_plane_top_em(2) <= art_plane_bot_em(2) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='art_plane_top_em(2) must be set > art_plane_bot_em(2).')
        end if
      elseif(ek_dir1(3)==1.0d0) then
        if( fe%coo(fe%inc_po_id(1,3),3) >= (fs%rlsize(3)/4.0d0) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to z-axis,', &
                       c2='source_loc1(3) must be set < al_em(3)/4.0d0.')
        end if
        if( fe%coo(fe%inc_po_id(1,3),3) >= art_plane_bot_em(3) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to z-axis,', &
                       c2='source_loc1(3) must be set < art_plane_bot_em(3).')
        end if
        if( art_plane_top_em(3) <= art_plane_bot_em(3) ) then
          call stop_em(   'When art_num_em > 0 with incident pulse propagating to x-axis,', &
                       c2='art_plane_top_em(3) must be set > art_plane_bot_em(3).')
        end if
      end if
      
      !set light speed for surrounding media
      fe%cspeed_smedia_art = cspeed_au / sqrt(fe%rep(art_smedia_id_em)*fe%rmu(art_smedia_id_em))
      
      !make sampling energy axis
      if    (art_ene_min_em>=0.0d0) then
        tmp_min = art_ene_min_em; tmp_max = art_ene_max_em;
      elseif(art_wav_min_em>=0.0d0) then
        tmp_min = art_wav_min_em; tmp_max = art_wav_max_em;
      end if
      allocate(fe%art_ene(art_num_em)); fe%art_ene(:)=0.0d0;
      do ii=1,art_num_em
        if (art_num_em==1) then
          fe%art_ene(ii) = tmp_min
        else
          fe%art_ene(ii) = tmp_min + ( (tmp_max-tmp_min)/(art_num_em-1) )*(dble(ii)-1)
        end if
        if(art_wav_min_em>=0.0d0) then
          !convert from wavelength to energy
          fe%art_ene(ii) = 2.0d0*pi*cspeed_au/fe%art_ene(ii)
        end if
      end do
      
      !allocate incident ex-z and hx-z
      allocate( fe%art_ex_inc_bt_ene(art_num_em,2,2),&
                fe%art_ey_inc_bt_ene(art_num_em,2,2),&
                fe%art_ez_inc_bt_ene(art_num_em,2,2),&
                fe%art_hx_inc_bt_ene(art_num_em,2,2),&
                fe%art_hy_inc_bt_ene(art_num_em,2,2),&
                fe%art_hz_inc_bt_ene(art_num_em,2,2) )
      fe%art_ex_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%art_ey_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%art_ez_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%art_hx_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%art_hy_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      fe%art_hz_inc_bt_ene(:,:,:) = (0.0d0,0.0d0);
      
      !set id & ipe on plane and allocate total e and h fields
      call find_point_id_only_em(fs%lg,itmp1,fe%Nd,fe%coo(:,:),art_plane_bot_em(:))
      call find_point_id_only_em(fs%lg,itmp2,fe%Nd,fe%coo(:,:),art_plane_top_em(:))
      if    (ek_dir1(1)==1.0d0) then !x-direction propagation, yz-plane----------------------!
        fe%id_on_plane(1) = itmp1(1) !bottom
        fe%id_on_plane(2) = itmp2(1) !top
        do ibt = 1,2
          if((fs%mg%is(1)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(1))) then
            fe%ipe_on_plane(ibt) = 1
          else
            fe%ipe_on_plane(ibt) = 0
          end if
        end do
        allocate( fe%art_ey_tot_bt_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_ez_tot_bt_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_hy_tot_bt_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_hz_tot_bt_ene(fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2) )
        fe%art_ey_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_ez_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
        fe%art_hy_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_hz_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      elseif(ek_dir1(2)==1.0d0) then !y-direction propagation, xz-plane----------------------!
        fe%id_on_plane(1) = itmp1(2) !bottom
        fe%id_on_plane(2) = itmp2(2) !top
        do ibt = 1,2
          if((fs%mg%is(2)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(2))) then
            fe%ipe_on_plane(ibt) = 1
          else
            fe%ipe_on_plane(ibt) = 0
          end if
        end do
        allocate( fe%art_ex_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_ez_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_hx_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2),&
                  fe%art_hz_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(3):fs%mg%ie(3),art_num_em,2,2) )
        fe%art_ex_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_ez_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
        fe%art_hx_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_hz_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      elseif(ek_dir1(3)==1.0d0) then !z-direction propagation, xy-plane----------------------!
        fe%id_on_plane(1) = itmp1(3) !bottom
        fe%id_on_plane(2) = itmp2(3) !top
        do ibt = 1,2
          if((fs%mg%is(3)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(3))) then
            fe%ipe_on_plane(ibt) = 1
          else
            fe%ipe_on_plane(ibt) = 0
          end if
        end do
        allocate( fe%art_ex_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),art_num_em,2,2),&
                  fe%art_ey_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),art_num_em,2,2),&
                  fe%art_hx_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),art_num_em,2,2),&
                  fe%art_hy_tot_bt_ene(fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),art_num_em,2,2) )
        fe%art_ex_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_ey_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
        fe%art_hx_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0); fe%art_hy_tot_bt_ene(:,:,:,:,:) = (0.0d0,0.0d0);
      end if
      
      !write information
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,*) "**************************"
        write(*,*) "art_num_em > 0:"
        write(*,*) "Absorption-, Reflection-, and Transmission-ratas"
        write(*,*) "will be outputed at the end of calculation."
        write(*,*) "NOTE: Those are normalized by the spectral distribution of the incident pulse."
        write(*,*)
        write(*,*) "Bottom and top planes are set at"
        if    (ek_dir1(1)==1.0d0) then !x-propagatin, yz-plane
          write(*,'(A,ES12.5,A,ES12.5)') " x = ",fe%coo(fe%id_on_plane(1),1)*ulength_from_au,&
                                         " and ",fe%coo(fe%id_on_plane(2),1)*ulength_from_au
        elseif(ek_dir1(2)==1.0d0) then !y-propagatin, xz-plane
          write(*,'(A,ES12.5,A,ES12.5)') " y = ",fe%coo(fe%id_on_plane(1),2)*ulength_from_au,&
                                         " and ",fe%coo(fe%id_on_plane(2),2)*ulength_from_au
        elseif(ek_dir1(3)==1.0d0) then !z-propagatin, xy-plane
          write(*,'(A,ES12.5,A,ES12.5)') " z = ",fe%coo(fe%id_on_plane(1),3)*ulength_from_au,&
                                         " and ",fe%coo(fe%id_on_plane(2),3)*ulength_from_au
        end if
        write(*,*) "in the unit system, ",trim(unit_system),"."
        if(obs_samp_em==1) then
          write(*,*)
          write(*, '(A)') " --- CAUTION ---"
          write(*, '(A)') " In this calculation, obs_samp_em is set to 1."
          write(*,'(2A)') " It is strongly recommended to increase obs_samp_em",&
                          " as far as your intending sampling period allows"
          write(*, '(A)') " because the calculation cost of this art-option can be reduced with larger obs_samp_em."
        end if
        write(*,*) "**************************"
      end if
    end if
    
    !*** deallocate unused variables **************************************************************************!
    deallocate(fs%imedia,fe%rmedia);
    
    !*** completely input restart date ************************************************************************!
    if(yn_restart=='y') then
      !start to calculate elapsed_time
      elapsed_time = get_wtime()
      
      !input iter_now(substituted to iter_sta with +1), e_max and h_max
      call input_r_txt_em(fe%ifn,1,1,'iter_now',i0d=fe%iter_sta, &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'
      fe%iter_sta = fe%iter_sta + 1;
      call input_r_txt_em(fe%ifn,1,1,'e_max',   r0d=fe%e_max   , &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_txt_em(fe%ifn,1,1,'h_max',   r0d=fe%h_max   , &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      
      !input binary data
      is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
      is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
      is5(3) = fs%mg%is(3); ie5(3) = fs%mg%ie(3);
      is5(4) = 1          ; ie5(4) = 1          ;
      is5(5) = 1          ; ie5(5) = 1          ;
      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ex_y',&
                          r3d=fe%ex_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ex_z',&
                          r3d=fe%ex_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ey_z',&
                          r3d=fe%ey_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ey_x',&
                          r3d=fe%ey_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ez_x',&
                          r3d=fe%ez_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ez_y',&
                          r3d=fe%ez_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hx_y',&
                          r3d=fe%hx_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hx_z',&
                          r3d=fe%hx_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hy_z',&
                          r3d=fe%hy_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hy_x',&
                          r3d=fe%hy_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hz_x',&
                          r3d=fe%hz_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','hz_y',&
                          r3d=fe%hz_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
      if(restart_io/=0)error stop 'Required EM restart read failed'

      call eh_sendrecv(fs,fe,'e'); call eh_sendrecv(fs,fe,'h');
      if(ld%enabled) then
        is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1)       ;
        is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2)       ;
        is5(3) = fs%mg%is(3); ie5(3) = fs%mg%ie(3)       ;
        is5(4) = 1          ; ie5(4) = ld%max_pole_num_ld;
        is5(5) = 1          ; ie5(5) = 1                 ;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjx_ld',&
                            r4d=ld%rjx_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjy_ld',&
                            r4d=ld%rjy_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjz_ld',&
                            r4d=ld%rjz_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','px_ld' ,&
                            r4d= ld%px_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','py_ld' ,&
                            r4d= ld%py_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','pz_ld' ,&
                            r4d= ld%pz_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        is5(4) = 1; ie5(4) = 1;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjx_fdtd_ld',&
                            r3d=ld%rjx_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjy_fdtd_ld',&
                            r3d=ld%rjy_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','rjz_fdtd_ld',&
                            r3d=ld%rjz_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

      end if
      if(fe%flag_lr) then
        is5(1) = 1; ie5(1) = nt_em;
        is5(2) = 1; ie5(2) = 1;
        is5(3) = 1; ie5(3) = 1;
        is5(4) = 1; ie5(4) = 1;
        is5(5) = 1; ie5(5) = 1;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','time_lr',r1d=fe%time_lr, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        is5(2) = 1; ie5(2) = 3;
        if(yn_periodic=='n') then
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','dip_lr' ,r2d=fe%dip_lr , &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        elseif(yn_periodic=='y') then
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','curr_lr',r2d=fe%curr_lr, &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','e_lr'   ,r2d=fe%e_lr   , &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end if
      end if
      if(fe%flag_obs) then
      if(sum(fe%iobs_num_ene(:))>0) then
        is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1)               ;
        is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2)               ;
        is5(3) = 1          ; ie5(3) = 1                         ;
        is5(4) = 1          ; ie5(4) = maxval(fe%iobs_num_ene(:));
        is5(5) = 1          ; ie5(5) = 2                         ;
        do ii=1,obs_num_em
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ex_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ex_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ey_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ey_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ez_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ez_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hx_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hx_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hy_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hy_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hz_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hz_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jx_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jx_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jy_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jy_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jz_xy_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jz_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
        is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
        is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
        do ii=1,obs_num_em
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ex_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ex_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ey_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ey_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ez_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ez_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hx_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hx_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hy_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hy_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hz_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hz_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jx_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jx_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jy_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jy_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jz_yz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jz_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
        is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
        is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
        do ii=1,obs_num_em
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ex_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ex_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ey_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ey_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_ez_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_ez_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hx_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hx_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hy_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hy_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_hz_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_hz_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jx_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jx_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jy_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jy_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','obs_jz_xz_ene_iobs'//trim(adjustl(save_name)),&
                              c5d=fe%obs_jz_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
      end if
      end if
      if(fe%flag_ase) then
        if    (ek_dir1(1)==1.0d0) then
          is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
        elseif(ek_dir1(2)==1.0d0) then
          is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
        elseif(ek_dir1(3)==1.0d0) then
          is5(1) = fs%mg%is(3); ie5(1) = fs%mg%ie(3);
        end if
        is5(2) = 1; ie5(2) = ase_num_em;
        is5(3) = 1; ie5(3) = 2;
        is5(4) = 1; ie5(4) = 1;
        is5(5) = 1; ie5(5) = 1;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ex_inc_pa_ene',c3d=fe%ase_ex_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ey_inc_pa_ene',c3d=fe%ase_ey_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ez_inc_pa_ene',c3d=fe%ase_ez_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hx_inc_pa_ene',c3d=fe%ase_hx_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hy_inc_pa_ene',c3d=fe%ase_hy_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hz_inc_pa_ene',c3d=fe%ase_hz_inc_pa_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        is5(1) = 1; ie5(1) = ase_num_em;
        is5(2) = 1; ie5(2) = 2;
        is5(3) = 1; ie5(3) = 2;
        is5(4) = 1; ie5(4) = 1;
        is5(5) = 1; ie5(5) = 1;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_ex_inc_bt_ene',c3d=fe%ase_ex_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_ey_inc_bt_ene',c3d=fe%ase_ey_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_ez_inc_bt_ene',c3d=fe%ase_ez_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_hx_inc_bt_ene',c3d=fe%ase_hx_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_hy_inc_bt_ene',c3d=fe%ase_hy_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','ase_hz_inc_bt_ene',c3d=fe%ase_hz_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
        is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
        is5(3) = 1;           ie5(3) = ase_num_em;
        is5(4) = 1;           ie5(4) = 1;
        is5(5) = 1;           ie5(5) = 2;
        do ii = 1,2
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ex_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ex_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ey_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ey_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hx_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hx_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hy_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hy_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
        is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
        is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
        do ii = 1,2
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ey_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ey_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ez_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ez_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hy_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hy_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hz_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hz_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
        is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
        is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
        do ii = 1,2
          write(save_name,*) ii
          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ex_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ex_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_ez_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_ez_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hx_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hx_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

          call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','ase_hz_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                              c5d=fe%ase_hz_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
          if(restart_io/=0)error stop 'Required EM restart read failed'

        end do
      end if
      if(fe%flag_art) then
        is5(1) = 1; ie5(1) = art_num_em;
        is5(2) = 1; ie5(2) = 2         ;
        is5(3) = 1; ie5(3) = 2         ;
        is5(4) = 1; ie5(4) = 1         ;
        is5(5) = 1; ie5(5) = 1         ;
        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_ex_inc_bt_ene',c3d=fe%art_ex_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_ey_inc_bt_ene',c3d=fe%art_ey_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_ez_inc_bt_ene',c3d=fe%art_ez_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_hx_inc_bt_ene',c3d=fe%art_hx_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_hy_inc_bt_ene',c3d=fe%art_hy_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'single','art_hz_inc_bt_ene',c3d=fe%art_hz_inc_bt_ene, &
            io_status=restart_io)
        if(restart_io/=0)error stop 'Required EM restart read failed'

        if    (ek_dir1(1)==1.0d0) then
          is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
          is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
          is5(3) = 1          ; ie5(3) = art_num_em ;
          is5(4) = 1          ; ie5(4) = 1          ;
          is5(5) = 1          ; ie5(5) = 2          ;
          do ii = 1,2
            write(save_name,*) ii
            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ey_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ey_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ez_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ez_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hy_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hy_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hz_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hz_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

          end do
        elseif(ek_dir1(2)==1.0d0) then
          is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
          is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
          is5(3) = 1          ; ie5(3) = art_num_em ;
          is5(4) = 1          ; ie5(4) = 1          ;
          is5(5) = 1          ; ie5(5) = 2          ;
          do ii = 1,2
            write(save_name,*) ii
            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ex_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ex_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ez_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ez_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hx_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hx_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hz_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hz_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

          end do
        elseif(ek_dir1(3)==1.0d0) then
          is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
          is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
          is5(3) = 1          ; ie5(3) = art_num_em ;
          is5(4) = 1          ; ie5(4) = 1          ;
          is5(5) = 1          ; ie5(5) = 2          ;
          do ii = 1,2
            write(save_name,*) ii
            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ex_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ex_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_ey_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_ey_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hx_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hx_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

            call input_r_bin_em(fe%ifn,is5,ie5,nsg_p,flag_same_p,'all','art_hy_tot_bt_ene_ibt'//&
                                trim(adjustl(save_name)),c5d=fe%art_hy_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
            if(restart_io/=0)error stop 'Required EM restart read failed'

          end do
        end if            
      end if
      
      !copy obsx_at_point_rt.data and obsx_plane_integral_rt.data
      if(fe%flag_obs .and. comm_is_root(nproc_id_global)) then
        do ii=1,obs_num_em
          write(tmp_c(2),*) ii
          tmp_c(1) = trim(adjustl(base_directory))     //'/obs'//trim(adjustl(tmp_c(2)))//'_at_point_rt.data'
          tmp_c(2) = trim(adjustl(directory_read_data))//'/obs'//trim(adjustl(tmp_c(2)))//'_at_point_rt.data'
          call txtfile_copy_em(trim(adjustl(tmp_c(2))),trim(adjustl(tmp_c(1))))
          if(yn_obs_plane_integral_em(ii)=='y') then
            write(tmp_c(2),*) ii
            tmp_c(1) = trim(adjustl(base_directory))     //'/obs'//trim(adjustl(tmp_c(2)))//'_xy_integral_rt.data'
            tmp_c(2) = trim(adjustl(directory_read_data))//'/obs'//trim(adjustl(tmp_c(2)))//'_xy_integral_rt.data'
            call txtfile_copy_em(trim(adjustl(tmp_c(2))),trim(adjustl(tmp_c(1))))
            write(tmp_c(2),*) ii
            tmp_c(1) = trim(adjustl(base_directory))     //'/obs'//trim(adjustl(tmp_c(2)))//'_yz_integral_rt.data'
            tmp_c(2) = trim(adjustl(directory_read_data))//'/obs'//trim(adjustl(tmp_c(2)))//'_yz_integral_rt.data'
            call txtfile_copy_em(trim(adjustl(tmp_c(2))),trim(adjustl(tmp_c(1))))
            write(tmp_c(2),*) ii
            tmp_c(1) = trim(adjustl(base_directory))     //'/obs'//trim(adjustl(tmp_c(2)))//'_xz_integral_rt.data'
            tmp_c(2) = trim(adjustl(directory_read_data))//'/obs'//trim(adjustl(tmp_c(2)))//'_xz_integral_rt.data'
            call txtfile_copy_em(trim(adjustl(tmp_c(2))),trim(adjustl(tmp_c(1))))
          end if
        end do
      end if
      
      !write elapsed time
      elapsed_time = get_wtime() - elapsed_time
      if(comm_is_root(nproc_id_global)) then
        write(*,*)
        write(*,'(A,f16.8)') " elapsed time for inputing restart data [s] = ", elapsed_time
      end if
    end if
    
    !*** write start ******************************************************************************************!
    ! Publish the new layout only after every required old payload was read.
    ! Input/output restart directories may be the same during repartition.
    if(fe%flag_o_restart)then
      restart_io=0
      !output nproc_size_global, nt_em, and obs_samp_em
      call output_r_txt_em(fe%ifn,1,1,'single','nproc_size_global',i0d=nproc_size_global, &
            io_status=restart_io)
      call output_r_txt_em(fe%ifn,1,1,'single','nt_em',            i0d=nt_em            , &
            io_status=restart_io)
      call output_r_txt_em(fe%ifn,1,1,'single','obs_samp_em',      i0d=obs_samp_em      , &
            io_status=restart_io)

      !output array size
      call output_r_txt_em(fe%ifn,1,3,'all','mg_is'      ,i1d=fs%mg%is      , &
            io_status=restart_io)
      call output_r_txt_em(fe%ifn,1,3,'all','mg_ie'      ,i1d=fs%mg%ie      , &
            io_status=restart_io)
      call output_r_txt_em(fe%ifn,1,3,'all','mg_is_array',i1d=fs%mg%is_array, &
            io_status=restart_io)
      call output_r_txt_em(fe%ifn,1,3,'all','mg_ie_array',i1d=fs%mg%ie_array, &
            io_status=restart_io)
      restart_bad=0
      if(restart_io/=0)restart_bad=1
      call comm_summation(restart_bad,restart_total,nproc_group_global)
      if(restart_total/=0)error stop 'Required EM restart metadata save failed; incomplete checkpoint'
    endif

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
    subroutine eh_allocate_vector(values,lo,hi)
      real(8),allocatable,intent(out) :: values(:)
      integer,intent(in) :: lo,hi
      integer :: ios
      if(present(status))then
        if(allocation_status/=0)return
        ios=0
#ifdef SALMONPY_RECOVERY_TESTING
        auxiliary_slot=auxiliary_slot+1
        if(present(allocation_failure_slot))then
          if(100+auxiliary_slot==allocation_failure_slot)ios=1
        endif
#endif
        if(ios==0)allocate(values(lo:hi),stat=ios)
        if(ios/=0)allocation_status=1120
      else
        allocate(values(lo:hi))
      endif
    end subroutine eh_allocate_vector

    subroutine eh_allocate(ista,iend,allocate_mode,num4d,num5d,i3d,r3d,i4d,r4d,i5d,r5d,allocation_status)
      implicit none
      integer,            intent(in)           :: ista(3),iend(3)
      character(3),       intent(in)           :: allocate_mode
      integer,            intent(in), optional :: num4d,num5d
      integer,allocatable,intent(out),optional :: i3d(:,:,:)
      real(8),allocatable,intent(out),optional :: r3d(:,:,:)
      integer,allocatable,intent(out),optional :: i4d(:,:,:,:)
      real(8),allocatable,intent(out),optional :: r4d(:,:,:,:)
      integer,allocatable,intent(out),optional :: i5d(:,:,:,:,:)
      real(8),allocatable,intent(out),optional :: r5d(:,:,:,:,:)
    
      integer,optional,intent(inout) :: allocation_status
      integer :: ios

      ! Native core arrays return STAT failures. Standalone and later optional
      ! arrays keep their existing path; never clear an earlier phase failure.
      if(present(allocation_status).and.present(status))then
        if(allocation_status/=0)return
        ios=0
#ifdef SALMONPY_RECOVERY_TESTING
        allocation_slot=allocation_slot+1
        if(present(allocation_failure_slot))then
          if(allocation_slot==allocation_failure_slot)ios=1
        endif
#endif
        if(ios==0)then
          select case(allocate_mode)
          case('i3d')
            allocate(i3d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3)),stat=ios)
          case('r3d')
            allocate(r3d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3)),stat=ios)
          case default
            ios=1
          end select
        endif
        if(ios/=0)then
          allocation_status=1;return
        endif
        select case(allocate_mode)
        case('i3d')
          i3d(:,:,:)=0
        case('r3d')
          r3d(:,:,:)=0d0
        end select
        return
      endif

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
      case('r4d')
        allocate(r4d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3),num4d))
        r4d(:,:,:,:)=0
      case('i5d')
        allocate(i5d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3),num4d,num5d))
        i5d(:,:,:,:,:)=0.0d0
      case('r5d')
        allocate(r5d(ista(1):iend(1),ista(2):iend(2),ista(3):iend(3),num4d,num5d))
        r5d(:,:,:,:,:)=0.0d0
      end select
      
      return
    end subroutine eh_allocate
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ set fdtd coefficient ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_material_allocation_start(slot,ios)
      integer,intent(in) :: slot
      integer,intent(out) :: ios
      ios=0
#ifdef SALMONPY_RECOVERY_TESTING
      if(present(allocation_failure_slot))then
        if(allocation_failure_slot==200+slot)ios=1
      endif
#endif
    end subroutine eh_material_allocation_start

    subroutine eh_material_allocation_result(name,ios,ierr)
      character(*),intent(in) :: name
      integer,intent(in) :: ios
      integer,intent(inout) :: ierr
      if(ios==0)return
      ierr=2
      write(*,'(A,I0,A,I0,2A)')'SALMONPY_MATERIAL_ALLOCATION_ERROR rank=',nproc_id_global, &
        ' stat=',ios,' array=',name
    end subroutine eh_material_allocation_result

    subroutine eh_read_diagonal(ierr)
      integer,intent(out)::ierr
      integer::unit,ios,nrecord,version,record,imedia,np,q,axis
      character(64)::magic
      logical::seen(0:200)
      ierr=0;seen=.false.
      if(ierr==0)then
        call eh_material_allocation_start(1,ios)
        if(ios==0)allocate(diagonal_eps(3,0:media_num),stat=ios)
        call eh_material_allocation_result('diagonal_eps',ios,ierr)
      endif
      if(ierr==0)then
        call eh_material_allocation_start(2,ios)
        if(ios==0)allocate(diagonal_f(3,0:media_num,100),stat=ios)
        call eh_material_allocation_result('diagonal_f',ios,ierr)
      endif
      if(ierr==0)then
        do axis=1,3
          diagonal_eps(axis,:)=epsilon_em(0:media_num)
          diagonal_f(axis,:,:)=f_ld(0:media_num,:)
        enddo
      endif
      if(comm_is_root(nproc_id_global).and.ierr==0)then
        open(newunit=unit,file=resolve_path('salmonpy-ld-diagonal.dat'),status='old',action='read',iostat=ios)
        if(ios/=0)then
          ierr=1
        else
          read(unit,*,iostat=ios)magic,version
          if(ios==0)then
            if(magic/='SALMONPY_LD_DIAGONAL'.or.version/=1)ios=1
          endif
          if(ios==0)read(unit,*,iostat=ios)nrecord
          if(ios==0)then
            if(nrecord/=count(media_type(0:media_num)=='ld-diagonal'))ios=1
          endif
          if(ios==0)then
            do record=1,nrecord
              read(unit,*,iostat=ios)imedia,np
              if(ios/=0)exit
              if(imedia<1.or.imedia>media_num.or.np<1.or.np>100)then
                ios=1;exit
              endif
              if(seen(imedia).or.media_type(imedia)/='ld-diagonal'.or.np/=pole_num_ld(imedia))then
                ios=1;exit
              endif
              seen(imedia)=.true.
              read(unit,*,iostat=ios)diagonal_eps(:,imedia)
              if(ios/=0)exit
              do q=1,np
                read(unit,*,iostat=ios)diagonal_f(:,imedia,q)
                if(ios/=0)exit
              enddo
              if(ios/=0)exit
            enddo
          endif
          if(ios/=0)ierr=1
          close(unit)
        endif
      endif
      ! Same startup gate as the former root status broadcast; no extra collective.
      call comm_get_max(ierr,nproc_group_global)
      if(ierr/=0)return
      call comm_bcast(diagonal_eps,nproc_group_global)
      call comm_bcast(diagonal_f,nproc_group_global)
      do imedia=0,media_num
        if(media_type(imedia)/='ld-diagonal')cycle
        do axis=1,3
          if(.not.ieee_is_finite(diagonal_eps(axis,imedia)))ierr=1
          if(diagonal_eps(axis,imedia)<=0d0)ierr=1
          do q=1,pole_num_ld(imedia)
            if(.not.ieee_is_finite(diagonal_f(axis,imedia,q)))ierr=1
            if(diagonal_f(axis,imedia,q)<0d0)ierr=1
          enddo
        enddo
      enddo
    end subroutine

    subroutine eh_read_subpixel(ierr)
      integer,intent(out)::ierr
      integer::dims(3)
      ierr=0
      if(comm_is_root(nproc_id_global)) &
        call read_subpixel_profile(resolve_path('salmonpy-subpixel.dat'),subpixel_profile,ierr)
      call comm_bcast(ierr,nproc_group_global)
      if(ierr/=0)return
      if(comm_is_root(nproc_id_global))dims=[subpixel_profile%nm,subpixel_profile%ns,subpixel_profile%nq]
      call comm_bcast(dims,nproc_group_global)
      if(.not.comm_is_root(nproc_id_global))then
        subpixel_profile%nm=dims(1);subpixel_profile%ns=dims(2);subpixel_profile%nq=dims(3)
        call allocate_subpixel_profile(subpixel_profile)
      endif
      if(dims(1)/=media_num)ierr=1
      call comm_bcast(subpixel_profile%omega,nproc_group_global)
      call comm_bcast(subpixel_profile%rtol,nproc_group_global)
      call comm_bcast(subpixel_profile%maxiter,nproc_group_global)
      call comm_bcast(subpixel_profile%background,nproc_group_global)
      call comm_bcast(subpixel_profile%strength,nproc_group_global)
      call comm_bcast(subpixel_profile%resonance,nproc_group_global)
      call comm_bcast(subpixel_profile%damping,nproc_group_global)
      call comm_bcast(subpixel_profile%owner,nproc_group_global)
      call comm_bcast(subpixel_profile%ld_phase,nproc_group_global)
      call comm_bcast(subpixel_profile%shape_size,nproc_group_global)
      call comm_bcast(subpixel_profile%center,nproc_group_global)
      call comm_bcast(subpixel_profile%kind,nproc_group_global)
      call comm_bcast(subpixel_profile%material,nproc_group_global)
      if(.not.any(subpixel_profile%ld_phase/=0))then
        if(comm_is_root(nproc_id_global))write(*,'(A)') &
          'WARNING: no LD phase specified; reference subpixel smoothing not applied.'
        subpixel_profile%nm=0
      endif
    end subroutine

    subroutine eh_subpixel_coefficients(ierr)
      integer,intent(out)::ierr
      real(8),allocatable::effective_f(:,:,:),c0(:),c1(:),c3(:)
      real(8)::background(3,3),point(3),matrix(3,3),m4(3),m5(3),drive,err,maxerr,bytes,global_error(1)
      integer::lo(3),hi(3),nq,ix,iy,iz,axis,q,a,b,flag,have_cross,m,next(3),active_lo(3),active_hi(3),cell_index(3)
      real(8)::counts(3),total_counts(3),local_bounds(2),global_bounds(2),dt_bound
      integer::unit,ios,close_status
      logical::live
      logical,allocatable::excluded(:,:,:)
      integer::p,ox,oy,oz,site(3),index(3),excluded_count
      ierr=0;counts=0;maxerr=0d0;local_bounds=1d0
      lo=fs%mg%is_array;hi=fs%mg%ie_array;nq=subpixel_profile%nq
      allocate(excluded(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      excluded=.false.;excluded_count=0
      if(present(excluded_response_points))then
        excluded_count=size(excluded_response_points,2)
        do p=1,excluded_count
          site=fs%lg%is+excluded_response_points(:,p)
          do oz=-1,1
          do oy=-1,1
          do ox=-1,1
            index=site+[ox,oy,oz]
            if(all(index>=lo).and.all(index<=hi))excluded(index(1),index(2),index(3))=.true.
          enddo
          enddo
          enddo
        enddo
      endif
      ! A box spanning a periodic direction represents a continuous extrusion.
      do q=1,subpixel_profile%ns
        if(subpixel_profile%kind(q)/=1)cycle
        do axis=1,3
          if(fs%a_bc(axis,1)/='periodic'.or.fs%a_bc(axis,2)/='periodic')cycle
          if(subpixel_profile%shape_size(axis,q)>=fs%rlsize(axis)) &
            subpixel_profile%shape_size(axis,q)=8d0*fs%rlsize(axis)+2d0*abs(subpixel_profile%center(axis,q))
        enddo
      enddo
      if(allocated(fe%tensor_raw))deallocate(fe%tensor_raw)
      if(allocated(fe%tensor_corr))deallocate(fe%tensor_corr)
      if(allocated(fe%tensor_rates))deallocate(fe%tensor_rates)
      allocate(fe%tensor_raw(3,3,0:subpixel_profile%nm+nq,0:0),fe%tensor_corr(3,3,0:nq,0:1), &
        fe%tensor_rates(3,nq,0:0),fe%tensor_weight(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0:nq), &
        fe%tensor_diag(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        fe%tensor_owner(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        fe%subpixel_corr(3,3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0:nq),fe%subpixel_cross(0:nq))
      fe%tensor_raw=0d0;fe%tensor_corr=0d0;fe%tensor_weight=0d0
      fe%tensor_diag=1d0;fe%tensor_owner=0;fe%subpixel_corr=0d0
      fe%subpixel_enabled=.true.;fe%tensor_enabled=.true.;fe%subpixel_reference=subpixel_profile%omega
      fe%tensor_rtol=subpixel_profile%rtol;fe%tensor_maxiter=subpixel_profile%maxiter
      do m=0,subpixel_profile%nm
        fe%tensor_raw(:,:,m,0)=subpixel_profile%background(:,:,m)
      enddo
      do q=1,nq
        fe%tensor_raw(:,:,subpixel_profile%nm+q,0)=subpixel_profile%strength(:,:,q)
        fe%tensor_rates(:,q,0)=[1d0,subpixel_profile%resonance(q),subpixel_profile%damping(q)]
      enddo
      ! Every oscillator slot has one rate pair across the whole domain.
      ! This keeps the existing ADE recurrence compatible with its tensor drive.
      call ld_finalize(ld)
      ld%enabled=.true.;ld%num_ld=1;ld%max_pole_num_ld=nq
      allocate(ld%media_ld(1));ld%media_ld(1)=min(1,subpixel_profile%nm)
      bytes=real(product(int(hi-lo+1,8)),8)*8d0*real(30*nq+12,8)
      if(comm_is_root(nproc_id_global))then
        write(*,'(A)')'SALMONPY_SUBPIXEL_REFERENCE_V2'
        write(*,'(A,ES16.8,A,I0,A,ES16.8)')'Subpixel reference omega=',subpixel_profile%omega, &
          '; global oscillator slots=',nq,'; approximate local coefficient/state bytes=',bytes
      endif
      call ld_allocate(ld,lo,hi)
      allocate(effective_f(3,3,nq),c0(nq),c1(nq),c3(nq))
      do q=1,nq
        c0(q)=1d0+subpixel_profile%damping(q)*dt_em/2d0+(subpixel_profile%resonance(q)*dt_em)**2/4d0
        c1(q)=(1d0-subpixel_profile%damping(q)*dt_em/2d0 &
          -(subpixel_profile%resonance(q)*dt_em)**2/4d0)/c0(q)
        c3(q)=subpixel_profile%resonance(q)**2*dt_em/c0(q)
        if(.not.ieee_is_finite(c0(q)).or.c0(q)<=0d0)ierr=1
        if(.not.ieee_is_finite(c1(q)).or..not.ieee_is_finite(c3(q)))ierr=1
      enddo
      call comm_get_max(ierr,nproc_group_global)
      if(ierr/=0)return
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        cell_index=[ix,iy,iz]
        point=[fe%coo(ix,1),fe%coo(iy,2),fe%coo(iz,3)]
        if(excluded(ix,iy,iz))then
          call subpixel_bulk(subpixel_profile,fs%imedia(ix,iy,iz),background,effective_f)
          flag=5;err=0d0
        else
          call subpixel_at(subpixel_profile,point,fs%hgs,background,effective_f,flag,err)
        endif
        call count_sample(flag,counts)
        maxerr=max(maxerr,err)
        do q=0,nq
          matrix=background
          if(q>0)matrix=effective_f(:,:,q)
          do b=1,3
            do a=1,3
              if(a==b.or.matrix(a,a)==0d0.or.matrix(b,b)==0d0)cycle
              fe%subpixel_corr(a,b,ix,iy,iz,q)=.5d0*(matrix(a,b)+matrix(b,a))/sqrt(matrix(a,a))/sqrt(matrix(b,b))
            enddo
          enddo
        enddo
        matrix=fe%subpixel_corr(:,:,ix,iy,iz,0)
        do a=1,3
          matrix(a,a)=1d0
        enddo
        local_bounds(1)=min(local_bounds(1),subpixel_spd_bound(matrix))
        do axis=1,3
          select case(axis)
          case(1)
            active_lo=max(fe%iex_y_is,fe%iex_z_is);active_hi=min(fe%iex_y_ie,fe%iex_z_ie)
          case(2)
            active_lo=max(fe%iey_z_is,fe%iey_x_is);active_hi=min(fe%iey_z_ie,fe%iey_x_ie)
          case(3)
            active_lo=max(fe%iez_x_is,fe%iez_y_is);active_hi=min(fe%iez_x_ie,fe%iez_y_ie)
          end select
          live=all([ix,iy,iz]>=active_lo).and.all([ix,iy,iz]<=active_hi)
          do a=1,3
            if(fs%a_bc(a,1)=='pml'.and.cell_index(a)<fs%lg%is(a)+9)live=.false.
            if(fs%a_bc(a,2)=='pml'.and.cell_index(a)>fs%lg%ie(a)-9)live=.false.
          enddo
          if(.not.live)cycle
          point=[fe%coo(ix,1),fe%coo(iy,2),fe%coo(iz,3)];point(axis)=point(axis)+.5d0*fs%hgs(axis)
          if(excluded(ix,iy,iz))then
            flag=5;err=0d0
          else
            call subpixel_at(subpixel_profile,point,fs%hgs,background,effective_f,flag,err)
          endif
          if(flag>=2)then
            m=fs%imedia(ix,iy,iz);next=[ix,iy,iz];next(axis)=next(axis)+1
            if(fs%imedia(next(1),next(2),next(3))==0)m=0
            call subpixel_bulk(subpixel_profile,m,background,effective_f)
          endif
          call count_sample(flag,counts);maxerr=max(maxerr,err)
          fe%tensor_owner(axis,ix,iy,iz)=1
          fe%tensor_weight(axis,ix,iy,iz,0)=sqrt(background(axis,axis))
          local_bounds(2)=min(local_bounds(2),background(axis,axis))
          m4(axis)=0d0
          do q=1,nq
            drive=effective_f(axis,axis,q)*dt_em/(8d0*pi*c0(q))
            if(.not.ieee_is_finite(drive).or.drive<0d0)then
              ierr=1;cycle
            endif
            fe%tensor_weight(axis,ix,iy,iz,q)=sqrt(drive);m4(axis)=m4(axis)+drive
            select case(axis)
            case(1)
              ld%c1_jx_ld(ix,iy,iz,q)=c1(q);ld%c2_jx_ld(ix,iy,iz,q)=drive;ld%c3_jx_ld(ix,iy,iz,q)=c3(q)
            case(2)
              ld%c1_jy_ld(ix,iy,iz,q)=c1(q);ld%c2_jy_ld(ix,iy,iz,q)=drive;ld%c3_jy_ld(ix,iy,iz,q)=c3(q)
            case(3)
              ld%c1_jz_ld(ix,iy,iz,q)=c1(q);ld%c2_jz_ld(ix,iy,iz,q)=drive;ld%c3_jz_ld(ix,iy,iz,q)=c3(q)
            end select
          enddo
          m5(axis)=background(axis,axis)-2d0*pi*dt_em*m4(axis)
          m4(axis)=background(axis,axis)+2d0*pi*dt_em*m4(axis)
          fe%tensor_diag(axis,ix,iy,iz)=m4(axis)
          if(.not.ieee_is_finite(m4(axis)).or.m4(axis)<=0d0)ierr=1
          if(.not.ieee_is_finite(m5(axis)))ierr=1
          select case(axis)
          case(1)
            fe%c1_ex_y(ix,iy,iz)=m5(axis)/m4(axis);fe%c1_ex_z(ix,iy,iz)=m5(axis)/m4(axis)
            fe%c2_ex_y(ix,iy,iz)=cspeed_au*dt_em/m4(axis)/fs%hgs(2)
            fe%c2_ex_z(ix,iy,iz)=-cspeed_au*dt_em/m4(axis)/fs%hgs(3)
            fe%c2_jx(ix,iy,iz)=-4d0*pi*dt_em/m4(axis)
          case(2)
            fe%c1_ey_z(ix,iy,iz)=m5(axis)/m4(axis);fe%c1_ey_x(ix,iy,iz)=m5(axis)/m4(axis)
            fe%c2_ey_z(ix,iy,iz)=cspeed_au*dt_em/m4(axis)/fs%hgs(3)
            fe%c2_ey_x(ix,iy,iz)=-cspeed_au*dt_em/m4(axis)/fs%hgs(1)
            fe%c2_jy(ix,iy,iz)=-4d0*pi*dt_em/m4(axis)
          case(3)
            fe%c1_ez_x(ix,iy,iz)=m5(axis)/m4(axis);fe%c1_ez_y(ix,iy,iz)=m5(axis)/m4(axis)
            fe%c2_ez_x(ix,iy,iz)=cspeed_au*dt_em/m4(axis)/fs%hgs(1)
            fe%c2_ez_y(ix,iy,iz)=-cspeed_au*dt_em/m4(axis)/fs%hgs(2)
            fe%c2_jz(ix,iy,iz)=-4d0*pi*dt_em/m4(axis)
          end select
        enddo
      enddo
      enddo
      enddo
      call comm_get_max(ierr,nproc_group_global)
      if(ierr/=0)return
      do q=0,nq
        have_cross=merge(1,0,any(fe%subpixel_corr(:,:,:,:,:,q)/=0d0))
        call comm_get_max(have_cross,nproc_group_global)
        fe%subpixel_cross(q)=have_cross/=0
      enddo
      call comm_summation(counts,total_counts,3,nproc_group_global)
      call comm_get_max([maxerr],global_error,1,nproc_group_global)
      maxerr=global_error(1)
      call comm_get_min(local_bounds,global_bounds,2,nproc_group_global)
      dt_bound=sqrt(max(0d0,product(global_bounds)))/(cspeed_au*sqrt(sum((1d0/fs%hgs)**2)))
      if(comm_is_root(nproc_id_global))then
        if(dt_bound<=0d0)then
          write(*,'(A)')'WARNING: subpixel background CFL bound unavailable; evolution continues.'
        else if(dt_em>=dt_bound)then
          write(*,'(A,ES16.8)')'WARNING: dt exceeds conservative subpixel background bound; continuing. bound=',dt_bound
        endif
        if(total_counts(2)>0)write(*,'(A,I0)') &
          'WARNING: subpixel corners/multiple interfaces use unsmoothed coefficients; samples=',int(total_counts(2),8)
        if(total_counts(3)>0)write(*,'(A,I0)') &
          'WARNING: singular reference uses unsmoothed coefficients; samples=',int(total_counts(3),8)
        open(newunit=unit,file=resolve_path('salmonpy-subpixel-summary.json'),status='replace',action='write',iostat=ios)
        if(ios==0)then
          if(ios==0)write(unit,'(A)',iostat=ios)'{'
          if(ios==0)write(unit,'(A,ES24.16E3,A)',iostat=ios)' "reference_omega": ',subpixel_profile%omega,','
          if(ios==0)write(unit,'(A,I0,A)',iostat=ios)' "global_poles": ',nq,','
          if(ios==0)write(unit,'(A,ES24.16E3,A)',iostat=ios)' "conservative_dt_bound": ',dt_bound,','
          if(ios==0)write(unit,'(A,I0,A)',iostat=ios)' "excluded_response_points": ',excluded_count,','
          if(ios==0)write(unit,'(A,I0,A)',iostat=ios)' "mixed_samples": ',int(total_counts(1),8),','
          if(ios==0)write(unit,'(A,I0,A)',iostat=ios)' "geometry_fallback_samples": ',int(total_counts(2),8),','
          if(ios==0)write(unit,'(A,I0,A)',iostat=ios)' "reference_fallback_samples": ',int(total_counts(3),8),','
          if(ios==0)write(unit,'(A,ES24.16E3)',iostat=ios)' "reference_match_max_abs_error": ',maxerr
          if(ios==0)write(unit,'(A)',iostat=ios)'}'
          close(unit,iostat=close_status)
          if(ios/=0.or.close_status/=0) &
            write(*,'(A)')'WARNING: subpixel summary write incomplete; evolution continues.'
        else
          write(*,'(A)')'WARNING: subpixel summary could not be written; evolution continues.'
        endif
      endif
    end subroutine

    subroutine count_sample(flag,counts)
      integer,intent(in)::flag
      real(8),intent(inout)::counts(3)
      if(flag==1)counts(1)=counts(1)+1
      if(flag==2.or.flag==3)counts(2)=counts(2)+1
      if(flag==4)counts(3)=counts(3)+1
    end subroutine

    subroutine eh_subpixel_impulse
      real(8),allocatable::kick(:,:,:,:),off(:,:,:,:)
      real(8)::factor,omega,gamma,value(3)
      integer::lo(3),hi(3),q,ix,iy,iz
      lo=fs%mg%is_array;hi=fs%mg%ie_array
      allocate(kick(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),off(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      kick(1,:,:,:)=tmp_x1+tmp_x2;kick(2,:,:,:)=tmp_y1+tmp_y2;kick(3,:,:,:)=tmp_z1+tmp_z2
      do q=1,ld%max_pole_num_ld
        call eh_tensor_off(fs,fe,q,kick,off)
        omega=fe%tensor_rates(2,q,0);gamma=fe%tensor_rates(3,q,0)
        factor=-2d0*(1d0+gamma*dt_em/2d0+(omega*dt_em)**2/4d0)/dt_em
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          value=factor*(fe%tensor_weight(:,ix,iy,iz,q)**2*kick(:,ix,iy,iz)+off(:,ix,iy,iz))
          ld%rjx_ld(ix,iy,iz,q)=value(1);ld%rjy_ld(ix,iy,iz,q)=value(2);ld%rjz_ld(ix,iy,iz,q)=value(3)
          ld%rjx_fdtd_ld(ix,iy,iz)=ld%rjx_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jx_ld(ix,iy,iz,q))*value(1)
          ld%rjy_fdtd_ld(ix,iy,iz)=ld%rjy_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jy_ld(ix,iy,iz,q))*value(2)
          ld%rjz_fdtd_ld(ix,iy,iz)=ld%rjz_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jz_ld(ix,iy,iz,q))*value(3)
        enddo
        enddo
        enddo
      enddo
    end subroutine

    subroutine eh_read_tensor(ierr)
      integer,intent(out)::ierr
      integer::unit,ios,nrecord,version,record,m,np,q,a,b,maxp
      real(8)::mat(3,3),inverse(3,3),scale,tol,det,minor,bound
      logical::seen(0:200)
      character(64)::magic
      maxp=maxval(pole_num_ld(0:media_num))
      ierr=0;seen=.false.
      if(maxp<1.or.maxp>100)then
        ierr=1;return
      endif
      if(.not.allocated(diagonal_eps))then
        if(ierr==0)then
          call eh_material_allocation_start(1,ios)
          if(ios==0)allocate(diagonal_eps(3,0:media_num),stat=ios)
          call eh_material_allocation_result('diagonal_eps',ios,ierr)
        endif
        if(ierr==0)then
          call eh_material_allocation_start(2,ios)
          if(ios==0)allocate(diagonal_f(3,0:media_num,100),stat=ios)
          call eh_material_allocation_result('diagonal_f',ios,ierr)
        endif
        if(ierr==0)then
          do a=1,3
            diagonal_eps(a,:)=epsilon_em(0:media_num)
            diagonal_f(a,:,:)=f_ld(0:media_num,:)
          enddo
        endif
      endif
      if(ierr==0)then
        call eh_material_allocation_start(3,ios)
        if(ios==0)allocate(fe%tensor_raw(3,3,0:maxp,0:media_num),stat=ios)
        call eh_material_allocation_result('fe%tensor_raw',ios,ierr)
      endif
      if(ierr==0)then
        call eh_material_allocation_start(4,ios)
        if(ios==0)allocate(fe%tensor_corr(3,3,0:maxp,0:media_num),stat=ios)
        call eh_material_allocation_result('fe%tensor_corr',ios,ierr)
      endif
      if(ierr==0)then
        call eh_material_allocation_start(5,ios)
        if(ios==0)allocate(fe%tensor_rates(3,maxp,0:media_num),stat=ios)
        call eh_material_allocation_result('fe%tensor_rates',ios,ierr)
      endif
      if(ierr==0)then
        fe%tensor_rates=0d0
        do m=0,media_num
          do q=1,min(maxp,pole_num_ld(m))
            fe%tensor_rates(:,q,m)=[omega_p_ld(m),omega_ld(m,q),gamma_ld(m,q)]
          enddo
        enddo
        fe%tensor_raw=0d0;fe%tensor_corr=0d0
      endif
      if(comm_is_root(nproc_id_global).and.ierr==0)then
        open(newunit=unit,file=resolve_path('salmonpy-ld-tensor.dat'),status='old',action='read',iostat=ios)
        if(ios==0)then
          read(unit,*,iostat=ios)magic,version
          if(ios==0)then
            if(magic/='SALMONPY_LD_TENSOR'.or.version/=1)ios=1
          endif
          if(ios==0)read(unit,*,iostat=ios)nrecord
          if(ios==0)then
            if(nrecord/=count(media_type(0:media_num)=='ld-tensor'))ios=1
          endif
          if(ios==0)read(unit,*,iostat=ios)fe%tensor_rtol,fe%tensor_maxiter
          if(ios==0)then
            do record=1,nrecord
              read(unit,*,iostat=ios)m,np
              if(ios/=0)exit
              if(m<1.or.m>media_num.or.np<1.or.np>maxp)then
                ios=1;exit
              endif
              if(seen(m).or.media_type(m)/='ld-tensor'.or.np/=pole_num_ld(m))then
                ios=1;exit
              endif
              seen(m)=.true.
              do q=0,np
                read(unit,*,iostat=ios)fe%tensor_raw(:,:,q,m)
                if(ios/=0)exit
              enddo
              if(ios/=0)exit
            enddo
          endif
          close(unit)
        endif
        if(ios/=0)ierr=1
      endif
      ! Agree local allocation and root input errors before broadcasting arrays.
      call comm_get_max(ierr,nproc_group_global)
      if(ierr/=0)return
      call comm_bcast(fe%tensor_raw,nproc_group_global)
      call comm_bcast(fe%tensor_rtol,nproc_group_global)
      call comm_bcast(fe%tensor_maxiter,nproc_group_global)
      if(.not.ieee_is_finite(fe%tensor_rtol))ierr=1
      if(fe%tensor_rtol<=0d0.or.fe%tensor_rtol>=1d0.or.fe%tensor_maxiter<1)ierr=1
      tol=64d0*epsilon(1d0)
      do m=1,media_num
        if(media_type(m)/='ld-tensor')cycle
        do q=0,pole_num_ld(m)
          do b=1,3
            do a=1,3
              if(.not.ieee_is_finite(fe%tensor_raw(a,b,q,m)))ierr=1
            enddo
          enddo
          if(ierr/=0)return
          scale=maxval(abs(fe%tensor_raw(:,:,q,m)))
          if(scale==0d0)then
            if(q==0)ierr=1
            cycle
          endif
          mat=fe%tensor_raw(:,:,q,m)/scale
          do a=1,3
            if(mat(a,a)<0d0)ierr=1
            if(q==0.and.mat(a,a)<=0d0)ierr=1
            do b=1,3
              if(abs(mat(a,b)-mat(b,a))>tol)ierr=1
              if(mat(a,a)==0d0.and.mat(a,b)/=0d0)ierr=1
            enddo
          enddo
          do a=1,2
            do b=a+1,3
              minor=mat(a,a)*mat(b,b)-mat(a,b)*mat(b,a)
              if(minor < -tol)ierr=1
            enddo
          enddo
          det=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)) &
             -mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1)) &
             +mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))
          if(det < -tol)ierr=1
          if(q==0.and.det<=0d0)ierr=1
          if(ierr/=0)return
          do a=1,3
            do b=1,3
              if(a==b.or.mat(a,a)==0d0.or.mat(b,b)==0d0)cycle
              fe%tensor_corr(a,b,q,m)=.5d0*(mat(a,b)+mat(b,a))/sqrt(mat(a,a))/sqrt(mat(b,b))
            enddo
          enddo
          if(q==0)then
            diagonal_eps(:,m)=[fe%tensor_raw(1,1,q,m),fe%tensor_raw(2,2,q,m),fe%tensor_raw(3,3,q,m)]
            inverse(1,:)=[mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2), &
              mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3),mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)]
            inverse(2,:)=[mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3), &
              mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1),mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)]
            inverse(3,:)=[mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1), &
              mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2),mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)]
            bound=scale*det/maxval(sum(abs(inverse),dim=2))
            if(.not.ieee_is_finite(bound).or.bound<=0d0)then
              ierr=1;return
            endif
            dt_cfl=min(dt_cfl,sqrt(bound)/(cspeed_au*sqrt(sum((1d0/fs%hgs)**2))))
          else
            diagonal_f(:,m,q)=[fe%tensor_raw(1,1,q,m),fe%tensor_raw(2,2,q,m),fe%tensor_raw(3,3,q,m)]
          endif
        enddo
      enddo
      fe%tensor_enabled=ierr==0
    end subroutine eh_read_tensor

    subroutine eh_tensor_impulse
      real(8),allocatable::kick(:,:,:,:),cross(:,:,:,:)
      real(8)::delta(3),factor
      integer::lo(3),hi(3),q,m,ix,iy,iz,k,a
      lo=fs%mg%is_array;hi=fs%mg%ie_array
      allocate(kick(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        cross(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      kick(1,:,:,:)=tmp_x1+tmp_x2;kick(2,:,:,:)=tmp_y1+tmp_y2;kick(3,:,:,:)=tmp_z1+tmp_z2
      do q=1,ld%max_pole_num_ld
        call eh_tensor_off(fs,fe,q,kick,cross)
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          delta=0d0
          do a=1,3
            m=fe%tensor_owner(a,ix,iy,iz)
            if(m==0)cycle
            do k=1,ld%num_ld
              if(ld%media_ld(k)==m)exit
            enddo
            if(q>pole_num_ld(m))cycle
            factor=-2d0*tmp_c0_ld(q,k)/dt_em
            delta(a)=factor*cross(a,ix,iy,iz)
          enddo
          if(ld%c1_jx_ld(ix,iy,iz,q)>0d0)then
            ld%rjx_ld(ix,iy,iz,q)=ld%rjx_ld(ix,iy,iz,q)+delta(1)
            ld%rjx_fdtd_ld(ix,iy,iz)=ld%rjx_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jx_ld(ix,iy,iz,q))*delta(1)
          endif
          if(ld%c1_jy_ld(ix,iy,iz,q)>0d0)then
            ld%rjy_ld(ix,iy,iz,q)=ld%rjy_ld(ix,iy,iz,q)+delta(2)
            ld%rjy_fdtd_ld(ix,iy,iz)=ld%rjy_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jy_ld(ix,iy,iz,q))*delta(2)
          endif
          if(ld%c1_jz_ld(ix,iy,iz,q)>0d0)then
            ld%rjz_ld(ix,iy,iz,q)=ld%rjz_ld(ix,iy,iz,q)+delta(3)
            ld%rjz_fdtd_ld(ix,iy,iz)=ld%rjz_fdtd_ld(ix,iy,iz)+.5d0*(1d0+ld%c1_jz_ld(ix,iy,iz,q))*delta(3)
          endif
        enddo
        enddo
        enddo
      enddo
    end subroutine eh_tensor_impulse

    subroutine eh_tensor_weights
      integer::lo(3),hi(3),ix,iy,iz,a,m,q,next(3),pmax,active_lo(3),active_hi(3)
      lo=fs%mg%is_array;hi=fs%mg%ie_array;pmax=ubound(fe%tensor_raw,3)
      allocate(fe%tensor_weight(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),0:pmax), &
        fe%tensor_diag(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        fe%tensor_owner(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      fe%tensor_weight=0d0;fe%tensor_diag=1d0;fe%tensor_owner=0
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        if(fe%c2_jx(ix,iy,iz)/=0d0)fe%tensor_diag(1,ix,iy,iz)=-4d0*pi*dt_em/fe%c2_jx(ix,iy,iz)
        if(fe%c2_jy(ix,iy,iz)/=0d0)fe%tensor_diag(2,ix,iy,iz)=-4d0*pi*dt_em/fe%c2_jy(ix,iy,iz)
        if(fe%c2_jz(ix,iy,iz)/=0d0)fe%tensor_diag(3,ix,iy,iz)=-4d0*pi*dt_em/fe%c2_jz(ix,iy,iz)
        m=fs%imedia(ix,iy,iz)
        if(media_type(m)/='ld-tensor')cycle
        do a=1,3
          select case(a)
          case(1)
            active_lo=max(fe%iex_y_is,fe%iex_z_is);active_hi=min(fe%iex_y_ie,fe%iex_z_ie)
          case(2)
            active_lo=max(fe%iey_z_is,fe%iey_x_is);active_hi=min(fe%iey_z_ie,fe%iey_x_ie)
          case(3)
            active_lo=max(fe%iez_x_is,fe%iez_y_is);active_hi=min(fe%iez_x_ie,fe%iez_y_ie)
          end select
          if(any([ix,iy,iz]<active_lo).or.any([ix,iy,iz]>active_hi))then
            select case(a)
            case(1)
              ld%c1_jx_ld(ix,iy,iz,:)=0d0;ld%c2_jx_ld(ix,iy,iz,:)=0d0;ld%c3_jx_ld(ix,iy,iz,:)=0d0
            case(2)
              ld%c1_jy_ld(ix,iy,iz,:)=0d0;ld%c2_jy_ld(ix,iy,iz,:)=0d0;ld%c3_jy_ld(ix,iy,iz,:)=0d0
            case(3)
              ld%c1_jz_ld(ix,iy,iz,:)=0d0;ld%c2_jz_ld(ix,iy,iz,:)=0d0;ld%c3_jz_ld(ix,iy,iz,:)=0d0
            end select
            cycle
          endif
          next=[ix,iy,iz];next(a)=next(a)+1
          if(fs%imedia(next(1),next(2),next(3))==0)cycle
          fe%tensor_owner(a,ix,iy,iz)=m
          fe%tensor_weight(a,ix,iy,iz,0)=sqrt(fe%tensor_raw(a,a,0,m))
          do q=1,pole_num_ld(m)
            select case(a)
            case(1)
              fe%tensor_weight(a,ix,iy,iz,q)=sqrt(ld%c2_jx_ld(ix,iy,iz,q))
            case(2)
              fe%tensor_weight(a,ix,iy,iz,q)=sqrt(ld%c2_jy_ld(ix,iy,iz,q))
            case(3)
              fe%tensor_weight(a,ix,iy,iz,q)=sqrt(ld%c2_jz_ld(ix,iy,iz,q))
            end select
          enddo
        enddo
      enddo
      enddo
      enddo
    end subroutine eh_tensor_weights

    real(8) function eh_axis_strength(axis,imedia,pole) result(value)
      integer,intent(in)::axis,imedia,pole
      value=f_ld(imedia,pole)
      if((media_type(imedia)=='ld-diagonal'.or.media_type(imedia)=='ld-tensor'))value=diagonal_f(axis,imedia,pole)
    end function

    subroutine eh_coeff
      implicit none
      real(8)  :: c1_e,c2_e_x,c2_e_y,c2_e_z,c2_j,c1_h,c2_h_x,c2_h_y,c2_h_z
      
      !prepare coefficient
      select case(media_type(ii))
      case('pec')
        c1_e=0.0d0; c2_e_x=0.0d0; c2_e_y=0.0d0; c2_e_z=0.0d0; c2_j=0.0d0;
      case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
        c1_e   = tmp_c5_ld(icount_ld) / tmp_c4_ld(icount_ld)
        c2_e_x = cspeed_au*dt_em      / tmp_c4_ld(icount_ld) / fs%hgs(1)
        c2_e_y = cspeed_au*dt_em      / tmp_c4_ld(icount_ld) / fs%hgs(2)
        c2_e_z = cspeed_au*dt_em      / tmp_c4_ld(icount_ld) / fs%hgs(3)
        c2_j   =  4.0d0*pi*dt_em      / tmp_c4_ld(icount_ld)
      case default
        c1_e   = (1.0d0-2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em) &
                /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)
        c2_e_x = (cspeed_au/fe%rep(ii)*dt_em) &
                /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(1)
        c2_e_y = (cspeed_au/fe%rep(ii)*dt_em) &
                /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(2)
        c2_e_z = (cspeed_au/fe%rep(ii)*dt_em) &
                /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)/fs%hgs(3)
        c2_j   = (4.0d0*pi/fe%rep(ii)*dt_em) &
                /(1.0d0+2.0d0*pi*fe%sig(ii)/fe%rep(ii)*dt_em)
      end select
      c1_h   = 1.0d0
      c2_h_x = cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(1)
      c2_h_y = cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(2)
      c2_h_z = cspeed_au/fe%rmu(ii)*dt_em/fs%hgs(3)
      call comm_bcast(c1_e  ,nproc_group_global)
      call comm_bcast(c2_e_x,nproc_group_global)
      call comm_bcast(c2_e_y,nproc_group_global)
      call comm_bcast(c2_e_z,nproc_group_global)
      call comm_bcast(c2_j  ,nproc_group_global)
      call comm_bcast(c1_h  ,nproc_group_global)
      call comm_bcast(c2_h_x,nproc_group_global)
      call comm_bcast(c2_h_y,nproc_group_global)
      call comm_bcast(c2_h_z,nproc_group_global)
      
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
            !------ ex and jx --------------------------------------------------!
            if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor'))then
              c1_e=diagonal_c5(icount_ld,1)/diagonal_c4(icount_ld,1)
              c2_e_x=cspeed_au*dt_em/diagonal_c4(icount_ld,1)/fs%hgs(1)
              c2_e_y=cspeed_au*dt_em/diagonal_c4(icount_ld,1)/fs%hgs(2)
              c2_e_z=cspeed_au*dt_em/diagonal_c4(icount_ld,1)/fs%hgs(3)
              c2_j=4d0*pi*dt_em/diagonal_c4(icount_ld,1)
            endif
            ! Geometric staircase ownership: lower-coordinate cell, never ID priority.
            ! A vacuum endpoint retains the legacy single-material boundary.
            if(fs%imedia(ix+1,iy,iz)/=0) then
              fe%c1_ex_y(ix,iy,iz) =  c1_e; fe%c2_ex_y(ix,iy,iz) =  c2_e_y;
              fe%c1_ex_z(ix,iy,iz) =  c1_e; fe%c2_ex_z(ix,iy,iz) = -c2_e_z;
              fe%c2_jx(ix,iy,iz)   = -c2_j;
              select case(media_type(ii))
              case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
                do ij=1,pole_num_ld(ii)
                  ld%c1_jx_ld(ix,iy,iz,ij) = tmp_c1_ld(ij,icount_ld)
                  ld%c2_jx_ld(ix,iy,iz,ij) = tmp_c2_ld(ij,icount_ld)
                  if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor')) &
                    ld%c2_jx_ld(ix,iy,iz,ij)=diagonal_c2(ij,icount_ld,1)
                  ld%c3_jx_ld(ix,iy,iz,ij) = tmp_c3_ld(ij,icount_ld)
                end do
              end select
            end if
            
            !------ ey and jy --------------------------------------------------!
            if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor'))then
              c1_e=diagonal_c5(icount_ld,2)/diagonal_c4(icount_ld,2)
              c2_e_x=cspeed_au*dt_em/diagonal_c4(icount_ld,2)/fs%hgs(1)
              c2_e_y=cspeed_au*dt_em/diagonal_c4(icount_ld,2)/fs%hgs(2)
              c2_e_z=cspeed_au*dt_em/diagonal_c4(icount_ld,2)/fs%hgs(3)
              c2_j=4d0*pi*dt_em/diagonal_c4(icount_ld,2)
            endif
            ! Geometric staircase ownership: lower-coordinate cell, never ID priority.
            ! A vacuum endpoint retains the legacy single-material boundary.
            if(fs%imedia(ix,iy+1,iz)/=0) then
              fe%c1_ey_z(ix,iy,iz) =  c1_e; fe%c2_ey_z(ix,iy,iz) =  c2_e_z;
              fe%c1_ey_x(ix,iy,iz) =  c1_e; fe%c2_ey_x(ix,iy,iz) = -c2_e_x;
              fe%c2_jy(ix,iy,iz)   = -c2_j;
              select case(media_type(ii))
              case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
                do ij=1,pole_num_ld(ii)
                  ld%c1_jy_ld(ix,iy,iz,ij) = tmp_c1_ld(ij,icount_ld)
                  ld%c2_jy_ld(ix,iy,iz,ij) = tmp_c2_ld(ij,icount_ld)
                  if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor')) &
                    ld%c2_jy_ld(ix,iy,iz,ij)=diagonal_c2(ij,icount_ld,2)
                  ld%c3_jy_ld(ix,iy,iz,ij) = tmp_c3_ld(ij,icount_ld)
                end do
              end select
            end if
            
            !------ ez and jz --------------------------------------------------!
            if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor'))then
              c1_e=diagonal_c5(icount_ld,3)/diagonal_c4(icount_ld,3)
              c2_e_x=cspeed_au*dt_em/diagonal_c4(icount_ld,3)/fs%hgs(1)
              c2_e_y=cspeed_au*dt_em/diagonal_c4(icount_ld,3)/fs%hgs(2)
              c2_e_z=cspeed_au*dt_em/diagonal_c4(icount_ld,3)/fs%hgs(3)
              c2_j=4d0*pi*dt_em/diagonal_c4(icount_ld,3)
            endif
            ! Geometric staircase ownership: lower-coordinate cell, never ID priority.
            ! A vacuum endpoint retains the legacy single-material boundary.
            if(fs%imedia(ix,iy,iz+1)/=0) then
              fe%c1_ez_x(ix,iy,iz) =  c1_e; fe%c2_ez_x(ix,iy,iz) =  c2_e_x;
              fe%c1_ez_y(ix,iy,iz) =  c1_e; fe%c2_ez_y(ix,iy,iz) = -c2_e_y;
              fe%c2_jz(ix,iy,iz)   = -c2_j;
              select case(media_type(ii))
              case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
                do ij=1,pole_num_ld(ii)
                  ld%c1_jz_ld(ix,iy,iz,ij) = tmp_c1_ld(ij,icount_ld)
                  ld%c2_jz_ld(ix,iy,iz,ij) = tmp_c2_ld(ij,icount_ld)
                  if((media_type(ii)=='ld-diagonal'.or.media_type(ii)=='ld-tensor')) &
                    ld%c2_jz_ld(ix,iy,iz,ij)=diagonal_c2(ij,icount_ld,3)
                  ld%c3_jz_ld(ix,iy,iz,ij) = tmp_c3_ld(ij,icount_ld)
                end do
              end select
            end if
            
          end if
        end do
        end do
        end do
      end if
      
      !update LD counter
      select case(media_type(ii))
      case('lorentz-drude','ld-diagonal','ld-tensor','ld-subpixel')
        icount_ld=icount_ld+1
      end select
      
      return
    end subroutine eh_coeff
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! H sites touch four material cells. Choose the first non-vacuum cell
    ! in a fixed geometrical order. Do not let media registration order overwrite mu.
    ! For a single material this assigns exactly the same legacy coefficients.
    subroutine eh_coeff_h_geometry
      implicit none
      integer :: i,j,k,n,m,axis,ids(4)
      real(8) :: cx,cy,cz
      do k=fs%mg%is(3),fs%mg%ie(3)
      do j=fs%mg%is(2),fs%mg%ie(2)
      do i=fs%mg%is(1),fs%mg%ie(1)
        do axis=1,3
          select case(axis)
          case(1)
            ids=(/fs%imedia(i,j,k),fs%imedia(i,j+1,k), &
                  fs%imedia(i,j,k+1),fs%imedia(i,j+1,k+1)/)
          case(2)
            ids=(/fs%imedia(i,j,k),fs%imedia(i+1,j,k), &
                  fs%imedia(i,j,k+1),fs%imedia(i+1,j,k+1)/)
          case(3)
            ids=(/fs%imedia(i,j,k),fs%imedia(i+1,j,k), &
                  fs%imedia(i,j+1,k),fs%imedia(i+1,j+1,k)/)
          end select
          m=0
          do n=1,4
            if(ids(n)/=0) then
              m=ids(n)
              exit
            end if
          end do
          if(m==0) cycle
          cx=cspeed_au/fe%rmu(m)*dt_em/fs%hgs(1)
          cy=cspeed_au/fe%rmu(m)*dt_em/fs%hgs(2)
          cz=cspeed_au/fe%rmu(m)*dt_em/fs%hgs(3)
          select case(axis)
          case(1)
            fe%c1_hx_y(i,j,k)=1.0d0; fe%c2_hx_y(i,j,k)=-cy
            fe%c1_hx_z(i,j,k)=1.0d0; fe%c2_hx_z(i,j,k)= cz
          case(2)
            fe%c1_hy_z(i,j,k)=1.0d0; fe%c2_hy_z(i,j,k)=-cz
            fe%c1_hy_x(i,j,k)=1.0d0; fe%c2_hy_x(i,j,k)= cx
          case(3)
            fe%c1_hz_x(i,j,k)=1.0d0; fe%c2_hz_x(i,j,k)=-cx
            fe%c1_hz_y(i,j,k)=1.0d0; fe%c2_hz_y(i,j,k)= cy
          end select
        end do
      end do
      end do
      end do
    end subroutine eh_coeff_h_geometry

    !+ check media +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_check_media_type(rep,rmu,sig,cm_type,cm_name)
      implicit none
      real(8),      intent(in)  :: rep,rmu,sig
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
          if     (ipulse==1) then
            write(tmp_c(1),*) 'For ek_dir1(1) = 1.0d0, epdir_re1(1) and epdir_im1(1) must be 0.0d0.'
          elseif (ipulse==2) then
            write(tmp_c(1),*) 'For ek_dir2(1) = 1.0d0, epdir_re2(1) and epdir_im2(1) must be 0.0d0.'
          end if
          call stop_em('Invalid input keywords:',c2=trim(adjustl(tmp_c(1))))
        else
          inc_dist='yz-plane'
        end if
      elseif(ek_dir(1)==0.0d0.and.ek_dir(2)==1.0d0.and.ek_dir(3)==0.0d0) then !y-direction propagation
        if(epdir_re(2)/=0.0d0.or.epdir_im(2)/=0.0d0) then
          if     (ipulse==1) then
            write(tmp_c(1),*) 'For ek_dir1(2) = 1.0d0, epdir_re1(2) and epdir_im1(2) must be 0.0d0.'
          elseif (ipulse==2) then
            write(tmp_c(1),*) 'For ek_dir2(2) = 1.0d0, epdir_re2(2) and epdir_im2(2) must be 0.0d0.'
          end if
          call stop_em('Invalid input keywords:',c2=trim(adjustl(tmp_c(1))))
        else
          inc_dist='xz-plane'
        end if
      elseif(ek_dir(1)==0.0d0.and.ek_dir(2)==0.0d0.and.ek_dir(3)==1.0d0) then !z-direction propagation
        if(epdir_re(3)/=0.0d0.or.epdir_im(3)/=0.0d0) then
          if     (ipulse==1) then
            write(tmp_c(1),*) 'For ek_dir1(3) = 1.0d0, epdir_re1(3) and epdir_im1(3) must be 0.0d0.'
          elseif (ipulse==2) then
            write(tmp_c(1),*) 'For ek_dir2(3) = 1.0d0, epdir_re2(3) and epdir_im2(3) must be 0.0d0.'
          end if
          call stop_em('Invalid input keywords:',c2=trim(adjustl(tmp_c(1))))
        else
          inc_dist='xy-plane'
        end if
      else
        if     (ipulse==1) then
          write(tmp_c(1),*) 'ek_dir1 is only allowed by (0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0).'
        elseif (ipulse==2) then
          write(tmp_c(1),*) 'ek_dir2 is only allowed by (0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0).'
        end if
        call stop_em('Invalid input keywords:',c2=trim(adjustl(tmp_c(1))))
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
        if     (ipulse==1) then
          call stop_em('set ae_shape1 to "Ecos2" or "Acos2".')
        elseif (ipulse==2) then
          call stop_em('set ae_shape2 to "Ecos2" or "Acos2".')
        end if
      end select
      
      !check phi_cep
      diff_cep=(phi_cep-0.25d0)*2.d0-int((phi_cep-0.25d0)*2.d0)
      if(ae_shape=='Ecos2'.and.abs(diff_cep)>=1.d-12)then
        if     (ipulse==1) then
          call stop_em('phi_cep1 must be equal to 0.25+0.5*i when "Ecos2" is specified for ae_shape1.')
        elseif (ipulse==2) then
          call stop_em('phi_cep2 must be equal to 0.25+0.5*i when "Ecos2" is specified for ae_shape2.')
        end if
      end if
      
      !set E_amplitude
      if(I_wcm2/=-1d0) &
        E_amplitude=sqrt(I_wcm2)*1.0d2*2.74492d1/(5.14223d11) !I[W/cm^2]->E[a.u.]
        
      return
    end subroutine eh_check_iw_parameter
    
    !+ CONTAINED IN eh_init ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ make gauss function used for beam spot ++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_make_gauss(rloc,s_p,s_l,g_xy,g_yz,g_xz,g_x,g_y,g_z)
      implicit none
      real(8), intent(in)  :: rloc(3),s_p(3),s_l(3)
      real(8), intent(out) :: g_xy(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(2):fs%mg%ie_array(2)),&
                              g_yz(fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3)),&
                              g_xz(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8), intent(out) :: g_x( fs%mg%is_array(1):fs%mg%ie_array(1) ),&
                              g_y( fs%mg%is_array(2):fs%mg%ie_array(2) ),&
                              g_z( fs%mg%is_array(3):fs%mg%ie_array(3) )
      integer :: i1s,i2s,i1,i2
      real(8) :: rtmp,gtmp
      
      !make 2D(plane) gauss function
      do ii=1,3
        if    (ii==1) then !xy
          i1s=1; i2s=2;
        elseif(ii==2) then !yz
          i1s=2; i2s=3;
        elseif(ii==3) then !xz
          i1s=1; i2s=3;
        end if
        if(s_p(ii)>0.0d0) then
          do i2=fs%mg%is_array(i2s),fs%mg%ie_array(i2s)
          do i1=fs%mg%is_array(i1s),fs%mg%ie_array(i1s)
            rtmp = sqrt( (fe%coo(i1,i1s)-rloc(i1s))**2.0d0 &
                        +(fe%coo(i2,i2s)-rloc(i2s))**2.0d0 )
            gtmp = exp( -0.5d0*((rtmp/s_p(ii))**2.0d0) )
            if    (ii==1) then !xy
              g_xy(i1,i2) = gtmp
            elseif(ii==2) then !yz
              g_yz(i1,i2) = gtmp
            elseif(ii==3) then !xz
              g_xz(i1,i2) = gtmp
            end if
          end do
          end do
        end if
      end do
      
      !make 1D(line) gauss function
      do ii=1,3
        if    (ii==1) then !x
          i1s=1;
        elseif(ii==2) then !y
          i1s=2;
        elseif(ii==3) then !z
          i1s=3;
        end if
        if(s_l(ii)>0.0d0) then
          do i1=fs%mg%is_array(i1s),fs%mg%ie_array(i1s)
            rtmp = sqrt( (fe%coo(i1,i1s)-rloc(i1s))**2.0d0 )
            gtmp = exp( -0.5d0*((rtmp/s_l(ii))**2.0d0) )
            if    (ii==1) then !x
              g_x(i1) = gtmp
            elseif(ii==2) then !y
              g_y(i1) = gtmp
            elseif(ii==3) then !z
              g_z(i1) = gtmp
            end if
          end do
        end if
      end do
      
      return
    end subroutine eh_make_gauss
    
  end subroutine eh_init
  
  !===========================================================================================
  != calculate eh-FDTD =======================================================================
  ! Adjoint paired edge->node and node->edge means. Projection is by material
  ! ownership, never by numeric material-ID priority. This leaves diagonal paths exact.
  subroutine eh_tensor_off(fs,fe,pole,v,out)
    use structures,only:s_fdtd_system
    use sendrecv_grid,only:update_overlap_real8
    implicit none
    type(s_fdtd_system),intent(inout)::fs
    type(ls_fdtd_eh),intent(in)::fe
    integer,intent(in)::pole
    real(8),intent(in)::v(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
      fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
    real(8),intent(out)::out(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
      fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
    real(8),allocatable::edge(:,:,:,:),node(:,:,:,:)
    real(8)::mean(3)
    integer::lo(3),hi(3),m,a,ix,iy,iz,next(3),nmodes
    lo=fs%mg%is_array;hi=fs%mg%ie_array;out=0d0
    if(fe%subpixel_enabled)then
      if(.not.fe%subpixel_cross(pole))return
    endif
    allocate(edge(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
      node(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    nmodes=ubound(fe%tensor_corr,4)
    if(fe%subpixel_enabled)nmodes=1
    do m=1,nmodes
      if(.not.fe%subpixel_enabled)then
        if(all(fe%tensor_corr(:,:,pole,m)==0d0))cycle
      endif
      edge=0d0;node=0d0
      where(fe%tensor_owner==m)edge=fe%tensor_weight(:,:,:,:,pole)*v
      do a=1,3
        call update_overlap_real8(fs%srg_ng,fs%mg,edge(a,:,:,:))
      enddo
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        mean=[.5d0*(edge(1,ix,iy,iz)+edge(1,ix-1,iy,iz)), &
              .5d0*(edge(2,ix,iy,iz)+edge(2,ix,iy-1,iz)), &
              .5d0*(edge(3,ix,iy,iz)+edge(3,ix,iy,iz-1))]
        if(fe%subpixel_enabled)then
          node(:,ix,iy,iz)=matmul(fe%subpixel_corr(:,:,ix,iy,iz,pole),mean)
        else
          node(:,ix,iy,iz)=matmul(fe%tensor_corr(:,:,pole,m),mean)
        endif
      enddo
      enddo
      enddo
      do a=1,3
        call update_overlap_real8(fs%srg_ng,fs%mg,node(a,:,:,:))
      enddo
      do iz=fs%mg%is(3),fs%mg%ie(3)
      do iy=fs%mg%is(2),fs%mg%ie(2)
      do ix=fs%mg%is(1),fs%mg%ie(1)
        do a=1,3
          if(fe%tensor_owner(a,ix,iy,iz)/=m)cycle
          next=[ix,iy,iz];next(a)=next(a)+1
          out(a,ix,iy,iz)=out(a,ix,iy,iz)+.5d0*fe%tensor_weight(a,ix,iy,iz,pole) &
            *(node(a,ix,iy,iz)+node(a,next(1),next(2),next(3)))
        enddo
      enddo
      enddo
      enddo
    enddo
  end subroutine eh_tensor_off

  subroutine eh_calc(fs,fe,ld,external_current,step_guard)
    use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
    use salmon_global,   only: dt_em,obs_num_em,obs_samp_em,yn_obs_plane_em,yn_obs_plane_integral_em,&
                               base_directory,nt_em,yn_periodic,ek_dir1,t1_t2,t1_start,              &
                               E_amplitude1,tw1,omega1,phi_cep1,epdir_re1,epdir_im1,ae_shape1,       &
                               E_amplitude2,tw2,omega2,phi_cep2,epdir_re2,epdir_im2,ae_shape2,       &
                               checkpoint_interval,ase_num_em,art_num_em
    use inputoutput,     only: utime_from_au,uenergy_from_au
    use parallelization, only: nproc_id_global,nproc_size_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_fdtd_system
    use misc_routines,   only: get_wtime
    use common_maxwell,  only: output_r_txt_em,output_r_bin_em,txtfile_copy_em,stop_em
    implicit none
    integer :: restart_io,restart_bad,restart_total
    type(s_ld_state), intent(inout) :: ld
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    ! Optional physical charge current on native E edges. No effect on CLI runs.
    real(8),optional,intent(in) :: external_current(3, &
      fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(2):fs%mg%ie_array(2), &
      fs%mg%is_array(3):fs%mg%ie_array(3))
    ! Optional native transaction gate. Batch calls retain the original path.
    interface
      logical function step_guard(stage)
        integer,intent(in) :: stage
      end function
    end interface
    optional :: step_guard
    integer :: observation_unit
    integer            :: iter,ii,ix,iy,iz
    integer            :: is5(5),ie5(5)
    real(8)            :: elapsed_time
    character(256)     :: save_name,tmp_name
    real(8),allocatable::tensor_drive(:,:,:,:,:)
    if(fe%tensor_enabled)then
      allocate(tensor_drive(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
        fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3),ld%max_pole_num_ld))
      tensor_drive=0d0
    endif
    !time-iteration
    do iter=fe%iter_sta,fe%iter_end
      !update iter_now
      fe%iter_now=iter
      if( comm_is_root(nproc_id_global) .and. mod(iter,10)==0 )then
        write(*,'(I12)') fe%iter_now
      end if
      
      !sotre old variables
      call ld_store_old(ld,fs%mg%is_array,fs%mg%ie_array,fs%mg%is,fs%mg%ie,dt_em, &
        fe%ex_y,fe%ex_z,fe%ey_z,fe%ey_x,fe%ez_x,fe%ez_y)
      
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
        if(fe%inc_dist1/='none') call eh_add_inc(1,E_amplitude1,tw1,omega1,phi_cep1,                       &
                                                 epdir_re1,epdir_im1,fe%c2_inc1_xyz,ae_shape1,fe%inc_dist1,&
                                                 fe%gbeam_width_xy1,fe%gbeam_width_yz1,fe%gbeam_width_xz1, &
                                                 fe%gbeam_width_x1 ,fe%gbeam_width_y1 ,fe%gbeam_width_z1 )
        if(fe%inc_dist2/='none') call eh_add_inc(2,E_amplitude2,tw2,omega2,phi_cep2,                       &
                                                 epdir_re2,epdir_im2,fe%c2_inc2_xyz,ae_shape2,fe%inc_dist2,&
                                                 fe%gbeam_width_xy2,fe%gbeam_width_yz2,fe%gbeam_width_xz2, &
                                                 fe%gbeam_width_x2 ,fe%gbeam_width_y2 ,fe%gbeam_width_z2 )
      end if
      if(ld%enabled) call eh_add_curr(ld%rjx_fdtd_ld(:,:,:),ld%rjy_fdtd_ld(:,:,:),ld%rjz_fdtd_ld(:,:,:))
      if(present(external_current)) then
        if(any(external_current/=0d0)) call eh_add_curr(external_current(1,:,:,:), &
          external_current(2,:,:,:),external_current(3,:,:,:))
      end if
      if(fe%tensor_enabled)then
        call eh_tensor_solve
        if(fe%tensor_failed)then
          if(present(step_guard))then
            if(.not.step_guard(20))return
          endif
          call stop_em('Nonfinite or non-positive tensor electric-field solve.')
        endif
      endif
      call eh_sendrecv(fs,fe,'e')
      if(present(step_guard))then
        if(.not.step_guard(20))return
      endif
      
      !update lorentz-drude
      if(fe%tensor_enabled)then
        call ld_update(ld,fs%mg%is_array,fs%mg%ie_array,fs%mg%is,fs%mg%ie,dt_em, &
          fe%ex_y,fe%ex_z,fe%ey_z,fe%ey_x,fe%ez_x,fe%ez_y,tensor_drive)
      else
      if(ld%enabled) call ld_update(ld,fs%mg%is_array,fs%mg%ie_array,fs%mg%is,fs%mg%ie,dt_em, &
        fe%ex_y,fe%ex_z,fe%ey_z,fe%ey_x,fe%ez_x,fe%ez_y)
      endif
      
      if(present(step_guard))then
        if(.not.step_guard(21))return
      endif

      !calculate linear response
      if(fe%flag_lr) call eh_calc_lr
      
      !store old h
      if( fe%flag_save .and. mod(iter,obs_samp_em)==0 )then
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2)
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
      if(present(step_guard))then
        if(.not.step_guard(22))return
        if(.not.step_guard(23))return
      endif
      


      !output
      if( fe%flag_save .and. mod(iter,obs_samp_em)==0 )then
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
        call eh_planes_write(fs,fe,iter)
        
        !observation point and plane data
        if(fe%flag_obs) then
          do ii=1,obs_num_em
            !point
            if(fe%iobs_po_pe(ii)==1) then
              write(save_name,*) ii
              save_name=trim(adjustl(base_directory))//'/obs'//trim(adjustl(save_name))//'_at_point_rt.data'
              call eh_point_output(fe,save_name,.false.,dble(iter)*dt_em*utime_from_au, &
                    [ &
                    fe%ex_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                    fe%ey_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                    fe%ez_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uVperm_from_au, &
                    fe%hx_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au, &
                    fe%hy_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au, &
                    fe%hz_s(fe%iobs_po_id(ii,1),fe%iobs_po_id(ii,2),fe%iobs_po_id(ii,3))*fe%uAperm_from_au ])

            end if
            
            !plane
            if(yn_obs_plane_em(ii)=='y') then
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ex_s,'ex', &
                                 fe%snapshot_io_warned)
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ey_s,'ey', &
                                 fe%snapshot_io_warned)
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%ez_s,'ez', &
                                 fe%snapshot_io_warned)
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hx_s,'hx', &
                                 fe%snapshot_io_warned)
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hy_s,'hy', &
                                 fe%snapshot_io_warned)
              call eh_save_plane(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uAperm_from_au,&
                                 fs%mg%is,fs%mg%ie,fs%lg%is,fs%lg%ie,fe%Nd,fe%ifn,ii,iter,fe%hz_s,'hz', &
                                 fe%snapshot_io_warned)
            end if
            
            !plane integral
            if(yn_obs_plane_integral_em(ii)=='y') then
              call eh_save_plane_integral(fe%iobs_po_id(ii,:),fe%iobs_pl_pe(ii,:),fe%uVperm_from_au,fe%uAperm_from_au, &
                                          dble(iter)*dt_em*utime_from_au,fs%mg%is,fs%mg%ie,fs%hgs,fe%Nd,fe%ifn,ii,     &
                                          fe%ex_s,fe%ey_s,fe%ez_s,fe%hx_s,fe%hy_s,fe%hz_s, &
                                          fe%plane_integral_io_warned)
            end if
            
            !obs_plane_ene_em option(spatial distribution at each energy point)
            if(fe%iobs_num_ene(ii)>0) call eh_calc_plane_ene(fs,fe,ii,iter)
          end do
          
          !check maximum
          call eh_update_max
        end if
        
        !ase option(calculate incident and scattering fields on closed surface[box shape])
        if(fe%flag_ase) call eh_calc_inc_sca_ase
        
        !art option(calculate incident and total fields on plane)
        if(fe%flag_art) call eh_calc_inc_tot_art
        
        !output data_for_restart
        if(fe%flag_o_restart .and. mod(iter,checkpoint_interval)==0) then
          !start to calculate elapsed_time
          elapsed_time = get_wtime()
          
          restart_io=0
          !iter_now, e_max and h_max
          call output_r_txt_em(fe%ifn,1,1,'single','iter_now',i0d=fe%iter_now, &
            io_status=restart_io)
          call output_r_txt_em(fe%ifn,1,1,'single',   'e_max',r0d=fe%e_max   , &
            io_status=restart_io)
          call output_r_txt_em(fe%ifn,1,1,'single',   'h_max',r0d=fe%h_max   , &
            io_status=restart_io)
          
          !bin data
          is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
          is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
          is5(3) = fs%mg%is(3); ie5(3) = fs%mg%ie(3);
          is5(4) = 1          ; ie5(4) = 1          ;
          is5(5) = 1          ; ie5(5) = 1          ;
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ex_y',r3d=fe%ex_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ex_z',r3d=fe%ex_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ey_z',r3d=fe%ey_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ey_x',r3d=fe%ey_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ez_x',r3d=fe%ez_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','ez_y',r3d=fe%ez_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hx_y',r3d=fe%hx_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hx_z',r3d=fe%hx_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hy_z',r3d=fe%hy_z(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hy_x',r3d=fe%hy_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hz_x',r3d=fe%hz_x(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          call output_r_bin_em(fe%ifn,is5,ie5,'all','hz_y',r3d=fe%hz_y(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          if(ld%enabled) then
            is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1)       ;
            is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2)       ;
            is5(3) = fs%mg%is(3); ie5(3) = fs%mg%ie(3)       ;
            is5(4) = 1          ; ie5(4) = ld%max_pole_num_ld;
            is5(5) = 1          ; ie5(5) = 1                 ;
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjx_ld',&
                                 r4d=ld%rjx_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjy_ld',&
                                 r4d=ld%rjy_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjz_ld',&
                                 r4d=ld%rjz_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','px_ld' ,&
                                 r4d= ld%px_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','py_ld' ,&
                                 r4d= ld%py_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','pz_ld' ,&
                                 r4d= ld%pz_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3),is5(4):ie5(4)), &
            io_status=restart_io)
            is5(4) = 1; ie5(4) = 1;
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjx_fdtd_ld',&
                                 r3d=ld%rjx_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjy_fdtd_ld',&
                                 r3d=ld%rjy_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','rjz_fdtd_ld',&
                                 r3d=ld%rjz_fdtd_ld(is5(1):ie5(1),is5(2):ie5(2),is5(3):ie5(3)), &
            io_status=restart_io)
          end if
          if(fe%flag_lr) then
            is5(1) = 1; ie5(1) = nt_em;
            is5(2) = 1; ie5(2) = 1    ;
            is5(3) = 1; ie5(3) = 1    ;
            is5(4) = 1; ie5(4) = 1    ;
            is5(5) = 1; ie5(5) = 1    ;
            call output_r_bin_em(fe%ifn,is5,ie5,'single','time_lr',r1d=fe%time_lr, &
            io_status=restart_io)
            is5(2) = 1; ie5(2) = 3    ;
            if(yn_periodic=='n') then
              call output_r_bin_em(fe%ifn,is5,ie5,'single','dip_lr' ,r2d=fe%dip_lr , &
            io_status=restart_io)
            elseif(yn_periodic=='y') then
              call output_r_bin_em(fe%ifn,is5,ie5,'single','curr_lr',r2d=fe%curr_lr, &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'single','e_lr'   ,r2d=fe%e_lr   , &
            io_status=restart_io)
            end if
          end if
          if(fe%flag_obs) then
          if(sum(fe%iobs_num_ene(:))>0) then
            is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1)               ;
            is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2)               ;
            is5(3) = 1          ; ie5(3) = 1                         ;
            is5(4) = 1          ; ie5(4) = maxval(fe%iobs_num_ene(:));
            is5(5) = 1          ; ie5(5) = 2                         ;
            do ii=1,obs_num_em
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ex_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ex_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ey_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ey_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ez_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ez_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hx_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hx_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hy_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hy_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hz_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hz_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jx_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jx_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jy_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jy_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jz_xy_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jz_xy_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,1), &
            io_status=restart_io)
            end do
            is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
            is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
            do ii=1,obs_num_em
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ex_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ex_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ey_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ey_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ez_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ez_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hx_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hx_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hy_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hy_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hz_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hz_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jx_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jx_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jy_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jy_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jz_yz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jz_yz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,2), &
            io_status=restart_io)
            end do
            is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
            is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
            do ii=1,obs_num_em
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ex_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ex_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ey_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ey_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_ez_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_ez_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hx_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hx_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hy_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hy_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_hz_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_hz_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jx_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jx_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jy_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jy_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','obs_jz_xz_ene_iobs'//trim(adjustl(save_name)),&
                                   c5d=fe%obs_jz_xz_ene(:,:,ii,:,:),ipe=fe%iobs_pl_pe(ii,3), &
            io_status=restart_io)
            end do
          end if
          end if
          if(fe%flag_ase) then
            if    (ek_dir1(1)==1.0d0) then
              is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
            elseif(ek_dir1(2)==1.0d0) then
              is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
            elseif(ek_dir1(3)==1.0d0) then
              is5(1) = fs%mg%is(3); ie5(1) = fs%mg%ie(3);
            end if
            is5(2) = 1; ie5(2) = ase_num_em;
            is5(3) = 1; ie5(3) = 2         ;
            is5(4) = 1; ie5(4) = 1         ;
            is5(5) = 1; ie5(5) = 1         ;
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ex_inc_pa_ene',c3d=fe%ase_ex_inc_pa_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ey_inc_pa_ene',c3d=fe%ase_ey_inc_pa_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ez_inc_pa_ene',c3d=fe%ase_ez_inc_pa_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hx_inc_pa_ene',c3d=fe%ase_hx_inc_pa_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hy_inc_pa_ene',c3d=fe%ase_hy_inc_pa_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hz_inc_pa_ene',c3d=fe%ase_hz_inc_pa_ene, &
            io_status=restart_io)
            is5(1) = 1; ie5(1) = ase_num_em;
            is5(2) = 1; ie5(2) = 2         ;
            is5(3) = 1; ie5(3) = 2         ;
            is5(4) = 1; ie5(4) = 1         ;
            is5(5) = 1; ie5(5) = 1         ;
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_ex_inc_bt_ene',c3d=fe%ase_ex_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_ey_inc_bt_ene',c3d=fe%ase_ey_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_ez_inc_bt_ene',c3d=fe%ase_ez_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_hx_inc_bt_ene',c3d=fe%ase_hx_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_hy_inc_bt_ene',c3d=fe%ase_hy_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','ase_hz_inc_bt_ene',c3d=fe%ase_hz_inc_bt_ene, &
            io_status=restart_io)
            is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
            is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
            is5(3) = 1          ; ie5(3) = ase_num_em ;
            is5(4) = 1          ; ie5(4) = 1          ;
            is5(5) = 1          ; ie5(5) = 2          ;
            do ii = 1,2
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ex_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ex_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ey_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ey_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hx_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hx_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hy_xy_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hy_xy_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(3,ii), &
            io_status=restart_io)
            end do
            is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
            is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
            do ii = 1,2
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ey_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ey_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ez_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ez_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hy_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hy_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hz_yz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hz_yz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(1,ii), &
            io_status=restart_io)
            end do
            is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
            is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
            do ii = 1,2
              write(save_name,*) ii
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ex_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ex_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_ez_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_ez_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hx_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hx_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
              call output_r_bin_em(fe%ifn,is5,ie5,'all','ase_hz_xz_sca_ene_ibt'//trim(adjustl(save_name)),&
                                   c5d=fe%ase_hz_xz_sca_ene(:,:,:,ii,:),ipe=fe%ipe_on_box_surf(2,ii), &
            io_status=restart_io)
            end do
          end if
          if(fe%flag_art) then
            is5(1) = 1; ie5(1) = art_num_em;
            is5(2) = 1; ie5(2) = 2         ;
            is5(3) = 1; ie5(3) = 2         ;
            is5(4) = 1; ie5(4) = 1         ;
            is5(5) = 1; ie5(5) = 1         ;
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_ex_inc_bt_ene',c3d=fe%art_ex_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_ey_inc_bt_ene',c3d=fe%art_ey_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_ez_inc_bt_ene',c3d=fe%art_ez_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_hx_inc_bt_ene',c3d=fe%art_hx_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_hy_inc_bt_ene',c3d=fe%art_hy_inc_bt_ene, &
            io_status=restart_io)
            call output_r_bin_em(fe%ifn,is5,ie5,'single','art_hz_inc_bt_ene',c3d=fe%art_hz_inc_bt_ene, &
            io_status=restart_io)
            if    (ek_dir1(1)==1.0d0) then
              is5(1) = fs%mg%is(2); ie5(1) = fs%mg%ie(2);
              is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
              is5(3) = 1          ; ie5(3) = art_num_em ;
              is5(4) = 1          ; ie5(4) = 1          ;
              is5(5) = 1          ; ie5(5) = 2          ;
              do ii = 1,2
                write(save_name,*) ii
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ey_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ey_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ez_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ez_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hy_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hy_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hz_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hz_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
              end do
            elseif(ek_dir1(2)==1.0d0) then
              is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
              is5(2) = fs%mg%is(3); ie5(2) = fs%mg%ie(3);
              is5(3) = 1          ; ie5(3) = art_num_em ;
              is5(4) = 1          ; ie5(4) = 1          ;
              is5(5) = 1          ; ie5(5) = 2          ;
              do ii = 1,2
                write(save_name,*) ii
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ex_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ex_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ez_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ez_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hx_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hx_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hz_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hz_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
              end do
            elseif(ek_dir1(3)==1.0d0) then
              is5(1) = fs%mg%is(1); ie5(1) = fs%mg%ie(1);
              is5(2) = fs%mg%is(2); ie5(2) = fs%mg%ie(2);
              is5(3) = 1          ; ie5(3) = art_num_em ;
              is5(4) = 1          ; ie5(4) = 1          ;
              is5(5) = 1          ; ie5(5) = 2          ;
              do ii = 1,2
                write(save_name,*) ii
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ex_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ex_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_ey_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_ey_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hx_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hx_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
                call output_r_bin_em(fe%ifn,is5,ie5,'all','art_hy_tot_bt_ene_ibt'//trim(adjustl(save_name)),&
                                     c5d=fe%art_hy_tot_bt_ene(:,:,:,ii,:),ipe=fe%ipe_on_plane(ii), &
            io_status=restart_io)
              end do
            end if            
          end if
          
          restart_bad=0
          if(restart_io/=0)restart_bad=1
          call comm_summation(restart_bad,restart_total,nproc_group_global)
          if(restart_total/=0)error stop 'Required EM restart state save failed; incomplete checkpoint'

          !copy obsx_at_point_rt.data and obsx_plane_integral_rt.data
          if(fe%flag_obs .and. comm_is_root(nproc_id_global)) then
            do ii=1,obs_num_em
              write(tmp_name,*) ii
              save_name = trim(adjustl(base_directory))//'/obs'//trim(adjustl(tmp_name))//'_at_point_rt.data'
              tmp_name  =                'data_for_restart/obs'//trim(adjustl(tmp_name))//'_at_point_rt.data'
              call txtfile_copy_em(trim(adjustl(save_name)),trim(adjustl(tmp_name)))
              if(yn_obs_plane_integral_em(ii)=='y') then
                write(tmp_name,*) ii
                save_name = trim(adjustl(base_directory))//'/obs'//trim(adjustl(tmp_name))//'_xy_integral_rt.data'
                tmp_name  =                'data_for_restart/obs'//trim(adjustl(tmp_name))//'_xy_integral_rt.data'
                call txtfile_copy_em(trim(adjustl(save_name)),trim(adjustl(tmp_name)))
                write(tmp_name,*) ii
                save_name = trim(adjustl(base_directory))//'/obs'//trim(adjustl(tmp_name))//'_yz_integral_rt.data'
                tmp_name  =                'data_for_restart/obs'//trim(adjustl(tmp_name))//'_yz_integral_rt.data'
                call txtfile_copy_em(trim(adjustl(save_name)),trim(adjustl(tmp_name)))
                write(tmp_name,*) ii
                save_name = trim(adjustl(base_directory))//'/obs'//trim(adjustl(tmp_name))//'_xz_integral_rt.data'
                tmp_name  =                'data_for_restart/obs'//trim(adjustl(tmp_name))//'_xz_integral_rt.data'
                call txtfile_copy_em(trim(adjustl(save_name)),trim(adjustl(tmp_name)))
              end if
            end do
          end if
          
          !write elapsed time
          elapsed_time = get_wtime() - elapsed_time
          if(comm_is_root(nproc_id_global)) then
            write(*,'(A,f16.8)') " elapsed time for outputing restart data [s] = ", elapsed_time
          end if
        end if
      end if
    end do
    
    return
  contains
    subroutine eh_tensor_solve
      real(8),allocatable::pred(:,:,:,:),old(:,:,:,:),b(:,:,:,:),x(:,:,:,:),res(:,:,:,:)
      real(8),allocatable::z(:,:,:,:),direction(:,:,:,:),ad(:,:,:,:),temp(:,:,:,:)
      real(8)::rho,next_rho,den,alpha,beta,normb,normr
      integer::lo(3),hi(3),q,k
      lo=fs%mg%is_array;hi=fs%mg%ie_array
      allocate(pred(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      pred(1,:,:,:)=fe%ex_y+fe%ex_z
      pred(2,:,:,:)=fe%ey_z+fe%ey_x
      pred(3,:,:,:)=fe%ez_x+fe%ez_y
      old=pred
      old(1,:,:,:)=ld%ex_old_ld;old(2,:,:,:)=ld%ey_old_ld;old(3,:,:,:)=ld%ez_old_ld
      b=fe%tensor_diag*pred;x=pred;res=pred;z=pred;direction=pred;ad=pred;temp=pred
      call eh_tensor_off(fs,fe,0,old,temp)
      b=b+temp
      do q=1,ld%max_pole_num_ld
        call eh_tensor_off(fs,fe,q,old,temp)
        b=b-2d0*acos(-1d0)*dt_em*temp
      enddo
      call tensor_product(x,ad)
      res=b-ad;z=res/fe%tensor_diag;direction=z
      rho=tensor_dot(res,z);normb=sqrt(tensor_dot(b,b))
      fe%tensor_failed=.false.
      do k=1,fe%tensor_maxiter
        normr=sqrt(tensor_dot(res,res))
        if(.not.ieee_is_finite(normr).or..not.ieee_is_finite(rho).or..not.ieee_is_finite(normb))then
          fe%tensor_failed=.true.;exit
        endif
        if(normr<=fe%tensor_rtol*normb)exit
        call tensor_product(direction,ad)
        den=tensor_dot(direction,ad)
        if(.not.ieee_is_finite(den).or.den<=0d0.or.rho<=0d0)then
          fe%tensor_failed=.true.;exit
        endif
        alpha=rho/den;x=x+alpha*direction;res=res-alpha*ad
        z=res/fe%tensor_diag;next_rho=tensor_dot(res,z)
        beta=next_rho/rho;direction=z+beta*direction;rho=next_rho
      enddo
      if(.not.fe%tensor_failed)then
        normr=sqrt(tensor_dot(res,res))
        if(normr>fe%tensor_rtol*normb.and.comm_is_root(nproc_id_global)) &
          write(*,'(A,I0,A,ES12.4)')'WARNING: tensor E solver budget at step ',iter, &
          '; relative residual=',normr/max(normb,tiny(1d0))
        fe%ex_y=fe%ex_y+.5d0*(x(1,:,:,:)-pred(1,:,:,:))
        fe%ex_z=fe%ex_z+.5d0*(x(1,:,:,:)-pred(1,:,:,:))
        fe%ey_z=fe%ey_z+.5d0*(x(2,:,:,:)-pred(2,:,:,:))
        fe%ey_x=fe%ey_x+.5d0*(x(2,:,:,:)-pred(2,:,:,:))
        fe%ez_x=fe%ez_x+.5d0*(x(3,:,:,:)-pred(3,:,:,:))
        fe%ez_y=fe%ez_y+.5d0*(x(3,:,:,:)-pred(3,:,:,:))
        old=old+x
        do q=1,ld%max_pole_num_ld
          call eh_tensor_off(fs,fe,q,old,temp)
          tensor_drive(:,:,:,:,q)=temp
        enddo
      endif
    end subroutine eh_tensor_solve

    subroutine tensor_product(v,out)
      real(8),intent(in)::v(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
        fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8),intent(out)::out(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
        fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8),allocatable::temp(:,:,:,:)
      integer::q
      allocate(temp,mold=v)
      call eh_tensor_off(fs,fe,0,v,temp)
      out=fe%tensor_diag*v+temp
      do q=1,ld%max_pole_num_ld
        call eh_tensor_off(fs,fe,q,v,temp)
        out=out+2d0*acos(-1d0)*dt_em*temp
      enddo
    end subroutine tensor_product

    real(8) function tensor_dot(a,b)result(value)
      real(8),intent(in)::a(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
        fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8),intent(in)::b(3,fs%mg%is_array(1):fs%mg%ie_array(1), &
        fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8)::local
      local=sum(a(:,fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3)) &
               *b(:,fs%mg%is(1):fs%mg%ie(1),fs%mg%is(2):fs%mg%ie(2),fs%mg%is(3):fs%mg%ie(3)))
      call comm_summation(local,value,nproc_group_global)
    end function tensor_dot


    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ store old variables +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ update lorentz-drude ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ calculate linear response +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_calc_lr
      use salmon_global, only: yn_periodic
      implicit none
      real(8) :: sum_lr(3),sum_lr2(3)
      
      !update time
      fe%time_lr(iter)=dble(iter)*dt_em
      
      if(yn_periodic=='n') then
        !calculate dipolemoment
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
        if(ld%enabled) then
          do ii=1,ld%max_pole_num_ld
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2) reduction( + : sum_lr )
            do iz=fs%mg%is(3),fs%mg%ie(3)
            do iy=fs%mg%is(2),fs%mg%ie(2)
            do ix=fs%mg%is(1),fs%mg%ie(1)
              sum_lr(1) = sum_lr(1) + ld%px_ld(ix,iy,iz,ii)
              sum_lr(2) = sum_lr(2) + ld%py_ld(ix,iy,iz,ii)
              sum_lr(3) = sum_lr(3) + ld%pz_ld(ix,iy,iz,ii)
            end do
            end do
            end do
!$omp end do
!$omp end parallel
          end do
        end if
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%dip_lr(iter,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3)
      elseif(yn_periodic=='y') then
        !calculate average current density
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
        if(ld%enabled) then
          do ii=1,ld%max_pole_num_ld
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2) reduction( + : sum_lr )
            do iz=fs%mg%is(3),fs%mg%ie(3)
            do iy=fs%mg%is(2),fs%mg%ie(2)
            do ix=fs%mg%is(1),fs%mg%ie(1)
              sum_lr(1) = sum_lr(1) + ld%rjx_ld(ix,iy,iz,ii)
              sum_lr(2) = sum_lr(2) + ld%rjy_ld(ix,iy,iz,ii)
              sum_lr(3) = sum_lr(3) + ld%rjz_ld(ix,iy,iz,ii)
            end do
            end do
            end do
!$omp end do
!$omp end parallel
          end do
        end if
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%curr_lr(iter,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3) &
                           /(fs%rlsize(1)*fs%rlsize(2)*fs%rlsize(3))
        
        !calculate average electric field
        sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) collapse(2) reduction( + : sum_lr )
        do iz=fs%mg%is(3),fs%mg%ie(3)
        do iy=fs%mg%is(2),fs%mg%ie(2)
        do ix=fs%mg%is(1),fs%mg%ie(1)
          sum_lr(1) = sum_lr(1) + ( fe%ex_y(ix,iy,iz)+fe%ex_z(ix,iy,iz) )
          sum_lr(2) = sum_lr(2) + ( fe%ey_z(ix,iy,iz)+fe%ey_x(ix,iy,iz) )
          sum_lr(3) = sum_lr(3) + ( fe%ez_x(ix,iy,iz)+fe%ez_y(ix,iy,iz) )
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
        fe%e_lr(iter,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3) &
                        /(fs%rlsize(1)*fs%rlsize(2)*fs%rlsize(3))
      end if
      
      return
    end subroutine eh_calc_lr
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ add incident current source +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine eh_add_inc(iord,amp,tw,omega,cep,ep_r,ep_i,c2_inc_xyz,aes,typ,&
                          g_xy,g_yz,g_xz,g_x,g_y,g_z)
      implicit none
      integer,      intent(in) :: iord
      real(8),      intent(in) :: amp,tw,omega,cep
      real(8),      intent(in) :: ep_r(3),ep_i(3),c2_inc_xyz(3)
      character(16),intent(in) :: aes,typ
      real(8),      intent(in) :: g_xy(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(2):fs%mg%ie_array(2)),&
                                  g_yz(fs%mg%is_array(2):fs%mg%ie_array(2),fs%mg%is_array(3):fs%mg%ie_array(3)),&
                                  g_xz(fs%mg%is_array(1):fs%mg%ie_array(1),fs%mg%is_array(3):fs%mg%ie_array(3))
      real(8),      intent(in) :: g_x( fs%mg%is_array(1):fs%mg%ie_array(1) ),&
                                  g_y( fs%mg%is_array(2):fs%mg%ie_array(2) ),&
                                  g_z( fs%mg%is_array(3):fs%mg%ie_array(3) )
      real(8)                  :: t_sta,t,e_inc_r,e_inc_i
      real(8)                  :: add_inc(3)
      
      !calculate incident pulse
      if(iord==1) then
        t_sta = t1_start
      elseif(iord==2) then
        t_sta = t1_start + t1_t2
      end if
      t = (dble(iter)-0.5d0)*dt_em - t_sta
      call eh_calc_e_inc(amp,t,tw,omega,cep,aes,e_inc_r,e_inc_i)
!      e_inc_r = amp*exp(-0.5d0*(( ((dble(iter)-0.5d0)*dt_em-10.0d0*tw1)/tw1 )**2.0d0) ) !test time factor
      add_inc(:) = e_inc_r*ep_r(:) + e_inc_i*ep_i(:)
      
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
            fe%ex_z(ix,iy,iz)=fe%ex_z(ix,iy,iz)+c2_inc_xyz(3)*add_inc(1)*g_xy(ix,iy)*g_x(ix)*g_y(iy)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy)
          do iy=fe%iey_z_is(2),fe%iey_z_ie(2)
          do ix=fe%iey_z_is(1),fe%iey_z_ie(1)
            fe%ey_z(ix,iy,iz)=fe%ey_z(ix,iy,iz)+c2_inc_xyz(3)*add_inc(2)*g_xy(ix,iy)*g_x(ix)*g_y(iy)
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
            fe%ey_x(ix,iy,iz)=fe%ey_x(ix,iy,iz)+c2_inc_xyz(1)*add_inc(2)*g_yz(iy,iz)*g_y(iy)*g_z(iz)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz)
          do iz=fe%iez_x_is(3),fe%iez_x_ie(3)
          do iy=fe%iez_x_is(2),fe%iez_x_ie(2)
            fe%ez_x(ix,iy,iz)=fe%ez_x(ix,iy,iz)+c2_inc_xyz(1)*add_inc(3)*g_yz(iy,iz)*g_y(iy)*g_z(iz)
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
            fe%ex_y(ix,iy,iz)=fe%ex_y(ix,iy,iz)+c2_inc_xyz(2)*add_inc(1)*g_xz(ix,iz)*g_x(ix)*g_z(iz)
          end do
          end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz)
          do iz=fe%iez_y_is(3),fe%iez_y_ie(3)
          do ix=fe%iez_y_is(1),fe%iez_y_ie(1)
            fe%ez_y(ix,iy,iz)=fe%ez_y(ix,iy,iz)+c2_inc_xyz(2)*add_inc(3)*g_xz(ix,iz)*g_x(ix)*g_z(iz)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
    !+ calculate incident and scattering fields for ase option +++++++++++++++++++++++++++++++
    subroutine eh_calc_inc_sca_ase
      use salmon_global,  only: nt_em,ek_dir1,epdir_re1,epdir_im1,ae_shape1,phi_cep1,E_amplitude1,&
                                ase_num_em,ase_smedia_id_em
      use math_constants, only: zi
      implicit none
      real(8),allocatable :: ex_inc_pa(:),ey_inc_pa(:),ez_inc_pa(:),hx_inc_pa(:),hy_inc_pa(:),hz_inc_pa(:)
      real(8)             :: ex_inc_bt(2),ey_inc_bt(2),ez_inc_bt(2),hx_inc_bt(2),hy_inc_bt(2),hz_inc_bt(2)
      integer             :: ip,ibt,ie
      real(8)             :: t,t_ori,t_inc,wf,e_inc_r,e_inc_i,ex_inc,ey_inc,ez_inc,hx_inc,hy_inc,hz_inc
      complex(8)          :: f_factor
      character(1)        :: p_axis
      
      !set time information and calcularate window function
      t     = dble(iter)*dt_em
      t_ori = t - t1_start
      call eh_calc_wf(t,dble(nt_em)*dt_em,wf)
      
      !speciy propagation axis
      if    (ek_dir1(1)==1.0d0) then
        p_axis = 'x'; ip = 1;
      elseif(ek_dir1(2)==1.0d0) then
        p_axis = 'y'; ip = 2;
      elseif(ek_dir1(3)==1.0d0) then
        p_axis = 'z'; ip = 3;
      end if
      
      !allocation for incident field on propagation axis
      allocate( ex_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)),&
                ey_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)),&
                ez_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)),&
                hx_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)),&
                hy_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)),&
                hz_inc_pa(fs%mg%is(ip):fs%mg%ie(ip)) )
      ex_inc_pa(:)=0.0d0; ey_inc_pa(:)=0.0d0; ez_inc_pa(:)=0.0d0;
      hx_inc_pa(:)=0.0d0; hy_inc_pa(:)=0.0d0; hz_inc_pa(:)=0.0d0;
      
      !calculate incident field on propagation axis
      do ii = fs%mg%is(ip),fs%mg%ie(ip)
        t_inc = t_ori - abs( fe%coo(ii,ip) - fe%coo(fe%inc_po_id(1,ip),ip) )/fe%cspeed_smedia_ase
        call eh_calc_e_inc(E_amplitude1,t_inc,tw1,omega1,phi_cep1,ae_shape1,e_inc_r,e_inc_i)
        call eh_calc_polarization_inc(epdir_re1,epdir_im1,e_inc_r,e_inc_i,   &
                                      ex_inc_pa(ii),ey_inc_pa(ii),ez_inc_pa(ii),&
                                      hx_inc_pa(ii),hy_inc_pa(ii),hz_inc_pa(ii),p_axis)
        hx_inc_pa(ii) = hx_inc_pa(ii)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
        hy_inc_pa(ii) = hy_inc_pa(ii)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
        hz_inc_pa(ii) = hz_inc_pa(ii)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
      end do
      
      !calculate incident field on bottom and top
      do ibt = 1,2
        t_inc = t_ori - abs( fe%coo(fe%id_on_box_surf(ip,ibt),ip) - fe%coo(fe%inc_po_id(1,ip),ip) )/fe%cspeed_smedia_ase
        call eh_calc_e_inc(E_amplitude1,t_inc,tw1,omega1,phi_cep1,ae_shape1,e_inc_r,e_inc_i)
        call eh_calc_polarization_inc(epdir_re1,epdir_im1,e_inc_r,e_inc_i,   &
                                      ex_inc_bt(ibt),ey_inc_bt(ibt),ez_inc_bt(ibt),&
                                      hx_inc_bt(ibt),hy_inc_bt(ibt),hz_inc_bt(ibt),p_axis)
        hx_inc_bt(ibt) = hx_inc_bt(ibt)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
        hy_inc_bt(ibt) = hy_inc_bt(ibt)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
        hz_inc_bt(ibt) = hz_inc_bt(ibt)*sqrt(fe%rep(ase_smedia_id_em)/fe%rmu(ase_smedia_id_em))
      end do
      
      !Fourier transformation for incident field
!$omp parallel
!$omp do private(ibt,ii,ie,f_factor) collapse(1)
      do ie = 1,ase_num_em
        f_factor = dt_em*dble(obs_samp_em)*exp(zi*fe%ase_ene(ie)*t)
        do ii = fs%mg%is(ip),fs%mg%ie(ip)
          fe%ase_ex_inc_pa_ene(ii,ie,1) = fe%ase_ex_inc_pa_ene(ii,ie,1) + f_factor*ex_inc_pa(ii)*wf
          fe%ase_ex_inc_pa_ene(ii,ie,2) = fe%ase_ex_inc_pa_ene(ii,ie,2) + f_factor*ex_inc_pa(ii)
          fe%ase_ey_inc_pa_ene(ii,ie,1) = fe%ase_ey_inc_pa_ene(ii,ie,1) + f_factor*ey_inc_pa(ii)*wf
          fe%ase_ey_inc_pa_ene(ii,ie,2) = fe%ase_ey_inc_pa_ene(ii,ie,2) + f_factor*ey_inc_pa(ii)
          fe%ase_ez_inc_pa_ene(ii,ie,1) = fe%ase_ez_inc_pa_ene(ii,ie,1) + f_factor*ez_inc_pa(ii)*wf
          fe%ase_ez_inc_pa_ene(ii,ie,2) = fe%ase_ez_inc_pa_ene(ii,ie,2) + f_factor*ez_inc_pa(ii)
          fe%ase_hx_inc_pa_ene(ii,ie,1) = fe%ase_hx_inc_pa_ene(ii,ie,1) + f_factor*hx_inc_pa(ii)*wf
          fe%ase_hx_inc_pa_ene(ii,ie,2) = fe%ase_hx_inc_pa_ene(ii,ie,2) + f_factor*hx_inc_pa(ii)
          fe%ase_hy_inc_pa_ene(ii,ie,1) = fe%ase_hy_inc_pa_ene(ii,ie,1) + f_factor*hy_inc_pa(ii)*wf
          fe%ase_hy_inc_pa_ene(ii,ie,2) = fe%ase_hy_inc_pa_ene(ii,ie,2) + f_factor*hy_inc_pa(ii)
          fe%ase_hz_inc_pa_ene(ii,ie,1) = fe%ase_hz_inc_pa_ene(ii,ie,1) + f_factor*hz_inc_pa(ii)*wf
          fe%ase_hz_inc_pa_ene(ii,ie,2) = fe%ase_hz_inc_pa_ene(ii,ie,2) + f_factor*hz_inc_pa(ii)
        end do
        do ibt = 1,2
          fe%ase_ex_inc_bt_ene(ie,ibt,1) = fe%ase_ex_inc_bt_ene(ie,ibt,1) + f_factor*ex_inc_bt(ibt)*wf
          fe%ase_ex_inc_bt_ene(ie,ibt,2) = fe%ase_ex_inc_bt_ene(ie,ibt,2) + f_factor*ex_inc_bt(ibt)
          fe%ase_ey_inc_bt_ene(ie,ibt,1) = fe%ase_ey_inc_bt_ene(ie,ibt,1) + f_factor*ey_inc_bt(ibt)*wf
          fe%ase_ey_inc_bt_ene(ie,ibt,2) = fe%ase_ey_inc_bt_ene(ie,ibt,2) + f_factor*ey_inc_bt(ibt)
          fe%ase_ez_inc_bt_ene(ie,ibt,1) = fe%ase_ez_inc_bt_ene(ie,ibt,1) + f_factor*ez_inc_bt(ibt)*wf
          fe%ase_ez_inc_bt_ene(ie,ibt,2) = fe%ase_ez_inc_bt_ene(ie,ibt,2) + f_factor*ez_inc_bt(ibt)
          fe%ase_hx_inc_bt_ene(ie,ibt,1) = fe%ase_hx_inc_bt_ene(ie,ibt,1) + f_factor*hx_inc_bt(ibt)*wf
          fe%ase_hx_inc_bt_ene(ie,ibt,2) = fe%ase_hx_inc_bt_ene(ie,ibt,2) + f_factor*hx_inc_bt(ibt)
          fe%ase_hy_inc_bt_ene(ie,ibt,1) = fe%ase_hy_inc_bt_ene(ie,ibt,1) + f_factor*hy_inc_bt(ibt)*wf
          fe%ase_hy_inc_bt_ene(ie,ibt,2) = fe%ase_hy_inc_bt_ene(ie,ibt,2) + f_factor*hy_inc_bt(ibt)
          fe%ase_hz_inc_bt_ene(ie,ibt,1) = fe%ase_hz_inc_bt_ene(ie,ibt,1) + f_factor*hz_inc_bt(ibt)*wf
          fe%ase_hz_inc_bt_ene(ie,ibt,2) = fe%ase_hz_inc_bt_ene(ie,ibt,2) + f_factor*hz_inc_bt(ibt)
        end do
      end do
!$omp end do
!$omp end parallel
      
      !xy plane: Fourier transformation for scattering field on closed surface[box shape]
!$omp parallel
!$omp do private(ibt,ix,iy,ie,f_factor,ex_inc,ey_inc,hx_inc,hy_inc) collapse(3)
      do ie  = 1,ase_num_em
      do iy  = fs%mg%is(2),fs%mg%ie(2)
      do ix  = fs%mg%is(1),fs%mg%ie(1)
      do ibt = 1,2
        if((fs%mg%is(3)<=fe%id_on_box_surf(3,ibt)).and.(fe%id_on_box_surf(3,ibt)<=fs%mg%ie(3))) then
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            ex_inc = ex_inc_pa(ix ); ey_inc = ey_inc_pa(ix );
            hx_inc = hx_inc_pa(ix ); hy_inc = hy_inc_pa(ix );
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            ex_inc = ex_inc_pa(iy ); ey_inc = ey_inc_pa(iy );
            hx_inc = hx_inc_pa(iy ); hy_inc = hy_inc_pa(iy );
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            ex_inc = ex_inc_bt(ibt); ey_inc = ey_inc_bt(ibt);
            hx_inc = hx_inc_bt(ibt); hy_inc = hy_inc_bt(ibt);
          end if
          f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%ase_ene(ie)*t)
          fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,1) = fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,1) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%ex_s(ix,iy,fe%id_on_box_surf(3,ibt))-ex_inc)*wf
          fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,2) = fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,2) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%ex_s(ix,iy,fe%id_on_box_surf(3,ibt))-ex_inc)
          fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,1) = fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,1) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%ey_s(ix,iy,fe%id_on_box_surf(3,ibt))-ey_inc)*wf
          fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,2) = fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,2) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%ey_s(ix,iy,fe%id_on_box_surf(3,ibt))-ey_inc)
          fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,1) = fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,1) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%hx_s(ix,iy,fe%id_on_box_surf(3,ibt))-hx_inc)*wf
          fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,2) = fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,2) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%hx_s(ix,iy,fe%id_on_box_surf(3,ibt))-hx_inc)
          fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,1) = fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,1) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%hy_s(ix,iy,fe%id_on_box_surf(3,ibt))-hy_inc)*wf
          fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,2) = fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,2) + dble(fe%iase_mask_xy(ix,iy)) &
                                                *f_factor*(fe%hy_s(ix,iy,fe%id_on_box_surf(3,ibt))-hy_inc)
        end if
      end do
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !yz plane: Fourier transformation for scattering field on closed surface[box shape]
!$omp parallel
!$omp do private(ibt,iy,iz,ie,f_factor,ey_inc,ez_inc,hy_inc,hz_inc) collapse(3)
      do ie  = 1,ase_num_em
      do iz  = fs%mg%is(3),fs%mg%ie(3)
      do iy  = fs%mg%is(2),fs%mg%ie(2)
      do ibt = 1,2
        if((fs%mg%is(1)<=fe%id_on_box_surf(1,ibt)).and.(fe%id_on_box_surf(1,ibt)<=fs%mg%ie(1))) then
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            ey_inc = ey_inc_bt(ibt); ez_inc = ez_inc_bt(ibt);
            hy_inc = hy_inc_bt(ibt); hz_inc = hz_inc_bt(ibt);
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            ey_inc = ey_inc_pa(iy ); ez_inc = ez_inc_pa(iy );
            hy_inc = hy_inc_pa(iy ); hz_inc = hz_inc_pa(iy );
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            ey_inc = ey_inc_pa(iz ); ez_inc = ez_inc_pa(iz );
            hy_inc = hy_inc_pa(iz ); hz_inc = hz_inc_pa(iz );
          end if
          f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%ase_ene(ie)*t)
          fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,1) = fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,1) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%ey_s(fe%id_on_box_surf(1,ibt),iy,iz)-ey_inc)*wf
          fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,2) = fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,2) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%ey_s(fe%id_on_box_surf(1,ibt),iy,iz)-ey_inc)
          fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,1) = fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,1) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%ez_s(fe%id_on_box_surf(1,ibt),iy,iz)-ez_inc)*wf
          fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,2) = fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,2) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%ez_s(fe%id_on_box_surf(1,ibt),iy,iz)-ez_inc)
          fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,1) = fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,1) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%hy_s(fe%id_on_box_surf(1,ibt),iy,iz)-hy_inc)*wf
          fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,2) = fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,2) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%hy_s(fe%id_on_box_surf(1,ibt),iy,iz)-hy_inc)
          fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,1) = fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,1) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%hz_s(fe%id_on_box_surf(1,ibt),iy,iz)-hz_inc)*wf
          fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,2) = fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,2) + dble(fe%iase_mask_yz(iy,iz)) &
                                                *f_factor*(fe%hz_s(fe%id_on_box_surf(1,ibt),iy,iz)-hz_inc)
        end if
      end do
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !xz plane: Fourier transformation for scattering field on closed surface[box shape]
!$omp parallel
!$omp do private(ibt,ix,iz,ie,f_factor,ex_inc,ez_inc,hx_inc,hz_inc) collapse(3)
      do ie  = 1,ase_num_em
      do iz  = fs%mg%is(3),fs%mg%ie(3)
      do ix  = fs%mg%is(1),fs%mg%ie(1)
      do ibt = 1,2
        if((fs%mg%is(2)<=fe%id_on_box_surf(2,ibt)).and.(fe%id_on_box_surf(2,ibt)<=fs%mg%ie(2))) then
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            ex_inc = ex_inc_pa(ix ); ez_inc = ez_inc_pa(ix );
            hx_inc = hx_inc_pa(ix ); hz_inc = hz_inc_pa(ix );
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            ex_inc = ex_inc_bt(ibt); ez_inc = ez_inc_bt(ibt);
            hx_inc = hx_inc_bt(ibt); hz_inc = hz_inc_bt(ibt);
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            ex_inc = ex_inc_pa(iz ); ez_inc = ez_inc_pa(iz );
            hx_inc = hx_inc_pa(iz ); hz_inc = hz_inc_pa(iz );
          end if
          f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%ase_ene(ie)*t)
          fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,1) = fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,1) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%ex_s(ix,fe%id_on_box_surf(2,ibt),iz)-ex_inc)*wf
          fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,2) = fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,2) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%ex_s(ix,fe%id_on_box_surf(2,ibt),iz)-ex_inc)
          fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,1) = fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,1) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%ez_s(ix,fe%id_on_box_surf(2,ibt),iz)-ez_inc)*wf
          fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,2) = fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,2) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%ez_s(ix,fe%id_on_box_surf(2,ibt),iz)-ez_inc)
          fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,1) = fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,1) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%hx_s(ix,fe%id_on_box_surf(2,ibt),iz)-hx_inc)*wf
          fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,2) = fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,2) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%hx_s(ix,fe%id_on_box_surf(2,ibt),iz)-hx_inc)
          fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,1) = fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,1) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%hz_s(ix,fe%id_on_box_surf(2,ibt),iz)-hz_inc)*wf
          fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,2) = fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,2) + dble(fe%iase_mask_xz(ix,iz)) &
                                                *f_factor*(fe%hz_s(ix,fe%id_on_box_surf(2,ibt),iz)-hz_inc)
        end if
      end do
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      
      !deallocation
      deallocate( ex_inc_pa,ey_inc_pa,ez_inc_pa,hx_inc_pa,hy_inc_pa,hz_inc_pa )
      
      return
    end subroutine eh_calc_inc_sca_ase
    
    !+ CONTAINED IN eh_calc ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ calculate incident and total fields for art option ++++++++++++++++++++++++++++++++++++
    subroutine eh_calc_inc_tot_art
      use salmon_global,  only: nt_em,ek_dir1,epdir_re1,epdir_im1,ae_shape1,phi_cep1,E_amplitude1,&
                                art_num_em,art_smedia_id_em
      use math_constants, only: zi
      implicit none
      real(8)             :: ex_inc_bt(2),ey_inc_bt(2),ez_inc_bt(2),hx_inc_bt(2),hy_inc_bt(2),hz_inc_bt(2)
      integer             :: ip,ibt,ie
      real(8)             :: t,t_ori,t_inc,wf,e_inc_r,e_inc_i
      complex(8)          :: f_factor
      character(1)        :: p_axis
      
      !set time information and calcularate window function
      t     = dble(iter)*dt_em
      t_ori = t - t1_start
      call eh_calc_wf(t,dble(nt_em)*dt_em,wf)
      
      !speciy propagation axis
      if    (ek_dir1(1)==1.0d0) then
        p_axis = 'x'; ip = 1;
      elseif(ek_dir1(2)==1.0d0) then
        p_axis = 'y'; ip = 2;
      elseif(ek_dir1(3)==1.0d0) then
        p_axis = 'z'; ip = 3;
      end if
      
      !calculate incident field on plane
      do ibt = 1,2
        t_inc = t_ori - abs( fe%coo(fe%id_on_plane(ibt),ip) - fe%coo(fe%inc_po_id(1,ip),ip) )/fe%cspeed_smedia_art
        call eh_calc_e_inc(E_amplitude1,t_inc,tw1,omega1,phi_cep1,ae_shape1,e_inc_r,e_inc_i)
        call eh_calc_polarization_inc(epdir_re1,epdir_im1,e_inc_r,e_inc_i,   &
                                      ex_inc_bt(ibt),ey_inc_bt(ibt),ez_inc_bt(ibt),&
                                      hx_inc_bt(ibt),hy_inc_bt(ibt),hz_inc_bt(ibt),p_axis)
        hx_inc_bt(ibt) = hx_inc_bt(ibt)*sqrt(fe%rep(art_smedia_id_em)/fe%rmu(art_smedia_id_em))
        hy_inc_bt(ibt) = hy_inc_bt(ibt)*sqrt(fe%rep(art_smedia_id_em)/fe%rmu(art_smedia_id_em))
        hz_inc_bt(ibt) = hz_inc_bt(ibt)*sqrt(fe%rep(art_smedia_id_em)/fe%rmu(art_smedia_id_em))
      end do
      
      !Fourier transformation for incident field on plane
!$omp parallel
!$omp do private(ibt,ie,f_factor) collapse(1)
      do ie = 1,art_num_em
        f_factor = dt_em*dble(obs_samp_em)*exp(zi*fe%art_ene(ie)*t)
        do ibt = 1,2
          fe%art_ex_inc_bt_ene(ie,ibt,1) = fe%art_ex_inc_bt_ene(ie,ibt,1) + f_factor*ex_inc_bt(ibt)*wf
          fe%art_ex_inc_bt_ene(ie,ibt,2) = fe%art_ex_inc_bt_ene(ie,ibt,2) + f_factor*ex_inc_bt(ibt)
          fe%art_ey_inc_bt_ene(ie,ibt,1) = fe%art_ey_inc_bt_ene(ie,ibt,1) + f_factor*ey_inc_bt(ibt)*wf
          fe%art_ey_inc_bt_ene(ie,ibt,2) = fe%art_ey_inc_bt_ene(ie,ibt,2) + f_factor*ey_inc_bt(ibt)
          fe%art_ez_inc_bt_ene(ie,ibt,1) = fe%art_ez_inc_bt_ene(ie,ibt,1) + f_factor*ez_inc_bt(ibt)*wf
          fe%art_ez_inc_bt_ene(ie,ibt,2) = fe%art_ez_inc_bt_ene(ie,ibt,2) + f_factor*ez_inc_bt(ibt)
          fe%art_hx_inc_bt_ene(ie,ibt,1) = fe%art_hx_inc_bt_ene(ie,ibt,1) + f_factor*hx_inc_bt(ibt)*wf
          fe%art_hx_inc_bt_ene(ie,ibt,2) = fe%art_hx_inc_bt_ene(ie,ibt,2) + f_factor*hx_inc_bt(ibt)
          fe%art_hy_inc_bt_ene(ie,ibt,1) = fe%art_hy_inc_bt_ene(ie,ibt,1) + f_factor*hy_inc_bt(ibt)*wf
          fe%art_hy_inc_bt_ene(ie,ibt,2) = fe%art_hy_inc_bt_ene(ie,ibt,2) + f_factor*hy_inc_bt(ibt)
          fe%art_hz_inc_bt_ene(ie,ibt,1) = fe%art_hz_inc_bt_ene(ie,ibt,1) + f_factor*hz_inc_bt(ibt)*wf
          fe%art_hz_inc_bt_ene(ie,ibt,2) = fe%art_hz_inc_bt_ene(ie,ibt,2) + f_factor*hz_inc_bt(ibt)
        end do
      end do
!$omp end do
!$omp end parallel
      
      !Fourier transformation for total field on plane
      if    (ek_dir1(1)==1.0d0) then !x-direction propagation, yz-plane----------------------!
!$omp parallel
!$omp do private(ibt,iy,iz,ie,f_factor) collapse(3)
        do ie  = 1,art_num_em
        do iz  = fs%mg%is(3),fs%mg%ie(3)
        do iy  = fs%mg%is(2),fs%mg%ie(2)
        do ibt = 1,2
          if((fs%mg%is(1)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(1))) then
            f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%art_ene(ie)*t)
            fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,1) = fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,1) &
                                                  +f_factor*fe%ey_s(fe%id_on_plane(ibt),iy,iz)*wf
            fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,2) = fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,2) &
                                                  +f_factor*fe%ey_s(fe%id_on_plane(ibt),iy,iz)
            fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,1) = fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,1) &
                                                  +f_factor*fe%ez_s(fe%id_on_plane(ibt),iy,iz)*wf
            fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,2) = fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,2) &
                                                  +f_factor*fe%ez_s(fe%id_on_plane(ibt),iy,iz)
            fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,1) = fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,1) &
                                                  +f_factor*fe%hy_s(fe%id_on_plane(ibt),iy,iz)*wf
            fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,2) = fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,2) &
                                                  +f_factor*fe%hy_s(fe%id_on_plane(ibt),iy,iz)
            fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,1) = fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,1) &
                                                  +f_factor*fe%hz_s(fe%id_on_plane(ibt),iy,iz)*wf
            fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,2) = fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,2) &
                                                  +f_factor*fe%hz_s(fe%id_on_plane(ibt),iy,iz)
          end if
        end do
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(ek_dir1(2)==1.0d0) then !y-direction propagation, xz-plane----------------------!
!$omp parallel
!$omp do private(ibt,ix,iz,ie,f_factor) collapse(3)
        do ie  = 1,art_num_em
        do iz  = fs%mg%is(3),fs%mg%ie(3)
        do ix  = fs%mg%is(1),fs%mg%ie(1)
        do ibt = 1,2
          if((fs%mg%is(2)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(2))) then
            f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%art_ene(ie)*t)
            fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,1) = fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,1) &
                                                  +f_factor*fe%ex_s(ix,fe%id_on_plane(ibt),iz)*wf
            fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,2) = fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,2) &
                                                  +f_factor*fe%ex_s(ix,fe%id_on_plane(ibt),iz)
            fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,1) = fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,1) &
                                                  +f_factor*fe%ez_s(ix,fe%id_on_plane(ibt),iz)*wf
            fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,2) = fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,2) &
                                                  +f_factor*fe%ez_s(ix,fe%id_on_plane(ibt),iz)
            fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,1) = fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,1) &
                                                  +f_factor*fe%hx_s(ix,fe%id_on_plane(ibt),iz)*wf
            fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,2) = fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,2) &
                                                  +f_factor*fe%hx_s(ix,fe%id_on_plane(ibt),iz)
            fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,1) = fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,1) &
                                                  +f_factor*fe%hz_s(ix,fe%id_on_plane(ibt),iz)*wf
            fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,2) = fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,2) &
                                                  +f_factor*fe%hz_s(ix,fe%id_on_plane(ibt),iz)
          end if
        end do
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      elseif(ek_dir1(3)==1.0d0) then !z-direction propagation, xy-plane----------------------!
!$omp parallel
!$omp do private(ibt,ix,iy,ie,f_factor) collapse(3)
        do ie  = 1,art_num_em
        do iy  = fs%mg%is(2),fs%mg%ie(2)
        do ix  = fs%mg%is(1),fs%mg%ie(1)
        do ibt = 1,2
          if((fs%mg%is(3)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(3))) then
            f_factor                             = dt_em*dble(obs_samp_em)*exp(zi*fe%art_ene(ie)*t)
            fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,1) = fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,1) &
                                                  +f_factor*fe%ex_s(ix,iy,fe%id_on_plane(ibt))*wf
            fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,2) = fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,2) &
                                                  +f_factor*fe%ex_s(ix,iy,fe%id_on_plane(ibt))
            fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,1) = fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,1) &
                                                  +f_factor*fe%ey_s(ix,iy,fe%id_on_plane(ibt))*wf
            fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,2) = fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,2) &
                                                  +f_factor*fe%ey_s(ix,iy,fe%id_on_plane(ibt))
            fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,1) = fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,1) &
                                                  +f_factor*fe%hx_s(ix,iy,fe%id_on_plane(ibt))*wf
            fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,2) = fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,2) &
                                                  +f_factor*fe%hx_s(ix,iy,fe%id_on_plane(ibt))
            fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,1) = fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,1) &
                                                  +f_factor*fe%hy_s(ix,iy,fe%id_on_plane(ibt))*wf
            fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,2) = fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,2) &
                                                  +f_factor*fe%hy_s(ix,iy,fe%id_on_plane(ibt))
          end if
        end do
        end do
        end do
        end do
!$omp end do
!$omp end parallel
      end if
      
      return
    end subroutine eh_calc_inc_tot_art
    
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
  subroutine eh_finalize(fs,fe,ld)
    use salmon_global,   only: dt_em,unit_system,yn_periodic,e_impulse,sysname,                          &
                               nt_em,nenergy,de,base_directory,obs_num_em,obs_samp_em,                   &
                               obs_plane_ene_em,yn_obs_plane_em,ase_num_em,ase_ene_min_em,ase_wav_min_em,&
                               art_num_em,art_ene_min_em,art_wav_min_em,ek_dir1,media_type
    use inputoutput,     only: utime_from_au,ulength_from_au,uenergy_from_au
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication,   only: comm_is_root,comm_summation
    use structures,      only: s_fdtd_system
    use phys_constants,  only: cspeed_au
    use math_constants,  only: pi,zi
    implicit none
    type(s_ld_state), intent(inout) :: ld
    type(s_fdtd_system),intent(in)    :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    complex(8),allocatable :: z_ex1(:,:),z_ey1(:,:),z_ez1(:,:),z_ex2(:,:),z_ey2(:,:),z_ez2(:,:),&
                              z_hx1(:,:),z_hy1(:,:),z_hz1(:,:),z_hx2(:,:),z_hy2(:,:),z_hz2(:,:),&
                              z_jx1(:,:),z_jy1(:,:),z_jz1(:,:),z_jx2(:,:),z_jy2(:,:),z_jz2(:,:)
    integer :: observation_unit,spectrum_unit,spectrum_units(2)
    integer :: info_status,info_close_status
    integer :: plane_map_status,plane_map_close_status
    logical :: plane_map_warned,optional_spectrum_warned,lr_opened,band_opened(2)
    integer :: lr_io_status,band_io_status(2)
    character(256) :: spectrum_names(2)
    integer                :: ii,ij,ik,il,i1,i1s,i2,i2s,ix,iy,iz,ie,ibt
    real(8)                :: rn_vec(2)
    real(8)                :: ss_ds(3,2),se_ds(3,2),ss_sum1,ss_sum2,se_sum1,se_sum2,s_inc
                              !ase: integrated pointing vector for sca.- and ext.-cross-sections
                              !     === 1st index ===   === 2nd index ===
                              !     1:ix on yz plane  | 1:bottom and 2:top
                              !     2:iy on xz plane  |
                              !     3:iz on xy plane  |
    real(8)                :: si_sum1(2),si_sum2(2),& !art: 1:bottom and 2:top
                              sr_sum1,sr_sum2,st_sum1,st_sum2
    complex(8)             :: z_tmp(3)
    complex(8)             :: zex_inc,zey_inc,zez_inc,zhx_inc,zhy_inc,zhz_inc,zex,zey,zez,zhx,zhy,zhz
    character(2)           :: plane_name
    character(256)         :: iobs_name,iene_name,wf_name,save_name
    
    plane_map_warned=.false.
    optional_spectrum_warned=.false.

    !*** output linear response(matter dipole pm and current jm are outputted: pm = -dip and jm = -curr) ******!
    if(fe%flag_lr) then
      if(any(fe%rep/=1d0).or.any(media_type=='ld-diagonal').or.any(media_type=='ld-tensor').or.fe%subpixel_enabled)then
        ! Keep the valid pole-current initial condition and raw trajectory.
        ! The legacy optical normalization omits instantaneous polarization.
        if(comm_is_root(nproc_id_global))then
          write(*,'(A)') 'WARNING: non-unit LD background: legacy impulse spectra omitted; raw pole data retained.'
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  '# Dispersive-pole data only; instantaneous background excluded. Time in a.u.'
          if(yn_periodic=='n')then
            if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  '# 1:Time[a.u.] 2:ddm_e_x[a.u.] 3:ddm_e_y[a.u.] 4:ddm_e_z[a.u.]'
            do ii=1,nt_em
              if(lr_io_status==0)write(spectrum_unit,'(F16.8,3(1X,E23.15E3))',iostat=lr_io_status) &
                  fe%time_lr(ii),fe%dip_lr(ii,1:3)
            enddo
          else
            if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  '# 1:Time[a.u.] 2:Jm_x[a.u.] 3:Jm_y[a.u.] 4:Jm_z[a.u.]'
            do ii=1,nt_em
              if(lr_io_status==0)write(spectrum_unit,'(F16.8,3(1X,E23.15E3))',iostat=lr_io_status) &
                  fe%time_lr(ii),fe%curr_lr(ii,1:3)
            enddo
          endif
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
        endif
      else
      if(yn_periodic=='n') then
        !output time-dependent dipole data
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Real time calculation:"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# ddm_e: Change of dipole moment (electrons/plus definition)"
          select case(unit_system)
          case('au','a.u.')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[a.u.]",                &
                  2, "ddm_e_x[a.u.]",             &
                  3, "ddm_e_y[a.u.]",             &
                  4, "ddm_e_z[a.u.]"
          case('A_eV_fs')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[fs]",                  &
                  2, "ddm_e_x[Angstrom]",         &
                  3, "ddm_e_y[Angstrom]",         &
                  4, "ddm_e_z[Angstrom]"
          end select
          do ii=1,nt_em
            if(lr_io_status==0)write(spectrum_unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=lr_io_status) &
                  fe%time_lr(ii)*utime_from_au,                  &
                  fe%dip_lr(ii,1:3)*ulength_from_au
            if(lr_io_status==0)write(spectrum_unit,*,iostat=lr_io_status)
          end do
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
        end if
        
        !output response data
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,1),fe%fr_lr(:,1),fe%fi_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,2),fe%fr_lr(:,2),fe%fi_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%dip_lr(:,3),fe%fr_lr(:,3),fe%fi_lr(:,3))
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Fourier-transform spectra:"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# alpha: Polarizability"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# df/dE: Strength function"
          select case(unit_system)
          case('au','a.u.')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
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
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
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
            if(lr_io_status==0)write(spectrum_unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=lr_io_status) &
                  dble(ii)*de*uenergy_from_au,                                        &
                  fe%fr_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,1)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fr_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,2)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fr_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  fe%fi_lr(ii,3)/(-e_impulse)*(ulength_from_au**2/fe%uVperm_from_au), &
                  2.0d0*dble(ii)*de/pi*fe%fi_lr(ii,1:3)/(-e_impulse)
            if(lr_io_status==0)write(spectrum_unit,*,iostat=lr_io_status)
          end do
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
        end if
      elseif(yn_periodic=='y') then
        !output time-dependent average matter current density data and average electric field data
        if(comm_is_root(nproc_id_global)) then
          !average matter current density data
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Real time calculation:"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Jm: Matter current density (electrons)"
          select case(unit_system)
          case('au','a.u.')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[a.u.]",                &
                  2, "Jm_x[a.u.]",                &
                  3, "Jm_y[a.u.]",                &
                  4, "Jm_z[a.u.]"
          case('A_eV_fs')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[fs]",                  &
                  2, "Jm_x[1/(fs*Angstrom^2)]",     &
                  3, "Jm_y[1/(fs*Angstrom^2)]",     &
                  4, "Jm_z[1/(fs*Angstrom^2)]"
          end select
          do ii=1,nt_em
            if(lr_io_status==0)write(spectrum_unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=lr_io_status) &
                  fe%time_lr(ii)*utime_from_au,                  &
                 -fe%curr_lr(ii,:)*((ulength_from_au/utime_from_au)/ulength_from_au**3)
            if(lr_io_status==0)write(spectrum_unit,*,iostat=lr_io_status)
          end do
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
          
          !average electric field data
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_rt_e.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Real time calculation:"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# E: Averaged electric field in the unit cell"
          select case(unit_system)
          case('au','a.u.')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[a.u.]",                &
                  2, "E_x[a.u.]",                 &
                  3, "E_y[a.u.]",                 &
                  4, "E_z[a.u.]"
          case('A_eV_fs')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
                  1, "Time[fs]",                  &
                  2, "E_x[V/Angstrom]",           &
                  3, "E_y[V/Angstrom]",           &
                  4, "E_z[V/Angstrom]"
          end select
          do ii=1,nt_em
            if(lr_io_status==0)write(spectrum_unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=lr_io_status) &
                  fe%time_lr(ii)*utime_from_au,                  &
                  fe%e_lr(ii,:)*fe%uVperm_from_au
            if(lr_io_status==0)write(spectrum_unit,*,iostat=lr_io_status)
          end do
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
        end if
        
        !output response data
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,1),fe%fr_lr(:,1),fe%fi_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,2),fe%fr_lr(:,2),fe%fi_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%curr_lr(:,3),fe%fr_lr(:,3),fe%fi_lr(:,3))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%e_lr(:,1)   ,fe%er_lr(:,1),fe%ei_lr(:,1))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%e_lr(:,2)   ,fe%er_lr(:,2),fe%ei_lr(:,2))
        call eh_fourier(nt_em,nenergy,dt_em,de,fe%time_lr,fe%e_lr(:,3)   ,fe%er_lr(:,3),fe%ei_lr(:,3))
        if(comm_is_root(nproc_id_global)) then
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_response.data'
          call spectrum_open(save_name,spectrum_unit,lr_opened,lr_io_status)
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# Fourier-transform spectra:"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# sigma: Conductivity"
          if(lr_io_status==0)write(spectrum_unit,'(A)',iostat=lr_io_status) &
                  "# eps: Dielectric constant"
          select case(unit_system)
          case('au','a.u.')
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
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
            if(lr_io_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=lr_io_status) &
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
            if(lr_io_status==0)write(spectrum_unit,"(F16.8,99(1X,E23.15E3))",advance='no',iostat=lr_io_status) &
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
            if(lr_io_status==0)write(spectrum_unit,*,iostat=lr_io_status)
          end do
          call spectrum_close(save_name,spectrum_unit,lr_opened,lr_io_status)
        end if
      end if
      end if
    end if
    
    !*** output for obs_plane_ene_em option *******************************************************************!
    if(fe%flag_obs) then
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
                         z_hz2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jx1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jx2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jy1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jy2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jz1(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)),&
                         z_jz2(fs%lg%is(i1s):fs%lg%ie(i1s),fs%lg%is(i2s):fs%lg%ie(i2s)) )
                z_ex1(:,:)=(0.0d0,0.0d0); z_ey1(:,:)=(0.0d0,0.0d0); z_ez1(:,:)=(0.0d0,0.0d0);
                z_ex2(:,:)=(0.0d0,0.0d0); z_ey2(:,:)=(0.0d0,0.0d0); z_ez2(:,:)=(0.0d0,0.0d0);
                z_hx1(:,:)=(0.0d0,0.0d0); z_hy1(:,:)=(0.0d0,0.0d0); z_hz1(:,:)=(0.0d0,0.0d0);
                z_hx2(:,:)=(0.0d0,0.0d0); z_hy2(:,:)=(0.0d0,0.0d0); z_hz2(:,:)=(0.0d0,0.0d0);
                z_jx1(:,:)=(0.0d0,0.0d0); z_jy1(:,:)=(0.0d0,0.0d0); z_jz1(:,:)=(0.0d0,0.0d0);
                z_jx2(:,:)=(0.0d0,0.0d0); z_jy2(:,:)=(0.0d0,0.0d0); z_jz2(:,:)=(0.0d0,0.0d0);
                         
                !collect lg-data
                do i2=fs%mg%is(i2s),fs%mg%ie(i2s)
                do i1=fs%mg%is(i1s),fs%mg%ie(i1s)
                  if(ik==1)     then !xy plane
                    z_ex1(i1,i2) = fe%obs_ex_xy_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2) = fe%obs_ey_xy_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2) = fe%obs_ez_xy_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2) = fe%obs_hx_xy_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2) = fe%obs_hy_xy_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2) = fe%obs_hz_xy_ene(i1,i2,ii,ij,il)
                    z_jx1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jx_xy_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ex_xy_ene(i1,i2,ii,ij,il)
                    z_jy1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jy_xy_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ey_xy_ene(i1,i2,ii,ij,il)
                    z_jz1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jz_xy_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ez_xy_ene(i1,i2,ii,ij,il)
                  elseif(ik==2) then !yz plane
                    z_ex1(i1,i2) = fe%obs_ex_yz_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2) = fe%obs_ey_yz_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2) = fe%obs_ez_yz_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2) = fe%obs_hx_yz_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2) = fe%obs_hy_yz_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2) = fe%obs_hz_yz_ene(i1,i2,ii,ij,il)
                    z_jx1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jx_yz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ex_yz_ene(i1,i2,ii,ij,il)
                    z_jy1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jy_yz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ey_yz_ene(i1,i2,ii,ij,il)
                    z_jz1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jz_yz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ez_yz_ene(i1,i2,ii,ij,il)
                  elseif(ik==3) then !xz plane
                    z_ex1(i1,i2) = fe%obs_ex_xz_ene(i1,i2,ii,ij,il)
                    z_ey1(i1,i2) = fe%obs_ey_xz_ene(i1,i2,ii,ij,il)
                    z_ez1(i1,i2) = fe%obs_ez_xz_ene(i1,i2,ii,ij,il)
                    z_hx1(i1,i2) = fe%obs_hx_xz_ene(i1,i2,ii,ij,il)                            
                    z_hy1(i1,i2) = fe%obs_hy_xz_ene(i1,i2,ii,ij,il)
                    z_hz1(i1,i2) = fe%obs_hz_xz_ene(i1,i2,ii,ij,il)
                    z_jx1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jx_xz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ex_xz_ene(i1,i2,ii,ij,il)
                    z_jy1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jy_xz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ey_xz_ene(i1,i2,ii,ij,il)
                    z_jz1(i1,i2) = (                 cspeed_au/(4.0d0*pi)) * fe%obs_jz_xz_ene(i1,i2,ii,ij,il) &
                                  +(zi*obs_plane_ene_em(ii,ij)/(4.0d0*pi)) * fe%obs_ez_xz_ene(i1,i2,ii,ij,il)
                  end if
                end do
                end do
                call comm_summation(z_ex1,z_ex2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_ey1,z_ey2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_ez1,z_ez2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hx1,z_hx2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hy1,z_hy2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_hz1,z_hz2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_jx1,z_jx2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_jy1,z_jy2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                call comm_summation(z_jz1,z_jz2,fs%lg%num(i1s)*fs%lg%num(i2s),nproc_group_global)
                
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
                  open(newunit=spectrum_unit,file=resolve_path(save_name),iostat=plane_map_status)
                  if(plane_map_status==0)then
                  if(il==1)     then
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# Fourier-transformed spatial distribution with window function:"
                  elseif(il==2) then
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# Fourier-transformed spatial distribution without window function:"
                  end if
                  if(ik==1)     then !xy plane
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_1: Grid ID along x axis"
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_2: Grid ID along y axis"
                  elseif(ik==2) then !yz plane
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_1: Grid ID along y axis"
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_2: Grid ID along z axis"
                  elseif(ik==3) then !xz plane
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_1: Grid ID along x axis"
                    if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# ID_2: Grid ID along z axis"
                  end if
                  if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# E: Electric field"
                  if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# H: Magnetic field"
                  if(plane_map_status==0)write(spectrum_unit,'(A)',iostat=plane_map_status) &
                        "# J: Electric current density defined by J = c/(4*pi)*rot(H) + i\omega/(4*pi)*E"
                  select case(unit_system)
                  case('au','a.u.')
                    if(plane_map_status==0)write(spectrum_unit,'(A,E23.15E3,A)',iostat=plane_map_status) &
                        "# Sampling energy: ",obs_plane_ene_em(ii,ij),' a.u.'
                    if(plane_map_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=plane_map_status) &
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
                          14,"Im(H_z)[a.u.]",             &
                          15,"Re(J_x)[a.u.]",             &
                          16,"Im(J_x)[a.u.]",             &
                          17,"Re(J_y)[a.u.]",             &
                          18,"Im(J_y)[a.u.]",             &
                          19,"Re(J_z)[a.u.]",             &
                          20,"Im(J_z)[a.u.]"
                  case('A_eV_fs')
                    if(plane_map_status==0)write(spectrum_unit,'(A,E23.15E3,A)',iostat=plane_map_status) &
                        "# Sampling energy: ",obs_plane_ene_em(ii,ij)*uenergy_from_au,' eV'
                    if(plane_map_status==0)write(spectrum_unit,'("#",99(1X,I0,":",A))',iostat=plane_map_status) &
                          1, "ID_1[none]",                   &
                          2, "ID_2[none]",                   &
                          3, "Re(E_x)[V/Angstrom*fs]",       &
                          4, "Im(E_x)[V/Angstrom*fs]",       &
                          5, "Re(E_y)[V/Angstrom*fs]",       &
                          6, "Im(E_y)[V/Angstrom*fs]",       &
                          7, "Re(E_z)[V/Angstrom*fs]",       &
                          8, "Im(E_z)[V/Angstrom*fs]",       &
                          9, "Re(H_x)[A/Angstrom*fs]",       &
                          10,"Im(H_x)[A/Angstrom*fs]",       &
                          11,"Re(H_y)[A/Angstrom*fs]",       &
                          12,"Im(H_y)[A/Angstrom*fs]",       &
                          13,"Re(H_z)[A/Angstrom*fs]",       &
                          14,"Im(H_z)[A/Angstrom*fs]",       &
                          15,"Re(J_x)[1/(fs*Angstrom^2)*fs]",  &
                          16,"Im(J_x)[1/(fs*Angstrom^2)*fs]",  &
                          17,"Re(J_y)[1/(fs*Angstrom^2)*fs]",  &
                          18,"Im(J_y)[1/(fs*Angstrom^2)*fs]",  &
                          19,"Re(J_z)[1/(fs*Angstrom^2)*fs]",  &
                          20,"Im(J_z)[1/(fs*Angstrom^2)*fs]"
                  end select
                  do i2=fs%lg%is(i2s),fs%lg%ie(i2s)
                  do i1=fs%lg%is(i1s),fs%lg%ie(i1s)
                    if(plane_map_status==0)write(spectrum_unit,'(I8,I8,99(1X,E23.15E3))',iostat=plane_map_status) &
                          i1,i2,                                               &
                          real( z_ex2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          aimag(z_ex2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          real( z_ey2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          aimag(z_ey2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          real( z_ez2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          aimag(z_ez2(i1,i2))*fe%uVperm_from_au*utime_from_au, &
                          real( z_hx2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          aimag(z_hx2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          real( z_hy2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          aimag(z_hy2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          real( z_hz2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          aimag(z_hz2(i1,i2))*fe%uAperm_from_au*utime_from_au, &
                          real( z_jx2(i1,i2))*(1.0d0/(ulength_from_au**2)),    &
                          aimag(z_jx2(i1,i2))*(1.0d0/(ulength_from_au**2)),    &
                          real( z_jy2(i1,i2))*(1.0d0/(ulength_from_au**2)),    &
                          aimag(z_jy2(i1,i2))*(1.0d0/(ulength_from_au**2)),    &
                          real( z_jz2(i1,i2))*(1.0d0/(ulength_from_au**2)),    &
                          aimag(z_jz2(i1,i2))*(1.0d0/(ulength_from_au**2))
                  end do
                  end do
                    close(spectrum_unit,iostat=plane_map_close_status)
                    if(plane_map_status==0)plane_map_status=plane_map_close_status
                  endif
                  if(plane_map_status/=0.and..not.plane_map_warned)then
                    write(*,'(A,I0,2A)')'SALMONPY_FREQUENCY_MAP_WARNING iostat=',plane_map_status, &
                      ' file=',trim(save_name)//'; calculation finished; frequency map may be incomplete'
                    plane_map_warned=.true.
                  endif
                end if
                
                !deallocation
                deallocate(z_ex1,z_ey1,z_ez1,z_hx1,z_hy1,z_hz1,z_jx1,z_jy1,z_jz1,&
                           z_ex2,z_ey2,z_ez2,z_hx2,z_hy2,z_hz2,z_jx2,z_jy2,z_jz2)
              end do !with or without window function-----------------------------------------------------------------!
            end do !plane---------------------------------------------------------------------------------------------!
          end do !energy----------------------------------------------------------------------------------------------!
        end if
      end do !observation point---------------------------------------------------------------------------------------!
    end if
    
    !*** ase option *******************************************************************************************!
    if(fe%flag_ase) then
      !set normal vector for bottom and top
      rn_vec(1) = -1.0d0; rn_vec(2) = 1.0d0;
      
      !calculate and output absorption-, scattering-, and extinction-cross-sections
      do ie = 1,ase_num_em
      do il = 1,2
        !open file and write header
        if((comm_is_root(nproc_id_global)).and.(ie==1)) then
          if    (il==1) then
            wf_name = 'with_wf'
          elseif(il==2) then
            wf_name = 'without_wf'
          end if
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_ase_'//trim(wf_name)//'.data'
          spectrum_names(il)=save_name
          call spectrum_open(spectrum_names(il),spectrum_units(il),band_opened(il),band_io_status(il))
          if    (il==1) then
            if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# Absorption-, Scattering-, and Extinction-cross-sections"     //&
                                   " normalized by the spectral distribution of the incident pulse"//&
                                   " (with window function):"
          elseif(il==2) then
            if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# Absorption-, Scattering-, and Extinction-cross-sections"     //&
                                   " normalized by the spectral distribution of the incident pulse"//&
                                   " (without window function):"
          end if
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# sigma_a: Absorption cross-section"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# sigma_s: Scattering cross-section"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# sigma_e: Extinction cross-section"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# S: Pointing vector of incident pulse along propagation direction"
          select case(unit_system)
          case('au','a.u.')
            if    (ase_ene_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Energy[a.u.]",                 &
                    2, "sigma_a[a.u.]",                &
                    3, "sigma_s[a.u.]",                &
                    4, "sigma_e[a.u.]",                &
                    5, "S[a.u.]"
            elseif(ase_wav_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Wavelength[a.u.]",             &
                    2, "sigma_a[a.u.]",                &
                    3, "sigma_s[a.u.]",                &
                    4, "sigma_e[a.u.]",                &
                    5, "S[a.u.]"
            end if
          case('A_eV_fs')
            if    (ase_ene_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Energy[eV]",                   &
                    2, "sigma_a[Angstrom^2]",          &
                    3, "sigma_s[Angstrom^2]",          &
                    4, "sigma_e[Angstrom^2]",          &
                    5, "S[VA/Angstrom^2*fs^2]"
            elseif(ase_wav_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Wavelength[Angstrom]",         &
                    2, "sigma_a[Angstrom^2]",          &
                    3, "sigma_s[Angstrom^2]",          &
                    4, "sigma_e[Angstrom^2]",          &
                    5, "S[VA/Angstrom^2*fs^2]"
            end if
          end select
        end if
        
        !initialize
        ss_ds(:,:) = 0.0d0; ss_sum1 = 0.0d0; ss_sum2 = 0.0d0;
        se_ds(:,:) = 0.0d0; se_sum1 = 0.0d0; se_sum2 = 0.0d0;
        
        !xy plane: integration(z component of pointing vector)
!$omp parallel
!$omp do private(ibt,ix,iy,zex_inc,zey_inc,zhx_inc,zhy_inc,s_inc) reduction(+:ss_ds,se_ds) collapse(2)
        do iy  = fs%mg%is(2),fs%mg%ie(2)
        do ix  = fs%mg%is(1),fs%mg%ie(1)
        do ibt = 1,2
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_pa_ene( ix,ie,il); zey_inc = fe%ase_ey_inc_pa_ene( ix,ie,il);
            zhx_inc = fe%ase_hx_inc_pa_ene( ix,ie,il); zhy_inc = fe%ase_hy_inc_pa_ene( ix,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ey_inc_pa_ene( ix,ie,il)*conjg(fe%ase_hz_inc_pa_ene( ix,ie,il)) &
                                 -fe%ase_ez_inc_pa_ene( ix,ie,il)*conjg(fe%ase_hy_inc_pa_ene( ix,ie,il)) )
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_pa_ene( iy,ie,il); zey_inc = fe%ase_ey_inc_pa_ene( iy,ie,il);
            zhx_inc = fe%ase_hx_inc_pa_ene( iy,ie,il); zhy_inc = fe%ase_hy_inc_pa_ene( iy,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ez_inc_pa_ene( iy,ie,il)*conjg(fe%ase_hx_inc_pa_ene( iy,ie,il)) &
                                 -fe%ase_ex_inc_pa_ene( iy,ie,il)*conjg(fe%ase_hz_inc_pa_ene( iy,ie,il)) )
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_bt_ene(ie,ibt,il); zey_inc = fe%ase_ey_inc_bt_ene(ie,ibt,il);
            zhx_inc = fe%ase_hx_inc_bt_ene(ie,ibt,il); zhy_inc = fe%ase_hy_inc_bt_ene(ie,ibt,il);
            s_inc   = 0.5d0*real( fe%ase_ex_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hy_inc_bt_ene(ie,ibt,il)) &
                                 -fe%ase_ey_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hx_inc_bt_ene(ie,ibt,il)) )
          end if
          ss_ds(1,ibt) = ss_ds(1,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,il)*conjg(fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,il))   &
                              -fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,il)*conjg(fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,il)) ) &
                        /s_inc
          se_ds(1,ibt) = se_ds(1,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( zex_inc                              *conjg(fe%ase_hy_xy_sca_ene(ix,iy,ie,ibt,il))   &
                              -zey_inc                              *conjg(fe%ase_hx_xy_sca_ene(ix,iy,ie,ibt,il))   &
                              +fe%ase_ex_xy_sca_ene(ix,iy,ie,ibt,il)*conjg(zhy_inc                              )   &
                              -fe%ase_ey_xy_sca_ene(ix,iy,ie,ibt,il)*conjg(zhx_inc                              ) ) &
                        /s_inc
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        ss_ds(1,:) = ss_ds(1,:)*fs%hgs(1)*fs%hgs(2)
        se_ds(1,:) = se_ds(1,:)*fs%hgs(1)*fs%hgs(2)
        
        !yz plane: integration(x component of pointing vector)
!$omp parallel
!$omp do private(ibt,iy,iz,zey_inc,zez_inc,zhy_inc,zhz_inc,s_inc) reduction(+:ss_ds,se_ds) collapse(2)
        do iz  = fs%mg%is(3),fs%mg%ie(3)
        do iy  = fs%mg%is(2),fs%mg%ie(2)
        do ibt = 1,2
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            zey_inc = fe%ase_ey_inc_bt_ene(ie,ibt,il); zez_inc = fe%ase_ez_inc_bt_ene(ie,ibt,il);
            zhy_inc = fe%ase_hy_inc_bt_ene(ie,ibt,il); zhz_inc = fe%ase_hz_inc_bt_ene(ie,ibt,il);
            s_inc   = 0.5d0*real( fe%ase_ey_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hz_inc_bt_ene(ie,ibt,il)) &
                                 -fe%ase_ez_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hy_inc_bt_ene(ie,ibt,il)) )
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            zey_inc = fe%ase_ey_inc_pa_ene( iy,ie,il); zez_inc = fe%ase_ez_inc_pa_ene( iy,ie,il);
            zhy_inc = fe%ase_hy_inc_pa_ene( iy,ie,il); zhz_inc = fe%ase_hz_inc_pa_ene( iy,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ez_inc_pa_ene( iy,ie,il)*conjg(fe%ase_hx_inc_pa_ene( iy,ie,il)) &
                                 -fe%ase_ex_inc_pa_ene( iy,ie,il)*conjg(fe%ase_hz_inc_pa_ene( iy,ie,il)) )
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            zey_inc = fe%ase_ey_inc_pa_ene( iz,ie,il); zez_inc = fe%ase_ez_inc_pa_ene( iz,ie,il);
            zhy_inc = fe%ase_hy_inc_pa_ene( iz,ie,il); zhz_inc = fe%ase_hz_inc_pa_ene( iz,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ex_inc_pa_ene( iz,ie,il)*conjg(fe%ase_hy_inc_pa_ene( iz,ie,il)) &
                                 -fe%ase_ey_inc_pa_ene( iz,ie,il)*conjg(fe%ase_hx_inc_pa_ene( iz,ie,il)) )
          end if
          ss_ds(2,ibt) = ss_ds(2,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,il)*conjg(fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,il))   &
                              -fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,il)*conjg(fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,il)) ) &
                        /s_inc
          se_ds(2,ibt) = se_ds(2,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( zey_inc                              *conjg(fe%ase_hz_yz_sca_ene(iy,iz,ie,ibt,il))   &
                              -zez_inc                              *conjg(fe%ase_hy_yz_sca_ene(iy,iz,ie,ibt,il))   &
                              +fe%ase_ey_yz_sca_ene(iy,iz,ie,ibt,il)*conjg(zhz_inc                              )   &
                              -fe%ase_ez_yz_sca_ene(iy,iz,ie,ibt,il)*conjg(zhy_inc                              ) ) &
                        /s_inc
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        ss_ds(2,:) = ss_ds(2,:)*fs%hgs(2)*fs%hgs(3)
        se_ds(2,:) = se_ds(2,:)*fs%hgs(2)*fs%hgs(3)
        
        !xz plane: integration(y component of pointing vector)
!$omp parallel
!$omp do private(ibt,ix,iz,zex_inc,zez_inc,zhx_inc,zhz_inc,s_inc) reduction(+:ss_ds,se_ds) collapse(2)
        do iz  = fs%mg%is(3),fs%mg%ie(3)
        do ix  = fs%mg%is(1),fs%mg%ie(1)
        do ibt = 1,2
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_pa_ene( ix,ie,il); zez_inc = fe%ase_ez_inc_pa_ene( ix,ie,il);
            zhx_inc = fe%ase_hx_inc_pa_ene( ix,ie,il); zhz_inc = fe%ase_hz_inc_pa_ene( ix,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ey_inc_pa_ene( ix,ie,il)*conjg(fe%ase_hz_inc_pa_ene( ix,ie,il)) &
                                 -fe%ase_ez_inc_pa_ene( ix,ie,il)*conjg(fe%ase_hy_inc_pa_ene( ix,ie,il)) )
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_bt_ene(ie,ibt,il); zez_inc = fe%ase_ez_inc_bt_ene(ie,ibt,il);
            zhx_inc = fe%ase_hx_inc_bt_ene(ie,ibt,il); zhz_inc = fe%ase_hz_inc_bt_ene(ie,ibt,il);
            s_inc   = 0.5d0*real( fe%ase_ez_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hx_inc_bt_ene(ie,ibt,il)) &
                                 -fe%ase_ex_inc_bt_ene(ie,ibt,il)*conjg(fe%ase_hz_inc_bt_ene(ie,ibt,il)) )
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            zex_inc = fe%ase_ex_inc_pa_ene( iz,ie,il); zez_inc = fe%ase_ez_inc_pa_ene( iz,ie,il);
            zhx_inc = fe%ase_hx_inc_pa_ene( iz,ie,il); zhz_inc = fe%ase_hz_inc_pa_ene( iz,ie,il);
            s_inc   = 0.5d0*real( fe%ase_ex_inc_pa_ene( iz,ie,il)*conjg(fe%ase_hy_inc_pa_ene( iz,ie,il)) &
                                 -fe%ase_ey_inc_pa_ene( iz,ie,il)*conjg(fe%ase_hx_inc_pa_ene( iz,ie,il)) )
          end if
          ss_ds(3,ibt) = ss_ds(3,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,il)*conjg(fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,il))   &
                              -fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,il)*conjg(fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,il)) ) &
                        /s_inc
          se_ds(3,ibt) = se_ds(3,ibt) + rn_vec(ibt)*0.5d0                                                           &
                        *real( zez_inc                              *conjg(fe%ase_hx_xz_sca_ene(ix,iz,ie,ibt,il))   &
                              -zex_inc                              *conjg(fe%ase_hz_xz_sca_ene(ix,iz,ie,ibt,il))   &
                              +fe%ase_ez_xz_sca_ene(ix,iz,ie,ibt,il)*conjg(zhx_inc                              )   &
                              -fe%ase_ex_xz_sca_ene(ix,iz,ie,ibt,il)*conjg(zhz_inc                              ) ) &
                        /s_inc
        end do
        end do
        end do
!$omp end do
!$omp end parallel
        ss_ds(3,:) = ss_ds(3,:)*fs%hgs(1)*fs%hgs(3)
        se_ds(3,:) = se_ds(3,:)*fs%hgs(1)*fs%hgs(3)
        
        !output
        ss_sum1=sum(ss_ds(:,:)); se_sum1=-sum(se_ds(:,:));
        call comm_summation(ss_sum1,ss_sum2,nproc_group_global)
        call comm_summation(se_sum1,se_sum2,nproc_group_global)
        if(comm_is_root(nproc_id_global)) then
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%ase_ey_inc_bt_ene(ie,1,il)*conjg(fe%ase_hz_inc_bt_ene(ie,1,il)) &
                               -fe%ase_ez_inc_bt_ene(ie,1,il)*conjg(fe%ase_hy_inc_bt_ene(ie,1,il)) )
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%ase_ez_inc_bt_ene(ie,1,il)*conjg(fe%ase_hx_inc_bt_ene(ie,1,il)) &
                               -fe%ase_ex_inc_bt_ene(ie,1,il)*conjg(fe%ase_hz_inc_bt_ene(ie,1,il)) )
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%ase_ex_inc_bt_ene(ie,1,il)*conjg(fe%ase_hy_inc_bt_ene(ie,1,il)) &
                               -fe%ase_ey_inc_bt_ene(ie,1,il)*conjg(fe%ase_hx_inc_bt_ene(ie,1,il)) )
          end if
          if    (ase_ene_min_em>=0.0d0) then
            if(band_io_status(il)==0)write(spectrum_units(il),"(F16.8,99(1X,E23.15E3))",iostat=band_io_status(il)) &
                  fe%ase_ene(ie)   *uenergy_from_au,          &
                  (se_sum2-ss_sum2)*(ulength_from_au**2.0d0), &
                  ss_sum2          *(ulength_from_au**2.0d0), &
                  se_sum2          *(ulength_from_au**2.0d0), &
                  s_inc            *fe%uVperm_from_au*fe%uAperm_from_au*(utime_from_au**2.0d0)
          elseif(ase_wav_min_em>=0.0d0) then
            if(band_io_status(il)==0)write(spectrum_units(il),"(F16.8,99(1X,E23.15E3))",iostat=band_io_status(il)) &
                  2.0d0*pi*cspeed_au/fe%ase_ene(ie)*ulength_from_au,          &
                  (se_sum2-ss_sum2)                *(ulength_from_au**2.0d0), &
                  ss_sum2                          *(ulength_from_au**2.0d0), &
                  se_sum2                          *(ulength_from_au**2.0d0), &
                  s_inc                            *fe%uVperm_from_au*fe%uAperm_from_au*(utime_from_au**2.0d0)
          end if
        end if
        
        !close file
        if((comm_is_root(nproc_id_global)).and.(ie==ase_num_em)) then
          call spectrum_close(spectrum_names(il),spectrum_units(il),band_opened(il),band_io_status(il))
        end if
      end do
      end do
    end if
    
    !*** art option *******************************************************************************************!
    if(fe%flag_art) then
      !set normal vector for bottom and top
      rn_vec(1) = -1.0d0; rn_vec(2) = 1.0d0;
      
      !calculate and output absorption-, reflection-, and transmission-rates
      do ie = 1,art_num_em
      do il = 1,2
        !open file and write header
        if((comm_is_root(nproc_id_global)).and.(ie==1)) then
          if    (il==1) then
            wf_name = 'with_wf'
          elseif(il==2) then
            wf_name = 'without_wf'
          end if
          save_name=trim(adjustl(base_directory))//'/'//trim(adjustl(sysname))//'_art_'//trim(wf_name)//'.data'
          spectrum_names(il)=save_name
          call spectrum_open(spectrum_names(il),spectrum_units(il),band_opened(il),band_io_status(il))
          if    (il==1) then
            if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# Absorption-, Reflection-, and Transmission-rates"            //&
                                   " normalized by the spectral distribution of the incident pulse"//&
                                   " (with window function):"
          elseif(il==2) then
            if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# Absorption-, Reflection-, and Transmission-rates"            //&
                                   " normalized by the spectral distribution of the incident pulse"//&
                                   " (without window function):"
          end if
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# A: Absorption rate"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# R: Reflection rate"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# T: Transmission rate"
          if(band_io_status(il)==0)write(spectrum_units(il),'(A)',iostat=band_io_status(il)) &
                  "# S: Pointing vector of incident pulse along propagation direction"
          select case(unit_system)
          case('au','a.u.')
            if    (art_ene_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Energy[a.u.]",                 &
                    2, "A %",                          &
                    3, "R %",                          &
                    4, "T %",                          &
                    5, "S[a.u.]"
            elseif(art_wav_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Wavelength[a.u.]",             &
                    2, "A %",                          &
                    3, "R %",                          &
                    4, "T %",                          &
                    5, "S[a.u.]"
            end if
          case('A_eV_fs')
            if    (art_ene_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Energy[eV]",                   &
                    2, "A %",                          &
                    3, "R %",                          &
                    4, "T %",                          &
                    5, "S[VA/Angstrom^2*fs^2]"
            elseif(art_wav_min_em>=0.0d0) then
              if(band_io_status(il)==0)write(spectrum_units(il),'("#",99(1X,I0,":",A))',iostat=band_io_status(il)) &
                    1, "Wavelength[Angstrom]",         &
                    2, "A %",                          &
                    3, "R %",                          &
                    4, "T %",                          &
                    5, "S[VA/Angstrom^2*fs^2]"
            end if
          end select
        end if
        
        !initialize
        si_sum1(:) = 0.0d0; sr_sum1 = 0.0d0; st_sum1 = 0.0d0;
        si_sum2(:) = 0.0d0; sr_sum2 = 0.0d0; st_sum2 = 0.0d0;
        
        !integration
        if    (ek_dir1(1)==1.0d0) then !x-direction propagation, yz-plane----------------------!
!$omp parallel
!$omp do private(ibt,iy,iz,zey_inc,zez_inc,zhy_inc,zhz_inc,zey,zez,zhy,zhz) &
!$omp reduction(+:si_sum1,sr_sum1,st_sum1) collapse(2)
          do iz  = fs%mg%is(3),fs%mg%ie(3)
          do iy  = fs%mg%is(2),fs%mg%ie(2)
          do ibt = 1,2
            if((fs%mg%is(1)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(1))) then
              zey_inc = fe%art_ey_inc_bt_ene(ie,ibt,il); zez_inc = fe%art_ez_inc_bt_ene(ie,ibt,il);
              zhy_inc = fe%art_hy_inc_bt_ene(ie,ibt,il); zhz_inc = fe%art_hz_inc_bt_ene(ie,ibt,il);
              si_sum1(ibt) = si_sum1(ibt) + rn_vec(2)*0.5d0*real( zey_inc*conjg(zhz_inc)-zez_inc*conjg(zhy_inc) )
              if    (ibt==1) then
                zey = fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,il) - zey_inc
                zez = fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,il) - zez_inc
                zhy = fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,il) - zhy_inc
                zhz = fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,il) - zhz_inc
                sr_sum1 = sr_sum1 + rn_vec(ibt)*0.5d0*real( zey*conjg(zhz)-zez*conjg(zhy) )
              elseif(ibt==2) then
                zey = fe%art_ey_tot_bt_ene(iy,iz,ie,ibt,il)
                zez = fe%art_ez_tot_bt_ene(iy,iz,ie,ibt,il)
                zhy = fe%art_hy_tot_bt_ene(iy,iz,ie,ibt,il)
                zhz = fe%art_hz_tot_bt_ene(iy,iz,ie,ibt,il)
                st_sum1 = st_sum1 + rn_vec(ibt)*0.5d0*real( zey*conjg(zhz)-zez*conjg(zhy) )
              end if
            end if
          end do
          end do
          end do
!$omp end do
!$omp end parallel
          si_sum1(:) = si_sum1(:)*fs%hgs(2)*fs%hgs(3)
          sr_sum1    = sr_sum1   *fs%hgs(2)*fs%hgs(3)
          st_sum1    = st_sum1   *fs%hgs(2)*fs%hgs(3)
        elseif(ek_dir1(2)==1.0d0) then !y-direction propagation, xz-plane----------------------!
!$omp parallel
!$omp do private(ibt,ix,iz,zex_inc,zez_inc,zhx_inc,zhz_inc,zex,zez,zhx,zhz) &
!$omp reduction(+:si_sum1,sr_sum1,st_sum1) collapse(2)
          do iz  = fs%mg%is(3),fs%mg%ie(3)
          do ix  = fs%mg%is(1),fs%mg%ie(1)
          do ibt = 1,2
            if((fs%mg%is(2)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(2))) then
              zex_inc = fe%art_ex_inc_bt_ene(ie,ibt,il); zez_inc = fe%art_ez_inc_bt_ene(ie,ibt,il);
              zhx_inc = fe%art_hx_inc_bt_ene(ie,ibt,il); zhz_inc = fe%art_hz_inc_bt_ene(ie,ibt,il);
              si_sum1(ibt) = si_sum1(ibt) + rn_vec(2)*0.5d0*real( zez_inc*conjg(zhx_inc)-zex_inc*conjg(zhz_inc) )
              if    (ibt==1) then
                zex = fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,il) - zex_inc
                zez = fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,il) - zez_inc
                zhx = fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,il) - zhx_inc
                zhz = fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,il) - zhz_inc
                sr_sum1 = sr_sum1 + rn_vec(ibt)*0.5d0*real( zez*conjg(zhx)-zex*conjg(zhz) )
              elseif(ibt==2) then
                zex = fe%art_ex_tot_bt_ene(ix,iz,ie,ibt,il)
                zez = fe%art_ez_tot_bt_ene(ix,iz,ie,ibt,il)
                zhx = fe%art_hx_tot_bt_ene(ix,iz,ie,ibt,il)
                zhz = fe%art_hz_tot_bt_ene(ix,iz,ie,ibt,il)
                st_sum1 = st_sum1 + rn_vec(ibt)*0.5d0*real( zez*conjg(zhx)-zex*conjg(zhz) )
              end if
            end if
          end do
          end do
          end do
!$omp end do
!$omp end parallel
          si_sum1(:) = si_sum1(:)*fs%hgs(1)*fs%hgs(3)
          sr_sum1    = sr_sum1   *fs%hgs(1)*fs%hgs(3)
          st_sum1    = st_sum1   *fs%hgs(1)*fs%hgs(3)
        elseif(ek_dir1(3)==1.0d0) then !z-direction propagation, xy-plane----------------------!
!$omp parallel
!$omp do private(ibt,ix,iy,zex_inc,zey_inc,zhx_inc,zhy_inc,zex,zey,zhx,zhy) &
!$omp reduction(+:si_sum1,sr_sum1,st_sum1) collapse(2)
          do iy  = fs%mg%is(2),fs%mg%ie(2)
          do ix  = fs%mg%is(1),fs%mg%ie(1)
          do ibt = 1,2
            if((fs%mg%is(3)<=fe%id_on_plane(ibt)).and.(fe%id_on_plane(ibt)<=fs%mg%ie(3))) then
              zex_inc = fe%art_ex_inc_bt_ene(ie,ibt,il); zey_inc = fe%art_ey_inc_bt_ene(ie,ibt,il);
              zhx_inc = fe%art_hx_inc_bt_ene(ie,ibt,il); zhy_inc = fe%art_hy_inc_bt_ene(ie,ibt,il);
              si_sum1(ibt) = si_sum1(ibt) + rn_vec(2)*0.5d0*real( zex_inc*conjg(zhy_inc)-zey_inc*conjg(zhx_inc) )
              if    (ibt==1) then
                zex = fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,il) - zex_inc
                zey = fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,il) - zey_inc
                zhx = fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,il) - zhx_inc
                zhy = fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,il) - zhy_inc
                sr_sum1 = sr_sum1 + rn_vec(ibt)*0.5d0*real( zex*conjg(zhy)-zey*conjg(zhx) )
              elseif(ibt==2) then
                zex = fe%art_ex_tot_bt_ene(ix,iy,ie,ibt,il)
                zey = fe%art_ey_tot_bt_ene(ix,iy,ie,ibt,il)
                zhx = fe%art_hx_tot_bt_ene(ix,iy,ie,ibt,il)
                zhy = fe%art_hy_tot_bt_ene(ix,iy,ie,ibt,il)
                st_sum1 = st_sum1 + rn_vec(ibt)*0.5d0*real( zex*conjg(zhy)-zey*conjg(zhx) )
              end if
            end if
          end do
          end do
          end do
!$omp end do
!$omp end parallel
          si_sum1(:) = si_sum1(:)*fs%hgs(1)*fs%hgs(2)
          sr_sum1    = sr_sum1   *fs%hgs(1)*fs%hgs(2)
          st_sum1    = st_sum1   *fs%hgs(1)*fs%hgs(2)
        end if
        
        !output
        call comm_summation(si_sum1(:),si_sum2(:),2,nproc_group_global)
        call comm_summation(sr_sum1   ,sr_sum2   ,  nproc_group_global)
        call comm_summation(st_sum1   ,st_sum2   ,  nproc_group_global)
        sr_sum2 = sr_sum2 / si_sum2(1)
        st_sum2 = st_sum2 / si_sum2(2)
        if(comm_is_root(nproc_id_global)) then
          if    (ek_dir1(1)==1.0d0) then !x-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%art_ey_inc_bt_ene(ie,1,il)*conjg(fe%art_hz_inc_bt_ene(ie,1,il)) &
                               -fe%art_ez_inc_bt_ene(ie,1,il)*conjg(fe%art_hy_inc_bt_ene(ie,1,il)) )
          elseif(ek_dir1(2)==1.0d0) then !y-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%art_ez_inc_bt_ene(ie,1,il)*conjg(fe%art_hx_inc_bt_ene(ie,1,il)) &
                               -fe%art_ex_inc_bt_ene(ie,1,il)*conjg(fe%art_hz_inc_bt_ene(ie,1,il)) )
          elseif(ek_dir1(3)==1.0d0) then !z-direction propagation--------------------------------!
            s_inc = 0.5d0*real( fe%art_ex_inc_bt_ene(ie,1,il)*conjg(fe%art_hy_inc_bt_ene(ie,1,il)) &
                               -fe%art_ey_inc_bt_ene(ie,1,il)*conjg(fe%art_hx_inc_bt_ene(ie,1,il)) )
          end if
          if    (art_ene_min_em>=0.0d0) then
            if(band_io_status(il)==0)write(spectrum_units(il),"(F16.8,99(1X,E23.15E3))",iostat=band_io_status(il)) &
                  fe%art_ene(ie)         *uenergy_from_au,&
                  (1.0d0-sr_sum2-st_sum2)*100.0d0,        &
                  sr_sum2                *100.0d0,        &
                  st_sum2                *100.0d0,        &
                  s_inc                  *fe%uVperm_from_au*fe%uAperm_from_au*(utime_from_au**2.0d0)
          elseif(art_wav_min_em>=0.0d0) then
            if(band_io_status(il)==0)write(spectrum_units(il),"(F16.8,99(1X,E23.15E3))",iostat=band_io_status(il)) &
                  2.0d0*pi*cspeed_au/fe%art_ene(ie)*ulength_from_au,&
                  (1.0d0-sr_sum2-st_sum2)          *100.0d0,        &
                  sr_sum2                          *100.0d0,        &
                  st_sum2                          *100.0d0,        &
                  s_inc                            *fe%uVperm_from_au*fe%uAperm_from_au*(utime_from_au**2.0d0)
          end if
        end if
        
        !close file
        if((comm_is_root(nproc_id_global)).and.(ie==art_num_em)) then
          call spectrum_close(spectrum_names(il),spectrum_units(il),band_opened(il),band_io_status(il))
        end if
      end do
      end do
    end if
    
    !*** make information file used for animation *************************************************************!
    if(fe%flag_obs) then
      if(comm_is_root(nproc_id_global)) then
        open(newunit=observation_unit,file=resolve_path(trim(base_directory)//"/obs0_info.data"),iostat=info_status)
        if(info_status==0)then
          if(info_status==0)write(observation_unit,'(A,A23)',iostat=info_status) &
            'unit_system          =',trim(unit_system)
          if(info_status==0)write(observation_unit,'(A,A23)',iostat=info_status) &
            'yn_periodic          =',yn_periodic
          if(info_status==0)write(observation_unit,'(A,E23.15E3)',iostat=info_status) &
            'dt_em                =',dt_em*utime_from_au
          if(info_status==0)write(observation_unit,'(A,I23)',iostat=info_status) &
            'nt_em                =',nt_em
          if(info_status==0)write(observation_unit,'(3(A,E23.15E3))',iostat=info_status) &
            'al_em                =',fs%rlsize(1)*ulength_from_au,', ',&
                                                                   fs%rlsize(2)*ulength_from_au,', ',&
                                                                   fs%rlsize(3)*ulength_from_au
          if(info_status==0)write(observation_unit,'(3(A,E23.15E3))',iostat=info_status) &
            'dl_em                =',fs%hgs(1)*ulength_from_au,', ',&
                                                                   fs%hgs(2)*ulength_from_au,', ',&
                                                                   fs%hgs(3)*ulength_from_au
          if(info_status==0)write(observation_unit,'(3(A,I23))',iostat=info_status) &
            'lg_sta               =',fs%lg%is(1),', ',fs%lg%is(2),', ',fs%lg%is(3)
          if(info_status==0)write(observation_unit,'(3(A,I23))',iostat=info_status) &
            'lg_end               =',fs%lg%ie(1),', ',fs%lg%ie(2),', ',fs%lg%ie(3)
          if(info_status==0)write(observation_unit,'(A,I23)',iostat=info_status) &
            'obs_num_em           =',obs_num_em
          if(info_status==0)write(observation_unit,'(A,I23)',iostat=info_status) &
            'obs_samp_em          =',obs_samp_em
          do ii=1,obs_num_em
            if(info_status==0)write(observation_unit,'(A,I3,A,A23)',iostat=info_status) &
            'yn_obs_plane_em(',&
                                                          ii,') =',yn_obs_plane_em(ii)
          end do
          if(info_status==0)write(observation_unit,'(A,E23.15E3)',iostat=info_status) &
            'e_max                =',fe%e_max
          if(info_status==0)write(observation_unit,'(A,E23.15E3)',iostat=info_status) &
            'h_max                =',fe%h_max
          close(observation_unit,iostat=info_close_status)
          if(info_status==0)info_status=info_close_status
        endif
        ! Optional animation metadata must not prevent completed-state cleanup.
        if(info_status/=0)write(*,'(A,I0,2A)')'SALMONPY_OBSERVATION_INFO_WARNING iostat=', &
          info_status,' file=',trim(base_directory)//'/obs0_info.data; calculation finished; info may be incomplete'
      end if
    end if
    
    !*** deallocate *******************************************************************************************!
    deallocate(fe%ex_y,fe%c1_ex_y,fe%c2_ex_y,fe%ex_z,fe%c1_ex_z,fe%c2_ex_z,&
               fe%ey_z,fe%c1_ey_z,fe%c2_ey_z,fe%ey_x,fe%c1_ey_x,fe%c2_ey_x,&
               fe%ez_x,fe%c1_ez_x,fe%c2_ez_x,fe%ez_y,fe%c1_ez_y,fe%c2_ez_y,&
               fe%hx_y,fe%c1_hx_y,fe%c2_hx_y,fe%hx_z,fe%c1_hx_z,fe%c2_hx_z,&
               fe%hy_z,fe%c1_hy_z,fe%c2_hy_z,fe%hy_x,fe%c1_hy_x,fe%c2_hy_x,&
               fe%hz_x,fe%c1_hz_x,fe%c2_hz_x,fe%hz_y,fe%c1_hz_y,fe%c2_hz_y)
    if(fe%flag_save) deallocate(fe%ex_s,fe%ey_s,fe%ez_s,fe%hx_s,fe%hy_s,fe%hz_s)
    
    !*** write end ********************************************************************************************!
    if(comm_is_root(nproc_id_global)) then
      write(*,*) "-------------------------------------------------------"
      write(*,*) "**************************"
      write(*,*) "FDTD end"
      write(*,*) "**************************"
    end if
    
    return
    call ld_finalize(ld)
  contains
    subroutine spectrum_open(name,unit,opened,ios)
      character(*),intent(in) :: name
      integer,intent(out) :: unit,ios
      logical,intent(out) :: opened
      unit=0
      open(newunit=unit,file=resolve_path(name),iostat=ios)
      opened=(ios==0)
    end subroutine spectrum_open
    subroutine spectrum_close(name,unit,opened,ios)
      character(*),intent(in) :: name
      integer,intent(in) :: unit
      logical,intent(inout) :: opened
      integer,intent(inout) :: ios
      integer :: close_status
      if(opened)then
        close(unit,iostat=close_status)
        opened=.false.
        if(ios==0)ios=close_status
      endif
      if(ios/=0.and..not.optional_spectrum_warned)then
        write(*,'(A,I0,2A)')'SALMONPY_OPTIONAL_SPECTRUM_WARNING iostat=',ios, &
          ' file=',trim(name)//'; calculation finished; optional spectrum may be incomplete'
        optional_spectrum_warned=.true.
      endif
    end subroutine spectrum_close
  end subroutine eh_finalize
  
  !===========================================================================================
  != prepare mpi, grid, and sendrecv enviroments =============================================
  subroutine eh_mpi_grid_sr(fs,fe,status)
    use salmon_global,     only: nproc_rgrid,nproc_k,nproc_ob
    use parallelization,   only: nproc_id_global,nproc_group_global
    use communication,     only: comm_bcast,comm_is_root
    use set_numcpu,        only: set_numcpu_general,iprefer_domain_distribution,check_numcpu
    use init_communicator, only: init_communicator_dft
    use common_maxwell,    only: stop_em
    use sendrecv_grid,     only: create_sendrecv_neig,init_sendrecv_grid
    use structures,        only: s_fdtd_system, s_parallel_info
    use initialization_sub
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh),   intent(inout) :: fe
    type(s_parallel_info)             :: info
    integer                           :: neig_ng_eh(1:2,1:3)
    integer                           :: ii,grid_status
    integer,optional,intent(out)      :: status
    logical                           :: if_stop
    
    if(present(status))status=0
    !set mpi condition
    if((nproc_k==0).and.(nproc_ob==0).and.(sum(nproc_rgrid(:))==0)) then
      if(present(status))then
        call set_numcpu_general(iprefer_domain_distribution,1,1,nproc_group_global,info,status=grid_status)
        if(grid_status/=0)then
          status=1117;return
        endif
      else
      call set_numcpu_general(iprefer_domain_distribution,1,1,nproc_group_global,info)
      endif
    else
      if((nproc_k>1).or.(nproc_ob>1))then
        if(present(status))then
          status=1105;return
        endif
        call stop_em('For theory = maxwell, nproc_k and nproc_ob must be 1.')
      endif
      info%npk       = nproc_k
      info%nporbital = nproc_ob
      info%nprgrid   = nproc_rgrid
    end if
    if (comm_is_root(nproc_id_global)) then
      if_stop = .not. check_numcpu(nproc_group_global,info)
    end if
    call comm_bcast(if_stop,nproc_group_global)
    if(if_stop)then
      if(present(status))then
        status=1105;return
      endif
      stop 'fail: check_numcpu'
    endif
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      write(*,*) "Parallelization information:"
      write(*,*) "nproc_k          =",info%npk
      write(*,*) "nproc_ob         =",info%nporbital
      write(*,*) "nproc_rgrid(1:3) =",info%nprgrid(1),info%nprgrid(2),info%nprgrid(3)
      write(*,*) "**************************"
    end if
    if(present(status))then
      call init_communicator_dft(nproc_group_global,info,grid_status)
      if(grid_status/=0)then
        status=merge(1107,1105,grid_status==1);return
      endif
    else
      call init_communicator_dft(nproc_group_global,info)
    endif
    
    !initialize r-grid
    call init_grid_whole(fs%rlsize,fs%hgs,fs%lg)
    if(present(status))then
      call init_grid_parallel(info,fs%lg,fs%mg,grid_status)
      if(grid_status/=0)then
        status=1106;return
      endif
    else
      call init_grid_parallel(info,fs%lg,fs%mg) ! legacy standalone behavior
    endif
    !### This process about mg is temporal. #################################!
    !### With modifying init_grid_parallel to be applied to arbitrary Nd, ###!
    !### this process will be removed.#######################################!
    fs%mg%is_overlap(1:3)=fs%mg%is(1:3)-fe%Nd
    fs%mg%ie_overlap(1:3)=fs%mg%ie(1:3)+fe%Nd
    fs%mg%is_array(1:3)  =fs%mg%is(1:3)-fe%Nd
    fs%mg%ie_array(1:3)  =fs%mg%ie(1:3)+fe%Nd
    !########################################################################!
    
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
    call create_sendrecv_neig(neig_ng_eh,info, &
      periodic=fs%a_bc(:,1)=='periodic'.and.fs%a_bc(:,2)=='periodic') ! neighboring node array
    call init_sendrecv_grid(fs%srg_ng,fs%mg,1,info%icomm_r,neig_ng_eh)
    
    return
  end subroutine eh_mpi_grid_sr
  
  !===========================================================================================
  != send and receive eh =====================================================================
  subroutine eh_sendrecv(fs,fe,var)
    use sendrecv_grid,  only: update_overlap_real8,init_sendrecv_grid,dealloc_cache,s_sendrecv_grid
    use structures,     only: s_fdtd_system
    implicit none
    type(s_fdtd_system),intent(inout) :: fs
    type(ls_fdtd_eh), intent(inout)   :: fe
    character(1),intent(in)           :: var
    integer                           :: ix,iy,iz,idir,iaxis,ipass
    type(s_sendrecv_grid)              :: material_halo
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
      ! Material ownership at H sites needs edge/corner ghosts, not just faces.
      ! Use a separate initialization-only cache; field communication is unchanged.
      call init_sendrecv_grid(material_halo,fs%mg,1,fs%srg_ng%icomm,fs%srg_ng%neig)
      do idir=1,3
        do iaxis=1,3
          if(iaxis==idir) cycle
          material_halo%is_block(iaxis,:,:,idir)=fs%mg%is_array(iaxis)
          material_halo%ie_block(iaxis,:,:,idir)=fs%mg%ie_array(iaxis)
        end do
      end do
      ! Each completed exchange advances information across one process face.
      ! Three passes fill faces, edges and corners, including periodic neighbors.
      do ipass=1,3
        call update_overlap_real8(material_halo,fs%mg,fe%rmedia)
      end do
      call dealloc_cache(material_halo)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
!$omp do private(ix,iy,iz) collapse(2)
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
      
      !for obs_plane_ene_em option
      if(fe%flag_obs) then
      if(sum(fe%iobs_num_ene(:))>0) then
        call update_overlap_real8(fs%srg_ng,fs%mg,fe%hx_s)
        call update_overlap_real8(fs%srg_ng,fs%mg,fe%hy_s)
        call update_overlap_real8(fs%srg_ng,fs%mg,fe%hz_s)
      end if
      end if
      
      !deallocate temporary variable
      deallocate(f1,f2,f3)
    end if
    
    return
  end subroutine eh_sendrecv
  
  !===========================================================================================
  != calculate finite difference in eh =======================================================

  
  !===========================================================================================
  != calculate incident pulse ================================================================
  subroutine eh_calc_e_inc(amp,t,tw,omega,cep,aes,e_inc_r,e_inc_i)
    use math_constants, only: pi
    implicit none
    real(8),      intent(in)  :: amp,t,tw,omega,cep
    character(16),intent(in)  :: aes
    real(8),      intent(out) :: e_inc_r,e_inc_i
    real(8) :: theta1,theta2_r,theta2_i,alpha,beta,gamma
    
    theta1   = pi/tw*(t-0.5d0*tw)                         !for cos(theta1)**2
    alpha    = pi/tw                                      !for cos(theta1)**2
    theta2_r = omega*(t-0.5d0*tw) + cep*2d0*pi            !for cos(theta2)
    theta2_i = omega*(t-0.5d0*tw) + cep*2d0*pi+3d0/2d0*pi !for cos(theta2), where this is translated to sin.
    beta     = omega                                      !for cos(theta2)
    if((t>=0.0d0).and.(t<=tw)) then
      gamma = 1.0d0
    else
      gamma = 0.0d0
    end if
    if    (aes=='Ecos2') then
      e_inc_r = amp*cos(theta1)**2*cos(theta2_r)*gamma
      e_inc_i = amp*cos(theta1)**2*cos(theta2_i)*gamma
    elseif(aes=='Acos2') then
      e_inc_r = -amp*(-alpha*sin(2.d0*theta1)*cos(theta2_r) &
                      -beta*cos(theta1)**2*sin(theta2_r))/beta*gamma
      e_inc_i = -amp*(-alpha*sin(2.d0*theta1)*cos(theta2_i) &
                      -beta*cos(theta1)**2*sin(theta2_i))/beta*gamma
    else
      e_inc_r=0.0d0; e_inc_i=0.0d0;
    end if
    
    return
  end subroutine eh_calc_e_inc
  
  !===========================================================================================
  != calculate polarizaiton of incident pulse ================================================
  subroutine eh_calc_polarization_inc(ep_r,ep_i,e_inc_r,e_inc_i,ex,ey,ez,hx,hy,hz,p_axis)
    implicit none
    real(8),      intent(in)    :: ep_r(3),ep_i(3)
    real(8),      intent(in)    :: e_inc_r,e_inc_i
    real(8),      intent(inout) :: ex,ey,ez,hx,hy,hz
    character(1), intent(in)    :: p_axis
    
    !---------!
    ! CAUTION !
    !----------------------------------------------------------------------------------------!
    ! hx-z calculated here assumes an incident pulse propagated from bottom to top direction !
    !----------------------------------------------------------------------------------------!
    
    !set initial value
    ex=0.0d0; ey=0.0d0; ez=0.0d0; hx=0.0d0; hy=0.0d0; hz=0.0d0; 
    
    !x-polarization
    if( (abs(ep_r(1))>=0.0d0) .or. (abs(ep_i(1))>=0.0d0) ) then
        ex = ex + e_inc_r*ep_r(1) + e_inc_i*ep_i(1)
      if     (p_axis=='y') then
        hz = hz - e_inc_r*ep_r(1) - e_inc_i*ep_i(1)
      elseif (p_axis=='z') then
        hy = hy + e_inc_r*ep_r(1) + e_inc_i*ep_i(1)
      end if
    end if
    
    !y-polarization
    if( (abs(ep_r(2))>=0.0d0) .or. (abs(ep_i(2))>=0.0d0) ) then
        ey = ey + e_inc_r*ep_r(2) + e_inc_i*ep_i(2)
      if     (p_axis=='z') then
        hx = hx - e_inc_r*ep_r(2) - e_inc_i*ep_i(2)
      elseif (p_axis=='x') then
        hz = hz + e_inc_r*ep_r(2) + e_inc_i*ep_i(2)
      end if
    end if
    
    !z-polarization
    if( (abs(ep_r(3))>=0.0d0) .or. (abs(ep_i(3))>=0.0d0) ) then
        ez = ez + e_inc_r*ep_r(3) + e_inc_i*ep_i(3)
      if     (p_axis=='x') then
        hy = hy - e_inc_r*ep_r(3) - e_inc_i*ep_i(3)
      elseif (p_axis=='y') then
        hx = hx + e_inc_r*ep_r(3) + e_inc_i*ep_i(3)
      end if
    end if
    
    return
  end subroutine eh_calc_polarization_inc
  
  !===========================================================================================
  != calculate window funciton ===============================================================
  subroutine eh_calc_wf(t,t_max,wf)
    implicit none
    real(8),      intent(in)  :: t,t_max
    real(8),      intent(out) :: wf
    
    wf = 1.0d0 -3.0d0*(t/t_max)**2.0d0 +2.0d0*(t/t_max)**3.0d0
    
    return
  end subroutine eh_calc_wf
  
  !===========================================================================================
  != save plane data =========================================================================
  subroutine eh_save_plane(id,ipl,conv,ng_is,ng_ie,lg_is,lg_ie,Nd,ifn,iobs,iter,f,var,io_warned)
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
    integer :: output_unit,ios,close_status
    logical,optional,intent(inout) :: io_warned
    logical :: report
    character(2),intent(in) :: var
    real(8),allocatable     :: save_pl(:,:),save_pl2(:,:)
    integer                 :: ii,inum,i1,i1s,i2,i2s
    character(2)            :: plane_name
    character(256)          :: iobs_name,iter_name,save_name
    
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
        open(newunit=output_unit,file=resolve_path(save_name),iostat=ios)
        if(ios==0)then
          do i2=lg_is(i2s),lg_ie(i2s)
          do i1=lg_is(i1s),lg_ie(i1s)
            if(ios==0)write(output_unit,'(I8,I8,1X,E23.15E3)',iostat=ios)i1,i2,save_pl2(i1,i2)*conv
          end do
          end do
          close(output_unit,iostat=close_status)
          if(ios==0)ios=close_status
        endif
        if(ios/=0)then
          report=.true.
          if(present(io_warned))report=.not.io_warned
          if(report)write(*,'(A,I0,2A)')'SALMONPY_PLANE_SNAPSHOT_WARNING iostat=',ios, &
            ' file=',trim(save_name)//'; calculation continues; snapshot may be incomplete'
          if(present(io_warned))io_warned=.true.
        endif
      end if
      
      !deallocate
      deallocate(save_pl,save_pl2)
    end do
    
    return
  end subroutine eh_save_plane
  
  !===========================================================================================
  != save plane data =========================================================================
  subroutine eh_save_plane_integral(id,ipl,conv_e,conv_h,ti,ng_is,ng_ie,dl,Nd,ifn,iobs,ex,ey,ez,hx,hy,hz,io_warned)
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
    logical,optional,intent(inout) :: io_warned
    character(2)     :: plane_name
    character(256)   :: save_name
    
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
        call eh_integral_output(save_name,.false.,ti, &
              [ex_sum2*conv_e,ey_sum2*conv_e,ez_sum2*conv_e,hx_sum2*conv_h,hy_sum2*conv_h,hz_sum2*conv_h, &
               px_sum2*(uenergy_from_au/utime_from_au), &
               py_sum2*(uenergy_from_au/utime_from_au), &
               pz_sum2*(uenergy_from_au/utime_from_au)],io_warned)
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
    integer    :: ii,ij,i1,i2,ix,iy,iz
    real(8)    :: t,wf,rtmp_x,rtmp_y,rtmp_z
    complex(8) :: f_factor
    
    !update time-information and window function
    t = dble(iter) * dt_em
    call eh_calc_wf(t,dble(nt_em)*dt_em,wf)
    
    !Fourier transformation
    do ii=1,3
      if(fe%iobs_pl_pe(iobs,ii)==1) then
        do ij=1,fe%iobs_num_ene(iobs)
          !update f_factor
          f_factor = dt_em*dble(obs_samp_em)*exp(zi*obs_plane_ene_em(iobs,ij)*t)
          
          !update time-integration
          if    (ii==1) then !xy plane
!$omp parallel
!$omp do private(i1,i2,ix,iy,iz,rtmp_x,rtmp_y,rtmp_z)
            do i2=fs%mg%is(2),fs%mg%ie(2)
            do i1=fs%mg%is(1),fs%mg%ie(1)
              ix = i1; iy = i2; iz = fe%iobs_po_id(iobs,3);
              rtmp_x = ( fe%hz_s(ix  ,iy+1,iz  ) - fe%hz_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0) &
                      -( fe%hy_s(ix  ,iy  ,iz+1) - fe%hy_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0)
              rtmp_y = ( fe%hx_s(ix  ,iy  ,iz+1) - fe%hx_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0) &
                      -( fe%hz_s(ix+1,iy  ,iz  ) - fe%hz_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0)
              rtmp_z = ( fe%hy_s(ix+1,iy  ,iz  ) - fe%hy_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0) &
                      -( fe%hx_s(ix  ,iy+1,iz  ) - fe%hx_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0)
              fe%obs_ex_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ex_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%ex_s(ix,iy,iz)*wf
              fe%obs_ex_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ex_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%ex_s(ix,iy,iz)
              fe%obs_ey_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ey_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%ey_s(ix,iy,iz)*wf
              fe%obs_ey_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ey_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%ey_s(ix,iy,iz)
              fe%obs_ez_xy_ene(i1,i2,iobs,ij,1) = fe%obs_ez_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%ez_s(ix,iy,iz)*wf
              fe%obs_ez_xy_ene(i1,i2,iobs,ij,2) = fe%obs_ez_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%ez_s(ix,iy,iz)
              fe%obs_hx_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hx_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%hx_s(ix,iy,iz)*wf
              fe%obs_hx_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hx_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%hx_s(ix,iy,iz)
              fe%obs_hy_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hy_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%hy_s(ix,iy,iz)*wf
              fe%obs_hy_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hy_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%hy_s(ix,iy,iz)
              fe%obs_hz_xy_ene(i1,i2,iobs,ij,1) = fe%obs_hz_xy_ene(i1,i2,iobs,ij,1) + f_factor*fe%hz_s(ix,iy,iz)*wf
              fe%obs_hz_xy_ene(i1,i2,iobs,ij,2) = fe%obs_hz_xy_ene(i1,i2,iobs,ij,2) + f_factor*fe%hz_s(ix,iy,iz)
              fe%obs_jx_xy_ene(i1,i2,iobs,ij,1) = fe%obs_jx_xy_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_x*wf
              fe%obs_jx_xy_ene(i1,i2,iobs,ij,2) = fe%obs_jx_xy_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_x
              fe%obs_jy_xy_ene(i1,i2,iobs,ij,1) = fe%obs_jy_xy_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_y*wf
              fe%obs_jy_xy_ene(i1,i2,iobs,ij,2) = fe%obs_jy_xy_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_y
              fe%obs_jz_xy_ene(i1,i2,iobs,ij,1) = fe%obs_jz_xy_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_z*wf
              fe%obs_jz_xy_ene(i1,i2,iobs,ij,2) = fe%obs_jz_xy_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_z
            end do
            end do
!$omp end do
!$omp end parallel
          elseif(ii==2) then !yz plane
!$omp parallel
!$omp do private(i1,i2,ix,iy,iz,rtmp_x,rtmp_y,rtmp_z)
            do i2=fs%mg%is(3),fs%mg%ie(3)
            do i1=fs%mg%is(2),fs%mg%ie(2)
              ix = fe%iobs_po_id(iobs,1); iy = i1; iz = i2;
              rtmp_x = ( fe%hz_s(ix  ,iy+1,iz  ) - fe%hz_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0) &
                      -( fe%hy_s(ix  ,iy  ,iz+1) - fe%hy_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0)
              rtmp_y = ( fe%hx_s(ix  ,iy  ,iz+1) - fe%hx_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0) &
                      -( fe%hz_s(ix+1,iy  ,iz  ) - fe%hz_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0)
              rtmp_z = ( fe%hy_s(ix+1,iy  ,iz  ) - fe%hy_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0) &
                      -( fe%hx_s(ix  ,iy+1,iz  ) - fe%hx_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0)
              fe%obs_ex_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ex_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ex_s(ix,iy,iz)*wf
              fe%obs_ex_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ex_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ex_s(ix,iy,iz)
              fe%obs_ey_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ey_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ey_s(ix,iy,iz)*wf
              fe%obs_ey_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ey_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ey_s(ix,iy,iz)
              fe%obs_ez_yz_ene(i1,i2,iobs,ij,1) = fe%obs_ez_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ez_s(ix,iy,iz)*wf
              fe%obs_ez_yz_ene(i1,i2,iobs,ij,2) = fe%obs_ez_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ez_s(ix,iy,iz)
              fe%obs_hx_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hx_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hx_s(ix,iy,iz)*wf
              fe%obs_hx_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hx_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hx_s(ix,iy,iz)
              fe%obs_hy_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hy_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hy_s(ix,iy,iz)*wf
              fe%obs_hy_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hy_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hy_s(ix,iy,iz)
              fe%obs_hz_yz_ene(i1,i2,iobs,ij,1) = fe%obs_hz_yz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hz_s(ix,iy,iz)*wf
              fe%obs_hz_yz_ene(i1,i2,iobs,ij,2) = fe%obs_hz_yz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hz_s(ix,iy,iz)
              fe%obs_jx_yz_ene(i1,i2,iobs,ij,1) = fe%obs_jx_yz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_x*wf
              fe%obs_jx_yz_ene(i1,i2,iobs,ij,2) = fe%obs_jx_yz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_x
              fe%obs_jy_yz_ene(i1,i2,iobs,ij,1) = fe%obs_jy_yz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_y*wf
              fe%obs_jy_yz_ene(i1,i2,iobs,ij,2) = fe%obs_jy_yz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_y
              fe%obs_jz_yz_ene(i1,i2,iobs,ij,1) = fe%obs_jz_yz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_z*wf
              fe%obs_jz_yz_ene(i1,i2,iobs,ij,2) = fe%obs_jz_yz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_z
            end do
            end do
!$omp end do
!$omp end parallel
          elseif(ii==3) then !xz plane
!$omp parallel
!$omp do private(i1,i2,ix,iy,iz,rtmp_x,rtmp_y,rtmp_z)
            do i2=fs%mg%is(3),fs%mg%ie(3)
            do i1=fs%mg%is(1),fs%mg%ie(1)
              ix = i1; iy = fe%iobs_po_id(iobs,2); iz = i2;
              rtmp_x = ( fe%hz_s(ix  ,iy+1,iz  ) - fe%hz_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0) &
                      -( fe%hy_s(ix  ,iy  ,iz+1) - fe%hy_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0)
              rtmp_y = ( fe%hx_s(ix  ,iy  ,iz+1) - fe%hx_s(ix  ,iy  ,iz-1) ) / (fs%hgs(3)*2.0d0) &
                      -( fe%hz_s(ix+1,iy  ,iz  ) - fe%hz_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0)
              rtmp_z = ( fe%hy_s(ix+1,iy  ,iz  ) - fe%hy_s(ix-1,iy  ,iz  ) ) / (fs%hgs(1)*2.0d0) &
                      -( fe%hx_s(ix  ,iy+1,iz  ) - fe%hx_s(ix  ,iy-1,iz  ) ) / (fs%hgs(2)*2.0d0)
              fe%obs_ex_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ex_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ex_s(ix,iy,iz)*wf
              fe%obs_ex_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ex_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ex_s(ix,iy,iz)
              fe%obs_ey_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ey_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ey_s(ix,iy,iz)*wf
              fe%obs_ey_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ey_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ey_s(ix,iy,iz)
              fe%obs_ez_xz_ene(i1,i2,iobs,ij,1) = fe%obs_ez_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%ez_s(ix,iy,iz)*wf
              fe%obs_ez_xz_ene(i1,i2,iobs,ij,2) = fe%obs_ez_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%ez_s(ix,iy,iz)
              fe%obs_hx_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hx_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hx_s(ix,iy,iz)*wf
              fe%obs_hx_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hx_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hx_s(ix,iy,iz)
              fe%obs_hy_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hy_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hy_s(ix,iy,iz)*wf
              fe%obs_hy_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hy_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hy_s(ix,iy,iz)
              fe%obs_hz_xz_ene(i1,i2,iobs,ij,1) = fe%obs_hz_xz_ene(i1,i2,iobs,ij,1) + f_factor*fe%hz_s(ix,iy,iz)*wf
              fe%obs_hz_xz_ene(i1,i2,iobs,ij,2) = fe%obs_hz_xz_ene(i1,i2,iobs,ij,2) + f_factor*fe%hz_s(ix,iy,iz)
              fe%obs_jx_xz_ene(i1,i2,iobs,ij,1) = fe%obs_jx_xz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_x*wf
              fe%obs_jx_xz_ene(i1,i2,iobs,ij,2) = fe%obs_jx_xz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_x
              fe%obs_jy_xz_ene(i1,i2,iobs,ij,1) = fe%obs_jy_xz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_y*wf
              fe%obs_jy_xz_ene(i1,i2,iobs,ij,2) = fe%obs_jy_xz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_y
              fe%obs_jz_xz_ene(i1,i2,iobs,ij,1) = fe%obs_jz_xz_ene(i1,i2,iobs,ij,1) + f_factor*rtmp_z*wf
              fe%obs_jz_xz_ene(i1,i2,iobs,ij,2) = fe%obs_jz_xz_ene(i1,i2,iobs,ij,2) + f_factor*rtmp_z
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
    real(8)              :: hw,wf
    complex(8)           :: zf
    
    !apply window function
    if(yn_wf_em=='y') then
      do it=1,nt
        call eh_calc_wf(ti(it),maxval(ti(:)),wf)
        ft_wf(it)=ft(it)*wf
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
  
end module fdtd_eh
