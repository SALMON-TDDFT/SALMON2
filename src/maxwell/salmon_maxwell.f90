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
!-----------------------------------------------------------------------------------------
module salmon_maxwell
  implicit none
  
  type s_fdtd_work
    !share
    
    !weyl
    
    !coulomb
    
    !eh
    real(8)             :: c_0             !light speed
    integer             :: Nd              !number of additional grid in mpi
    integer             :: iter_sta        !start of time-iteration
    integer             :: iter_end        !end of time-iteration
    integer             :: ifn             !file number for inputing or saving data
    integer             :: ioddeven(3)     !odd or even grid paterns
    integer             :: ipml_l          !pml parameter
    real(8)             :: pml_m           !pml parameter
    real(8)             :: pml_r           !pml parameter
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
    integer             :: inc_num         !number of incident current source
    integer,allocatable :: inc_po_pe(:)    !processor element at incident current source point
    integer,allocatable :: inc_li_pe(:,:)  !processor element at incident current source line
    integer,allocatable :: inc_pl_pe(:,:)  !processor element at incident current source plane
    integer,allocatable :: inc_po_id(:,:)  !id at incident current source point
    character(16)       :: inc_dist1       !spatial distribution type of inc.
    character(16)       :: inc_dist2       !spatial distribution type of inc.
    real(8)             :: c2_inc_xyz(3)   !coeff. for inc., xyz(1:3) means propa. direc. of the inc.
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
    real(8),allocatable :: ex_s(:,:,:),ey_s(:,:,:),ez_s(:,:,:)                                        !e for save
    real(8),allocatable :: hx_s(:,:,:),hy_s(:,:,:),hz_s(:,:,:)                                        !h for save
    real(8),allocatable :: c2_jx(:,:,:),c2_jy(:,:,:),c2_jz(:,:,:)       !coeff. for general curr. dens.
    integer             :: inum_d                                       !Drude: number of media
    integer,allocatable :: imedia_d(:)                                  !Drude: imedia number for drude model
    integer,allocatable :: idx_d(:,:,:,:),idy_d(:,:,:,:),idz_d(:,:,:,:) !Drude: id for each component
    real(8),allocatable :: rjx_d(:,:,:,:),rjy_d(:,:,:,:),rjz_d(:,:,:,:) !Drude: poparization current density
    real(8),allocatable :: rjx_sum_d(:,:,:),rjy_sum_d(:,:,:),&          !Drude: sum of poparization current density
                           rjz_sum_d(:,:,:)                             
    real(8),allocatable :: wex_d(:,:,:,:),wey_d(:,:,:,:),wez_d(:,:,:,:) !Drude: weight function for smoothing
    real(8),allocatable :: c1_j_d(:),c2_j_d(:)                          !Drude: coefficient for j
    real(8),allocatable :: rmedia(:,:,:)                                !Material information for tmp.
    real(8),allocatable :: time_lr(:)                                   !LR: time
    integer             :: iter_lr                                      !LR: time iteration for save
    real(8),allocatable :: fr_lr(:,:)                                   !LR: Re[f]
    real(8),allocatable :: fi_lr(:,:)                                   !LR: Im[f]
    real(8),allocatable :: rjx_lr(:,:,:),rjy_lr(:,:,:),rjz_lr(:,:,:)    !LR: poparization current density
    real(8),allocatable :: px_lr(:,:,:), py_lr(:,:,:), pz_lr(:,:,:)     !LR: poparization vector
    real(8),allocatable :: curr_lr(:,:)                                 !LR: current
    real(8),allocatable :: dip_lr(:,:)                                  !LR: dipolemoment
  end type s_fdtd_work
  
  contains
  
  subroutine init_maxwell(fs,ff,fw)
    use inputoutput, only: theory, use_ms_maxwell
    use structures,  only: s_fdtd_system, s_fdtd_field
    implicit none
    type(s_fdtd_system) :: fs
    type(s_fdtd_field)  :: ff
    type(s_fdtd_work)   :: fw
    
    select case(theory)
    case('Maxwell+TDDFT')
      !this selection is temporary. 
      !After removing use_ms_maxwell, this selection is revised.
      select case(use_ms_maxwell) 
      case('y')
        fs%gauge = 'weyl' 
      case('s')
        fs%gauge = 'coulomb' 
      end select
    case('Maxwell')
      fs%gauge = 'eh' 
    case default
      stop 'invalid theory'
    end select
    
    select case(fs%gauge)
    case('weyl')
      call weyl_init(fs,ff,fw)
    case('coulomb')
      call coulomb_init(fs,ff,fw)
    case('eh')
      call eh_init(fs,fw)
    end select
  end subroutine init_maxwell
  
  subroutine finalize_maxwell(fs,ff,fw)
    use structures,  only: s_fdtd_system, s_fdtd_field
    implicit none
    type(s_fdtd_system) :: fs
    type(s_fdtd_field)  :: ff
    type(s_fdtd_work)   :: fw
    
    select case(fs%gauge)
    case('weyl')
      call weyl_finalize(fs,ff,fw)
    case('coulomb')
      call coulomb_finalize(fs,ff,fw)
    case('eh')
      call eh_finalize(fs,fw)
    end select
  end subroutine finalize_maxwell
  
  subroutine calc_maxwell(fs,ff,fw)
    use structures,  only: s_fdtd_system, s_fdtd_field
    implicit none
    type(s_fdtd_system) :: fs
    type(s_fdtd_field)  :: ff
    type(s_fdtd_work)   :: fw
    
    select case(fs%gauge)
    case('weyl')
      call weyl_calc(fs,ff,fw)
    case('coulomb')
      call coulomb_calc(fs,ff,fw)
    case('eh')
      call eh_calc(fs,fw)
    end select
  end subroutine calc_maxwell
  
end module salmon_maxwell
