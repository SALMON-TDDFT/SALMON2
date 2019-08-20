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
subroutine eh_init(fs,fw)
  use inputoutput,          only: nt_em,al_em,dl_em,dt_em,boundary_em,&
                                  utime_from_au,ulength_from_au,uenergy_from_au,unit_system,&
                                  uenergy_to_au,ulength_to_au,ucharge_to_au,iperiodic,directory,&
                                  imedia_num,shape_file,epsilon,rmu,sigma,type_media,&
                                  pole_num_ld,omega_p_ld,f_ld,gamma_ld,omega_ld,&
                                  iobs_num_em,obs_loc_em,wave_input,trans_longi,e_impulse,nenergy,&
                                  source_loc1,ek_dir1,epdir_re1,epdir_im1,ae_shape1,&
                                  phi_cep1,rlaser_int_wcm2_1,amplitude1,&
                                  source_loc2,ek_dir2,epdir_re2,epdir_im2,ae_shape2,&
                                  phi_cep2,rlaser_int_wcm2_2,amplitude2
  use salmon_parallel,      only: nproc_id_global, nproc_group_global
  use salmon_communication, only: comm_is_root, comm_bcast
  use structures,           only: s_fdtd_system
  use salmon_maxwell,       only: ls_fdtd_work
  use math_constants,       only: pi
  implicit none
  type(s_fdtd_system),intent(inout) :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  integer                           :: ii,ij,ix,iy,iz,icount,icount_ld,iflag_lr,iflag_pml
  real(8)                           :: dt_cfl,diff_cep
  character(1)                      :: dir
  character(128)                    :: save_name
  
  !set initial parameter and value
  fw%c_0        = 1.370359991378353d2
  fw%Nd         = 1
  fw%iter_sta   = 1
  fw%iter_end   = nt_em
  fs%rlsize(:)  = al_em(:)
  fs%hgs(:)     = dl_em(:)
  fw%ifn        = 600
  fw%ipml_l     = 8
  fw%pml_m      = 4.0d0
  fw%pml_r      = 1.0d-7
  do ii=1,3
  do ij=1,2
    select case(boundary_em(ii,ij))
    case('default')
      if(iperiodic==0) then
        fs%a_bc(ii,ij)='pml'
        iflag_pml   =1
      elseif(iperiodic==3) then
        fs%a_bc(ii,ij)='periodic'
      end if
    case('pml')
      fs%a_bc(ii,ij)='pml'
      iflag_pml   =1
    case('pec')
      if(comm_is_root(nproc_id_global).and.(iperiodic==3)) &
      write(*,*) "For iperiodic = 3, boundary_em must be default or pml."
      stop
      fs%a_bc(ii,ij)='pec'
    case('periodic')
      if(comm_is_root(nproc_id_global).and.(iperiodic==0)) &
      write(*,*) "For iperiodic = 0, boundary_em must be default, pml, or pec."
      stop
    end select
  end do
  end do
  select case(unit_system)
  case('au','a.u.')
    fw%uVperm_from_au=1.0d0
    fw%uAperm_from_au=1.0d0
  case('A_eV_fs')
    !see amplitude1 or amplitude2 in src/io/iunputoutput.f90
    fw%uVperm_from_au=1/(uenergy_to_au/ulength_to_au/ucharge_to_au)
    fw%uAperm_from_au=fw%uVperm_from_au
  end select
  
  !prepare GCEED(set mpi condition, gird, and sendrecv environment)
  call eh_prep_GCEED(fs,fw)
  
  !set coordinate
  allocate(fw%coo(minval(fs%lg%is(:))-fw%Nd:maxval(fs%lg%ie(:))+fw%Nd,3))
  call set_coo_em(iperiodic,fw%Nd,fw%ioddeven(:),fs%lg%is(:),fs%lg%ie(:),fs%hgs(:),fw%coo(:,:))
  
  !set and check dt
  dt_cfl=1.0d0/( &
         fw%c_0*sqrt( (1.0d0/fs%hgs(1))**2.0d0+(1.0d0/fs%hgs(2))**2.0d0+(1.0d0/fs%hgs(3))**2.0d0 ) &
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
  
  !basic allocation in eh-FDTD
  call eh_allocate
  
  !input fdtd shape
  allocate(fs%imedia(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
  fs%imedia(:,:,:)=0
  allocate(fw%rmedia(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
  fw%rmedia(:,:,:)=0.0d0
  if(imedia_num>0) then
    !check file format and input shape file
    if(comm_is_root(nproc_id_global)) write(*,*)
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
    if(index(shape_file,".cube", back=.true.)/=0) then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "shape file is inputed by .cube format."
      end if
      call eh_input_shape(fw%ifn,fs%ng%is,fs%ng%ie,fs%lg%is,fs%lg%ie,fw%Nd,fs%imedia,'cu')
      fw%rmedia(:,:,:)=dble(fs%imedia(:,:,:))
      call eh_sendrecv(fs,fw,'r')
      fs%imedia(:,:,:)=int(fw%rmedia(:,:,:)+1d-3)
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
    if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
  end if
  
  !prepare Lorentz-Drude
  fw%num_ld=0
  do ii=0,imedia_num
    select case(type_media(ii))
    case('lorentz-drude')
      fw%num_ld=fw%num_ld+1
    end select
  end do
  if(fw%num_ld>0) then
    !set counter, make media_ld, make max_pole_num_ld, and check pole_num_ld condition
    icount_ld=1; fw%max_pole_num_ld=0;
    allocate(fw%media_ld(fw%num_ld))
    fw%media_ld(:)=0;
    do ii=0,imedia_num
      select case(type_media(ii))
      case('lorentz-drude')
        fw%media_ld(icount_ld)=ii
        icount_ld=icount_ld+1;
        if(fw%max_pole_num_ld<pole_num_ld(ii)) fw%max_pole_num_ld=pole_num_ld(ii)
        if(pole_num_ld(ii)<=0) then
          if(comm_is_root(nproc_id_global)) &
            write(*,*) "For type_media = lorentz-drude, pole_num_ld must be equal to or larger than 1."
          stop
        end if
      end select
    end do
    
    !reset counter
    icount_ld=1
    
    !allocate drude variable
    allocate(fw%idx_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%num_ld),&
             fw%idy_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%num_ld),&
             fw%idz_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%num_ld) )
    fw%idx_ld(:,:,:,:)=0; fw%idy_ld(:,:,:,:)=0; fw%idz_ld(:,:,:,:)=0;
    allocate( fw%rjx_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld),&
              fw%rjy_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld),&
              fw%rjz_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld) )
    fw%rjx_ld(:,:,:,:,:)=0.0d0; fw%rjy_ld(:,:,:,:,:)=0.0d0; fw%rjz_ld(:,:,:,:,:)=0.0d0;
    allocate( fw%rjx_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                            fs%ng%is_array(2):fs%ng%ie_array(2),&
                            fs%ng%is_array(3):fs%ng%ie_array(3)),&
              fw%rjy_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                            fs%ng%is_array(2):fs%ng%ie_array(2),&
                            fs%ng%is_array(3):fs%ng%ie_array(3)),&
              fw%rjz_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                            fs%ng%is_array(2):fs%ng%ie_array(2),&
                            fs%ng%is_array(3):fs%ng%ie_array(3)) )
    fw%rjx_sum_ld(:,:,:)=0.0d0; fw%rjy_sum_ld(:,:,:)=0.0d0; fw%rjz_sum_ld(:,:,:)=0.0d0;
    allocate( fw%px_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld),&
              fw%py_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld),&
              fw%pz_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                       fs%ng%is_array(2):fs%ng%ie_array(2),&
                       fs%ng%is_array(3):fs%ng%ie_array(3),fw%max_pole_num_ld,fw%num_ld) )
    fw%px_ld(:,:,:,:,:)=0.0d0; fw%py_ld(:,:,:,:,:)=0.0d0; fw%pz_ld(:,:,:,:,:)=0.0d0;
    allocate( fw%px_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                           fs%ng%is_array(2):fs%ng%ie_array(2),&
                           fs%ng%is_array(3):fs%ng%ie_array(3)),&
              fw%py_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                           fs%ng%is_array(2):fs%ng%ie_array(2),&
                           fs%ng%is_array(3):fs%ng%ie_array(3)),&
              fw%pz_sum_ld(fs%ng%is_array(1):fs%ng%ie_array(1),&
                           fs%ng%is_array(2):fs%ng%ie_array(2),&
                           fs%ng%is_array(3):fs%ng%ie_array(3)) )
    fw%px_sum_ld(:,:,:)=0.0d0; fw%py_sum_ld(:,:,:)=0.0d0; fw%pz_sum_ld(:,:,:)=0.0d0;
    allocate(fw%c1_j_ld(fw%max_pole_num_ld,fw%num_ld),&
             fw%c2_j_ld(fw%max_pole_num_ld,fw%num_ld),&
             fw%c3_j_ld(fw%max_pole_num_ld,fw%num_ld))
    fw%c1_j_ld(:,:)=0.0d0; fw%c2_j_ld(:,:)=0.0d0; fw%c3_j_ld(:,:)=0.0d0;
  end if
  
  !set fdtd coeffient and write media information
  allocate(fw%rep(0:imedia_num),fw%rmu(0:imedia_num),fw%sig(0:imedia_num))
  fw%rep(:)=1.0d0; fw%rmu(:)=1.0d0; fw%sig(:)=0.0d0;
  do ii=0,imedia_num
    call eh_coeff
  end do
  deallocate(fs%imedia); deallocate(fw%rmedia);
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,'(A,I3)')          ' imedia_num = ',imedia_num
    do ii=0,imedia_num
      write(*,*) "=========================="
      write(*,'(A,I3,A)')      ' id ',ii, ':'
      select case(type_media(ii))
      case ('vacuum')
        if(epsilon(ii)/=1d0 .or. rmu(ii)/=1d0 .or. sigma(ii)/=0d0) then
          write(*,'(A)'  )     ' type_media  =  constant media'
        else
          write(*,'(A,A)')     ' type_media  =  ', trim(type_media(ii))
        end if
      case('lorentz-drude')
        write(*,'(A,A)')       ' type_media  =  ', trim(type_media(ii))
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
      case default
        write(*,'(A,A)')       ' type_media  =  ', trim(type_media(ii))
      end select
      write(*,'(A,ES12.5)')    ' epsilon     = ', epsilon(ii)
      write(*,'(A,ES12.5)')    ' rmu         = ', rmu(ii)
      write(*,'(A,ES12.5)'   ) ' sigma       = ', sigma(ii)
    end do
    write(*,*) "**************************"
  end if
  
  !set calculation area
  fw%iex_y_is(:)=fs%ng%is(:); fw%iex_y_ie(:)=fs%ng%ie(:);
  fw%iex_z_is(:)=fs%ng%is(:); fw%iex_z_ie(:)=fs%ng%ie(:);
  fw%iey_z_is(:)=fs%ng%is(:); fw%iey_z_ie(:)=fs%ng%ie(:);
  fw%iey_x_is(:)=fs%ng%is(:); fw%iey_x_ie(:)=fs%ng%ie(:);
  fw%iez_x_is(:)=fs%ng%is(:); fw%iez_x_ie(:)=fs%ng%ie(:);
  fw%iez_y_is(:)=fs%ng%is(:); fw%iez_y_ie(:)=fs%ng%ie(:);
  fw%ihx_y_is(:)=fs%ng%is(:); fw%ihx_y_ie(:)=fs%ng%ie(:);
  fw%ihx_z_is(:)=fs%ng%is(:); fw%ihx_z_ie(:)=fs%ng%ie(:);
  fw%ihy_z_is(:)=fs%ng%is(:); fw%ihy_z_ie(:)=fs%ng%ie(:);
  fw%ihy_x_is(:)=fs%ng%is(:); fw%ihy_x_ie(:)=fs%ng%ie(:);
  fw%ihz_x_is(:)=fs%ng%is(:); fw%ihz_x_ie(:)=fs%ng%ie(:);
  fw%ihz_y_is(:)=fs%ng%is(:); fw%ihz_y_ie(:)=fs%ng%ie(:);
  if((fs%a_bc(1,1)=='pml').and.(fs%ng%is(1)==fs%lg%is(1))) then !x, bottom
    fw%iey_x_is(1)=fs%ng%is(1)+1; fw%iez_x_is(1)=fs%ng%is(1)+1;
  end if
  if((fs%a_bc(1,2)=='pml').and.(fs%ng%ie(1)==fs%lg%ie(1))) then !x, top
    fw%iex_y_ie(1)=fs%ng%ie(1)-1; fw%iex_z_ie(1)=fs%ng%ie(1)-1;
    fw%iey_x_ie(1)=fs%ng%ie(1)-1; fw%iez_x_ie(1)=fs%ng%ie(1)-1;
    fw%ihy_z_ie(1)=fs%ng%ie(1)-1; fw%ihy_x_ie(1)=fs%ng%ie(1)-1;
    fw%ihz_x_ie(1)=fs%ng%ie(1)-1; fw%ihz_y_ie(1)=fs%ng%ie(1)-1;
  end if
  if((fs%a_bc(2,1)=='pml').and.(fs%ng%is(2)==fs%lg%is(2))) then !y, bottom
    fw%iex_y_is(2)=fs%ng%is(2)+1; fw%iez_y_is(2)=fs%ng%is(2)+1;
  end if
  if((fs%a_bc(2,2)=='pml').and.(fs%ng%ie(2)==fs%lg%ie(2))) then !y, top
    fw%iex_y_ie(2)=fs%ng%ie(2)-1; fw%iey_z_ie(2)=fs%ng%ie(2)-1;
    fw%iey_x_ie(2)=fs%ng%ie(2)-1; fw%iez_y_ie(2)=fs%ng%ie(2)-1;
    fw%ihx_y_ie(2)=fs%ng%ie(2)-1; fw%ihx_z_ie(2)=fs%ng%ie(2)-1;
    fw%ihz_x_ie(2)=fs%ng%ie(2)-1; fw%ihz_y_ie(2)=fs%ng%ie(2)-1;
  end if
  if((fs%a_bc(3,1)=='pml').and.(fs%ng%is(3)==fs%lg%is(3))) then !z, bottom
    fw%iex_z_is(3)=fs%ng%is(3)+1; fw%iey_z_is(3)=fs%ng%is(3)+1;
  end if
  if((fs%a_bc(3,2)=='pml').and.(fs%ng%ie(3)==fs%lg%ie(3))) then !z, top
    fw%iex_z_ie(3)=fs%ng%ie(3)-1; fw%iey_z_ie(3)=fs%ng%ie(3)-1;
    fw%iez_x_ie(3)=fs%ng%ie(3)-1; fw%iez_y_ie(3)=fs%ng%ie(3)-1;
    fw%ihx_y_ie(3)=fs%ng%ie(3)-1; fw%ihx_z_ie(3)=fs%ng%ie(3)-1;
    fw%ihy_z_ie(3)=fs%ng%ie(3)-1; fw%ihy_x_ie(3)=fs%ng%ie(3)-1;
  end if
  
  !set pml
  call eh_set_pml(1,fw%c1_ey_x,fw%c2_ey_x,fw%c1_ez_x,fw%c2_ez_x,&
                    fw%c1_hy_x,fw%c2_hy_x,fw%c1_hz_x,fw%c2_hz_x) !x direction
  call eh_set_pml(2,fw%c1_ez_y,fw%c2_ez_y,fw%c1_ex_y,fw%c2_ex_y,&
                    fw%c1_hz_y,fw%c2_hz_y,fw%c1_hx_y,fw%c2_hx_y) !y direction
  call eh_set_pml(3,fw%c1_ex_z,fw%c2_ex_z,fw%c1_ey_z,fw%c2_ey_z,&
                    fw%c1_hx_z,fw%c2_hx_z,fw%c1_hy_z,fw%c2_hy_z) !z direction
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
                                 fw%coo(fs%lg%is(ii),ii)*ulength_from_au,' to ',&
                                 fw%coo(fs%lg%is(ii)+fw%ipml_l,ii)*ulength_from_au,'.'
        if(fs%a_bc(ii,2)=='pml') write(*,'(A,A,A,ES12.5,A,ES12.5,A)') &
                                 ' PML has been set for ',dir,'-direction: ',&
                                 fw%coo(fs%lg%ie(ii)-fw%ipml_l,ii)*ulength_from_au,' to ',&
                                 fw%coo(fs%lg%ie(ii),ii)*ulength_from_au,'.'
      end do
      write(*,*) "**************************"
    end if
  end if
  
  !prepare observation
  if(iobs_num_em>0) then
    !set initial
    allocate(fw%ex_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%ey_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%ez_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%hx_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%hy_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%hz_s(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ex_s(:,:,:)=0; fw%ey_s(:,:,:)=0; fw%ez_s(:,:,:)=0; 
    fw%hx_s(:,:,:)=0; fw%hy_s(:,:,:)=0; fw%hz_s(:,:,:)=0; 
    allocate(fw%iobs_po_id(iobs_num_em,3)) !1:x,        2:y,        3:z
    allocate(fw%iobs_po_pe(iobs_num_em))
    allocate(fw%iobs_li_pe(iobs_num_em,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(fw%iobs_pl_pe(iobs_num_em,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    fw%iobs_po_id(:,:)=0; fw%iobs_po_pe(:)=0; fw%iobs_li_pe(:,:)=0; fw%iobs_pl_pe(:,:)=0; 
    fw%e_max=0.0d0; fw%h_max=0.0d0;
    
    !search observation point
    do ii=1,iobs_num_em
      call find_point_em(obs_loc_em(ii,:),fw%iobs_po_id(ii,:),&
                         fw%iobs_po_pe(ii),fw%iobs_li_pe(ii,:),fw%iobs_pl_pe(ii,:),fs%ng%is(:),fs%ng%ie(:),&
                         minval(fs%lg%is)-fw%Nd,maxval(fs%lg%ie)+fw%Nd,fw%coo(:,:))
    end do
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if(iobs_num_em==1) then
        write(*,*) "Observation point is placed at"
      else
        write(*,*) "Observation points are placed at"
      end if
      do ii=1,iobs_num_em
        write(*,'(I3,A,3ES14.5)') ii,":",(fw%coo(fw%iobs_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
      end do
      write(*,*) "**************************"
      do ii=1,iobs_num_em
        write(save_name,*) ii
        save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(save_name))//'_at_point.data'
        open(fw%ifn,file=save_name)
        select case(unit_system)
        case('au','a.u.')
          write(fw%ifn,'(A)') "# time[a.u.], Ex[a.u.], Ey[a.u.], Ez[a.u.], Hx[a.u.], Hy[a.u.], Hz[a.u.]" 
        case('A_eV_fs')
          write(fw%ifn,'(A)') "# time[fs], Ex[V/Ang.], Ey[V/Ang.], Ez[V/Ang.], Hx[A/Ang.], Hy[A/Ang.], Hz[A/Ang.]" 
        end select
        close(fw%ifn)
      end do
    end if
  end if
  
  !check incident current source condition
  select case(wave_input)
  case('source')
    !linear response
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ae_shape1/2:"
        write(*,*) "For ae_shape1/2 = impulse, wave_input must be default(do not set source)."
      end if
      stop
    end if
    
    !source1
    if    (ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then 
      fw%inc_dist1='none'
    elseif(ek_dir1(1)==1.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(1)/=0.0d0.or.epdir_im1(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(1) and epdir_im1(1):"
          write(*,*) "For theory = Maxwell and ek_dir1(1) = 1.0d0, epdir_re1(1) and epdir_im1(1) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='yz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==1.0d0.and.ek_dir1(3)==0.0d0) then
      if(epdir_re1(2)/=0.0d0.or.epdir_im1(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(2) and epdir_im1(2):"
          write(*,*) "For theory = Maxwell and ek_dir1(2) = 1.0d0, epdir_re1(2) and epdir_im1(2) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='xz-plane'
      end if
    elseif(ek_dir1(1)==0.0d0.and.ek_dir1(2)==0.0d0.and.ek_dir1(3)==1.0d0) then
      if(epdir_re1(3)/=0.0d0.or.epdir_im1(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re1(3) and epdir_im1(3):"
          write(*,*) "For theory = Maxwell and ek_dir1(3) = 1.0d0, epdir_re1(3) and epdir_im1(3) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist1='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir1:"
        write(*,*) "For theory = Maxwell, ek_dir1 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
    
    !source2
    if    (ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then 
      fw%inc_dist2='none'
    elseif(ek_dir2(1)==1.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(1)/=0.0d0.or.epdir_im2(1)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(1) and epdir_im2(1):"
          write(*,*) "For theory = Maxwell and ek_dir2(1) = 1.0d0, epdir_re2(1) and epdir_im2(1) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='yz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==1.0d0.and.ek_dir2(3)==0.0d0) then
      if(epdir_re2(2)/=0.0d0.or.epdir_im2(2)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(2) and epdir_im2(2):"
          write(*,*) "For theory = Maxwell and ek_dir2(2) = 1.0d0, epdir_re2(2) and epdir_im2(2) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='xz-plane'
      end if
    elseif(ek_dir2(1)==0.0d0.and.ek_dir2(2)==0.0d0.and.ek_dir2(3)==1.0d0) then
      if(epdir_re2(3)/=0.0d0.or.epdir_im2(3)/=0.0d0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "invalid epdir_re2(3) and epdir_im2(3):"
          write(*,*) "For theory = Maxwell and ek_dir2(3) = 1.0d0, epdir_re2(3) and epdir_im2(3) must be 0.0d0."
        end if
        stop
      else
        fw%inc_dist2='xy-plane'
      end if
    else
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid ek_dir2:"
        write(*,*) "For theory = Maxwell, ek_dir2 is only allowed by"
        write(*,*) "(0d0,0d0,0d0),(1d0,0d0,0d0),(0d0,1d0,0d0),or (0d0,0d0,1d0)."
      end if
      stop
    end if
  case('point','x-line','y-line','z-line')
    !these selection are for debug
    fw%inc_dist1=wave_input; fw%inc_dist2='none';
    if(comm_is_root(nproc_id_global)) write(*,*) trim(wave_input), " source is used."
  case default
    fw%inc_dist1='none'; fw%inc_dist2='none';
    if(ae_shape1/='impulse'.and.ae_shape2/='impulse') then
      if(comm_is_root(nproc_id_global)) then
        write(*,*) "invalid wave_input:"
        write(*,*) "For theory = Maxwell, wave_input must be source"
        write(*,*) "or ae_shape1 and/or ae_shape2 must be impulse."
      end if
      stop
    end if
  end select
  
  !prepare incident current source
  if((fw%inc_dist1=='none').and.(fw%inc_dist2=='none')) then
    fw%inc_num=0
  else
    fw%inc_num=2
  end if
  if(fw%inc_num>0) then
    !set initial
    allocate(fw%inc_po_id(fw%inc_num,3)) !1:x,        2:y,        3:z
    allocate(fw%inc_po_pe(fw%inc_num))
    allocate(fw%inc_li_pe(fw%inc_num,3)) !1:x-line,   2:y-line,   3:z-line
    allocate(fw%inc_pl_pe(fw%inc_num,3)) !1:xy-plane, 2:yz-plane, 3:xz-plane
    fw%inc_po_id(:,:)=0; fw%inc_po_pe(:)=0; fw%inc_li_pe(:,:)=0; fw%inc_pl_pe(:,:)=0; 
    do ii=1,3
      fw%c2_inc_xyz(ii)=(fw%c_0/fw%rep(0)*dt_em) &
                         /(1.0d0+2.0d0*pi*fw%sig(0)/fw%rep(0)*dt_em) &
                         *2.0d0/( fs%hgs(ii)*sqrt(fw%rmu(0)/fw%rep(0)) )
    end do
    
    !search incident current source point and check others
    if(fw%inc_dist1/='none') then
      ii=1
      call find_point_em(source_loc1(:),fw%inc_po_id(ii,:),&
                         fw%inc_po_pe(ii),fw%inc_li_pe(ii,:),fw%inc_pl_pe(ii,:),fs%ng%is(:),fs%ng%ie(:),&
                         minval(fs%lg%is(:))-fw%Nd,maxval(fs%lg%ie(:))+fw%Nd,fw%coo(:,:))
      select case(ae_shape1)
      case("Ecos2","Acos2")
        continue
      case default
        if(comm_is_root(nproc_id_global)) write(*,*) 'set ae_shape1 to "Ecos2" or "Acos2".'
        stop
      end select
      diff_cep=(phi_cep1-0.25d0)*2.d0-int((phi_cep1-0.25d0)*2.d0)
      if(ae_shape1=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        if(comm_is_root(nproc_id_global)) write(*,*) &
          "phi_cep1 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape1."
        stop
      end if
      if(rlaser_int_wcm2_1/=-1d0) &
        amplitude1=sqrt(rlaser_int_wcm2_1)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    if(fw%inc_dist2/='none') then
      ii=2
      call find_point_em(source_loc2(:),fw%inc_po_id(ii,:),&
                         fw%inc_po_pe(ii),fw%inc_li_pe(ii,:),fw%inc_pl_pe(ii,:),fs%ng%is(:),fs%ng%ie(:),&
                         minval(fs%lg%is(:))-fw%Nd,maxval(fs%lg%ie(:))+fw%Nd,fw%coo(:,:))
      select case(ae_shape2)
      case("Ecos2","Acos2")
        continue
      case default
        if(comm_is_root(nproc_id_global)) write(*,*) 'set ae_shape2 to "Ecos2" or "Acos2".'
        stop
      end select
      diff_cep=(phi_cep2-0.25d0)*2.d0-int((phi_cep2-0.25d0)*2.d0)
      if(ae_shape2=="Ecos2".and.abs(diff_cep)>=1.d-12)then
        if(comm_is_root(nproc_id_global)) write(*,*) &
          "phi_cep2 must be equal to 0.25+0.5*i when Ecos2 is specified for ae_shape2."
        stop
      end if
      if(rlaser_int_wcm2_2/=-1d0) &
        amplitude2=sqrt(rlaser_int_wcm2_2)*1.0d2*2.74492d1/(5.14223d11)!I[W/cm^2]->E[a.u.]
    end if
    
    !write information
    if(comm_is_root(nproc_id_global)) then
      write(*,*)
      write(*,*) "**************************"
      if((fw%inc_dist1=='none').or.(fw%inc_dist2=='none')) then
        write(*,*) "Incident current source is placed at"
      else
        write(*,*) "Incident current sources are placed at"
      end if
      if(fw%inc_dist1/='none') then
        ii=1
        write(*,'(I8,A,3ES14.5,A)') ii,":",(fw%coo(fw%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir1:",ek_dir1
      end if
      if(fw%inc_dist2/='none') then
        ii=2
        write(*,'(I8,A,3ES14.5,A)') ii,":",(fw%coo(fw%inc_po_id(ii,ix),ix)*ulength_from_au,ix=1,3)
        write(*,'(A,3ES14.5)') " ek_dir2:",ek_dir2
      end if
      write(*,*) "**************************"
    end if
  end if
  
  !prepare linear response
  if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
    !check condition
    iflag_lr=0
    if(iperiodic==3.and.trans_longi/='tr') iflag_lr=1
    do ii=0,imedia_num
      if(fw%rep(ii)/=1.0d0.or.fw%rmu(ii)/=1.0d0.or.fw%sig(ii)/=0.0d0) iflag_lr=1
      if(ii==0) then
        select case(type_media(ii))
        case('vacuum')
          continue
        case default
          iflag_lr=1
        end select
      else
        select case(type_media(ii))
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
        write(*,*) "epsilon and rmu must be 1.0d0."
        write(*,*) "sigma must be 0.0d0."
        write(*,*) "type_media(i) must be drude, where i > 0."
        if(iperiodic==3) write(*,*) "trans_longi must be tr."
      end if
      stop
    end if
    
    !set initial current density
    if(fw%num_ld>0) then
      do ii=1,fw%num_ld
      do ij=1,pole_num_ld(fw%media_ld(ii))
        do iz=fs%ng%is(3),fs%ng%ie(3)
        do iy=fs%ng%is(2),fs%ng%ie(2)
        do ix=fs%ng%is(1),fs%ng%ie(1)
          if(fw%idx_ld(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjx_ld(ix,iy,iz,ij,ii)=fw%rjx_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
                                       /(4.0d0*pi)*e_impulse*(epdir_re1(1)+epdir_im1(1))
            end if
            if(ae_shape2=='impulse') then
              fw%rjx_ld(ix,iy,iz,ij,ii)=fw%rjx_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
                                       /(4.0d0*pi)*e_impulse*(epdir_re2(1)+epdir_im2(1))
            end if
          end if
          if(fw%idy_ld(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjy_ld(ix,iy,iz,ij,ii)=fw%rjy_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
                                       /(4.0d0*pi)*e_impulse*(epdir_re1(2)+epdir_im1(2))
            end if
            if(ae_shape2=='impulse') then
              fw%rjy_ld(ix,iy,iz,ij,ii)=fw%rjy_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
                                       /(4.0d0*pi)*e_impulse*(epdir_re2(2)+epdir_im2(2))
            end if
          end if
          if(fw%idz_ld(ix,iy,iz,ii)==1) then
            if(ae_shape1=='impulse') then
              fw%rjz_ld(ix,iy,iz,ij,ii)=fw%rjz_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
                                       /(4.0d0*pi)*e_impulse*(epdir_re1(3)+epdir_im1(3))
            end if
            if(ae_shape2=='impulse') then
              fw%rjz_ld(ix,iy,iz,ij,ii)=fw%rjz_ld(ix,iy,iz,ij,ii) &
                                       -(f_ld(fw%media_ld(ii),ij)*omega_p_ld(fw%media_ld(ii))**2.0d0) &
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
    allocate(fw%time_lr(nt_em))
    fw%time_lr(:)=0.0d0
    fw%iter_lr=1
    allocate(fw%fr_lr(0:nenergy,3),fw%fi_lr(0:nenergy,3))
    fw%fr_lr(:,:)=0.0d0; fw%fi_lr(:,:)=0.0d0;
    if(iperiodic==0) then
      allocate(fw%px_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)),&
               fw%py_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)),&
               fw%pz_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)) )
      fw%px_lr(:,:,:)=0.0d0; fw%py_lr(:,:,:)=0.0d0; fw%pz_lr(:,:,:)=0.0d0;
      allocate(fw%dip_lr(nt_em,3))
      fw%dip_lr(:,:)=0.0d0
    elseif(iperiodic==3) then
      allocate(fw%rjx_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                         fs%ng%is_array(2):fs%ng%ie_array(2),&
                         fs%ng%is_array(3):fs%ng%ie_array(3)),&
               fw%rjy_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                         fs%ng%is_array(2):fs%ng%ie_array(2),&
                         fs%ng%is_array(3):fs%ng%ie_array(3)),&
               fw%rjz_lr(fs%ng%is_array(1):fs%ng%ie_array(1),&
                         fs%ng%is_array(2):fs%ng%ie_array(2),&
                         fs%ng%is_array(3):fs%ng%ie_array(3)) )
      fw%rjx_lr(:,:,:)=0.0d0; fw%rjy_lr(:,:,:)=0.0d0; fw%rjz_lr(:,:,:)=0.0d0;
      allocate(fw%curr_lr(nt_em,3),fw%e_lr(nt_em,3))
      fw%curr_lr(:,:)=0.0d0; fw%e_lr(:,:)=0.0d0;
      allocate(fw%er_lr(0:nenergy,3),fw%ei_lr(0:nenergy,3))
      fw%er_lr(:,:)=0.0d0; fw%ei_lr(:,:)=0.0d0;
    end if
  end if
  
  !write strat
  if(comm_is_root(nproc_id_global)) then
    write(*,*)
    write(*,*) "**************************"
    write(*,*) "FDTD start"
    write(*,*) "**************************"
    write(*,*) "timestep"
    write(*,*) "-------------------------------------------------------"
  end if
  
contains
  
  !=========================================================================================
  != e and h allocation ====================================================================
  subroutine eh_allocate
    implicit none
    
    !e
    allocate(fw%ex_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ex_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ex_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ex_y(:,:,:)=0.0d0; fw%c1_ex_y(:,:,:)=0.0d0; fw%c2_ex_y(:,:,:)=0.0d0;
    allocate(fw%ex_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ex_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ex_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ex_z(:,:,:)=0.0d0; fw%c1_ex_z(:,:,:)=0.0d0; fw%c2_ex_z(:,:,:)=0.0d0;
    allocate(fw%ey_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ey_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ey_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ey_z(:,:,:)=0.0d0; fw%c1_ey_z(:,:,:)=0.0d0; fw%c2_ey_z(:,:,:)=0.0d0;
    allocate(fw%ey_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ey_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ey_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ey_x(:,:,:)=0.0d0; fw%c1_ey_x(:,:,:)=0.0d0; fw%c2_ey_x(:,:,:)=0.0d0;
    allocate(fw%ez_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ez_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ez_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ez_x(:,:,:)=0.0d0; fw%c1_ez_x(:,:,:)=0.0d0; fw%c2_ez_x(:,:,:)=0.0d0;
    allocate(fw%ez_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_ez_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_ez_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%ez_y(:,:,:)=0.0d0; fw%c1_ez_y(:,:,:)=0.0d0; fw%c2_ez_y(:,:,:)=0.0d0;
    
    !h
    allocate(fw%hx_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hx_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hx_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hx_y(:,:,:)=0.0d0; fw%c1_hx_y(:,:,:)=0.0d0; fw%c2_hx_y(:,:,:)=0.0d0;
    allocate(fw%hx_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hx_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hx_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hx_z(:,:,:)=0.0d0; fw%c1_hx_z(:,:,:)=0.0d0; fw%c2_hx_z(:,:,:)=0.0d0;
    allocate(fw%hy_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hy_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hy_z(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hy_z(:,:,:)=0.0d0; fw%c1_hy_z(:,:,:)=0.0d0; fw%c2_hy_z(:,:,:)=0.0d0;
    allocate(fw%hy_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hy_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hy_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hy_x(:,:,:)=0.0d0; fw%c1_hy_x(:,:,:)=0.0d0; fw%c2_hy_x(:,:,:)=0.0d0;
    allocate(fw%hz_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hz_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hz_x(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hz_x(:,:,:)=0.0d0; fw%c1_hz_x(:,:,:)=0.0d0; fw%c2_hz_x(:,:,:)=0.0d0;
    allocate(fw%hz_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                     fs%ng%is_array(2):fs%ng%ie_array(2),&
                     fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c1_hz_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    allocate(fw%c2_hz_y(fs%ng%is_array(1):fs%ng%ie_array(1),&
                        fs%ng%is_array(2):fs%ng%ie_array(2),&
                        fs%ng%is_array(3):fs%ng%ie_array(3)))
    fw%hz_y(:,:,:)=0.0d0; fw%c1_hz_y(:,:,:)=0.0d0; fw%c2_hz_y(:,:,:)=0.0d0;
    
    allocate(fw%c2_jx(fs%ng%is_array(1):fs%ng%ie_array(1),&
                      fs%ng%is_array(2):fs%ng%ie_array(2),&
                      fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%c2_jy(fs%ng%is_array(1):fs%ng%ie_array(1),&
                      fs%ng%is_array(2):fs%ng%ie_array(2),&
                      fs%ng%is_array(3):fs%ng%ie_array(3)),&
             fw%c2_jz(fs%ng%is_array(1):fs%ng%ie_array(1),&
                      fs%ng%is_array(2):fs%ng%ie_array(2),&
                      fs%ng%is_array(3):fs%ng%ie_array(3)) )
    fw%c2_jx(:,:,:)=0.0d0; fw%c2_jy(:,:,:)=0.0d0; fw%c2_jz(:,:,:)=0.0d0;
    
  end subroutine eh_allocate
  
  !=========================================================================================
  != set fdtd coefficient ==================================================================
  subroutine eh_coeff
    implicit none
    real(8)  :: c1_e,c2_e_x,c2_e_y,c2_e_z,c1_h,c2_h_x,c2_h_y,c2_h_z,c2_j,&
                c1_e_mid,c2_e_x_mid,c2_e_y_mid,c2_e_z_mid,c2_j_mid
    
    !set constant parameter
    fw%rep(ii)=epsilon(ii); fw%rmu(ii)=rmu(ii); fw%sig(ii)=sigma(ii);
    
    !prepare coefficient
    c1_e  =(1.0d0-2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em)
    c2_e_x=(fw%c_0/fw%rep(ii)*dt_em) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em)/fs%hgs(1)
    c2_e_y=(fw%c_0/fw%rep(ii)*dt_em) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em)/fs%hgs(2)
    c2_e_z=(fw%c_0/fw%rep(ii)*dt_em) &
           /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em)/fs%hgs(3)
    call comm_bcast(c1_e,  nproc_group_global)
    call comm_bcast(c2_e_x,nproc_group_global)
    call comm_bcast(c2_e_y,nproc_group_global)
    call comm_bcast(c2_e_z,nproc_group_global)
    c1_h=1.0d0
    c2_h_x=fw%c_0/fw%rmu(ii)*dt_em/fs%hgs(1)
    c2_h_y=fw%c_0/fw%rmu(ii)*dt_em/fs%hgs(2)
    c2_h_z=fw%c_0/fw%rmu(ii)*dt_em/fs%hgs(3)
    call comm_bcast(c1_h,  nproc_group_global)
    call comm_bcast(c2_h_x,nproc_group_global)
    call comm_bcast(c2_h_y,nproc_group_global)
    call comm_bcast(c2_h_z,nproc_group_global)
    c2_j=(4.0d0*pi/fw%rep(ii)*dt_em) &
         /(1.0d0+2.0d0*pi*fw%sig(ii)/fw%rep(ii)*dt_em)
    call comm_bcast(c2_j,nproc_group_global)
    
    !check type_media
    select case(type_media(ii))
    case('pec')
      c1_e=0.0d0; c2_e_x=0.0d0; c2_e_y=0.0d0; c2_e_z=0.0d0;
    case('lorentz-drude')
      do iz=fs%ng%is(3),fs%ng%ie(3)
      do iy=fs%ng%is(2),fs%ng%ie(2)
      do ix=fs%ng%is(1),fs%ng%ie(1)
        if(fs%imedia(ix,iy,iz)==ii) then
          if(fs%imedia(ix+1,iy,iz)==ii) then !x
            fw%idx_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
            fw%idx_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
            do ij=1,fw%num_ld
              if(fw%media_ld(ij)==fs%imedia(ix+1,iy,iz)) then
                fw%idx_ld(ix,iy,iz,ij)=1;
              end if
            end do
          end if
          if(fs%imedia(ix,iy+1,iz)==ii) then !y
            fw%idy_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
            fw%idy_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
            do ij=1,fw%num_ld
              if(fw%media_ld(ij)==fs%imedia(ix,iy+1,iz)) then
                fw%idy_ld(ix,iy,iz,ij)=1;
              end if
            end do
          end if
          if(fs%imedia(ix,iy,iz+1)==ii) then !z
            fw%idz_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
            fw%idz_ld(ix,iy,iz,icount_ld)=1;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
            do ij=1,fw%num_ld
              if(fw%media_ld(ij)==fs%imedia(ix,iy,iz+1)) then
                fw%idz_ld(ix,iy,iz,ij)=1;
              end if
            end do
          end if
        end if
      end do
      end do
      end do
      do ij=1,pole_num_ld(ii)
        fw%c1_j_ld(ij,icount_ld)=(1.0d0-gamma_ld(ii,ij)*dt_em/2.0d0) &
                                 / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
        fw%c2_j_ld(ij,icount_ld)=(f_ld(ii,ij)*(omega_p_ld(ii)**2.0d0)*dt_em/(4.0d0*pi)) &
                                 / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
        fw%c3_j_ld(ij,icount_ld)=((omega_ld(ii,ij)**2.0d0)*dt_em) &
                                 / (1.0d0+gamma_ld(ii,ij)*dt_em/2.0d0);
      end do
      icount_ld=icount_ld+1
    end select
    
    !set coefficient
    if(ii==0) then
      fw%c1_ex_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ex_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_e_y
      fw%c1_ex_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ex_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_e_z
          
      fw%c1_ey_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ey_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_e_z
      fw%c1_ey_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ey_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_e_x
        
      fw%c1_ez_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ez_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_e_x
      fw%c1_ez_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_e
      fw%c2_ez_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_e_y
        
      fw%c1_hx_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hx_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_h_y
      fw%c1_hx_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hx_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_h_z
        
      fw%c1_hy_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hy_z(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_h_z
      fw%c1_hy_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hy_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_h_x
      
      fw%c1_hz_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hz_x(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=-c2_h_x
      fw%c1_hz_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c1_h
      fw%c2_hz_y(fs%ng%is(1):fs%ng%ie(1),&
                 fs%ng%is(2):fs%ng%ie(2),&
                 fs%ng%is(3):fs%ng%ie(3))=c2_h_y
                  
      fw%c2_jx(fs%ng%is(1):fs%ng%ie(1),&
               fs%ng%is(2):fs%ng%ie(2),&
               fs%ng%is(3):fs%ng%ie(3))=-c2_j
      fw%c2_jy(fs%ng%is(1):fs%ng%ie(1),&
               fs%ng%is(2):fs%ng%ie(2),&
               fs%ng%is(3):fs%ng%ie(3))=-c2_j
      fw%c2_jz(fs%ng%is(1):fs%ng%ie(1),&
               fs%ng%is(2):fs%ng%ie(2),&
               fs%ng%is(3):fs%ng%ie(3))=-c2_j
    else
      do iz=fs%ng%is(3),fs%ng%ie(3)
      do iy=fs%ng%is(2),fs%ng%ie(2)
      do ix=fs%ng%is(1),fs%ng%ie(1)
        if(fs%imedia(ix,iy,iz)==ii) then
          !ex and jx
          if(fs%imedia(ix+1,iy,iz)==ii) then
            fw%c1_ex_y(ix,iy,iz)=c1_e; fw%c2_ex_y(ix,iy,iz)= c2_e_y;
            fw%c1_ex_z(ix,iy,iz)=c1_e; fw%c2_ex_z(ix,iy,iz)=-c2_e_z;
            fw%c2_jx(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)<ii) then
            fw%c1_ex_y(ix,iy,iz)=c1_e; fw%c2_ex_y(ix,iy,iz)= c2_e_y;
            fw%c1_ex_z(ix,iy,iz)=c1_e; fw%c2_ex_z(ix,iy,iz)=-c2_e_z;
            fw%c2_jx(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix+1,iy,iz)/=0.and.fs%imedia(ix+1,iy,iz)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em)
            c2_e_y_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /fs%hgs(2)
            c2_e_z_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /fs%hgs(3)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em)
            fw%c1_ex_y(ix,iy,iz)=c1_e_mid; fw%c2_ex_y(ix,iy,iz)= c2_e_y_mid;
            fw%c1_ex_z(ix,iy,iz)=c1_e_mid; fw%c2_ex_z(ix,iy,iz)=-c2_e_z_mid;
            fw%c2_jx(ix,iy,iz)=-c2_j_mid;
          end if
          
          !ey and jy
          if(fs%imedia(ix,iy+1,iz)==ii) then
            fw%c1_ey_z(ix,iy,iz)=c1_e; fw%c2_ey_z(ix,iy,iz)= c2_e_z;
            fw%c1_ey_x(ix,iy,iz)=c1_e; fw%c2_ey_x(ix,iy,iz)=-c2_e_x;
            fw%c2_jy(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)<ii) then
            fw%c1_ey_z(ix,iy,iz)=c1_e; fw%c2_ey_z(ix,iy,iz)= c2_e_z;
            fw%c1_ey_x(ix,iy,iz)=c1_e; fw%c2_ey_x(ix,iy,iz)=-c2_e_x;
            fw%c2_jy(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy+1,iz)/=0.and.fs%imedia(ix,iy+1,iz)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*dt_em)
            c2_e_z_mid=(fw%c_0/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /fs%hgs(3)
            c2_e_x_mid=(fw%c_0/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /fs%hgs(1)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix,iy+1,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy+1,iz))/epsilon(fs%imedia(ix,iy+1,iz))*dt_em)
            fw%c1_ey_z(ix,iy,iz)=c1_e_mid; fw%c2_ey_z(ix,iy,iz)= c2_e_z_mid;
            fw%c1_ey_x(ix,iy,iz)=c1_e_mid; fw%c2_ey_x(ix,iy,iz)=-c2_e_x_mid;
            fw%c2_jy(ix,iy,iz)=-c2_j_mid;
          end if
          
          !ez and jz
          if(fs%imedia(ix,iy,iz+1)==ii) then
            fw%c1_ez_x(ix,iy,iz)=c1_e; fw%c2_ez_x(ix,iy,iz)= c2_e_x;
            fw%c1_ez_y(ix,iy,iz)=c1_e; fw%c2_ez_y(ix,iy,iz)=-c2_e_y;
            fw%c2_jz(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)<ii) then
            fw%c1_ez_x(ix,iy,iz)=c1_e; fw%c2_ez_x(ix,iy,iz)= c2_e_x;
            fw%c1_ez_y(ix,iy,iz)=c1_e; fw%c2_ez_y(ix,iy,iz)=-c2_e_y;
            fw%c2_jz(ix,iy,iz)=-c2_j;
          elseif(fs%imedia(ix,iy,iz+1)/=0.and.fs%imedia(ix,iy,iz+1)>ii) then
            c1_e_mid  =(1.0d0-2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*dt_em)
            c2_e_x_mid=(fw%c_0/epsilon(fs%imedia(ix,iy,iz+1))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*dt_em) &
                       /fs%hgs(1)
            c2_e_y_mid=(fw%c_0/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix+1,iy,iz))/epsilon(fs%imedia(ix+1,iy,iz))*dt_em) &
                       /fs%hgs(2)
            c2_j_mid  =(4.0d0*pi/epsilon(fs%imedia(ix,iy,iz+1))*dt_em) &
                       /(1.0d0+2.0d0*pi*sigma(fs%imedia(ix,iy,iz+1))/epsilon(fs%imedia(ix,iy,iz+1))*dt_em)
            fw%c1_ez_x(ix,iy,iz)=c1_e_mid; fw%c2_ez_x(ix,iy,iz)= c2_e_x_mid;
            fw%c1_ez_y(ix,iy,iz)=c1_e_mid; fw%c2_ez_y(ix,iy,iz)=-c2_e_y_mid;
            fw%c2_jz(ix,iy,iz)=-c2_j_mid;
          end if
          
          !hx
          fw%c1_hx_y(ix,iy-1:iy,iz-1:iz)=c1_h; fw%c2_hx_y(ix,iy-1:iy,iz-1:iz)=-c2_h_y;
          fw%c1_hx_z(ix,iy-1:iy,iz-1:iz)=c1_h; fw%c2_hx_z(ix,iy-1:iy,iz-1:iz)= c2_h_z;
          
          !hy
          fw%c1_hy_z(ix-1:ix,iy,iz-1:iz)=c1_h; fw%c2_hy_z(ix-1:ix,iy,iz-1:iz)=-c2_h_z;
          fw%c1_hy_x(ix-1:ix,iy,iz-1:iz)=c1_h; fw%c2_hy_x(ix-1:ix,iy,iz-1:iz)= c2_h_x;
          
          !hz
          fw%c1_hz_x(ix-1:ix,iy-1:iy,iz)=c1_h; fw%c2_hz_x(ix-1:ix,iy-1:iy,iz)=-c2_h_x;
          fw%c1_hz_y(ix-1:ix,iy-1:iy,iz)=c1_h; fw%c2_hz_y(ix-1:ix,iy-1:iy,iz)= c2_h_y;
        end if
      end do
      end do
      end do
    end if
    
  end subroutine eh_coeff
  
  !=========================================================================================
  != set pml ===============================================================================
  subroutine eh_set_pml(idir,c1_e1,c2_e1,c1_e2,c2_e2,c1_h1,c2_h1,c1_h2,c2_h2)
    implicit none
    integer,intent(in)  :: idir
    real(8),intent(out) :: c1_e1(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c2_e1(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c1_e2(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c2_e2(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c1_h1(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c2_h1(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c1_h2(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3)),&
                           c2_h2(fs%ng%is_array(1):fs%ng%ie_array(1),&
                                 fs%ng%is_array(2):fs%ng%ie_array(2),&
                                 fs%ng%is_array(3):fs%ng%ie_array(3))
    integer :: ista,iend
    real(8) :: pml_del,s_max
    real(8) :: s_l(fw%ipml_l+1),sh_l(fw%ipml_l), &
               c1_pml(fw%ipml_l+1),c2_pml(fw%ipml_l+1),c1_pml_h(fw%ipml_l),c2_pml_h(fw%ipml_l)

    !set pml conductivity
    pml_del=fs%hgs(idir)
    s_max=-(fw%pml_m+1.0d0)*log(fw%pml_r)/(2.0d0*dble(fw%ipml_l)*pml_del) &
          *fw%c_0/(4.0d0*pi)*sqrt(fw%rep(0)/fw%rmu(0));
    do ii=1,(fw%ipml_l+1)
      s_l(ii)=s_max*(&
                    (dble(fw%ipml_l)*pml_del-(dble(ii)-1.0d0)*pml_del)/(dble(fw%ipml_l)*pml_del)&
                    )**fw%pml_m;
    end do
    do ii=1,fw%ipml_l
      sh_l(ii)=(fw%rmu(0)/fw%rep(0)) &
               *s_max*(&
                      (dble(fw%ipml_l)*pml_del-(dble(ii)-0.5d0)*pml_del)/(dble(fw%ipml_l)*pml_del)&
                      )**fw%pml_m;
    end do
    
    !set pml coefficient
    do ii=1,(fw%ipml_l+1)
      c1_pml(ii)=(1.0d0-2.0d0*pi*s_l(ii)/fw%rep(0)*dt_em) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/fw%rep(0)*dt_em);
      c2_pml(ii)=(fw%c_0/fw%rep(0)*dt_em) &
                 /(1.0d0+2.0d0*pi*s_l(ii)/fw%rep(0)*dt_em)/pml_del
    end do
    call comm_bcast(c1_pml,nproc_group_global)
    call comm_bcast(c2_pml,nproc_group_global)
    do ii=1,fw%ipml_l
      c1_pml_h(ii)=(1.0d0-2.0d0*pi*sh_l(ii)/fw%rmu(0)*dt_em) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/fw%rmu(0)*dt_em);
      c2_pml_h(ii)=(fw%c_0/fw%rmu(0)*dt_em) &
                   /(1.0d0+2.0d0*pi*sh_l(ii)/fw%rmu(0)*dt_em)/pml_del
    end do
    call comm_bcast(c1_pml_h,nproc_group_global)
    call comm_bcast(c2_pml_h,nproc_group_global)
    
    !set pml(bottom)
    if((fs%a_bc(idir,1)=='pml').and.(fs%ng%is(idir)<=(fs%lg%is(idir)+fw%ipml_l))) then
      !e
      iend=fs%lg%is(idir)+fw%ipml_l
      if(fs%ng%ie(idir)<iend) then
        iend=fs%ng%ie(idir)
      end if
      icount=1
      do ii=fs%ng%is(idir),iend
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e1(ii,:,:)=-c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_e2(ii,:,:)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e2(ii,:,:)= c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e1(:,ii,:)=-c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_e2(:,ii,:)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e2(:,ii,:)= c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e1(:,:,ii)=-c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_e2(:,:,ii)= c1_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_e2(:,:,ii)= c2_pml(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        end if
        icount=icount+1
      end do
      
      !h
      if(iend==(fs%lg%is(idir)+fw%ipml_l)) then
        iend=iend-1
      end if
      icount=1
      do ii=fs%ng%is(idir),iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h1(ii,:,:)= c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_h2(ii,:,:)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h2(ii,:,:)=-c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h1(:,ii,:)= c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_h2(:,ii,:)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h2(:,ii,:)=-c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h1(:,:,ii)= c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c1_h2(:,:,ii)= c1_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
          c2_h2(:,:,ii)=-c2_pml_h(fs%ng%is(idir)-fs%lg%is(idir)+icount)
        end if
        icount=icount+1
      end do
    end if
    
    !set pml(top)
    if((fs%a_bc(idir,2)=='pml').and.(fs%ng%ie(idir)>=(fs%lg%ie(idir)-fw%ipml_l))) then
      !e
      ista=fs%lg%ie(idir)-fw%ipml_l
      if(fs%ng%is(idir)>ista) then
        ista=fs%ng%is(idir)
      end if
      icount=1
      do ii=ista,fs%ng%ie(idir)
        if(idir==1) then
          c1_e1(ii,:,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(ii,:,:)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(ii,:,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(ii,:,:)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_e1(:,ii,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(:,ii,:)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(:,ii,:)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(:,ii,:)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_e1(:,:,ii)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e1(:,:,ii)=-c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_e2(:,:,ii)= c1_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_e2(:,:,ii)= c2_pml((fw%ipml_l+1)-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
      
      !h
      if(fs%ng%ie(idir)==fs%lg%ie(idir)) then
        iend=fs%ng%ie(idir)-1
      else
        iend=fs%ng%ie(idir)
      end if
      icount=1
      do ii=ista,iend
        if(idir==1) then
          c1_h1(ii,:,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(ii,:,:)= c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(ii,:,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(ii,:,:)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==2) then
          c1_h1(:,ii,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(:,ii,:)= c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(:,ii,:)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(:,ii,:)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        elseif(idir==3) then
          c1_h1(:,:,ii)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h1(:,:,ii)= c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c1_h2(:,:,ii)= c1_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
          c2_h2(:,:,ii)=-c2_pml_h(fw%ipml_l-(ista-(fs%lg%ie(idir)-fw%ipml_l)+(icount-1)))
        end if
        icount=icount+1
      end do
    end if
    
  end subroutine eh_set_pml
  
end subroutine eh_init

!=========================================================================================
!= input fdtd shape data =================================================================
subroutine eh_input_shape(ifn,ng_is,ng_ie,lg_is,lg_ie,Nd,imat,format)
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
  use inputoutput,          only: shape_file
  implicit none
  integer,intent(in)      :: ifn,Nd
  integer,intent(in)      :: ng_is(3),ng_ie(3),lg_is(3),lg_ie(3)
  integer,intent(out)     :: imat(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                  ng_is(2)-Nd:ng_ie(2)+Nd,&
                                  ng_is(3)-Nd:ng_ie(3)+Nd)
  character(2),intent(in) :: format
  real(8),allocatable     :: rtmp1d(:)
  integer                 :: inum(3),inum_check(3)
  integer                 :: ii,ij,ix,iy,iz,iflag_x,iflag_y,iflag_z
  
  !open file
  open(ifn,file=trim(shape_file), status='old')
  
  if(trim(format)=='cu') then
    !check grid information
    inum(:)=lg_ie(:)-lg_is(:)+1
    read(ifn,*); read (ifn,*); read (ifn,*); !skip
    allocate(rtmp1d(4))
    read (ifn,*) rtmp1d; inum_check(1)=int(rtmp1d(1)+1d-3);
    read (ifn,*) rtmp1d; inum_check(2)=int(rtmp1d(1)+1d-3);
    read (ifn,*) rtmp1d; inum_check(3)=int(rtmp1d(1)+1d-3);
    deallocate(rtmp1d)
    do ii=1,3
      if(inum(ii)/=inum_check(ii)) then
        if(comm_is_root(nproc_id_global)) write(*,*) "al_em or dl_em does not mutch shape file."
        stop
      end if
    end do
    read (ifn,*); !skip
    
    !input shape(general case)
    allocate(rtmp1d(6))
    ix=lg_is(1); iy=lg_is(2); iz=lg_is(3);
    do ii=1,int(inum(1)*inum(2)*inum(3)/6)
      read (ifn,*) rtmp1d
      do ij=1,6
        !check flag and write imat
        iflag_x=0; iflag_y=0; iflag_z=0;
        if(ix>=ng_is(1) .and. ix<=ng_ie(1)) iflag_x=1
        if(iy>=ng_is(2) .and. iy<=ng_ie(2)) iflag_y=1
        if(iz>=ng_is(3) .and. iz<=ng_ie(3)) iflag_z=1
        if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
          imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
        end if        
        
        !update iz, iy, ix 
        iz=iz+1                                          !iz
        if(iz>lg_ie(3))                      iz=lg_is(3) !iz
        if(iz==lg_is(3))                     iy=iy+1     !iy
        if(iy>lg_ie(2))                      iy=lg_is(2) !iy
        if(iz==lg_is(3) .and. iy==lg_is(2)) ix=ix+1      !ix
      end do
    end do
    deallocate(rtmp1d)
    
    !input shape(special case)
    if(mod(inum(1)*inum(2)*inum(3),6)>0) then
      allocate(rtmp1d(mod(inum(1)*inum(2)*inum(3),6)))
      read (ifn,*) rtmp1d
      do ij=1,mod(inum(1)*inum(2)*inum(3),6)
        !check flag and write imat
        iflag_x=0; iflag_y=0; iflag_z=0;
        if(ix>=ng_is(1) .and. ix<=ng_ie(1)) iflag_x=1
        if(iy>=ng_is(2) .and. iy<=ng_ie(2)) iflag_y=1
        if(iz>=ng_is(3) .and. iz<=ng_ie(3)) iflag_z=1
        if(iflag_x==1 .and. iflag_y==1 .and. iflag_z==1) then
          imat(ix,iy,iz)=int(rtmp1d(ij)+1d-3)
        end if        
        
        !update iz, iy, ix 
        iz=iz+1                                          !iz
        if(iz>lg_ie(3))                      iz=lg_is(3) !iz
        if(iz==lg_is(3))                     iy=iy+1     !iy
        if(iy>lg_ie(2))                      iy=lg_is(2) !iy
        if(iz==lg_is(3) .and. iy==lg_is(2)) ix=ix+1      !ix
      end do      
      deallocate(rtmp1d)
    end if
  elseif(trim(format)=='mp') then
  end if
  
  !close file
  close(ifn)

end subroutine eh_input_shape

!=========================================================================================
!= prepare GCEED =========================================================================
!= (This routine is temporary) ===========================================================
!= (With unifying ARTED and GCEED, this routine will be removed) =========================
subroutine eh_prep_GCEED(fs,fw)
  use inputoutput,       only: nproc_domain_orbital,nproc_domain_general,num_kgrid,nproc_k,nproc_ob,isequential,iperiodic
  use salmon_parallel,   only: nproc_id_orbitalgrid,nproc_id_global,nproc_size_global,nproc_group_global
  use set_numcpu,        only: set_numcpu_gs
  use scf_data,          only: nproc_d_o,nproc_d_g,nproc_d_o_mul,nproc_d_g_mul_dm,nproc_d_g_dm,&
                               k_sta,k_end,k_num,num_kpoints_3d,num_kpoints_rd,&
                               rLsize,Harray,Hgs,Hvol,imesh_oddeven,&
                               lg_sta,lg_end,lg_num, &
                               mg_sta,mg_end,mg_num, &
                               ng_sta,ng_end,ng_num,&
                               ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,Nd, &
                               ista_Mxin,iend_Mxin,inum_Mxin,&
                               ista_Mxin_s,iend_Mxin_s,inum_Mxin_s
  use new_world_sub,     only: make_new_world
  use init_sendrecv_sub, only: init_updown,iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
  use sendrecv_grid,     only: init_sendrecv_grid
  use structures,        only: s_fdtd_system, s_orbital_parallel
  use salmon_maxwell,    only: ls_fdtd_work
  implicit none
  type(s_fdtd_system),intent(inout) :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  integer                           :: neig_ng_eh(1:3,1:2)
  integer                           :: ii
  type(s_orbital_parallel)          :: info
  
  !set mpi condition
  num_kpoints_3d(1:3)=num_kgrid(1:3)
  num_kpoints_rd=num_kpoints_3d(1)*num_kpoints_3d(2)*num_kpoints_3d(3)
  nproc_d_o=nproc_domain_orbital
  nproc_d_g=nproc_domain_general
  call set_numcpu_gs(nproc_d_o,nproc_d_g,nproc_d_g_dm)
  nproc_d_o_mul=nproc_d_o(1)*nproc_d_o(2)*nproc_d_o(3)
  nproc_d_g_mul_dm=nproc_d_g_dm(1)*nproc_d_g_dm(2)*nproc_d_g_dm(3)
  call make_new_world(info)
  call setk(k_sta,k_end,k_num,num_kpoints_rd,nproc_k,nproc_id_orbitalgrid)
  
  !set grid and odd or even grid paterns
  rLsize(:,1)=fs%rlsize(:); Harray(:,1)=fs%hgs(:);
  Hgs(:)=Harray(:,1); Hvol=Hgs(1)*Hgs(2)*Hgs(3);
  call set_imesh_oddeven(1)
  call setlg(fs%lg,lg_sta,lg_end,lg_num,ista_Mx_ori,iend_Mx_ori,inum_Mx_ori,    &
             Hgs,Nd,rLsize(:,1),imesh_oddeven,iperiodic)
  allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1), &
           inum_Mxin(3,0:nproc_size_global-1))
  call setmg(fs%mg,mg_sta,mg_end,mg_num,ista_Mxin,iend_Mxin,inum_Mxin,  &
             lg_sta,lg_num,nproc_size_global,nproc_id_global,nproc_d_o,nproc_k,nproc_ob,isequential,1)
  allocate(ista_Mxin_s(3,0:nproc_size_global-1),iend_Mxin_s(3,0:nproc_size_global-1))
  allocate(inum_Mxin_s(3,0:nproc_size_global-1))
  call setng(fs%ng,ng_sta,ng_end,ng_num,ista_Mxin_s,iend_Mxin_s,inum_Mxin_s, &
             nproc_size_global,nproc_id_global,nproc_d_o,nproc_d_g_dm,ista_Mxin,iend_Mxin,isequential)
  fw%ioddeven(:)=imesh_oddeven(:);
  
  !set sendrecv environment
  call init_updown
  neig_ng_eh(1,1)=iup_array(2); neig_ng_eh(1,2)=idw_array(2);
  neig_ng_eh(2,1)=jup_array(2); neig_ng_eh(2,2)=jdw_array(2);
  neig_ng_eh(3,1)=kup_array(2); neig_ng_eh(3,2)=kdw_array(2);
  !This process about ng is temporal. 
  !With modifying set_ng to be applied to arbitrary Nd, this process will be removed.
  fs%ng%is_overlap(1:3)=fs%ng%is(1:3)-fw%Nd
  fs%ng%ie_overlap(1:3)=fs%ng%ie(1:3)+fw%Nd
  fs%ng%is_array(1:3)  =fs%ng%is(1:3)-fw%Nd
  fs%ng%ie_array(1:3)  =fs%ng%ie(1:3)+fw%Nd
  if(allocated(fs%ng%idx)) deallocate(fs%ng%idx)
  if(allocated(fs%ng%idy)) deallocate(fs%ng%idy)
  if(allocated(fs%ng%idz)) deallocate(fs%ng%idz)
  allocate(fs%ng%idx(fs%ng%is_overlap(1):fs%ng%ie_overlap(1)), &
           fs%ng%idy(fs%ng%is_overlap(2):fs%ng%ie_overlap(2)), &
           fs%ng%idz(fs%ng%is_overlap(3):fs%ng%ie_overlap(3)))
  do ii=fs%ng%is_overlap(1),fs%ng%ie_overlap(1)
    fs%ng%idx(ii)=ii
  end do
  do ii=fs%ng%is_overlap(2),fs%ng%ie_overlap(2)
    fs%ng%idy(ii)=ii
  end do
  do ii=fs%ng%is_overlap(3),fs%ng%ie_overlap(3)
    fs%ng%idz(ii)=ii
  end do
  fs%ng%Nd=fw%Nd
  call init_sendrecv_grid(fs%srg_ng,fs%ng,1,nproc_group_global,neig_ng_eh)
  
end subroutine eh_prep_GCEED
