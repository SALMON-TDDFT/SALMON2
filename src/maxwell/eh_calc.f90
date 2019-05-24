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
subroutine eh_calc(fs,fw)
  use inputoutput,          only: iobs_num_em,iobs_samp_em,obs_plane_em,directory,utime_from_au,t1_t2,t1_delay,&
                                  amplitude1,pulse_tw1,omega1,phi_cep1,epdir_re1,epdir_im1,ae_shape1,&
                                  amplitude2,pulse_tw2,omega2,phi_cep2,epdir_re2,epdir_im2,ae_shape2
  use salmon_parallel,      only: nproc_id_global,nproc_size_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_summation
  use structures,           only: s_fdtd_system
  use salmon_maxwell,       only: s_fdtd_work
  implicit none
  type(s_fdtd_system) :: fs
  type(s_fdtd_work)   :: fw
  integer             :: iter,ii,ix,iy,iz
  real(8),parameter   :: pi=3.141592653589793d0
  character(128)      :: save_name
  
  !time-iteration
  do iter=fw%iter_sta,fw%iter_end
    !update iter_now
    fs%iter_now=iter
    if(comm_is_root(nproc_id_global))then
      write(*,*) fs%iter_now
    end if
    
    !update drude
    if(fw%inum_d>0) then
      call eh_update_drude
    end if
    
    !calculate linear response
    if(ae_shape1=='impulse'.or.ae_shape2=='impulse') then
      call eh_calc_lr
    end if
    
    !update e
    call eh_fd(fw%iex_y_sta,fw%iex_y_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ex_y,fw%c2_ex_y,fw%ex_y,fw%hz_x,fw%hz_y,      'e','y') !ex_y
    call eh_fd(fw%iex_z_sta,fw%iex_z_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ex_z,fw%c2_ex_z,fw%ex_z,fw%hy_z,fw%hy_x,      'e','z') !ex_z
    call eh_fd(fw%iey_z_sta,fw%iey_z_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ey_z,fw%c2_ey_z,fw%ey_z,fw%hx_y,fw%hx_z,      'e','z') !ey_z
    call eh_fd(fw%iey_x_sta,fw%iey_x_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ey_x,fw%c2_ey_x,fw%ey_x,fw%hz_x,fw%hz_y,      'e','x') !ey_x
    call eh_fd(fw%iez_x_sta,fw%iez_x_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ez_x,fw%c2_ez_x,fw%ez_x,fw%hy_z,fw%hy_x,      'e','x') !ez_x
    call eh_fd(fw%iez_y_sta,fw%iez_y_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_ez_y,fw%c2_ez_y,fw%ez_y,fw%hx_y,fw%hx_z,      'e','y') !ez_y
    if(fw%inc_num>0) then                                !add incident current source
      if(fw%inc_dist1/='none') call eh_add_inc(1,amplitude1,pulse_tw1,omega1,phi_cep1,&
                                                  epdir_re1,epdir_im1,ae_shape1,fw%inc_dist1)
      if(fw%inc_dist2/='none') call eh_add_inc(2,amplitude2,pulse_tw2,omega2,phi_cep2,&
                                                  epdir_re2,epdir_im2,ae_shape2,fw%inc_dist2)
    end if
    if(fw%inum_d>0) then
      call eh_add_curr(fw%rjx_sum_d(:,:,:),fw%rjy_sum_d(:,:,:),fw%rjz_sum_d(:,:,:))
    end if
    call eh_sendrecv(fs,fw,'e')
    
    !store old h
    if( (iobs_num_em>0).and.(mod(iter,iobs_samp_em)==0) )then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=(fs%ng_sta(3)-fw%Nd),(fs%ng_end(3)+fw%Nd)
      do iy=(fs%ng_sta(2)-fw%Nd),(fs%ng_end(2)+fw%Nd)
      do ix=(fs%ng_sta(1)-fw%Nd),(fs%ng_end(1)+fw%Nd)
        fw%hx_s(ix,iy,iz)=fw%hx_y(ix,iy,iz)+fw%hx_z(ix,iy,iz)
        fw%hy_s(ix,iy,iz)=fw%hy_z(ix,iy,iz)+fw%hy_x(ix,iy,iz)
        fw%hz_s(ix,iy,iz)=fw%hz_x(ix,iy,iz)+fw%hz_y(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
    
    !update h
    call eh_fd(fw%ihx_y_sta,fw%ihx_y_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hx_y,fw%c2_hx_y,fw%hx_y,fw%ez_x,fw%ez_y,      'h','y') !hx_y
    call eh_fd(fw%ihx_z_sta,fw%ihx_z_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hx_z,fw%c2_hx_z,fw%hx_z,fw%ey_z,fw%ey_x,      'h','z') !hx_z
    call eh_fd(fw%ihy_z_sta,fw%ihy_z_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hy_z,fw%c2_hy_z,fw%hy_z,fw%ex_y,fw%ex_z,      'h','z') !hy_z
    call eh_fd(fw%ihy_x_sta,fw%ihy_x_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hy_x,fw%c2_hy_x,fw%hy_x,fw%ez_x,fw%ez_y,      'h','x') !hy_x
    call eh_fd(fw%ihz_x_sta,fw%ihz_x_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hz_x,fw%c2_hz_x,fw%hz_x,fw%ey_z,fw%ey_x,      'h','x') !hz_x
    call eh_fd(fw%ihz_y_sta,fw%ihz_y_end,      fs%ng_sta,fs%ng_end,fw%Nd,&
               fw%c1_hz_y,fw%c2_hz_y,fw%hz_y,fw%ex_y,fw%ex_z,      'h','y') !hz_y
    call eh_sendrecv(fs,fw,'h')
    
    !observation
    if( (iobs_num_em>0).and.(mod(iter,iobs_samp_em)==0) )then
      !prepare e and h for save
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=(fs%ng_sta(3)-fw%Nd),(fs%ng_end(3)+fw%Nd)
      do iy=(fs%ng_sta(2)-fw%Nd),(fs%ng_end(2)+fw%Nd)
      do ix=(fs%ng_sta(1)-fw%Nd),(fs%ng_end(1)+fw%Nd)
        fw%ex_s(ix,iy,iz)=fw%ex_y(ix,iy,iz)+fw%ex_z(ix,iy,iz)
        fw%ey_s(ix,iy,iz)=fw%ey_z(ix,iy,iz)+fw%ey_x(ix,iy,iz)
        fw%ez_s(ix,iy,iz)=fw%ez_x(ix,iy,iz)+fw%ez_y(ix,iy,iz)
        fw%hx_s(ix,iy,iz)=( fw%hx_s(ix,iy,iz)+(fw%hx_y(ix,iy,iz)+fw%hx_z(ix,iy,iz)) )/2.0d0
        fw%hy_s(ix,iy,iz)=( fw%hy_s(ix,iy,iz)+(fw%hy_z(ix,iy,iz)+fw%hy_x(ix,iy,iz)) )/2.0d0
        fw%hz_s(ix,iy,iz)=( fw%hz_s(ix,iy,iz)+(fw%hz_x(ix,iy,iz)+fw%hz_y(ix,iy,iz)) )/2.0d0
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      call eh_sendrecv(fs,fw,'s')
      
      !save data
      do ii=1,iobs_num_em
        !point
        if(fw%iobs_po_pe(ii)==1) then
          write(save_name,*) ii
          save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(save_name))//'_at_point.data'
          open(fw%ifn,file=save_name,status='old',position='append')
          write(fw%ifn, '(E13.5)',advance="no") dble(iter)*fs%dt*utime_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%ex_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uVperm_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%ey_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uVperm_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%ez_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uVperm_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%hx_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uAperm_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%hy_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uAperm_from_au
          write(fw%ifn,'(E16.6e3)',advance="no") &
                fw%hz_s(fw%iobs_po_id(ii,1),fw%iobs_po_id(ii,2),fw%iobs_po_id(ii,3))*fw%uAperm_from_au
          close(fw%ifn)
        end if
        
        !plane
        if(obs_plane_em(ii)=='y') then
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uVperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%ex_s,'ex')
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uVperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%ey_s,'ey')
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uVperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%ez_s,'ez')
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uAperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%hx_s,'hx')
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uAperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%hy_s,'hy')
          call eh_save_plane(fw%iobs_po_id(ii,:),fw%iobs_pl_pe(ii,:),fw%uAperm_from_au,&
                             fs%ng_sta,fs%ng_end,fs%lg_sta,fs%lg_end,fw%Nd,fw%ifn,ii,iter,fw%hz_s,'hz')
        end if
      end do
      
      !check maximum
      call eh_update_max
    end if
    
  end do
  
contains
  
  !=========================================================================================
  != update drude ==========================================================================
  subroutine eh_update_drude
    implicit none
    
    !initialize
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      fw%rjx_sum_d(ix,iy,iz)=0.0d0; fw%rjy_sum_d(ix,iy,iz)=0.0d0; fw%rjz_sum_d(ix,iy,iz)=0.0d0;
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !update drude current
    do ii=1,fw%inum_d
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        fw%rjx_d(ix,iy,iz,ii)= fw%c1_j_d(ii)*fw%rjx_d(ix,iy,iz,ii) &
                               +fw%c2_j_d(ii)*( fw%ex_y(ix,iy,iz)+fw%ex_z(ix,iy,iz) )&
                               *dble(fw%idx_d(ix,iy,iz,ii))
        fw%rjy_d(ix,iy,iz,ii)= fw%c1_j_d(ii)*fw%rjy_d(ix,iy,iz,ii) &
                               +fw%c2_j_d(ii)*( fw%ey_z(ix,iy,iz)+fw%ey_x(ix,iy,iz) )&
                               *dble(fw%idy_d(ix,iy,iz,ii))
        fw%rjz_d(ix,iy,iz,ii)= fw%c1_j_d(ii)*fw%rjz_d(ix,iy,iz,ii) &
                               +fw%c2_j_d(ii)*( fw%ez_x(ix,iy,iz)+fw%ez_y(ix,iy,iz) )&
                               *dble(fw%idz_d(ix,iy,iz,ii))
        fw%rjx_sum_d(ix,iy,iz)=fw%rjx_sum_d(ix,iy,iz)+fw%wex_d(ix,iy,iz,ii)*fw%rjx_d(ix,iy,iz,ii)
        fw%rjy_sum_d(ix,iy,iz)=fw%rjy_sum_d(ix,iy,iz)+fw%wey_d(ix,iy,iz,ii)*fw%rjy_d(ix,iy,iz,ii)
        fw%rjz_sum_d(ix,iy,iz)=fw%rjz_sum_d(ix,iy,iz)+fw%wez_d(ix,iy,iz,ii)*fw%rjz_d(ix,iy,iz,ii)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end do
    
  end subroutine eh_update_drude
  
  !=========================================================================================
  != calculate linear response =============================================================
  subroutine eh_calc_lr
    use inputoutput,          only: iperiodic
    use salmon_parallel,      only: nproc_group_global
    use salmon_communication, only: comm_summation
    implicit none
    real(8) :: sum_lr_x,sum_lr_y,sum_lr_z
    real(8) :: sum_lr(3),sum_lr2(3)
    
    !update time
    fw%time_lr(fw%iter_lr)=dble(fw%iter_lr)*fs%dt
    
    !initialize current density
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      fw%rjx_lr(ix,iy,iz)=0.0d0; fw%rjy_lr(ix,iy,iz)=0.0d0; fw%rjz_lr(ix,iy,iz)=0.0d0;
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !add all current density
    if(fw%inum_d>0) then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        fw%rjx_lr(ix,iy,iz)=fw%rjx_lr(ix,iy,iz)+fw%rjx_sum_d(ix,iy,iz)
        fw%rjy_lr(ix,iy,iz)=fw%rjy_lr(ix,iy,iz)+fw%rjy_sum_d(ix,iy,iz)
        fw%rjz_lr(ix,iy,iz)=fw%rjz_lr(ix,iy,iz)+fw%rjz_sum_d(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
    
    !calculate dip or curr
    if(iperiodic==0) then
!$omp parallel
!$omp do private(ix,iy,iz)
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        fw%px_lr(ix,iy,iz)=fw%px_lr(ix,iy,iz)+fw%rjx_lr(ix,iy,iz)*fs%dt
        fw%py_lr(ix,iy,iz)=fw%py_lr(ix,iy,iz)+fw%rjy_lr(ix,iy,iz)*fs%dt
        fw%pz_lr(ix,iy,iz)=fw%pz_lr(ix,iy,iz)+fw%rjz_lr(ix,iy,iz)*fs%dt
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
      sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        sum_lr_x=sum_lr_x+fw%px_lr(ix,iy,iz)
        sum_lr_y=sum_lr_y+fw%py_lr(ix,iy,iz)
        sum_lr_z=sum_lr_z+fw%pz_lr(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
      call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
      fw%dip_lr(fw%iter_lr,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3)
    elseif(iperiodic==3) then
      sum_lr_x=0.0d0;  sum_lr_y=0.0d0;  sum_lr_z=0.0d0;
      sum_lr(:)=0.0d0; sum_lr2(:)=0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_lr_x,sum_lr_y,sum_lr_z )
      do iz=fs%ng_sta(3),fs%ng_end(3)
      do iy=fs%ng_sta(2),fs%ng_end(2)
      do ix=fs%ng_sta(1),fs%ng_end(1)
        sum_lr_x=sum_lr_x+fw%rjx_lr(ix,iy,iz)
        sum_lr_y=sum_lr_y+fw%rjy_lr(ix,iy,iz)
        sum_lr_z=sum_lr_z+fw%rjz_lr(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      sum_lr(1)=sum_lr_x; sum_lr(2)=sum_lr_y; sum_lr(3)=sum_lr_z;
      call comm_summation(sum_lr,sum_lr2,3,nproc_group_global)
      fw%curr_lr(fw%iter_lr,:)=sum_lr2(:)*fs%hgs(1)*fs%hgs(2)*fs%hgs(3) &
                                 /(fs%rlsize(1)*fs%rlsize(2)*fs%rlsize(3))
    end if
    
    !update time iteration
    fw%iter_lr=fw%iter_lr+1
    
  end subroutine eh_calc_lr
  
  !=========================================================================================
  != add incident current source ===========================================================
  subroutine eh_add_inc(iord,amp,tw,omega,cep,ep_r,ep_i,aes,typ)
    implicit none
    integer,intent(in)       :: iord
    real(8),intent(in)       :: amp,tw,omega,cep
    real(8),intent(in)       :: ep_r(3),ep_i(3)
    character(16),intent(in) :: aes,typ
    real(8)                  :: t_sta,t,theta1,theta2_r,theta2_i,alpha,beta,gamma,tf_r,tf_i
    real(8)                  :: add_inc(3)
    
    !calculate time factor and adding current
    if(iord==1) then
      t_sta=t1_delay
    elseif(iord==2) then
      t_sta=t1_delay+t1_t2
    end if
    t=(dble(iter)-0.5d0)*fs%dt-t_sta
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
             -beta*cos(theta1)**2*sin(theta2_r))*gamma
      tf_i=-(-alpha*sin(2.d0*theta1)*cos(theta2_i)   &
             -beta*cos(theta1)**2*sin(theta2_i))*gamma
    end if
!    tf_r=exp(-0.5d0*(( ((dble(iter)-0.5d0)*fs%dt-10.0d0*pulse_tw1)/pulse_tw1 )**2.0d0) ) !test time factor
    add_inc(:)=amp*(tf_r*ep_r(:)+tf_i*ep_i(:))
    
    if(typ=='point') then
      if(fw%inc_po_pe(iord)==1) then
        ix=fw%inc_po_id(iord,1); iy=fw%inc_po_id(iord,2); iz=fw%inc_po_id(iord,3);
        fw%ex_y(ix,iy,iz)=add_inc(1)/2.0d0
        fw%ex_z(ix,iy,iz)=add_inc(1)/2.0d0
        fw%ey_z(ix,iy,iz)=add_inc(2)/2.0d0
        fw%ey_x(ix,iy,iz)=add_inc(2)/2.0d0
        fw%ez_x(ix,iy,iz)=add_inc(3)/2.0d0
        fw%ez_y(ix,iy,iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='x-line') then
      if(fw%inc_li_pe(iord,1)==1) then
        iy=fw%inc_po_id(iord,2); iz=fw%inc_po_id(iord,3);
        fw%ex_y(fw%iex_y_sta(1):fw%iex_y_end(1),iy,iz)=add_inc(1)/2.0d0
        fw%ex_z(fw%iex_z_sta(1):fw%iex_z_end(1),iy,iz)=add_inc(1)/2.0d0
        fw%ey_z(fw%iey_z_sta(1):fw%iey_z_end(1),iy,iz)=add_inc(2)/2.0d0
        fw%ey_x(fw%iey_x_sta(1):fw%iey_x_end(1),iy,iz)=add_inc(2)/2.0d0
        fw%ez_x(fw%iez_x_sta(1):fw%iez_x_end(1),iy,iz)=add_inc(3)/2.0d0
        fw%ez_y(fw%iez_y_sta(1):fw%iez_y_end(1),iy,iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='y-line') then
      if(fw%inc_li_pe(iord,2)==1) then
        ix=fw%inc_po_id(iord,1); iz=fw%inc_po_id(iord,3);
        fw%ex_y(ix,fw%iex_y_sta(2):fw%iex_y_end(2),iz)=add_inc(1)/2.0d0
        fw%ex_z(ix,fw%iex_z_sta(2):fw%iex_z_end(2),iz)=add_inc(1)/2.0d0
        fw%ey_z(ix,fw%iey_z_sta(2):fw%iey_z_end(2),iz)=add_inc(2)/2.0d0
        fw%ey_x(ix,fw%iey_x_sta(2):fw%iey_x_end(2),iz)=add_inc(2)/2.0d0
        fw%ez_x(ix,fw%iez_x_sta(2):fw%iez_x_end(2),iz)=add_inc(3)/2.0d0
        fw%ez_y(ix,fw%iez_y_sta(2):fw%iez_y_end(2),iz)=add_inc(3)/2.0d0
      end if
    elseif(typ=='z-line') then
      if(fw%inc_li_pe(iord,3)==1) then
        ix=fw%inc_po_id(iord,1); iy=fw%inc_po_id(iord,2);
        fw%ex_y(ix,iy,fw%iex_y_sta(3):fw%iex_y_end(3))=add_inc(1)/2.0d0
        fw%ex_z(ix,iy,fw%iex_z_sta(3):fw%iex_z_end(3))=add_inc(1)/2.0d0
        fw%ey_z(ix,iy,fw%iey_z_sta(3):fw%iey_z_end(3))=add_inc(2)/2.0d0
        fw%ey_x(ix,iy,fw%iey_x_sta(3):fw%iey_x_end(3))=add_inc(2)/2.0d0
        fw%ez_x(ix,iy,fw%iez_x_sta(3):fw%iez_x_end(3))=add_inc(3)/2.0d0
        fw%ez_y(ix,iy,fw%iez_y_sta(3):fw%iez_y_end(3))=add_inc(3)/2.0d0
      end if
    elseif(typ=='xy-plane') then !z propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(fw%inc_pl_pe(iord,1)==1) then
        iz=fw%inc_po_id(iord,3)
!$omp parallel
!$omp do private(ix,iy)
        do iy=fw%iex_z_sta(2),fw%iex_z_end(2)
        do ix=fw%iex_z_sta(1),fw%iex_z_end(1)
          fw%ex_z(ix,iy,iz)=fw%ex_z(ix,iy,iz)+fw%c2_inc_xyz(3)*add_inc(1)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy)
        do iy=fw%iey_z_sta(2),fw%iey_z_end(2)
        do ix=fw%iey_z_sta(1),fw%iey_z_end(1)
          fw%ey_z(ix,iy,iz)=fw%ey_z(ix,iy,iz)+fw%c2_inc_xyz(3)*add_inc(2)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(typ=='yz-plane') then !x propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(fw%inc_pl_pe(iord,2)==1) then
        ix=fw%inc_po_id(iord,1)
!$omp parallel
!$omp do private(iy,iz)
        do iz=fw%iey_x_sta(3),fw%iey_x_end(3)
        do iy=fw%iey_x_sta(2),fw%iey_x_end(2)
          fw%ey_x(ix,iy,iz)=fw%ey_x(ix,iy,iz)+fw%c2_inc_xyz(1)*add_inc(2)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(iy,iz)
        do iz=fw%iez_x_sta(3),fw%iez_x_end(3)
        do iy=fw%iez_x_sta(2),fw%iez_x_end(2)
          fw%ez_x(ix,iy,iz)=fw%ez_x(ix,iy,iz)+fw%c2_inc_xyz(1)*add_inc(3)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    elseif(typ=='xz-plane') then !y propagation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(fw%inc_pl_pe(iord,3)==1) then
        iy=fw%inc_po_id(iord,2)
!$omp parallel
!$omp do private(ix,iz)
        do iz=fw%iex_y_sta(3),fw%iex_y_end(3)
        do ix=fw%iex_y_sta(1),fw%iex_y_end(1)
          fw%ex_y(ix,iy,iz)=fw%ex_y(ix,iy,iz)+fw%c2_inc_xyz(2)*add_inc(1)
        end do
        end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iz)
        do iz=fw%iez_y_sta(3),fw%iez_y_end(3)
        do ix=fw%iez_y_sta(1),fw%iez_y_end(1)
          fw%ez_y(ix,iy,iz)=fw%ez_y(ix,iy,iz)+fw%c2_inc_xyz(2)*add_inc(3)
        end do
        end do
!$omp end do
!$omp end parallel
      end if
    end if
    
  end subroutine eh_add_inc
  
  !=========================================================================================
  != add current ===========================================================================
  subroutine eh_add_curr(rjx,rjy,rjz)
    implicit none
    real(8),intent(in) :: rjx(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                              fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                              fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                          rjy(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                              fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                              fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
                          rjz(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                              fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                              fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd)
    
    !ex
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iex_y_sta(3),fw%iex_y_end(3)
    do iy=fw%iex_y_sta(2),fw%iex_y_end(2)
    do ix=fw%iex_y_sta(1),fw%iex_y_end(1)
      fw%ex_y(ix,iy,iz)=fw%ex_y(ix,iy,iz)+fw%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iex_z_sta(3),fw%iex_z_end(3)
    do iy=fw%iex_z_sta(2),fw%iex_z_end(2)
    do ix=fw%iex_z_sta(1),fw%iex_z_end(1)
      fw%ex_z(ix,iy,iz)=fw%ex_z(ix,iy,iz)+fw%c2_jx(ix,iy,iz)*rjx(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !ey
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iey_z_sta(3),fw%iey_z_end(3)
    do iy=fw%iey_z_sta(2),fw%iey_z_end(2)
    do ix=fw%iey_z_sta(1),fw%iey_z_end(1)
      fw%ey_z(ix,iy,iz)=fw%ey_z(ix,iy,iz)+fw%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iey_x_sta(3),fw%iey_x_end(3)
    do iy=fw%iey_x_sta(2),fw%iey_x_end(2)
    do ix=fw%iey_x_sta(1),fw%iey_x_end(1)
      fw%ey_x(ix,iy,iz)=fw%ey_x(ix,iy,iz)+fw%c2_jy(ix,iy,iz)*rjy(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
    !ez
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iez_x_sta(3),fw%iez_x_end(3)
    do iy=fw%iez_x_sta(2),fw%iez_x_end(2)
    do ix=fw%iez_x_sta(1),fw%iez_x_end(1)
      fw%ez_x(ix,iy,iz)=fw%ez_x(ix,iy,iz)+fw%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fw%iez_y_sta(3),fw%iez_y_end(3)
    do iy=fw%iez_y_sta(2),fw%iez_y_end(2)
    do ix=fw%iez_y_sta(1),fw%iez_y_end(1)
      fw%ez_y(ix,iy,iz)=fw%ez_y(ix,iy,iz)+fw%c2_jz(ix,iy,iz)*rjz(ix,iy,iz)/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    
  end subroutine eh_add_curr
  
  !=========================================================================================
  != check and update maximum of e and h ===================================================
  subroutine eh_update_max
    implicit none
    real(8) :: fe(fs%ng_sta(1):fs%ng_end(1),fs%ng_sta(2):fs%ng_end(2),fs%ng_sta(3):fs%ng_end(3)),&
               fh(fs%ng_sta(1):fs%ng_end(1),fs%ng_sta(2):fs%ng_end(2),fs%ng_sta(3):fs%ng_end(3))
    real(8) :: e_max_tmp(0:nproc_size_global-1), h_max_tmp(0:nproc_size_global-1),&
               e_max_tmp2(0:nproc_size_global-1),h_max_tmp2(0:nproc_size_global-1)
    
    e_max_tmp(:)=0.0d0; h_max_tmp(:)=0.0d0;
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      fe(ix,iy,iz)=sqrt( fw%ex_s(ix,iy,iz)**2.0d0 + fw%ey_s(ix,iy,iz)**2.0d0 + fw%ez_s(ix,iy,iz)**2.0d0 )
      fh(ix,iy,iz)=sqrt( fw%hx_s(ix,iy,iz)**2.0d0 + fw%hy_s(ix,iy,iz)**2.0d0 + fw%hz_s(ix,iy,iz)**2.0d0 )
      if(e_max_tmp(nproc_id_global)<fe(ix,iy,iz)) e_max_tmp(nproc_id_global)=fe(ix,iy,iz)
      if(h_max_tmp(nproc_id_global)<fh(ix,iy,iz)) h_max_tmp(nproc_id_global)=fh(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(e_max_tmp,e_max_tmp2,nproc_size_global,nproc_group_global)
    call comm_summation(h_max_tmp,h_max_tmp2,nproc_size_global,nproc_group_global)
    e_max_tmp2(:)=e_max_tmp2(:)*fw%uVperm_from_au
    h_max_tmp2(:)=h_max_tmp2(:)*fw%uAperm_from_au
    if(fw%e_max<maxval(e_max_tmp2(:))) fw%e_max=maxval(e_max_tmp2(:))
    if(fw%h_max<maxval(h_max_tmp2(:))) fw%h_max=maxval(h_max_tmp2(:))
    
  end subroutine eh_update_max
  
end subroutine eh_calc

!=========================================================================================
!= calculate finite difference in eh =====================================================
subroutine eh_fd(ista,iend,ng_sta,ng_end,Nd,c1,c2,f1,f2,f3,var,dir)
  implicit none
  integer,intent(in)      :: ista(3),iend(3),ng_sta(3),ng_end(3)
  integer,intent(in)      :: Nd
  real(8),intent(in)      :: c1(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd),&
                             c2(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  real(8),intent(inout)   :: f1(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  real(8),intent(in)      :: f2(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd),&
                             f3(ng_sta(1)-Nd:ng_end(1)+Nd,&
                                ng_sta(2)-Nd:ng_end(2)+Nd,&
                                ng_sta(3)-Nd:ng_end(3)+Nd)
  character(1),intent(in) :: var,dir
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
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix-1,iy,iz)+f3(ix-1,iy,iz)) )
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
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix,iy-1,iz)+f3(ix,iy-1,iz)) )
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
                     (f2(ix,iy,iz)+f3(ix,iy,iz))-(f2(ix,iy,iz-1)+f3(ix,iy,iz-1)) )
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
                     (f2(ix+1,iy,iz)+f3(ix+1,iy,iz))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
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
                     (f2(ix,iy+1,iz)+f3(ix,iy+1,iz))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
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
                     (f2(ix,iy,iz+1)+f3(ix,iy,iz+1))-(f2(ix,iy,iz)+f3(ix,iy,iz)) )
      end do
      end do
      end do
!$omp end do
!$omp end parallel
    end if
  end if
  
end subroutine eh_fd

!=========================================================================================
!= save plane data =======================================================================
subroutine eh_save_plane(id,ipl,conv,ng_sta,ng_end,lg_sta,lg_end,Nd,ifn,iobs,iter,f,var)
  use inputoutput,          only: directory
  use salmon_parallel,      only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_summation
  implicit none
  integer,intent(in)      :: id(3),ipl(3)
  real(8),intent(in)      :: conv
  integer,intent(in)      :: ng_sta(3),ng_end(3),lg_sta(3),lg_end(3)
  integer,intent(in)      :: Nd,ifn,iobs,iter
  real(8),intent(in)      :: f(ng_sta(1)-Nd:ng_end(1)+Nd,&
                               ng_sta(2)-Nd:ng_end(2)+Nd,&
                               ng_sta(3)-Nd:ng_end(3)+Nd)
  character(2),intent(in) :: var
  real(8),allocatable     :: save_pl(:,:),save_pl2(:,:)
  integer          :: ii,inum,i1,i1s,i2,i2s
  character(2)     :: plane_name
  character(128)   :: iobs_name,iter_name,save_name
  
  do ii=1,3
    !allocate
    if(ii==1) then     !xy
      i1s=1; i2s=2; plane_name='xy';
    elseif(ii==2) then !yz
      i1s=2; i2s=3; plane_name='yz';
    elseif(ii==3) then !xz
      i1s=1; i2s=3; plane_name='xz';
    end if
    allocate(save_pl(lg_sta(i1s):lg_end(i1s),lg_sta(i2s):lg_end(i2s)),&
             save_pl2(lg_sta(i1s):lg_end(i1s),lg_sta(i2s):lg_end(i2s)))
    save_pl(:,:)=0.0d0; save_pl2(:,:)=0.0d0
    inum=(lg_end(i1s)-lg_sta(i1s)+1)*(lg_end(i2s)-lg_sta(i2s)+1)
    
    !prepare save data
    if(ipl(ii)==1) then
      do i2=ng_sta(i2s),ng_end(i2s)
      do i1=ng_sta(i1s),ng_end(i1s)
        if(ii==1) then     !xy
          save_pl(i1,i2)=f(i1,i2,id(3))
        elseif(ii==2) then !yz
          save_pl(i1,i2)=f(id(1),i1,i2)
        elseif(ii==3) then !xz
          save_pl(i1,i2)=f(i1,id(2),i2)
        end if
      end do
      end do
    end if
    call comm_summation(save_pl,save_pl2,inum,nproc_group_global)
    
    !make save data
    if(comm_is_root(nproc_id_global)) then
      write(iobs_name,*) iobs
      write(iter_name,*) iter
      save_name=trim(adjustl(directory))//'/obs'//trim(adjustl(iobs_name))//'_'//var//&
                '_'//plane_name//'_'//trim(adjustl(iter_name))//'.data'
      open(ifn,file=save_name)
      do i2=lg_sta(i2s),lg_end(i2s)
      do i1=lg_sta(i1s),lg_end(i1s)
        write(ifn,'(I8,I8,E16.6e3)') i1,i2,save_pl2(i1,i2)*conv
      end do
      end do
      close(ifn)
    end if
    
    !deallocate
    deallocate(save_pl,save_pl2)
  end do
  
end subroutine eh_save_plane

!=========================================================================================
!= send and receive eh ===================================================================
!= (This routine is temporary) ===========================================================
!= (With unifying ARTED and GCEED, this routine will be removed) =========================
subroutine eh_sendrecv(fs,fw,var)
  use sendrecv_grid,  only: update_overlap_real8
  use structures,     only: s_fdtd_system
  use salmon_maxwell, only: s_fdtd_work
  implicit none
  type(s_fdtd_system)     :: fs
  type(s_fdtd_work)       :: fw
  character(1),intent(in) :: var
  integer                 :: ix,iy,iz
  real(8),allocatable     :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
  
  if(var=='e') then
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ex_y)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ex_z)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ey_z)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ey_x)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ez_x)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%ez_y)
  elseif(var=='h') then
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hx_y)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hx_z)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hy_z)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hy_x)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hz_x)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hz_y)
  elseif(var=='r') then
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%rmedia)
  elseif(var=='s') then
    !allocate temporary variable
    allocate(f1(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             f2(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd),&
             f3(fs%ng_sta(1)-fw%Nd:fs%ng_end(1)+fw%Nd,&
                fs%ng_sta(2)-fw%Nd:fs%ng_end(2)+fw%Nd,&
                fs%ng_sta(3)-fw%Nd:fs%ng_end(3)+fw%Nd))
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    
    !spatially adjust e for save
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      f1(ix,iy,iz)=( fw%ex_s(ix,iy,iz)+fw%ex_s(ix-1,iy,iz) )/2.0d0
      f2(ix,iy,iz)=( fw%ey_s(ix,iy,iz)+fw%ey_s(ix,iy-1,iz) )/2.0d0
      f3(ix,iy,iz)=( fw%ez_s(ix,iy,iz)+fw%ez_s(ix,iy,iz-1) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    fw%ex_s(:,:,:)=f1(:,:,:); fw%ey_s(:,:,:)=f2(:,:,:); fw%ez_s(:,:,:)=f3(:,:,:);
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    
    !spatially adjust h for save
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      f1(ix,iy,iz)=( fw%hx_s(ix,iy,iz)+fw%hx_s(ix,iy-1,iz) )/2.0d0
      f2(ix,iy,iz)=( fw%hy_s(ix,iy,iz)+fw%hy_s(ix,iy,iz-1) )/2.0d0
      f3(ix,iy,iz)=( fw%hz_s(ix,iy,iz)+fw%hz_s(ix-1,iy,iz) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    fw%hx_s(:,:,:)=f1(:,:,:); fw%hy_s(:,:,:)=f2(:,:,:); fw%hz_s(:,:,:)=f3(:,:,:);
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hx_s)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hy_s)
    call update_overlap_real8(fs%srg_ng,fs%ng,fw%hz_s)
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng_sta(3),fs%ng_end(3)
    do iy=fs%ng_sta(2),fs%ng_end(2)
    do ix=fs%ng_sta(1),fs%ng_end(1)
      f1(ix,iy,iz)=( fw%hx_s(ix,iy,iz)+fw%hx_s(ix,iy,iz-1) )/2.0d0
      f2(ix,iy,iz)=( fw%hy_s(ix,iy,iz)+fw%hy_s(ix-1,iy,iz) )/2.0d0
      f3(ix,iy,iz)=( fw%hz_s(ix,iy,iz)+fw%hz_s(ix,iy-1,iz) )/2.0d0
    end do
    end do
    end do
!$omp end do
!$omp end parallel
    fw%hx_s(:,:,:)=f1(:,:,:); fw%hy_s(:,:,:)=f2(:,:,:); fw%hz_s(:,:,:)=f3(:,:,:);
    
    !deallocate temporary variable
    deallocate(f1,f2,f3)
  end if
  
end subroutine eh_sendrecv
