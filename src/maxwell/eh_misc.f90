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

!=========================================================================================
!= prepare mpi, grid, and sendrecv enviroments============================================
subroutine eh_mpi_grid_sr(fs,fw)
  use salmon_global,     only: nproc_domain_orbital,nproc_domain_general,yn_periodic, &
                               nproc_k,nproc_ob
  use salmon_parallel,   only: nproc_group_global
  use set_numcpu,        only: set_numcpu_general,iprefer_domain_distribution
  use init_communicator, only: init_communicator_dft
  use sendrecv_grid,     only: create_sendrecv_neig_ng,init_sendrecv_grid
  use structures,        only: s_fdtd_system, s_orbital_parallel, s_field_parallel, s_process_info
  use salmon_maxwell,    only: ls_fdtd_work
  use initialization_sub
  implicit none
  type(s_fdtd_system),intent(inout) :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  type(s_process_info)              :: pinfo
  type(s_orbital_parallel)          :: info
  type(s_field_parallel)            :: info_field
  integer                           :: neig_ng_eh(1:2,1:3)
  integer                           :: ii,iperi
  
  !set mpi condition
  pinfo%npk              = nproc_k
  pinfo%nporbital        = nproc_ob
  pinfo%npdomain_orbital = nproc_domain_orbital
  pinfo%npdomain_general = nproc_domain_general
  call set_numcpu_general(iprefer_domain_distribution,1,1,nproc_group_global,pinfo)
  call init_communicator_dft(nproc_group_global,pinfo,info,info_field)
  
  !initialize r-grid
  call init_grid_whole(fs%rlsize,fs%hgs,fs%lg)
  call init_grid_parallel(info%id_rko,info%isize_rko,pinfo,fs%lg,fs%mg,fs%ng) ! lg --> mg & ng
  !### This process about ng is temporal. #####################!
  !### With modifying set_ng to be applied to arbitrary Nd, ###!
  !### this process will be removed.###########################!
  fs%ng%is_overlap(1:3)=fs%ng%is(1:3)-fw%Nd
  fs%ng%ie_overlap(1:3)=fs%ng%ie(1:3)+fw%Nd
  fs%ng%is_array(1:3)  =fs%ng%is(1:3)-fw%Nd
  fs%ng%ie_array(1:3)  =fs%ng%ie(1:3)+fw%Nd
  !############################################################!
  
  !prepare for setting sendrecv environment
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
  
  !set sendrecv environment
  if    (yn_periodic=='n') then
    iperi=0
  elseif(yn_periodic=='y') then
    iperi=3
  end if
  call create_sendrecv_neig_ng(neig_ng_eh,pinfo,info_field,iperi) ! neighboring node array
  call init_sendrecv_grid(fs%srg_ng,fs%ng,1,info_field%icomm_all,neig_ng_eh)
  
end subroutine eh_mpi_grid_sr

!=========================================================================================
!= input fdtd shape data =================================================================
subroutine eh_input_shape(ifn,ng_is,ng_ie,lg_is,lg_ie,Nd,imat,format)
  use salmon_global,        only: shape_file
  use salmon_parallel,      only: nproc_id_global
  use salmon_communication, only: comm_is_root
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
        if(iz==lg_is(3) .and. iy==lg_is(2))  ix=ix+1     !ix
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
!= send and receive eh ===================================================================
subroutine eh_sendrecv(fs,fw,var)
  use sendrecv_grid,  only: update_overlap_real8
  use structures,     only: s_fdtd_system
  use salmon_maxwell, only: ls_fdtd_work
  implicit none
  type(s_fdtd_system),intent(inout) :: fs
  type(ls_fdtd_work), intent(inout) :: fw
  character(1),intent(in)           :: var
  integer                           :: ix,iy,iz
  real(8),allocatable               :: f1(:,:,:),f2(:,:,:),f3(:,:,:)
  
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
    allocate(f1(fs%ng%is_array(1):fs%ng%ie_array(1),&
                fs%ng%is_array(2):fs%ng%ie_array(2),&
                fs%ng%is_array(3):fs%ng%ie_array(3)),&
             f2(fs%ng%is_array(1):fs%ng%ie_array(1),&
                fs%ng%is_array(2):fs%ng%ie_array(2),&
                fs%ng%is_array(3):fs%ng%ie_array(3)),&
             f3(fs%ng%is_array(1):fs%ng%ie_array(1),&
                fs%ng%is_array(2):fs%ng%ie_array(2),&
                fs%ng%is_array(3):fs%ng%ie_array(3)))
    f1(:,:,:)=0.0d0; f2(:,:,:)=0.0d0; f3(:,:,:)=0.0d0;
    
    !spatially adjust e for save
!$omp parallel
!$omp do private(ix,iy,iz)
    do iz=fs%ng%is(3),fs%ng%ie(3)
    do iy=fs%ng%is(2),fs%ng%ie(2)
    do ix=fs%ng%is(1),fs%ng%ie(1)
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
    do iz=fs%ng%is(3),fs%ng%ie(3)
    do iy=fs%ng%is(2),fs%ng%ie(2)
    do ix=fs%ng%is(1),fs%ng%ie(1)
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
    do iz=fs%ng%is(3),fs%ng%ie(3)
    do iy=fs%ng%is(2),fs%ng%ie(2)
    do ix=fs%ng%is(1),fs%ng%ie(1)
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

!=========================================================================================
!= calculate finite difference in eh =====================================================
subroutine eh_fd(ista,iend,ng_is,ng_ie,Nd,c1,c2,f1,f2,f3,var,dir)
  implicit none
  integer,intent(in)      :: ista(3),iend(3),ng_is(3),ng_ie(3)
  integer,intent(in)      :: Nd
  real(8),intent(in)      :: c1(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                ng_is(2)-Nd:ng_ie(2)+Nd,&
                                ng_is(3)-Nd:ng_ie(3)+Nd),&
                             c2(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                ng_is(2)-Nd:ng_ie(2)+Nd,&
                                ng_is(3)-Nd:ng_ie(3)+Nd)
  real(8),intent(inout)   :: f1(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                ng_is(2)-Nd:ng_ie(2)+Nd,&
                                ng_is(3)-Nd:ng_ie(3)+Nd)
  real(8),intent(in)      :: f2(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                ng_is(2)-Nd:ng_ie(2)+Nd,&
                                ng_is(3)-Nd:ng_ie(3)+Nd),&
                             f3(ng_is(1)-Nd:ng_ie(1)+Nd,&
                                ng_is(2)-Nd:ng_ie(2)+Nd,&
                                ng_is(3)-Nd:ng_ie(3)+Nd)
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
subroutine eh_save_plane(id,ipl,conv,ng_is,ng_ie,lg_is,lg_ie,Nd,ifn,iobs,iter,f,var)
  use salmon_global,        only: base_directory
  use salmon_parallel,      only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root,comm_summation
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
    if(ii==1) then     !xy
      i1s=1; i2s=2; plane_name='xy';
    elseif(ii==2) then !yz
      i1s=2; i2s=3; plane_name='yz';
    elseif(ii==3) then !xz
      i1s=1; i2s=3; plane_name='xz';
    end if
    allocate(save_pl(lg_is(i1s):lg_ie(i1s),lg_is(i2s):lg_ie(i2s)),&
             save_pl2(lg_is(i1s):lg_ie(i1s),lg_is(i2s):lg_ie(i2s)))
    save_pl(:,:)=0.0d0; save_pl2(:,:)=0.0d0
    inum=(lg_ie(i1s)-lg_is(i1s)+1)*(lg_ie(i2s)-lg_is(i2s)+1)
    
    !prepare save data
    if(ipl(ii)==1) then
      do i2=ng_is(i2s),ng_ie(i2s)
      do i1=ng_is(i1s),ng_ie(i1s)
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
  
end subroutine eh_save_plane

!=========================================================================================
!= Fourier transformation in eh ==========================================================
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
  
end subroutine eh_fourier
