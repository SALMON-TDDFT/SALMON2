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
!=======================================================================

#include "config.h"

module checkpoint_restart_sub
  implicit none

  integer,parameter,private :: write_mode = 1
  integer,parameter,private :: read_mode  = 2

contains

!===================================================================================================================================

subroutine init_dir_out_restart(ofl)
  use structures,  only: s_ofile
  use filesystem,  only: atomic_create_directory
  use salmon_global,   only: theory,write_rt_wfn_k
  use parallelization, only: nproc_id_global,nproc_group_global
  implicit none
  type(s_ofile), intent(inout) :: ofl

  select case(theory)
    case('dft','dft_md','dft_band')
      ofl%dir_out_restart = 'data_for_restart/'
      call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
    case('tddft_response','tddft_pulse','single_scale_maxwell_tddft')
      if (write_rt_wfn_k == 'y') then
        ofl%dir_out_restart = 'data_for_restart_rt/'
        call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
      end if
    case('multi_scale_maxwell_tddft')
      if (write_rt_wfn_k == 'y') then
        ofl%dir_out_restart = 'data_for_restart_ms'
        call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
      end if
    case('dft2tddft')
      ofl%dir_out_restart = 'data_for_restart_rt/'
      call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
  end select

end subroutine init_dir_out_restart


subroutine generate_checkpoint_directory_name(header,iter,gdir,pdir)
  use parallelization, only: nproc_id_global
  implicit none
  character(*),  intent(in)  :: header
  integer,       intent(in)  :: iter
  character(256),intent(out) :: gdir
  character(256),intent(out) :: pdir

  ! global directory
  write(gdir,'(A,A,A,I6.6,A)') "checkpoint_",trim(header),"_",iter,"/"
  ! process private directory
  write(pdir,'(A,A,I6.6,A)')   trim(gdir),'rank_',nproc_id_global,'/'
end subroutine generate_checkpoint_directory_name

subroutine generate_restart_directory_name(basedir,gdir,pdir)
  use parallelization, only: nproc_id_global
  implicit none
  character(*),  intent(in)  :: basedir
  character(256),intent(out) :: gdir
  character(256),intent(out) :: pdir

  ! global directory
  write(gdir,'(A,I6.6,A)')   trim(basedir)
  ! process private directory
  write(pdir,'(A,A,I6.6,A)') trim(gdir),'rank_',nproc_id_global,'/'
end subroutine generate_restart_directory_name


subroutine checkpoint_gs(lg,mg,system,info,spsi,iter,mixing,odir)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital, s_mixing
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint
  use parallelization, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_mixing)          ,intent(in) :: mixing
  character(*),optional   ,intent(in) :: odir

  character(256) :: gdir,wdir
  logical :: iself

  if (present(odir)) then
    ! save wavefunction
    call generate_restart_directory_name(odir,gdir,wdir)
    iself = .false.
  else
    call generate_checkpoint_directory_name('gs',iter,gdir,wdir)
    call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)
    iself = (yn_self_checkpoint == 'y')
  end if

  if (iself) then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_Rion(wdir,system)
  call write_Velocity(wdir,system)
  call write_bin(wdir,lg,mg,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine checkpoint_gs

subroutine restart_gs(lg,mg,system,info,spsi,iter,mixing)
  use structures, only: s_rgrid, s_dft_system,s_parallel_info, s_orbital, s_mixing, s_mixing
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint,&
                           read_gs_restart_data
  implicit none
  type(s_rgrid)             ,intent(in) :: lg, mg
  type(s_dft_system)     ,intent(inout) :: system
  type(s_parallel_info)  ,intent(in) :: info
  type(s_orbital)        ,intent(inout) :: spsi
  integer                  ,intent(out) :: iter
  type(s_mixing)         ,intent(inout) :: mixing

  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (.not. iself) then
    wdir = gdir
  end if
  if (read_gs_restart_data=='rho' .or. &
      read_gs_restart_data=='rho_inout:single' .or. &
      read_gs_restart_data=='all:single'  ) then
    wdir = gdir
    iself = .false. !turn off self_checkpoint format in reading and behave like single file format
  end if

  call read_bin(wdir,lg,mg,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine restart_gs

subroutine checkpoint_opt(Miopt,opt,odir)
  use structures, only: s_opt
  use salmon_global, only: natom
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root
  use filesystem, only: get_filehandle, atomic_create_directory
  implicit none
  type(s_opt) :: opt
  integer :: Miopt,fh_opt_bin, NA3
  character(256) :: dir_file_out, gdir,wdir
  character(*),optional   ,intent(in) :: odir

  if (present(odir)) then
    call generate_restart_directory_name(odir,gdir,wdir)
  else
    call generate_checkpoint_directory_name('gs',Miopt,gdir,wdir)
    call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)
  end if

  dir_file_out = trim(gdir)//"opt.bin"
  NA3 = 3*natom

  if (comm_is_root(nproc_id_global))then
     fh_opt_bin = get_filehandle()
     open(fh_opt_bin,file=dir_file_out,form='unformatted')
     write(fh_opt_bin) Miopt
     write(fh_opt_bin) opt%a_dRion(1:NA3)
     write(fh_opt_bin) opt%dFion(1:NA3)
     write(fh_opt_bin) opt%Hess_mat(1:NA3,1:NA3)
     write(fh_opt_bin) opt%Hess_mat_last(1:NA3,1:NA3)
     close(fh_opt_bin)
  endif

end subroutine checkpoint_opt

subroutine restart_opt(Miopt,opt)
  use structures, only: s_opt
  use salmon_global, only: directory_read_data,natom
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root,comm_bcast
  use filesystem, only: get_filehandle
  implicit none
  type(s_opt) :: opt
  integer :: access
  integer :: Miopt,istat, fh_opt_bin, NA3, comm
  character(256) :: dir_file_in, gdir,wdir
  
  call generate_restart_directory_name(directory_read_data,gdir,wdir)
  dir_file_in = trim(gdir)//"opt.bin"

  NA3 = 3*natom

  if (comm_is_root(nproc_id_global))then
     istat = access( dir_file_in, " ")
     if(istat==0) then  !if file exists
        fh_opt_bin = get_filehandle()
        open(fh_opt_bin,file=dir_file_in,form='unformatted')
        read(fh_opt_bin) Miopt
        read(fh_opt_bin) opt%a_dRion(1:NA3)
        read(fh_opt_bin) opt%dFion(1:NA3)
        read(fh_opt_bin) opt%Hess_mat(1:NA3,1:NA3)
        read(fh_opt_bin) opt%Hess_mat_last(1:NA3,1:NA3)
        close(fh_opt_bin)
     endif
  endif

  comm = nproc_group_global
  call comm_bcast(Miopt             ,comm)
  call comm_bcast(opt%a_dRion       ,comm)
  call comm_bcast(opt%dFion         ,comm)
  call comm_bcast(opt%Hess_mat      ,comm)
  call comm_bcast(opt%Hess_mat_last ,comm)

end subroutine restart_opt

subroutine checkpoint_rt(lg,mg,system,info,spsi,iter,rt,Vh_stock1,Vh_stock2,singlescale,idir)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital, s_scalar, s_singlescale, s_rt
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint
  use parallelization, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info)   ,intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_rt)              ,intent(in) :: rt
  type(s_scalar)          ,intent(in) :: Vh_stock1,Vh_stock2
  type(s_singlescale)     ,intent(in) :: singlescale
  character(*),optional   ,intent(in) :: idir

  character(256) :: gdir,wdir
  logical :: iself

  if (present(idir)) then
    ! save wavefunction
    call generate_restart_directory_name(idir,gdir,wdir)
    iself = .false.
  else
    call generate_checkpoint_directory_name('rt',iter,gdir,wdir)
    call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)
    iself = (yn_self_checkpoint == 'y')
  end if

  if (iself) then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_Rion(wdir,system)
  call write_Velocity(wdir,system)
  call write_bin(wdir,lg,mg,system,info,spsi,iter &
                ,Vh_stock1=Vh_stock1,Vh_stock2=Vh_stock2,is_self_checkpoint=iself)
  if(singlescale%flag_use) then
    call write_singlescale(wdir,lg,mg,info,singlescale,system%Ac_micro,system%div_Ac,is_self_checkpoint=iself)
  end if

  call write_rtdata(wdir,iter,lg,mg,system,info,iself,rt)
  
end subroutine checkpoint_rt

subroutine restart_rt(lg,mg,system,info,spsi,iter,rt,Vh_stock1,Vh_stock2)
  use structures, only: s_rgrid, s_dft_system,s_parallel_info, s_orbital, s_mixing, s_scalar, s_rt
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint
  implicit none
  type(s_rgrid)          ,intent(in)    :: lg, mg
  type(s_dft_system)     ,intent(inout) :: system
  type(s_parallel_info)  ,intent(in)    :: info
  type(s_orbital)        ,intent(inout) :: spsi
  integer                ,intent(out)   :: iter
  type(s_rt)             ,intent(inout) :: rt
  type(s_scalar)         ,intent(inout) :: Vh_stock1,Vh_stock2

  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = yn_restart =='y' .and. yn_self_checkpoint == 'y'
  if (.not. iself) then
    wdir = gdir
  end if

  call read_bin(wdir,lg,mg,system,info,spsi,iter &
               ,Vh_stock1=Vh_stock1,Vh_stock2=Vh_stock2,is_self_checkpoint=iself)
       
  if(yn_restart =='y') then
    call read_rtdata(wdir,iter,lg,mg,system,info,iself,rt)
  end if
  
end subroutine restart_rt

!===================================================================================================================================

subroutine write_bin(odir,lg,mg,system,info,spsi,iter,mixing,Vh_stock1,Vh_stock2,is_self_checkpoint)
  use salmon_global, only: theory,calc_mode,write_gs_restart_data
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital, s_mixing, s_scalar
  use parallelization, only: nproc_id_global, nproc_size_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)            ,intent(in) :: odir
  type(s_rgrid)           ,intent(in) :: lg, mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_mixing) ,optional,intent(in) :: mixing
  type(s_scalar) ,optional,intent(in) :: Vh_stock1,Vh_stock2
  logical        ,optional,intent(in) :: is_self_checkpoint

  logical :: flag_GS, flag_RT
  integer :: iu1_w
  character(100) :: dir_file_out
  logical :: iself

  flag_GS = (theory=='dft'.or.theory=='dft_md'.or.calc_mode=='GS')
  flag_RT = (theory=='tddft_response'.or.theory=='tddft_pulse'.or.calc_mode=='RT')

  if (present(is_self_checkpoint)) then
    iself = is_self_checkpoint
  else
    iself = .false.
  end if

  iu1_w = 97

  !information
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"info.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')

    write(iu1_w) system%nk
    write(iu1_w) system%no

    write(iu1_w) iter               ! iteration number (same format for gs and rt calculations)
    write(iu1_w) nproc_size_global  ! number of process (to debug)
    
    write(iu1_w) system%if_real_orbital

    close(iu1_w)
  end if

  !occupation
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"occupation.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')
    write(iu1_w) system%rocc(1:system%no,1:system%nk,1:system%nspin)
    close(iu1_w)
  end if

  !wave fucntion
  if( flag_GS .and. &
      (write_gs_restart_data=='rho_inout'.or.write_gs_restart_data=='rho'))then
     ! this case => do not write wavefunction
  else
     call write_wavefunction(odir,lg,mg,system,info,spsi,iself)
  endif

  !rho_inout
  if( flag_GS.and. &
     (write_gs_restart_data=='rho'.or.write_gs_restart_data=='wfn'))then
     ! this case => do not write rho_inout
  else
     if (present(mixing)) then
        call write_rho_inout(odir,lg,mg,system,info,mixing,iself)
     end if
  end if

  !rho (only for GS)
  if( flag_GS .and. write_gs_restart_data=='rho' )then
     call write_rho(odir,lg,mg,system,info,mixing)
  endif

  !Vh_stock
  if( flag_RT )then
    if (present(Vh_stock1) .and. present(Vh_stock2)) then
      call write_Vh_stock(odir,lg,mg,info,Vh_stock1,Vh_stock2,iself)
    end if
  end if

end subroutine write_bin

!===================================================================================================================================

subroutine read_bin(idir,lg,mg,system,info,spsi,iter,mixing,Vh_stock1,Vh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system,s_parallel_info, s_orbital, s_mixing, s_scalar
  use parallelization, only: nproc_id_global,nproc_group_global,nproc_size_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: yn_restart, theory,calc_mode,read_gs_restart_data, yn_reset_step_restart
  implicit none
  character(*)              ,intent(in) :: idir
  type(s_rgrid)             ,intent(in) :: lg, mg
  type(s_dft_system)     ,intent(inout) :: system
  type(s_parallel_info)  ,intent(in) :: info
  type(s_orbital)        ,intent(inout) :: spsi
  integer                  ,intent(out) :: iter
  type(s_mixing),optional,intent(inout) :: mixing
  type(s_scalar),optional,intent(inout) :: Vh_stock1,Vh_stock2
  logical       ,optional,intent(in)    :: is_self_checkpoint

  logical :: flag_GS, flag_RT
  logical :: flag_read_info, flag_read_occ
  integer :: mk,mo
  real(8),allocatable :: roccbox(:,:,:)
  integer :: iu1_r
  character(256) :: dir_file_in
  integer :: comm,itt,nprocs
  logical :: iself,if_real_orbital

  flag_GS = (theory=='dft'.or.theory=='dft_md'.or.theory=='dft_band'.or.calc_mode=='GS')
  flag_RT = (theory=='tddft_response'.or.theory=='tddft_pulse'.or.calc_mode=='RT')

  flag_read_info = .true.
  flag_read_occ  = .true.
  if( flag_GS ) then
     if(read_gs_restart_data=='rho'.or.read_gs_restart_data(1:9)=='rho_inout')then
        flag_read_info = .false.
        flag_read_occ  = .false.
     endif
  endif

  if (present(is_self_checkpoint)) then
    iself = is_self_checkpoint
  else
    iself = .false.
  end if

  comm = nproc_group_global

  iu1_r = 96

  !information
  !first to be read
  if(flag_read_info) then
     if(comm_is_root(nproc_id_global))then
        dir_file_in = trim(idir)//"info.bin"
        open(iu1_r,file=dir_file_in,form='unformatted')

        read(iu1_r) mk
        read(iu1_r) mo

        read(iu1_r) itt
        read(iu1_r) nprocs
        
        read(iu1_r) if_real_orbital

        close(iu1_r)
     end if
     call comm_bcast(mk,comm)
     call comm_bcast(mo,comm)
     call comm_bcast(itt,comm)
     call comm_bcast(nprocs,comm)
     call comm_bcast(if_real_orbital,comm)

     if((theory=='dft'.or.theory=='dft_band'.or.calc_mode=='GS').or.  &
        ((theory=='tddft_response'.or.theory=='tddft_pulse'.or.calc_mode=='RT').and.yn_restart=='y'))then
        iter = itt
        if ( theory == 'dft_band' ) iter=0
     end if

     !debug check
     if (yn_restart == 'y') then
        if (nprocs /= nproc_size_global) then
           if(comm_is_root(nproc_id_global)) &
           write(*,*) 'Warning: number of processes do not match!'
          !stop 'number of processes do not match!'
        end if
     end if
  else
     mk = system%nk
     mo = system%no
  endif

  !occupation
  if(flag_read_occ) then
  if(comm_is_root(nproc_id_global))then
     dir_file_in = trim(idir)//"occupation.bin"
     open(iu1_r,file=dir_file_in,form='unformatted')

     allocate(roccbox(mo,mk,system%nspin))
     read(iu1_r) roccbox(1:mo,1:mk,1:system%nspin)
     system%rocc(1:system%no,1:system%nk,1:system%nspin) = &
     roccbox    (1:system%no,1:system%nk,1:system%nspin)
     deallocate(roccbox)

     close(iu1_r)
  end if
  call comm_bcast(system%rocc,comm)
  endif

  !wave function
  if( flag_GS .and. &
      (read_gs_restart_data(1:9)=='rho_inout'.or.read_gs_restart_data=='rho'))then
     ! this case => do not read wavefunction
  else
     call read_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,if_real_orbital,iself)
  endif
  !rho_inout
  if( flag_GS.and. &
     (read_gs_restart_data(1:3)=='all'.or.read_gs_restart_data(1:9)=='rho_inout'))then
    if (present(mixing)) then
      call read_rho_inout(idir,lg,mg,system,info,mixing,iself)
    end if
  end if

  !rho (only for GS)
  if( flag_GS .and. read_gs_restart_data=='rho' )then
     call read_rho(idir,lg,mg,system,info,mixing)
  endif

  !Vh_stock
 !if( flag_RT .and. yn_restart=='y')then
  if( flag_RT .and. yn_restart=='y' .and. yn_reset_step_restart=='n' )then
    if (present(Vh_stock1) .and. present(Vh_stock2)) then
      call read_Vh_stock(idir,lg,mg,info,Vh_stock1,Vh_stock2,iself)
    end if
  end if

end subroutine read_bin


!===================================================================================================================================

subroutine write_wavefunction(odir,lg,mg,system,info,spsi,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital
  use communication, only: comm_is_root, comm_summation, comm_bcast
#ifdef USE_MPI
#else
  use salmon_global, only: method_wf_distributor
#endif
  implicit none
  character(*),            intent(in) :: odir
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital)                     :: spsi
  logical,                 intent(in) :: is_self_checkpoint

  integer :: iu2_w
  character(256) :: dir_file_out
  logical :: is_written

#ifdef USE_MPI
#else
  integer :: im,ik,io,is
#endif

  iu2_w = 87

  if(is_self_checkpoint) then
    ! write all processes (each process dump data)
    dir_file_out = trim(odir)//"wfn.bin"
    open(iu2_w,file=dir_file_out,form='unformatted',access='stream')

    if(allocated(spsi%rwf))then
      write (iu2_w) spsi%rwf(mg%is(1):mg%ie(1),   &
                             mg%is(2):mg%ie(2),   &
                             mg%is(3):mg%ie(3),   &
                             1:system%nspin,      &
                             info%io_s:info%io_e, &
                             info%ik_s:info%ik_e, &
                             info%im_s:info%im_e)
    else if(allocated(spsi%zwf))then
      write (iu2_w) spsi%zwf(mg%is(1):mg%ie(1),   &
                             mg%is(2):mg%ie(2),   &
                             mg%is(3):mg%ie(3),   &
                             1:system%nspin,      &
                             info%io_s:info%io_e, &
                             info%ik_s:info%ik_e, &
                             info%im_s:info%im_e)
    end if
  else
#ifdef USE_MPI
    call distributed_rw_wavefunction(odir,lg,mg,system,info,spsi,system%nk,system%no,system%if_real_orbital,write_mode)
#else
    ! single process execution
    select case(method_wf_distributor)
    case('single') ; call write_all    ! create single shared file
    case('slice')  ; call write_sliced ! sliced shared file
    case default   ; stop 'write_wavefunction: fatal error'
    end select
#endif
  end if

  !close file iu2_w
  inquire(iu2_w, opened=is_written)
  if (is_written) close(iu2_w)

#ifdef USE_MPI
#else
contains
  subroutine write_all
    implicit none
    dir_file_out = trim(odir)//"wfn.bin"
    open(iu2_w,file=dir_file_out,form='unformatted',access='stream')

    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
    do is=1,system%nspin
      if(allocated(spsi%rwf))then
        write (iu2_w) spsi%rwf(lg%is(1):lg%ie(1),   &
                               lg%is(2):lg%ie(2),   &
                               lg%is(3):lg%ie(3),   &
                               is,io,ik,im)
      else if(allocated(spsi%zwf))then
        write (iu2_w) spsi%zwf(lg%is(1):lg%ie(1),   &
                               lg%is(2):lg%ie(2),   &
                               lg%is(3):lg%ie(3),   &
                               is,io,ik,im)
      end if
    end do
    end do
    end do
    end do
  end subroutine

  subroutine write_sliced
    use salmon_global, only: nblock_wf_distribute
    use filesystem
    implicit none
    integer :: nblock_orbital
    integer :: ik,io,nb
    integer :: mo

    mo = info%io_e - info%io_s + 1
    nblock_orbital = min(mo,nblock_wf_distribute)

    ! create all directory
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e,nblock_orbital
      write (dir_file_out,'(A,I6.6,A,I6.6)') trim(odir)//'k_',ik,'_ob_',io
      call create_directory(dir_file_out)
    end do
    end do

    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
      nb = ((io - 1) / nblock_orbital) * nblock_orbital + 1
      write (dir_file_out,'(A,I6.6,A,I6.6,A,I6.6,A)') trim(odir)//'k_',ik,'_ob_',nb,'/wfn_ob_',io,'.dat'
      open(iu2_w,file=dir_file_out,form='unformatted',access='stream')

      do is=1,system%nspin
        if(allocated(spsi%rwf))then
          write (iu2_w) spsi%rwf(lg%is(1):lg%ie(1),   &
                                 lg%is(2):lg%ie(2),   &
                                 lg%is(3):lg%ie(3),   &
                                 is,io,ik,1)
        else if(allocated(spsi%zwf))then
          write (iu2_w) spsi%zwf(lg%is(1):lg%ie(1),   &
                                 lg%is(2):lg%ie(2),   &
                                 lg%is(3):lg%ie(3),   &
                                 is,io,ik,1)
        end if
      end do

      close(iu2_w)
    end do
    end do
  end subroutine
#endif
end subroutine write_wavefunction

!===================================================================================================================================

subroutine write_rho_inout(odir,lg,mg,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)                        :: odir
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_mixing)          ,intent(in) :: mixing
  logical                 ,intent(in) :: is_self_checkpoint
  !
  character(100) ::  dir_file_out
  integer :: i,ix,iy,iz,is
  integer :: iu1_w
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)

  iu1_w = 97
  dir_file_out = trim(odir)//"rho_inout.bin"

  if(is_self_checkpoint) then
    ! write all processes
    open(iu1_w,file=dir_file_out,form='unformatted',access='stream')
    do i=1,mixing%num_rho_stock+1
      write(iu1_w) mixing%rho_in (i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    end do
    do i=1,mixing%num_rho_stock
      write(iu1_w) mixing%rho_out(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    end do
    if(system%nspin==2) then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
          write(iu1_w) mixing%rho_s_in (i,is)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        end do
        do i=1,mixing%num_rho_stock
          write(iu1_w) mixing%rho_s_out(i,is)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        end do
      end do
    end if
    close(iu1_w)
  else
    ! write root process
    if(comm_is_root(nproc_id_global))then
      open(iu1_w,file=dir_file_out,form='unformatted')
    end if

    allocate(matbox (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

#ifdef USE_OPENACC
!$acc parallel loop collapse(2)
#else
!$omp parallel do collapse(2) private(iz,iy,ix)
#endif
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      matbox2(ix,iy,iz)=0.d0
    end do
    end do
    end do

    do i=1,mixing%num_rho_stock+1
#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix)
#else
!$omp parallel do collapse(2) private(iz,iy,ix)
#endif
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        matbox2(ix,iy,iz) = mixing%rho_in(i)%f(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
      if(comm_is_root(nproc_id_global))then
        write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
    end do

    do i=1,mixing%num_rho_stock
#ifdef USE_OPENACC
!$acc parallel loop private(iz,iy,ix)
#else
!$omp parallel do collapse(2) private(iz,iy,ix)
#endif
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        matbox2(ix,iy,iz) = mixing%rho_out(i)%f(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
      if(comm_is_root(nproc_id_global))then
        write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
    end do

    if(system%nspin == 2)then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            matbox2(ix,iy,iz) = mixing%rho_s_in(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
          call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
          if(comm_is_root(nproc_id_global))then
            write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
          end if
        end do

        do i=1,mixing%num_rho_stock
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            matbox2(ix,iy,iz) = mixing%rho_s_out(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
          call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
          if(comm_is_root(nproc_id_global))then
            write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
          end if
        end do
      end do
    end if

    if(comm_is_root(nproc_id_global))then
      close(iu1_w)
    end if

    deallocate(matbox,matbox2)
  end if
end subroutine write_rho_inout

subroutine write_rho(odir,lg,mg,system,info,mixing)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)                        :: odir
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_dft_system)      ,intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_mixing)          ,intent(in) :: mixing
  !
  character(100) ::  dir_file_out
  integer :: i,ix,iy,iz,is
  integer :: iu1_w
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)

  iu1_w = 97
  dir_file_out = trim(odir)//"rho.bin"

  ! write root process
  if(comm_is_root(nproc_id_global)) open(iu1_w,file=dir_file_out,form='unformatted')

  allocate(matbox (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

  !$omp parallel do collapse(2) private(iz,iy,ix)
  do iz=lg%is(3),lg%ie(3)
  do iy=lg%is(2),lg%ie(2)
  do ix=lg%is(1),lg%ie(1)
     matbox2(ix,iy,iz)=0.d0
  end do
  end do
  end do

  i=mixing%num_rho_stock+1
  !$omp parallel do collapse(2) private(iz,iy,ix)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
     matbox2(ix,iy,iz) = mixing%rho_in(i)%f(ix,iy,iz)
  end do
  end do
  end do
  call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
  if(comm_is_root(nproc_id_global))then
     write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  end if

  if(system%nspin == 2)then
     do is=1,2
     i=mixing%num_rho_stock+1
     !$omp parallel do collapse(2) private(iz,iy,ix)
     do iz=mg%is(3),mg%ie(3)
     do iy=mg%is(2),mg%ie(2)
     do ix=mg%is(1),mg%ie(1)
        matbox2(ix,iy,iz) = mixing%rho_s_in(i,is)%f(ix,iy,iz)
     end do
     end do
     end do
     call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)
     if(comm_is_root(nproc_id_global))then
        write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
     end if
     end do
  endif

  if(comm_is_root(nproc_id_global)) close(iu1_w)

  deallocate(matbox,matbox2)

end subroutine write_rho

!===================================================================================================================================

subroutine write_Vh_stock(odir,lg,mg,info,Vh_stock1,Vh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_parallel_info, s_scalar
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)   :: odir
  type(s_rgrid), intent(in)    :: lg,mg
  type(s_parallel_info),intent(in) :: info
  type(s_scalar),intent(in) :: Vh_stock1,Vh_stock2
  logical,intent(in) :: is_self_checkpoint

  character(256) ::  dir_file_out
  integer :: iu1_w
  real(8),allocatable :: matbox0(:,:,:),matbox1(:,:,:),matbox2(:,:,:)
  integer :: ix,iy,iz

  iu1_w = 97
  dir_file_out = trim(odir)//"Vh_stock.bin"

  if (is_self_checkpoint) then
    ! write all processes
    open(iu1_w,file=dir_file_out,form='unformatted',access='stream')
    write(iu1_w) Vh_stock1%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) Vh_stock2%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    close(iu1_w)
  else
    ! write root process
    allocate(matbox0(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      matbox0(ix,iy,iz) = 0d0
    end do
    end do
    end do

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      matbox0(ix,iy,iz) = Vh_stock1%f(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(matbox0,matbox1,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      matbox0(ix,iy,iz) = Vh_stock2%f(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(matbox0,matbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)

    if(comm_is_root(nproc_id_global))then
      open(iu1_w,file=dir_file_out,form='unformatted')
      write(iu1_w) matbox1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      write(iu1_w) matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      close(iu1_w)
    end if

    deallocate(matbox0,matbox1,matbox2)
  end if

end subroutine write_Vh_stock

!===================================================================================================================================

subroutine write_singlescale(odir,lg,mg,info,singlescale,Ac,div_Ac,is_self_checkpoint)
  use structures
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: method_singlescale
  implicit none
  character(*)            ,intent(in) :: odir
  type(s_rgrid)           ,intent(in) :: lg,mg
  type(s_parallel_info)   ,intent(in) :: info
  type(s_singlescale)     ,intent(in) :: singlescale
  type(s_vector)          ,intent(in) :: Ac
  type(s_scalar)          ,intent(in) :: div_Ac
  logical                 ,intent(in) :: is_self_checkpoint
  !
  character(256) ::  dir_file_out
  integer :: iu1_w
  real(8),allocatable :: matbox0(:,:,:,:,:),matbox1(:,:,:,:,:)
  real(8),allocatable :: v0(:,:,:,:,:),v1(:,:,:,:,:)
  real(8),allocatable :: b0(:,:,:,:),b1(:,:,:,:)
  real(8),allocatable :: d0(:,:,:,:),d1(:,:,:,:)
  complex(8),allocatable :: z0(:,:,:,:,:),z1(:,:,:,:,:)
  integer :: ix,iy,iz

  iu1_w = 97
  dir_file_out = trim(odir)//"singlescale.bin"

  if (is_self_checkpoint) then
    ! write all processes
    open(iu1_w,file=dir_file_out,form='unformatted',access='stream')
    write(iu1_w) singlescale%vec_Ac_m(-1:1,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:3)
    write(iu1_w) Ac%v(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) singlescale%vec_je_old(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) singlescale%vec_Ac_old(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) singlescale%vec_Ac_boundary_bottom_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    write(iu1_w) singlescale%vec_Ac_boundary_top_old   (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    write(iu1_w) singlescale%vec_Ac_boundary_bottom    (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    write(iu1_w) singlescale%vec_Ac_boundary_top       (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    write(iu1_w) div_Ac%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) singlescale%div_Ac_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    write(iu1_w) singlescale%Energy_joule
    write(iu1_w) singlescale%Energy_poynting
    if(method_singlescale/='3d') then
      write(iu1_w) singlescale%Ac_zt_m(lg%is(3)-1:lg%ie(3)+1,-1:1,1:3)
      write(iu1_w) singlescale%zf_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      write(iu1_w) singlescale%zc_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      write(iu1_w) singlescale%zs_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      write(iu1_w) singlescale%Ac_zt_boundary_bottom
      write(iu1_w) singlescale%Ac_zt_boundary_top
      write(iu1_w) singlescale%Ac_zt_boundary_bottom_old
      write(iu1_w) singlescale%Ac_zt_boundary_top_old
    end if
    close(iu1_w)
  else
    ! write root process
    allocate(matbox0(-1:1,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))
    allocate(matbox1(-1:1,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))

!$omp workshare
    matbox0 = 0d0
!$omp end workshare

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      matbox0(-1:1,ix,iy,iz,1:3) = singlescale%vec_Ac_m(-1:1,ix,iy,iz,1:3)
    end do
    end do
    end do
    call comm_summation(matbox0,matbox1,9*lg%num(1)*lg%num(2)*lg%num(3),info%icomm_r)

    allocate(v0(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))
    allocate(v1(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))

!$omp workshare
    v0 = 0d0
!$omp end workshare

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      v0(1:3,ix,iy,iz,1) = Ac%v(1:3,ix,iy,iz)
      v0(1:3,ix,iy,iz,2) = singlescale%vec_je_old(1:3,ix,iy,iz)
      v0(1:3,ix,iy,iz,3) = singlescale%vec_Ac_old(1:3,ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(v0,v1,3*lg%num(1)*lg%num(2)*lg%num(3)*3,info%icomm_r)
    
    allocate(b0(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),1:3,1:4))
    allocate(b1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),1:3,1:4))

!$omp workshare
    b0 = 0d0
!$omp end workshare

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      b0(ix,iy,1:3,1) = singlescale%vec_Ac_boundary_bottom_old(ix,iy,1:3)
      b0(ix,iy,1:3,2) = singlescale%vec_Ac_boundary_top_old   (ix,iy,1:3)
      b0(ix,iy,1:3,3) = singlescale%vec_Ac_boundary_bottom    (ix,iy,1:3)
      b0(ix,iy,1:3,4) = singlescale%vec_Ac_boundary_top       (ix,iy,1:3)
    end do
    end do
    call comm_summation(b0,b1,3*lg%num(1)*lg%num(2)*4,info%icomm_r)
    
    allocate(d0(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:2))
    allocate(d1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:2))

!$omp workshare
    d0 = 0d0
!$omp end workshare

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      d0(ix,iy,iz,1) = div_Ac%f(ix,iy,iz)
      d0(ix,iy,iz,2) = singlescale%div_Ac_old(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(d0,d1,lg%num(1)*lg%num(2)*lg%num(3)*2,info%icomm_r)
    
    if(method_singlescale/='3d') then
    
      allocate(z0(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,1:3))
      allocate(z1(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,1:3))

  !$omp workshare
      z0 = 0d0
  !$omp end workshare

  !$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        z0(ix,iy,iz,0:3,1) = singlescale%zf_old(ix,iy,iz,0:3)
        z0(ix,iy,iz,0:3,2) = singlescale%zc_old(ix,iy,iz,0:3)
        z0(ix,iy,iz,0:3,3) = singlescale%zs_old(ix,iy,iz,0:3)
      end do
      end do
      end do
      call comm_summation(z0,z1,lg%num(1)*lg%num(2)*lg%num(3)*4*3,info%icomm_r)
      
    end if

    if(comm_is_root(nproc_id_global))then
      open(iu1_w,file=dir_file_out,form='unformatted')
      write(iu1_w) matbox1(-1:1,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3)
      write(iu1_w) v1(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3)
      write(iu1_w) b1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),1:3,1:4)
      write(iu1_w) d1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:2)
      write(iu1_w) singlescale%Energy_joule
      write(iu1_w) singlescale%Energy_poynting
      if(method_singlescale/='3d') then
        write(iu1_w) singlescale%Ac_zt_m(lg%is(3)-1:lg%ie(3)+1,-1:1,1:3)
        write(iu1_w) z1(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,1)
        write(iu1_w) z1(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,2)
        write(iu1_w) z1(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,3)
        write(iu1_w) singlescale%Ac_zt_boundary_bottom
        write(iu1_w) singlescale%Ac_zt_boundary_top
        write(iu1_w) singlescale%Ac_zt_boundary_bottom_old
        write(iu1_w) singlescale%Ac_zt_boundary_top_old
      end if
      close(iu1_w)
    end if

    deallocate(matbox0,matbox1,v0,v1,b0,b1,d0,d1)
    if(method_singlescale/='3d') deallocate(z0,z1)
  end if

end subroutine write_singlescale

!===================================================================================================================================

subroutine read_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,if_real_orbital,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital, &
  &                     allocate_orbital_real, deallocate_orbital
  use communication, only: comm_is_root, comm_summation, comm_bcast
#ifdef USE_MPI
#else
  use salmon_global, only: method_wf_distributor
#endif
  implicit none
  character(*),            intent(in) :: idir
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_orbital)                     :: spsi
  integer,                 intent(in) :: mk, mo
  logical,                 intent(in) :: if_real_orbital,is_self_checkpoint

  integer :: iu2_r
  character(256) :: dir_file_in
  logical :: is_read
#ifdef USE_MPI
#else
  real(8),allocatable :: ddummy(:,:,:)
  complex(8),allocatable :: zdummy(:,:,:)
#endif

  iu2_r = 86

  if(is_self_checkpoint) then
    ! read all processes (each process load dumped data)
    dir_file_in = trim(idir)//"wfn.bin"
    open(iu2_r,file=dir_file_in,form='unformatted',access='stream')

    if (allocated(spsi%rwf)) then
      read (iu2_r) spsi%rwf(mg%is(1):mg%ie(1),   &
                            mg%is(2):mg%ie(2),   &
                            mg%is(3):mg%ie(3),   &
                            1:system%nspin,      &
                            info%io_s:info%io_e, &
                            info%ik_s:info%ik_e, &
                            info%im_s:info%im_e)
    else if (allocated(spsi%zwf)) then
      read (iu2_r) spsi%zwf(mg%is(1):mg%ie(1),   &
                            mg%is(2):mg%ie(2),   &
                            mg%is(3):mg%ie(3),   &
                            1:system%nspin,      &
                            info%io_s:info%io_e, &
                            info%ik_s:info%ik_e, &
                            info%im_s:info%im_e)
    end if
  else
#ifdef USE_MPI
    call distributed_rw_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,if_real_orbital,read_mode)
#else
    ! single process execution
    select case(method_wf_distributor)
    case('single') ; call read_all    ! create single shared file
    case('slice')  ; call read_sliced ! sliced shared file
    case default   ; stop 'read_wavefunction: fatal error'
    end select
#endif
  end if

  !close file iu2_r
  inquire(iu2_r, opened=is_read)
  if (is_read) close(iu2_r)

#ifdef USE_MPI
#else
contains
  subroutine read_all
    implicit none
    integer :: im,ik,io,is

    dir_file_in = trim(idir)//"wfn.bin"
    open(iu2_r,file=dir_file_in,form='unformatted',access='stream')

    if (if_real_orbital) then
      allocate(ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    else
      allocate(zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    end if

    do im=info%im_s,info%im_e
    do ik=1,mk
    do io=1,mo
    do is=1,system%nspin
      if (if_real_orbital) then
        read (iu2_r) ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      else
        read (iu2_r) zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      endif

      if (info%ik_s <= ik .and. ik <= info%ik_e .and. &
          info%io_s <= io .and. io <= info%io_e) then
        if (if_real_orbital) then
          if (allocated(spsi%rwf)) then
            spsi%rwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,im) = ddummy(:,:,:)
          else
            spsi%zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,im) = dcmplx(ddummy(:,:,:))
          end if
        else
          spsi%zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,im) = zdummy(:,:,:)
        end if
      end if
    end do
    end do
    end do
    end do

    if (allocated(ddummy)) deallocate(ddummy)
    if (allocated(zdummy)) deallocate(zdummy)
  end subroutine

  subroutine read_sliced
    use salmon_global, only: nblock_wf_distribute
    use filesystem
    implicit none
    integer :: nblock_orbital
    integer :: ik,io,nb,io_e,is

    nblock_orbital = min(mo,nblock_wf_distribute)

    if (if_real_orbital) then
      allocate(ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    else
      allocate(zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    end if

    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
      nb = ((io - 1) / nblock_orbital) * nblock_orbital + 1
      write (dir_file_in,'(A,I6.6,A,I6.6,A,I6.6,A)') trim(idir)//'k_',ik,'_ob_',nb,'/wfn_ob_',io,'.dat'
      open(iu2_r,file=dir_file_in,form='unformatted',access='stream')

      do is=1,system%nspin
        if (if_real_orbital) then
          read (iu2_r) ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
        else
          read (iu2_r) zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
        endif

        if (if_real_orbital) then
          if (allocated(spsi%rwf)) then
            spsi%rwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,1) = ddummy(:,:,:)
          else
            spsi%zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,1) = dcmplx(ddummy(:,:,:))
          end if
        else
          spsi%zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,1) = zdummy(:,:,:)
        end if
      end do

      close(iu2_r)
    end do
    end do
  end subroutine
#endif
end subroutine read_wavefunction

!===================================================================================================================================

subroutine read_rho_inout(idir,lg,mg,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*), intent(in) :: idir
  type(s_rgrid), intent(in)    :: lg,mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_mixing),intent(inout) :: mixing
  logical,intent(in) :: is_self_checkpoint

  integer :: iu1_r
  integer :: i,ix,iy,iz,is
  real(8),allocatable :: matbox(:,:,:)
  character(100) :: dir_file_in

  iu1_r = 96
  dir_file_in = trim(idir)//"rho_inout.bin"

  if(is_self_checkpoint) then
    ! read all processses
    open(iu1_r,file=dir_file_in,form='unformatted',access='stream')
    do i=1,mixing%num_rho_stock+1
      read(iu1_r) mixing%rho_in (i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    end do
    do i=1,mixing%num_rho_stock
      read(iu1_r) mixing%rho_out(i)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    end do
    if(system%nspin==2) then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
          read(iu1_r) mixing%rho_s_in (i,is)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        end do
        do i=1,mixing%num_rho_stock
          read(iu1_r) mixing%rho_s_out(i,is)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
        end do
      end do
    end if
    close(iu1_r)
  else
    ! read root process
    if(comm_is_root(nproc_id_global))then
      open(iu1_r,file=dir_file_in,form='unformatted')
    end if

    allocate(matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

    do i=1,mixing%num_rho_stock+1
      if(comm_is_root(nproc_id_global))then
        read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
      call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_in(i)%f(ix,iy,iz)=matbox(ix,iy,iz)
      end do
      end do
      end do
    end do

    do i=1,mixing%num_rho_stock
      if(comm_is_root(nproc_id_global))then
        read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
      call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        mixing%rho_out(i)%f(ix,iy,iz)=matbox(ix,iy,iz)
      end do
      end do
      end do
    end do

    if(system%nspin == 2)then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
          if(comm_is_root(nproc_id_global))then
            read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
          end if
          call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            mixing%rho_s_in(i,is)%f(ix,iy,iz)=matbox(ix,iy,iz)
          end do
          end do
          end do
        end do

        do i=1,mixing%num_rho_stock
          if(comm_is_root(nproc_id_global))then
            read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
          end if
          call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            mixing%rho_s_out(i,is)%f(ix,iy,iz)=matbox(ix,iy,iz)
          end do
          end do
          end do
        end do
      end do
    end if

    if(comm_is_root(nproc_id_global))then
      close(iu1_r)
    end if

    deallocate(matbox)
  end if

end subroutine read_rho_inout

subroutine read_rho(idir,lg,mg,system,info,mixing)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*), intent(in) :: idir
  type(s_rgrid), intent(in)    :: lg,mg
  type(s_dft_system),intent(in) :: system
  type(s_parallel_info),intent(in) :: info
  type(s_mixing),intent(inout) :: mixing

  integer :: iu1_r
  integer :: i,ix,iy,iz,is
  real(8),allocatable :: matbox(:,:,:)
  character(100) :: dir_file_in

  iu1_r = 96
  dir_file_in = trim(idir)//"rho.bin"

  ! read root process
  if(comm_is_root(nproc_id_global)) open(iu1_r,file=dir_file_in,form='unformatted')

  allocate(matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

  i = mixing%num_rho_stock+1
  if(comm_is_root(nproc_id_global))then
     read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
  end if
  call comm_bcast(matbox,info%icomm_rko)

  !$omp parallel do collapse(2)
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
     mixing%rho_in(i)%f(ix,iy,iz)=matbox(ix,iy,iz)
  end do
  end do
  end do

  if(system%nspin == 2)then
     do is=1,2
        i = mixing%num_rho_stock+1
        if(comm_is_root(nproc_id_global))then
           read(iu1_r) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
        end if
        call comm_bcast(matbox,info%icomm_rko)

        !$omp parallel do collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
           mixing%rho_s_in(i,is)%f(ix,iy,iz)=matbox(ix,iy,iz)
        end do
        end do
        end do

     end do
  end if

  if(comm_is_root(nproc_id_global)) close(iu1_r)

  deallocate(matbox)

end subroutine read_rho

!===================================================================================================================================

subroutine read_Vh_stock(idir,lg,mg,info,Vh_stock1,Vh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_parallel_info, s_scalar
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*), intent(in) :: idir
  type(s_rgrid), intent(in)    :: lg,mg
  type(s_parallel_info),intent(in) :: info
  type(s_scalar),intent(inout) :: Vh_stock1,Vh_stock2
  logical,intent(in) :: is_self_checkpoint

  integer :: iu1_r
  integer :: ix,iy,iz
  real(8),allocatable :: matbox1(:,:,:),matbox2(:,:,:)
  character(100) :: dir_file_in

  iu1_r = 96
  dir_file_in = trim(idir)//"Vh_stock.bin"

  if (is_self_checkpoint) then
    ! read all processes
    open(iu1_r,file=dir_file_in,form='unformatted',access='stream')
    read(iu1_r) Vh_stock1%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) Vh_stock2%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    close(iu1_r)
  else
    ! read root process
    allocate(matbox1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

    if(comm_is_root(nproc_id_global))then
      open(iu1_r,file=dir_file_in,form='unformatted')
      read(iu1_r) matbox1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      read(iu1_r) matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      close(iu1_r)
    end if

    call comm_bcast(matbox1,info%icomm_rko)
    call comm_bcast(matbox2,info%icomm_rko)

!$omp parallel do collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vh_stock1%f(ix,iy,iz)=matbox1(ix,iy,iz)
    end do
    end do
    end do

!$omp parallel do collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      Vh_stock2%f(ix,iy,iz)=matbox2(ix,iy,iz)
    end do
    end do
    end do

    deallocate(matbox1,matbox2)
  end if
end subroutine read_Vh_stock

!===================================================================================================================================

subroutine restart_singlescale(comm,lg,mg,singlescale,Ac,div_Ac)
  use structures
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: directory_read_data,yn_self_checkpoint,method_singlescale
  implicit none
  integer      ,intent(in) :: comm
  type(s_rgrid),intent(in) :: lg,mg
  type(s_singlescale)      :: singlescale
  type(s_vector)           :: Ac
  type(s_scalar)           :: div_Ac
  !
  integer :: iu1_r, ierr
  integer :: ix,iy,iz
  real(8),allocatable :: matbox1(:,:,:,:,:),matbox2(:,:,:,:,:),matbox3(:,:,:,:),matbox4(:,:,:,:)
  complex(8),allocatable :: zbox(:,:,:,:,:)
  character(100) :: dir_file_in
  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_self_checkpoint == 'y')
  if (.not. iself) then
    wdir = gdir
  end if

  if (comm_is_root(nproc_id_global)) then
     write(*,*) "The data for the single-scale Maxwell-TDDFT method is read from restart directory"
  endif
  iu1_r = 96
  dir_file_in = trim(wdir)//"singlescale.bin"

  if (iself) then
    ! read all processes
    open(iu1_r,file=dir_file_in,form='unformatted',access='stream',status="old",err=20)
    read(iu1_r) singlescale%vec_Ac_m(-1:1,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1:3)
    read(iu1_r) Ac%v(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) singlescale%vec_je_old(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) singlescale%vec_Ac_old(1:3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) singlescale%vec_Ac_boundary_bottom_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    read(iu1_r) singlescale%vec_Ac_boundary_top_old   (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    read(iu1_r) singlescale%vec_Ac_boundary_bottom    (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    read(iu1_r) singlescale%vec_Ac_boundary_top       (mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),1:3)
    read(iu1_r) div_Ac%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) singlescale%div_Ac_old(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
    read(iu1_r) singlescale%Energy_joule
    read(iu1_r) singlescale%Energy_poynting
    if(method_singlescale/='3d') then
      read(iu1_r) singlescale%Ac_zt_m(lg%is(3)-1:lg%ie(3)+1,-1:1,1:3)
      read(iu1_r) singlescale%zf_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      read(iu1_r) singlescale%zc_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      read(iu1_r) singlescale%zs_old(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3)
      read(iu1_r) singlescale%Ac_zt_boundary_bottom
      read(iu1_r) singlescale%Ac_zt_boundary_top
      read(iu1_r) singlescale%Ac_zt_boundary_bottom_old
      read(iu1_r) singlescale%Ac_zt_boundary_top_old
    end if
    close(iu1_r)
  else
    ! read root process
    allocate(matbox1(-1:1,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))
    allocate(matbox2(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3))
    allocate(matbox3(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),1:3,1:4))
    allocate(matbox4(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:2))
    if(method_singlescale/='3d') allocate(zbox(1:lg%num(1),1:mg%num(2),1:mg%num(3),0:3,1:3))

    if(comm_is_root(nproc_id_global))then
      open(iu1_r,file=dir_file_in,form='unformatted',status='old',iostat=ierr)
      if(ierr==0) then
      read(iu1_r) matbox1(-1:1,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3)
      read(iu1_r) matbox2(1:3,lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:3)
      read(iu1_r) matbox3(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),1:3,1:4)
      read(iu1_r) matbox4(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),1:2)
      read(iu1_r) singlescale%Energy_joule
      read(iu1_r) singlescale%Energy_poynting
      if(method_singlescale/='3d') then
        read(iu1_r) singlescale%Ac_zt_m(lg%is(3)-1:lg%ie(3)+1,-1:1,1:3)
        read(iu1_r) zbox(1:lg%num(1),1:lg%num(2),1:lg%num(3),0:3,1:3)
        read(iu1_r) singlescale%Ac_zt_boundary_bottom
        read(iu1_r) singlescale%Ac_zt_boundary_top
        read(iu1_r) singlescale%Ac_zt_boundary_bottom_old
        read(iu1_r) singlescale%Ac_zt_boundary_top_old
      end if
      close(iu1_r)
      endif
    end if

    call comm_bcast(ierr,comm)
    if(ierr.ne.0) goto 20

    call comm_bcast(matbox1,comm)
    call comm_bcast(matbox2,comm)
    call comm_bcast(matbox3,comm)
    call comm_bcast(matbox4,comm)
    call comm_bcast(singlescale%Energy_joule,comm)
    call comm_bcast(singlescale%Energy_poynting,comm)
    if(method_singlescale/='3d') then
      call comm_bcast(singlescale%Ac_zt_m,comm)
      call comm_bcast(zbox,comm)
      call comm_bcast(singlescale%Ac_zt_boundary_bottom,comm)
      call comm_bcast(singlescale%Ac_zt_boundary_top,comm)
      call comm_bcast(singlescale%Ac_zt_boundary_bottom_old,comm)
      call comm_bcast(singlescale%Ac_zt_boundary_top_old,comm)
    end if

!$omp parallel do collapse(2)
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      singlescale%vec_Ac_m(-1:1,ix,iy,iz,1:3) = matbox1(-1:1,ix,iy,iz,1:3)
      Ac%v(1:3,ix,iy,iz)                   = matbox2(1:3,ix,iy,iz,1)
      singlescale%vec_je_old(1:3,ix,iy,iz) = matbox2(1:3,ix,iy,iz,2)
      singlescale%vec_Ac_old(1:3,ix,iy,iz) = matbox2(1:3,ix,iy,iz,3)
      div_Ac%f(ix,iy,iz)               = matbox4(ix,iy,iz,1)
      singlescale%div_Ac_old(ix,iy,iz) = matbox4(ix,iy,iz,2)
    end do
    end do
    end do
    
    !$omp parallel do collapse(2)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      singlescale%vec_Ac_boundary_bottom_old(ix,iy,1:3) = matbox3(ix,iy,1:3,1)
      singlescale%vec_Ac_boundary_top_old   (ix,iy,1:3) = matbox3(ix,iy,1:3,2)
      singlescale%vec_Ac_boundary_bottom    (ix,iy,1:3) = matbox3(ix,iy,1:3,3)
      singlescale%vec_Ac_boundary_top       (ix,iy,1:3) = matbox3(ix,iy,1:3,4)
    end do
    end do

    if(method_singlescale/='3d') then
      !$omp parallel do collapse(2)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=1,lg%num(1)
        singlescale%zf_old(ix,iy,iz,0:3) = zbox(ix,iy,iz,0:3,1)
        singlescale%zc_old(ix,iy,iz,0:3) = zbox(ix,iy,iz,0:3,2)
        singlescale%zs_old(ix,iy,iz,0:3) = zbox(ix,iy,iz,0:3,3)
      end do
      end do
      end do
    end if

    deallocate(matbox1,matbox2,matbox3,matbox4)
    if(method_singlescale/='3d') deallocate(zbox)
  end if

  return

20 continue  !restart from e.x. theory=tddft_pulse or dft
  if (comm_is_root(nproc_id_global)) then
     write(*,*) "no singlescale.bin was found: restart from other theory is assumed"
  endif
  return
end subroutine restart_singlescale

!===================================================================================================================================
!(currently not used: see subroutine "read_atomic_coordinates")
subroutine restart_Rion(system)
  use structures, only: s_dft_system
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint
  implicit none
  type(s_dft_system)     ,intent(inout) :: system
  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (.not. iself) then
    wdir = gdir
  end if

  call read_Rion(wdir,system)
end subroutine restart_Rion

subroutine restart_Velocity(system)
  use structures, only: s_dft_system
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint
  implicit none
  type(s_dft_system), intent(inout) :: system
  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (.not. iself) then
    wdir = gdir
  end if

  if (comm_is_root(nproc_id_global)) then
     write(*,*) "  Initial velocities is read from restart directory"
  endif


  call read_Velocity(wdir,system)
end subroutine restart_Velocity


subroutine write_Rion(odir,system)
  use structures, only: s_dft_system
  use salmon_global, only: natom, atom_name, kion, unit_length,yn_opt,flag_opt_atom
!  use inputoutput, only: au_length_aa   !?? why error??
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_dft_system),intent(in) :: system
  integer :: iu1_w, ia
  real(8) :: uconv
  character(*)   :: odir
  character(256) :: dir_file_out

  real(8),parameter :: au_length_aa = 0.52917721067d0

  iu1_w = 87

  ! atomic coordinate
  if(unit_length=='AA')then ; uconv = au_length_aa
  else                      ; uconv = 1d0   !au
  endif

  if(comm_is_root(nproc_id_global)) then
     dir_file_out = trim(odir)//"atomic_coor.txt"
     open(iu1_w, file=dir_file_out, status="unknown")
     if(yn_opt == 'y')then
        do ia = 1,natom
           write(iu1_w,7000) trim(atom_name(ia)), system%Rion(1:3,ia)*uconv, kion(ia), flag_opt_atom(ia)
        enddo
     else
        do ia = 1,natom
           write(iu1_w,7100) trim(atom_name(ia)), system%Rion(1:3,ia)*uconv, kion(ia)
        enddo
     endif
     close(iu1_w)
  end if

7000 format("'",a,"'  ",3f24.18,i4,'   ',a)
7100 format("'",a,"'  ",3f24.18,i4)
9000 format(a)

end subroutine write_Rion

subroutine write_Velocity(odir,system)
  use structures, only: s_dft_system
  use salmon_global, only: natom
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_dft_system),intent(in) :: system
  integer :: iu1_w, ia
  character(*)   :: odir
  character(256) :: dir_file_out

  iu1_w = 87

  ! atomic velocity [au]
  if(comm_is_root(nproc_id_global)) then
     dir_file_out = trim(odir)//"atomic_vel.txt"
     open(iu1_w, file=dir_file_out, status="unknown")
     do ia = 1,natom
        write(iu1_w,8000) system%Velocity(1:3,ia)
     enddo
     close(iu1_w)
  end if

8000 format(3f24.14)
9000 format(a)

end subroutine write_Velocity

!(currently not used: see subroutine "read_atomic_coordinates")
subroutine read_Rion(idir,system)
  use structures, only: s_dft_system
  use salmon_global, only: natom, atom_name, kion, unit_length !, Rion
!  use inputoutput, only: au_length_aa  !?? why error??
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_dft_system) :: system
  character(*),intent(in) :: idir
  real(8) :: uconv
  integer :: iu1_w, ia, comm
  character(256) :: dir_file_out

  real(8),parameter :: au_length_aa = 0.52917721067d0

  iu1_w = 87
  comm = nproc_group_global

  ! atomic coordinate
  if(unit_length=='AA')then ; uconv = au_length_aa
  else                      ; uconv = 1d0   !au
  endif

  if(comm_is_root(nproc_id_global)) then
     dir_file_out = trim(idir)//"atomic_coor.txt"
     open(iu1_w, file=dir_file_out, status="old",err=10)
    !read(iu1_w,'(a)') line
    !if(index(line,"&atomic_coor").ne.0) then
    !   stop 'must be &atomic_coor in atomic_coor.txt'
    !endif
     do ia = 1,natom
        read(iu1_w,*) atom_name(ia), system%Rion(1:3,ia), kion(ia)
        system%Rion(1:3,ia) = system%Rion(1:3,ia)/uconv
     enddo
     close(iu1_w)
     write(*,*) "  read atomic coordinates from restart data"
  end if

  call comm_bcast(system%Rion,comm)
  call comm_bcast(atom_name,comm)
  call comm_bcast(kion,comm)
  !Rion(:,:) = system%Rion(:,:)  !remove later 

10 continue

end subroutine read_Rion

subroutine read_Velocity(idir,system)
  use structures, only: s_dft_system
  use salmon_global, only: natom
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_dft_system) :: system
  character(*),intent(in) :: idir
  integer :: iu1_w, ia, comm
  character(256) :: dir_file_out

  iu1_w = 87
  comm = nproc_group_global

  ! atomic velocity [au]
  if(comm_is_root(nproc_id_global)) then
     dir_file_out = trim(idir)//"atomic_vel.txt"
     open(iu1_w, file=dir_file_out, status="old",err=20)
     do ia = 1,natom
        read(iu1_w,*) system%Velocity(1:3,ia)
     enddo
!     if(ensemble=="NVT" .and. thermostat=="nose-hoover")then
!        read(iu1_w,*,err=100,end=100) md%xi_nh !if no value, skip reading (xi_nh=0)
!     endif
!100  continue
     close(iu1_w)
     write(*,*) "  read atomic velocities from restart data"
  end if
  call comm_bcast(system%Velocity,comm)
!  call comm_bcast(md%xi_nh,comm)

20 continue

end subroutine read_Velocity

!===================================================================================================================================

subroutine write_rtdata(wdir,itt,lg,mg,system,info,iself,rt)
  use structures
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use filesystem, only: create_directory
  use salmon_global, only: trans_longi
  implicit none
  character(*),            intent(in) :: wdir
  integer,                 intent(in) :: itt
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_parallel_info),   intent(in) :: info
  logical,                 intent(in) :: iself
  type(s_rt),              intent(in) :: rt
  !
  integer,parameter :: iunit = 333
  integer :: ierr,comm
  character(256) :: tdir,filename
  
  comm = info%icomm_rko
  
  if(trans_longi /= 'tr') then
    if(comm_is_root(nproc_id_global)) then
      filename = trim(wdir)//"rtdata.bin"
      open(iunit,file=filename,form='unformatted',iostat=ierr)
      write(iunit) rt%Ac_tot(:,itt),rt%Ac_tot(:,itt+1)
      write(iunit) rt%Ac_ind(:,itt),rt%Ac_ind(:,itt+1)
      write(iunit) rt%Ac_ind(:,0) ! for the 2D approximation
      close(iunit)
    end if
  end if
  
end subroutine write_rtdata

subroutine read_rtdata(wdir,itt,lg,mg,system,info,iself,rt)
  use structures
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: trans_longi
  implicit none
  character(*),            intent(in) :: wdir
  integer,                 intent(in) :: itt
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_parallel_info),   intent(in) :: info
  logical,                 intent(in) :: iself
  type(s_rt)                          :: rt
  !
  integer,parameter :: iunit = 333
  integer :: ierr,comm,i1,i2
  character(256) :: tdir,filename
  
  comm = info%icomm_rko
  
  if(trans_longi /= 'tr') then
    if(comm_is_root(nproc_id_global)) then
      filename = trim(wdir)//"rtdata.bin"
      open(iunit,file=filename,form='unformatted',status='old',iostat=ierr)
      read(iunit) rt%Ac_tot(:,itt),rt%Ac_tot(:,itt+1)
      read(iunit) rt%Ac_ind(:,itt),rt%Ac_ind(:,itt+1)
      read(iunit) rt%Ac_ind(:,0) ! for the 2D approximation
      close(iunit)
      write(*,*) "  read rtdata from restart data"
    end if
    call comm_bcast(rt%Ac_tot,comm)
    call comm_bcast(rt%Ac_ind,comm)
  end if

end subroutine read_rtdata

!===================================================================================================================================

#ifdef USE_MPI
#define MPI_CHECK(X) call X; call errcheck(ierr)
subroutine distributed_rw_wavefunction(iodir,lg,mg,system,info,spsi,mk,mo,if_real_orbital,rw_mode)
  use structures, only: s_rgrid, s_dft_system, s_parallel_info, s_orbital, &
  &                     allocate_orbital_real, deallocate_orbital
  use salmon_global, only: yn_periodic,method_wf_distributor,nblock_wf_distribute
  use mpi
  implicit none
  character(*),            intent(in)    :: iodir
  type(s_rgrid),           intent(in)    :: lg, mg
  type(s_dft_system),      intent(in)    :: system
  type(s_parallel_info),intent(in)    :: info
  type(s_orbital),         intent(inout) :: spsi
  integer,                 intent(in)    :: mk, mo
  logical,                 intent(in)    :: if_real_orbital
  integer,                 intent(in)    :: rw_mode

  character(256) :: iofile
  integer :: icomm
  integer :: gsize(7), lsize(7), lstart(7)
  integer :: source_type, local_type, global_type
  integer :: iopen_flag, minfo, mfile
  integer :: ierr

  type(s_orbital) :: dummy

  ! determine source type
  if (allocated(spsi%rwf)) then
    source_type = MPI_DOUBLE
  else if (allocated(spsi%zwf)) then
    if (rw_mode == read_mode .and. if_real_orbital) then
      source_type = MPI_DOUBLE
    else
      source_type = MPI_DOUBLE_COMPLEX
    end if
  end if

  call set_mpi_info

  select case(method_wf_distributor)
  case('single') ; call rw_all    ! create single shared file
  case('slice')  ; call rw_sliced ! sliced shared file
  case default   ; stop 'rw_wavefunction: fatal error'
  end select

  call free_mpi_info
  call deallocate_orbital(dummy)

contains
  subroutine rw_all
    implicit none

    iofile = trim(iodir)//"wfn.bin"
    icomm = info%icomm_rko

    ! requires data conversion from double to double complex
    if (source_type == MPI_DOUBLE) then
      call allocate_orbital_real(system%nspin,mg,info,dummy)
    end if

    ! create MPI_Type (Window) of process-local wave function
    gsize  = [mg%ie_array(1:3) - mg%is_array(1:3) + 1, system%nspin, info%numo, info%numk, 1]
    lsize  = [mg%ie(1:3)       - mg%is(1:3)       + 1, system%nspin, info%numo, info%numk, 1]
    lstart = [mg%is(1:3)       - mg%is_array(1:3) + 1, 1,            1,         1,         1] - 1

    MPI_CHECK(MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, local_type, ierr))
    MPI_CHECK(MPI_Type_commit(local_type, ierr))

    ! create MPI_Type (Window) of global wave function
    gsize  = [lg%ie(1:3) - lg%is(1:3) + 1, system%nspin, mo,        mk,        1]
    lstart = [mg%is(1:3)                 , 1,            info%io_s, info%ik_s, 1] - 1
    if (yn_periodic == 'n') then
      lstart(1:3) = lstart(1:3) + lg%num(1:3)/2
    end if

    MPI_CHECK(MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, global_type, ierr))
    MPI_CHECK(MPI_Type_commit(global_type, ierr))

    ! write/read file
    MPI_CHECK(MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr))
    MPI_CHECK(MPI_File_set_view(mfile, 0_MPI_OFFSET_KIND, local_type, global_type, 'native', MPI_INFO_NULL, ierr))

    select case(rw_mode)
      case (write_mode)
        if (allocated(spsi%rwf)) then
          MPI_CHECK(MPI_File_write_all(mfile, spsi%rwf, 1, local_type, MPI_STATUS_IGNORE, ierr))
        else if (allocated(spsi%zwf)) then
          MPI_CHECK(MPI_File_write_all(mfile, spsi%zwf, 1, local_type, MPI_STATUS_IGNORE, ierr))
        end if
      case (read_mode)
        if (allocated(spsi%rwf)) then
          if (source_type /= MPI_DOUBLE) stop 'source_type /= MPI_DOUBLE'
          MPI_CHECK(MPI_File_read_all(mfile, spsi%rwf, 1, local_type, MPI_STATUS_IGNORE, ierr))
        else if (allocated(spsi%zwf)) then
          if (source_type == MPI_DOUBLE) then
            ! convert from double to double complex
            ! NOTE: When simulating large-scale isolated system, it's possible that SALMON hangs by failing memory allocation.
            MPI_CHECK(MPI_File_read_all(mfile, dummy%rwf, 1, local_type, MPI_STATUS_IGNORE, ierr))
            spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,:,:,:) &
              = dcmplx(dummy%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,:,:,:))
          else
            MPI_CHECK(MPI_File_read_all(mfile, spsi%zwf, 1, local_type, MPI_STATUS_IGNORE, ierr))
          end if
        end if
    end select

    MPI_CHECK(MPI_File_close(mfile, ierr))
    MPI_CHECK(MPI_Type_free(global_type, ierr))
    MPI_CHECK(MPI_Type_free( local_type, ierr))
  end subroutine rw_all

  subroutine rw_sliced
    use filesystem
    use communication, only: comm_is_root
    implicit none
    integer :: nblock_orbital
    integer :: ik,io,nb,io_e,iret
    logical :: check
    type(s_parallel_info) :: dummy_info

    icomm = info%icomm_r
    nblock_orbital = min(mo,nblock_wf_distribute)

    ! requires data conversion from double to double complex
    dummy_info%io_s = 1
    dummy_info%io_e = 1 !nblock_orbital
    dummy_info%ik_s = 1
    dummy_info%ik_e = 1
    dummy_info%im_s = 1
    dummy_info%im_e = 1
    if (rw_mode == read_mode .and. if_real_orbital) then
      call allocate_orbital_real(system%nspin,mg,dummy_info,dummy)
    end if

    ! create MPI_Type (Window) of process-local wave function
    gsize  = [mg%ie_array(1:3) - mg%is_array(1:3) + 1, system%nspin, 1, 1, 1]
    lsize  = [mg%ie(1:3)       - mg%is(1:3)       + 1, system%nspin, 1, 1, 1]
    lstart = [mg%is(1:3)       - mg%is_array(1:3) + 1, 1,            1, 1, 1] - 1

    MPI_CHECK(MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, local_type, ierr))
    MPI_CHECK(MPI_Type_commit(local_type, ierr))

    ! create MPI_Type (Window) of global wave function
    gsize  = [lg%ie(1:3) - lg%is(1:3) + 1, system%nspin, 1, 1, 1]
    lstart = [mg%is(1:3)                 , 1,            1, 1, 1] - 1
    if (yn_periodic == 'n') then
      lstart(1:3) = lstart(1:3) + lg%num(1:3)/2
    end if

    MPI_CHECK(MPI_Type_create_subarray(7, gsize, lsize, lstart, MPI_ORDER_FORTRAN, source_type, global_type, ierr))
    MPI_CHECK(MPI_Type_commit(global_type, ierr))

    ! create all directory...
    if (rw_mode == write_mode) then
      if (mod(info%io_e,nblock_orbital) > 0) then
        io_e = (info%io_e / nblock_orbital + 1) * nblock_orbital
      else
        io_e = info%io_e
      end if

      do ik=info%ik_s,info%ik_e
      do io=1,mo+nblock_orbital-1,nblock_orbital
        if (mod((io - 1) / nblock_orbital, info%isize_o) == info%id_o .and. comm_is_root(info%id_r)) then
          nb = ((io - 1) / nblock_orbital) * nblock_orbital + 1
          write (iofile,'(A,I6.6,A,I6.6)') trim(iodir)//'k_',ik,'_ob_',nb
          call create_directory(iofile)
        end if
      end do
      end do

      ! check own path
      do while(.true.)
        check = .true.
        do ik=info%ik_s,info%ik_e
        do io=info%io_s,io_e,nblock_orbital
          nb = ((io - 1) / nblock_orbital) * nblock_orbital + 1
          write (iofile,'(A,I6.6,A,I6.6)') trim(iodir)//'k_',ik,'_ob_',nb
          check = check .and. directory_exists(iofile)
        end do
        end do
        if (check) exit
      end do
    end if

    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
      nb = ((io - 1) / nblock_orbital) * nblock_orbital + 1
      write (iofile,'(A,I6.6,A,I6.6,A,I6.6,A)') trim(iodir)//'k_',ik,'_ob_',nb,'/wfn_ob_',io,'.dat'

      ! write/read file
      MPI_CHECK(MPI_File_open(icomm, iofile, iopen_flag, minfo, mfile, ierr))
      MPI_CHECK(MPI_File_set_view(mfile, 0_MPI_OFFSET_KIND, local_type, global_type, 'native', MPI_INFO_NULL, ierr))

      select case(rw_mode)
        case (write_mode)
          if (allocated(spsi%rwf)) then
            MPI_CHECK(MPI_File_write_all(mfile, spsi%rwf(:,:,:,:,io,ik,1), 1, local_type, MPI_STATUS_IGNORE, ierr))
          else if (allocated(spsi%zwf)) then
            MPI_CHECK(MPI_File_write_all(mfile, spsi%zwf(:,:,:,:,io,ik,1), 1, local_type, MPI_STATUS_IGNORE, ierr))
          end if
        case (read_mode)
          if (allocated(spsi%rwf)) then
            if (source_type /= MPI_DOUBLE) stop 'source_type /= MPI_DOUBLE'
            MPI_CHECK(MPI_File_read_all(mfile, spsi%rwf(:,:,:,:,io,ik,1), 1, local_type, MPI_STATUS_IGNORE, ierr))
          else if (allocated(spsi%zwf)) then
            if (source_type == MPI_DOUBLE) then
              ! convert from double to double complex
              ! NOTE: When simulating large-scale isolated system, it's possible that
              !       SALMON hangs by failing memory allocation.
              MPI_CHECK(MPI_File_read_all(mfile, dummy%rwf(:,:,:,:,1,1,1), 1, local_type, MPI_STATUS_IGNORE, ierr))
              spsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,io,ik,1) &
                = dcmplx(dummy%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,1,1,1))
            else
              MPI_CHECK(MPI_File_read_all(mfile, spsi%zwf(:,:,:,:,io,ik,1), 1, local_type, MPI_STATUS_IGNORE, ierr))
            end if
          end if
      end select

      MPI_CHECK(MPI_File_close(mfile, ierr))
    end do
    end do

    MPI_CHECK(MPI_Type_free(global_type, ierr))
    MPI_CHECK(MPI_Type_free( local_type, ierr))
  end subroutine rw_sliced

  subroutine set_mpi_info
    implicit none

    minfo = MPI_INFO_NULL
    select case(rw_mode)
      case (write_mode)
        iopen_flag = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)

        ! NOTE: Tuning parameter for distributed file system (eg. lustre)
        !       use all OSS (Object Storage Server)
        MPI_CHECK(MPI_Info_create(minfo, ierr))
        MPI_CHECK(MPI_Info_set(minfo, 'striping_factor', '-1', ierr))
      case (read_mode)
        iopen_flag = MPI_MODE_RDONLY
      case default
        stop 'iopen_flag'
    end select
  end subroutine set_mpi_info

  subroutine free_mpi_info
    implicit none

    select case(rw_mode)
      case (write_mode)
        MPI_CHECK(MPI_Info_free(minfo, ierr))
    end select
  end subroutine

  subroutine errcheck(errcode)
    use mpi, only: MPI_MAX_ERROR_STRING, MPI_SUCCESS
    use communication, only: comm_finalize
    implicit none
    integer, intent(in) :: errcode
    character(MPI_MAX_ERROR_STRING) :: errstr
    integer                         :: retlen, ierr
    if (errcode /= MPI_SUCCESS) then
      call MPI_Error_string(errcode, errstr, retlen, ierr)
      print *, 'MPI Error:', errstr
      call comm_finalize
      stop
    end if
  end subroutine
end subroutine
#endif

end module checkpoint_restart_sub
