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
!=======================================================================

#include "config.h"

module checkpoint_restart_sub
  implicit none

  integer,parameter,private :: write_mode = 1
  integer,parameter,private :: read_mode  = 2

contains

!===================================================================================================================================

subroutine init_dir_out_restart(ofile)
  use structures,  only: s_ofile
  use filesystem,  only: atomic_create_directory
  use salmon_global,   only: theory,write_rt_wfn_k
  use salmon_parallel, only: nproc_id_global,nproc_group_global
  implicit none
  type(s_ofile), intent(inout) :: ofile

  select case(theory)
    case('DFT','DFT_MD')
      ofile%dir_out_restart = 'data_for_restart/'
      call atomic_create_directory(ofile%dir_out_restart,nproc_group_global,nproc_id_global)
    case('TDDFT_response','TDDFT_pulse','Single_scale_Maxwell_TDDFT')
      if (write_rt_wfn_k == 'y') then
        ofile%dir_out_restart = 'data_for_restart_rt/'
        call atomic_create_directory(ofile%dir_out_restart,nproc_group_global,nproc_id_global)
      end if
  end select

end subroutine init_dir_out_restart


subroutine generate_checkpoint_directory_name(header,iter,gdir,pdir)
  use salmon_parallel, only: nproc_id_global
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
  use salmon_parallel, only: nproc_id_global
  implicit none
  character(*),  intent(in)  :: basedir
  character(256),intent(out) :: gdir
  character(256),intent(out) :: pdir

  ! global directory
  write(gdir,'(A,I6.6,A)')   trim(basedir)
  ! process private directory
  write(pdir,'(A,A,I6.6,A)') trim(gdir),'rank_',nproc_id_global,'/'
end subroutine generate_restart_directory_name


subroutine checkpoint_gs(lg,mg,ng,system,info,spsi,iter,mixing)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint
  use salmon_parallel, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg, ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_mixing)          ,intent(in) :: mixing

  character(256) :: gdir,wdir
  logical :: iself

  call generate_checkpoint_directory_name('gs',iter,gdir,wdir)
  call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)

  iself = (yn_self_checkpoint == 'y')
  if (iself) then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_bin(wdir,lg,mg,ng,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine checkpoint_gs

subroutine restart_gs(lg,mg,ng,system,info,spsi,iter,mixing)
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_mixing
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint
  implicit none
  type(s_rgrid)             ,intent(in) :: lg, mg, ng
  type(s_dft_system)     ,intent(inout) :: system
  type(s_orbital_parallel)  ,intent(in) :: info
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
  call read_bin(wdir,lg,mg,ng,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine restart_gs


subroutine checkpoint_rt(lg,mg,ng,system,info,spsi,iter,sVh_stock1,sVh_stock2)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_scalar
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint
  use salmon_parallel, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg, ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_scalar)          ,intent(in) :: sVh_stock1,sVh_stock2

  character(256) :: gdir,wdir
  logical :: iself

  call generate_checkpoint_directory_name('rt',iter,gdir,wdir)
  call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)

  iself = (yn_self_checkpoint == 'y')
  if (iself) then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_bin(wdir,lg,mg,ng,system,info,spsi,iter &
                ,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2,is_self_checkpoint=iself)
end subroutine checkpoint_rt

subroutine restart_rt(lg,mg,ng,system,info,spsi,iter,sVh_stock1,sVh_stock2)
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint
  implicit none
  type(s_rgrid)             ,intent(in) :: lg, mg, ng
  type(s_dft_system)     ,intent(inout) :: system
  type(s_orbital_parallel)  ,intent(in) :: info
  type(s_orbital)        ,intent(inout) :: spsi
  integer                  ,intent(out) :: iter
  type(s_scalar)         ,intent(inout) :: sVh_stock1,sVh_stock2

  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (.not. iself) then
    wdir = gdir
  end if
  call read_bin(wdir,lg,mg,ng,system,info,spsi,iter &
               ,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2,is_self_checkpoint=iself)
end subroutine restart_rt

!===================================================================================================================================

subroutine write_bin(odir,lg,mg,ng,system,info,spsi,iter,mixing,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use inputoutput, only: theory,calc_mode
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use salmon_parallel, only: nproc_id_global, nproc_size_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)            ,intent(in) :: odir
  type(s_rgrid)           ,intent(in) :: lg, mg, ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_mixing) ,optional,intent(in) :: mixing
  type(s_scalar) ,optional,intent(in) :: sVh_stock1,sVh_stock2
  logical        ,optional,intent(in) :: is_self_checkpoint

  integer :: iu1_w
  character(100) :: dir_file_out
  logical :: iself

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
  call write_wavefunction(odir,lg,mg,system,info,spsi,iself)

  !rho_inout
  if(theory=='DFT'.or.calc_mode=='GS')then
    call write_rho_inout(odir,lg,ng,system,info,mixing,iself)
  end if

  !Vh_stock
  if(theory=='TDDFT_response'.or.theory=='TDDFT_pulse'.or.calc_mode=='RT')then
    call write_Vh_stock(odir,lg,ng,info,sVh_stock1,sVh_stock2,iself)
  end if

end subroutine write_bin

!===================================================================================================================================

subroutine read_bin(idir,lg,mg,ng,system,info,spsi,iter,mixing,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use inputoutput, only: theory,calc_mode
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use salmon_parallel, only: nproc_id_global,nproc_group_global,nproc_size_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: yn_restart
  implicit none
  character(*)              ,intent(in) :: idir
  type(s_rgrid)             ,intent(in) :: lg, mg, ng
  type(s_dft_system)     ,intent(inout) :: system
  type(s_orbital_parallel)  ,intent(in) :: info
  type(s_orbital)        ,intent(inout) :: spsi
  integer                  ,intent(out) :: iter
  type(s_mixing),optional,intent(inout) :: mixing
  type(s_scalar),optional,intent(inout) :: sVh_stock1,sVh_stock2
  logical       ,optional,intent(in)    :: is_self_checkpoint

  integer :: mk,mo

  real(8),allocatable :: roccbox(:,:,:)
  integer :: iu1_r
  character(256) :: dir_file_in
  integer :: comm,itt,nprocs
  logical :: iself

  if (present(is_self_checkpoint)) then
    iself = is_self_checkpoint
  else
    iself = .false.
  end if

  comm = nproc_group_global

  iu1_r = 96

  !information
  !first to be read
  if(comm_is_root(nproc_id_global))then
    dir_file_in = trim(idir)//"info.bin"
    open(iu1_r,file=dir_file_in,form='unformatted')

    read(iu1_r) mk
    read(iu1_r) mo

    read(iu1_r) itt
    read(iu1_r) nprocs

    close(iu1_r)
  end if
  call comm_bcast(mk,comm)
  call comm_bcast(mo,comm)
  call comm_bcast(itt,comm)
  call comm_bcast(nprocs,comm)

  if((theory=='DFT'.or.calc_mode=='GS').or.  &
     ((theory=='TDDFT_response'.or.theory=='TDDFT_pulse'.or.calc_mode=='RT').and.yn_restart=='y'))then
    iter = itt
  end if

  !debug check
  if (yn_restart == 'y') then
    if (nprocs /= nproc_size_global) then
      stop 'number of processes do not match!'
    end if
  end if

  !occupation
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

  !wave function
  call read_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,is_self_checkpoint)

  !rho_inout
  if(theory=='DFT'.or.calc_mode=='GS')then
    call read_rho_inout(idir,lg,ng,system,info,mixing,is_self_checkpoint)
  end if

  !Vh_stock
  if((theory=='TDDFT_response'.or.theory=='TDDFT_pulse'.or.calc_mode=='RT').and.yn_restart=='y')then
    call read_Vh_stock(idir,lg,ng,info,sVh_stock1,sVh_stock2,is_self_checkpoint)
  end if

end subroutine read_bin


!===================================================================================================================================

subroutine write_wavefunction(odir,lg,mg,system,info,spsi,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*),            intent(in) :: odir
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
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
    call distributed_rw_wavefunction(odir,lg,mg,system,info,spsi,system%nk,system%no,write_mode)
#else
    ! single process execution
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
#endif
  end if

  !close file iu2_w
  inquire(iu2_w, opened=is_written)
  if (is_written) close(iu2_w)
end subroutine write_wavefunction

!===================================================================================================================================

subroutine write_rho_inout(odir,lg,ng,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)                        :: odir
  type(s_rgrid)           ,intent(in) :: lg,ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
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
      write(iu1_w) mixing%srho_in (i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    end do
    do i=1,mixing%num_rho_stock
      write(iu1_w) mixing%srho_out(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    end do
    if(system%nspin==2) then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
          write(iu1_w) mixing%srho_s_in (i,is)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
        end do
        do i=1,mixing%num_rho_stock
          write(iu1_w) mixing%srho_s_out(i,is)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
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

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=lg%is(3),lg%ie(3)
    do iy=lg%is(2),lg%ie(2)
    do ix=lg%is(1),lg%ie(1)
      matbox2(ix,iy,iz)=0.d0
    end do
    end do
    end do

    do i=1,mixing%num_rho_stock+1
!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        matbox2(ix,iy,iz) = mixing%srho_in(i)%f(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
      if(comm_is_root(nproc_id_global))then
        write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
    end do

    do i=1,mixing%num_rho_stock
!$omp parallel do collapse(2) private(iz,iy,ix)
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        matbox2(ix,iy,iz) = mixing%srho_out(i)%f(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
      if(comm_is_root(nproc_id_global))then
        write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      end if
    end do

    if(system%nspin == 2)then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            matbox2(ix,iy,iz) = mixing%srho_s_in(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
          call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
          if(comm_is_root(nproc_id_global))then
            write(iu1_w) matbox(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
          end if
        end do

        do i=1,mixing%num_rho_stock
!$omp parallel do collapse(2) private(iz,iy,ix)
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            matbox2(ix,iy,iz) = mixing%srho_s_out(i,is)%f(ix,iy,iz)
          end do
          end do
          end do
          call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
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

!===================================================================================================================================

subroutine write_Vh_stock(odir,lg,ng,info,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_orbital_parallel, s_scalar
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)   :: odir
  type(s_rgrid), intent(in)    :: lg,ng
  type(s_orbital_parallel),intent(in) :: info
  type(s_scalar),intent(in) :: sVh_stock1,sVh_stock2
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
    write(iu1_w) sVh_stock1%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    write(iu1_w) sVh_stock2%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    close(iu1_w)
  else
    ! write root process
    allocate(matbox0(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox1(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      matbox0(ix,iy,iz) = 0d0
    end do
    end do
    end do

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      matbox0(ix,iy,iz) = sVh_stock1%f(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(matbox0,matbox1,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

!$omp parallel do collapse(2) private(iz,iy,ix)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      matbox0(ix,iy,iz) = sVh_stock2%f(ix,iy,iz)
    end do
    end do
    end do
    call comm_summation(matbox0,matbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

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

subroutine read_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, &
  &                     allocate_orbital_real, deallocate_orbital
#ifdef USE_MPI
#else
  use salmon_global, only: yn_periodic
#endif
  implicit none
  character(*),            intent(in) :: idir
  type(s_rgrid),           intent(in) :: lg, mg
  type(s_dft_system),      intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)                     :: spsi
  integer,                 intent(in) :: mk, mo
  logical,                 intent(in) :: is_self_checkpoint

  integer :: iu2_r
  character(256) :: dir_file_in
  logical :: is_read

#ifdef USE_MPI
#else
  real(8),allocatable :: ddummy(:,:,:)
  complex(8),allocatable :: zdummy(:,:,:)
  integer :: im,ik,io,is
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
    call distributed_rw_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,read_mode)
#else
    ! single process execution
    dir_file_in = trim(idir)//"wfn.bin"
    open(iu2_r,file=dir_file_in,form='unformatted',access='stream')

    if (yn_periodic == 'n') then
      allocate(ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    else
      allocate(zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    end if

    do im=info%im_s,info%im_e
    do ik=1,mk
    do io=1,mo
    do is=1,system%nspin
      if (yn_periodic == 'n') then
        read (iu2_r) ddummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      else
        read (iu2_r) zdummy(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3))
      endif

      if (info%ik_s <= ik .and. ik <= info%ik_e .and. &
          info%io_s <= io .and. io <= info%io_e) then
        if (yn_periodic == 'n') then
          if (allocated(spsi%rwf)) then
            spsi%rwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,im) = ddummy(:,:,:)
          else
            spsi%zwf(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3),is,io,ik,im) = cmplx(ddummy(:,:,:))
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
#endif
  end if

  !close file iu2_r
  inquire(iu2_r, opened=is_read)
  if (is_read) close(iu2_r)
end subroutine read_wavefunction

!===================================================================================================================================

subroutine read_rho_inout(idir,lg,ng,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*), intent(in) :: idir
  type(s_rgrid), intent(in)    :: lg,ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
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
      read(iu1_r) mixing%srho_in (i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    end do
    do i=1,mixing%num_rho_stock
      read(iu1_r) mixing%srho_out(i)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    end do
    if(system%nspin==2) then
      do is=1,2
        do i=1,mixing%num_rho_stock+1
          read(iu1_r) mixing%srho_s_in (i,is)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
        end do
        do i=1,mixing%num_rho_stock
          read(iu1_r) mixing%srho_s_out(i,is)%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
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
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        mixing%srho_in(i)%f(ix,iy,iz)=matbox(ix,iy,iz)
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
      do iz=ng%is(3),ng%ie(3)
      do iy=ng%is(2),ng%ie(2)
      do ix=ng%is(1),ng%ie(1)
        mixing%srho_out(i)%f(ix,iy,iz)=matbox(ix,iy,iz)
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
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            mixing%srho_s_in(i,is)%f(ix,iy,iz)=matbox(ix,iy,iz)
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
          do iz=ng%is(3),ng%ie(3)
          do iy=ng%is(2),ng%ie(2)
          do ix=ng%is(1),ng%ie(1)
            mixing%srho_s_out(i,is)%f(ix,iy,iz)=matbox(ix,iy,iz)
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

!===================================================================================================================================

subroutine read_Vh_stock(idir,lg,ng,info,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_orbital_parallel, s_scalar
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*), intent(in) :: idir
  type(s_rgrid), intent(in)    :: lg,ng
  type(s_orbital_parallel),intent(in) :: info
  type(s_scalar),intent(inout) :: sVh_stock1,sVh_stock2
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
    read(iu1_r) sVh_stock1%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
    read(iu1_r) sVh_stock2%f(ng%is(1):ng%ie(1),ng%is(2):ng%ie(2),ng%is(3):ng%ie(3))
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
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      sVh_stock1%f(ix,iy,iz)=matbox1(ix,iy,iz)
    end do
    end do
    end do

!$omp parallel do collapse(2)
    do iz=ng%is(3),ng%ie(3)
    do iy=ng%is(2),ng%ie(2)
    do ix=ng%is(1),ng%ie(1)
      sVh_stock2%f(ix,iy,iz)=matbox2(ix,iy,iz)
    end do
    end do
    end do

    deallocate(matbox1,matbox2)
  end if
end subroutine read_Vh_stock

!===================================================================================================================================

#ifdef USE_MPI
#define MPI_CHECK(X) call X; call errcheck(ierr)
subroutine distributed_rw_wavefunction(iodir,lg,mg,system,info,spsi,mk,mo,rw_mode)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, &
  &                     allocate_orbital_real, deallocate_orbital
  use salmon_global, only: yn_periodic
  use mpi
  implicit none
  character(*),            intent(in)    :: iodir
  type(s_rgrid),           intent(in)    :: lg, mg
  type(s_dft_system),      intent(in)    :: system
  type(s_orbital_parallel),intent(in)    :: info
  type(s_orbital),         intent(inout) :: spsi
  integer,                 intent(in)    :: mk, mo
  integer,                 intent(in)    :: rw_mode

  character(256) :: iofile
  integer :: gsize(3), lsize(3), lstart(3)
  integer :: wfn_gsize(7), wfn_lsize(7), wfn_lstart(7)
  integer :: stype, rtype, ftype
  integer :: iopen_flag, minfo, mfile
  integer :: ierr
  type(s_orbital) :: dummy

  iofile = trim(iodir)//"wfn.bin"

  select case(rw_mode)
    case (write_mode)
      if (allocated(spsi%rwf)) then
        stype = MPI_DOUBLE
      else if (allocated(spsi%zwf)) then
        stype = MPI_DOUBLE_COMPLEX
      end if

      iopen_flag = ior(MPI_MODE_CREATE,MPI_MODE_WRONLY)

      ! Tuning parameter for distributed file system (eg. lustre)
      MPI_CHECK(MPI_Info_create(minfo, ierr))
      ! use all OSS (Object Storage Server)
      MPI_CHECK(MPI_Info_set(minfo, 'striping_factor', '-1', ierr))
    case (read_mode)
      if (yn_periodic == 'n') then
        stype = MPI_DOUBLE
      else
        stype = MPI_DOUBLE_COMPLEX
      end if

      iopen_flag = MPI_MODE_RDONLY

      minfo = MPI_INFO_NULL
  end select

  gsize(:)  = mg%ie_array(:) - mg%is_array(:) + 1
  lsize(:)  = mg%ie(:)       - mg%is(:)       + 1
  lstart(:) = mg%is(:)       - mg%is_array(:) + 1

  wfn_gsize  = [gsize (1:3), system%nspin, info%numo, info%numk, 1]
  wfn_lsize  = [lsize (1:3), system%nspin, info%numo, info%numk, 1]
  wfn_lstart = [lstart(1:3), 1,            1,         1,         1] - 1

  MPI_CHECK(MPI_Type_create_subarray(7, wfn_gsize, wfn_lsize, wfn_lstart, MPI_ORDER_FORTRAN, stype, rtype, ierr))
  MPI_CHECK(MPI_Type_commit(rtype, ierr))

  gsize(:)  = lg%ie(:) - lg%is(:) + 1
  lstart(:) = mg%is(:)
  if (yn_periodic == 'n') then
    lstart(:) = lstart(:) + lg%num(:)/2
  end if

  wfn_gsize  = [gsize (1:3), system%nspin, mo,        mk,        1]
  wfn_lstart = [lstart(1:3), 1,            info%io_s, info%ik_s, 1] - 1

  MPI_CHECK(MPI_Type_create_subarray(7, wfn_gsize, wfn_lsize, wfn_lstart, MPI_ORDER_FORTRAN, stype, ftype, ierr))
  MPI_CHECK(MPI_Type_commit(ftype, ierr))

  MPI_CHECK(MPI_File_open(MPI_COMM_WORLD, iofile, iopen_flag, minfo, mfile, ierr))
  MPI_CHECK(MPI_File_set_view(mfile, 0_MPI_OFFSET_KIND, rtype, ftype, 'native', MPI_INFO_NULL, ierr))

  select case(rw_mode)
    case (write_mode)
      if (allocated(spsi%rwf)) then
        MPI_CHECK(MPI_File_write_all(mfile, spsi%rwf, 1, rtype, MPI_STATUS_IGNORE, ierr))
      else if (allocated(spsi%zwf)) then
        MPI_CHECK(MPI_File_write_all(mfile, spsi%zwf, 1, rtype, MPI_STATUS_IGNORE, ierr))
      end if
    case (read_mode)
      if (allocated(spsi%rwf)) then
        if (stype /= MPI_DOUBLE) stop 'unsupported: stype /= MPI_DOUBLE'
        MPI_CHECK(MPI_File_read_all(mfile, spsi%rwf, 1, rtype, MPI_STATUS_IGNORE, ierr))
      else if (allocated(spsi%zwf)) then
        if (stype == MPI_DOUBLE) then
          ! read double, convert to double complex
          ! NOTE: When simulating large-scale isolated system, it's possible that
          !       SALMON hangs by failing memory allocation.
          call allocate_orbital_real(system%nspin, mg, info, dummy)
          MPI_CHECK(MPI_File_read_all(mfile, dummy%rwf, 1, rtype, MPI_STATUS_IGNORE, ierr))
          spsi%zwf = cmplx(dummy%rwf)
          call deallocate_orbital(dummy)
        else
          ! read double complex
          MPI_CHECK(MPI_File_read_all(mfile, spsi%zwf, 1, rtype, MPI_STATUS_IGNORE, ierr))
        end if
      end if
  end select

  MPI_CHECK(MPI_File_close(mfile, ierr))

  MPI_CHECK(MPI_Type_free(ftype, ierr))
  MPI_CHECK(MPI_Type_free(rtype, ierr))

  select case(rw_mode)
    case (write_mode)
      MPI_CHECK(MPI_Info_free(minfo, ierr))
  end select

contains
  subroutine errcheck(errcode)
    use mpi, only: MPI_MAX_ERROR_STRING, MPI_SUCCESS
    use salmon_communication, only: comm_finalize
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

!===================================================================================================================================

subroutine set_ndfiles(datafiles_dist, num_datafiles)
  use salmon_communication, only: comm_get_globalinfo
  implicit none
  character(*), intent(in)  :: datafiles_dist
  integer,      intent(out) :: num_datafiles
  integer :: comm, irank, nprocs

  call comm_get_globalinfo(comm,irank,nprocs)
  select case(datafiles_dist)
    case('none')
      num_datafiles = 1
    case('orbital')
      num_datafiles = nprocs
    case default
      stop 'datafiles_dist'
  end select
end subroutine

!===================================================================================================================================

end module checkpoint_restart_sub
