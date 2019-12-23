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

module checkpoint_restart_sub
  implicit none

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
    case('DFT','DFT_MD')
      ofl%dir_out_restart = 'data_for_restart/'
      call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
    case('TDDFT_response','TDDFT_pulse','Single_scale_Maxwell_TDDFT')
      if (write_rt_wfn_k == 'y') then
        ofl%dir_out_restart = 'data_for_restart_rt/'
        call atomic_create_directory(ofl%dir_out_restart,nproc_group_global,nproc_id_global)
      end if
    case('DFT2TDDFT')
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


subroutine checkpoint_gs(lg,mg,ng,system,info,spsi,iter,mixing,odir)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint,yn_datafiles_dump
  use parallelization, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg, ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_mixing)          ,intent(in) :: mixing
  character(*),optional   ,intent(in) :: odir

  character(256) :: gdir,wdir
  logical :: iself

  if (present(odir)) then
    call generate_restart_directory_name(odir,gdir,wdir)
  else
    call generate_checkpoint_directory_name('gs',iter,gdir,wdir)
    call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)
  end if

  iself = (yn_self_checkpoint == 'y')
  if (iself .or. yn_datafiles_dump == 'y') then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_Rion(wdir,system)
  call write_Velocity(wdir,system)
  call write_bin(wdir,lg,mg,ng,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine checkpoint_gs

subroutine restart_gs(lg,mg,ng,system,info,spsi,iter,mixing)
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_mixing
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint,yn_datafiles_dump
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
  if (yn_datafiles_dump /= 'y' .and. .not. iself) then
    wdir = gdir
  end if

  call read_bin(wdir,lg,mg,ng,system,info,spsi,iter,mixing=mixing,is_self_checkpoint=iself)
end subroutine restart_gs


subroutine checkpoint_rt(lg,mg,ng,system,info,spsi,iter,sVh_stock1,sVh_stock2,idir)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_scalar
  use filesystem, only: atomic_create_directory,create_directory
  use salmon_global, only: yn_self_checkpoint,yn_datafiles_dump
  use parallelization, only: nproc_group_global,nproc_id_global
  implicit none
  type(s_rgrid)           ,intent(in) :: lg, mg, ng
  type(s_dft_system)      ,intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital)         ,intent(in) :: spsi
  integer                 ,intent(in) :: iter
  type(s_scalar)          ,intent(in) :: sVh_stock1,sVh_stock2
  character(*),optional   ,intent(in) :: idir

  character(256) :: gdir,wdir
  logical :: iself

  if (present(idir)) then
    call generate_restart_directory_name(idir,gdir,wdir)
  else
    call generate_checkpoint_directory_name('rt',iter,gdir,wdir)
    call atomic_create_directory(gdir,nproc_group_global,nproc_id_global)
  end if

  iself = (yn_self_checkpoint == 'y')
  if (iself .or. yn_datafiles_dump == 'y') then
    call create_directory(wdir)
  else
    wdir = gdir
  end if
  call write_Rion(wdir,system)
  call write_Velocity(wdir,system)
  call write_bin(wdir,lg,mg,ng,system,info,spsi,iter &
                ,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2,is_self_checkpoint=iself)
end subroutine checkpoint_rt

subroutine restart_rt(lg,mg,ng,system,info,spsi,iter,sVh_stock1,sVh_stock2)
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint,yn_datafiles_dump
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

  iself = yn_restart =='y' .and. yn_self_checkpoint == 'y'
  if (yn_datafiles_dump /= 'y' .and. .not. iself) then
    wdir = gdir
  end if

  call read_bin(wdir,lg,mg,ng,system,info,spsi,iter &
               ,sVh_stock1=sVh_stock1,sVh_stock2=sVh_stock2,is_self_checkpoint=iself)
end subroutine restart_rt


!===================================================================================================================================

subroutine write_bin(odir,lg,mg,ng,system,info,spsi,iter,mixing,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use salmon_global, only: theory,calc_mode
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use parallelization, only: nproc_id_global, nproc_size_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
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
  if((theory=='DFT'.or.calc_mode=='GS'))then
    if (present(mixing)) then
      call write_rho_inout(odir,lg,ng,system,info,mixing,iself)
    end if
  end if

  !Vh_stock
  if(theory=='TDDFT_response'.or.theory=='TDDFT_pulse'.or.calc_mode=='RT')then
    if (present(sVh_stock1) .and. present(sVh_stock2)) then
      call write_Vh_stock(odir,lg,ng,info,sVh_stock1,sVh_stock2,iself)
    end if
  end if

end subroutine write_bin

!===================================================================================================================================

subroutine read_bin(idir,lg,mg,ng,system,info,spsi,iter,mixing,sVh_stock1,sVh_stock2,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing, s_scalar
  use parallelization, only: nproc_id_global,nproc_group_global,nproc_size_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  use salmon_global, only: yn_restart, theory,calc_mode,yn_datafiles_dump
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
  if (yn_restart == 'y' .or. yn_datafiles_dump == 'y') then
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
  call read_wavefunction(idir,lg,mg,system,info,spsi,mk,mo,iself)

  !rho_inout
  if(theory=='DFT'.or.calc_mode=='GS')then
    if (present(mixing)) then
      call read_rho_inout(idir,lg,ng,system,info,mixing,iself)
    end if
  end if

  !Vh_stock
  if((theory=='TDDFT_response'.or.theory=='TDDFT_pulse'.or.calc_mode=='RT').and.yn_restart=='y')then
    if (present(sVh_stock1) .and. present(sVh_stock2)) then
      call read_Vh_stock(idir,lg,ng,info,sVh_stock1,sVh_stock2,iself)
    end if
  end if

end subroutine read_bin

!===================================================================================================================================

subroutine write_wavefunction(odir,lg,mg,system,info,spsi,is_self_checkpoint)
  use salmon_global, only: yn_datafiles_dump
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)   :: odir
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(in) :: spsi
  logical,intent(in) :: is_self_checkpoint

  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
  integer :: iu2_w
  type(s_rgrid) :: dg
  integer :: is,iob,ik
  integer :: ix,iy,iz
  character(256) ::  dir_file_out
  logical :: is_written

  iu2_w = 87

  call set_dg(lg,mg,dg,is_self_checkpoint .or. yn_datafiles_dump == 'y')

  if(is_self_checkpoint .or. yn_datafiles_dump == 'y') then
    ! write all processes (each process dump data)
    dir_file_out = trim(odir)//"wfn.bin"
    open(iu2_w,file=dir_file_out,form='unformatted',access='stream')
  else if(comm_is_root(nproc_id_global)) then
    ! write root process
    dir_file_out = trim(odir)//"wfn.bin"
    open(iu2_w,file=dir_file_out,form='unformatted')
  end if

  !write wavefunction
  if(is_self_checkpoint .or. yn_datafiles_dump == 'y')then
    if(allocated(spsi%rwf))then
      write (iu2_w) spsi%rwf(dg%is(1):dg%ie(1),   &
                             dg%is(2):dg%ie(2),   &
                             dg%is(3):dg%ie(3),   &
                             1:system%nspin,      &
                             info%io_s:info%io_e, &
                             info%ik_s:info%ik_e, &
                             1)
    else if(allocated(spsi%zwf))then
      write (iu2_w) spsi%zwf(dg%is(1):dg%ie(1),   &
                             dg%is(2):dg%ie(2),   &
                             dg%is(3):dg%ie(3),   &
                             1:system%nspin,      &
                             info%io_s:info%io_e, &
                             info%ik_s:info%ik_e, &
                             1)
    end if
  else
    if(allocated(spsi%rwf))then
      allocate(matbox (lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
      allocate(matbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

!$omp parallel do collapse(2)
      do iz=lg%is(3),lg%ie(3)
      do iy=lg%is(2),lg%ie(2)
      do ix=lg%is(1),lg%ie(1)
        matbox(ix,iy,iz) = 0d0
      end do
      end do
      end do

      do ik=1,system%nk
      do iob=1,system%no
      do is=1,system%nspin
!$omp parallel do collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          matbox(ix,iy,iz)=0.d0
        end do
        end do
        end do

        if(info%ik_s <= ik  .and. ik  <= info%ik_e .and.   &
           info%io_s <= iob .and. iob <= info%io_e) then
  !$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            matbox(ix,iy,iz) = spsi%rwf(ix,iy,iz,is,iob,ik,1)
          end do
          end do
          end do
        end if
        call comm_summation(matbox,matbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

        if(comm_is_root(nproc_id_global))then
          write(iu2_w) matbox2(dg%is(1):dg%ie(1),dg%is(2):dg%ie(2),dg%is(3):dg%ie(3))
        end if
      end do
      end do
      end do

      deallocate(matbox,matbox2)
    else if(allocated(spsi%zwf))then
      allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
      allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

!$omp parallel do collapse(2)
      do iz=lg%is(3),lg%ie(3)
      do iy=lg%is(2),lg%ie(2)
      do ix=lg%is(1),lg%ie(1)
        cmatbox(ix,iy,iz) = 0d0
      end do
      end do
      end do

      do ik=1,system%nk
      do iob=1,system%no
      do is=1,system%nspin
!$omp parallel do collapse(2)
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          cmatbox(ix,iy,iz)=0.d0
        end do
        end do
        end do

        if(info%ik_s <= ik  .and. ik  <= info%ik_e .and.   &
           info%io_s <= iob .and. iob <= info%io_e) then
  !$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            cmatbox(ix,iy,iz) = spsi%zwf(ix,iy,iz,is,iob,ik,1)
          end do
          end do
          end do
        end if
        call comm_summation(cmatbox,cmatbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)

        if(comm_is_root(nproc_id_global))then
          write(iu2_w) cmatbox2(dg%is(1):dg%ie(1),dg%is(2):dg%ie(2),dg%is(3):dg%ie(3))
        end if
      end do
      end do
      end do

      deallocate(cmatbox,cmatbox2)
    end if
  end if

  !close file iu2_w
  inquire(iu2_w, opened=is_written)
  if (is_written) close(iu2_w)
end subroutine write_wavefunction

!===================================================================================================================================

subroutine write_rho_inout(odir,lg,ng,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
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
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
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
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital
  use salmon_global, only: iperiodic,yn_datafiles_dump
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*),intent(in) :: idir
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(inout) :: spsi
  integer,intent(in) :: mk,mo
  logical,intent(in) :: is_self_checkpoint

  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
  type(s_rgrid)                :: dg
  integer :: iu2_r
  integer :: is,iob,ik
  integer :: ix,iy,iz
  character(256) :: dir_file_in
  integer :: comm
  logical :: is_read

  iu2_r = 86
  comm = nproc_group_global

  call set_dg(lg,mg,dg,is_self_checkpoint .or. yn_datafiles_dump == 'y')

  if(is_self_checkpoint .or. yn_datafiles_dump == 'y') then
    ! read all processes (each process load dumped data)
    dir_file_in = trim(idir)//"wfn.bin"
    open(iu2_r,file=dir_file_in,form='unformatted',access='stream')
  else if(comm_is_root(nproc_id_global))then
    ! read root process
    dir_file_in = trim(idir)//"wfn.bin"
    open(iu2_r,file=dir_file_in,form='unformatted')
  end if

  if(is_self_checkpoint .or. yn_datafiles_dump == 'y')then
    if (allocated(spsi%rwf)) then
      read (iu2_r) spsi%rwf(dg%is(1):dg%ie(1),   &
                            dg%is(2):dg%ie(2),   &
                            dg%is(3):dg%ie(3),   &
                            1:system%nspin,      &
                            info%io_s:info%io_e, &
                            info%ik_s:info%ik_e, &
                            1)
    else if (allocated(spsi%zwf)) then
      read (iu2_r) spsi%zwf(dg%is(1):dg%ie(1),   &
                            dg%is(2):dg%ie(2),   &
                            dg%is(3):dg%ie(3),   &
                            1:system%nspin,      &
                            info%io_s:info%io_e, &
                            info%ik_s:info%ik_e, &
                            1)
    end if
  else
    allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
    allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

    do ik=1,mk
    do iob=1,mo
    do is=1,system%nspin
      if(iperiodic==0)then
        if(comm_is_root(nproc_id_global))then
          read(iu2_r) matbox2(dg%is(1):dg%ie(1),dg%is(2):dg%ie(2),dg%is(3):dg%ie(3))
        end if
        call comm_bcast(matbox2,comm)
        if(info%ik_s <= ik  .and. ik  <= info%ik_e .and.   &
           info%io_s <= iob .and. iob <= info%io_e) then
          if(allocated(spsi%rwf))then
  !$omp parallel do collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              spsi%rwf(ix,iy,iz,is,iob,ik,1) = matbox2(ix,iy,iz)
            end do
            end do
            end do
          else
  !$omp parallel do collapse(2)
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              spsi%zwf(ix,iy,iz,is,iob,ik,1) = cmplx(matbox2(ix,iy,iz))
            end do
            end do
            end do
          end if
        end if
      else if(iperiodic==3)then
        if(comm_is_root(nproc_id_global))then
          read(iu2_r) cmatbox2(dg%is(1):dg%ie(1),dg%is(2):dg%ie(2),dg%is(3):dg%ie(3))
        end if
        call comm_bcast(cmatbox2,comm)
        if(info%ik_s <= ik  .and. ik  <= info%ik_e .and.   &
           info%io_s <= iob .and. iob <= info%io_e) then
  !$omp parallel do collapse(2)
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,is,iob,ik,1) = cmatbox2(ix,iy,iz)
          end do
          end do
          end do
        end if
      end if
    end do
    end do
    end do

    deallocate(matbox,matbox2,cmatbox,cmatbox2)
  end if

  !close file iu2_r
  inquire(iu2_r, opened=is_read)
  if (is_read) close(iu2_r)
end subroutine read_wavefunction

!===================================================================================================================================

subroutine read_rho_inout(idir,lg,ng,system,info,mixing,is_self_checkpoint)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
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
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
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
!(currently not used: see subroutine "read_atomic_coordinates")
subroutine restart_Rion(system)
  use structures, only: s_dft_system
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint,yn_datafiles_dump
  implicit none
  type(s_dft_system)     ,intent(inout) :: system
  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (yn_datafiles_dump /= 'y' .and. .not. iself) then
    wdir = gdir
  end if

  call read_Rion(wdir,system)
end subroutine restart_Rion

subroutine restart_Velocity(system)
  use structures, only: s_dft_system
  use communication, only: comm_is_root
  use parallelization, only: nproc_id_global
  use salmon_global, only: directory_read_data,yn_restart,yn_self_checkpoint,yn_datafiles_dump
  implicit none
  type(s_dft_system), intent(inout) :: system
  character(256) :: gdir,wdir
  logical :: iself

  call generate_restart_directory_name(directory_read_data,gdir,wdir)

  iself = (yn_restart =='y' .and. yn_self_checkpoint == 'y')
  if (yn_datafiles_dump /= 'y' .and. .not. iself) then
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
  use inputoutput, only: au_length_aa
  use parallelization, only: nproc_id_global
  use communication, only: comm_is_root
  implicit none
  type(s_dft_system),intent(in) :: system
  integer :: iu1_w, ia
  real(8) :: uconv
  character(*)   :: odir
  character(256) :: dir_file_out

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

7000 format("'",a,"'  ",3f18.10,i4,'   ',a)
7100 format("'",a,"'  ",3f18.10,i4)
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
  use salmon_global, only: natom, atom_name, kion, rion, unit_length
  use inputoutput, only: au_length_aa
  use parallelization, only: nproc_id_global,nproc_group_global
  use communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_dft_system) :: system
  character(*),intent(in) :: idir
  real(8) :: uconv
  integer :: iu1_w, ia, comm
  character(256) :: dir_file_out

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
  Rion(:,:) = system%Rion(:,:)  !remove later 

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

subroutine set_dg(lg,mg,dg,is_self_checkpoint)
  use structures, only: s_rgrid 
  implicit none 
  type(s_rgrid),intent(in)     :: lg,mg
  type(s_rgrid),intent(inout)  :: dg
  logical,intent(in)           :: is_self_checkpoint

  if(is_self_checkpoint)then
    dg%is(1:3) =mg%is(1:3)
    dg%ie(1:3) =mg%ie(1:3)
    dg%num(1:3)=mg%num(1:3)
  else 
    dg%is(1:3) =lg%is(1:3)
    dg%ie(1:3) =lg%ie(1:3)
    dg%num(1:3)=lg%num(1:3)
  end if

end subroutine set_dg

!===================================================================================================================================

end module checkpoint_restart_sub
