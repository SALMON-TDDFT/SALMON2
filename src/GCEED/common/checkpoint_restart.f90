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

subroutine write_gs_bin(odir,lg,mg,ng,system,info,spsi,mixing,miter)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_rgrid), intent(in)    :: lg, mg, ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(in)  :: spsi
  type(s_mixing),intent(inout) :: mixing
  integer, intent(in)          :: miter
  
  integer :: is,iob,ik
  integer :: iu1_w
  character(100) :: dir_file_out
  character(*)   :: odir
 
  iu1_w = 97
  
  !system
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"system.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')
    write(iu1_w) system%nk
    write(iu1_w) system%no

    close(iu1_w)
  end if

  !iteration number
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"iteration.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')
    write(iu1_w) miter

    close(iu1_w)
  end if

  !wave fucntion
  call write_wavefunction(odir,lg,mg,system,info,spsi)

  !occupation
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"occupation.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')
    do is=1,system%nspin
    do ik=1,system%nk
    do iob=1,system%no
      write(iu1_w) system%rocc(iob,ik,is)
    end do
    end do
    end do
    close(iu1_w)
  end if
  
  !rho_inout
  call write_rho_inout(odir,lg,ng,system,info,mixing)

end subroutine write_gs_bin

!=======================================================================

subroutine read_gs_bin(lg,mg,ng,system,info,spsi,mixing,miter)
  use inputoutput, only: theory,calc_mode,directory_read_data
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing
  use salmon_parallel, only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(inout) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(inout)  :: spsi
  type(s_mixing),intent(inout) :: mixing
  integer, intent(out) :: miter

  integer :: mk,mo

  real(8),allocatable :: roccbox(:,:,:)
  integer :: iu1_r
  integer :: ik,iob,is 
  character(100) :: dir_file_in
  
  integer :: comm
  
  comm = nproc_group_global
 
  iu1_r = 96

  !system
  !first to be read 
  if(comm_is_root(nproc_id_global))then
    dir_file_in =trim(directory_read_data)//"system.bin"
    open(iu1_r,file=dir_file_in,form='unformatted')
    read(iu1_r) mk
    read(iu1_r) mo

    close(iu1_r)
  end if
  call comm_bcast(mk,comm)
  call comm_bcast(mo,comm)

  !iteration
  if(comm_is_root(nproc_id_global))then
    dir_file_in =trim(directory_read_data)//"iteration.bin"
    open(iu1_r,file=dir_file_in,form='unformatted')
    read(iu1_r) miter

    close(iu1_r)
  end if
  call comm_bcast(miter,comm)
  
  !wave function
  call read_wavefunction(lg,mg,system,info,spsi,mk,mo)
  
  !occupation 
  allocate(roccbox(mo,mk,system%nspin))

  if(comm_is_root(nproc_id_global))then
    dir_file_in =trim(directory_read_data)//"occupation.bin"
    open(iu1_r,file=dir_file_in,form='unformatted')

    do is=1,system%nspin
    do ik=1,mk
    do iob=1,mo
      read(iu1_r) roccbox(iob,ik,is)
    end do
    end do
    end do
    system%rocc(1:system%no,1:system%nk,1:system%nspin) = roccbox(1:system%no,1:system%nk,1:system%nspin)
  end if
  call comm_bcast(system%rocc,comm)
  deallocate(roccbox)  
 
  !rho_inout
  if(theory=='DFT'.or.calc_mode=='GS')then
    call read_rho_inout(lg,ng,system,info,mixing)
  end if
 
  return

end subroutine read_gs_bin

!=======================================================================
subroutine write_wavefunction(odir,lg,mg,system,info,spsi)
  use inputoutput, only: num_datafiles_out
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)   :: odir
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(in) :: spsi
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
  integer :: iu2_w
  integer :: ifilenum_data
  type(s_rgrid) :: dg
  integer :: is,iob,ik
  integer :: ix,iy,iz
  character(100) ::  dir_file_out
  character(8) :: filenumber_data
 
  iu2_w = 87

  allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  
  matbox=0.d0
  cmatbox=0.d0
 
  !filenumber for writing wave function
  call set_dg(lg,mg,dg,num_datafiles_out)
  
  !filenumber for writing wave function
  ifilenum_data = iu2_w
  
  !open file iu2_w
  if(num_datafiles_out==1.and.comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"wfn.bin"
    open(iu2_w,file=dir_file_out,form='unformatted')
  else if(num_datafiles_out==-1.or.   &
          (num_datafiles_out>1.and.nproc_id_global<num_datafiles_out))then
    write(filenumber_data, '(i6.6)') nproc_id_global
    dir_file_out = trim(odir)//"wfn"//trim(adjustl(filenumber_data))//".bin"
    open(iu2_w,file=dir_file_out,form='unformatted')
  end if
  
  !write wavefunction
  if(num_datafiles_out==-1)then
    if(allocated(spsi%rwf))then
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
      do is=1,system%nspin
        write(ifilenum_data) ((( spsi%rwf(ix,iy,iz,is,iob,ik,1),ix=dg%is(1),dg%ie(1)), &
                                                                iy=dg%is(2),dg%ie(2)), &
                                                                iz=dg%is(3),dg%ie(3))
      end do
      end do
      end do
    else if(allocated(spsi%zwf))then
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
      do is=1,system%nspin
        write(ifilenum_data) ((( spsi%zwf(ix,iy,iz,is,iob,ik,1),ix=dg%is(1),dg%ie(1)), &
                                                                iy=dg%is(2),dg%ie(2)), &
                                                                iz=dg%is(3),dg%ie(3))
      end do
      end do
      end do
    end if
  else
    if(allocated(spsi%rwf))then
      do ik=1,1
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
        if((num_datafiles_out==1.and.comm_is_root(nproc_id_global)).or.   &
           (num_datafiles_out>1.and.nproc_id_global<num_datafiles_out))then
            write(ifilenum_data) ((( matbox2(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                       iy=dg%is(2),dg%ie(2)), &
                                                       iz=dg%is(3),dg%ie(3))
        end if 
      end do
      end do
      end do
    else if(allocated(spsi%zwf))then
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
        if((num_datafiles_out==1.and.comm_is_root(nproc_id_global)).or.   &
           (num_datafiles_out>1.and.nproc_id_global<num_datafiles_out))then
            write(ifilenum_data) ((( cmatbox2(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                        iy=dg%is(2),dg%ie(2)), &
                                                        iz=dg%is(3),dg%ie(3))
        end if 
      end do
      end do
      end do
    end if 
  end if 
  
  !close file iu2_w
  if((num_datafiles_out==1.and.comm_is_root(nproc_id_global)).or.  &
     num_datafiles_out==-1.or.   &
     (num_datafiles_out>1.and.nproc_id_global<num_datafiles_out))then
    close(iu2_w)
  end if
  
  deallocate(matbox,matbox2,cmatbox,cmatbox2)

end subroutine write_wavefunction

!=======================================================================
subroutine write_rho_inout(odir,lg,ng,system,info,mixing)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  character(*)   :: odir
  type(s_rgrid), intent(in)    :: lg,ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_mixing),intent(inout) :: mixing
  character(100) ::  dir_file_out
  integer :: i,ix,iy,iz,is
  integer :: iu1_w
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)

  iu1_w = 97
  
  if(comm_is_root(nproc_id_global))then
    dir_file_out = trim(odir)//"rho_inout.bin"
    open(iu1_w,file=dir_file_out,form='unformatted')
  end if

  allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))

  do i=1,mixing%num_rho_stock+1
    matbox2=0.d0
    matbox2(ng%is(1):ng%ie(1),   &
            ng%is(2):ng%ie(2),   &
            ng%is(3):ng%ie(3))   &
       = mixing%srho_in(i)%f(ng%is(1):ng%ie(1), &
                             ng%is(2):ng%ie(2), &
                             ng%is(3):ng%ie(3))
  
    call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
  
    if(comm_is_root(nproc_id_global))then
      write(iu1_w) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
    end if
  end do
  
  do i=1,mixing%num_rho_stock
    matbox2=0.d0
    matbox2(ng%is(1):ng%ie(1),   &
            ng%is(2):ng%ie(2),   &
            ng%is(3):ng%ie(3))   &
       = mixing%srho_out(i)%f(ng%is(1):ng%ie(1), &
                              ng%is(2):ng%ie(2), &
                              ng%is(3):ng%ie(3))
  
    call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
    if(comm_is_root(nproc_id_global))then
      write(iu1_w) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
    end if
  end do
  
  if(system%nspin == 2)then
    do is=1,2
      do i=1,mixing%num_rho_stock+1
        matbox2=0.d0
        matbox2(ng%is(1):ng%ie(1),   &
                ng%is(2):ng%ie(2),   &
                ng%is(3):ng%ie(3))   &
          = mixing%srho_s_in(i,is)%f(ng%is(1):ng%ie(1), &
                                      ng%is(2):ng%ie(2), &
                                      ng%is(3):ng%ie(3))
  
        call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
  
        if(comm_is_root(nproc_id_global))then
          write(iu1_w) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
        end if
      end do
  
      do i=1,mixing%num_rho_stock
        matbox2=0.d0
        matbox2(ng%is(1):ng%ie(1),   &
                ng%is(2):ng%ie(2),   &
                ng%is(3):ng%ie(3))   &
          = mixing%srho_s_out(i,is)%f(ng%is(1):ng%ie(1),   &
                  ng%is(2):ng%ie(2),   &
                  ng%is(3):ng%ie(3))
  
        call comm_summation(matbox2,matbox,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
  
        if(comm_is_root(nproc_id_global))then
          write(iu1_w) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
        end if
      end do
      
    end do
  
  end if
  
  if(comm_is_root(nproc_id_global))then
    close(iu1_w)
  end if

  deallocate(matbox,matbox2)
  
end subroutine write_rho_inout

!=======================================================================
subroutine read_wavefunction(lg,mg,system,info,spsi,mk,mo)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital
  use inputoutput, only: theory,calc_mode,iperiodic,num_datafiles_in,directory_read_data
  use salmon_parallel, only: nproc_id_global,nproc_group_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_rgrid), intent(in) :: lg, mg
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(inout) :: spsi
  integer,intent(in) :: mk,mo
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
  type(s_rgrid)                :: dg
  integer :: iu2_r
  integer :: ifilenum_data
  integer :: is,iob,ik
  integer :: ix,iy,iz
  character(100) :: dir_file_in
  character(8) :: filenumber_data
  integer :: comm

  iu2_r = 86

  comm = nproc_group_global
 
  !set dg
  call set_dg(lg,mg,dg,num_datafiles_in)
  
  !filenumber for reading wave function
  ifilenum_data = iu2_r
  
  !read wavefunction
  allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  
  !open file iu2_r
  if(num_datafiles_in==1.and.comm_is_root(nproc_id_global))then
    dir_file_in = trim(directory_read_data)//"wfn.bin"
    open(ifilenum_data,file=dir_file_in,form='unformatted')
  else if(num_datafiles_in==-1.or.   &
     (num_datafiles_in>1.and.nproc_id_global<num_datafiles_in))then
    write(filenumber_data, '(i6.6)') nproc_id_global
    dir_file_in = trim(directory_read_data)//"wfn"//trim(adjustl(filenumber_data))//".bin"
    open(ifilenum_data,file=dir_file_in,form='unformatted')
  end if
  
  if(num_datafiles_in==-1)then
    if(iperiodic==0)then
      if(theory=='DFT'.or.calc_mode=='GS')then
        do ik=info%ik_s,info%ik_e
        do iob=info%io_s,info%io_e
        do is=1,system%nspin
          read(ifilenum_data) (((spsi%rwf(ix,iy,iz,is,iob,ik,1),ix=dg%is(1),dg%ie(1)), &
                                                                iy=dg%is(2),dg%ie(2)), &
                                                                iz=dg%is(3),dg%ie(3))
        end do
        end do
        end do
      else
        do ik=info%ik_s,info%ik_e
        do iob=info%io_s,info%io_e
        do is=1,system%nspin
          read(ifilenum_data) (((matbox(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                  iy=dg%is(2),dg%ie(2)), &
                                                  iz=dg%is(3),dg%ie(3))
  !$omp parallel do collapse(2)
          do iz=dg%is(3),dg%ie(3)
          do iy=dg%is(3),dg%ie(3)
          do ix=dg%is(3),dg%ie(3)
            spsi%zwf(ix,iy,iz,is,iob,ik,1)=cmplx(matbox(ix,iy,iz))
          end do
          end do
          end do
        end do
        end do
        end do
      end if
    else if(iperiodic==3)then
      do ik=info%ik_s,info%ik_e
      do iob=info%io_s,info%io_e
      do is=1,system%nspin
        read(ifilenum_data) (((spsi%zwf(ix,iy,iz,is,iob,ik,1),ix=dg%is(1),dg%ie(1)), &
                                                              iy=dg%is(2),dg%ie(2)), &
                                                              iz=dg%is(3),dg%ie(3))
      end do
      end do
      end do
    end if
  else
    do ik=1,mk
    do iob=1,mo
    do is=1,system%nspin
      if(iperiodic==0)then
        if(num_datafiles_in==1)then
          if(comm_is_root(nproc_id_global))then
            read(ifilenum_data) (((matbox2(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                     iy=dg%is(2),dg%ie(2)), &
                                                     iz=dg%is(3),dg%ie(3))
          end if
          call comm_bcast(matbox2,comm)
        else if(num_datafiles_in>1)then
          if(nproc_id_global<num_datafiles_in)then
            read(ifilenum_data) (((matbox(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                    iy=dg%is(2),dg%ie(2)), &
                                                    iz=dg%is(3),dg%ie(3))
          end if
          call comm_summation(matbox,matbox2,lg%num(1)*lg%num(2)*lg%num(3),comm)
        end if
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
        if(num_datafiles_in==1)then
          if(comm_is_root(nproc_id_global))then
            read(ifilenum_data) (((cmatbox2(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                      iy=dg%is(2),dg%ie(2)), &
                                                      iz=dg%is(3),dg%ie(3))
          end if
          call comm_bcast(cmatbox2,comm)
        else if(num_datafiles_in>1)then
          if(nproc_id_global<num_datafiles_in)then
            read(ifilenum_data) (((cmatbox(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                     iy=dg%is(2),dg%ie(2)), &
                                                     iz=dg%is(3),dg%ie(3))
          end if
          call comm_summation(cmatbox,cmatbox2,lg%num(1)*lg%num(2)*lg%num(3),comm)
        end if
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
  end if
  
  !close file iu2_r
  if((num_datafiles_in==1.and.comm_is_root(nproc_id_global)).or.  &
     num_datafiles_in==-1.or.   &
     (num_datafiles_in>1.and.nproc_id_global<num_datafiles_in))then
    close(iu2_r)
  end if

end subroutine read_wavefunction

!=======================================================================
subroutine read_rho_inout(lg,ng,system,info,mixing)
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  implicit none
  type(s_rgrid), intent(in)    :: lg,ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_mixing),intent(inout) :: mixing
  integer :: iu1_r
  integer :: i,ix,iy,iz
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
 

  iu1_r = 96

  if(comm_is_root(nproc_id_global))then
    read(iu1_r) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
  end if
  call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)  
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_in(mixing%num_rho_stock+1)%f(ix,iy,iz)=matbox(ix,iy,iz)
  end do
  end do
  end do
  
  if(comm_is_root(nproc_id_global))then
    read(iu1_r) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
  end if
  call comm_bcast(matbox,info%icomm_rko)

!$omp parallel do collapse(2)  
  do iz=ng%is(3),ng%ie(3)
  do iy=ng%is(2),ng%ie(2)
  do ix=ng%is(1),ng%ie(1)
    mixing%srho_out(mixing%num_rho_stock)%f(ix,iy,iz)=matbox(ix,iy,iz)
  end do
  end do
  end do

  
  if(system%nspin == 2)then
    do i=1,mixing%num_rho_stock+1
      if(comm_is_root(nproc_id_global))then
        read(iu1_r) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
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
        read(iu1_r) ((( matbox(ix,iy,iz),ix=lg%is(1),lg%ie(1)),iy=lg%is(2),lg%ie(2)),iz=lg%is(3),lg%ie(3))
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
  end if
 
  deallocate(matbox,matbox2)

end subroutine read_rho_inout

!=======================================================================
subroutine set_dg(lg,mg,dg,num_datafiles)
  use structures, only: s_rgrid 
  use salmon_parallel, only: nproc_id_global
  implicit none 
  type(s_rgrid),intent(in)     :: lg,mg
  type(s_rgrid),intent(inout)  :: dg
  integer,intent(in)           :: num_datafiles
  integer :: i,j,j1,j2,j3
  integer :: ibox
  integer :: nproc_xyz_datafile(3)

  if(num_datafiles==1)then
    dg%is(1:3) =lg%is(1:3)
    dg%ie(1:3) =lg%ie(1:3)
    dg%num(1:3)=lg%num(1:3)
  else if(num_datafiles==-1)then
    dg%is(1:3) =mg%is(1:3)
    dg%ie(1:3) =mg%ie(1:3)
    dg%num(1:3)=mg%num(1:3)
  else
    if(nproc_id_global<num_datafiles)then
      ibox=1
      nproc_xyz_datafile=1
      do i=1,19
      do j=3,1,-1
        if(ibox<num_datafiles)then
          nproc_xyz_datafile(j)=nproc_xyz_datafile(j)*2
          ibox=ibox*2
        end if
      end do
      end do
  
      do j3=0,nproc_xyz_datafile(3)-1
      do j2=0,nproc_xyz_datafile(2)-1
      do j1=0,nproc_xyz_datafile(1)-1
        ibox = j1 + nproc_xyz_datafile(1)*j2 + nproc_xyz_datafile(1)*nproc_xyz_datafile(2)*j3 
        if(ibox==nproc_id_global)then
          dg%is(1)=j1*lg%num(1)/nproc_xyz_datafile(1)+lg%is(1)
          dg%ie(1)=(j1+1)*lg%num(1)/nproc_xyz_datafile(1)+lg%is(1)-1
          dg%is(2)=j2*lg%num(2)/nproc_xyz_datafile(2)+lg%is(2)
          dg%ie(2)=(j2+1)*lg%num(2)/nproc_xyz_datafile(2)+lg%is(2)-1
          dg%is(3)=j3*lg%num(3)/nproc_xyz_datafile(3)+lg%is(3)
          dg%ie(3)=(j3+1)*lg%num(3)/nproc_xyz_datafile(3)+lg%is(3)-1
        end if
      end do
      end do
      end do
      dg%num(:)=dg%ie(:)-dg%is(:)+1
    end if
  end if

end subroutine set_dg

!===================================================================================================================================

end module checkpoint_restart_sub
