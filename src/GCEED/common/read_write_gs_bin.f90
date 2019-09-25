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

module read_write_gs_bin_sub
  implicit none

contains

!===================================================================================================================================

subroutine write_gs_bin(odir,lg,mg,ng,system,info,spsi,mixing,miter)
  use inputoutput, only: sysname,num_datafiles_out
  use structures, only: s_rgrid, s_dft_system, s_orbital_parallel, s_orbital, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use scf_data, only: file_out_gs_bin
  implicit none
  type(s_rgrid), intent(in)    :: lg, mg, ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(in)  :: spsi
  type(s_mixing),intent(inout) :: mixing
  integer, intent(in)          :: miter
  
  type(s_rgrid)                :: dg
  integer :: is,iob,ik,ix,iy,iz
  integer :: i,iu1_w,iu2_w
  integer :: myrank_datafiles
  real(8),allocatable :: matbox(:,:,:),matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:),cmatbox2(:,:,:)
  character(8) :: fileNumber_data
  character(100) :: file_out_gs_num_bin, dir_file_out
  character(*)   :: odir
  integer :: ifilenum_data
  integer :: version_num(2)
 
  iu1_w = 97
  iu2_w = 87
  
  if(comm_is_root(nproc_id_global))then
  
  !open file iu1_w
     dir_file_out = trim(odir)//file_out_gs_bin
     open(iu1_w,file=dir_file_out,form='unformatted')
  
  !version number
     version_num(1)=42
     version_num(2)=1
     write(iu1_w) version_num(1),version_num(2)
  
  !iteration number
     write(iu1_w) miter
  
  end if
  
  !!!!!!!!!!!!!!!!!
  ! wave function !
  !!!!!!!!!!!!!!!!!
  
  allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  
  matbox=0.d0
  cmatbox=0.d0
  
  !filenumber for writing wave function
  call set_dg(lg,mg,dg,num_datafiles_out)
  
  !filenumber for writing wave function
  if(num_datafiles_out==1)then
    ifilenum_data = iu1_w
  else
    ifilenum_data = iu2_w
  end if
  
  !open file iu2_w
  if(num_datafiles_out==1)then
  !
  else if(num_datafiles_out==-1)then
    write(fileNumber_data, '(i6.6)') myrank_datafiles
    file_out_gs_num_bin = trim(adjustl(sysname))//"_gs_"//trim(adjustl(fileNumber_data))//".bin"
    dir_file_out = trim(odir)//file_out_gs_num_bin
    open(ifilenum_data,file=dir_file_out,form='unformatted')
  else
    if(num_datafiles_out>1.and.nproc_id_global<num_datafiles_out)then
      write(fileNumber_data, '(i6.6)') myrank_datafiles
      file_out_gs_num_bin = trim(adjustl(sysname))//"_gs_"//trim(adjustl(fileNumber_data))//".bin"
      dir_file_out = trim(odir)//file_out_gs_num_bin
      open(ifilenum_data,file=dir_file_out,form='unformatted')
    end if
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
  if(num_datafiles_out==-1.or.   &
     (num_datafiles_out>1.and.nproc_id_global<num_datafiles_out))then
    close(iu2_w)
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!
  ! rho_in and rho_out !
  !!!!!!!!!!!!!!!!!!!!!!
  
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
  
  !close file iu1_w
  if(comm_is_root(nproc_id_global))then
    close(iu1_w)
  end if
  
  deallocate(matbox,matbox2)
  deallocate(cmatbox,cmatbox2)
  
end subroutine write_gs_bin

!=======================================================================

subroutine read_gs_bin(lg,mg,ng,system,info,spsi,mixing,miter)
  use inputoutput, only: theory,calc_mode,iperiodic,num_datafiles_in
  use structures, only: s_rgrid, s_dft_system,s_orbital_parallel, s_orbital, s_mixing
  use salmon_parallel, only: nproc_id_global
  use salmon_communication, only: comm_is_root, comm_summation, comm_bcast
  use scf_data, only: file_in_gs_bin
  implicit none
  type(s_rgrid),intent(in) :: lg
  type(s_rgrid),intent(in) :: mg
  type(s_rgrid),intent(in) :: ng
  type(s_dft_system),intent(in) :: system
  type(s_orbital_parallel),intent(in) :: info
  type(s_orbital), intent(inout)  :: spsi
  type(s_mixing),intent(inout) :: mixing
  integer, intent(out) :: miter
  
  type(s_rgrid)                :: dg
  real(8),allocatable :: matbox(:,:,:)
  real(8),allocatable :: matbox2(:,:,:)
  complex(8),allocatable :: cmatbox(:,:,:)
  complex(8),allocatable :: cmatbox2(:,:,:)
  integer :: version_num_box(2)
  integer :: iu1_r, iu2_r
  integer :: ifilenum_data
  integer :: ik,iob,is 
  integer :: i,ix,iy,iz
 
  iu1_r = 96
  iu2_r = 86

  if(comm_is_root(nproc_id_global))then
  
    write(*,*) file_in_gs_bin
    open(iu1_r,file=file_in_gs_bin,form='unformatted')

    read(iu1_r) version_num_box(1),version_num_box(2)
  
  end if
  
  call comm_bcast(version_num_box,info%icomm_rko)
  
  if(version_num_box(1)<=41)then
    stop 'You cannot use old restart files.'
  end if
  
  if(comm_is_root(nproc_id_global)) then
    read(iu1_r) miter
  end if
  
  call comm_bcast(miter,info%icomm_rko)
  
  !set dg
  call set_dg(lg,mg,dg,num_datafiles_in)
  
  !filenumber for reading wave function
  if(num_datafiles_in==1)then
    ifilenum_data = iu1_r
  else
    ifilenum_data = iu2_r
  end if
  
  !read wavefunction
  allocate(matbox(  lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(matbox2( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox( lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  allocate(cmatbox2(lg%is(1):lg%ie(1),lg%is(2):lg%ie(2),lg%is(3):lg%ie(3)))
  
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
            spsi%zwf(ix,iy,iz,is,iob,ik,1)=matbox(ix,iy,iz)
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
    do ik=1,system%nk
    do iob=1,system%no
    do is=1,system%nspin
      if(iperiodic==0)then
        if(num_datafiles_in==1)then
          write(*,*) nproc_id_global,iob
          if(comm_is_root(nproc_id_global))then
            read(ifilenum_data) (((matbox2(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                     iy=dg%is(2),dg%ie(2)), &
                                                     iz=dg%is(3),dg%ie(3))
          end if
          call comm_bcast(matbox2,info%icomm_rko)
        else if(num_datafiles_in>1)then
          if(nproc_id_global<num_datafiles_in)then
            read(ifilenum_data) (((matbox(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                    iy=dg%is(2),dg%ie(2)), &
                                                    iz=dg%is(3),dg%ie(3))
          end if
          call comm_summation(matbox,matbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
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
              spsi%zwf(ix,iy,iz,is,iob,ik,1) = matbox2(ix,iy,iz)
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
          call comm_bcast(cmatbox2,info%icomm_rko)
        else if(num_datafiles_in>1)then
          if(nproc_id_global<num_datafiles_in)then
            read(ifilenum_data) (((cmatbox(ix,iy,iz),ix=dg%is(1),dg%ie(1)), &
                                                     iy=dg%is(2),dg%ie(2)), &
                                                     iz=dg%is(3),dg%ie(3))
          end if
          call comm_summation(cmatbox,cmatbox2,lg%num(1)*lg%num(2)*lg%num(3),info%icomm_rko)
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
  
  !!!!!!!!!!!!!!!!!!!!!!
  ! rho_in and rho_out !
  !!!!!!!!!!!!!!!!!!!!!!
  
  if(theory=='DFT'.or.calc_mode=='GS')then
  
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
  end if
 
  deallocate(matbox,matbox2,cmatbox,cmatbox2)
 
  return

end subroutine read_gs_bin

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
  integer :: myrank_datafiles

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
      myrank_datafiles=nproc_id_global
  
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
        if(ibox==myrank_datafiles)then
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

end module read_write_gs_bin_sub
