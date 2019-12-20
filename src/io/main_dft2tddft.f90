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

#include "config.h"

subroutine main_dft2tddft
use structures
use salmon_global, only: directory_read_data,calc_mode, &
                         nproc_k,nproc_ob,nproc_domain_orbital,nproc_domain_general, &
                         target_nproc_k,target_nproc_ob,target_nproc_domain_orbital, &
                         target_nproc_domain_general, &
                         yn_periodic
use communication, only: comm_get_globalinfo,comm_bcast,comm_is_root,comm_sync_all
use salmon_xc
use timer
use initialization_sub
use init_communicator
use checkpoint_restart_sub
use set_numcpu
use filesystem, only: create_directory
implicit none
integer :: Miter

type(s_rgrid)            :: lg_scf,lg_rt
type(s_rgrid)            :: mg_scf,mg_rt
type(s_rgrid)            :: ng_scf,ng_rt
type(s_process_info)     :: pinfo_scf,pinfo_rt
type(s_orbital_parallel) :: info_scf,info_rt
type(s_field_parallel)   :: info_field_scf,info_field_rt
type(s_orbital)          :: spsi,shpsi
type(s_dft_system)       :: system_scf,system_rt
type(s_stencil)          :: stencil
type(s_ofile)            :: ofl

integer :: icomm,irank,nprocs
character(256) :: dir_file_out,gdir,pdir
integer,parameter :: fh = 41

call comm_get_globalinfo(icomm,irank,nprocs)

if(comm_is_root(irank))then
  print *, '=========================================='
  print *, 'DFT2TDDFT'
  print *, '  data redisribution for prep. TDDFT calc.'
  print *, '=========================================='
end if

call timer_begin(LOG_TOTAL)


! Read DFT calculation results
! ---------------------------------------------------------
call timer_begin(LOG_INIT_GS)
call init_dft_system(lg_scf,system_scf,stencil)

pinfo_scf%npk                 = nproc_k
pinfo_scf%nporbital           = nproc_ob
pinfo_scf%npdomain_orbital    = nproc_domain_orbital
pinfo_scf%npdomain_general    = nproc_domain_general
call init_process_distribution(system_scf,icomm,pinfo_scf)

call init_communicator_dft(icomm,pinfo_scf,info_scf,info_field_scf)
call init_grid_parallel(irank,nprocs,pinfo_scf,lg_scf,mg_scf,ng_scf)
call init_orbital_parallel_singlecell(system_scf,info_scf,pinfo_scf)

if (yn_periodic == 'y') then
  call allocate_orbital_complex(system_scf%nspin,mg_scf,info_scf,spsi)
else
  call allocate_orbital_real(system_scf%nspin,mg_scf,info_scf,spsi)
end if
call timer_end(LOG_INIT_GS)


call timer_begin(LOG_INIT_GS_RESTART)
call generate_restart_directory_name(directory_read_data,gdir,pdir)
call read_bin(pdir,lg_scf,mg_scf,ng_scf,system_scf,info_scf,spsi,Miter)
call timer_end(LOG_INIT_GS_RESTART)


! Redistribute and write to use TDDFT calculation.
! ---------------------------------------------------------

call timer_begin(LOG_INIT_RT)
calc_mode = 'RT' ! FIXME
call init_dft_system(lg_rt,system_rt,stencil)

pinfo_rt%npk                 = target_nproc_k
pinfo_rt%nporbital           = target_nproc_ob
pinfo_rt%npdomain_orbital    = target_nproc_domain_orbital
pinfo_rt%npdomain_general    = target_nproc_domain_general
call init_process_distribution(system_rt,icomm,pinfo_rt)

call init_communicator_dft(icomm,pinfo_rt,info_rt,info_field_rt)
call init_grid_parallel(irank,nprocs,pinfo_rt,lg_rt,mg_rt,ng_rt)
call init_orbital_parallel_singlecell(system_rt,info_rt,pinfo_rt)

call allocate_orbital_complex(system_rt%nspin,mg_rt,info_rt,shpsi)
call timer_end(LOG_INIT_RT)


call timer_begin(LOG_WRITE_RT_DATA)
call convert_wave_function
call timer_end(LOG_WRITE_RT_DATA)


call timer_begin(LOG_WRITE_RT_RESULTS)
call init_dir_out_restart(ofl)
call generate_restart_directory_name(ofl%dir_out_restart,gdir,pdir)
call create_directory(pdir)

if(comm_is_root(irank))then
!information
  dir_file_out = trim(pdir)//"info.bin"
  open(fh,file=dir_file_out,form='unformatted')
  write(fh) system_scf%nk
  write(fh) system_scf%no
  write(fh) Miter
  write(fh) nprocs
  close(fh)

!occupation
  dir_file_out = trim(pdir)//"occupation.bin"
  open(fh,file=dir_file_out,form='unformatted')
  write(fh) system_scf%rocc(1:system_scf%no,1:system_scf%nk,1:system_scf%nspin)
  close(fh)
end if

call write_wavefunction(pdir,lg_rt,mg_rt,system_rt,info_rt,shpsi,.false.)
call timer_end(LOG_WRITE_RT_RESULTS)

call timer_end(LOG_TOTAL)

contains

subroutine convert_wave_function
use communication, only: comm_summation, comm_create_group_byid, &
                                comm_proc_null, comm_send, comm_recv, comm_bcast
implicit none
integer,allocatable    :: is_ref(:,:,:),is_ref_srank(:,:,:),is_ref_rrank(:,:,:)
integer,allocatable    :: irank_src(:,:),irank_dst(:,:)
real(8),allocatable    :: dbuf1(:,:,:,:),dbuf2(:,:,:,:)
complex(8),allocatable :: zbuf1(:,:,:,:),zbuf2(:,:,:,:)
integer :: ik,io,is,jrank,krank
integer :: ix,iy,iz

! get referrer
allocate(is_ref      (system_scf%no,system_scf%nk,0:nprocs-1))
allocate(is_ref_srank(system_scf%no,system_scf%nk,0:nprocs-1)) ! SCF
allocate(is_ref_rrank(system_scf%no,system_scf%nk,0:nprocs-1)) ! RT

is_ref = 0
is_ref(info_scf%io_s:info_scf%io_e, &
       info_scf%ik_s:info_scf%ik_e, irank) = 1
call comm_summation(is_ref,is_ref_srank,size(is_ref),icomm)

is_ref = 0
is_ref(info_rt%io_s:info_rt%io_e, &
       info_rt%ik_s:info_rt%ik_e, irank) = 1
call comm_summation(is_ref,is_ref_rrank,size(is_ref),icomm)


! find src rank and dst rank
allocate(irank_src(system_scf%no,system_scf%nk))
allocate(irank_dst(system_scf%no,system_scf%nk))

irank_src = -1
irank_dst = -1

!$omp parallel do collapse(2) default(none) &
!$omp             private(ik,io,jrank) &
!$omp             shared(system_scf,nprocs) &
!$omp             shared(is_ref_srank,irank_src) &
!$omp             shared(is_ref_rrank,irank_dst)
do ik=1,system_scf%nk
do io=1,system_scf%no
  ! src (SCF) rank
  do jrank=0,nprocs-1
    if (is_ref_srank(io,ik,jrank) >= 1) then
      irank_src(io,ik) = jrank
      exit
    end if
  end do

  ! dst (RT) rank
  do jrank=0,nprocs-1
    if (is_ref_rrank(io,ik,jrank) >= 1) then
      irank_dst(io,ik) = jrank
      exit
    end if
  end do
end do
end do
!$omp end parallel do


! conversion
if (allocated(spsi%rwf)) then
  allocate(dbuf1(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
  allocate(dbuf2(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
end if
if (allocated(spsi%zwf)) then
  allocate(zbuf1(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
  allocate(zbuf2(lg_scf%is(1):lg_scf%ie(1), &
                 lg_scf%is(2):lg_scf%ie(2), &
                 lg_scf%is(3):lg_scf%ie(3),system_scf%nspin))
end if

do ik=1,system_scf%nk
do io=1,system_scf%no
  ! gather rgrid data
  if (info_scf%ik_s <= ik .and. ik <= info_scf%ik_e .and. &
      info_scf%io_s <= io .and. io <= info_scf%io_e) then
    if (allocated(dbuf1)) then
      dbuf1 = 0d0
!$omp parallel do collapse(3) default(none) &
!$omp             private(ix,iy,iz,is) &
!$omp             shared(system_scf,mg_scf,dbuf1,spsi) &
!$omp             shared(ik,io)
      do is=1,system_scf%nspin
      do iz=mg_scf%is(3),mg_scf%ie(3)
      do iy=mg_scf%is(2),mg_scf%ie(2)
      do ix=mg_scf%is(1),mg_scf%ie(1)
        dbuf1(ix,iy,iz,is) = spsi%rwf(ix,iy,iz,is,io,ik,1)
      end do
      end do
      end do
      end do
!$omp end parallel do
      call comm_summation(dbuf1,dbuf2,size(dbuf1),info_scf%icomm_r)
    end if
    if (allocated(zbuf1)) then
      zbuf1 = 0d0
!$omp parallel do collapse(3) default(none) &
!$omp             private(ix,iy,iz,is) &
!$omp             shared(system_scf,mg_scf,zbuf1,spsi) &
!$omp             shared(ik,io)
      do is=1,system_scf%nspin
      do iz=mg_scf%is(3),mg_scf%ie(3)
      do iy=mg_scf%is(2),mg_scf%ie(2)
      do ix=mg_scf%is(1),mg_scf%ie(1)
        zbuf1(ix,iy,iz,is) = spsi%zwf(ix,iy,iz,is,io,ik,1)
      end do
      end do
      end do
      end do
!$omp end parallel do
      call comm_summation(zbuf1,zbuf2,size(zbuf1),info_scf%icomm_r)
    end if
  end if

  ! send data from representive SCF rank to RT rank
  jrank = irank_src(io,ik)
  krank = irank_dst(io,ik)
  if (jrank >= 0 .and. krank >= 0) then
    if (jrank == krank) then
      ! sender == receiver
    else if (jrank == irank) then
      if (allocated(dbuf2)) call comm_send(dbuf2, krank, io*system_scf%nk+ik, icomm)
      if (allocated(zbuf2)) call comm_send(zbuf2, krank, io*system_scf%nk+ik, icomm)
    else if (krank == irank) then
      if (allocated(dbuf2)) call comm_recv(dbuf2, jrank, io*system_scf%nk+ik, icomm)
      if (allocated(zbuf2)) call comm_recv(zbuf2, jrank, io*system_scf%nk+ik, icomm)
    end if
  end if

  ! scatter rgrid data
  if (info_rt%ik_s <= ik .and. ik <= info_rt%ik_e .and. &
      info_rt%io_s <= io .and. io <= info_rt%io_e) then
    if (allocated(dbuf2)) then
      call comm_bcast(dbuf2, info_rt%icomm_r)
!$omp parallel do collapse(3) default(none) &
!$omp             private(ix,iy,iz,is) &
!$omp             shared(system_rt,mg_rt,dbuf2,shpsi) &
!$omp             shared(ik,io)
      do is=1,system_rt%nspin
      do iz=mg_rt%is(3),mg_rt%ie(3)
      do iy=mg_rt%is(2),mg_rt%ie(2)
      do ix=mg_rt%is(1),mg_rt%ie(1)
        shpsi%zwf(ix,iy,iz,is,io,ik,1) = cmplx(dbuf2(ix,iy,iz,is))
      end do
      end do
      end do
      end do
!$omp end parallel do
    end if
    if (allocated(zbuf2)) then
      call comm_bcast(zbuf2, info_rt%icomm_r)
!$omp parallel do collapse(3) default(none) &
!$omp             private(ix,iy,iz,is) &
!$omp             shared(system_rt,mg_rt,zbuf2,shpsi) &
!$omp             shared(ik,io)
      do is=1,system_rt%nspin
      do iz=mg_rt%is(3),mg_rt%ie(3)
      do iy=mg_rt%is(2),mg_rt%ie(2)
      do ix=mg_rt%is(1),mg_rt%ie(1)
        shpsi%zwf(ix,iy,iz,is,io,ik,1) = zbuf2(ix,iy,iz,is)
      end do
      end do
      end do
      end do
!$omp end parallel do
    end if
  end if
end do
end do

end subroutine convert_wave_function

end subroutine main_dft2tddft
