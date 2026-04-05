!
!  Copyright 2019-2026 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90-------100-------110-------120-------130

module lcfo_soi_init
  implicit none
  private
  public :: init_conventional_from_dcdft_soi

  character(32), parameter :: bdir_frag='./data_dcdft/fragments/'
  character(32), parameter :: binfile_wf_soi='wavefunctions_soi.bin'
  character(32), parameter :: binfile_bf_soi='basis_functions_soi.bin'
  character(32), parameter :: binfile_rg_soi='rgrid_index_soi.bin'

  type s_dc_fragment_cache_soi
    integer :: nxyz_domain(3) = 0
    integer :: lgnum_frag(3) = 0
    integer,allocatable :: n_basis_frag(:)
    integer,allocatable :: index_basis(:,:)
    integer,allocatable :: jxyz_tot(:,:)
    complex(8),allocatable :: coef_wf(:,:,:)
    complex(8),allocatable :: f_basis(:,:,:,:,:)
  end type s_dc_fragment_cache_soi

contains

  ! yn_conventional_from_dcdft==y (SOI):
  ! conventional RT calculation with wavefunctions reconstructed from DC-LCFO(SOI) data
  subroutine init_conventional_from_dcdft_soi(lg,mg,system,info,spsi)
    use communication, only: comm_is_root, comm_summation, comm_bcast
    use filesystem, only: get_filehandle
    use salmon_global, only: num_fragment, lcfo_frag_cache_size
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)                  :: spsi
    character(256) :: filename
    integer :: iunit, n_frag, nspin, nstate_tot
    integer :: n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    integer,parameter :: io_block_size = 2
    integer :: i,j,jfrag,ispin,io,jo,ix,iy,iz,ix_tot,iy_tot,iz_tot,n
    integer :: io_blk_s, io_blk_e, io_blk, nblk
    integer :: n_frag_local, ifrag_local, ifrag_blk_s, ifrag_blk_e, ncache
    integer :: frag_cache_batch_size
    integer,dimension(3) :: lgnum_tmp
    integer,allocatable :: n_mat(:),n_basis(:,:),index_basis(:,:,:),jxyz_tot(:,:)
    complex(8),allocatable :: wrk1(:,:,:,:),wrk2(:,:,:,:)
    type(s_dc_fragment_cache_soi),allocatable :: frag_cache(:)
    integer,allocatable :: local_frag_ids(:)

    nspin = system%nspin
    n_frag = product(num_fragment)
    nstate_tot = 0

    if(comm_is_root(info%id_rko)) then
      write(*,*) 'start init_conventional_from_dcdft_soi'
      write(*,*) "yn_conventional_from_dcdft==y (SOI): reconstruct from DC-LCFO(SOI) data"
      write(*,*) "read from ./data_dcdft/fragments/*/*_soi.bin"
    end if

    if(info%isize_rko < 1) stop "yn_conventional_from_dcdft=y (SOI): invalid MPI size."

    if(info%id_rko == mod(1-1, info%isize_rko)) then
      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), 1, '/', binfile_wf_soi
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      close(iunit)
      if(n_frag_file /= n_frag .or. nspin_file /= nspin) stop "data_dcdft(SOI): input mismatch"
      nstate_tot = nstate_tot_file
    end if
    call comm_bcast(nstate_tot, info%icomm_rko)
    frag_cache_batch_size = max(1, lcfo_frag_cache_size)

    n_frag_local = 0
    do jfrag=1,n_frag
      if(mod(jfrag-1, info%isize_rko) == info%id_rko) n_frag_local = n_frag_local + 1
    end do
    allocate(local_frag_ids(n_frag_local))
    ifrag_local = 0
    do jfrag=1,n_frag
      if(mod(jfrag-1, info%isize_rko) /= info%id_rko) cycle
      ifrag_local = ifrag_local + 1
      local_frag_ids(ifrag_local) = jfrag
    end do

    allocate(wrk1(lg%num(1),lg%num(2),lg%num(3),io_block_size))
    allocate(wrk2(lg%num(1),lg%num(2),lg%num(3),io_block_size))
    do ispin=1,nspin
      do io=info%io_s,info%io_e
        if(allocated(spsi%zwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,ispin,io,1,1) = (0d0,0d0)
          end do
          end do
          end do
        else if(allocated(spsi%rwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%rwf(ix,iy,iz,ispin,io,1,1) = 0d0
          end do
          end do
          end do
        end if
      end do

      do ifrag_blk_s=1,n_frag_local,frag_cache_batch_size
        ifrag_blk_e = min(ifrag_blk_s + frag_cache_batch_size - 1, n_frag_local)
        ncache = ifrag_blk_e - ifrag_blk_s + 1
        allocate(frag_cache(ncache))
        call load_fragment_cache_batch_soi(local_frag_ids, ifrag_blk_s, ifrag_blk_e, n_frag, nspin, nstate_tot, lg%num, frag_cache)
        do io_blk_s=1,system%no,io_block_size
          io_blk_e = min(io_blk_s + io_block_size - 1, system%no)
          nblk = io_blk_e - io_blk_s + 1
          wrk1(:,:,:,1:nblk) = (0d0,0d0)
          if(io_blk_s <= nstate_tot) then
            do ifrag_local=1,ncache
              do io_blk=1,nblk
                io = io_blk_s + io_blk - 1
                if(io > nstate_tot) cycle
                do jo=1,frag_cache(ifrag_local)%n_basis_frag(ispin) ; j = frag_cache(ifrag_local)%index_basis(jo,ispin)
                do iz=1,frag_cache(ifrag_local)%nxyz_domain(3); iz_tot = frag_cache(ifrag_local)%jxyz_tot(iz,3)
                do iy=1,frag_cache(ifrag_local)%nxyz_domain(2); iy_tot = frag_cache(ifrag_local)%jxyz_tot(iy,2)
                do ix=1,frag_cache(ifrag_local)%nxyz_domain(1); ix_tot = frag_cache(ifrag_local)%jxyz_tot(ix,1)
                  wrk1(ix_tot,iy_tot,iz_tot,io_blk) = wrk1(ix_tot,iy_tot,iz_tot,io_blk) &
                  & + frag_cache(ifrag_local)%f_basis(ix,iy,iz,ispin,jo) * frag_cache(ifrag_local)%coef_wf(jo,io,ispin)
                end do
                end do
                end do
                end do
              end do
            end do
          end if
          call comm_summation(wrk1,wrk2,product(lg%num(1:3))*nblk,info%icomm_rko)
          do io_blk=1,nblk
            io = io_blk_s + io_blk - 1
            if(info%io_s <= io .and. io <= info%io_e) then
              if(allocated(spsi%zwf)) then
                do iz=mg%is(3),mg%ie(3)
                do iy=mg%is(2),mg%ie(2)
                do ix=mg%is(1),mg%ie(1)
                  spsi%zwf(ix,iy,iz,ispin,io,1,1) = spsi%zwf(ix,iy,iz,ispin,io,1,1) + wrk2(ix,iy,iz,io_blk)
                end do
                end do
                end do
              else if(allocated(spsi%rwf)) then
                do iz=mg%is(3),mg%ie(3)
                do iy=mg%is(2),mg%ie(2)
                do ix=mg%is(1),mg%ie(1)
                  spsi%rwf(ix,iy,iz,ispin,io,1,1) = spsi%rwf(ix,iy,iz,ispin,io,1,1) + real(wrk2(ix,iy,iz,io_blk))
                end do
                end do
                end do
              end if
            end if
          end do
        end do
        call free_fragment_cache_batch_soi(frag_cache)
        deallocate(frag_cache)
      end do
    end do

    if(comm_is_root(info%id_rko)) then
      write(*,*) "end init_conventional_from_dcdft_soi"
    end if

    deallocate(wrk1,wrk2)
    if(allocated(local_frag_ids)) deallocate(local_frag_ids)
  end subroutine init_conventional_from_dcdft_soi

  subroutine load_fragment_cache_batch_soi(local_frag_ids, ifrag_blk_s, ifrag_blk_e, n_frag, nspin, nstate_tot, lgnum_tot, frag_cache)
    use filesystem, only: get_filehandle
    implicit none
    integer,intent(in) :: local_frag_ids(:)
    integer,intent(in) :: ifrag_blk_s, ifrag_blk_e
    integer,intent(in) :: n_frag, nspin, nstate_tot, lgnum_tot(3)
    type(s_dc_fragment_cache_soi),intent(inout) :: frag_cache(:)
    integer :: iunit, n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    integer :: ifrag_local, jfrag, i, j, n
    integer,allocatable :: n_mat(:),n_basis(:,:),index_basis(:,:,:)
    integer :: lgnum_tmp(3)
    character(256) :: filename

    do ifrag_local=1,size(frag_cache)
      jfrag = local_frag_ids(ifrag_blk_s + ifrag_local - 1)

      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_wf_soi
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
      if(n_frag_file /= n_frag .or. nspin_file /= nspin .or. nstate_tot_file /= nstate_tot) &
        stop "data_dcdft(SOI): input mismatch"
      allocate(n_mat(nspin))
      allocate(n_basis(n_frag,nspin))
      allocate(index_basis(nstate_frag_file,n_frag,nspin))
      allocate(frag_cache(ifrag_local)%coef_wf(nstate_frag_file,nstate_tot,nspin))
      read(iunit) n_mat(1:nspin)
      read(iunit) n_basis(1:n_frag,1:nspin)
      read(iunit) index_basis(1:nstate_frag_file,1:n_frag,1:nspin)
      read(iunit) frag_cache(ifrag_local)%coef_wf(1:nstate_frag_file,1:nstate_tot,1:nspin)
      close(iunit)

      allocate(frag_cache(ifrag_local)%n_basis_frag(nspin))
      allocate(frag_cache(ifrag_local)%index_basis(nstate_frag_file,nspin))
      frag_cache(ifrag_local)%n_basis_frag(1:nspin) = n_basis(jfrag,1:nspin)
      frag_cache(ifrag_local)%index_basis(1:nstate_frag_file,1:nspin) = index_basis(1:nstate_frag_file,jfrag,1:nspin)
      deallocate(n_mat,n_basis,index_basis)

      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_rg_soi
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) frag_cache(ifrag_local)%lgnum_frag(1:3), lgnum_tmp(1:3)
      if(any(lgnum_tmp /= lgnum_tot)) stop "data_dcdft(SOI): input mismatch (lg)"
      allocate(frag_cache(ifrag_local)%jxyz_tot(maxval(frag_cache(ifrag_local)%lgnum_frag),3))
      do n=1,3
        read(iunit) frag_cache(ifrag_local)%jxyz_tot(1:frag_cache(ifrag_local)%lgnum_frag(n),n)
      end do
      close(iunit)

      iunit = get_filehandle()
      write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_bf_soi
      open(iunit,file=filename,form='unformatted',access='stream')
      read(iunit) frag_cache(ifrag_local)%nxyz_domain(1:3),i,j
      read(iunit) lgnum_tmp(1:nspin)
      if(i /= nspin .or. j /= nstate_frag_file .or. any(lgnum_tmp(1:nspin) /= frag_cache(ifrag_local)%n_basis_frag(1:nspin))) then
        stop "data_dcdft(SOI): input mismatch (basis_functions_soi.bin)"
      end if
      allocate(frag_cache(ifrag_local)%f_basis(1:frag_cache(ifrag_local)%nxyz_domain(1), &
      & 1:frag_cache(ifrag_local)%nxyz_domain(2),1:frag_cache(ifrag_local)%nxyz_domain(3),1:nspin,1:nstate_frag_file))
      read(iunit) frag_cache(ifrag_local)%f_basis(1:frag_cache(ifrag_local)%nxyz_domain(1), &
      & 1:frag_cache(ifrag_local)%nxyz_domain(2),1:frag_cache(ifrag_local)%nxyz_domain(3),1:nspin,1:nstate_frag_file)
      close(iunit)
    end do
  end subroutine load_fragment_cache_batch_soi

  subroutine free_fragment_cache_batch_soi(frag_cache)
    implicit none
    type(s_dc_fragment_cache_soi),intent(inout) :: frag_cache(:)
    integer :: ifrag_local

    do ifrag_local=1,size(frag_cache)
      if(allocated(frag_cache(ifrag_local)%n_basis_frag)) deallocate(frag_cache(ifrag_local)%n_basis_frag)
      if(allocated(frag_cache(ifrag_local)%index_basis)) deallocate(frag_cache(ifrag_local)%index_basis)
      if(allocated(frag_cache(ifrag_local)%jxyz_tot)) deallocate(frag_cache(ifrag_local)%jxyz_tot)
      if(allocated(frag_cache(ifrag_local)%coef_wf)) deallocate(frag_cache(ifrag_local)%coef_wf)
      if(allocated(frag_cache(ifrag_local)%f_basis)) deallocate(frag_cache(ifrag_local)%f_basis)
    end do
  end subroutine free_fragment_cache_batch_soi

end module lcfo_soi_init
