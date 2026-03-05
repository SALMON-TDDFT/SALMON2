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

contains

  ! yn_conventional_from_dcdft==y (SOI):
  ! conventional RT calculation with wavefunctions reconstructed from DC-LCFO(SOI) data
  subroutine init_conventional_from_dcdft_soi(lg,mg,system,info,spsi)
    use communication, only: comm_is_root, comm_summation, comm_bcast
    use filesystem, only: get_filehandle
    use salmon_global, only: num_fragment
    use structures
    implicit none
    type(s_rgrid),        intent(in) :: lg,mg
    type(s_dft_system),   intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_orbital)                  :: spsi
    character(32), parameter :: bdir_frag='./data_dcdft/fragments/'
    character(32), parameter :: binfile_wf_soi='wavefunctions_soi.bin'
    character(32), parameter :: binfile_bf_soi='basis_functions_soi.bin'
    character(32), parameter :: binfile_rg_soi='rgrid_index_soi.bin'
    character(256) :: filename
    integer :: iunit, n_frag, nspin, nstate_frag, nstate_tot
    integer :: n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
    integer :: i,j,jfrag,ispin,io,jo,ix,iy,iz,ix_tot,iy_tot,iz_tot,n
    integer,dimension(3) :: lgnum_frag,lgnum_tmp,nxyz_domain
    integer,allocatable :: n_mat(:),n_basis(:,:),index_basis(:,:,:),jxyz_tot(:,:)
    complex(8),allocatable :: f_basis(:,:,:,:,:),coef_wf(:,:,:),wrk1(:,:,:),wrk2(:,:,:)

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

    allocate(wrk1(lg%num(1),lg%num(2),lg%num(3)))
    allocate(wrk2(lg%num(1),lg%num(2),lg%num(3)))
    do ispin=1,nspin
    do io=1,system%no
      wrk1 = (0d0,0d0)
      if(io <= nstate_tot) then
        do jfrag=1,n_frag
          if(mod(jfrag-1, info%isize_rko) /= info%id_rko) cycle

          iunit = get_filehandle()
          write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_wf_soi
          open(iunit,file=filename,form='unformatted',access='stream')
          read(iunit) n_frag_file, nspin_file, nstate_frag_file, nstate_tot_file
          if(n_frag_file /= n_frag .or. nspin_file /= nspin .or. nstate_tot_file /= nstate_tot) &
            stop "data_dcdft(SOI): input mismatch"
          nstate_frag = nstate_frag_file
          allocate(n_mat(nspin))
          allocate(n_basis(n_frag,nspin))
          allocate(index_basis(nstate_frag,n_frag,nspin))
          allocate(coef_wf(nstate_frag,nstate_tot,nspin))
          read(iunit) n_mat(1:nspin)
          read(iunit) n_basis(1:n_frag,1:nspin)
          read(iunit) index_basis(1:nstate_frag,1:n_frag,1:nspin)
          read(iunit) coef_wf(1:nstate_frag,1:nstate_tot,1:nspin)
          close(iunit)

          iunit = get_filehandle()
          write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_rg_soi
          open(iunit,file=filename,form='unformatted',access='stream')
          read(iunit) lgnum_frag(1:3), lgnum_tmp(1:3)
          if(any(lgnum_tmp /= lg%num)) stop "data_dcdft(SOI): input mismatch (lg)"
          allocate(jxyz_tot(maxval(lgnum_frag),3))
          do n=1,3
            read(iunit) jxyz_tot(1:lgnum_frag(n),n)
          end do
          close(iunit)

          iunit = get_filehandle()
          write(filename, '(a, i6.6, a, a)') trim(bdir_frag), jfrag, '/', binfile_bf_soi
          open(iunit,file=filename,form='unformatted',access='stream')
          read(iunit) nxyz_domain(1:3),i,j
          read(iunit) lgnum_tmp(1:nspin)
          if(i /= nspin .or. j /= nstate_frag .or. any(lgnum_tmp(1:nspin) /= n_basis(jfrag,1:nspin))) then
            stop "data_dcdft(SOI): input mismatch (basis_functions_soi.bin)"
          end if
          allocate(f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag))
          read(iunit) f_basis(1:nxyz_domain(1),1:nxyz_domain(2),1:nxyz_domain(3),1:nspin,1:nstate_frag)
          close(iunit)

          do jo=1,n_basis(jfrag,ispin) ; j = index_basis(jo,jfrag,ispin)
          do iz=1,nxyz_domain(3); iz_tot = jxyz_tot(iz,3)
          do iy=1,nxyz_domain(2); iy_tot = jxyz_tot(iy,2)
          do ix=1,nxyz_domain(1); ix_tot = jxyz_tot(ix,1)
            wrk1(ix_tot,iy_tot,iz_tot) = wrk1(ix_tot,iy_tot,iz_tot) &
            & + f_basis(ix,iy,iz,ispin,jo) * coef_wf(jo,io,ispin)
          end do
          end do
          end do
          end do

          deallocate(n_mat,n_basis,index_basis,jxyz_tot,coef_wf,f_basis)
        end do
      end if

      call comm_summation(wrk1,wrk2,product(lg%num(1:3)),info%icomm_rko)
      if(info%io_s <= io .and. io <= info%io_e) then
        if(allocated(spsi%zwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%zwf(ix,iy,iz,ispin,io,1,1) = wrk2(ix,iy,iz)
          end do
          end do
          end do
        else if(allocated(spsi%rwf)) then
          do iz=mg%is(3),mg%ie(3)
          do iy=mg%is(2),mg%ie(2)
          do ix=mg%is(1),mg%ie(1)
            spsi%rwf(ix,iy,iz,ispin,io,1,1) = real(wrk2(ix,iy,iz))
          end do
          end do
          end do
        end if
      end if
    end do
    end do

    if(comm_is_root(info%id_rko)) then
      write(*,*) "end init_conventional_from_dcdft_soi"
    end if

    deallocate(wrk1,wrk2)
  end subroutine init_conventional_from_dcdft_soi

end module lcfo_soi_init
