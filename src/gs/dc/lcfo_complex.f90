!
!  Copyright 2026 SALMON developers
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
module lcfo_complex
  implicit none

  private
  public :: check_lcfo_complex_options, dc_lcfo_complex
  public :: init_conventional_from_dcdft_complex

  real(8), parameter :: lcfo_tol = 1d-10

  type :: s_lcfo_complex_halo
    integer :: id_src = -1
    integer :: id_dst = -1
    integer :: ifrag_src = 0
    integer :: dvec(3) = 0
    integer :: length(3) = 0
    integer :: dsp_send(3) = 0
    integer :: dsp_recv(3) = 0
    complex(8), allocatable :: buf_send(:,:,:,:,:)
    complex(8), allocatable :: buf_recv(:,:,:,:,:)
    complex(8), allocatable :: mat_H_local(:,:,:)
  end type s_lcfo_complex_halo

  type :: s_complex_lcfo_writer
    integer :: unit_basis = -1
    integer :: unit_coeff = -1
    integer(8) :: bytes_basis = 0
    integer(8) :: bytes_coeff = 0
    character(96) :: run_id = ''
    character(256) :: file_basis = ''
    character(256) :: file_coeff = ''
  end type s_complex_lcfo_writer

  type :: s_complex_lcfo_fragment
    integer, allocatable :: n_basis(:), jxyz(:,:)
    complex(8), allocatable :: basis(:,:,:,:,:), coef(:,:,:)
  end type s_complex_lcfo_fragment

contains

  subroutine open_complex_lcfo_files(lg,dc,writer,local_status)
    use communication, only: comm_bcast
    use filesystem, only: get_filehandle
    use ieee_arithmetic, only: ieee_support_datatype
    use iso_fortran_env, only: int32,int64,real64,file_storage_size
    use salmon_global, only: base_directory,num_fragment
    use structures, only: s_dcdft,s_rgrid
    implicit none
    type(s_rgrid), intent(in) :: lg
    type(s_dcdft), intent(in) :: dc
    type(s_complex_lcfo_writer), intent(out) :: writer
    integer, intent(out) :: local_status
    integer :: clock_values(8),ios,iu_rgrid,n,close_status
    integer(int64) :: clock_count,clock_rate
    logical :: rgrid_open

    local_status = 0
    if (storage_size(0) /= 32 .or. storage_size(0d0) /= 64 .or. &
        storage_size(cmplx(0d0,0d0,8)) /= 128 .or. file_storage_size /= 8) then
      local_status = 1
      return
    end if
    if (.not.ieee_support_datatype(0d0)) then
      local_status = 1
      return
    end if

    if (dc%id_tot == 0) then
      call date_and_time(values=clock_values)
      call system_clock(count=clock_count,count_rate=clock_rate)
      write(writer%run_id,'(8(i4.4),":",i0,":",i0)') &
           clock_values,clock_count,clock_rate
    end if
    call comm_bcast(writer%run_id,dc%icomm_tot,0)

    if (dc%id_frag /= 0) return
    writer%file_basis = trim(base_directory)//'basis_functions.bin'
    writer%file_coeff = trim(base_directory)//'wavefunctions.bin'
    writer%unit_basis = get_filehandle()
    open(writer%unit_basis,file=trim(writer%file_basis),status='replace', &
         form='unformatted',access='stream',action='write',iostat=ios)
    if (ios /= 0) then
      local_status = 1
      return
    end if
    call write_complex_lcfo_header(writer%unit_basis,1_int32,writer%run_id,lg,dc,ios)
    if (ios /= 0) then
      local_status = 1
      close(writer%unit_basis,iostat=ios)
      writer%unit_basis = -1
      return
    end if
    writer%unit_coeff = get_filehandle()
    open(writer%unit_coeff,file=trim(writer%file_coeff),status='replace', &
         form='unformatted',access='stream',action='write',iostat=ios)
    if (ios /= 0) then
      local_status = 1
      close(writer%unit_basis,iostat=ios)
      writer%unit_basis = -1
      return
    end if
    call write_complex_lcfo_header(writer%unit_coeff,2_int32,writer%run_id,lg,dc,ios)
    if (ios /= 0) then
      local_status = 1
      close(writer%unit_basis,iostat=ios)
      close(writer%unit_coeff,iostat=ios)
      writer%unit_basis = -1
      writer%unit_coeff = -1
      return
    end if
    writer%bytes_basis = 304_int64 + 32_int64*int(dc%system_tot%nk,int64)
    writer%bytes_coeff = writer%bytes_basis

    iu_rgrid = get_filehandle()
    rgrid_open = .false.
    open(iu_rgrid,file=trim(base_directory)//'rgrid_index.bin', &
         status='replace',form='unformatted',access='stream',action='write',iostat=ios)
    if (ios == 0) rgrid_open = .true.
    if (ios == 0) write(iu_rgrid,iostat=ios) lg%num(1:3),dc%lg_tot%num(1:3)
    if (ios == 0) then
      do n=1,3
        write(iu_rgrid,iostat=ios) dc%jxyz_tot(1:lg%num(n),n)
        if (ios /= 0) exit
      end do
    end if
    close_status = 0
    if (rgrid_open) close(iu_rgrid,iostat=close_status)
    if (ios /= 0 .or. close_status /= 0) local_status = 1
  end subroutine open_complex_lcfo_files

  subroutine write_complex_lcfo_header(iu,file_kind,run_id,lg,dc,ios)
    use iso_fortran_env, only: int32,int64,real64
    use salmon_global, only: num_fragment
    use structures, only: s_dcdft,s_rgrid
    implicit none
    integer, intent(in) :: iu
    integer(int32), intent(in) :: file_kind
    character(96), intent(in) :: run_id
    type(s_rgrid), intent(in) :: lg
    type(s_dcdft), intent(in) :: dc
    integer, intent(out) :: ios
    character(16), parameter :: magic='SLCFO_COMPLEX_V1'
    integer(int32) :: control(6),meta(20)
    integer(int64) :: header_bytes
    real(real64) :: real_meta(10)

    control = [1_int32,int(z'01020304',int32),32_int32,64_int32,file_kind,0_int32]
    meta(1:3) = int(dc%lg_tot%num,int32)
    meta(4:6) = int(lg%num,int32)
    meta(7:9) = int(dc%nxyz_domain,int32)
    meta(10:12) = int(num_fragment,int32)
    meta(13:15) = int(dc%jxyz_tot(1,1:3),int32)
    meta(16:20) = [int(dc%i_frag,int32),int(dc%system_tot%nspin,int32), &
         int(dc%system_tot%nk,int32),int(dc%nstate_frag,int32),int(dc%nstate_tot,int32)]
    real_meta(1:9) = reshape(real(dc%system_tot%primitive_a,real64),[9])
    real_meta(10) = real(dc%system_tot%hvol,real64)
    header_bytes = 304_int64 + 32_int64*int(dc%system_tot%nk,int64)
    write(iu,iostat=ios) magic
    if (ios == 0) write(iu,iostat=ios) control
    if (ios == 0) write(iu,iostat=ios) header_bytes
    if (ios == 0) write(iu,iostat=ios) run_id
    if (ios == 0) write(iu,iostat=ios) meta
    if (ios == 0) write(iu,iostat=ios) real_meta
    if (ios == 0) write(iu,iostat=ios) real(dc%system_tot%vec_k,real64)
    if (ios == 0) write(iu,iostat=ios) real(dc%system_tot%wtk,real64)
  end subroutine write_complex_lcfo_header

  subroutine write_complex_lcfo_k_record(writer,ik,dc,n_basis,n_mat,index_basis,basis,coef,local_status)
    use iso_fortran_env, only: int32,int64
    use structures, only: s_dcdft
    implicit none
    type(s_complex_lcfo_writer), intent(inout) :: writer
    integer, intent(in) :: ik
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: n_basis(:,:,:),n_mat(:,:),index_basis(:,:,:,:)
    complex(8), intent(in) :: basis(:,:,:,:,:),coef(:,:,:)
    integer, intent(out) :: local_status
    integer(int32) :: ik_wire,ispin_wire,nb_wire,nmat_wire
    integer(int32), allocatable :: meta_basis(:),rows(:)
    integer(int64) :: basis_payload,coeff_payload
    integer :: ispin,nb,ios

    local_status = 0
    if (dc%id_frag /= 0) return
    call complex_lcfo_payload_sizes(dc,ik,n_basis,n_mat,basis_payload,coeff_payload,local_status)
    if (local_status /= 0) return
    ik_wire = int(ik,int32)
    write(writer%unit_basis,iostat=ios) ik_wire,basis_payload
    if (ios == 0) write(writer%unit_coeff,iostat=ios) ik_wire,coeff_payload
    do ispin=1,dc%system_tot%nspin
      nb = n_basis(dc%i_frag,ispin,ik)
      nb_wire = int(nb,int32)
      ispin_wire = int(ispin,int32)
      if (ios == 0) write(writer%unit_basis,iostat=ios) ispin_wire,nb_wire
      if (ios == 0 .and. nb > 0) then
        call write_wire_complex_array(writer%unit_basis, &
             reshape(basis(:,:,:,ispin,1:nb),[size(basis,1)*size(basis,2)*size(basis,3)*nb]),ios)
      end if
      nmat_wire = int(n_mat(ispin,ik),int32)
      if (ios == 0) write(writer%unit_coeff,iostat=ios) ispin_wire,nb_wire,nmat_wire
      if (ios == 0) write(writer%unit_coeff,iostat=ios) int(n_basis(:,ispin,ik),int32)
      if (ios == 0 .and. nb > 0) then
        allocate(rows(nb))
        rows = int(index_basis(1:nb,dc%i_frag,ispin,ik),int32)
        write(writer%unit_coeff,iostat=ios) rows
        deallocate(rows)
      end if
      if (ios == 0 .and. nb > 0) then
        call write_wire_complex_array(writer%unit_coeff, &
             reshape(coef(1:nb,1:dc%nstate_tot,ispin),[nb*dc%nstate_tot]),ios)
      end if
      if (ios /= 0) exit
    end do
    if (ios /= 0) then
      local_status = 1
      return
    end if
    if (writer%bytes_basis > huge(writer%bytes_basis)-12_int64-basis_payload .or. &
        writer%bytes_coeff > huge(writer%bytes_coeff)-12_int64-coeff_payload) then
      local_status = 1
      return
    end if
    writer%bytes_basis = writer%bytes_basis+12_int64+basis_payload
    writer%bytes_coeff = writer%bytes_coeff+12_int64+coeff_payload
  end subroutine write_complex_lcfo_k_record

  subroutine complex_lcfo_payload_sizes(dc,ik,n_basis,n_mat,basis_payload,coeff_payload,status)
    use iso_fortran_env, only: int64
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    integer, intent(in) :: ik,n_basis(:,:,:),n_mat(:,:)
    integer(int64), intent(out) :: basis_payload,coeff_payload
    integer, intent(out) :: status
    integer(int64) :: factors(5)
    integer :: ispin,nb

    basis_payload = 0_int64
    coeff_payload = 0_int64
    status = 0
    do ispin=1,dc%system_tot%nspin
      nb = n_basis(dc%i_frag,ispin,ik)
      if (nb < 0 .or. nb > dc%nstate_frag .or. n_mat(ispin,ik) < dc%nstate_tot) then
        status = 1
        return
      end if
      call safe_add_factors(basis_payload,[8_int64],status)
      factors = [int(dc%nxyz_domain(1),int64),int(dc%nxyz_domain(2),int64), &
           int(dc%nxyz_domain(3),int64),int(nb,int64),16_int64]
      call safe_add_factors(basis_payload,factors,status)
      call safe_add_factors(coeff_payload,[12_int64],status)
      call safe_add_factors(coeff_payload,[4_int64,int(dc%n_frag,int64)],status)
      factors = [int(nb,int64),int(dc%nstate_tot,int64),16_int64,1_int64,1_int64]
      call safe_add_factors(coeff_payload,factors,status)
      call safe_add_factors(coeff_payload,[4_int64,int(nb,int64)],status)
      if (status /= 0) return
    end do
    if (basis_payload > huge(basis_payload)-12_int64 .or. &
        coeff_payload > huge(coeff_payload)-12_int64) status = 1
  end subroutine complex_lcfo_payload_sizes

  subroutine safe_add_factors(total,factors,status)
    use iso_fortran_env, only: int64
    implicit none
    integer(int64), intent(inout) :: total
    integer(int64), intent(in) :: factors(:)
    integer, intent(inout) :: status
    integer(int64) :: value
    integer :: i
    value = 1_int64
    do i=1,size(factors)
      if (factors(i) < 0_int64) then
        status = 1
        return
      end if
      if (factors(i) /= 0_int64) then
        if (value > huge(value)/factors(i)) then
          status = 1
          return
        end if
      end if
      value = value*factors(i)
    end do
    if (total > huge(total)-value) then
      status = 1
      return
    end if
    total = total+value
  end subroutine safe_add_factors

  subroutine write_wire_complex_array(iu,values,ios)
    use iso_fortran_env, only: real64
    implicit none
    integer, intent(in) :: iu
    complex(8), intent(in) :: values(:)
    integer, intent(inout) :: ios
    integer, parameter :: chunk=65536
    real(real64), allocatable :: wire(:)
    integer :: first,last,j,n
    if (ios /= 0) return
    allocate(wire(2*min(chunk,max(1,size(values)))))
    first = 1
    do while(first <= size(values))
      last = min(size(values),first+chunk-1)
      n = last-first+1
      do j=1,n
        wire(2*j-1) = real(values(first+j-1),real64)
        wire(2*j) = aimag(values(first+j-1))
      end do
      write(iu,iostat=ios) wire(1:2*n)
      if (ios /= 0) exit
      first = last+1
    end do
    deallocate(wire)
  end subroutine write_wire_complex_array

  subroutine finish_complex_lcfo_files(dc,writer,local_status)
    use iso_fortran_env, only: int32,int64
    use structures, only: s_dcdft
    implicit none
    type(s_dcdft), intent(in) :: dc
    type(s_complex_lcfo_writer), intent(inout) :: writer
    integer, intent(out) :: local_status
    character(16), parameter :: footer_magic='SLCFO_DONE_V1'
    integer(int32) :: completed_k
    integer(int64) :: file_bytes,actual_bytes
    integer :: ios,close_status

    local_status = 0
    if (dc%id_frag /= 0) return
    completed_k = int(dc%system_tot%nk,int32)
    file_bytes = writer%bytes_basis+124_int64
    write(writer%unit_basis,iostat=ios) footer_magic,writer%run_id,completed_k,file_bytes
    if (ios == 0) then
      writer%bytes_basis = file_bytes
      close(writer%unit_basis,iostat=close_status)
      if (close_status /= 0) ios = close_status
      writer%unit_basis = -1
    end if
    if (ios == 0) then
      inquire(file=trim(writer%file_basis),size=actual_bytes,iostat=ios)
      if (ios == 0 .and. actual_bytes /= writer%bytes_basis) ios = 1
    end if
    file_bytes = writer%bytes_coeff+124_int64
    if (ios == 0) write(writer%unit_coeff,iostat=ios) &
         footer_magic,writer%run_id,completed_k,file_bytes
    if (ios == 0) then
      writer%bytes_coeff = file_bytes
      close(writer%unit_coeff,iostat=close_status)
      if (close_status /= 0) ios = close_status
      writer%unit_coeff = -1
    end if
    if (ios == 0) then
      inquire(file=trim(writer%file_coeff),size=actual_bytes,iostat=ios)
      if (ios == 0 .and. actual_bytes /= writer%bytes_coeff) ios = 1
    end if
    if (ios /= 0) local_status = 1
    if (writer%unit_basis >= 0) then
      close(writer%unit_basis,iostat=close_status)
      writer%unit_basis = -1
      if (close_status /= 0) local_status = 1
    end if
    if (writer%unit_coeff >= 0) then
      close(writer%unit_coeff,iostat=close_status)
      writer%unit_coeff = -1
      if (close_status /= 0) local_status = 1
    end if
  end subroutine finish_complex_lcfo_files

  subroutine check_lcfo_complex_options(system,info,dc)
    use ieee_arithmetic, only: ieee_is_finite
    use salmon_global, only: energy_cut,lambda_cut,lcfo_eigensolver, &
         num_fragment,yn_dc_lcfo_diag,yn_spinorbit,theory
    use structures, only: s_dcdft,s_dft_system,s_parallel_info
    implicit none
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_dcdft), intent(in) :: dc
    integer :: ik,n

#if defined(USE_OPENACC) || defined(USE_CUDA)
    stop "DC-LCFO complex: GPU/OpenACC/CUDA is unsupported."
#endif
    if (yn_spinorbit == 'y') stop "DC-LCFO complex: spin-orbit and noncollinear calculations are unsupported."
    if (trim(theory) /= 'dft') stop "DC-LCFO complex: theory must be dft."
    select case(trim(lcfo_eigensolver))
    case('lapack')
    case('chefsi')
#ifndef USE_SCALAPACK
      stop "DC-LCFO complex: lcfo_eigensolver='chefsi' requires ScaLAPACK."
#endif
    case default
      stop "DC-LCFO complex: supported eigensolvers are 'lapack' and 'chefsi'."
    end select
    if (yn_dc_lcfo_diag /= 'y') &
      stop "DC-LCFO complex: yn_dc_lcfo_diag='y' is required."
    if (system%if_real_orbital) stop "DC-LCFO complex: complex orbitals are required."
    if (system%nspin < 1 .or. system%nspin > 2) stop "DC-LCFO complex: unsupported spin count."
    if (.not.ieee_is_finite(system%hvol) .or. system%hvol <= 0d0) &
      stop "DC-LCFO complex: invalid grid volume."
    if (system%nk < 1 .or. dc%nstate_tot < 1 .or. dc%nstate_frag < 1) &
      stop "DC-LCFO complex: invalid k-point or state count."
    if (info%npk < 1 .or. info%npk > system%nk .or. info%nporbital < 1 .or. &
        info%nporbital > system%no) stop "DC-LCFO complex: invalid MPI distribution."
    if (info%numk < 1) stop "DC-LCFO complex: every k-point group must own a k point."
    if (dc%system_tot%nk /= system%nk) stop "DC-LCFO complex: total and fragment nk differ."
    if (.not.allocated(system%vec_k) .or. .not.allocated(system%wtk)) &
      stop "DC-LCFO complex: k-point data is not allocated."
    if (.not.allocated(dc%system_tot%vec_k) .or. .not.allocated(dc%system_tot%wtk)) &
      stop "DC-LCFO complex: total k-point data is not allocated."
    if (maxval(abs(system%vec_k-dc%system_tot%vec_k)) > lcfo_tol .or. &
        maxval(abs(system%wtk-dc%system_tot%wtk)) > lcfo_tol) &
      stop "DC-LCFO complex: fragment and total k-point lists differ."
    if (.not.ieee_is_finite(energy_cut) .or. .not.ieee_is_finite(lambda_cut) .or. &
        lambda_cut <= 0d0) stop "DC-LCFO complex: invalid LCFO cutoff."
    if (any(dc%nxyz_domain < 1) .or. any(dc%nxyz_buffer < 0) .or. &
        any(dc%nxyz_buffer > dc%nxyz_domain)) stop "DC-LCFO complex: invalid domain/buffer."
    if ((dc%id_frag == 0) .neqv. (info%id_rko == 0)) &
      stop "DC-LCFO complex: fragment root and rko root do not agree."
    if (dc%n_frag /= product(num_fragment)) stop "DC-LCFO complex: fragment count mismatch."
    if (dc%nstate_frag /= system%no) stop "DC-LCFO complex: nstate_frag and system%no differ."
    do n=1,3
      if (num_fragment(n) > 1) then
        if (any(abs(system%vec_k(n,1:system%nk)) > lcfo_tol)) &
          stop "DC-LCFO complex: k component along a fragment-split direction must be zero."
      else if (dc%nxyz_buffer(n) /= 0) then
        stop "DC-LCFO complex: buffer is allowed only along split directions."
      end if
    end do
    do ik=1,system%nk
      if (.not.ieee_is_finite(system%wtk(ik))) stop "DC-LCFO complex: non-finite k weight."
    end do
  end subroutine check_lcfo_complex_options

  subroutine dc_lcfo_complex(lg,mg,system,info,stencil,ppg,energy,v_local,spsi,shpsi,sttpsi,srg,dc)
    use communication, only: comm_bcast,comm_irecv,comm_isend,comm_summation,comm_wait_all
    use eigen_subdiag_sub, only: eigen_zheev
    use hamiltonian, only: hpsi
    use ieee_arithmetic, only: ieee_is_finite
    use filesystem, only: get_filehandle
    use math_constants, only: pi
    use salmon_global, only: energy_cut,lambda_cut,sysname, &
         lcfo_eigensolver,lcfo_diag_chefsi_filter_degree, &
         lcfo_diag_chefsi_filter_chunk_size,lcfo_diag_chefsi_max_cycle, &
         lcfo_diag_chefsi_residual_tolerance,yn_dc_lcfo_diag
    use structures
    implicit none
    type(s_rgrid),         intent(in) :: lg,mg
    type(s_dft_system),    intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_stencil),       intent(in) :: stencil
    type(s_pp_grid),       intent(in) :: ppg
    type(s_dft_energy),    intent(in) :: energy
    type(s_scalar),        intent(in) :: v_local(system%nspin)
    type(s_orbital),       intent(in) :: spsi
    type(s_orbital)                  :: shpsi,sttpsi
    type(s_sendrecv_grid)             :: srg
    type(s_dcdft)                     :: dc

    type(s_lcfo_complex_halo) :: halo(26)
    type(s_complex_lcfo_writer) :: writer
    integer, allocatable :: id_array(:),n_basis(:,:,:),n_mat(:,:),index_basis(:,:,:,:)
    integer, allocatable :: req_send(:),req_recv(:)
    real(8), allocatable :: esp_tot(:,:,:)
    complex(8), allocatable :: f_basis(:,:,:,:,:),work_basis(:,:,:,:,:)
    complex(8), allocatable :: hf(:,:,:,:,:),work_hf(:,:,:,:,:)
    complex(8), allocatable :: mat_h_local(:,:,:)
    complex(8), allocatable :: hsend(:,:),hmat(:,:),vmat(:,:)
    complex(8), allocatable :: coef_frag(:,:,:)
    real(8), allocatable :: eval(:)
    integer :: n_halo,ik,ispin,io,ix,iy,iz,n,nactive,istat,global_status
    integer :: nspin,nk,m,pass,j,ifrag,i,jj,src_frag,tag,ios
    integer, allocatable :: nb(:)
    real(8) :: hvol,herm_err,basis_err,ortho_err,resid_err,norm_h
    logical :: owns_ik

    call check_lcfo_complex_options(system,info,dc)
    if (.not.allocated(spsi%zwf) .or. .not.allocated(sttpsi%zwf) .or. &
        .not.allocated(shpsi%zwf)) stop "DC-LCFO complex: complex orbital scratch is not allocated."
    if (.not.allocated(energy%esp)) stop "DC-LCFO complex: SCF eigenvalues are not allocated."
    if (.not.ieee_is_finite(system%mu)) stop "DC-LCFO complex: chemical potential is not finite."
    if (any(shape(energy%esp) /= [system%no,system%nk,system%nspin])) &
      stop "DC-LCFO complex: SCF eigenvalue shape is inconsistent."
    if (any(.not.ieee_is_finite(energy%esp))) &
      stop "DC-LCFO complex: SCF eigenvalues are not finite."

    if (dc%id_tot == 0) write(*,*) "start DC-LCFO complex"
    hvol = system%hvol
    nspin = system%nspin
    nk = system%nk
    m = dc%nstate_frag

    allocate(id_array(dc%n_frag))
    call init_fragment_roots(id_array)
    call init_halo(n_halo,halo,id_array)
    allocate(n_basis(dc%n_frag,nspin,nk),n_mat(nspin,nk))
    allocate(index_basis(m,dc%n_frag,nspin,nk),esp_tot(dc%nstate_tot,nspin,nk))
    allocate(coef_frag(m,dc%nstate_tot,nspin))
    n_basis = 0
    n_mat = 0
    index_basis = 0
    esp_tot = 0d0
    coef_frag = (0d0,0d0)
    if (yn_dc_lcfo_diag == 'y') then
      call open_complex_lcfo_files(lg,dc,writer,istat)
      call check_collective_status(istat,"open complex LCFO output",0,0)
    end if

    do ik=1,nk
      owns_ik = (info%ik_s <= ik .and. ik <= info%ik_e)
      allocate(f_basis(dc%nxyz_domain(1),dc%nxyz_domain(2),dc%nxyz_domain(3),nspin,m))
      allocate(nb(nspin))
      call build_basis(ik,f_basis,nb,basis_err)
      call collect_basis_dimensions(ik,nb)
      if (any(n_mat(:,ik) < dc%nstate_tot)) &
        stop "DC-LCFO complex: nstate_tot exceeds a k-dependent Hamiltonian dimension."

      sttpsi%zwf = (0d0,0d0)
      shpsi%zwf = (0d0,0d0)
      if (owns_ik) then
        do io=info%io_s,info%io_e
          if (io <= m) then
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              if (ix <= dc%nxyz_domain(1) .and. iy <= dc%nxyz_domain(2) .and. &
                  iz <= dc%nxyz_domain(3)) &
                sttpsi%zwf(ix,iy,iz,1:nspin,io,ik,1) = f_basis(ix,iy,iz,1:nspin,io)
            end do
            end do
            end do
          end if
        end do
      end if

      allocate(hf(lg%num(1),lg%num(2),lg%num(3),nspin,m))
      allocate(work_hf(lg%num(1),lg%num(2),lg%num(3),nspin,m))
      call hpsi(sttpsi,shpsi,info,mg,v_local,system,stencil,srg,ppg)
      work_hf = (0d0,0d0)
      if (owns_ik) then
        do io=info%io_s,info%io_e
          if (io <= m) then
            do ispin=1,nspin
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              work_hf(ix,iy,iz,ispin,io) = shpsi%zwf(ix,iy,iz,ispin,io,ik,1)
            end do
            end do
            end do
            end do
          end if
        end do
      end if
      call comm_summation(work_hf,hf,size(hf),info%icomm_rko)
      deallocate(work_hf)

      allocate(mat_h_local(m,m,nspin))
      call assemble_halo_hamiltonian(ik,f_basis,hf,mat_h_local,n_halo,halo,nactive, &
           req_send,req_recv)
      deallocate(hf)
      coef_frag = (0d0,0d0)

      if (trim(lcfo_eigensolver) == 'chefsi') then
#ifdef USE_SCALAPACK
        call diag_chefsi_complex_driver(ik)
#else
        stop "DC-LCFO complex: CheFSI requires a ScaLAPACK build."
#endif
      else
      do ispin=1,nspin
        n = n_mat(ispin,ik)
        allocate(hsend(n,n),hmat(n,n),vmat(n,n),eval(n))
        hsend = (0d0,0d0)
        if (dc%id_frag == 0) then
          ifrag = dc%i_frag
          do io=1,n_basis(ifrag,ispin,ik)
            i = index_basis(io,ifrag,ispin,ik)
            do jj=1,n_basis(ifrag,ispin,ik)
              j = index_basis(jj,ifrag,ispin,ik)
              hsend(i,j) = mat_h_local(io,jj,ispin)
            end do
          end do
          do j=1,n_halo
            src_frag = halo(j)%ifrag_src
            if (.not.allocated(halo(j)%mat_h_local)) cycle
            do io=1,n_basis(ifrag,ispin,ik)
              i = index_basis(io,ifrag,ispin,ik)
              do jj=1,n_basis(src_frag,ispin,ik)
                n = index_basis(jj,src_frag,ispin,ik)
                hsend(n,i) = hsend(n,i) + 0.5d0*halo(j)%mat_h_local(jj,io,ispin)
                hsend(i,n) = hsend(i,n) + &
                     0.5d0*conjg(halo(j)%mat_h_local(jj,io,ispin))
              end do
            end do
          end do
        end if
        call comm_summation(hsend,hmat,size(hmat),dc%icomm_tot)
        herm_err = maxval(abs(hmat-conjg(transpose(hmat)))) / &
             max(1d0,maxval(abs(hmat)))
        istat = 0
        if (.not.ieee_is_finite(herm_err) .or. herm_err > lcfo_tol) istat = 1
        call check_collective_status(istat,"Hermitian Hamiltonian",ik,ispin)
        hmat = 0.5d0*(hmat+conjg(transpose(hmat)))
        do i=1,n
          hmat(i,i) = dcmplx(real(hmat(i,i),8),0d0)
        end do
        vmat = (0d0,0d0)
        istat = 0
        if (dc%id_tot == 0) then
          call eigen_zheev(hmat,eval,vmat,lapack_info=istat)
          if (istat /= 0) istat = 1
          if (istat == 0) call check_eigensystem(hmat,vmat,eval,dc%nstate_tot, &
               ortho_err,resid_err)
          if (istat == 0 .and. (ortho_err > lcfo_tol .or. resid_err > lcfo_tol)) istat = 1
          if (dc%id_tot == 0) write(*,'(a,2i6,4(1x,es12.4))') &
               "complex LCFO k/spin, hermitian/basis/eigen/residual:",ik,ispin, &
               herm_err,basis_err,ortho_err,resid_err
        end if
        call check_collective_status(istat,"ZHEEV",ik,ispin)
        call comm_bcast(eval,dc%icomm_tot,0)
        call comm_bcast(vmat,dc%icomm_tot,0)
        esp_tot(:,ispin,ik) = eval(1:dc%nstate_tot)
        do io=1,n_basis(dc%i_frag,ispin,ik)
          i = index_basis(io,dc%i_frag,ispin,ik)
          coef_frag(io,1:dc%nstate_tot,ispin) = vmat(i,1:dc%nstate_tot)
        end do
        deallocate(hsend,hmat,vmat,eval)
      end do
      end if
      if (yn_dc_lcfo_diag == 'y') then
        call write_complex_lcfo_k_record(writer,ik,dc,n_basis,n_mat,index_basis, &
             f_basis,coef_frag,istat)
        call check_collective_status(istat,"write complex LCFO k record",ik,0)
      end if
      deallocate(f_basis)
      deallocate(mat_h_local)
      call deallocate_halo_buffers(n_halo,halo)
      deallocate(nb)
    end do

    if (yn_dc_lcfo_diag == 'y') then
      call finish_complex_lcfo_files(dc,writer,istat)
      call check_collective_status(istat,"finish complex LCFO output",nk,0)
    end if

    call write_complex_eigenvalues(esp_tot,n_basis,n_mat)
    deallocate(coef_frag,esp_tot,index_basis,n_mat,n_basis,id_array)
    if (dc%id_tot == 0) write(*,*) "end DC-LCFO complex"

  contains

    subroutine init_fragment_roots(ids)
      use communication, only: comm_summation
      implicit none
      integer, intent(out) :: ids(:)
      integer :: local(size(ids))
      local = 0
      if (dc%id_frag == 0) local(dc%i_frag) = dc%id_tot + 1
      call comm_summation(local,ids,size(ids),dc%icomm_tot)
      ids = ids - 1
      if (any(ids < 0)) stop "DC-LCFO complex: fragment root map is incomplete."
    end subroutine init_fragment_roots

    subroutine init_halo(nh,halos,ids)
      use salmon_global, only: num_fragment
      implicit none
      integer, intent(out) :: nh
      type(s_lcfo_complex_halo), intent(out) :: halos(:)
      integer, intent(in) :: ids(:)
      integer :: lx,ly,lz,h,ifg,n,d(3),ir1(3),ir2(3),delta(3)
      integer :: nhdir(3)
      nhdir = 0
      do n=1,3
        if (num_fragment(n) > 1) nhdir(n) = 1
      end do
      h = 0
      do lx=-nhdir(1),nhdir(1)
      do ly=-nhdir(2),nhdir(2)
      do lz=-nhdir(3),nhdir(3)
        if (lx == 0 .and. ly == 0 .and. lz == 0) cycle
        h = h + 1
        halos(h)%dvec = [lx,ly,lz]
        halos(h)%id_dst = -1
        halos(h)%id_src = -1
        do ifg=1,dc%n_frag
          ir1 = dc%ixyz_frag(:,ifg)
          ir2 = dc%ixyz_frag(:,dc%i_frag) + halos(h)%dvec*dc%nxyz_domain
          delta = mod(ir1-ir2,dc%lg_tot%num)
          if (all(delta == 0) .and. halos(h)%id_dst < 0) halos(h)%id_dst = ids(ifg)
          ir2 = dc%ixyz_frag(:,dc%i_frag) - halos(h)%dvec*dc%nxyz_domain
          delta = mod(ir1-ir2,dc%lg_tot%num)
          if (all(delta == 0) .and. halos(h)%id_src < 0) then
            halos(h)%id_src = ids(ifg)
            halos(h)%ifrag_src = ifg
          end if
        end do
        if (halos(h)%id_dst < 0 .or. halos(h)%id_src < 0) &
          stop "DC-LCFO complex: halo neighbor was not found."
        do n=1,3
          select case(halos(h)%dvec(n))
          case(0)
            halos(h)%length(n) = dc%nxyz_domain(n)
            halos(h)%dsp_send(n) = 0
            halos(h)%dsp_recv(n) = 0
          case(1)
            halos(h)%length(n) = dc%nxyz_buffer(n)
            halos(h)%dsp_send(n) = dc%nxyz_domain(n)-dc%nxyz_buffer(n)
            halos(h)%dsp_recv(n) = dc%nxyz_domain(n)+dc%nxyz_buffer(n)
          case(-1)
            halos(h)%length(n) = dc%nxyz_buffer(n)
            halos(h)%dsp_send(n) = 0
            halos(h)%dsp_recv(n) = dc%nxyz_domain(n)
          end select
        end do
      end do
      end do
      end do
      nh = h
    end subroutine init_halo

    subroutine build_basis(ik0,basis,nb0,basis_err0)
      use communication, only: comm_bcast,comm_summation
      use eigen_subdiag_sub, only: eigen_zheev
      implicit none
      integer, intent(in) :: ik0
      complex(8), intent(out) :: basis(:,:,:,:,:)
      integer, intent(out) :: nb0(:)
      real(8), intent(out) :: basis_err0
      complex(8), allocatable :: local(:,:,:,:,:),work(:,:,:,:,:),smat(:,:,:),umat(:,:,:)
      real(8), allocatable :: lambda(:,:)
      integer :: io0,jo0,isp0,ix0,iy0,iz0,ieig,ib,pass0,status0
      complex(8) :: coeff
      real(8) :: norm2

      allocate(local(size(basis,1),size(basis,2),size(basis,3),nspin,m))
      allocate(work(size(basis,1),size(basis,2),size(basis,3),nspin,m))
      allocate(smat(m,m,nspin),umat(m,m,nspin),lambda(m,nspin))
      local = (0d0,0d0)
      do io0=info%io_s,info%io_e
        if (io0 <= m .and. info%ik_s <= ik0 .and. ik0 <= info%ik_e) then
          do isp0=1,nspin
          do iz0=mg%is(3),mg%ie(3)
          do iy0=mg%is(2),mg%ie(2)
          do ix0=mg%is(1),mg%ie(1)
            if (ix0 <= dc%nxyz_domain(1) .and. iy0 <= dc%nxyz_domain(2) .and. &
                iz0 <= dc%nxyz_domain(3)) then
              if (energy%esp(io0,ik0,isp0)-system%mu < energy_cut) &
                local(ix0,iy0,iz0,isp0,io0) = spsi%zwf(ix0,iy0,iz0,isp0,io0,ik0,1)
            end if
          end do
          end do
          end do
          end do
        end if
      end do
      call comm_summation(local,basis,size(basis),info%icomm_rko)
      deallocate(local)
      smat = (0d0,0d0)
      umat = (0d0,0d0)
      lambda = 0d0
      status0 = 0
      if (dc%id_frag == 0) then
        status0 = 0
        do isp0=1,nspin
        do io0=1,m
        do jo0=1,m
          smat(io0,jo0,isp0) = hvol*sum(conjg(basis(:,:,:,isp0,io0))*basis(:,:,:,isp0,jo0))
        end do
        end do
        end do
        status0 = 0
        do isp0=1,nspin
          call eigen_zheev(smat(:,:,isp0),lambda(:,isp0),umat(:,:,isp0),lapack_info=status0)
          if (status0 /= 0) exit
          if (any(.not.ieee_is_finite(lambda(:,isp0))) .or. &
              minval(lambda(:,isp0)) < -lcfo_tol*max(1d0,maxval(abs(lambda(:,isp0))))) then
            status0 = 1
            exit
          end if
        end do
      end if
      call check_collective_status(status0,"fragment overlap ZHEEV",ik0,0)
      work = basis
      basis = (0d0,0d0)
      nb0 = 0
      if (dc%id_frag == 0) then
        do isp0=1,nspin
          ib = 0
          do ieig=m,1,-1
            if (lambda(ieig,isp0) > lambda_cut) then
              ib = ib + 1
              do jo0=1,m
                basis(:,:,:,isp0,ib) = basis(:,:,:,isp0,ib) + &
                     work(:,:,:,isp0,jo0)*umat(jo0,ieig,isp0)/sqrt(lambda(ieig,isp0))
              end do
            end if
          end do
          nb0(isp0) = ib
        end do
        do pass0=1,2
          do isp0=1,nspin
            ib = count(lambda(:,isp0) > lambda_cut)
            do io0=1,ib
              do jo0=1,io0-1
                coeff = hvol*sum(conjg(basis(:,:,:,isp0,jo0))*basis(:,:,:,isp0,io0))
                basis(:,:,:,isp0,io0) = basis(:,:,:,isp0,io0) - &
                     basis(:,:,:,isp0,jo0)*coeff
              end do
              norm2 = hvol*sum(abs(basis(:,:,:,isp0,io0))**2)
              if (.not.ieee_is_finite(norm2) .or. norm2 <= 100d0*epsilon(1d0)) then
                status0 = 1
                exit
              end if
              basis(:,:,:,isp0,io0) = basis(:,:,:,isp0,io0)/sqrt(norm2)
            end do
          end do
        end do
      end if
      call check_collective_status(status0,"Gram-Schmidt",ik0,0)
      basis_err0 = 0d0
      if (dc%id_frag == 0) then
        do isp0=1,nspin
          do io0=1,nb0(isp0)
          do jo0=1,nb0(isp0)
            coeff = hvol*sum(conjg(basis(:,:,:,isp0,io0))*basis(:,:,:,isp0,jo0))
            if (io0 == jo0) coeff = coeff-(1d0,0d0)
            basis_err0 = max(basis_err0,abs(coeff))
          end do
          end do
        end do
      end if
      call comm_bcast(basis_err0,info%icomm_rko,0)
      status0 = 0
      if (basis_err0 > lcfo_tol .or. .not.ieee_is_finite(basis_err0)) status0 = 1
      call check_collective_status(status0,"fragment basis orthogonality",ik0,0)
      call comm_bcast(nb0,info%icomm_rko,0)
      call comm_bcast(basis,info%icomm_rko,0)
      deallocate(lambda,umat,smat,work)
    end subroutine build_basis

    subroutine collect_basis_dimensions(ik0,nb0)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: ik0,nb0(:)
      integer :: local_n(dc%n_frag,nspin),isp0,ifg0,io0,idx0
      local_n = 0
      if (dc%id_frag == 0) local_n(dc%i_frag,:) = nb0
      call comm_summation(local_n,n_basis(:,:,ik0),dc%n_frag*nspin,dc%icomm_tot)
      do isp0=1,nspin
        idx0 = 0
        do ifg0=1,dc%n_frag
          do io0=1,n_basis(ifg0,isp0,ik0)
            idx0 = idx0 + 1
            index_basis(io0,ifg0,isp0,ik0) = idx0
          end do
        end do
        n_mat(isp0,ik0) = idx0
      end do
    end subroutine collect_basis_dimensions

    subroutine assemble_halo_hamiltonian(ik0,basis,hf0,diag_h,nh,halos,nact,rs,rr)
      use communication, only: comm_irecv,comm_isend,comm_wait_all
      implicit none
      integer, intent(in) :: ik0,nh
      complex(8), intent(in) :: basis(:,:,:,:,:),hf0(:,:,:,:,:)
      complex(8), intent(out) :: diag_h(:,:,:)
      type(s_lcfo_complex_halo), intent(inout) :: halos(:)
      integer, intent(out) :: nact
      integer, allocatable, intent(out) :: rs(:),rr(:)
      integer :: h,ia,ib,ic,isp0,io0,jo0,jj0,k0,tag_send,tag_recv
      integer :: l(3),d(3),nreq

      diag_h = (0d0,0d0)
      if (dc%id_frag == 0) then
        do isp0=1,nspin
        do io0=1,n_basis(dc%i_frag,isp0,ik0)
        do jo0=1,n_basis(dc%i_frag,isp0,ik0)
          diag_h(io0,jo0,isp0) = hvol*sum(conjg(basis(:,:,:,isp0,io0))* &
               hf0(1:dc%nxyz_domain(1),1:dc%nxyz_domain(2),1:dc%nxyz_domain(3),isp0,jo0))
        end do
        end do
        end do
      end if

      nact = 0
      do h=1,nh
        if (all(halos(h)%length > 0)) nact = nact + 1
      end do
      allocate(rs(max(1,nact)),rr(max(1,nact)))
      nreq = 0
      if (dc%id_frag == 0) then
        do h=1,nh
          if (.not.all(halos(h)%length > 0)) cycle
          nreq = nreq + 1
          l = halos(h)%length
          d = halos(h)%dsp_send
          allocate(halos(h)%buf_send(l(1),l(2),l(3),nspin,m), &
               halos(h)%buf_recv(l(1),l(2),l(3),nspin,m))
          do ic=1,l(3)
          do ib=1,l(2)
          do ia=1,l(1)
            halos(h)%buf_send(ia,ib,ic,:,:) = basis(d(1)+ia,d(2)+ib,d(3)+ic,:,:)
          end do
          end do
          end do
          tag_send = halo_tag(dc%i_frag,halos(h)%dvec)
          tag_recv = halo_tag(halos(h)%ifrag_src,halos(h)%dvec)
          rr(nreq) = comm_irecv(halos(h)%buf_recv,halos(h)%id_src,tag_recv,dc%icomm_tot)
          rs(nreq) = comm_isend(halos(h)%buf_send,halos(h)%id_dst,tag_send,dc%icomm_tot)
        end do
        if (nreq > 0) then
          call comm_wait_all(rr(1:nreq))
          call comm_wait_all(rs(1:nreq))
        end if
        nreq = 0
        do h=1,nh
          if (.not.all(halos(h)%length > 0)) cycle
          nreq = nreq + 1
          l = halos(h)%length
          d = halos(h)%dsp_recv
          allocate(halos(h)%mat_h_local(m,m,nspin))
          halos(h)%mat_h_local = (0d0,0d0)
          do isp0=1,nspin
          do jo0=1,m
          do io0=1,m
            do ic=1,l(3)
            do ib=1,l(2)
            do ia=1,l(1)
              halos(h)%mat_h_local(jo0,io0,isp0) = halos(h)%mat_h_local(jo0,io0,isp0) + &
                   hvol*conjg(halos(h)%buf_recv(ia,ib,ic,isp0,jo0))* &
                   hf0(d(1)+ia,d(2)+ib,d(3)+ic,isp0,io0)
            end do
            end do
            end do
          end do
          end do
          end do
          deallocate(halos(h)%buf_send,halos(h)%buf_recv)
        end do
      end if
      deallocate(rs,rr)
    end subroutine assemble_halo_hamiltonian

    integer function halo_tag(ifg,dvec) result(tag0)
      integer, intent(in) :: ifg,dvec(3)
      tag0 = (ifg-1)*27 + (dvec(1)+1)*9 + (dvec(2)+1)*3 + dvec(3) + 2
    end function halo_tag

    subroutine deallocate_halo_buffers(nh,halos)
      implicit none
      integer, intent(in) :: nh
      type(s_lcfo_complex_halo), intent(inout) :: halos(:)
      integer :: h
      do h=1,nh
        if (allocated(halos(h)%buf_send)) deallocate(halos(h)%buf_send)
        if (allocated(halos(h)%buf_recv)) deallocate(halos(h)%buf_recv)
        if (allocated(halos(h)%mat_h_local)) deallocate(halos(h)%mat_h_local)
      end do
    end subroutine deallocate_halo_buffers

    subroutine check_eigensystem(h0,v0,e0,nt,oe,re)
      implicit none
      complex(8), intent(in) :: h0(:,:),v0(:,:)
      real(8), intent(in) :: e0(:)
      integer, intent(in) :: nt
      real(8), intent(out) :: oe,re
      complex(8), allocatable :: gram(:,:),res(:)
      integer :: row,col,k
      complex(8) :: dot
      real(8) :: denom
      allocate(gram(size(v0,2),size(v0,2)))
      gram = (0d0,0d0)
      do row=1,size(v0,2)
      do col=1,size(v0,2)
        dot = (0d0,0d0)
        do k=1,size(v0,1)
          dot = dot + conjg(v0(k,row))*v0(k,col)
        end do
        gram(row,col) = dot
      end do
      end do
      do col=1,size(gram,1)
        gram(col,col) = gram(col,col)-(1d0,0d0)
      end do
      oe = maxval(abs(gram))
      re = 0d0
      denom = max(1d0,sqrt(sum(abs(h0)**2)))
      allocate(res(size(h0,1)))
      do col=1,nt
        do row=1,size(h0,1)
          dot = (0d0,0d0)
          do k=1,size(h0,2)
            dot = dot + h0(row,k)*v0(k,col)
          end do
          res(row) = dot-e0(col)*v0(row,col)
        end do
        re = max(re,sqrt(sum(abs(res)**2))/max(denom,abs(e0(col)),1d0))
      end do
      deallocate(gram,res)
      if (.not.ieee_is_finite(oe) .or. .not.ieee_is_finite(re) .or. &
          oe > lcfo_tol .or. re > lcfo_tol) then
        ! The caller converts the failed diagnostic into a collective status.
        re = max(re,2d0*lcfo_tol)
      end if
    end subroutine check_eigensystem

    subroutine check_collective_status(local_status,stage,ik0,isp0)
      use communication, only: comm_summation
      implicit none
      integer, intent(in) :: local_status,ik0,isp0
      character(*), intent(in) :: stage
      integer :: total_status
      call comm_summation(local_status,total_status,dc%icomm_tot)
      if (total_status /= 0) then
        if (dc%id_tot == 0) write(*,'(a,2i6,1x,a)') &
             "DC-LCFO complex failure at k/spin:",ik0,isp0,trim(stage)
        stop "DC-LCFO complex: collective failure."
      end if
    end subroutine check_collective_status

    subroutine write_complex_eigenvalues(values,nbasis,nmat0)
      use communication, only: comm_summation
      use filesystem, only: get_filehandle
      implicit none
      real(8), intent(in) :: values(:,:,:)
      integer, intent(in) :: nbasis(:,:,:),nmat0(:,:)
      integer :: iu,ik0,isp0,ib0,ifg0,local_status,total_status
      character(256) :: filename
      real(8) :: qred(3)
      local_status = 0
      if (dc%id_tot == 0) then
        iu = get_filehandle()
        filename = trim(dc%base_directory)//trim(sysname)//'_lcfo_complex_eigen.data'
        open(iu,file=filename,status='replace',form='formatted',iostat=local_status)
        if (local_status == 0) then
          write(iu,'(a)') '# DC-LCFO complex scalar eigenvalues; format_version=1'
          write(iu,'(a,1x,i0,1x,a,1x,i0,1x,a,1x,i0)') '# nk=',nk,'nspin=',nspin, &
               'nstate=',dc%nstate_tot
          write(iu,'(a)') '# columns: ik ispin iband eigenvalue_Ha'
          do ik0=1,nk
            qred = [dot_product(dc%system_tot%primitive_a(:,1),dc%system_tot%vec_k(:,ik0)), &
                 dot_product(dc%system_tot%primitive_a(:,2),dc%system_tot%vec_k(:,ik0)), &
                 dot_product(dc%system_tot%primitive_a(:,3),dc%system_tot%vec_k(:,ik0))]/(2d0*pi)
            write(iu,'(a,i0,a,3(1x,es26.16e3),a,3(1x,es26.16e3),a,es26.16e3)') &
                 '# k=',ik0,' reduced=',qred,' cartesian=',dc%system_tot%vec_k(:,ik0), &
                 ' weight=',dc%system_tot%wtk(ik0)
            write(iu,'(a,i0,a,*(1x,i0))') '# k=',ik0,' nmat=',nmat0(:,ik0)
            do ifg0=1,dc%n_frag
              write(iu,'(a,i0,a,i0,a,*(1x,i0))') '# k=',ik0,' fragment=',ifg0, &
                   ' n_basis=',nbasis(ifg0,:,ik0)
            end do
            do isp0=1,nspin
              do ib0=1,dc%nstate_tot
                write(iu,'(3(i0,1x),es26.16e3)') ik0,isp0,ib0,values(ib0,isp0,ik0)
              end do
            end do
          end do
          write(iu,'(a)') '# end DC-LCFO complex eigenvalues'
          close(iu)
        end if
      end if
      call comm_summation(local_status,total_status,dc%icomm_tot)
      if (total_status /= 0) stop "DC-LCFO complex: failed to write eigenvalue output."
    end subroutine write_complex_eigenvalues

#ifdef USE_SCALAPACK
    subroutine diag_chefsi_complex_driver(ik0)
      use lcfo_diag_chefsi_complex, only: diag_chefsi_complex
      implicit none
      integer, intent(in) :: ik0
      integer, allocatable :: halo_src(:),halo_dst(:),halo_root_src(:),halo_dvec(:,:)
      complex(8), allocatable :: h_halo(:,:,:,:)
      real(8), allocatable :: ortho(:),residual(:)
      integer :: h,frag0,solve_status

      allocate(halo_src(n_halo),halo_dst(n_halo),halo_root_src(n_halo))
      allocate(halo_dvec(3,n_halo),h_halo(m,m,nspin,n_halo))
      allocate(ortho(nspin),residual(nspin))
      h_halo=(0d0,0d0)
      do h=1,n_halo
        halo_src(h)=halo(h)%ifrag_src
        halo_root_src(h)=halo(h)%id_src
        halo_dvec(:,h)=halo(h)%dvec
        halo_dst(h)=0
        do frag0=1,dc%n_frag
          if(id_array(frag0)==halo(h)%id_dst) halo_dst(h)=frag0
        end do
        if(halo_dst(h)==0) stop "DC-LCFO complex CheFSI: destination fragment not found."
        if(dc%id_frag==0 .and. allocated(halo(h)%mat_H_local)) &
          h_halo(:,:,:,h)=halo(h)%mat_H_local
      end do
      call diag_chefsi_complex(dc,ik0,nspin,lcfo_diag_chefsi_filter_degree, &
        lcfo_diag_chefsi_filter_chunk_size,lcfo_diag_chefsi_max_cycle, &
        lcfo_diag_chefsi_residual_tolerance,n_basis(:,:,ik0),n_mat(:,ik0), &
        n_halo,halo_src,halo_dst,halo_root_src,halo_dvec,mat_h_local,h_halo, &
        esp_tot(:,:,ik0),coef_frag,ortho,residual,solve_status)
      call check_collective_status(solve_status,"complex CheFSI",ik0,0)
      if(dc%id_tot==0) write(*,'(a,2i6,2(1x,es12.4))') &
        "complex CheFSI orthogonality/residual:",ik0,nspin,maxval(ortho),maxval(residual)
      deallocate(residual,ortho,h_halo,halo_dvec,halo_root_src,halo_dst,halo_src)
    end subroutine diag_chefsi_complex_driver

#endif

  end subroutine dc_lcfo_complex

  subroutine init_conventional_from_dcdft_complex(lg,mg,system,info,spsi)
    use communication, only: comm_summation,comm_bcast
    use salmon_global, only: num_fragment
    use structures, only: s_dft_system,s_orbital,s_parallel_info,s_rgrid
    implicit none
    type(s_rgrid), intent(in) :: lg,mg
    type(s_dft_system), intent(in) :: system
    type(s_parallel_info), intent(in) :: info
    type(s_orbital), intent(inout) :: spsi
    integer, allocatable :: n_basis_all(:,:,:),n_basis_wire(:),coverage(:,:,:),coverage_sum(:,:,:)
    integer :: meta(20),local_status,total_status,f,ik,ispin,io,j,ix,iy,iz
    integer :: nfrag,nspin,nk,ngrid
    real(8) :: geom(10)
    real(8), allocatable :: vec_k(:,:),wtk(:)
    character(96) :: run_id
    character(32), parameter :: bdir='./data_dcdft/fragments/'
    type(s_complex_lcfo_fragment), allocatable :: frag(:)
    complex(8), allocatable :: wrk_local(:,:,:),wrk_sum(:,:,:)

#if defined(USE_OPENACC) || defined(USE_CUDA)
    stop "DC-LCFO complex reconstruction: GPU/OpenACC/CUDA is unsupported."
#endif
    if (system%if_real_orbital) stop "DC-LCFO complex reconstruction requires complex orbitals."
    if (.not.allocated(spsi%zwf)) stop "DC-LCFO complex reconstruction: complex orbital is not allocated."
    if (.not.allocated(system%vec_k) .or. .not.allocated(system%wtk)) &
      stop "DC-LCFO complex reconstruction: k-point data is not allocated."
    if (info%isize_ro < 1 .or. info%id_ro < 0 .or. info%id_ro >= info%isize_ro) &
      stop "DC-LCFO complex reconstruction: invalid r/o communicator."
    nfrag = product(num_fragment)
    nspin = system%nspin
    nk = system%nk
    if (nfrag < 1 .or. nspin < 1 .or. nk < 1 .or. system%no < 1) &
      stop "DC-LCFO complex reconstruction: invalid system dimensions."
    if (info%ik_s < 1 .or. info%ik_e > nk .or. info%ik_s > info%ik_e) &
      stop "DC-LCFO complex reconstruction: invalid local k range."
    ngrid = product(lg%num)
    if (ngrid < 1) stop "DC-LCFO complex reconstruction: invalid total grid."
    allocate(vec_k(3,nk),wtk(nk),n_basis_all(nfrag,nspin,nk), &
         n_basis_wire(nfrag*nspin*nk))
    meta = 0
    geom = 0d0
    vec_k = 0d0
    wtk = 0d0
    run_id = ''
    n_basis_all = 0
    local_status = 0
    if (info%id_ro == 0) then
      call read_complex_lcfo_reference(system,lg,nfrag,meta,geom,vec_k,wtk,run_id, &
           n_basis_all,local_status)
    end if
    call comm_summation(local_status,total_status,info%icomm_ro)
    if (total_status /= 0) stop "DC-LCFO complex reconstruction: reference metadata is invalid."
    call comm_bcast(meta,info%icomm_ro,0)
    call comm_bcast(geom,info%icomm_ro,0)
    call comm_bcast(vec_k,info%icomm_ro,0)
    call comm_bcast(wtk,info%icomm_ro,0)
    call comm_bcast(run_id,info%icomm_ro,0)
    n_basis_wire = reshape(n_basis_all,[size(n_basis_all)])
    call comm_bcast(n_basis_wire,info%icomm_ro,0)
    n_basis_all = reshape(n_basis_wire,shape(n_basis_all))

    allocate(frag(nfrag),coverage(lg%num(1),lg%num(2),lg%num(3)))
    allocate(coverage_sum(lg%num(1),lg%num(2),lg%num(3)))
    coverage = 0
    local_status = 0
    do f=1,nfrag
      if (mod(f-1,info%isize_ro) /= info%id_ro) cycle
      call validate_complex_lcfo_fragment(f,bdir,system,lg,nfrag,meta,geom,vec_k,wtk, &
           run_id,n_basis_all,frag(f)%jxyz,coverage,local_status)
      if (local_status /= 0) exit
    end do
    call comm_summation(local_status,total_status,info%icomm_ro)
    if (total_status /= 0) stop "DC-LCFO complex reconstruction: fragment preflight failed."
    call comm_summation(coverage,coverage_sum,ngrid,info%icomm_ro)
    if (any(coverage_sum /= 1)) then
      if (info%id_ro == 0) write(*,*) &
           "DC-LCFO complex reconstruction: fragment cores do not cover the total grid exactly once."
      stop "DC-LCFO complex reconstruction: invalid rgrid coverage."
    end if

    allocate(wrk_local(lg%num(1),lg%num(2),lg%num(3)))
    allocate(wrk_sum(lg%num(1),lg%num(2),lg%num(3)))
    if (info%id_ro == 0) then
      write(*,*) "start complex DC-LCFO wavefunction reconstruction"
      write(*,*) "complex LCFO format v1; fragments/k/spins:",nfrag,nk,nspin
    end if
    do ik=info%ik_s,info%ik_e
      do f=1,nfrag
        if (mod(f-1,info%isize_ro) /= info%id_ro) cycle
        call load_complex_lcfo_fragment_k(f,ik,bdir,system,lg,nfrag,meta,geom,vec_k,wtk, &
             run_id,n_basis_all,frag(f),local_status)
        if (local_status /= 0) exit
      end do
      call comm_summation(local_status,total_status,info%icomm_ro)
      if (total_status /= 0) stop "DC-LCFO complex reconstruction: k-record read failed."
      do ispin=1,nspin
        do io=1,system%no
          wrk_local = (0d0,0d0)
          do f=1,nfrag
            if (.not.allocated(frag(f)%basis)) cycle
            do j=1,frag(f)%n_basis(ispin)
              do iz=1,meta(9)
              do iy=1,meta(8)
              do ix=1,meta(7)
                wrk_local(frag(f)%jxyz(ix,1),frag(f)%jxyz(iy,2),frag(f)%jxyz(iz,3)) = &
                     wrk_local(frag(f)%jxyz(ix,1),frag(f)%jxyz(iy,2),frag(f)%jxyz(iz,3)) + &
                     frag(f)%basis(ix,iy,iz,ispin,j)*frag(f)%coef(j,io,ispin)
              end do
              end do
              end do
            end do
          end do
          call comm_summation(wrk_local,wrk_sum,ngrid,info%icomm_ro)
          if (info%io_s <= io .and. io <= info%io_e) then
            do iz=mg%is(3),mg%ie(3)
            do iy=mg%is(2),mg%ie(2)
            do ix=mg%is(1),mg%ie(1)
              spsi%zwf(ix,iy,iz,ispin,io,ik,1) = wrk_sum(ix,iy,iz)
            end do
            end do
            end do
          end if
        end do
      end do
      do f=1,nfrag
        if (allocated(frag(f)%basis)) deallocate(frag(f)%basis,frag(f)%coef,frag(f)%n_basis)
      end do
    end do
    if (info%id_ro == 0) write(*,*) "end complex DC-LCFO wavefunction reconstruction"
    deallocate(wrk_sum,wrk_local,coverage_sum,coverage,frag,n_basis_wire,n_basis_all,wtk,vec_k)
  end subroutine init_conventional_from_dcdft_complex

  subroutine open_complex_lcfo_read(path,file_kind,system,lg,unit,meta,geom,vec_k,wtk,run_id,file_size,status)
    use filesystem, only: get_filehandle
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: int32,int64,real64,file_storage_size
    use structures, only: s_dft_system,s_rgrid
    implicit none
    character(*), intent(in) :: path
    integer, intent(in) :: file_kind
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg
    integer, intent(out) :: unit,status
    integer, intent(out) :: meta(20)
    real(real64), intent(out) :: geom(10),vec_k(:,:),wtk(:)
    character(96), intent(out) :: run_id
    integer(int64), intent(out) :: file_size
    character(16) :: magic
    integer(int32) :: control(6)
    integer(int64) :: header_bytes
    integer :: ios

    status = 1
    unit = -1
    file_size = 0_int64
    if (storage_size(0) /= 32 .or. storage_size(0d0) /= 64 .or. &
        storage_size(cmplx(0d0,0d0,8)) /= 128 .or. file_storage_size /= 8) return
    unit = get_filehandle()
    open(unit,file=trim(path),status='old',form='unformatted', &
         access='stream',action='read',iostat=ios)
    if (ios /= 0) then
      unit = -1
      return
    end if
    read(unit,iostat=ios) magic,control,header_bytes,run_id,meta,geom,vec_k,wtk
    if (ios /= 0) goto 900
    if (magic /= 'SLCFO_COMPLEX_V1' .or. any(control /= &
        [1_int32,int(z'01020304',int32),32_int32,64_int32,int(file_kind,int32),0_int32])) goto 900
    if (header_bytes /= 304_int64+32_int64*int(system%nk,int64)) goto 900
    if (size(vec_k,1) /= 3 .or. size(vec_k,2) /= system%nk .or. size(wtk) /= system%nk) goto 900
    if (any(meta(1:3) /= lg%num) .or. any(meta(10:12) < 1) .or. &
        meta(16) < 1 .or. meta(17) /= system%nspin .or. meta(18) /= system%nk .or. &
        meta(19) < 1 .or. meta(20) < system%no) goto 900
    if (any(meta(7:9) < 1) .or. any(meta(7:9) > meta(4:6))) goto 900
    if (maxval(abs(geom(1:9)-reshape(system%primitive_a,[9]))) > &
        1d-12*max(1d0,maxval(abs(geom(1:9))))) goto 900
    if (abs(geom(10)-system%hvol) > 1d-12*max(1d0,abs(geom(10)))) goto 900
    if (any(abs(vec_k-system%vec_k) > 1d-12*max(1d0,maxval(abs(vec_k)))) .or. &
        any(abs(wtk-system%wtk) > 1d-12*max(1d0,maxval(abs(wtk))))) goto 900
    if (any(.not.ieee_is_finite(geom)) .or. any(.not.ieee_is_finite(vec_k)) .or. &
        any(.not.ieee_is_finite(wtk))) goto 900
    inquire(unit=unit,size=file_size,iostat=ios)
    if (ios /= 0 .or. file_size < header_bytes+124_int64) goto 900
    status = 0
    return
900 continue
    close(unit,iostat=ios)
    unit = -1
  end subroutine open_complex_lcfo_read

  subroutine skip_complex_lcfo_reals(unit,nreal,status)
    use iso_fortran_env, only: real64,int64
    implicit none
    integer, intent(in) :: unit
    integer(int64), intent(in) :: nreal
    integer, intent(inout) :: status
    real(real64) :: buffer(8192)
    integer(int64) :: remain
    integer :: n,ios
    if (status /= 0) return
    if (nreal < 0_int64) then
      status = 1
      return
    end if
    remain = nreal
    do while(remain > 0_int64)
      n = int(min(remain,int(size(buffer),int64)))
      read(unit,iostat=ios) buffer(1:n)
      if (ios /= 0) then
        status = 1
        return
      end if
      remain = remain-int(n,int64)
    end do
  end subroutine skip_complex_lcfo_reals

  subroutine skip_complex_lcfo_bytes(unit,nbytes,status)
    use iso_fortran_env, only: int8,int64
    implicit none
    integer, intent(in) :: unit
    integer(int64), intent(in) :: nbytes
    integer, intent(inout) :: status
    integer(int8) :: buffer(65536)
    integer(int64) :: remain
    integer :: n,ios
    if (status /= 0) return
    if (nbytes < 0_int64) then
      status=1
      return
    end if
    remain=nbytes
    do while(remain > 0_int64)
      n=int(min(remain,int(size(buffer),int64)))
      read(unit,iostat=ios) buffer(1:n)
      if (ios /= 0) then
        status=1
        return
      end if
      remain=remain-int(n,int64)
    end do
  end subroutine skip_complex_lcfo_bytes

  subroutine read_complex_lcfo_reference(system,lg,nfrag,meta,geom,vec_k,wtk,run_id,n_basis_all,status)
    use iso_fortran_env, only: int32,int64,real64
    use salmon_global, only: num_fragment
    use structures, only: s_dft_system,s_rgrid
    implicit none
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg
    integer, intent(in) :: nfrag
    integer, intent(out) :: meta(20),n_basis_all(:,:,:),status
    real(real64), intent(out) :: geom(10),vec_k(:,:),wtk(:)
    character(96), intent(out) :: run_id
    character(32), parameter :: bdir='./data_dcdft/fragments/'
    character(256) :: filename
    integer :: iu_basis,iu_coeff,ios,fcount,ik,ispin,nb,nmat,j,prefix
    integer(int32) :: ik_wire,spin_wire,nb_wire,nmat_wire,counts(size(n_basis_all,1)),row
    integer(int32) :: footer_k
    integer(int64) :: payload,payload_start,payload_end,file_size,footer_size
    character(16) :: footer_magic
    character(96) :: footer_run
    integer :: meta_basis(20)
    real(real64) :: geom_basis(10),vec_basis(3,size(vec_k,2)),wt_basis(size(wtk))
    character(96) :: run_basis
    integer :: row_index

    status = 1
    meta = 0
    geom = 0d0
    vec_k = 0d0
    wtk = 0d0
    run_id = ''
    n_basis_all = 0
    fcount = product(num_fragment)
    if (fcount /= nfrag .or. nfrag < 1) then
      write(*,*) "DC-LCFO complex reconstruction: input/file fragment count mismatch:",fcount,nfrag
      return
    end if
    write(filename,'(a,i6.6,a,a)') trim(bdir),1,'/','basis_functions.bin'
    call open_complex_lcfo_read(filename,1,system,lg,iu_basis,meta,geom,vec_k,wtk,run_id,file_size,status)
    if (status /= 0) then
      write(*,*) "DC-LCFO complex reconstruction: invalid reference basis header:",trim(filename)
      return
    end if
    close(iu_basis,iostat=ios)
    if (ios /= 0) then
      write(*,*) "DC-LCFO complex reconstruction: failed closing reference basis file."
      status = 1
      return
    end if
    write(filename,'(a,i6.6,a,a)') trim(bdir),1,'/','wavefunctions.bin'
    call open_complex_lcfo_read(filename,2,system,lg,iu_coeff,meta_basis,geom_basis, &
         vec_basis,wt_basis,run_basis,file_size,status)
    if (status /= 0) then
      write(*,*) "DC-LCFO complex reconstruction: invalid reference coefficient header:",trim(filename)
      return
    end if
    if (any(meta_basis(1:15) /= meta(1:15)) .or. any(meta_basis(16:20) /= meta(16:20)) .or. &
        run_basis /= run_id .or. any(abs(geom_basis-geom) > 0d0) .or. &
        any(abs(vec_basis-vec_k) > 0d0) .or. any(abs(wt_basis-wtk) > 0d0)) then
      close(iu_coeff)
      write(*,*) "DC-LCFO complex reconstruction: inconsistent reference basis/coefficient headers."
      status = 1
      return
    end if
    do ik=1,system%nk
      read(iu_coeff,iostat=ios) ik_wire,payload
      if (ios /= 0) exit
      inquire(unit=iu_coeff,pos=payload_start,iostat=ios)
      if (ios /= 0 .or. ik_wire /= ik .or. payload < 0_int64) exit
      payload_end = payload_start
      do ispin=1,system%nspin
        read(iu_coeff,iostat=ios) spin_wire,nb_wire,nmat_wire
        if (ios /= 0) exit
        if (spin_wire /= ispin .or. nb_wire < 0 .or. nb_wire > meta(19)) then
          ios = 1
          exit
        end if
        read(iu_coeff,iostat=ios) counts
        if (ios /= 0) exit
        if (any(counts < 0) .or. any(counts > meta(19)) .or. &
            int(sum(counts),int64) /= int(nmat_wire,int64) .or. counts(1) /= nb_wire) then
          ios = 1
          exit
        end if
        n_basis_all(:,ispin,ik) = counts
        prefix = 0
        do j=1,int(nb_wire)
          read(iu_coeff,iostat=ios) row
          if (ios /= 0) exit
          row_index = prefix+j
          if (row /= row_index) ios = 1
        end do
        if (ios /= 0) exit
        call skip_complex_lcfo_reals(iu_coeff,2_int64*int(nb_wire,int64)* &
             int(meta(20),int64),ios)
        if (ios /= 0) exit
      end do
      if (ios /= 0) exit
      inquire(unit=iu_coeff,pos=payload_end,iostat=ios)
      if (ios /= 0 .or. payload_end-payload_start /= payload) then
        write(*,*) "DC-LCFO complex reconstruction: coefficient payload length/IO status:", &
             ik,payload,payload_end-payload_start,ios
        ios = 1
        exit
      end if
    end do
    if (ios /= 0) write(*,*) &
         "DC-LCFO complex reconstruction: invalid reference coefficient record k/spin:",ik,ispin
    if (ios == 0) then
      read(iu_coeff,iostat=ios) footer_magic,footer_run,footer_k,footer_size
      if (ios == 0) then
        inquire(unit=iu_coeff,pos=payload_end,iostat=ios)
        if (ios == 0) then
          if (footer_magic /= 'SLCFO_DONE_V1' .or. footer_run /= run_id .or. &
              footer_k /= system%nk .or. footer_size /= file_size .or. payload_end-1_int64 /= file_size) ios = 1
        end if
      end if
    end if
    if (ios /= 0) write(*,*) "DC-LCFO complex reconstruction: invalid reference coefficient footer."
    close(iu_coeff,iostat=fcount)
    if (ios /= 0 .or. fcount /= 0) then
      status = 1
      return
    end if
    if (any(n_basis_all < 0)) then
      status = 1
      return
    end if
    status = 0
  end subroutine read_complex_lcfo_reference

  subroutine validate_complex_lcfo_fragment(f,bdir,system,lg,nfrag,ref_meta,ref_geom,ref_k, &
       ref_wtk,ref_run,n_basis_all,jxyz,coverage,status)
    use filesystem, only: get_filehandle
    use iso_fortran_env, only: int32,int64,real64
    use structures, only: s_dft_system,s_rgrid
    implicit none
    integer, intent(in) :: f,nfrag,ref_meta(20),n_basis_all(:,:,:)
    character(*), intent(in) :: bdir
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg
    real(real64), intent(in) :: ref_geom(10),ref_k(:,:),ref_wtk(:)
    character(96), intent(in) :: ref_run
    integer, allocatable, intent(out) :: jxyz(:,:)
    integer, intent(inout) :: coverage(:,:,:)
    integer, intent(out) :: status
    character(256) :: filename
    character(16) :: footer_basis,footer_coeff
    character(96) :: run_basis,run_coeff
    integer :: ub,uc,ur,ios,close_ios,ik,ispin,d,nb,nmat,j,prefix,ix,iy,iz
    integer :: meta_b(20),meta_c(20),gfrag(3),gtot(3),counts(size(n_basis_all,1))
    integer(int32) :: ikb,ikc,spin_wire,nb_wire,nmat_wire,footer_kb,footer_kc,row
    integer(int64) :: size_b,size_c,payload_b,payload_c,start_b,start_c,end_b,end_c
    integer(int64) :: nreal
    real(real64) :: geom_b(10),geom_c(10),vk_b(3,size(ref_k,2)),vk_c(3,size(ref_k,2))
    real(real64) :: w_b(size(ref_wtk)),w_c(size(ref_wtk))
    integer(int64) :: footer_bytes_b,footer_bytes_c

    status = 1
    write(filename,'(a,i6.6,2a)') trim(bdir),f,'/','basis_functions.bin'
    call open_complex_lcfo_read(filename,1,system,lg,ub,meta_b,geom_b,vk_b,w_b, &
         run_basis,size_b,ios)
    if (ios /= 0) return
    write(filename,'(a,i6.6,2a)') trim(bdir),f,'/','wavefunctions.bin'
    call open_complex_lcfo_read(filename,2,system,lg,uc,meta_c,geom_c,vk_c,w_c, &
         run_coeff,size_c,ios)
    if (ios /= 0) then
      close(ub)
      return
    end if
    if (.not.complex_lcfo_headers_match(meta_b,geom_b,vk_b,w_b,run_basis, &
        meta_c,geom_c,vk_c,w_c,run_coeff,ref_meta,ref_geom,ref_k,ref_wtk,ref_run)) goto 900
    if (meta_b(16) /= f .or. any(meta_b(4:6) /= ref_meta(4:6))) goto 900
    write(filename,'(a,i6.6,2a)') trim(bdir),f,'/','rgrid_index.bin'
    ur = get_filehandle()
    open(ur,file=trim(filename),status='old',form='unformatted',access='stream', &
         action='read',iostat=ios)
    if (ios /= 0) goto 900
    read(ur,iostat=ios) gfrag,gtot
    if (ios /= 0 .or. any(gfrag /= meta_b(4:6)) .or. any(gtot /= lg%num) .or. &
        any(meta_b(7:9) > gfrag)) then
      ios = 1
      goto 890
    end if
    allocate(jxyz(maxval(gfrag),3))
    jxyz = 0
    do d=1,3
      read(ur,iostat=ios) jxyz(1:gfrag(d),d)
      if (ios /= 0) exit
      if (any(jxyz(1:gfrag(d),d) < 1) .or. any(jxyz(1:gfrag(d),d) > lg%num(d))) then
        ios = 1
        exit
      end if
      if (jxyz(1,d) /= meta_b(12+d)) then
        ios = 1
        exit
      end if
    end do
890 continue
    close(ur,iostat=close_ios)
    if (ios /= 0 .or. close_ios /= 0) goto 900

    do ik=1,system%nk
      read(ub,iostat=ios) ikb,payload_b
      if (ios /= 0) exit
      inquire(unit=ub,pos=start_b,iostat=ios)
      if (ios /= 0 .or. ikb /= ik .or. payload_b < 0_int64) exit
      do ispin=1,system%nspin
        read(ub,iostat=ios) spin_wire,nb_wire
        if (ios /= 0) exit
        if (spin_wire /= ispin .or. nb_wire /= n_basis_all(f,ispin,ik)) then
          ios = 1
          exit
        end if
        nreal = 2_int64
        do d=7,9
          if (int(meta_b(d),int64) > huge(nreal)/nreal) then
            ios = 1
            exit
          end if
          nreal = nreal*int(meta_b(d),int64)
        end do
        if (ios /= 0) exit
        if (int(nb_wire,int64) > huge(nreal)/max(1_int64,nreal)) then
          ios = 1
          exit
        end if
        nreal = nreal*int(nb_wire,int64)
        call skip_complex_lcfo_reals(ub,nreal,ios)
        if (ios /= 0) exit
      end do
      if (ios /= 0) exit
      inquire(unit=ub,pos=end_b,iostat=ios)
      if (ios /= 0 .or. end_b-start_b /= payload_b) then
        ios = 1
        exit
      end if

      read(uc,iostat=ios) ikc,payload_c
      if (ios /= 0) exit
      inquire(unit=uc,pos=start_c,iostat=ios)
      if (ios /= 0 .or. ikc /= ik .or. payload_c < 0_int64) exit
      do ispin=1,system%nspin
        read(uc,iostat=ios) spin_wire,nb_wire,nmat_wire
        if (ios /= 0) exit
        read(uc,iostat=ios) counts
        if (ios /= 0) exit
        if (spin_wire /= ispin .or. any(counts /= n_basis_all(:,ispin,ik)) .or. &
            nb_wire /= counts(f) .or. nmat_wire /= sum(counts) .or. &
            nb_wire < 0 .or. nb_wire > meta_c(19)) then
          ios = 1
          exit
        end if
        prefix = sum(counts(1:f-1))
        do j=1,nb_wire
          read(uc,iostat=ios) row
          if (ios /= 0) exit
          if (row /= prefix+j) then
            ios = 1
            exit
          end if
        end do
        if (ios /= 0) exit
        call skip_complex_lcfo_reals(uc,2_int64*int(nb_wire,int64)* &
             int(meta_c(20),int64),ios)
        if (ios /= 0) exit
      end do
      if (ios /= 0) exit
      inquire(unit=uc,pos=end_c,iostat=ios)
      if (ios /= 0 .or. end_c-start_c /= payload_c) then
        ios = 1
        exit
      end if
    end do
    if (ios /= 0) goto 900
    read(ub,iostat=ios) footer_basis,run_basis,footer_kb,footer_bytes_b
    if (ios /= 0) goto 900
    inquire(unit=ub,pos=end_b,iostat=ios)
    if (ios /= 0 .or. footer_basis /= 'SLCFO_DONE_V1' .or. run_basis /= ref_run .or. &
        footer_kb /= system%nk .or. footer_bytes_b /= size_b .or. end_b-1_int64 /= size_b) goto 900
    read(uc,iostat=ios) footer_coeff,run_coeff,footer_kc,footer_bytes_c
    if (ios /= 0) goto 900
    inquire(unit=uc,pos=end_c,iostat=ios)
    if (ios /= 0 .or. footer_coeff /= 'SLCFO_DONE_V1' .or. run_coeff /= ref_run .or. &
        footer_kc /= system%nk .or. footer_bytes_c /= size_c .or. end_c-1_int64 /= size_c) goto 900
    do iz=1,meta_b(9)
    do iy=1,meta_b(8)
    do ix=1,meta_b(7)
      coverage(jxyz(ix,1),jxyz(iy,2),jxyz(iz,3)) = &
           coverage(jxyz(ix,1),jxyz(iy,2),jxyz(iz,3))+1
    end do
    end do
    end do
    status = 0
900 continue
    close(ub,iostat=close_ios)
    if (close_ios /= 0) status = 1
    close(uc,iostat=close_ios)
    if (close_ios /= 0) status = 1
    if (status /= 0 .and. allocated(jxyz)) deallocate(jxyz)
  end subroutine validate_complex_lcfo_fragment

  logical function complex_lcfo_headers_match(meta_a,geom_a,k_a,wtk_a,run_a, &
       meta_b,geom_b,k_b,wtk_b,run_b,ref_meta,ref_geom,ref_k,ref_wtk,ref_run)
    use iso_fortran_env, only: real64
    implicit none
    integer, intent(in) :: meta_a(20),meta_b(20),ref_meta(20)
    real(real64), intent(in) :: geom_a(10),geom_b(10),k_a(:,:),k_b(:,:),wtk_a(:),wtk_b(:)
    real(real64), intent(in) :: ref_geom(10),ref_k(:,:),ref_wtk(:)
    character(96), intent(in) :: run_a,run_b,ref_run
    integer :: ma(20),mb(20),mr(20)
    real(real64) :: scale
    ma=meta_a
    mb=meta_b
    mr=ref_meta
    ma(12:16)=0
    mb(12:16)=0
    mr(12:16)=0
    complex_lcfo_headers_match = .false.
    if (any(ma /= mr) .or. any(mb /= mr)) return
    if (run_a /= run_b .or. run_a /= ref_run) return
    scale=max(1d0,maxval(abs(ref_geom)))
    if (maxval(abs(geom_a-ref_geom)) > 1d-12*scale .or. &
        maxval(abs(geom_b-ref_geom)) > 1d-12*scale) return
    scale=max(1d0,maxval(abs(ref_k)))
    if (maxval(abs(k_a-ref_k)) > 1d-12*scale .or. maxval(abs(k_b-ref_k)) > 1d-12*scale) return
    scale=max(1d0,maxval(abs(ref_wtk)))
    if (maxval(abs(wtk_a-ref_wtk)) > 1d-12*scale .or. maxval(abs(wtk_b-ref_wtk)) > 1d-12*scale) return
    complex_lcfo_headers_match = .true.
  end function complex_lcfo_headers_match

  subroutine read_wire_complex_values(unit,values,status)
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: real64
    implicit none
    integer, intent(in) :: unit
    complex(8), intent(out) :: values(:)
    integer, intent(inout) :: status
    integer, parameter :: chunk=4096
    real(real64) :: wire(2*chunk)
    integer :: first,n,j,ios
    if (status /= 0) return
    first=1
    do while(first <= size(values))
      n=min(chunk,size(values)-first+1)
      read(unit,iostat=ios) wire(1:2*n)
      if (ios /= 0) then
        status=1
        return
      end if
      if (any(.not.ieee_is_finite(wire(1:2*n)))) then
        status=1
        return
      end if
      do j=1,n
        values(first+j-1)=cmplx(wire(2*j-1),wire(2*j),kind=8)
      end do
      first=first+n
    end do
  end subroutine read_wire_complex_values

  subroutine load_complex_lcfo_fragment_k(f,ik,bdir,system,lg,nfrag,ref_meta,ref_geom, &
       ref_k,ref_wtk,ref_run,n_basis_all,fragment,status)
    use iso_fortran_env, only: int32,int64,real64
    use structures, only: s_dft_system,s_rgrid
    implicit none
    integer, intent(in) :: f,ik,nfrag,ref_meta(20),n_basis_all(:,:,:)
    character(*), intent(in) :: bdir
    type(s_dft_system), intent(in) :: system
    type(s_rgrid), intent(in) :: lg
    real(real64), intent(in) :: ref_geom(10),ref_k(:,:),ref_wtk(:)
    character(96), intent(in) :: ref_run
    type(s_complex_lcfo_fragment), intent(inout) :: fragment
    integer, intent(out) :: status
    character(256) :: filename
    character(96) :: run_b,run_c
    integer :: ub,uc,ios,ispin,nb,ik0,jrow,prefix,meta_b(20),meta_c(20),stage
    integer(int32) :: ikwire,spinwire,nbwire,nmatwire,counts(size(n_basis_all,1)),row
    integer(int64) :: size_b,size_c,payload_b,payload_c,nreal
    real(real64) :: geom_b(10),geom_c(10),vk_b(3,size(ref_k,2)),vk_c(3,size(ref_k,2))
    real(real64) :: w_b(size(ref_wtk)),w_c(size(ref_wtk))
    complex(8), allocatable :: values(:),matrix_tmp(:,:)

    status=1
    stage=1
    if (ik < 1 .or. ik > system%nk .or. .not.allocated(fragment%jxyz)) return
    write(filename,'(a,i6.6,2a)') trim(bdir),f,'/','basis_functions.bin'
    call open_complex_lcfo_read(filename,1,system,lg,ub,meta_b,geom_b,vk_b,w_b,run_b,size_b,ios)
    if (ios /= 0) return
    write(filename,'(a,i6.6,2a)') trim(bdir),f,'/','wavefunctions.bin'
    call open_complex_lcfo_read(filename,2,system,lg,uc,meta_c,geom_c,vk_c,w_c,run_c,size_c,ios)
    if (ios /= 0) then
      close(ub)
      return
    end if
    if (.not.complex_lcfo_headers_match(meta_b,geom_b,vk_b,w_b,run_b, &
        meta_c,geom_c,vk_c,w_c,run_c,ref_meta,ref_geom,ref_k,ref_wtk,ref_run)) goto 900
    if (meta_b(16) /= f) goto 900
    stage=2
    do ik0=1,ik
      stage=3
      read(ub,iostat=ios) ikwire,payload_b
      if (ios /= 0 .or. ikwire /= ik0 .or. payload_b < 0_int64) goto 900
      if (ik0 < ik) then
        call skip_complex_lcfo_bytes(ub,payload_b,ios)
        if (ios /= 0) goto 900
      else
        allocate(fragment%n_basis(system%nspin))
        allocate(fragment%basis(meta_b(7),meta_b(8),meta_b(9),system%nspin,meta_b(19)))
        fragment%n_basis=0
        fragment%basis=(0d0,0d0)
        do ispin=1,system%nspin
          read(ub,iostat=ios) spinwire,nbwire
          if (ios /= 0 .or. spinwire /= ispin .or. nbwire /= n_basis_all(f,ispin,ik)) goto 900
          nb=int(nbwire)
          fragment%n_basis(ispin)=nb
          nreal=2_int64*int(meta_b(7),int64)*int(meta_b(8),int64)*int(meta_b(9),int64)*int(nb,int64)
          if (nb > 0) then
            if (nreal/2_int64 > int(huge(1),int64)) goto 900
            allocate(values(int(nreal/2_int64)))
            call read_wire_complex_values(ub,values,ios)
            if (ios /= 0) goto 900
            fragment%basis(:,:,:,ispin,1:nb)=reshape(values, &
                 [meta_b(7),meta_b(8),meta_b(9),nb])
            deallocate(values)
          end if
        end do
      end if
      stage=4
      read(uc,iostat=ios) ikwire,payload_c
      if (ios /= 0 .or. ikwire /= ik0 .or. payload_c < 0_int64) goto 900
      if (ik0 < ik) then
        call skip_complex_lcfo_bytes(uc,payload_c,ios)
        if (ios /= 0) goto 900
      else
        allocate(fragment%coef(meta_c(19),meta_c(20),system%nspin))
        fragment%coef=(0d0,0d0)
        do ispin=1,system%nspin
          stage=5
          read(uc,iostat=ios) spinwire,nbwire,nmatwire
          if (ios /= 0) goto 900
          read(uc,iostat=ios) counts
          if (ios /= 0) goto 900
          nb=int(nbwire)
          if (spinwire /= ispin .or. nb /= n_basis_all(f,ispin,ik) .or. &
              any(counts /= n_basis_all(:,ispin,ik)) .or. nmatwire /= sum(counts)) goto 900
          prefix=sum(counts(1:f-1))
          do jrow=1,nb
            read(uc,iostat=ios) row
            if (ios /= 0 .or. row /= prefix+jrow) goto 900
          end do
          if (nb > 0) then
            if (int(nb,int64)*int(meta_c(20),int64) > int(huge(1),int64)) goto 900
            allocate(values(nb*meta_c(20)),matrix_tmp(nb,meta_c(20)))
            call read_wire_complex_values(uc,values,ios)
            if (ios /= 0) goto 900
            matrix_tmp=reshape(values,[nb,meta_c(20)])
            fragment%coef(1:nb,1:system%no,ispin)=matrix_tmp(:,1:system%no)
            deallocate(matrix_tmp,values)
          end if
        end do
      end if
    end do
    status=0
900 continue
    if (status /= 0) write(*,*) &
         "DC-LCFO complex reconstruction: k record load failed (fragment,k,stage,io):",f,ik,stage,ios
    close(ub,iostat=ios)
    if (ios /= 0) status=1
    close(uc,iostat=ios)
    if (ios /= 0) status=1
    if (status /= 0) then
      if (allocated(fragment%n_basis)) deallocate(fragment%n_basis)
      if (allocated(fragment%basis)) deallocate(fragment%basis)
      if (allocated(fragment%coef)) deallocate(fragment%coef)
    end if
  end subroutine load_complex_lcfo_fragment_k

end module lcfo_complex
