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

module dc_lcfo_flux_operator
  implicit none
  private

  public :: apply_dc_lcfo_flux_hpsi_rwf, apply_dc_lcfo_flux_hpsi_zwf

  type halo_cg_info
    integer :: id_src,id_dst,ifrag_src,ifrag_dst,dvec(3),length(3),dsp_send(3),axis
  end type halo_cg_info

  type(halo_cg_info), allocatable, target, save :: flux_halo_cache(:)
  integer, save :: flux_halo_cache_n = -1
  integer, save :: flux_halo_cache_frag = -1
  integer, save :: flux_halo_cache_nfrag = -1
  integer, save :: flux_halo_cache_isize = -1
  integer, save :: flux_halo_cache_ncore(3) = -1
  integer, save :: flux_halo_cache_nbuf(3) = -1
  integer, save :: flux_halo_cache_radius(3) = -1
  integer, save :: flux_halo_cache_num_fragment(3) = -1

contains

subroutine prepare_dc_lcfo_flux_halo_cache(dc,ncore,stencil_radius,npg_base)
  use structures
  use salmon_global, only: num_fragment
  implicit none
  type(s_dcdft),intent(in) :: dc
  integer,intent(in) :: ncore(3),stencil_radius(3),npg_base
  integer :: nh(3),ir1(3),ir2(3),d(3)
  integer :: lx,ly,lz,ifrag,n,idir,nonzero_dirs,rank_same_subgroup
  logical :: cache_valid

  cache_valid = allocated(flux_halo_cache)
  cache_valid = cache_valid .and. flux_halo_cache_frag == dc%i_frag
  cache_valid = cache_valid .and. flux_halo_cache_nfrag == dc%n_frag
  cache_valid = cache_valid .and. flux_halo_cache_isize == dc%isize_tot
  cache_valid = cache_valid .and. all(flux_halo_cache_ncore == ncore)
  cache_valid = cache_valid .and. all(flux_halo_cache_nbuf == dc%nxyz_buffer)
  cache_valid = cache_valid .and. all(flux_halo_cache_radius == stencil_radius)
  cache_valid = cache_valid .and. all(flux_halo_cache_num_fragment == num_fragment)
  if(cache_valid) return

  nh = 0
  do n=1,3
    if(dc%nxyz_buffer(n) > ncore(n)) stop "DC-LCFO-Flux CG: buffer > domain"
    if(num_fragment(n) > 1 .and. dc%nxyz_buffer(n) < stencil_radius(n)) &
    & stop "DC-LCFO-Flux CG: buffer is smaller than the active stencil radius"
    if(num_fragment(n) > 1) nh(n) = 1
  end do

  if(allocated(flux_halo_cache)) deallocate(flux_halo_cache)
  allocate(flux_halo_cache(26))
  flux_halo_cache_n = 0
  do lx=-nh(1),nh(1)
  do ly=-nh(2),nh(2)
  do lz=-nh(3),nh(3)
    if(lx==0 .and. ly==0 .and. lz==0) cycle
    nonzero_dirs = count([lx,ly,lz] /= 0)
    if(nonzero_dirs /= 1) cycle

    flux_halo_cache_n = flux_halo_cache_n + 1
    flux_halo_cache(flux_halo_cache_n)%dvec = [lx,ly,lz]
    flux_halo_cache(flux_halo_cache_n)%axis = 0
    do idir=1,3
      if(flux_halo_cache(flux_halo_cache_n)%dvec(idir) /= 0) &
      & flux_halo_cache(flux_halo_cache_n)%axis = idir
    end do
    flux_halo_cache(flux_halo_cache_n)%id_dst = -1
    flux_halo_cache(flux_halo_cache_n)%id_src = -1
    flux_halo_cache(flux_halo_cache_n)%ifrag_dst = -1
    flux_halo_cache(flux_halo_cache_n)%ifrag_src = -1
    do ifrag=1,dc%n_frag
      ir1 = dc%ixyz_frag(:,ifrag)
      ir2 = dc%ixyz_frag(:,dc%i_frag) + flux_halo_cache(flux_halo_cache_n)%dvec*ncore
      d = mod(ir1 - ir2, dc%lg_tot%num(1:3))
      if(all(d == 0) .and. flux_halo_cache(flux_halo_cache_n)%id_dst < 0) then
        rank_same_subgroup = -1
        if(dc%id_frag < npg_base) rank_same_subgroup = (ifrag-1)*npg_base + dc%id_frag
        flux_halo_cache(flux_halo_cache_n)%id_dst = rank_same_subgroup
        flux_halo_cache(flux_halo_cache_n)%ifrag_dst = ifrag
      end if

      ir2 = dc%ixyz_frag(:,dc%i_frag) - flux_halo_cache(flux_halo_cache_n)%dvec*ncore
      d = mod(ir1 - ir2, dc%lg_tot%num(1:3))
      if(all(d == 0) .and. flux_halo_cache(flux_halo_cache_n)%id_src < 0) then
        rank_same_subgroup = -1
        if(dc%id_frag < npg_base) rank_same_subgroup = (ifrag-1)*npg_base + dc%id_frag
        flux_halo_cache(flux_halo_cache_n)%id_src = rank_same_subgroup
        flux_halo_cache(flux_halo_cache_n)%ifrag_src = ifrag
      end if
    end do
    if(flux_halo_cache(flux_halo_cache_n)%id_dst < 0 .or. &
    &  flux_halo_cache(flux_halo_cache_n)%id_src < 0) &
    & stop "DC-LCFO-Flux CG: missing halo rank"

    do n=1,3
      select case(flux_halo_cache(flux_halo_cache_n)%dvec(n))
      case(0)
        flux_halo_cache(flux_halo_cache_n)%length(n) = ncore(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = 0
      case(1)
        flux_halo_cache(flux_halo_cache_n)%length(n) = dc%nxyz_buffer(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = ncore(n) - dc%nxyz_buffer(n)
      case(-1)
        flux_halo_cache(flux_halo_cache_n)%length(n) = dc%nxyz_buffer(n)
        flux_halo_cache(flux_halo_cache_n)%dsp_send(n) = 0
      end select
    end do
  end do
  end do
  end do

  flux_halo_cache_frag = dc%i_frag
  flux_halo_cache_nfrag = dc%n_frag
  flux_halo_cache_isize = dc%isize_tot
  flux_halo_cache_ncore = ncore
  flux_halo_cache_nbuf = dc%nxyz_buffer
  flux_halo_cache_radius = stencil_radius
  flux_halo_cache_num_fragment = num_fragment
end subroutine prepare_dc_lcfo_flux_halo_cache

subroutine apply_dc_lcfo_flux_hpsi_rwf(mg,system,info,stencil,dc,psi,hpsi)
  use structures
  use communication, only: comm_isend, comm_irecv, comm_wait_all, comm_proc_null
  use dc_fragment_geometry, only: get_fragment_domain
  implicit none
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_dcdft),          intent(in)    :: dc
  type(s_orbital),        intent(in)    :: psi
  type(s_orbital),        intent(inout) :: hpsi
  !
  type(halo_cg_info), pointer :: halo(:)
  integer :: ncore(3),l(3),lo(3),hi(3)
  integer :: i_halo,n_halo,n,axis,dist
  integer :: ix,iy,iz,ixg,iyg,izg,il(3),ih(3)
  integer :: io,ilo,ispin,nio,itag_send,itag_recv,itag_dir
  integer :: npg_base,rank_send,rank_recv
  integer :: dist_max,stencil_radius(3)
  integer :: ireq_send(1),ireq_recv(1)
  real(8) :: flux_coef
  real(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:)
  logical :: owns_send_face,owns_target_face

  if(.not. allocated(psi%rwf)) return
  if(.not. allocated(hpsi%rwf)) return
  if(dc%id_frag < 0) return
  if(info%ik_s /= 1 .or. info%ik_e /= 1) stop "DC-LCFO-Flux CG: ik/=1"
  if(info%im_s /= 1 .or. info%im_e /= 1) stop "DC-LCFO-Flux CG: im/=1"
  if(.not. stencil%if_orthogonal) stop "DC-LCFO-Flux CG: nonorthogonal lattice is not supported"
  if(dc%optimized_fragment_geometry) stop "DC-LCFO-Flux CG: optimized fragment geometry is not supported"
  if(mod(dc%isize_tot,dc%n_frag) /= 0) stop "DC-LCFO-Flux CG: MPI size must be divisible by number of fragments"

  call get_fragment_domain(dc,dc%i_frag,ncore)
  do n=1,3
    stencil_radius(n) = active_laplacian_radius(stencil,n)
  end do
  nio = info%io_e - info%io_s + 1
  if(nio <= 0) return
  if(info%io_s < lbound(psi%rwf,5) .or. info%io_e > ubound(psi%rwf,5)) &
  & stop "DC-LCFO-Flux CG: psi orbital range mismatch"
  if(info%io_s < lbound(hpsi%rwf,5) .or. info%io_e > ubound(hpsi%rwf,5)) &
  & stop "DC-LCFO-Flux CG: hpsi orbital range mismatch"
  if(system%nspin > ubound(psi%rwf,4) .or. system%nspin > ubound(hpsi%rwf,4)) &
  & stop "DC-LCFO-Flux CG: spin range mismatch"

  npg_base = dc%isize_tot / dc%n_frag

  call prepare_dc_lcfo_flux_halo_cache(dc,ncore,stencil_radius,npg_base)
  n_halo = flux_halo_cache_n
  halo => flux_halo_cache

  do i_halo=1,n_halo
    l = halo(i_halo)%length
    allocate(buf_send(l(1),l(2),l(3),system%nspin,nio))
    allocate(buf_recv(l(1),l(2),l(3),system%nspin,nio))
    buf_send = 0d0
    buf_recv = 0d0
    lo(1) = max(1,lbound(psi%rwf,1) - halo(i_halo)%dsp_send(1))
    lo(2) = max(1,lbound(psi%rwf,2) - halo(i_halo)%dsp_send(2))
    lo(3) = max(1,lbound(psi%rwf,3) - halo(i_halo)%dsp_send(3))
    hi(1) = min(l(1),ubound(psi%rwf,1) - halo(i_halo)%dsp_send(1))
    hi(2) = min(l(2),ubound(psi%rwf,2) - halo(i_halo)%dsp_send(2))
    hi(3) = min(l(3),ubound(psi%rwf,3) - halo(i_halo)%dsp_send(3))

    if(.not. any(lo > hi)) then
!$omp parallel do private(io,ilo,ispin,iz,iy,ix,ixg,iyg,izg) schedule(static)
      do io=info%io_s,info%io_e
        ilo = io - info%io_s + 1
        do ispin=1,system%nspin
        do iz=lo(3),hi(3)
        do iy=lo(2),hi(2)
        do ix=lo(1),hi(1)
          ixg = halo(i_halo)%dsp_send(1) + ix
          iyg = halo(i_halo)%dsp_send(2) + iy
          izg = halo(i_halo)%dsp_send(3) + iz
          buf_send(ix,iy,iz,ispin,ilo) = &
          & psi%rwf(ixg,iyg,izg,ispin,io,1,1)
        end do
        end do
        end do
        end do
      end do
!$omp end parallel do
    end if

    owns_send_face = owns_flux_send_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)
    owns_target_face = owns_flux_target_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)

    ! Tag convention follows the receiver-side halo label:
    ! dvec points from the source fragment to the destination fragment.
    ! A receiver with the same dvec waits for source=current-dvec, so both
    ! endpoints use source_fragment + offset(dvec)*n_frag.
    itag_dir = halo_tag_offset(halo(i_halo)%dvec)
    itag_send = dc%i_frag + itag_dir*dc%n_frag
    itag_recv = halo(i_halo)%ifrag_src + itag_dir*dc%n_frag
    rank_send = comm_proc_null
    rank_recv = comm_proc_null
    if(owns_send_face) rank_send = flux_send_destination_rank(info,dc,halo(i_halo),npg_base)
    if(owns_target_face) rank_recv = flux_recv_source_rank(info,dc,halo(i_halo),npg_base)
    ireq_send(1) = comm_isend(buf_send,rank_send,itag_send,dc%icomm_tot)
    ireq_recv(1) = comm_irecv(buf_recv,rank_recv,itag_recv,dc%icomm_tot)
    call comm_wait_all(ireq_recv)
    call comm_wait_all(ireq_send)

    if(.not. owns_target_face) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    axis = halo(i_halo)%axis
    dist_max = min(stencil_radius(axis),min(ncore(axis),l(axis)))
    lo = 1
    hi = l
    do n=1,3
      if(n /= axis) then
        lo(n) = max(lo(n),mg%is(n))
        hi(n) = min(hi(n),mg%ie(n))
      end if
    end do
    if(halo(i_halo)%dvec(axis) > 0) then
      lo(axis) = max(1,mg%is(axis))
      hi(axis) = min(dist_max,mg%ie(axis))
    else
      lo(axis) = max(1,ncore(axis) - mg%ie(axis) + 1)
      hi(axis) = min(dist_max,ncore(axis) - mg%is(axis) + 1)
    end if
    if(any(lo > hi)) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    if(halo(i_halo)%dvec(axis) > 0) then
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = dist
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    else
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = ncore(axis) - dist + 1
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%rwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    end if

    deallocate(buf_send,buf_recv)
  end do

end subroutine apply_dc_lcfo_flux_hpsi_rwf

subroutine apply_dc_lcfo_flux_hpsi_zwf(mg,system,info,stencil,dc,psi,hpsi)
  use structures
  use communication, only: comm_isend, comm_irecv, comm_wait_all, comm_proc_null
  use dc_fragment_geometry, only: get_fragment_domain
  implicit none
  type(s_rgrid),          intent(in)    :: mg
  type(s_dft_system),     intent(in)    :: system
  type(s_parallel_info),  intent(in)    :: info
  type(s_stencil),        intent(in)    :: stencil
  type(s_dcdft),          intent(in)    :: dc
  type(s_orbital),        intent(in)    :: psi
  type(s_orbital),        intent(inout) :: hpsi
  !
  type(halo_cg_info), pointer :: halo(:)
  integer :: ncore(3),l(3),lo(3),hi(3)
  integer :: i_halo,n_halo,n,axis,dist
  integer :: ix,iy,iz,ixg,iyg,izg,il(3),ih(3)
  integer :: io,ilo,ispin,nio,itag_send,itag_recv,itag_dir
  integer :: npg_base,rank_send,rank_recv
  integer :: dist_max,stencil_radius(3)
  integer :: ireq_send(1),ireq_recv(1)
  real(8) :: flux_coef
  complex(8),allocatable :: buf_send(:,:,:,:,:),buf_recv(:,:,:,:,:)
  logical :: owns_send_face,owns_target_face

  if(.not. allocated(psi%zwf)) return
  if(.not. allocated(hpsi%zwf)) return
  if(dc%id_frag < 0) return
  if(info%ik_s /= 1 .or. info%ik_e /= 1) stop "DC-LCFO-Flux: ik/=1"
  if(info%im_s /= 1 .or. info%im_e /= 1) stop "DC-LCFO-Flux: im/=1"
  if(.not. stencil%if_orthogonal) stop "DC-LCFO-Flux: nonorthogonal lattice is not supported"
  if(dc%optimized_fragment_geometry) stop "DC-LCFO-Flux: optimized fragment geometry is not supported"
  if(mod(dc%isize_tot,dc%n_frag) /= 0) stop "DC-LCFO-Flux: MPI size must be divisible by number of fragments"

  call get_fragment_domain(dc,dc%i_frag,ncore)
  do n=1,3
    stencil_radius(n) = active_laplacian_radius(stencil,n)
  end do
  nio = info%io_e - info%io_s + 1
  if(nio <= 0) return
  if(info%io_s < lbound(psi%zwf,5) .or. info%io_e > ubound(psi%zwf,5)) &
  & stop "DC-LCFO-Flux: psi orbital range mismatch"
  if(info%io_s < lbound(hpsi%zwf,5) .or. info%io_e > ubound(hpsi%zwf,5)) &
  & stop "DC-LCFO-Flux: hpsi orbital range mismatch"
  if(system%nspin > ubound(psi%zwf,4) .or. system%nspin > ubound(hpsi%zwf,4)) &
  & stop "DC-LCFO-Flux: spin range mismatch"

  npg_base = dc%isize_tot / dc%n_frag

  call prepare_dc_lcfo_flux_halo_cache(dc,ncore,stencil_radius,npg_base)
  n_halo = flux_halo_cache_n
  halo => flux_halo_cache

  do i_halo=1,n_halo
    l = halo(i_halo)%length
    allocate(buf_send(l(1),l(2),l(3),system%nspin,nio))
    allocate(buf_recv(l(1),l(2),l(3),system%nspin,nio))
    buf_send = (0d0,0d0)
    buf_recv = (0d0,0d0)
    lo(1) = max(1,lbound(psi%zwf,1) - halo(i_halo)%dsp_send(1))
    lo(2) = max(1,lbound(psi%zwf,2) - halo(i_halo)%dsp_send(2))
    lo(3) = max(1,lbound(psi%zwf,3) - halo(i_halo)%dsp_send(3))
    hi(1) = min(l(1),ubound(psi%zwf,1) - halo(i_halo)%dsp_send(1))
    hi(2) = min(l(2),ubound(psi%zwf,2) - halo(i_halo)%dsp_send(2))
    hi(3) = min(l(3),ubound(psi%zwf,3) - halo(i_halo)%dsp_send(3))

    if(.not. any(lo > hi)) then
!$omp parallel do private(io,ilo,ispin,iz,iy,ix,ixg,iyg,izg) schedule(static)
      do io=info%io_s,info%io_e
        ilo = io - info%io_s + 1
        do ispin=1,system%nspin
        do iz=lo(3),hi(3)
        do iy=lo(2),hi(2)
        do ix=lo(1),hi(1)
          ixg = halo(i_halo)%dsp_send(1) + ix
          iyg = halo(i_halo)%dsp_send(2) + iy
          izg = halo(i_halo)%dsp_send(3) + iz
          buf_send(ix,iy,iz,ispin,ilo) = &
          & psi%zwf(ixg,iyg,izg,ispin,io,1,1)
        end do
        end do
        end do
        end do
      end do
!$omp end parallel do
    end if

    owns_send_face = owns_flux_send_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)
    owns_target_face = owns_flux_target_face(mg,ncore,halo(i_halo)%dvec,stencil_radius)

    itag_dir = halo_tag_offset(halo(i_halo)%dvec)
    itag_send = dc%i_frag + itag_dir*dc%n_frag
    itag_recv = halo(i_halo)%ifrag_src + itag_dir*dc%n_frag
    rank_send = comm_proc_null
    rank_recv = comm_proc_null
    if(owns_send_face) rank_send = flux_send_destination_rank(info,dc,halo(i_halo),npg_base)
    if(owns_target_face) rank_recv = flux_recv_source_rank(info,dc,halo(i_halo),npg_base)
    ireq_send(1) = comm_isend(buf_send,rank_send,itag_send,dc%icomm_tot)
    ireq_recv(1) = comm_irecv(buf_recv,rank_recv,itag_recv,dc%icomm_tot)
    call comm_wait_all(ireq_recv)
    call comm_wait_all(ireq_send)

    if(.not. owns_target_face) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    axis = halo(i_halo)%axis
    dist_max = min(stencil_radius(axis),min(ncore(axis),l(axis)))
    lo = 1
    hi = l
    do n=1,3
      if(n /= axis) then
        lo(n) = max(lo(n),mg%is(n))
        hi(n) = min(hi(n),mg%ie(n))
      end if
    end do
    if(halo(i_halo)%dvec(axis) > 0) then
      lo(axis) = max(1,mg%is(axis))
      hi(axis) = min(dist_max,mg%ie(axis))
    else
      lo(axis) = max(1,ncore(axis) - mg%ie(axis) + 1)
      hi(axis) = min(dist_max,ncore(axis) - mg%is(axis) + 1)
    end if
    if(any(lo > hi)) then
      deallocate(buf_send,buf_recv)
      cycle
    end if

    if(halo(i_halo)%dvec(axis) > 0) then
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = dist
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    else
!$omp parallel do collapse(3) private(iz,iy,ix,ih,il,dist,flux_coef,io,ilo,ispin) schedule(static)
      do iz=lo(3),hi(3)
      do iy=lo(2),hi(2)
      do ix=lo(1),hi(1)
        ih = [ix,iy,iz]
        dist = ih(axis)
        il = ih
        il(axis) = ncore(axis) - dist + 1
        flux_coef = -0.5d0 * stencil%coef_lap(dist,axis)
        do io=info%io_s,info%io_e
          ilo = io - info%io_s + 1
          do ispin=1,system%nspin
            hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) = &
            & hpsi%zwf(il(1),il(2),il(3),ispin,io,1,1) &
            & + flux_coef * buf_recv(ix,iy,iz,ispin,ilo)
          end do
        end do
      end do
      end do
      end do
!$omp end parallel do
    end if

    deallocate(buf_send,buf_recv)
  end do

end subroutine apply_dc_lcfo_flux_hpsi_zwf

integer function flux_send_destination_rank(info,dc,halo,npg_base) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_dcdft),intent(in) :: dc
  type(halo_cg_info),intent(in) :: halo
  integer,intent(in) :: npg_base
  integer :: local_rank,face_side

  if(halo%dvec(halo%axis) > 0) then
    face_side = 0
  else
    face_side = info%nprgrid(halo%axis) - 1
  end if
  local_rank = rank_at_flux_face(info,halo%axis,face_side)
  rank = (halo%ifrag_dst - 1)*npg_base + local_rank
  if(rank < 0 .or. rank >= dc%isize_tot) stop "DC-LCFO-Flux CG: invalid send rank"
end function flux_send_destination_rank

integer function flux_recv_source_rank(info,dc,halo,npg_base) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  type(s_dcdft),intent(in) :: dc
  type(halo_cg_info),intent(in) :: halo
  integer,intent(in) :: npg_base
  integer :: local_rank,face_side

  if(halo%dvec(halo%axis) > 0) then
    face_side = info%nprgrid(halo%axis) - 1
  else
    face_side = 0
  end if
  local_rank = rank_at_flux_face(info,halo%axis,face_side)
  rank = (halo%ifrag_src - 1)*npg_base + local_rank
  if(rank < 0 .or. rank >= dc%isize_tot) stop "DC-LCFO-Flux CG: invalid recv rank"
end function flux_recv_source_rank

integer function rank_at_flux_face(info,axis,face_side) result(rank)
  use structures
  implicit none
  type(s_parallel_info),intent(in) :: info
  integer,intent(in) :: axis,face_side
  integer :: iaddr(5)

  iaddr = info%iaddress
  iaddr(axis) = face_side
  rank = info%imap(iaddr(1),iaddr(2),iaddr(3),iaddr(4),iaddr(5))
end function rank_at_flux_face

integer function active_laplacian_radius(stencil,axis) result(radius)
  use structures
  implicit none
  type(s_stencil),intent(in) :: stencil
  integer,intent(in) :: axis
  integer :: dist

  radius = 0
  do dist=1,size(stencil%coef_lap,1)
    if(abs(stencil%coef_lap(dist,axis)) > 0d0) radius = dist
  end do
end function active_laplacian_radius

logical function owns_flux_target_face(mg,ncore,dvec,stencil_radius) result(owns)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: ncore(3),dvec(3),stencil_radius(3)
  integer :: axis,idir,dist,ix(3)

  owns = .false.
  axis = 0
  do idir=1,3
    if(dvec(idir) /= 0) axis = idir
  end do
  if(axis == 0) return

  do dist=1,min(stencil_radius(axis),ncore(axis))
    ix = mg%is
    if(dvec(axis) > 0) then
      ix(axis) = dist
    else
      ix(axis) = ncore(axis) - dist + 1
    end if
    if(all(ix >= mg%is) .and. all(ix <= mg%ie)) then
      owns = .true.
      return
    end if
  end do
end function owns_flux_target_face

logical function owns_flux_send_face(mg,ncore,dvec,stencil_radius) result(owns)
  use structures
  implicit none
  type(s_rgrid),intent(in) :: mg
  integer,intent(in) :: ncore(3),dvec(3),stencil_radius(3)
  integer :: axis,idir,dist,ix(3)

  owns = .false.
  axis = 0
  do idir=1,3
    if(dvec(idir) /= 0) axis = idir
  end do
  if(axis == 0) return

  do dist=1,min(stencil_radius(axis),ncore(axis))
    ix = mg%is
    if(dvec(axis) > 0) then
      ix(axis) = ncore(axis) - dist + 1
    else
      ix(axis) = dist
    end if
    if(all(ix >= mg%is) .and. all(ix <= mg%ie)) then
      owns = .true.
      return
    end if
  end do
end function owns_flux_send_face

integer function halo_tag_offset(dvec) result(offset)
  implicit none
  integer,intent(in) :: dvec(3)

  offset = (dvec(1) + 1)*9 + (dvec(2) + 1)*3 + (dvec(3) + 1)
end function halo_tag_offset

end module dc_lcfo_flux_operator
