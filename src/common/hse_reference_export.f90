! Optional serial export for an independent HSE reference driver. No hybrid
! potential is inserted into SALMON by this module. The destination must exist.
#include "config.h"
module hse_reference_export
  use structures
  use iso_fortran_env, only: int8,int32
  implicit none
  private
  public :: export_hse_reference
contains
  subroutine export_hse_reference(lg,mg,system,info,stencil,srg,ppg,ppn,spsi,shpsi, &
                                  rho,vlocal,vh,vxc,vpsl,energy,iteration,converged)
    use hamiltonian, only: hpsi
    use salmon_global, only: yn_periodic,yn_spinorbit,yn_jm,yn_dc,yn_symmetrized_stencil,xc
    use pseudo_pt_plusU_sub, only: PLUS_U_ON
    type(s_rgrid),intent(in) :: lg,mg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    type(s_pp_nlcc),intent(in) :: ppn
    type(s_orbital),intent(inout) :: spsi,shpsi
    type(s_scalar),intent(in) :: rho,vlocal(:),vh,vxc(:),vpsl
    type(s_dft_energy),intent(in) :: energy
    integer,intent(in) :: iteration
    logical,intent(in) :: converged
    character(:),allocatable :: directory
    integer :: n,status,u,ios,ik,io,ilma,ia,j,ix,iy,iz,a,b,c,d,e,f
    complex(8),allocatable :: projectors(:,:,:,:)
    real(8),allocatable :: raw_uv(:,:),raw_xyz(:,:,:)
    integer,allocatable :: raw_idx(:,:,:),raw_count(:)
    integer(int8) :: endian_bytes(4)
    logical :: marker_exists
    call get_environment_variable('SALMON_HSE_REFERENCE_EXPORT',length=n,status=status)
    if(status/=0.or.n==0) return
    allocate(character(n)::directory)
    call get_environment_variable('SALMON_HSE_REFERENCE_EXPORT',value=directory,status=status)
    if(status/=0) error stop 'HSE reference export: cannot read directory'
    ! Invalidate an old snapshot before any validation or output. Only rank zero
    ! touches the marker; unsupported MPI export is rejected below on all ranks.
    if(info%id_rko==0) then
      inquire(file=directory//'/complete.txt',exist=marker_exists)
      if(marker_exists) then
        open(newunit=u,file=directory//'/complete.txt',status='old',iostat=ios)
        if(ios/=0) error stop 'HSE reference export: cannot remove stale completion marker'
        close(u,status='delete',iostat=ios)
        if(ios/=0) error stop 'HSE reference export: cannot delete stale completion marker'
      endif
    endif
#ifdef USE_OPENACC
    error stop 'HSE reference export: CPU only'
#endif
    if(info%isize_rko/=1.or.info%numm/=1) error stop 'HSE reference export: serial only'
    if(yn_periodic/='y'.or.system%nspin/=1.or..not.stencil%if_orthogonal) &
      error stop 'HSE reference export: periodic unpolarized orthogonal system required'
    if(yn_spinorbit/='n'.or.yn_jm/='n'.or.yn_dc/='n'.or.PLUS_U_ON.or. &
       yn_symmetrized_stencil=='y'.or.system%xc_payload%use_tau_operator) &
      error stop 'HSE reference export: unsupported Hamiltonian extension'
    if(.not.allocated(spsi%zwf).or..not.allocated(ppg%zekr_uv)) &
      error stop 'HSE reference export: complex periodic orbitals/projectors required'
    if(allocated(system%Ac_micro%v).or.any(abs(system%vec_Ac)>1d-14)) &
      error stop 'HSE reference export: zero field required'
    if(storage_size(1d0)/=64.or.storage_size((1d0,0d0))/=128) &
      error stop 'HSE reference export: unsupported real/complex storage size'
    a=mg%is(1);b=mg%ie(1);c=mg%is(2);d=mg%ie(2);e=mg%is(3);f=mg%ie(3)
    call hpsi(spsi,shpsi,info,mg,vlocal,system,stencil,srg,ppg)
    open(newunit=u,file=directory//'/metadata.txt',status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'HSE reference export: destination must exist and be writable'
    write(u,'(a)') 'SALMON_HSE_REFERENCE_V1'
    write(u,'(a)') '# Binary files: native-endian real64/complex128, unformatted stream, Fortran order.'
    write(u,'(a)') '# psi/hpsi: periodic u_nk, shape (nx,ny,nz,no,nk), physical grid without halos.'
    write(u,'(a)') '# projectors: shape (nx,ny,nz,nlma,nk), duplicate grid entries accumulated.'
    write(u,'(a)') '# Hnl psi = B diag(rinv_uvu) B^H psi, WITHOUT extra hvol in this native convention.'
    write(u,'(a)') '# Fields rho/vlocal/vh/vxc/vpsl: shape (nx,ny,nz); all quantities in atomic units.'
    write(u,'(a)') '# rho is the SCF mixed density used for potentials, not necessarily the orbital density.'
    write(u,'(a)') '# This is an instantaneous fixed-geometry snapshot, including when exported during optimization.'
    write(u,'(a)') '# k.bin (3,nk), weights.bin (nk), occupations.bin (no,nk), eigenvalues.bin (no,nk).'
    write(u,'(a)') '# geometry.bin: hgs(3), primitive_a(3,3), primitive_b(3,3), Rion(3,nion).'
    write(u,'(a)') '# coordinates.bin: x(nx),y(ny),z(nz); stencil.bin: lap0,lap(4,3),nab(4,3).'
    write(u,'(a)') '# kinetic=-.5 Laplacian - i k.grad + k^2/2; lap coefficients are unscaled derivative.'
    write(u,'(a)') '# Raw projectors: uv(nps,nlma); positions(3,nps,nlma); indices int32(3,nps,nlma).'
    write(u,'(a)') '# Counts int32(nlma); indices 1-based physical grid; inactive entries zero.'
    write(u,'(a)') '# B(k+A) accumulated as uv*exp(-i(k+A).position) at each raw index.'
    write(u,'(a)') '# energies.bin: total,kinetic,Hartree,xc,ion-ion,ion-local,ion-nonlocal.'
    endian_bytes=transfer(1_int32,endian_bytes)
    write(u,'(a,1x,l1)') 'endian_little',endian_bytes(1)==1_int8
    write(u,'(a)') 'rho_semantics scf_mixed_potential_density'
    write(u,'(a,3(1x,i0))') 'grid',mg%num
    write(u,'(a,1x,i0)') 'no',system%no
    write(u,'(a,1x,i0)') 'nk',system%nk
    write(u,'(a,1x,i0)') 'nion',system%nion
    write(u,'(a,1x,i0)') 'nlma',ppg%nlma
    write(u,'(a,1x,i0)') 'nps',ppg%nps
    write(u,'(a,1x,es25.16)') 'hvol',system%hvol
    write(u,'(a,1x,i0)') 'iteration',iteration
    write(u,'(a,1x,l1)') 'converged',converged
    write(u,'(a,1x,a)') 'xc',trim(xc)
    write(u,'(a,1x,l1)') 'nlcc_available',allocated(ppn%rho_nlcc)
    close(u)
    if(allocated(ppn%rho_nlcc)) then
      call open_binary('rho_nlcc.bin',u);write(u) ppn%rho_nlcc(a:b,c:d,e:f);close(u)
    endif
    call open_binary('psi.bin',u)
    do ik=1,system%nk
      do io=1,system%no
        write(u) spsi%zwf(a:b,c:d,e:f,1,io,ik,info%im_s)
      enddo
    enddo
    close(u)
    call open_binary('hpsi.bin',u)
    do ik=1,system%nk
      do io=1,system%no
        write(u) shpsi%zwf(a:b,c:d,e:f,1,io,ik,info%im_s)
      enddo
    enddo
    close(u)
    allocate(projectors(a:b,c:d,e:f,ppg%nlma))
    call open_binary('projectors.bin',u)
    do ik=1,system%nk
      projectors=0d0
      do ilma=1,ppg%nlma
        ia=ppg%ia_tbl(ilma)
        do j=1,ppg%mps(ia)
          ix=ppg%jxyz(1,j,ia);iy=ppg%jxyz(2,j,ia);iz=ppg%jxyz(3,j,ia)
          if(ix<a.or.ix>b.or.iy<c.or.iy>d.or.iz<e.or.iz>f) &
            error stop 'HSE reference export: projector outside physical serial grid'
          projectors(ix,iy,iz,ilma)=projectors(ix,iy,iz,ilma)+ppg%zekr_uv(j,ilma,ik)
        enddo
      enddo
      write(u) projectors
    enddo
    close(u)
    allocate(raw_uv(ppg%nps,ppg%nlma),raw_xyz(3,ppg%nps,ppg%nlma), &
             raw_idx(3,ppg%nps,ppg%nlma),raw_count(ppg%nlma))
    raw_uv=0d0;raw_xyz=0d0;raw_idx=0
    do ilma=1,ppg%nlma
      ia=ppg%ia_tbl(ilma);n=ppg%mps(ia);raw_count(ilma)=n
      raw_uv(1:n,ilma)=ppg%uv(1:n,ilma)
      raw_xyz(:,1:n,ilma)=ppg%rxyz(:,1:n,ia)
      do j=1,n
        raw_idx(:,j,ilma)=ppg%jxyz(:,j,ia)-mg%is+1
      enddo
    enddo
    if(storage_size(1)/=32) error stop 'HSE reference export: expected int32'
    call open_binary('raw_projectors.bin',u);write(u) raw_uv;close(u)
    call open_binary('projector_positions.bin',u);write(u) raw_xyz;close(u)
    call open_binary('projector_indices.bin',u);write(u) raw_idx;close(u)
    call open_binary('projector_counts.bin',u);write(u) raw_count;close(u)
    call open_binary('rinv_uvu.bin',u);write(u) ppg%rinv_uvu(1:ppg%nlma);close(u)
    call field('rho.bin',rho);call field('vlocal.bin',vlocal(1))
    call field('vh.bin',vh);call field('vxc.bin',vxc(1));call field('vpsl.bin',vpsl)
    call open_binary('k.bin',u);write(u) system%vec_k;close(u)
    call open_binary('weights.bin',u);write(u) system%wtk;close(u)
    call open_binary('occupations.bin',u);write(u) system%rocc(:,:,1);close(u)
    call open_binary('eigenvalues.bin',u);write(u) energy%esp(:,:,1);close(u)
    call open_binary('geometry.bin',u)
    write(u) system%hgs,system%primitive_a,system%primitive_b,system%Rion
    close(u)
    call open_binary('coordinates.bin',u)
    write(u) lg%coordinate(a:b,1),lg%coordinate(c:d,2),lg%coordinate(e:f,3)
    close(u)
    call open_binary('stencil.bin',u);write(u) stencil%coef_lap0,stencil%coef_lap,stencil%coef_nab;close(u)
    call open_binary('energies.bin',u)
    write(u) energy%E_tot,energy%E_kin,energy%E_h,energy%E_xc, &
             energy%E_ion_ion,energy%E_ion_loc,energy%E_ion_nloc
    close(u)
    ! A reader must require this marker. Write it only after all arrays and
    ! metadata have been successfully closed; an interrupted export is invalid.
    open(newunit=u,file=directory//'/complete.txt',status='replace',action='write',iostat=ios)
    if(ios/=0) error stop 'HSE reference export: cannot write completion marker'
    write(u,'(a)',iostat=ios) 'SALMON_HSE_REFERENCE_V1_COMPLETE'
    if(ios/=0) error stop 'HSE reference export: failed writing completion marker'
    close(u,iostat=ios)
    if(ios/=0) error stop 'HSE reference export: failed closing completion marker'
    write(*,'(a)') 'HSE reference export written to '//directory
  contains
    subroutine open_binary(name,unit)
      character(*),intent(in) :: name
      integer,intent(out) :: unit
      integer :: ierr
      open(newunit=unit,file=directory//'/'//name,status='replace',access='stream', &
           form='unformatted',action='write',iostat=ierr)
      if(ierr/=0) error stop 'HSE reference export: failed opening binary output'
    end subroutine
    subroutine field(name,value)
      character(*),intent(in) :: name
      type(s_scalar),intent(in) :: value
      integer :: unit
      call open_binary(name,unit)
      write(unit) value%f(a:b,c:d,e:f)
      close(unit)
    end subroutine
  end subroutine
end module
