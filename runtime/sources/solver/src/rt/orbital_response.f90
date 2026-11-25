! Spinless magnetic response of the existing molecular grid Hamiltonian.
! The result is ONLY -d<E>/dH; electric polarization current is not included.
module salmonpy_orbital_response
 use structures
 use salmonpy_grid_field,only:grid_magnetic_add,grid_phase_magnetic_add
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 private
 public::orbital_grid_magnetization
contains
 subroutine orbital_grid_magnetization(mg,lg,system,info,stencil,ppg,rt,psi,M,status)
  type(s_rgrid),intent(in)::mg,lg
  type(s_dft_system),intent(in)::system
  type(s_parallel_info),intent(in)::info
  type(s_stencil),intent(in)::stencil
  type(s_pp_grid),intent(in)::ppg
  type(s_rt),intent(in)::rt
  type(s_orbital),intent(in)::psi
  real(8),allocatable,intent(out)::M(:,:,:,:)
  integer,intent(out)::status
  integer::ix,iy,iz,io,ik,im,d,j,a,ilma,atom,code
  real(8)::r(3),q(3),kAc(3),occupation,point(3),finish(3),phase_gradient
  complex(8)::psi0,w(3),overlap,projection
  status=34
  if(.not.allocated(rt%near_grid_values).or..not.allocated(psi%zwf))return
  if(system%nspin/=1.or.info%im_s/=info%im_e.or..not.stencil%if_orthogonal)return
  if(any(mg%is/=lg%is).or.any(mg%ie/=lg%ie))return
  if(info%io_s/=1.or.info%io_e/=system%no.or.info%ik_s/=1.or.info%ik_e/=system%nk)return
  allocate(M(3,rt%near_grid_counts(1),rt%near_grid_counts(2),rt%near_grid_counts(3)),stat=code)
  status=32;if(code/=0)return
  M=0d0
  do iz=mg%is(3),mg%ie(3)
  do iy=mg%is(2),mg%ie(2)
  do ix=mg%is(1),mg%ie(1)
   point=[lg%coordinate(ix,1),lg%coordinate(iy,2),lg%coordinate(iz,3)]
   r=point-rt%near_geometry(1:3);q=0d0
   do im=info%im_s,info%im_e
   do ik=info%ik_s,info%ik_e
   do io=info%io_s,info%io_e
    occupation=system%rocc(io,ik,1)*system%wtk(ik)*system%Hvol
    psi0=psi%zwf(ix,iy,iz,1,io,ik,im);w=cmplx(0d0,0d0,8)
    do j=1,4
     w(1)=w(1)+stencil%coef_nab(j,1)*(psi%zwf(mg%idx(ix+j),iy,iz,1,io,ik,im)- &
                                             psi%zwf(mg%idx(ix-j),iy,iz,1,io,ik,im))
     w(2)=w(2)+stencil%coef_nab(j,2)*(psi%zwf(ix,mg%idy(iy+j),iz,1,io,ik,im)- &
                                             psi%zwf(ix,mg%idy(iy-j),iz,1,io,ik,im))
     w(3)=w(3)+stencil%coef_nab(j,3)*(psi%zwf(ix,iy,mg%idz(iz+j),1,io,ik,im)- &
                                             psi%zwf(ix,iy,mg%idz(iz-j),1,io,ik,im))
    enddo
    kAc=system%vec_k(:,ik)
    if(allocated(system%Ac_micro%v))kAc=kAc+system%Ac_micro%v(:,ix,iy,iz)
    q=q+occupation*(aimag(conjg(psi0)*w)+kAc*abs(psi0)**2)
   enddo
   enddo
   enddo
   call grid_magnetic_add(rt%near_grid_counts,rt%near_grid_origins(:,:,2),rt%near_grid_spacing, &
                         rt%near_grid_center,r,q,M,code)
   status=35;if(code/=0)return
  enddo
  enddo
  enddo
  ! The finite-rank nonlocal operator is |beta><beta|/uvu.
  ! beta_j=exp(-i phase_j)*uv_j, so dE/dphase_j=-2 Im(S* beta_j* psi_j)/uvu.
  do ilma=1,ppg%nlma
   atom=ppg%ia_tbl(ilma)
   do ik=info%ik_s,info%ik_e
   do im=info%im_s,info%im_e
   do io=info%io_s,info%io_e
    overlap=cmplx(0d0,0d0,8)
    do j=1,ppg%mps(atom)
     ix=ppg%jxyz(1,j,atom);iy=ppg%jxyz(2,j,atom);iz=ppg%jxyz(3,j,atom)
     overlap=overlap+conjg(ppg%zekr_uV(j,ilma,ik))*psi%zwf(ix,iy,iz,1,io,ik,im)
    enddo
    occupation=system%rocc(io,ik,1)*system%wtk(ik)*system%Hvol*ppg%rinv_uvu(ilma)
    do j=1,ppg%mps(atom)
     ix=ppg%jxyz(1,j,atom);iy=ppg%jxyz(2,j,atom);iz=ppg%jxyz(3,j,atom)
     projection=conjg(ppg%zekr_uV(j,ilma,ik))*psi%zwf(ix,iy,iz,1,io,ik,im)
     phase_gradient=-2d0*occupation*aimag(conjg(overlap)*projection)
     if(phase_gradient==0d0)cycle
     finish=system%Rion(:,atom)+ppg%rxyz(:,j,atom)
     call grid_phase_magnetic_add(rt%near_grid_counts,rt%near_grid_origins(:,:,2),rt%near_grid_spacing, &
      rt%near_grid_center,rt%near_geometry(1:3),system%Rion(:,atom),finish,phase_gradient,M,code)
     if(code/=0)return
    enddo
   enddo
   enddo
   enddo
  enddo
  status=35
  do iz=1,rt%near_grid_counts(3)
  do iy=1,rt%near_grid_counts(2)
  do ix=1,rt%near_grid_counts(1)
  do d=1,3
   if(.not.ieee_is_finite(M(d,ix,iy,iz)))return
  enddo
  enddo
  enddo
  enddo
  status=0
 end subroutine
end module
