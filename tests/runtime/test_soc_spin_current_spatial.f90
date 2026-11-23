! Fixed spinors/projectors: serial versus spatial nonlocal Pauli current.
! Includes a projector absent from one rank and independent spatial replicas.
#include "config.h"
program test_soc_spin_current_spatial
 use mpi
 use communication,only:comm_summation
 use structures,only:s_pp_grid,s_dft_system,s_rgrid,s_stencil,s_parallel_info,s_orbital
 use noncollinear_module,only:calc_spin_current
 use pseudo_pt_current_so,only:calc_spin_current_nonlocal
#ifndef SOC_SERIAL_ONLY
 use pseudo_pt_current_so,only:calc_spin_current_nonlocal_rdivided
#endif
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 type(s_pp_grid)::full,part
 type(s_dft_system)::sys
 type(s_rgrid)::grid,local_grid
 type(s_stencil)::stencil
 type(s_parallel_info)::info,local_info
 type(s_orbital)::orb
 real(8)::micro(3,0:3,8,2,2),band(3,0:3,1,1),local_band(3,0:3,1,1)
 real(8),allocatable::local_micro(:,:,:,:,:)
 complex(8)::psi(8,2,2,2),u0,r0(3,0:3),v0,pauli(2,2,0:3),spinor(2),projector(2)
 real(8)::manual(3,0:3)
 real(8)::reference(3,0:3),local(3,0:3),total(3,0:3),err,scale,largest
 integer::rank,np,comm_r,rr,nr,ierr,i,j,k,s,a,l,n,source_j(4,2),owner(3),lo,hi,ch,d
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 if(np/=1.and.np/=2.and.np/=4)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
 call MPI_Comm_split(MPI_COMM_WORLD,rank/2,rank,comm_r,ierr)
 call MPI_Comm_rank(comm_r,rr,ierr);call MPI_Comm_size(comm_r,nr,ierr)
 lo=1+rr*(8/nr);hi=lo+8/nr-1
 scale=1d0+.2d0*(rank/2)
 do s=1,2
 do k=1,2
 do j=1,2
 do i=1,8
  psi(i,j,k,s)=scale*cmplx(.2d0*i+.3d0*j-.1d0*k+.07d0*s,sin(.21d0*i+.33d0*j+.12d0*k+s),8)
 enddo
 enddo
 enddo
 enddo
 owner=[1,1,2]
 allocate(full%mps(2),full%ia_tbl_so(3),full%rinv_uvu_so(3))
 allocate(full%jxyz(3,4,2),full%rxyz(3,4,2),full%zekr_uv_so(4,3,1,2,1))
 full%mps=4;full%ia_tbl_so=owner;full%rinv_uvu_so=[.7d0,-.4d0,1.1d0]
 do a=1,2
 do j=1,4
  i=2*j-1
  if(a==2)i=j
  full%jxyz(:,j,a)=[i,1+mod(j,2),1+mod(j/2,2)]
  full%rxyz(:,j,a)=[.2d0*j+.03d0*a,-.1d0*j+.04d0*a,.3d0*j-.02d0*a]
 enddo
 enddo
 do l=1,3
 do s=1,2
 do j=1,4
  full%zekr_uv_so(j,l,1,s,1)=cmplx(.11d0*j+.03d0*l+.01d0*s,.05d0*j*s-.02d0*l,8)
 enddo
 enddo
 enddo
 part=full;part%mps=0
 do a=1,2
 do j=1,4
  i=full%jxyz(1,j,a)
  if(i<lo.or.i>hi)cycle
  n=part%mps(a)+1;part%mps(a)=n;source_j(n,a)=j
  part%jxyz(:,n,a)=full%jxyz(:,j,a);part%rxyz(:,n,a)=full%rxyz(:,j,a)
 enddo
 enddo
 do l=1,3
 a=owner(l)
 do j=1,part%mps(a)
  part%zekr_uv_so(j,l,1,:,1)=full%zekr_uv_so(source_j(j,a),l,1,:,1)
 enddo
 enddo
 pauli=(0d0,0d0)
 pauli(1,1,0)=1d0;pauli(2,2,0)=1d0
 pauli(1,2,1)=1d0;pauli(2,1,1)=1d0
 pauli(1,2,2)=(0d0,-1d0);pauli(2,1,2)=(0d0,1d0)
 pauli(1,1,3)=1d0;pauli(2,2,3)=-1d0
 manual=0d0
 do l=1,3
  a=owner(l);u0=(0d0,0d0);r0=(0d0,0d0)
  do j=1,4
   i=full%jxyz(1,j,a);k=full%jxyz(2,j,a);n=full%jxyz(3,j,a)
   spinor=psi(i,k,n,:);projector=full%zekr_uv_so(j,l,1,:,1)
   u0=u0+dot_product(projector,spinor)
   do ch=0,3
    v0=dot_product(projector,matmul(pauli(:,:,ch),spinor))
    r0(:,ch)=r0(:,ch)+full%rxyz(:,j,a)*v0
   enddo
  enddo
  manual=manual+2d0*aimag(conjg(r0)*u0*full%rinv_uvu_so(l))
 enddo
 call calc_spin_current_nonlocal(reference,psi,full,[1,1,1],[8,2,2],1)
#ifdef SOC_SERIAL_ONLY
 total=reference
#else
 call calc_spin_current_nonlocal_rdivided(local,psi,part,[1,1,1],[8,2,2],1,comm_r)
 call comm_summation(local,total,12,comm_r)
#endif
 do s=0,3
 do i=1,3
  if(.not.ieee_is_finite(total(i,s)))call MPI_Abort(MPI_COMM_WORLD,2,ierr)
 enddo
 enddo
 err=max(maxval(abs(total-reference)),maxval(abs(reference-manual)));largest=max(1d0,maxval(abs(reference)))
 if(rank==0)print *,'SOC_NONLOCAL_REFERENCE_ERROR',np,err
 if(err>5d-13*largest)call MPI_Abort(MPI_COMM_WORLD,3,ierr)
 if(rank==0)print *,'SOC_SPATIAL_NONLOCAL_PASS',np,err
#ifdef SOC_SERIAL_ONLY
 call MPI_Comm_free(comm_r,ierr)
 call MPI_Finalize(ierr)
 stop
#endif
 ! Full local gradient + nonlocal Pauli current with prescribed valid halos.
 sys%no=1;sys%nk=1;sys%nspin=2;sys%ngrid=32
 allocate(sys%rocc(1,1,1),sys%wtk(1),sys%vec_k(3,1))
 sys%rocc=.8d0;sys%wtk=1d0;sys%vec_k(:,1)=[.12d0,-.05d0,.08d0]
 sys%vec_Ac=[-.01d0,.02d0,.03d0];sys%rmatrix_B=0d0
 do i=1,3
  sys%rmatrix_B(i,i)=1d0
 enddo
 stencil%coef_nab=0d0;stencil%coef_nab(1,:)=.5d0
 grid%is=1;grid%ie=[8,2,2];grid%num=grid%ie;grid%is_array=1;grid%ie_array=grid%ie
 allocate(grid%idx(-3:12),grid%idy(-3:6),grid%idz(-3:6))
 do i=-3,12
  grid%idx(i)=1+modulo(i-1,8)
 enddo
 do i=-3,6
  grid%idy(i)=1+modulo(i-1,2);grid%idz(i)=1+modulo(i-1,2)
 enddo
 local_grid=grid;local_grid%is(1)=lo;local_grid%ie(1)=hi;local_grid%num(1)=hi-lo+1
 deallocate(local_grid%idx);allocate(local_grid%idx(lo-4:hi+4))
 do i=lo-4,hi+4
  local_grid%idx(i)=1+modulo(i-1,8)
 enddo
 allocate(orb%zwf(8,2,2,2,1,1,1));orb%zwf(:,:,:,:,1,1,1)=psi
 info%ik_s=1;info%ik_e=1;info%io_s=1;info%io_e=1
 info%if_divide_rspace=.false.;info%icomm_ko=MPI_COMM_SELF;info%icomm_rko=MPI_COMM_SELF
 call calc_spin_current(sys,grid,stencil,info,orb,full,micro,band)
 local_info=info;local_info%if_divide_rspace=nr>1
 local_info%icomm_r=comm_r;local_info%icomm_rko=comm_r
 allocate(local_micro(3,0:3,lo:hi,1:2,1:2))
 call calc_spin_current(sys,local_grid,stencil,local_info,orb,part,local_micro,local_band)
 do ch=0,3
 do d=1,3
  if(.not.ieee_is_finite(local_band(d,ch,1,1)))call MPI_Abort(MPI_COMM_WORLD,5,ierr)
  do k=1,2
  do j=1,2
  do i=lo,hi
   if(.not.ieee_is_finite(local_micro(d,ch,i,j,k)))call MPI_Abort(MPI_COMM_WORLD,6,ierr)
  enddo
  enddo
  enddo
 enddo
 enddo
 err=max(maxval(abs(local_band-band)),maxval(abs(local_micro-micro(:,:,lo:hi,:,:))))
 if(err>5d-13*max(1d0,maxval(abs(band)),maxval(abs(micro))))call MPI_Abort(MPI_COMM_WORLD,4,ierr)
 if(rank==0)print *,'SOC_SPATIAL_TOTAL_CURRENT_PASS',np,err
 call MPI_Comm_free(comm_r,ierr)
 call MPI_Finalize(ierr)
end program
