! Compensated products/sums for original-H refinement. VOLATILE temporaries
! preserve error-free transformations even in builds enabling reassociation.
module hsex_compensated
  implicit none
  private
  public :: hsex_original_residual,hsex_precise_gram,hsex_split_matrix,hsex_refinement_residual
contains
  subroutine add_product(a,b,s,c)
    implicit none
    real(8),intent(in)::a,b
    real(8),intent(inout)::s,c
    real(8),volatile::p,ta,tb,ah,al,bh,bl,t,z,e,v,w
    real(8),parameter::split=134217729d0
    p=a*b
    ta=split*a;tb=split*b
    v=ta-a;ah=ta-v;al=a-ah
    v=tb-b;bh=tb-v;bl=b-bh
    v=ah*bh;v=v-p
    w=ah*bl;v=v+w
    w=al*bh;v=v+w
    w=al*bl;e=v+w
    t=s+p;z=t-s
    v=t-z;v=s-v;w=p-z
    v=v+w;e=e+v
    c=c+e;s=t
  end subroutine

  subroutine hsex_split_matrix(n,m,a,rows,bits,parts,tail)
    implicit none
    integer,intent(in)::n,m,bits
    complex(8),intent(in)::a(n,m)
    logical,intent(in)::rows
    complex(8),intent(out)::parts(:,:,:)
    complex(8),intent(in),optional::tail(n,m)
    complex(8),allocatable::extra(:,:)
    complex(8),volatile::t,z,v,w,previous
    complex(8),allocatable::rest(:,:)
    real(8)::largest,ar,ai
    integer::i,j,k,expo,group,count,entry
    allocate(rest(n,m));rest=a
    if(present(tail))then
      allocate(extra(n,m));extra=tail
    endif
    count=merge(n,m,rows)
    do k=1,size(parts,3)
      do group=1,count
        largest=0d0
        do entry=1,merge(m,n,rows)
          i=merge(group,entry,rows);j=merge(entry,group,rows)
          largest=max(largest,abs(real(rest(i,j),8)),abs(aimag(rest(i,j))))
        enddo
        expo=exponent(largest)
        do entry=1,merge(m,n,rows)
          i=merge(group,entry,rows);j=merge(entry,group,rows)
          ar=scale(anint(scale(real(rest(i,j),8),bits-expo)),expo-bits)
          ai=scale(anint(scale(aimag(rest(i,j)),bits-expo)),expo-bits)
          parts(i,j,k)=cmplx(ar,ai,8);rest(i,j)=rest(i,j)-parts(i,j,k)
          if(present(tail))then
            previous=rest(i,j);t=previous+extra(i,j);z=t-previous
            v=t-z;v=previous-v;w=extra(i,j)-z;v=v+w
            rest(i,j)=t;extra(i,j)=v
          endif
        enddo
      enddo
    enddo
  end subroutine

  subroutine hsex_refinement_residual(n,g,hparts,hdiag,nh,bits,energy,eta,x,xlow,b,r)
    implicit none
    integer,intent(in)::n,g,nh,bits
    real(8),intent(in)::energy,eta,hdiag(n)
    complex(8),intent(in)::hparts(:,:,:),x(n,g),xlow(n,g),b(n,g)
    complex(8),intent(out)::r(n,g)
    complex(8),allocatable::xp(:,:,:),term(:,:),low(:,:)
    complex(8),volatile::t,z,v,w,previous
    real(8)::sr,si,cr,ci,xr,xi
    real(8),volatile::shift,shift_low,zz,vv,ww
    integer::order,ia,ib,i,j,nslices
    nslices=size(hparts,3)
    allocate(xp(n,g,nslices),term(n,g),low(n,g))
    if(nh>0)call hsex_split_matrix(n,g,x,.false.,bits,xp,xlow)
    r=b;low=0d0
    ! Each fixed-point dot product is exact in binary64: reserve bits for
    ! n terms, the two real products per complex component, and carries.
    ! Additional slice orders offset eta^-2 amplification for smaller eta.
    do order=0,nslices-1
      do ia=1,min(order+1,nh)
        ib=order+2-ia
        call zgemm('N','N',n,g,n,(-1d0,0d0),hparts(:,:,ia),n,xp(:,:,ib),n,(0d0,0d0),term,n)
        do j=1,g
          do i=1,n
            previous=r(i,j);t=previous+term(i,j);z=t-previous
            v=t-z;v=previous-v;w=term(i,j)-z;v=v+w
            low(i,j)=low(i,j)+v;r(i,j)=t
          enddo
        enddo
      enddo
    enddo
    do j=1,g
      do i=1,n
        ! Exact two-part E-H_ii avoids multiplying two equal resonant terms
        ! separately. Keeping the low part preserves the original operator.
        shift=energy-hdiag(i);zz=shift-energy
        vv=shift-zz;vv=energy-vv;ww=-hdiag(i)-zz;shift_low=vv+ww
        sr=real(r(i,j),8);si=aimag(r(i,j));cr=real(low(i,j),8);ci=aimag(low(i,j))
        xr=real(x(i,j),8);xi=aimag(x(i,j))
        call add_product(shift,xr,sr,cr)
        call add_product(shift_low,xr,sr,cr)
        call add_product(-eta,xi,sr,cr)
        call add_product(shift,xi,si,ci)
        call add_product(shift_low,xi,si,ci)
        call add_product(eta,xr,si,ci)
        xr=real(xlow(i,j),8);xi=aimag(xlow(i,j))
        call add_product(shift,xr,sr,cr)
        call add_product(shift_low,xr,sr,cr)
        call add_product(-eta,xi,sr,cr)
        call add_product(shift,xi,si,ci)
        call add_product(shift_low,xi,si,ci)
        call add_product(eta,xr,si,ci)
        r(i,j)=cmplx(sr+cr,si+ci,8)
      enddo
    enddo
  end subroutine

  subroutine hsex_original_residual(n,g,h,energy,eta,x,b,r) bind(C)
    use iso_c_binding, only: c_int,c_double,c_double_complex
    implicit none
    integer(c_int),value,intent(in)::n,g
    real(c_double),value,intent(in)::energy,eta
    complex(c_double_complex),intent(in)::h(n,n),x(n,g),b(n,g)
    complex(c_double_complex),intent(out)::r(n,g)
    real(8)::sr,si,cr,ci,hr,hi,xr,xi
    integer::i,j,k
    do j=1,g
      do i=1,n
        sr=real(b(i,j),8);si=aimag(b(i,j));cr=0d0;ci=0d0
        do k=1,n
          hr=real(h(i,k),8);hi=aimag(h(i,k));xr=real(x(k,j),8);xi=aimag(x(k,j))
          call add_product(-hr,xr,sr,cr)
          call add_product(hi,xi,sr,cr)
          call add_product(-hr,xi,si,ci)
          call add_product(-hi,xr,si,ci)
        enddo
        ! Form the shift in compensated arithmetic, not by rounding H_ii-E.
        xr=real(x(i,j),8);xi=aimag(x(i,j))
        call add_product(energy,xr,sr,cr)
        call add_product(-eta,xi,sr,cr)
        call add_product(energy,xi,si,ci)
        call add_product(eta,xr,si,ci)
        r(i,j)=cmplx(sr+cr,si+ci,8)
      enddo
    enddo
  end subroutine

  subroutine hsex_precise_gram(n,g,ld,b,x,xlow,gram) bind(C)
    use iso_c_binding, only: c_int,c_double_complex
    implicit none
    integer(c_int),value,intent(in)::n,g,ld
    complex(c_double_complex),intent(in)::b(n,g),x(n,g),xlow(n,g)
    complex(c_double_complex),intent(out)::gram(ld,g)
    real(8)::sr,si,cr,ci,br,bi,xr,xi
    integer::i,j,k
    do j=1,g
      do i=1,j
        sr=0d0;si=0d0;cr=0d0;ci=0d0
        do k=1,n
          br=.5d0*real(b(k,i),8);bi=.5d0*aimag(b(k,i));xr=real(x(k,j),8);xi=aimag(x(k,j))
          call add_product(br,xr,sr,cr)
          call add_product(bi,xi,sr,cr)
          call add_product(br,xi,si,ci)
          call add_product(-bi,xr,si,ci)
          br=.5d0*real(x(k,i),8);bi=.5d0*aimag(x(k,i));xr=real(b(k,j),8);xi=aimag(b(k,j))
          call add_product(br,xr,sr,cr)
          call add_product(bi,xi,sr,cr)
          call add_product(br,xi,si,ci)
          call add_product(-bi,xr,si,ci)
          br=.5d0*real(b(k,i),8);bi=.5d0*aimag(b(k,i));xr=real(xlow(k,j),8);xi=aimag(xlow(k,j))
          call add_product(br,xr,sr,cr)
          call add_product(bi,xi,sr,cr)
          call add_product(br,xi,si,ci)
          call add_product(-bi,xr,si,ci)
          br=.5d0*real(xlow(k,i),8);bi=.5d0*aimag(xlow(k,i));xr=real(b(k,j),8);xi=aimag(b(k,j))
          call add_product(br,xr,sr,cr)
          call add_product(bi,xi,sr,cr)
          call add_product(br,xi,si,ci)
          call add_product(-bi,xr,si,ci)
        enddo
        gram(i,j)=cmplx(sr+cr,si+ci,8)
        if(i==j)gram(i,j)=cmplx(sr+cr,0d0,8)
        gram(j,i)=conjg(gram(i,j))
      enddo
    enddo
  end subroutine
end module hsex_compensated
#ifdef HSEX_SALMON
#include "config.h"
#endif
! Shared target rotation, private q workspaces. Called by SALMON and the reference harness.
module hsex_q_response
  implicit none
  private
  public :: hsex_q_solve,hsex_q_abi
contains
  function hsex_q_abi() result(version) bind(C)
    use iso_c_binding, only: c_int
    implicit none
    integer(c_int)::version
    version=2
  end function

  subroutine hsex_q_solve(n,nb,nq,mg,nu,nc,nthreads,pump,adj,d,e,u,psi,phase,counts,levels,occ,eta, &
                         response,residual,status,team,refine,original_h) bind(C)
    use iso_c_binding, only: c_int,c_double,c_double_complex
    use ieee_arithmetic, only: ieee_is_finite
    use hsex_compensated, only: hsex_precise_gram,hsex_split_matrix,hsex_refinement_residual
#ifdef HSEX_SALMON
    use exx_blas_threads, only: exx_blas_thread_control
#endif
!$  use omp_lib, only: omp_get_num_threads
    implicit none
    integer(c_int),value,intent(in) :: n,nb,nq,mg,nu,nc,nthreads,pump,refine
    complex(c_double_complex),intent(in) :: adj(n,n),u(n,nu),psi(n,nb,nq),phase(n,mg,nq)
    complex(c_double_complex),intent(in) :: original_h(merge(n,1,refine/=0),merge(n,1,refine/=0))
    real(c_double),intent(in) :: d(n),e(n-1),levels(nb,nq),occ(nb,nq),eta(nc)
    integer(c_int),intent(in) :: counts(nq)
    complex(c_double_complex),intent(out) :: response(mg,mg,nc,nq)
    real(c_double),intent(out) :: residual(nq)
    integer(c_int),intent(out) :: status(nq),team
    complex(c_double_complex),allocatable :: b(:,:),rot(:,:),x(:,:),proj(:,:),gram(:,:)
    complex(c_double_complex),allocatable :: real_x(:,:),defect(:,:),correction(:,:),real_low(:,:),old_gram(:,:)
    complex(c_double_complex),allocatable :: diag(:),dl(:),du(:),du2(:)
    complex(c_double_complex),allocatable :: hparts(:,:,:),hoff(:,:)
    real(c_double),allocatable :: hdiag(:)
    real(c_double),allocatable :: rd(:)
    integer,allocatable :: piv(:)
    complex(c_double_complex) :: r
    complex(c_double_complex),volatile :: previous,t,z,vv,ww
    logical :: converged
    real(c_double) :: rn,bn,err,floor_norm
    integer :: q,v,c,g,i,j,info,old_mode,api_status,blas_worker(3),blas_parent(3),iteration,bits,nslices,nh
#ifdef HSEX_ACCELERATE
    interface
      integer(c_int) function blas_set(mode) bind(C,name='BLASSetThreading')
        import c_int
        integer(c_int),value :: mode
      end function
      integer(c_int) function blas_get() bind(C,name='BLASGetThreading')
        import c_int
      end function
    end interface
#endif
    response=0d0;residual=0d0;status=0;team=1
    bits=(50-ceiling(log(dble(max(n,2)))/log(2d0)))/2
    if(refine/=0)then
      nslices=5+ceiling(max(0d0,2d0*(log(5d-4)-log(minval(eta)))/log(2d0))/bits)
      if(nslices>12)then
        status=5
        return
      endif
      allocate(hparts(n,n,nslices),hdiag(n))
      allocate(hoff(n,n));hoff=original_h
      do i=1,n
        hdiag(i)=real(original_h(i,i),c_double)
        hoff(i,i)=cmplx(0d0,aimag(original_h(i,i)),c_double)
      enddo
      call hsex_split_matrix(n,n,hoff,.true.,bits,hparts)
      deallocate(hoff)
      ! Skip only identically zero matrix slices, never small nonzero values.
      nh=0
      do i=1,nslices
        if(any(hparts(:,:,i)/=cmplx(0d0,0d0,c_double)))nh=i
      enddo
    endif
    blas_parent=0
#if defined(HSEX_SALMON) && defined(HAVE_EXX_OPENBLAS_THREADS)
    call exx_blas_thread_control(1,blas_parent)
#endif
    !$omp parallel num_threads(nthreads) default(shared) &
    !$omp& private(b,rot,x,proj,gram,diag,dl,du,du2,rd,piv,q,v,c,g,i,j,info,r,rn,bn,err,floor_norm, &
    !$omp& old_mode,api_status,blas_worker,real_x,defect,correction,iteration,real_low,old_gram, &
    !$omp& previous,t,z,vv,ww,converged)
    !$omp single
!$  team=omp_get_num_threads()
    !$omp end single
    old_mode=0;api_status=0;blas_worker=0
#ifdef HSEX_SALMON
    call exx_blas_thread_control(1,blas_worker)
#endif
#ifdef HSEX_ACCELERATE
    old_mode=blas_get()
    api_status=blas_set(1_c_int)
#endif
    allocate(b(n,mg),rot(n,mg),x(n,mg),proj(max(1,nu),mg),gram(mg,mg))
    if(refine/=0)allocate(real_x(n,mg),defect(n,mg),correction(n,mg),real_low(n,mg),old_gram(mg,mg))
    allocate(diag(n),dl(n-1),du(n-1),du2(max(1,n-2)),rd(n),piv(n))
    !$omp do schedule(dynamic,1)
    do q=1,nq
      if(api_status/=0)then
        status(q)=-1
        cycle
      endif
#if defined(HSEX_SALMON) && defined(HAVE_HSEX_APPLE_GNU) && defined(HAVE_EXX_OPENBLAS_THREADS)
      ! Match the established SALMON guard for mixed GNU/LLVM OpenMP OpenBLAS.
      ! Independent q work remains parallel on compatible production backends.
      !$omp critical(hsex_q_vendor_serial)
#endif
      g=counts(q)
      do v=1,nb
        if(occ(v,q)==0d0)cycle
        do j=1,g
          do i=1,n
            b(i,j)=psi(i,v,q)*phase(i,j,q)
          enddo
        enddo
        if(nu>0)then
          call zgemm('C','N',nu,g,n,(1d0,0d0),u,n,b,n,(0d0,0d0),proj,max(1,nu))
          call zgemm('N','N',n,g,nu,(-1d0,0d0),u,n,proj,max(1,nu),(1d0,0d0),b,n)
        endif
        call zgemm('N','N',n,g,n,(1d0,0d0),adj,n,b,n,(0d0,0d0),rot,n)
        do c=1,nc
          x(:,:g)=rot(:,:g)
          if(pump==0)then
            rd=d-levels(v,q);dl=cmplx(e,0d0,c_double)
            call zptsv(n,g,rd,dl,x,n,info)
            floor_norm=1d-12
          else
            diag=cmplx(d-levels(v,q),-eta(c),c_double)
            dl=cmplx(e,0d0,c_double);du=dl
            call zgttrf(n,dl,diag,du,du2,piv,info)
            if(info==0)call zgttrs('N',n,g,dl,diag,du,du2,piv,x,n,info)
            floor_norm=1d-14
          endif
          if(info/=0)then
            status(q)=1
            cycle
          endif
          do j=1,g
            rn=0d0;bn=0d0
            do i=1,n
              r=cmplx(d(i)-levels(v,q),-eta(c),c_double)*x(i,j)-rot(i,j)
              if(i>1)r=r+e(i-1)*x(i-1,j)
              if(i<n)r=r+e(i)*x(i+1,j)
              rn=rn+abs(r)**2;bn=bn+abs(rot(i,j))**2
            enddo
            err=sqrt(rn)/max(sqrt(bn),floor_norm)
            if(.not.ieee_is_finite(err))then
              status(q)=2
            else if(refine==0.or.pump==0)then
              residual(q)=max(residual(q),err)
              if(err>1d-8)status(q)=2
            endif
          enddo
          ! With refinement, the final ORIGINAL-H residual below is decisive;
          ! the initial transformed solve is only a preconditioned estimate.
          if(refine/=0.and.pump/=0)then
            call zgemm('C','N',n,g,n,(1d0,0d0),adj,n,x,n,(0d0,0d0),real_x,n)
            real_low=0d0;old_gram=0d0;converged=.false.
            do iteration=1,10
              call hsex_refinement_residual(n,g,hparts,hdiag,nh,bits,levels(v,q),eta(c),real_x,real_low,b,defect)
              call zgemm('N','N',n,g,n,(1d0,0d0),adj,n,defect,n,(0d0,0d0),correction,n)
              call zgttrs('N',n,g,dl,diag,du,du2,piv,correction,n,info)
              if(info/=0)then
                status(q)=1
                exit
              endif
              call zgemm('C','N',n,g,n,(1d0,0d0),adj,n,correction,n,(0d0,0d0),defect,n)
              do j=1,g
                do i=1,n
                  previous=real_x(i,j);t=previous+defect(i,j);z=t-previous
                  vv=t-z;vv=previous-vv;ww=defect(i,j)-z;vv=vv+ww
                  real_low(i,j)=real_low(i,j)+vv;real_x(i,j)=t
                  previous=real_x(i,j);t=previous+real_low(i,j);z=t-previous
                  vv=t-z;vv=previous-vv;ww=real_low(i,j)-z;vv=vv+ww
                  real_x(i,j)=t;real_low(i,j)=vv
                enddo
              enddo
              call hsex_precise_gram(n,g,mg,b,real_x,real_low,gram)
              if(iteration>=2)then
                converged=.true.
                do j=1,g
                  do i=1,g
                    err=abs(gram(i,j)-old_gram(i,j))/(1d0+abs(gram(i,j)))
                    if(.not.ieee_is_finite(err))then
                      converged=.false.
                    else if(err>1d-12)then
                      converged=.false.
                    endif
                  enddo
                enddo
                if(converged)exit
              endif
              old_gram(:g,:g)=gram(:g,:g)
            enddo
            if(info/=0)cycle
            if(.not.converged.and.status(q)==0)status(q)=4
            call hsex_refinement_residual(n,g,hparts,hdiag,nh,bits,levels(v,q),eta(c),real_x,real_low,b,defect)
            do j=1,g
              rn=sum(abs(defect(:,j))**2);bn=sum(abs(b(:,j))**2)
              err=sqrt(rn)/max(sqrt(bn),floor_norm)
              if(.not.ieee_is_finite(err))then
                status(q)=2
              else
                residual(q)=max(residual(q),err)
                if(err>1d-8)status(q)=2
              endif
            enddo
          else
            call zgemm('C','N',g,g,n,(1d0,0d0),rot,n,x,n,(0d0,0d0),gram,mg)
          endif
          do j=1,g
            do i=1,g
              r=gram(i,j)
              if(pump/=0)r=.5d0*(r+conjg(gram(j,i)))
              response(i,j,c,q)=response(i,j,c,q)+.5d0*occ(v,q)*r
              if(.not.ieee_is_finite(real(response(i,j,c,q),c_double)))status(q)=3
              if(.not.ieee_is_finite(aimag(response(i,j,c,q))))status(q)=3
            enddo
          enddo
        enddo
      enddo
#if defined(HSEX_SALMON) && defined(HAVE_HSEX_APPLE_GNU) && defined(HAVE_EXX_OPENBLAS_THREADS)
      !$omp end critical(hsex_q_vendor_serial)
#endif
    enddo
    !$omp end do
    deallocate(b,rot,x,proj,gram,diag,dl,du,du2,rd,piv)
    if(allocated(real_x))deallocate(real_x,defect,correction,real_low,old_gram)
#ifdef HSEX_ACCELERATE
    api_status=blas_set(old_mode)
    if(api_status/=0)then
      !$omp critical
      status=-1
      !$omp end critical
    endif
#endif
#ifdef HSEX_SALMON
    call exx_blas_thread_control(0,blas_worker)
#endif
    !$omp end parallel
#ifdef HSEX_SALMON
    call exx_blas_thread_control(0,blas_parent)
#endif
  end subroutine hsex_q_solve
end module hsex_q_response
