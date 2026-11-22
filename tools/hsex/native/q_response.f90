! Compensated products/sums for original-H refinement. VOLATILE temporaries
! preserve error-free transformations even in builds enabling reassociation.
module hsex_compensated
  implicit none
  private
  public :: hsex_original_residual,hsex_precise_gram,hsex_split_matrix,hsex_refinement_residual
contains
  recursive subroutine add_product(a,b,s,c)
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

  recursive subroutine hsex_split_matrix(n,m,a,rows,bits,parts,tail)
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

  recursive subroutine hsex_refinement_residual(n,g,hparts,hdiag,nh,bits,energy,eta,x,xlow,b,r)
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

  recursive subroutine hsex_precise_gram(n,g,ld,b,x,xlow,gram) bind(C)
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
  public :: hsex_q_solve,hsex_q_abi,hsex_q_shared,hsex_q_workspace
  public :: hsex_q_prepare,hsex_q_release,hsex_q_worker_prepare,hsex_q_worker_release,hsex_q_one
  type :: hsex_q_shared
    private
    integer :: bits=0,nh=0,panel_bands=1,blas_parent(3)=0
    complex(8),allocatable :: hparts(:,:,:)
    real(8),allocatable :: hdiag(:)
  end type hsex_q_shared
  type :: hsex_q_workspace
    private
    integer :: old_mode=0,api_status=0,blas_worker(3)=0
    complex(8),allocatable :: b(:,:),rot(:,:),x(:,:),proj(:,:),gram(:,:)
    complex(8),allocatable :: wide_b(:,:),wide_rot(:,:),wide_proj(:,:)
    integer,allocatable :: source_bands(:),piv(:)
    complex(8),allocatable :: real_x(:,:),defect(:,:),correction(:,:),real_low(:,:),old_gram(:,:)
    complex(8),allocatable :: diag(:),dl(:),du(:),du2(:)
    real(8),allocatable :: rd(:)
  end type hsex_q_workspace
#ifdef HSEX_ACCELERATE
  interface
    integer(c_int) function blas_set(mode) bind(C,name='BLASSetThreading')
      use iso_c_binding, only: c_int
      implicit none
      integer(c_int),value :: mode
    end function blas_set
    integer(c_int) function blas_get() bind(C,name='BLASGetThreading')
      use iso_c_binding, only: c_int
      implicit none
    end function blas_get
  end interface
#endif

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
    type(hsex_q_shared) :: shared
    type(hsex_q_workspace) :: workspace
    integer :: q,code
    response=0d0;residual=0d0;status=0;team=1
    call hsex_q_prepare(n,nb,mg,nu,pump,refine,eta,original_h,shared,code)
    if(code/=0)then
      status=code
      call hsex_q_release(shared)
      return
    endif
    !$omp parallel num_threads(nthreads) default(shared) private(workspace,q,code)
    !$omp single
!$  team=omp_get_num_threads()
    !$omp end single
    call hsex_q_worker_prepare(n,nb,mg,nu,pump,refine,shared,workspace)
    !$omp do schedule(dynamic,1)
    do q=1,nq
      call hsex_q_one(n,nb,mg,nu,nc,pump,refine,adj,d,e,u,psi(:,:,q),phase(:,:,q),counts(q), &
        levels(:,q),occ(:,q),eta,shared,workspace,response(:,:,:,q),residual(q),status(q))
    enddo
    !$omp end do
    call hsex_q_worker_release(workspace,code)
    if(code/=0)then
      !$omp critical(hsex_q_policy_error)
      status=-1
      !$omp end critical(hsex_q_policy_error)
    endif
    !$omp end parallel
    call hsex_q_release(shared)
  end subroutine hsex_q_solve

  subroutine hsex_q_prepare(n,nb,mg,nu,pump,refine,eta,original_h,shared,status)
    use hsex_compensated, only: hsex_split_matrix
#ifdef HSEX_SALMON
    use exx_blas_threads, only: exx_blas_thread_control
#endif
    implicit none
    integer,intent(in) :: n,nb,mg,nu,pump,refine
    real(8),intent(in) :: eta(:)
    complex(8),intent(in) :: original_h(:,:)
    type(hsex_q_shared),intent(out) :: shared
    integer,intent(out) :: status
    complex(8),allocatable :: hoff(:,:)
    integer :: i,nslices
    status=0
    shared%bits=(50-ceiling(log(dble(max(n,2)))/log(2d0)))/2
    shared%panel_bands=max(1,min(nb,4,int(64d0*1024**2/(48d0*max(1,n)*max(1,mg)))))
    if(refine/=0)then
      nslices=5+ceiling(max(0d0,2d0*(log(5d-4)-log(minval(eta)))/log(2d0))/shared%bits)
      if(nslices>12)then
        status=5
        return
      endif
      allocate(shared%hparts(n,n,nslices),shared%hdiag(n),hoff(n,n))
      hoff=original_h
      do i=1,n
        shared%hdiag(i)=real(original_h(i,i),8)
        hoff(i,i)=cmplx(0d0,aimag(original_h(i,i)),8)
      enddo
      call hsex_split_matrix(n,n,hoff,.true.,shared%bits,shared%hparts)
      deallocate(hoff)
      do i=1,nslices
        if(any(shared%hparts(:,:,i)/=cmplx(0d0,0d0,8)))shared%nh=i
      enddo
    endif
#if defined(HSEX_SALMON) && defined(HAVE_EXX_OPENBLAS_THREADS)
    call exx_blas_thread_control(1,shared%blas_parent)
#endif
  end subroutine hsex_q_prepare

  subroutine hsex_q_release(shared)
#ifdef HSEX_SALMON
    use exx_blas_threads, only: exx_blas_thread_control
#endif
    implicit none
    type(hsex_q_shared),intent(inout) :: shared
#ifdef HSEX_SALMON
    call exx_blas_thread_control(0,shared%blas_parent)
#endif
    if(allocated(shared%hparts))deallocate(shared%hparts,shared%hdiag)
  end subroutine hsex_q_release

  subroutine hsex_q_worker_prepare(n,nb,mg,nu,pump,refine,shared,workspace)
#ifdef HSEX_SALMON
    use exx_blas_threads, only: exx_blas_thread_control
#endif
    implicit none
    integer,intent(in) :: n,nb,mg,nu,pump,refine
    type(hsex_q_shared),intent(in) :: shared
    type(hsex_q_workspace),intent(out) :: workspace
    integer :: panel_bands
    panel_bands=shared%panel_bands
#ifdef HSEX_SALMON
    call exx_blas_thread_control(1,workspace%blas_worker)
#endif
#ifdef HSEX_ACCELERATE
    workspace%old_mode=blas_get()
    workspace%api_status=blas_set(1)
#endif
    allocate(workspace%b(n,mg),workspace%rot(n,mg),workspace%x(n,mg),workspace%proj(max(1,nu),mg), &
      workspace%gram(mg,mg),workspace%source_bands(nb))
    if(pump==0)allocate(workspace%wide_b(n,mg*panel_bands),workspace%wide_rot(n,mg*panel_bands), &
      workspace%wide_proj(max(1,nu),mg*panel_bands))
    if(refine/=0)allocate(workspace%real_x(n,mg),workspace%defect(n,mg),workspace%correction(n,mg), &
      workspace%real_low(n,mg),workspace%old_gram(mg,mg))
    allocate(workspace%diag(n),workspace%dl(n-1),workspace%du(n-1),workspace%du2(max(1,n-2)), &
      workspace%rd(n),workspace%piv(n))
  end subroutine hsex_q_worker_prepare

  subroutine hsex_q_worker_release(workspace,status)
#ifdef HSEX_SALMON
    use exx_blas_threads, only: exx_blas_thread_control
#endif
    implicit none
    type(hsex_q_workspace),intent(inout) :: workspace
    integer,intent(out) :: status
    status=0
#ifdef HSEX_ACCELERATE
    if(blas_set(workspace%old_mode)/=0)status=-1
#endif
#ifdef HSEX_SALMON
    call exx_blas_thread_control(0,workspace%blas_worker)
#endif
    deallocate(workspace%b,workspace%rot,workspace%x,workspace%proj,workspace%gram, &
      workspace%source_bands,workspace%diag,workspace%dl,workspace%du,workspace%du2,workspace%rd,workspace%piv)
    if(allocated(workspace%wide_b))deallocate(workspace%wide_b,workspace%wide_rot,workspace%wide_proj)
    if(allocated(workspace%real_x))deallocate(workspace%real_x,workspace%defect,workspace%correction, &
      workspace%real_low,workspace%old_gram)
  end subroutine hsex_q_worker_release

  recursive subroutine hsex_q_one(n,nb,mg,nu,nc,pump,refine,adj,d,e,u,psi,phase,g,levels,occ,eta, &
      shared,workspace,response,residual,status)
    use ieee_arithmetic, only: ieee_is_finite
    use iso_c_binding, only: c_double
    use hsex_compensated, only: hsex_precise_gram,hsex_refinement_residual
    implicit none
    integer,intent(in) :: n,nb,mg,nu,nc,pump,refine,g
    complex(8),intent(in) :: adj(n,n),u(n,nu),psi(n,nb),phase(n,mg)
    real(8),intent(in) :: d(n),e(n-1),levels(nb),occ(nb),eta(nc)
    type(hsex_q_shared),intent(in) :: shared
    type(hsex_q_workspace),intent(inout) :: workspace
    complex(8),intent(out) :: response(mg,mg,nc)
    real(8),intent(out) :: residual
    integer,intent(out) :: status
    integer :: v,c,i,j,info,iteration,active_bands,band_index,panel_first,panel_last,width,column,slot,ldg,panel_bands
    real(8) :: rn,bn,err,floor_norm
    complex(8) :: r
    complex(8),volatile :: previous,t,z,vv,ww
    logical :: converged
    response=0d0;residual=0d0;status=0
    if(workspace%api_status/=0)then
      status=-1
      return
    endif
    ldg=size(workspace%gram,1);panel_bands=shared%panel_bands
    if(g<1.or.g>mg.or.mg>ldg)then
      status=-1
      return
    endif
    associate(b=>workspace%b,rot=>workspace%rot,x=>workspace%x,proj=>workspace%proj,gram=>workspace%gram, &
      source_bands=>workspace%source_bands, &
      diag=>workspace%diag,dl=>workspace%dl,du=>workspace%du,du2=>workspace%du2,rd=>workspace%rd,piv=>workspace%piv)
#if defined(HSEX_SALMON) && defined(HAVE_HSEX_APPLE_GNU) && defined(HAVE_EXX_OPENBLAS_THREADS)
    !$omp critical(hsex_q_vendor_serial)
#endif
      active_bands=0
      do v=1,nb
        if(occ(v)==0d0)cycle
        active_bands=active_bands+1;source_bands(active_bands)=v
      enddo
      do band_index=1,active_bands
        v=source_bands(band_index)
        if(pump==0)then
          associate(wide_b=>workspace%wide_b,wide_rot=>workspace%wide_rot,wide_proj=>workspace%wide_proj)
          if(mod(band_index-1,panel_bands)==0)then
            panel_first=band_index;panel_last=min(active_bands,panel_first+panel_bands-1)
            width=g*(panel_last-panel_first+1)
            do slot=panel_first,panel_last
              column=g*(slot-panel_first)
              do j=1,g
                do i=1,n
                  wide_b(i,column+j)=psi(i,source_bands(slot))*phase(i,j)
                enddo
              enddo
            enddo
            if(nu>0)then
              call zgemm('C','N',nu,width,n,(1d0,0d0),u,n,wide_b,n,(0d0,0d0),wide_proj,max(1,nu))
              call zgemm('N','N',n,width,nu,(-1d0,0d0),u,n,wide_proj,max(1,nu),(1d0,0d0),wide_b,n)
            endif
            call zgemm('N','N',n,width,n,(1d0,0d0),adj,n,wide_b,n,(0d0,0d0),wide_rot,n)
          endif
          column=g*(band_index-panel_first)
          rot(:,:g)=wide_rot(:,column+1:column+g)
          end associate
        else
          ! Keep the original-H refinement's projected RHS and arithmetic.
          do j=1,g
            do i=1,n
              b(i,j)=psi(i,v)*phase(i,j)
            enddo
          enddo
          if(nu>0)then
            call zgemm('C','N',nu,g,n,(1d0,0d0),u,n,b,n,(0d0,0d0),proj,max(1,nu))
            call zgemm('N','N',n,g,nu,(-1d0,0d0),u,n,proj,max(1,nu),(1d0,0d0),b,n)
          endif
          call zgemm('N','N',n,g,n,(1d0,0d0),adj,n,b,n,(0d0,0d0),rot,n)
        endif
        do c=1,nc
          x(:,:g)=rot(:,:g)
          if(pump==0)then
            rd=d-levels(v);dl=cmplx(e,0d0,c_double)
            call zptsv(n,g,rd,dl,x,n,info)
            floor_norm=1d-12
          else
            diag=cmplx(d-levels(v),-eta(c),c_double)
            dl=cmplx(e,0d0,c_double);du=dl
            call zgttrf(n,dl,diag,du,du2,piv,info)
            if(info==0)call zgttrs('N',n,g,dl,diag,du,du2,piv,x,n,info)
            floor_norm=1d-14
          endif
          if(info/=0)then
            status=1
            cycle
          endif
          do j=1,g
            rn=0d0;bn=0d0
            do i=1,n
              r=cmplx(d(i)-levels(v),-eta(c),c_double)*x(i,j)-rot(i,j)
              if(i>1)r=r+e(i-1)*x(i-1,j)
              if(i<n)r=r+e(i)*x(i+1,j)
              rn=rn+abs(r)**2;bn=bn+abs(rot(i,j))**2
            enddo
            err=sqrt(rn)/max(sqrt(bn),floor_norm)
            if(.not.ieee_is_finite(err))then
              status=2
            else if(refine==0.or.pump==0)then
              residual=max(residual,err)
              if(err>1d-8)status=2
            endif
          enddo
          ! With refinement, the final ORIGINAL-H residual below is decisive;
          ! the initial transformed solve is only a preconditioned estimate.
          if(refine/=0.and.pump/=0)then
            associate(real_x=>workspace%real_x,defect=>workspace%defect,correction=>workspace%correction, &
              real_low=>workspace%real_low,old_gram=>workspace%old_gram)
            call zgemm('C','N',n,g,n,(1d0,0d0),adj,n,x,n,(0d0,0d0),real_x,n)
            real_low=0d0;old_gram=0d0;converged=.false.
            do iteration=1,10
              call hsex_refinement_residual(n,g,shared%hparts,shared%hdiag,shared%nh,shared%bits, &
                levels(v),eta(c),real_x,real_low,b,defect)
              call zgemm('N','N',n,g,n,(1d0,0d0),adj,n,defect,n,(0d0,0d0),correction,n)
              call zgttrs('N',n,g,dl,diag,du,du2,piv,correction,n,info)
              if(info/=0)then
                status=1
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
              call hsex_precise_gram(n,g,ldg,b,real_x,real_low,gram)
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
            if(.not.converged.and.status==0)status=4
            call hsex_refinement_residual(n,g,shared%hparts,shared%hdiag,shared%nh,shared%bits, &
                levels(v),eta(c),real_x,real_low,b,defect)
            do j=1,g
              rn=sum(abs(defect(:,j))**2);bn=sum(abs(b(:,j))**2)
              err=sqrt(rn)/max(sqrt(bn),floor_norm)
              if(.not.ieee_is_finite(err))then
                status=2
              else
                residual=max(residual,err)
                if(err>1d-8)status=2
              endif
            enddo
            end associate
          else
            call zgemm('C','N',g,g,n,(1d0,0d0),rot,n,x,n,(0d0,0d0),gram,ldg)
          endif
          do j=1,g
            do i=1,g
              r=gram(i,j)
              if(pump/=0)r=.5d0*(r+conjg(gram(j,i)))
              response(i,j,c)=response(i,j,c)+.5d0*occ(v)*r
              if(.not.ieee_is_finite(real(response(i,j,c),c_double)))status=3
              if(.not.ieee_is_finite(aimag(response(i,j,c))))status=3
            enddo
          enddo
        enddo
      enddo
#if defined(HSEX_SALMON) && defined(HAVE_HSEX_APPLE_GNU) && defined(HAVE_EXX_OPENBLAS_THREADS)
    !$omp end critical(hsex_q_vendor_serial)
#endif
    end associate
  end subroutine hsex_q_one

end module hsex_q_response
