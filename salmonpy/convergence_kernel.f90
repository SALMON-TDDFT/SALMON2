! SALMONpy convergence diagnostics. All IEEE inquiries are scalar.
module salmonpy_convergence_kernel
  use iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
contains
  integer(c_int) function sp_convergence_check(n,values,positive) bind(c)
    integer(c_int),value :: n,positive
    real(c_double),intent(in) :: values(*)
    integer :: i
    sp_convergence_check=4
    if(n<1)return
    do i=1,n
      if(.not.ieee_is_finite(values(i)))return
      if(positive/=0.and.values(i)<=0d0)return
    enddo
    sp_convergence_check=0
  end function

  integer(c_int) function sp_convergence_norm(n,values,reference,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: values(*)
    real(c_double),value :: reference,atol,rtol
    real(c_double),intent(out) :: stats(6)
    integer(c_int),intent(out) :: accepted
    real(c_double) :: loss,gain,limit,delta
    integer :: i
    sp_convergence_norm=4;accepted=0;stats=0d0
    if(.not.ieee_is_finite(reference))return
    if(reference<=0d0)return
    if(sp_convergence_check(n,values,1_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    loss=0d0;gain=0d0
    do i=1,n
      delta=reference-values(i)
      if(.not.ieee_is_finite(delta))return
      loss=max(loss,delta);gain=max(gain,-delta)
    enddo
    stats(1)=loss;stats(2)=gain;stats(3)=loss/reference;stats(4)=gain/reference
    stats(5)=max(loss,gain);stats(6)=stats(5)/reference
    do i=1,6
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    limit=atol+rtol*reference
    if(.not.ieee_is_finite(limit))return
    if(stats(5)<=limit)accepted=1
    sp_convergence_norm=0
  end function

  integer(c_int) function sp_convergence_compare(n,coarse,fine,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: coarse(*),fine(*)
    real(c_double),value :: atol,rtol
    real(c_double),intent(out) :: stats(2)
    integer(c_int),intent(out) :: accepted
    integer :: i
    real(c_double) :: delta,scale,limit
    sp_convergence_compare=4;stats=0d0;accepted=1;scale=0d0
    if(sp_convergence_check(n,coarse,0_c_int)/=0)return
    if(sp_convergence_check(n,fine,0_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    do i=1,n
      delta=abs(coarse(i)-fine(i));limit=atol+rtol*abs(fine(i))
      if(.not.ieee_is_finite(delta).or..not.ieee_is_finite(limit))return
      stats(1)=max(stats(1),delta);scale=max(scale,abs(fine(i)))
      if(delta>limit)accepted=0
    enddo
    if(scale>0d0)then
      stats(2)=stats(1)/scale
      if(.not.ieee_is_finite(stats(2)))stats(2)=huge(1d0)
    else if(stats(1)>0d0)then
      stats(2)=huge(1d0)
    endif
    sp_convergence_compare=0
  end function
  subroutine scaled_square(value,scale,ssq)
    real(c_double),intent(in) :: value
    real(c_double),intent(inout) :: scale,ssq
    real(c_double) :: magnitude
    magnitude=abs(value)
    if(magnitude==0d0)return
    if(scale<magnitude)then
      ssq=1d0+ssq*(scale/magnitude)**2
      scale=magnitude
    else
      ssq=ssq+(magnitude/scale)**2
    endif
  end subroutine

  integer(c_int) function sp_convergence_waveforms(n,coarse,fine,atol,rtol,stats,accepted) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: coarse(*),fine(*)
    real(c_double),value :: atol,rtol
    real(c_double),intent(out) :: stats(4)
    integer(c_int),intent(out) :: accepted
    real(c_double) :: difference,es,eq,rs,rq,limit
    integer :: i
    sp_convergence_waveforms=4;stats=0d0;accepted=0
    if(sp_convergence_check(n,coarse,0_c_int)/=0)return
    if(sp_convergence_check(n,fine,0_c_int)/=0)return
    if(.not.ieee_is_finite(atol).or..not.ieee_is_finite(rtol))return
    if(atol<0d0.or.rtol<0d0)return
    es=0d0;eq=1d0;rs=0d0;rq=1d0
    do i=1,n
      difference=coarse(i)-fine(i)
      if(.not.ieee_is_finite(difference))return
      call scaled_square(difference,es,eq)
      call scaled_square(fine(i),rs,rq)
      stats(4)=max(stats(4),abs(difference))
    enddo
    stats(1)=es*sqrt(eq/dble(n))
    stats(3)=rs*sqrt(rq/dble(n))
    if(stats(3)>0d0)then
      stats(2)=stats(1)/stats(3)
      if(.not.ieee_is_finite(stats(2)))stats(2)=huge(1d0)
    endif
    do i=1,4
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    limit=atol+rtol*stats(3)
    ! Inclusive boundary, allowing only machine rounding, not a physical margin.
    limit=limit+8d0*epsilon(1d0)*max(stats(1),stats(3))
    if(.not.ieee_is_finite(limit))return
    if(stats(3)>0d0.or.atol>0d0)then
      if(stats(1)<=limit)accepted=1
    endif
    sp_convergence_waveforms=0
  end function
  integer(c_int) function sp_optical_ratio(n,omega,j,e,threshold,sigma,eps,mask) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: omega(*)
    complex(c_double_complex),intent(in) :: j(*),e(*)
    real(c_double),value :: threshold
    complex(c_double_complex),intent(out) :: sigma(*),eps(*)
    integer(c_int),intent(out) :: mask(*)
    real(c_double) :: peak,cutoff,nan
    complex(c_double_complex) :: z,d
    integer :: i
    sp_optical_ratio=4
    if(n<1)return
    if(.not.ieee_is_finite(threshold))return
    if(threshold<0d0.or.threshold>=1d0)return
    peak=0d0
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      if(.not.ieee_is_finite(real(j(i))).or..not.ieee_is_finite(aimag(j(i))))return
      if(.not.ieee_is_finite(real(e(i))).or..not.ieee_is_finite(aimag(e(i))))return
      if(.not.ieee_is_finite(abs(e(i))))return
      peak=max(peak,abs(e(i)))
    enddo
    nan=ieee_value(0d0,ieee_quiet_nan)
    cutoff=threshold*peak
    do i=1,n
      sigma(i)=cmplx(nan,nan,kind=c_double);eps(i)=sigma(i);mask(i)=0
      if(abs(e(i))<=cutoff.or.abs(e(i))==0d0)cycle
      z=j(i)/e(i)
      if(.not.ieee_is_finite(real(z)).or..not.ieee_is_finite(aimag(z)))cycle
      sigma(i)=z;mask(i)=1
      if(omega(i)==0d0)cycle
      d=1d0+cmplx(0d0,4d0*acos(-1d0),kind=c_double)*z/omega(i)
      if(.not.ieee_is_finite(real(d)).or..not.ieee_is_finite(aimag(d)))cycle
      eps(i)=d;mask(i)=3
    enddo
    sp_optical_ratio=0
  end function
  integer(c_int) function sp_optical_diagnostic(n,omega,a,b,ma,mb,stats) bind(c)
    integer(c_int),value :: n
    real(c_double),intent(in) :: omega(*)
    complex(c_double_complex),intent(in) :: a(*),b(*)
    integer(c_int),intent(in) :: ma(*),mb(*)
    real(c_double),intent(out) :: stats(4)
    real(c_double) :: es,eq,rs,rq,delta
    integer :: i,count
    sp_optical_diagnostic=4;stats=0d0;es=0d0;eq=1d0;rs=0d0;rq=1d0;count=0
    if(n<1)return
    do i=1,n
      if(.not.ieee_is_finite(omega(i)))return
      if(i>1)then
        if(omega(i)<=omega(i-1))return
      endif
      if(ma(i)/=0.and.mb(i)/=0)then
        if(.not.ieee_is_finite(real(a(i))).or..not.ieee_is_finite(aimag(a(i))))return
        if(.not.ieee_is_finite(real(b(i))).or..not.ieee_is_finite(aimag(b(i))))return
        delta=abs(a(i)-b(i))
        if(.not.ieee_is_finite(delta))return
        call scaled_square(delta,es,eq)
        call scaled_square(abs(a(i)),rs,rq)
        count=count+1
      endif
      if(i>1)then
        if(ma(i)/=0.and.ma(i-1)/=0)then
          stats(4)=stats(4)+.5d0*(real(a(i))+real(a(i-1)))*(omega(i)-omega(i-1))
        endif
      endif
    enddo
    stats(1)=dble(count)
    if(count>0)then
      stats(2)=es*sqrt(eq/dble(count));stats(3)=rs*sqrt(rq/dble(count))
    endif
    do i=1,4
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    sp_optical_diagnostic=0
  end function
  ! Observation-only PRB85,045134 Appendix A4 regularization. No solver state.
  integer(c_int) function sp_optical_a4(n,m,kind,absorption,omega,input,valid,diagonal,low,high,cutoff, &
      raw_sigma,raw_eps,sigma,eps,raw_bits,bits,applied,coeff,stats,counts)bind(c)result(status)
    integer(c_int),value::n,m,kind,absorption
    real(c_double),intent(in)::omega(n)
    complex(c_double_complex),intent(in)::input(m,n)
    integer(c_int),intent(in)::valid(m,n),diagonal(m)
    real(c_double),value::low,high,cutoff
    complex(c_double_complex),intent(out)::raw_sigma(m,n),raw_eps(m,n),sigma(m,n),eps(m,n)
    integer(c_int),intent(out)::raw_bits(m,n),bits(m,n),applied(m,n),counts(3,m)
    real(c_double),intent(out)::coeff(4,m),stats(4,m)
    real(c_double)::nan,pi4,w,x,y,z,meanx,meany,meanz,sxx,sxy,sxz,dx,dy,dz
    real(c_double)::slope,intercept,abs_slope,abs_intercept,er,sr,si,ei,delta,ss,scale
    complex(c_double_complex)::s,e
    integer::i,j,k
    logical::fitpoint(n),finite_value
    status=1
    if(n<1.or.m<1.or.(kind/=0.and.kind/=1).or.(absorption/=0.and.absorption/=1))return
    if(.not.ieee_is_finite(low).or..not.ieee_is_finite(high).or..not.ieee_is_finite(cutoff))return
    if(low<=0d0.or.high<=low.or.cutoff<=0d0)return
    do i=1,n
      if(.not.ieee_is_finite(omega(i)).or.omega(i)<0d0)return
      if(i>1)then
        if(omega(i)<=omega(i-1))return
      endif
    enddo
    do j=1,m
      if(diagonal(j)/=0.and.diagonal(j)/=1)return
    enddo
    nan=ieee_value(0d0,ieee_quiet_nan);pi4=4d0*acos(-1d0)
    raw_sigma=cmplx(nan,nan,c_double);raw_eps=raw_sigma
    ! Keep supplied raw values even at masked bins; never treat them as data.
    if(kind==0)then
      raw_sigma=input
    else
      raw_eps=input
    endif
    raw_bits=0;bits=0;applied=0;coeff=0d0;stats=0d0;counts=0
    do i=1,n
      w=omega(i)
      do j=1,m
        if(valid(j,i)==0)cycle
        if(.not.ieee_is_finite(real(input(j,i),c_double)))return
        if(.not.ieee_is_finite(aimag(input(j,i))))return
        if(kind==0)then
          raw_sigma(j,i)=input(j,i);raw_bits(j,i)=1
          if(w==0d0)cycle
          e=real(diagonal(j),c_double)+cmplx(0d0,pi4,c_double)*input(j,i)/w
          if(ieee_is_finite(real(e,c_double)).and.ieee_is_finite(aimag(e)))then
            raw_eps(j,i)=e;raw_bits(j,i)=3
          endif
        else
          raw_eps(j,i)=input(j,i);raw_bits(j,i)=2
          if(w==0d0)cycle
          s=cmplx(0d0,-w/pi4,c_double)*(input(j,i)-real(diagonal(j),c_double))
          if(ieee_is_finite(real(s,c_double)).and.ieee_is_finite(aimag(s)))then
            raw_sigma(j,i)=s;raw_bits(j,i)=3
          endif
        endif
      enddo
    enddo
    sigma=raw_sigma;eps=raw_eps;bits=raw_bits
    do j=1,m
      k=0;meanx=0d0;meany=0d0;meanz=0d0;sxx=0d0;sxy=0d0;sxz=0d0;fitpoint=.false.
      do i=1,n
        w=omega(i)
        if(w<low.or.w>high.or.iand(raw_bits(j,i),2)==0)cycle
        x=(w/high)**2;y=real(raw_eps(j,i),c_double);z=0d0
        if(absorption==1)then
          if(iand(raw_bits(j,i),1)==0)cycle
          z=(real(raw_sigma(j,i),c_double)/w)/w
          if(.not.ieee_is_finite(z))cycle
        endif
        fitpoint(i)=.true.;k=k+1
        dx=x-meanx;dy=y-meany;dz=z-meanz
        meanx=meanx+dx/real(k,c_double);meany=meany+dy/real(k,c_double);meanz=meanz+dz/real(k,c_double)
        sxx=sxx+dx*(x-meanx);sxy=sxy+dx*(y-meany);sxz=sxz+dx*(z-meanz)
      enddo
      counts(1,j)=k
      if(k<2.or.sxx<=tiny(1d0))then
        counts(2,j)=1;cycle
      endif
      slope=sxy/sxx;intercept=meany-slope*meanx
      abs_slope=sxz/sxx;abs_intercept=meanz-abs_slope*meanx
      coeff(:,j)=[intercept,(slope/high)/high,abs_intercept,(abs_slope/high)/high]
      stats(4,j)=(real(diagonal(j),c_double)-intercept)/pi4
      finite_value=.true.
      do k=1,4
        if(.not.ieee_is_finite(coeff(k,j)))finite_value=.false.
      enddo
      if(.not.finite_value)then
        counts(2,j)=2;cycle
      endif
      ss=0d0;scale=0d0
      do i=1,n
        if(.not.fitpoint(i))cycle
        delta=real(raw_eps(j,i),c_double)-(intercept+slope*(omega(i)/high)**2)
        call scaled_square(delta,scale,ss)
      enddo
      stats(1,j)=scale*sqrt(ss/real(counts(1,j),c_double))
      do i=1,n
        w=omega(i)
        if(w>=cutoff)cycle
        ! Unmeasured positive-frequency data/components remain unmeasured.
        if(w>0d0.and.raw_bits(j,i)==0)cycle
        x=(w/high)**2
        er=intercept+slope*x
        si=-w*(er-real(diagonal(j),c_double))/pi4
        if(absorption==1)then
          sr=w*(w*(abs_intercept+abs_slope*x))
          ei=pi4*w*(abs_intercept+abs_slope*x)
        else
          if(w==0d0)then
            ! Paper replaces real epsilon only; DC imaginary limit is not inferred.
            if(iand(raw_bits(j,i),2)/=0)then
              if(ieee_is_finite(er))then
                eps(j,i)=cmplx(er,aimag(raw_eps(j,i)),c_double);applied(j,i)=1
                counts(3,j)=counts(3,j)+1
              endif
            endif
            cycle
          endif
          if(iand(raw_bits(j,i),1)==0)cycle
          sr=real(raw_sigma(j,i),c_double);ei=aimag(raw_eps(j,i))
        endif
        if(.not.ieee_is_finite(er).or..not.ieee_is_finite(si))then
          counts(2,j)=2;cycle
        endif
        if(.not.ieee_is_finite(sr).or..not.ieee_is_finite(ei))then
          counts(2,j)=2;cycle
        endif
        sigma(j,i)=cmplx(sr,si,c_double);eps(j,i)=cmplx(er,ei,c_double);bits(j,i)=3
        applied(j,i)=1;counts(3,j)=counts(3,j)+1
        if(iand(raw_bits(j,i),2)/=0)stats(2,j)=max(stats(2,j),abs(er-real(raw_eps(j,i),c_double)))
        if(iand(raw_bits(j,i),1)/=0)stats(3,j)=max(stats(3,j),abs(sr-real(raw_sigma(j,i),c_double)))
      enddo
    enddo
    status=0
  end function sp_optical_a4
  ! Direct requested-frequency observation transform; never a solver timestep.
  integer(c_int) function sp_observation_spectrum(n,m,nw,time,data,omega,window,out,valid,stats) bind(c) result(status)
    integer(c_int),value::n,m,nw,window
    real(c_double),intent(in)::time(n),data(m,n),omega(nw)
    complex(c_double_complex),intent(out)::out(m,nw)
    integer(c_int),intent(out)::valid(nw)
    real(c_double),intent(out)::stats(4)
    real(c_double)::dt,pi,u,weight,nan
    complex(c_double_complex)::phase,advance,term,y,acc(m),comp(m),next(m)
    integer::i,j,k
    status=1;stats=0d0
    if(n<2.or.m<1.or.nw<1.or.window<0.or.window>2)return
    if(window==1.and.n<3)return
    dt=time(2)-time(1);pi=acos(-1d0)
    if(.not.ieee_is_finite(dt).or.dt<=0d0)return
    do i=1,n
      if(.not.ieee_is_finite(time(i)))return
      if(i>1)then
        if(.not.ieee_is_finite(time(i)-time(i-1)))return
        if(abs((time(i)-time(i-1))-dt)>1d-8*dt)return
      endif
      do j=1,m
        if(.not.ieee_is_finite(data(j,i)))return
      enddo
    enddo
    do k=1,nw
      if(.not.ieee_is_finite(omega(k)).or.omega(k)<0d0)return
      if(k>1)then
        if(omega(k)<=omega(k-1))return
      endif
    enddo
    stats=[dt,time(n)-time(1),2d0*pi/(dt*real(n,c_double)),pi/dt]
    do j=1,4
      if(.not.ieee_is_finite(stats(j)))return
    enddo
    nan=ieee_value(0d0,ieee_quiet_nan);out=cmplx(nan,nan,c_double);valid=0
    do k=1,nw
      if(omega(k)>stats(4))cycle
      if(.not.ieee_is_finite(omega(k)*time(n)))return
      if(.not.ieee_is_finite(omega(k)*time(1)))return
      phase=exp(cmplx(0d0,omega(k)*time(1),c_double))
      advance=exp(cmplx(0d0,omega(k)*dt,c_double));acc=0d0;comp=0d0
      do i=1,n
        if(mod(i-1,256)==0)phase=exp(cmplx(0d0,omega(k)*time(i),c_double))
        u=real(i-1,c_double)/real(n-1,c_double)
        select case(window)
        case(0);weight=1d0
        case(1);weight=.5d0-.5d0*cos(2d0*pi*u)
        case(2);weight=1d0-3d0*u*u+2d0*u*u*u
        end select
        do j=1,m
          term=(dt*weight)*data(j,i)*phase
          y=term-comp(j);next(j)=acc(j)+y;comp(j)=(next(j)-acc(j))-y;acc(j)=next(j)
        enddo
        phase=phase*advance
      enddo
      do j=1,m
        if(.not.ieee_is_finite(real(acc(j),c_double)).or..not.ieee_is_finite(aimag(acc(j))))return
      enddo
      out(:,k)=acc;valid(k)=1
    enddo
    status=0
  end function

  integer(c_int) function sp_k_current(n,ne,nk,steps,groups,raw,dt,time0,sign,time,byk,total,stats) bind(c) result(status)
    integer(c_int),value::n,ne,nk,sign
    integer(c_int64_t),intent(in)::steps(n)
    integer(c_int),intent(in)::groups(ne)
    real(c_double),intent(in)::raw(3,ne,n)
    real(c_double),value::dt,time0
    real(c_double),intent(out)::time(n),byk(3,nk,n),total(3,n),stats(2,nk)
    real(c_double)::comp(3,nk),v,y,next,es,eq,rs,rq,den,nan,ts(3)
    integer::i,e,k,a
    status=1;stats=0d0
    if(n<2.or.ne<1.or.nk<1.or.(sign/=1.and.sign/=-1))return
    if(.not.ieee_is_finite(dt).or.dt<=0d0.or..not.ieee_is_finite(time0))return
    do e=1,ne
      if(groups(e)<1.or.groups(e)>nk)return
    enddo
    byk=0d0;total=0d0
    do i=1,n
      if(i>1)then
        if(steps(i)<=steps(i-1))return
      endif
      time(i)=time0+real(steps(i),c_double)*dt
      if(.not.ieee_is_finite(time(i)))return
      if(i>1)then
        if(time(i)<=time(i-1))return
      endif
      comp=0d0
      do e=1,ne
        k=groups(e)
        do a=1,3
          v=raw(a,e,i)
          if(.not.ieee_is_finite(v))return
          y=real(sign,c_double)*v-comp(a,k);next=byk(a,k,i)+y
          comp(a,k)=(next-byk(a,k,i))-y;byk(a,k,i)=next
        enddo
      enddo
      ts=0d0
      do k=1,nk
        do a=1,3
          if(.not.ieee_is_finite(byk(a,k,i)))return
          y=byk(a,k,i)-ts(a);next=total(a,i)+y;ts(a)=(next-total(a,i))-y;total(a,i)=next
          if(.not.ieee_is_finite(total(a,i)))return
        enddo
      enddo
    enddo
    rs=0d0;rq=0d0
    do i=1,n
      do a=1,3
        call scaled_square(total(a,i),rs,rq)
      enddo
    enddo
    den=rs*sqrt(rq/(3d0*real(n,c_double)));nan=ieee_value(0d0,ieee_quiet_nan)
    do k=1,nk
      es=0d0;eq=0d0
      do i=1,n
        do a=1,3
          call scaled_square(byk(a,k,i),es,eq)
        enddo
      enddo
      stats(1,k)=es*sqrt(eq/(3d0*real(n,c_double)));stats(2,k)=nan
      if(den>0d0)then
        v=stats(1,k)/den
        if(ieee_is_finite(v))stats(2,k)=v
      endif
    enddo
    status=0
  end function
  integer(c_int) function sp_observation_difference(n,omega,a,b,ma,mb,stats) bind(c) result(status)
    integer(c_int),value::n
    real(c_double),intent(in)::omega(n)
    complex(c_double_complex),intent(in)::a(n),b(n)
    integer(c_int),intent(in)::ma(n),mb(n)
    real(c_double),intent(out)::stats(5)
    real(c_double)::ratio
    stats=0d0
    status=sp_optical_diagnostic(n,omega,a,b,ma,mb,stats(1:4))
    stats(5)=ieee_value(0d0,ieee_quiet_nan)
    if(status/=0)return
    if(stats(3)>0d0)then
      ratio=stats(2)/stats(3)
      if(ieee_is_finite(ratio))stats(5)=ratio
    endif
  end function

  ! Collocated total E and physical charge j, time-major logical arrays.
  integer(c_int) function sp_absorbed_work(nt,np,time,e,j,work) bind(c)
    integer(c_int),value::nt,np
    real(c_double),intent(in)::time(nt),e(3,np,nt),j(3,np,nt)
    real(c_double),intent(out)::work(np)
    real(c_double)::dt,previous,power,increment
    integer::t,p,k
    sp_absorbed_work=4;work=0d0
    if(nt<2.or.np<1)return
    do t=1,nt
      if(.not.ieee_is_finite(time(t)))return
      if(t>1)then
        dt=time(t)-time(t-1)
        if(.not.ieee_is_finite(dt).or.dt<=0d0)return
      endif
      do p=1,np
        do k=1,3
          if(.not.ieee_is_finite(e(k,p,t)).or..not.ieee_is_finite(j(k,p,t)))return
        enddo
      enddo
    enddo
    do p=1,np
      previous=dot_product(e(:,p,1),j(:,p,1))
      if(.not.ieee_is_finite(previous))return
      do t=2,nt
        power=dot_product(e(:,p,t),j(:,p,t))
        if(.not.ieee_is_finite(power))return
        increment=(time(t)-time(t-1))*(0.5d0*previous+0.5d0*power)
        work(p)=work(p)+increment
        if(.not.ieee_is_finite(work(p)))return
        previous=power
      enddo
    enddo
    sp_absorbed_work=0
  end function

  ! All six port spectra use the same FFT/window/time convention.
  integer(c_int) function sp_port_spectra(n,ei,hi,er,hr,et,ht,normal,rtol,r,t,valid) bind(c)
    integer(c_int),value::n
    complex(c_double_complex),intent(in)::ei(3,n),hi(3,n),er(3,n),hr(3,n),et(3,n),ht(3,n)
    real(c_double),intent(in)::normal(3)
    real(c_double),value::rtol
    real(c_double),intent(out)::r(n),t(n)
    integer(c_int),intent(out)::valid(n)
    real(c_double)::flux(n,3),peak,v(3),length
    complex(c_double_complex)::a(3),b(3)
    integer::i,k,port
    sp_port_spectra=4;valid=0
    r=ieee_value(0d0,ieee_quiet_nan);t=r
    if(n<1.or..not.ieee_is_finite(rtol))return
    if(rtol<0d0.or.rtol>=1d0)return
    do k=1,3
      if(.not.ieee_is_finite(normal(k)))return
    enddo
    length=sqrt(dot_product(normal,normal))
    if(.not.ieee_is_finite(length).or.length<=0d0)return
    do i=1,n
      do port=1,3
        select case(port)
        case(1);a=ei(:,i);b=hi(:,i)
        case(2);a=er(:,i);b=hr(:,i)
        case(3);a=et(:,i);b=ht(:,i)
        end select
        do k=1,3
          if(.not.ieee_is_finite(real(a(k))).or..not.ieee_is_finite(aimag(a(k))))return
          if(.not.ieee_is_finite(real(b(k))).or..not.ieee_is_finite(aimag(b(k))))return
        enddo
        v=[real(a(2)*conjg(b(3))-a(3)*conjg(b(2))), &
           real(a(3)*conjg(b(1))-a(1)*conjg(b(3))), &
           real(a(1)*conjg(b(2))-a(2)*conjg(b(1)))]
        flux(i,port)=dot_product(v,normal)/length
        if(.not.ieee_is_finite(flux(i,port)))return
      enddo
    enddo
    peak=maxval(flux(:,1))
    do i=1,n
      if(flux(i,1)<=0d0.or.flux(i,1)<=rtol*peak)cycle
      r(i)=-flux(i,2)/flux(i,1);t(i)=flux(i,3)/flux(i,1)
      if(.not.ieee_is_finite(r(i)).or..not.ieee_is_finite(t(i)))return
      valid(i)=1
    enddo
    sp_port_spectra=0
  end function

  integer(c_int) function sp_port_fluence(nt,time,e,h,normal,fluence) bind(c)
    integer(c_int),value::nt
    real(c_double),intent(in)::time(nt),e(3,nt),h(3,nt),normal(3)
    real(c_double),intent(out)::fluence
    real(c_double)::length,power,previous,v(3),dt
    integer::i,k
    sp_port_fluence=4;fluence=0d0
    if(nt<2)return
    do k=1,3
      if(.not.ieee_is_finite(normal(k)))return
    enddo
    length=sqrt(dot_product(normal,normal))
    if(.not.ieee_is_finite(length).or.length<=0d0)return
    previous=0d0
    do i=1,nt
      if(.not.ieee_is_finite(time(i)))return
      do k=1,3
        if(.not.ieee_is_finite(e(k,i)).or..not.ieee_is_finite(h(k,i)))return
      enddo
      v=[e(2,i)*h(3,i)-e(3,i)*h(2,i),e(3,i)*h(1,i)-e(1,i)*h(3,i), &
         e(1,i)*h(2,i)-e(2,i)*h(1,i)]
      power=137.035999084d0/(4d0*acos(-1d0))*dot_product(v,normal)/length
      if(.not.ieee_is_finite(power))return
      if(i>1)then
        dt=time(i)-time(i-1)
        if(.not.ieee_is_finite(dt).or.dt<=0d0)return
        fluence=fluence+dt*(0.5d0*previous+0.5d0*power)
        if(.not.ieee_is_finite(fluence))return
      endif
      previous=power
    enddo
    sp_port_fluence=0
  end function

  integer(c_int) function sp_energy_balance(n,work,volumes,terms,area,rtol,atol,stats,accepted) bind(c)
    integer(c_int),value::n
    real(c_double),intent(in)::work(n),volumes(n),terms(5)
    real(c_double),value::area,rtol,atol
    real(c_double),intent(out)::stats(6)
    integer(c_int),intent(out)::accepted
    real(c_double)::total,residual,limit,energy(3)
    integer::i
    sp_energy_balance=4;accepted=0;stats=0d0
    if(n<1)return
    if(.not.ieee_is_finite(rtol).or..not.ieee_is_finite(atol))return
    if(rtol<0d0.or.atol<0d0)return
    do i=1,5
      if(.not.ieee_is_finite(terms(i)))return
    enddo
    if(.not.ieee_is_finite(area).or.area<=0d0)return
    energy=terms(1:3)*area
    do i=1,3
      if(.not.ieee_is_finite(energy(i)))return
    enddo
    if(energy(1)<=0d0)return
    total=0d0
    do i=1,n
      if(.not.ieee_is_finite(work(i)).or..not.ieee_is_finite(volumes(i)))return
      if(volumes(i)<=0d0)return
      total=total+work(i)*volumes(i)
      if(.not.ieee_is_finite(total))return
    enddo
    residual=energy(1)-energy(2)-energy(3)-total-terms(4)-terms(5)
    limit=atol+rtol*energy(1)
    stats=[energy,total,residual,residual/energy(1)]
    do i=1,6
      if(.not.ieee_is_finite(stats(i)))return
    enddo
    if(.not.ieee_is_finite(limit))return
    accepted=merge(1,0,abs(residual)<=limit)
    sp_energy_balance=0
  end function

  ! Local E cross H first, area integration second: no coherent field averaging.
  integer(c_int) function surface_geometry(np,normal,weights,unit,area) result(status)
    integer,intent(in)::np
    real(c_double),intent(in)::normal(3),weights(np)
    real(c_double),intent(out)::unit(3),area
    real(c_double)::length
    integer::k
    status=4;area=0d0
    if(np<1)return
    do k=1,3
      if(.not.ieee_is_finite(normal(k)))return
    enddo
    length=sqrt(dot_product(normal,normal))
    if(.not.ieee_is_finite(length).or.length<=0d0)return
    unit=normal/length
    do k=1,np
      if(.not.ieee_is_finite(weights(k)).or.weights(k)<0d0)return
      area=area+weights(k)
    enddo
    if(.not.ieee_is_finite(area).or.area<=0d0)return
    status=0
  end function

  integer(c_int) function sp_surface_flux(nt,np,time,e,h,normal,weights,power,energy,stats) bind(c)
    integer(c_int),value::nt,np
    real(c_double),intent(in)::time(nt),e(3,np,nt),h(3,np,nt),normal(3),weights(np)
    real(c_double),intent(out)::power(nt),energy(nt),stats(3)
    real(c_double)::unit(3),area,v(3),flux,dt
    integer::t,p,k
    sp_surface_flux=4;power=0d0;energy=0d0;stats=0d0
    if(nt<2)return
    if(surface_geometry(np,normal,weights,unit,area)/=0)return
    do t=1,nt
      if(.not.ieee_is_finite(time(t)))return
      flux=0d0
      do p=1,np
        do k=1,3
          if(.not.ieee_is_finite(e(k,p,t)).or..not.ieee_is_finite(h(k,p,t)))return
        enddo
        v=[e(2,p,t)*h(3,p,t)-e(3,p,t)*h(2,p,t), &
           e(3,p,t)*h(1,p,t)-e(1,p,t)*h(3,p,t), &
           e(1,p,t)*h(2,p,t)-e(2,p,t)*h(1,p,t)]
        flux=flux+weights(p)*dot_product(v,unit)
        if(.not.ieee_is_finite(flux))return
      enddo
      power(t)=137.035999084d0/(4d0*acos(-1d0))*flux
      if(.not.ieee_is_finite(power(t)))return
      if(t>1)then
        dt=time(t)-time(t-1)
        if(.not.ieee_is_finite(dt).or.dt<=0d0)return
        energy(t)=energy(t-1)+dt*(.5d0*power(t-1)+.5d0*power(t))
        if(.not.ieee_is_finite(energy(t)))return
      endif
    enddo
    stats=[area,energy(nt),energy(nt)/area]
    do k=1,3
      if(.not.ieee_is_finite(stats(k)))return
    enddo
    sp_surface_flux=0
  end function

  integer(c_int) function sp_surface_port_spectra(n,np,ei,hi,er,hr,et,ht,normal,weights, &
                                                 rtol,r,t,valid,area) bind(c)
    integer(c_int),value::n,np
    complex(c_double_complex),intent(in)::ei(3,np,n),hi(3,np,n),er(3,np,n),hr(3,np,n)
    complex(c_double_complex),intent(in)::et(3,np,n),ht(3,np,n)
    real(c_double),intent(in)::normal(3),weights(np)
    real(c_double),value::rtol
    real(c_double),intent(out)::r(n),t(n),area
    integer(c_int),intent(out)::valid(n)
    real(c_double)::unit(3),flux(n,3),v(3),peak
    complex(c_double_complex)::a(3),b(3)
    integer::i,p,k,port
    sp_surface_port_spectra=4;valid=0
    r=ieee_value(0d0,ieee_quiet_nan);t=r
    if(n<1.or..not.ieee_is_finite(rtol))return
    if(rtol<0d0.or.rtol>=1d0)return
    if(surface_geometry(np,normal,weights,unit,area)/=0)return
    flux=0d0
    do i=1,n
      do port=1,3
        do p=1,np
          select case(port)
          case(1);a=ei(:,p,i);b=hi(:,p,i)
          case(2);a=er(:,p,i);b=hr(:,p,i)
          case(3);a=et(:,p,i);b=ht(:,p,i)
          end select
          do k=1,3
            if(.not.ieee_is_finite(real(a(k))).or..not.ieee_is_finite(aimag(a(k))))return
            if(.not.ieee_is_finite(real(b(k))).or..not.ieee_is_finite(aimag(b(k))))return
          enddo
          v=[real(a(2)*conjg(b(3))-a(3)*conjg(b(2))), &
             real(a(3)*conjg(b(1))-a(1)*conjg(b(3))), &
             real(a(1)*conjg(b(2))-a(2)*conjg(b(1)))]
          flux(i,port)=flux(i,port)+weights(p)*dot_product(v,unit)
          if(.not.ieee_is_finite(flux(i,port)))return
        enddo
      enddo
    enddo
    peak=maxval(flux(:,1))
    do i=1,n
      if(flux(i,1)<=0d0.or.flux(i,1)<=rtol*peak)cycle
      r(i)=-flux(i,2)/flux(i,1);t(i)=flux(i,3)/flux(i,1)
      if(.not.ieee_is_finite(r(i)).or..not.ieee_is_finite(t(i)))return
      valid(i)=1
    enddo
    sp_surface_port_spectra=0
  end function
end module
