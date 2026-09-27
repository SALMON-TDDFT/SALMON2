program probe
  use dc_thermal
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(8) :: e(5),w(5),f(5),fp(5),fm(5),de(5),dw(5),df(5),mu,mup,mum,ts,tsp,tsm,dmu,dts,t,n,g,h
  integer :: status,mode
  e=[-.8d0,-.1d0,.2d0,.2d0,1.1d0];w=[.3d0,.8d0,.4d0,.7d0,0d0]
  t=.17d0;g=2d0;n=2.3d0;h=1d-5
  call solve_dc_thermal(e,w,t,g,n,mu,f,ts,status)
  if(status/=0)error stop 'weighted solve failed'
  if(abs(g*sum(w*f)-n)>1d-10)error stop 'charge constraint'
  if(abs(ts+g*t*sum(w*(f*log(f)+(1-f)*log(1-f))))>1d-13)error stop 'entropy'
  if(abs(f(3)-f(4))>1d-14)error stop 'degenerate occupations'
  do mode=1,3
    de=[.21d0,-.32d0,.43d0,.17d0,-.25d0];dw=[-.1d0,.2d0,.3d0,-.2d0,0d0]
    if(mode==1)dw=0d0
    if(mode==2)de=0d0
    call response_dc_thermal(e,w,t,g,mu,de,dw,dmu,df,dts,status)
    if(status/=0)error stop 'response failed'
    if(abs(g*sum(w*df+f*dw))>1d-12)error stop 'differentiated charge'
    call solve_dc_thermal(e+h*de,w+h*dw,t,g,n,mup,fp,tsp,status)
    if(status/=0)error stop 'plus solve'
    call solve_dc_thermal(e-h*de,w-h*dw,t,g,n,mum,fm,tsm,status)
    if(status/=0)error stop 'minus solve'
    if(abs(dmu-(mup-mum)/(2*h))>2d-7)error stop 'mu derivative'
    if(maxval(abs(df-(fp-fm)/(2*h)))>2d-7)error stop 'occupation derivative'
    if(abs(dts-(tsp-tsm)/(2*h))>2d-7)error stop 'entropy derivative'
  enddo
  call response_dc_thermal(e,w,t,g,mu,spread(1d0,1,5),spread(0d0,1,5),dmu,df,dts,status)
  if(status/=0.or.abs(dmu-1)>1d-13.or.maxval(abs(df))>1d-13.or.abs(dts)>1d-13) &
    error stop 'common energy shift response'
  call solve_dc_thermal(e+3d0,w,t,g,n,mup,fp,tsp,status)
  if(status/=0.or.abs(mup-mu-3)>1d-10.or.maxval(abs(f-fp))>1d-10)error stop 'energy gauge'
  call solve_dc_thermal(e,w,t,g,g*sum(w)+.1d0,mu,f,ts,status)
  if(status/=dc_thermal_capacity)error stop 'insufficient capacity accepted'
  call solve_dc_thermal(e,w,0d0,g,n,mu,f,ts,status)
  if(status/=dc_thermal_invalid)error stop 'zero temperature accepted'
  call solve_dc_thermal(e,-w,t,g,n,mu,f,ts,status)
  if(status/=dc_thermal_invalid)error stop 'negative weights accepted'
  call solve_dc_thermal(e,w,ieee_value(t,ieee_quiet_nan),g,n,mu,f,ts,status)
  if(status/=dc_thermal_invalid)error stop 'NaN temperature accepted'
  call solve_dc_thermal(e,w,t,g,g*sum(w),mu,f,ts,status)
  if(status/=0.or.abs(g*sum(w*f)-g*sum(w))>1d-10)error stop 'full capacity charge'
  call response_dc_thermal(e,w,t,g,mu,de,dw,dmu,df,dts,status)
  if(status/=dc_thermal_singular)error stop 'saturated response accepted'
  call solve_dc_thermal(e,w,t,g,0d0,mu,f,ts,status)
  if(status/=0.or.abs(g*sum(w*f))>1d-10)error stop 'empty charge'
  ! A small but resolvable hole density must not be snapped to full capacity.
  e=0d0;w=[1d0,0d0,0d0,0d0,0d0];n=2d0-5d-11
  call solve_dc_thermal(e,w,.1d0,g,n,mu,f,ts,status)
  if(status/=0.or.abs(g*sum(w*f)-n)>1d-13)error stop 'near-full interior charge'
  if(abs(mu-.1d0*log(n/(g-n)))>1d-4)error stop 'near-full interior mu'
  de=1d0;dw=0d0
  call response_dc_thermal(e,w,.1d0,g,mu,de,dw,dmu,df,dts,status)
  if(status/=0.or.abs(dmu-1)>1d-12)error stop 'near-full response incorrectly singular'
  e=[-1000d0,-1d0,0d0,1d0,1000d0];w=1d0
  call solve_dc_thermal(e,w,.001d0,g,5d0,mu,f,ts,status)
  if(status/=0.or.abs(mu)>1d-12.or.abs(f(3)-.5d0)>1d-12)error stop 'extreme tails'
  if(abs(ts-2d0*.001d0*log(2d0))>1d-14)error stop 'tail entropy'
  if(fermi_entropy(0d0)/=0d0.or.fermi_entropy(1d0)/=0d0)error stop 'endpoint entropy'
  print *, 'DC thermal solve and fixed-charge directional response passed'
end program
