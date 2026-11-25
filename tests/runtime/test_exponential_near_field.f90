! Diagnostic only. All sampling, fitting, rate selection and errors are Fortran.
! Uses separate regular-grid fit points, off-grid tuning points and final test points.
program test_exponential_near_field
 use iso_c_binding,only:c_double
 use salmonpy_solid_harmonics,only:near_basis
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_value,ieee_quiet_nan,ieee_positive_inf
 implicit none
 integer,parameter::nt=729,ns=129,nv=257
 integer,parameter::degrees(7)=[0,1,2,3,4,6,8]
 real(8),parameter::rates(5)=[.5d0,1.7d0,2d0,7.3d0,8d0]
 real(8)::x(3,nt),s(3,ns),q(3,nv),f(3,nt),t(3,ns),truth(3,nv),pred(3,nv),tune(3,ns)
 real(8)::rate,chosen,score,best,rel,absolute,worst,randoms(3,max(ns,nv)),bad(3,nt)
 real(8)::sentinel(3,nv)
 integer::i,j,k,n,profile,degree,degree_index,rindex,candidate,status,u,seedsize
 integer,allocatable::seed(:)
 character(2048)::output
 call get_command_argument(1,output)
 if(len_trim(output)==0)error stop 'provide a new CSV output path'
 open(newunit=u,file=trim(output),status='new',action='write')
 write(u,'(a)')'profile,decay_rate,total_degree,coefficients_per_component,method,used_rate,relative_rms,absolute_rms,max_absolute'
 n=0
 do k=0,8
 do j=0,8
 do i=0,8
  n=n+1;x(:,n)=[real(i,8),real(j,8),real(k,8)]/8d0-.5d0
 enddo
 enddo
 enddo
 call random_seed(size=seedsize);allocate(seed(seedsize));seed=1952;call random_seed(put=seed)
 call random_number(randoms);s=randoms(:,1:ns)-.5d0
 call random_number(randoms);q=randoms(:,1:nv)-.5d0
 ! Guard tests include scalar NaN/Inf rejection, rank loss and unchanged output.
 f=1d0;sentinel=7d0
 do i=1,2
  bad=f
  if(i==1)bad(1,1)=ieee_value(0d0,ieee_quiet_nan)
  if(i==2)bad(1,1)=ieee_value(0d0,ieee_positive_inf)
  call fit(x,bad,q,2,0d0,sentinel,status)
  if(status==0.or.any(sentinel/=7d0))error stop 'invalid field changed predictions'
 enddo
 bad=0d0
 call fit(bad,f,q,2,0d0,sentinel,status)
 if(status==0.or.any(sentinel/=7d0))error stop 'rank-deficient fit accepted'
 ! Uniform-field preservation and exact exponential recovery independent of rate tuning.
 call fit(x,f,q,0,0d0,pred,status)
 if(status/=0.or.maxval(abs(pred-1d0))>1d-12)error stop 'uniform fit changed'
 do profile=1,4
 do rindex=1,size(rates)
  if(profile==4.and.rindex>1)cycle
  rate=rates(rindex)
  if(profile==4)rate=0d0
  call field(x,profile,rate,f);call field(s,profile,rate,t);call field(q,profile,rate,truth)
  do degree_index=1,size(degrees)
   degree=degrees(degree_index)
   ! Select the rate separately at this degree; final test points are never consulted.
   best=huge(1d0);chosen=0d0
   do candidate=0,96
    call fit(x,f,s,degree,candidate*.25d0,tune,status)
    if(status/=0)cycle
    call errors(tune,t,score,absolute,worst)
    if(score<best)then
     best=score;chosen=candidate*.25d0
    endif
   enddo
   if(best==huge(1d0))error stop 'no valid exponential candidate'
   do j=1,3
    select case(j)
    case(1);score=0d0
    case(2);score=rate ! Oracle: supplied true rate, not a fitted parameter.
    case(3);score=chosen ! One additional selected parameter and tuning-set cost.
    end select
    call fit(x,f,q,degree,score,pred,status)
    if(status/=0)error stop 'experiment fit failed'
    call errors(pred,truth,rel,absolute,worst)
    if(profile==1.and.j==2.and.worst>1d-10)error stop 'known exponential not recovered'
    write(u,'(i1,a,f5.2,a,i1,a,i3,a,i1,a,4(es24.15,a))')profile,',',rate,',',degree,',', &
      (degree+1)*(degree+2)*(degree+3)/6,',',j,',',score,',',rel,',',absolute,',',worst,''
   enddo
  enddo
 enddo
 enddo
 close(u)
 print *, 'PASS: uniform, NaN/Inf/rank rejection, known exponential and independent held-out experiment'
contains
 subroutine field(points,kind,decay,values)
  real(8),intent(in)::points(:,:),decay
  integer,intent(in)::kind
  real(8),intent(out)::values(:,:)
  real(8)::z,v
  integer::ii
  values=0d0
  do ii=1,size(points,2)
   z=points(3,ii)+.5d0;v=exp(-decay*z)
   select case(kind)
   case(1) ! Pure exponential profile, not by itself a general Maxwell solution.
    values(:,ii)=v*[.3d0,-.2d0,1d0]
   case(2) ! Exact local electrostatic Laplace field; E=-grad[exp(-kz)cos(kx)/k].
    values(1,ii)=v*sin(decay*points(1,ii));values(3,ii)=v*cos(decay*points(1,ii))
   case(3) ! Two decay lengths: a stress test against the single-exponential assumption.
    values(:,ii)=(v+.3d0*exp(-3d0*decay*z))*[.3d0,-.2d0,1d0]
   case(4) ! Non-exponential control, exactly quadratic; rate hint0, not a decay model.
    values(:,ii)=(1d0+.3d0*sum(points(:,ii)**2))*[.3d0,-.2d0,1d0]
   end select
  enddo
 end subroutine
 subroutine errors(values,reference,relative,rms,maximum)
  real(8),intent(in)::values(:,:),reference(:,:)
  real(8),intent(out)::relative,rms,maximum
  real(8)::error,denom
  integer::ii,jj
  error=0d0;denom=0d0;maximum=0d0
  do ii=1,size(values,2)
  do jj=1,3
   error=error+(values(jj,ii)-reference(jj,ii))**2
   denom=denom+reference(jj,ii)**2
   maximum=max(maximum,abs(values(jj,ii)-reference(jj,ii)))
  enddo
  enddo
  relative=sqrt(error/max(denom,tiny(1d0)));rms=sqrt(error/size(values))
 end subroutine
 subroutine fit(points,values,queries,L,decay,prediction,code)
  real(8),intent(in)::points(:,:),values(:,:),queries(:,:),decay
  integer,intent(in)::L
  real(8),intent(inout)::prediction(:,:)
  integer,intent(out)::code
  real(8),allocatable::a(:,:),rhs(:,:),sv(:),work(:),buffer(:,:)
  real(8)::basis(165),querywork(1),weight
  integer::nn,p,ii,jj,rank,info,lwork,allocation_status
  code=1;nn=size(points,2);p=(L+1)*(L+2)*(L+3)/6
  if(L<0.or.L>8.or.nn<p.or..not.ieee_is_finite(decay))return
  allocate(a(nn,p),rhs(nn,3),sv(p),buffer(3,size(queries,2)),stat=allocation_status)
  if(allocation_status/=0)return
  do ii=1,nn
  do jj=1,3
   if(.not.ieee_is_finite(points(jj,ii)).or..not.ieee_is_finite(values(jj,ii)))return
  enddo
  enddo
  do ii=1,size(queries,2)
  do jj=1,3
   if(.not.ieee_is_finite(queries(jj,ii)))return
  enddo
  enddo
  call matrix(points,values,L,decay,a,rhs)
  call dgelss(nn,p,3,a,nn,rhs,nn,sv,1d-12,rank,querywork,-1,info)
  if(info/=0)return
  lwork=ceiling(querywork(1));allocate(work(lwork),stat=allocation_status)
  if(allocation_status/=0)return
  call matrix(points,values,L,decay,a,rhs)
  call dgelss(nn,p,3,a,nn,rhs,nn,sv,1d-12,rank,work,lwork,info)
  if(info/=0.or.rank/=p)return
  do ii=1,size(queries,2)
   call near_basis(queries(:,ii),L,basis);weight=exp(-decay*queries(3,ii))
   buffer(:,ii)=matmul(transpose(rhs(1:p,:)),basis(1:p))*weight
   do jj=1,3
    if(.not.ieee_is_finite(buffer(jj,ii)))return
   enddo
  enddo
  prediction=buffer;code=0
 end subroutine
 subroutine matrix(points,values,L,decay,a,rhs)
  real(8),intent(in)::points(:,:),values(:,:),decay
  integer,intent(in)::L
  real(8),intent(out)::a(:,:),rhs(:,:)
  real(8)::basis(165)
  integer::kk
  do kk=1,size(points,2)
   call near_basis(points(:,kk),L,basis)
   a(kk,:)=basis(1:size(a,2))*exp(-decay*points(3,kk))
  enddo
  rhs=transpose(values)
 end subroutine
end program
