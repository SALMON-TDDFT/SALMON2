program probe
 use mpi
 use hse_spatial
 use hse_ace
 use exx_orbitals, only: orbital_hermitian_action
 implicit none
 type(spatial_exx_state) :: op
 type(hse_ace_state) :: ace
 integer :: ierr,np,rank,n(3)=[32,16,16],m(3),g,x,y,z,st,k,mode,compact
 real(8) :: h(3)=[.5d0,.5d0,.5d0],omega,err,global_err,budget,bound,correction,error_scale
 complex(8),allocatable :: target(:,:,:),reference(:,:,:),value(:,:,:),compressed(:,:,:)
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 m=[n(1),n(2)/np,n(3)]
 allocate(target(product(m),2,1),reference(product(m),2,1),value(product(m),2,1),op%source(product(m),2))
 target=0;g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
  g=g+1
  if(z/=3.or.y+rank*m(2)/=3)cycle
  if(x==2)target(g,:,1)=[(1d0,0d0),(0d0,1d-6)]/sqrt(product(h))
  if(x==20)target(g,:,1)=[(0d0,1d-6),(1d0,0d0)]/sqrt(product(h))
 enddo;enddo;enddo
 allocate(compressed(product(m),2,1))
 op%source=target(:,:,1)
 do k=0,2
  omega=.11d0*k
  op%screen_mode=0;op%compact=.false.
  call apply(reference)
  do compact=0,1
   op%compact=compact==1
   do mode=1,2
    op%screen_mode=mode;op%screen_tolerance=1d-2
    call apply(value)
    if(op%screen_candidates/=2)error stop 'expected two small off-diagonal pairs'
    bound=op%screen_bound
    if(bound<=0.or.bound>op%screen_tolerance)error stop 'invalid bound'
    err=sum(abs(value-reference)**2)*product(h)
    call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(mode==1.and.sqrt(global_err)>1d-11)error stop 'diagnostics changed action'
    if(mode==2.and.(op%screen_skipped/=2.or.sqrt(global_err)>bound+1d-11))error stop 'omission exceeds bound'
    if(mode==2)then
      compressed=value
      call orbital_hermitian_action(target,value,product(h),MPI_COMM_WORLD,MPI_COMM_SELF,0d0,correction,st)
      if(st==0.or.any(value/=compressed))error stop 'over-budget correction accepted or modified action'
      call orbital_hermitian_action(target,value,product(h),MPI_COMM_WORLD,MPI_COMM_SELF, &
        op%screen_tolerance-bound,correction,st)
      if(st/=0)error stop 'Hermitian completion failed'
      bound=bound+correction
      err=sum(abs(value-reference)**2)*product(h)
      call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
      if(sqrt(global_err)>bound+1d-11)error stop 'completed action exceeds bound'
    endif
    call hse_ace_build(ace,target,value,product(h),st,sum_grid)
    if(st/=0)error stop 'screened ACE metric invalid'
    call hse_ace_apply(ace,target,compressed,st,sum_grid)
    if(st/=0.or.maxval(abs(compressed-value))>1d-10)error stop 'ACE interpolation failed'
    if(rank==0)print *,'PAIR mode/compact/omega/error/bound: ',mode,compact,omega,sqrt(global_err),bound
   enddo
   op%screen_tolerance=0;call apply(value)
   if(op%screen_skipped/=0)error stop 'zero tolerance omitted finite pair'
   if(maxval(abs(value-reference))>1d-11)error stop 'zero tolerance changed action'
  enddo
 enddo
 ! Squared densities underflow here, but the finite L1 density does not.
 target=target*1d-200
 op%screen_mode=2;op%screen_tolerance=1d-250
 call apply(value)
 if(op%screen_skipped/=0)error stop 'underflow falsely certified finite pair'
 op%screen_mode=0;call apply(reference)
 op%screen_mode=2;op%screen_tolerance=1d-190;call apply(value)
 if(op%screen_skipped/=4.or.op%screen_bound<=0d0)error stop 'finite omitted bound underflowed'
 err=maxval(abs(value-reference))
 call MPI_Allreduce(err,error_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierr)
 err=sum((abs(value-reference)/error_scale)**2)*product(h)
 call MPI_Allreduce(err,global_err,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
 if(error_scale*sqrt(global_err)>op%screen_bound)error stop 'tiny omitted bound invalid'
 ! Completion must also measure tiny finite changes without squaring to zero.
 value=reference;value(:,1,1)=value(:,1,1)*(1d0,1d0)
 call orbital_hermitian_action(target*1d200,value,product(h),MPI_COMM_WORLD,MPI_COMM_SELF, &
   1d-190,correction,st)
 if(st/=0.or.correction<=0d0)error stop 'tiny completion norm underflowed'
 call localized_chain()
 call MPI_Finalize(ierr)
contains
 subroutine localized_chain()
  type(spatial_exx_state) :: chain
  complex(8),allocatable :: t(:,:,:),exact(:,:,:),fast(:,:,:)
  integer :: nn(3),mm(3),gg,xx,yy,zz,jj,ii,nnorb,ss,cc,kk
  real(8) :: ee,ge,ww
  nnorb=8;nn=[64,16,16];mm=[64,16/np,16]
  allocate(t(product(mm),nnorb,1),exact(product(mm),nnorb,1),fast(product(mm),nnorb,1))
  allocate(chain%source(product(mm),nnorb));chain%source=0d0
  gg=0
  do zz=0,mm(3)-1;do yy=0,mm(2)-1;do xx=0,mm(1)-1
   gg=gg+1
   if(yy+rank*mm(2)/=2.or.zz/=2)cycle
   do jj=1,nnorb
    if(xx==8*(jj-1)+1)chain%source(gg,jj)=cmplx(1d0,0d0,8)/sqrt(product(h))
   enddo
  enddo;enddo;enddo
  t(:,:,1)=chain%source+cmplx(1d-12,-1d-12,8)
  do kk=0,1
   ww=.11d0*kk;chain%screen_mode=0;chain%compact=.false.
   call spatial_exx_apply(chain,nn,h,[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
       MPI_COMM_WORLD,4d0,t,exact,ss,omega=ww)
   if(ss/=0)error stop 'chain reference'
   do cc=0,1
    chain%screen_mode=2;chain%screen_tolerance=1d-7;chain%compact=cc==1
    call spatial_exx_apply(chain,nn,h,[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
        MPI_COMM_WORLD,4d0,t,fast,ss,omega=ww)
    if(ss/=0)error stop 'chain apply'
    if(chain%pair_products/=nnorb)error stop 'distant pairs generated full products'
    if(chain%screen_skipped/=nnorb*(nnorb-1))error stop 'chain pair reduction'
    ee=sum(abs(fast-exact)**2)*product(h)
    call MPI_Allreduce(ee,ge,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
    if(sqrt(ge)>chain%screen_bound+1d-13.or.chain%screen_bound>1d-7)error stop 'chain bound'
    if(rank==0)print *,'CHAIN ACTION omega/products/skipped/error/bound',ww,chain%pair_products, &
        chain%screen_skipped,sqrt(ge),chain%screen_bound
   enddo
  enddo
 end subroutine
 subroutine sum_grid(a)
  complex(8),intent(inout) :: a(:,:)
  complex(8) :: b(size(a,1),size(a,2))
  call MPI_Allreduce(a,b,size(a),MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  a=b
 end subroutine
 subroutine apply(a)
  complex(8),intent(out) :: a(:,:,:)
  call spatial_exx_apply(op,n,h,[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF],MPI_COMM_WORLD, &
    4d0,target,a,st,omega=omega)
  if(st/=0)error stop 'apply failed'
 end subroutine
end program
