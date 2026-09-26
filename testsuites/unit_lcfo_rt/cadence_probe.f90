program cadence_probe
 use lcfo_rt_wannier
 use lcfo_rt_basis
 implicit none
 complex(8) :: c(4,2),changed(4,2),stored_c(4,2),u(2,2),expected(32,2)
 complex(8),allocatable :: source(:,:,:),first(:,:,:)
 integer :: header(7),iu,j,step,g
 real(8) :: radius,delta(3),position(3)
 character(32) :: value
 real(8) :: h(3),centers(3,2)
 allocate(lcfo_basis(32,4),lcfo_counts(1),lcfo_offsets(2),lcfo_origins(3,1))
 lcfo_basis=0d0;lcfo_counts=4;lcfo_offsets=[0,4];lcfo_origins=0
 lcfo_basis(1,1)=1d0;lcfo_basis(2,2)=1d0;lcfo_basis(5,3)=1d0;lcfo_basis(6,4)=1d0
 c(:,1)=.5d0;c(:,2)=[.5d0,.5d0,-.5d0,-.5d0]
 call lcfo_mlwf_configure()
 call lcfo_mlwf_source(c,lcfo_basis,[1,2,3,4],first)
 open(newunit=iu,file='lcfo_mlwf_initial.bin',form='unformatted',access='stream')
 read(iu)header,h,stored_c,u,centers
 close(iu)
 call get_environment_variable('SALMON_LCFO_RT_RADIUS',value)
 read(value,*)radius
 ! Step1 always transports, both predictor and accepted endpoint.
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(c,lcfo_basis,[1,2,3,4],source)
 call lcfo_mlwf_stage(1)
 call lcfo_mlwf_stage(2)
 call lcfo_mlwf_accept_cached()
 ! Step2 holds U even when the occupied orbitals acquire relative phases.
 do j=1,2;changed(:,j)=c(:,j)*exp(cmplx(0d0,.7d0*j*j,8));enddo
 expected=matmul(lcfo_basis,matmul(changed,u))
 if(radius>0d0)then
  do j=1,2;do g=1,32
   position=real([mod(g-1,8),mod((g-1)/8,2),(g-1)/16],8)
   delta=modulo(position-centers(:,j)+.5d0*lcfo_grid,real(lcfo_grid,8))-.5d0*lcfo_grid
   if(sum(delta**2)>radius**2)expected(g,j)=0d0
  enddo;enddo
 endif
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(changed,lcfo_basis,[1,2,3,4],source)
 if(maxval(abs(source(:,:,1)-expected))>1d-12)error stop 'U not held on step2'
 call lcfo_mlwf_stage(1)
 call lcfo_mlwf_stage(2)
 call lcfo_mlwf_accept_cached()
 ! Step3 must reset the held-step dephasing against the last transported frame.
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(changed,lcfo_basis,[1,2,3,4],source)
 if(maxval(abs(source(:,:,1)-first(:,:,1)))>1d-12)error stop 'held-step dephasing accumulated'
 call lcfo_mlwf_stage(2)
 call lcfo_mlwf_accept_cached()
 ! Cached acceptance must restore the step3 U after predictor rollback.
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(changed,lcfo_basis,[1,2,3,4],source)
 if(maxval(abs(source(:,:,1)-first(:,:,1)))>1d-12)error stop 'cached U lost on step4'
 call lcfo_mlwf_stage(2)
 call lcfo_mlwf_accept_cached()
 do step=5,101
  do j=1,2;changed(:,j)=c(:,j)*exp(cmplx(0d0,.23d0*step*j*j,8));enddo
  call lcfo_mlwf_stage(0)
  call lcfo_mlwf_source(changed,lcfo_basis,[1,2,3,4],source)
  if(mod(step,2)==1)then
   if(maxval(abs(source-first))>1d-11)error stop 'stationary-density secular spreading'
  endif
  call lcfo_mlwf_stage(1)
  call lcfo_mlwf_stage(2)
  call lcfo_mlwf_accept_cached()
 enddo
 write(*,*) 'Physical-step U cadence and rollback passed'
end program
