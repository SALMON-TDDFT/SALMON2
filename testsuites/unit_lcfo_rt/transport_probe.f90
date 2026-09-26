program transport_probe
 use lcfo_rt_wannier
 use lcfo_rt_basis
 use hse_wannier_gauge, only: gauge_transport
 implicit none
 complex(8) :: c(4,2),changed(4,2)
 complex(8),allocatable :: first(:,:,:),next(:,:,:)
 integer :: j,g,iu,header(7),status
 real(8) :: h(3),centers(3,2),minimum_overlap,distance,delta(3),position(3)
 complex(8) :: stored_c(4,2),initial_u(2,2),reference(4,2,1),up(2,2,1),predicted(4,2),corrected(4,2)
 complex(8) :: v1(4),v2(4),expected(32,2,1)
 allocate(lcfo_basis(32,4),lcfo_counts(1),lcfo_offsets(2),lcfo_origins(3,1))
 lcfo_basis=0d0;lcfo_counts=4;lcfo_offsets=[0,4];lcfo_origins=0
 lcfo_basis(1,1)=1d0;lcfo_basis(2,2)=1d0;lcfo_basis(5,3)=1d0;lcfo_basis(6,4)=1d0
 c(:,1)=.5d0;c(:,2)=[.5d0,.5d0,-.5d0,-.5d0]
 call lcfo_mlwf_configure()
 call lcfo_mlwf_source(c,lcfo_basis,[1,2,3,4],first)
 do j=1,2
   changed(:,j)=c(:,j)*exp(cmplx(0d0,.7d0*j*j,8))
 enddo
 call lcfo_mlwf_source(changed,lcfo_basis,[1,2,3,4],next)
 write(*,*) 'stationary-density orbital-phase source error=',maxval(abs(next-first))
 if(maxval(abs(next-first))>1d-12)error stop 'Spurious WF spreading under occupied orbital phases'
 ! Regression: predicted frame -> rollback -> exact corrected cache hit -> next step.
 open(newunit=iu,file='lcfo_mlwf_initial.bin',form='unformatted',access='stream')
 read(iu)header,h,stored_c,initial_u,centers
 close(iu)
 if(header(2)/=2)error stop 'expected xyz center format version 2'
 reference(:,:,1)=matmul(c,initial_u)
 v1=[.5d0,-.5d0,.5d0,-.5d0];v2=[.5d0,-.5d0,-.5d0,.5d0]
 predicted(:,1)=cos(.2d0)*c(:,1)+sin(.2d0)*v1;predicted(:,2)=c(:,2)
 call gauge_transport(reshape(predicted,[4,2,1]),reference,1d0,up,minimum_overlap,status)
 if(status/=0)error stop 'expected predicted transport'
 reference(:,:,1)=matmul(predicted,up(:,:,1))
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(predicted,lcfo_basis,[1,2,3,4],next)
 call lcfo_mlwf_stage(1)
 call lcfo_mlwf_stage(2)
 call lcfo_mlwf_accept_cached()
 corrected(:,1)=cos(.2d0)*c(:,1)+sin(.2d0)*(cos(.7d0)*v1+sin(.7d0)*v2)
 corrected(:,2)=cos(.15d0)*c(:,2)+sin(.15d0)*(-sin(.7d0)*v1+cos(.7d0)*v2)
 call gauge_transport(reshape(corrected,[4,2,1]),reference,1d0,up,minimum_overlap,status)
 if(status/=0)error stop 'expected corrected transport'
 expected(:,:,1)=matmul(lcfo_basis,matmul(corrected,up(:,:,1)))
 do j=1,2;do g=1,32
   position=real([mod(g-1,8),mod((g-1)/8,2),(g-1)/16],8)
   delta=modulo(position-centers(:,j)+.5d0*lcfo_grid,real(lcfo_grid,8))-.5d0*lcfo_grid
   distance=sqrt(sum(delta**2))
   if(distance>1d0)expected(g,j,1)=0d0
 enddo;enddo
 call lcfo_mlwf_stage(0)
 call lcfo_mlwf_source(corrected,lcfo_basis,[1,2,3,4],next)
 if(maxval(abs(next-expected))>1d-12)error stop 'Predictor rollback/cache reference contamination'
 call lcfo_mlwf_stage(2)
 write(*,*) 'Predictor rollback and cached acceptance passed'
 ! No real-space source is required on a retained-ACE step. Transport through
 ! a nontrivial occupied subspace, then compare the next actual source.
 call lcfo_mlwf_track(corrected)
 reference(:,:,1)=matmul(corrected,up(:,:,1))
 call gauge_transport(reshape(predicted,[4,2,1]),reference,1d0,up,minimum_overlap,status)
 expected(:,:,1)=matmul(lcfo_basis,matmul(predicted,up(:,:,1)))
 do j=1,2;do g=1,32
   position=real([mod(g-1,8),mod((g-1)/8,2),(g-1)/16],8)
   delta=modulo(position-centers(:,j)+.5d0*lcfo_grid,real(lcfo_grid,8))-.5d0*lcfo_grid
   distance=sqrt(sum(delta**2))
   if(distance>1d0)expected(g,j,1)=0d0
 enddo;enddo
 call lcfo_mlwf_source(predicted,lcfo_basis,[1,2,3,4],next)
 if(maxval(abs(next-expected))>1d-12)error stop 'Transport-only frame lost'
 write(*,*) 'Transport-only retained-ACE frame passed'
end program
