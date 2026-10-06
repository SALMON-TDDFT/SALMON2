program active_wf_probe
 use lcfo_wf_support
 implicit none
 type(s_lcfo_wf_plan) :: plan
 type(s_lcfo_wf_kernel) :: kernel
 complex(8) :: basis(12,3),frame(3,4),dense(12,4),gram(3,3)
 complex(8),allocatable :: compact(:,:),selected(:,:)
 real(8) :: positions(3,12),centers(3,4),length(3),radius,d(3),total,inside,reference,loss
 logical :: protected(4),keep(12,4)
 integer :: i,j,k,trial,g
 complex(8),allocatable :: local_result(:,:)
 length=[20d0,12d0,8d0];g=0
 do k=0,1;do j=0,1;do i=0,2
  g=g+1;positions(:,g)=[dble(i),dble(j),dble(k)]
 enddo;enddo;enddo
 centers(:,1)=[19.8d0,.2d0,.3d0];centers(:,2)=[10d0,6d0,4d0]
 centers(:,3)=[1.2d0,.1d0,.7d0];centers(:,4)=[11d0,6d0,4d0]
 do j=1,3;do i=1,12;basis(i,j)=cmplx(sin(.19d0*i*j),cos(.13d0*(i+j)),8);enddo;enddo
 do j=1,4;do i=1,3;frame(i,j)=cmplx(cos(.17d0*i*j),sin(.23d0*(i-j)),8);enddo;enddo
 basis(1:6,3)=0d0;basis(7:12,1)=0d0
 basis(12,:)=0d0
 basis(1,3)=cmplx(1d-25,-2d-25,8)
 dense=matmul(basis,frame);gram=matmul(conjg(transpose(basis)),basis)
 do trial=1,4
  radius=1.1d0;protected=.false.
  if(trial==1)protected(4)=.true.
  if(trial==2)radius=.001d0
  if(trial==3)radius=0d0
  if(trial==4)protected=.true.
  do j=1,4;do i=1,12
   d=modulo(positions(:,i)-centers(:,j)+.5d0*length,length)-.5d0*length
   keep(i,j)=radius==0d0.or.protected(j).or.sum(d*d)<=radius*radius
  enddo;enddo
  call lcfo_wf_plan_init(plan,positions,centers,length,radius,protected)
  call lcfo_wf_reconstruct(plan,basis,frame,compact)
  call lcfo_wf_kernel_init(kernel,plan,basis)
  call lcfo_wf_kernel_apply(kernel,frame(:,plan%columns),local_result)
  if(any(shape(local_result)/=shape(compact)))error stop 'local kernel shape'
  if(any(abs(local_result-compact)>1d-12))error stop 'local kernel differs'
  if(trial==1.and.kernel%products>=12*3*size(plan%columns))error stop 'no reconstruction work saved'
  if(size(compact,2)/=count(any(keep,dim=1)))error stop 'inactive columns reconstructed'
  if(trial==1.and.size(compact,2)/=3)error stop 'periodic/protected selection'
  if(trial==2.and.size(compact,2)/=0)error stop 'empty support'
  do j=1,size(plan%columns)
   k=plan%columns(j)
   do i=1,12
    if(abs(compact(i,j)-merge(dense(i,k),(0d0,0d0),keep(i,k)))>1d-12)error stop 'masked WF differs'
   enddo
  enddo
  call lcfo_wf_reconstruct(plan,basis,frame(:,plan%columns),selected,compact=.true.)
  if(any(shape(selected)/=shape(compact)))error stop 'compact input shape'
  if(any(abs(selected-compact)>1d-12))error stop 'compact input reconstruction'
  total=lcfo_wf_total_norm(gram,frame);reference=sum(abs(dense)**2)
  if(abs(total-reference)>1d-12*reference)error stop 'nonorthogonal total norm'
  inside=sum(abs(compact)**2);loss=sum(abs(dense)**2,mask=.not.keep)
  if(abs((total-inside)-loss)>1d-12*reference)error stop 'discarded norm changed'
 enddo
 ! Isolated tiny nonzero term: an absolute tolerance on O(1) data would miss its removal.
 call lcfo_wf_plan_init(plan,positions(:,1:1),centers(:,1:1),length,0d0,protected(1:1))
 call lcfo_wf_kernel_init(kernel,plan,reshape([cmplx(1d-25,-2d-25,8)],[1,1]))
 call lcfo_wf_kernel_apply(kernel,reshape([cmplx(3d0,4d0,8)],[1,1]),local_result)
 if(abs(local_result(1,1)/cmplx(11d-25,-2d-25,8)-1d0)>1d-14)error stop 'tiny term dropped'
 if(kernel%products/=1)error stop 'tiny term work count'
 write(*,*)'Active WF reconstruction and Gram diagnostic passed'
end program
