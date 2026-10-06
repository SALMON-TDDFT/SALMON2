program probe
 use lcfo_projection
 implicit none
 type(s_lcfo_projection) :: plan
 complex(8) :: basis(13,7),action(13,7),reference(7,7)
 complex(8),allocatable :: result(:,:)
 integer :: i,j,trial
 do j=1,7;do i=1,13
  basis(i,j)=cmplx(sin(.3d0*i*j),cos(.2d0*(i+j)),8)
  action(i,j)=cmplx(cos(.2d0*i*j),sin(.4d0*(i-j)),8)
 enddo;enddo
 do trial=1,3
  if(trial==2)then
   basis([1,3,6,10],:)=0d0;basis(:,[2,4,5])=0d0
  endif
  if(trial==3)basis=0d0
  reference=.07d0*matmul(conjg(transpose(basis)),action)
  reference=.5d0*(reference+conjg(transpose(reference)))
  call lcfo_projection_init(plan,basis,.07d0)
  call lcfo_projection_apply(plan,action,result)
  if(maxval(abs(result-reference))>1d-12)error stop 'projection mismatch'
  if(maxval(abs(result-conjg(transpose(result))))>1d-14)error stop 'Hermiticity'
  if(trial==2)then
   if(size(plan%rows)/=9.or.size(plan%columns)/=4)error stop 'zero support not removed'
  endif
  if(trial==3.and.size(plan%rows)/=0)error stop 'empty support'
 enddo
 print *, 'Compact complex projection passed'
end program
