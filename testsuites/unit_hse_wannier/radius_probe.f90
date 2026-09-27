program probe
 use hse_wannier
 implicit none
 type(s_hse_wannier) :: op
 integer,parameter :: n(3)=[8,6,4],ng=192
 integer :: nk,ik,g,j,a,status,kind,p(3)
 real(8) :: h(3)=[.7d0,.8d0,.9d0],length(3),center(3),delta(3),radius,omega,before,expected_loss
 real(8),allocatable :: k(:,:)
 complex(8),allocatable :: saved(:,:),target(:,:,:),action(:,:,:),dense(:,:,:),full(:,:,:)
 complex(8) :: metric(3,3)
 do nk=1,2
 allocate(k(3,nk));k=0
 if(nk==2)k(1,2)=acos(-1d0)/(n(1)*h(1))
 allocate(target(ng,3,nk),action(ng,3,nk),dense(ng,3,nk),full(ng,3,nk))
 do kind=1,2
 omega=0d0;if(kind==1)omega=.11d0
 call wannier_init(op,n,[nk,1,1],h,k,omega,status)
 if(status/=0)error stop 'init'
 length=op%ns*h
 allocate(op%source(op%ngs,3),saved(op%ngs,3))
 do g=1,op%ngs
  do j=1,2
   center=0;if(j==2)center(1)=length(1)-h(1)
   delta=modulo(op%point(:,g)*h-center+length/2,length)-length/2
   op%source(g,j)=exp(-sum(delta**2)/.6d0)
  enddo
  op%source(g,3)=1d0/sqrt(real(op%ngs,8)) ! undefined center, stays uncut
 enddo
 saved=op%source
 do ik=1,nk;do j=1,3;do g=1,ng
 target(g,j,ik)=cmplx(sin(real(g*j+ik,8)),cos(real(g+7*j*ik,8)),8)
 enddo;enddo;enddo
 call wannier_apply(op,target,full,status)
 if(status/=0)error stop 'full action'
 call wannier_truncate_source(op,0d0,status)
 if(status/=0.or.any(op%source/=saved))error stop 'zero radius'
 call wannier_truncate_source(op,sum(length),status)
 if(status/=0.or.any(op%source/=saved))error stop 'large radius'
 radius=1.1d0
 call wannier_truncate_source(op,radius,status)
 if(status/=0)error stop 'mask'
 if(op%protected_sources/=1)error stop 'protected count'
 if(any(op%source(:,3)/=saved(:,3)))error stop 'protected source'
 do g=1,op%ngs;do j=1,2
 center=0;if(j==2)center(1)=length(1)-h(1)
 delta=modulo(op%point(:,g)*h-center+length/2,length)-length/2
 if(sum(delta**2)<=radius**2)then
  if(op%source(g,j)/=saved(g,j))error stop 'inside sphere'
 else
  if(op%source(g,j)/=(0d0,0d0))error stop 'outside sphere'
 endif
 enddo;enddo
 before=sum(abs(saved)**2)
 expected_loss=(before-sum(abs(op%source)**2))/before
 if(abs(op%discarded_norm_fraction-expected_loss)>1d-13.or.expected_loss<=0)error stop 'norm loss'
 op%compact_source_support=.true.
 call wannier_apply(op,target,action,status)
 if(status/=0)error stop 'compact apply'
 op%compact_source_support=.false.
 call wannier_apply(op,target,dense,status)
 if(status/=0.or.maxval(abs(action-dense))>1d-11)error stop 'compact/dense'
 if(maxval(abs(action-full))<1d-6)error stop 'radius has no effect'
 metric=0
 do ik=1,nk
 metric=metric+matmul(conjg(transpose(target(:,:,ik))),action(:,:,ik))
 enddo
 if(maxval(abs(metric-conjg(transpose(metric))))>1d-10)error stop 'Hermiticity'
 ! A lost overlap cannot inherit a past converged localization status.
 allocate(op%gauge(3,3,nk),op%previous(ng,3,nk))
 op%gauge=0d0;op%previous=0d0
 do ik=1,nk;do j=1,3
 op%gauge(j,j,ik)=1d0
 enddo;enddo
 op%last_localization_status=0
 call wannier_localize(op,target,0,1d-7,status)
 if(status/=0.or.op%last_localization_status==0)error stop 'stale localization after overlap loss'
 call wannier_destroy(op);deallocate(saved)
 enddo
 deallocate(k,target,action,dense,full)
 enddo
 print *,'radius checks passed'
end program
