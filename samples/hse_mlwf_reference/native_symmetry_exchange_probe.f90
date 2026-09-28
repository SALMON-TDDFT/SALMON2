program symmetry_exchange_probe
 use hse_symmetry
 use hse_exchange
 implicit none
 type(hse_symmetry_map)::map
 type(hse_kernel)::kernel
 real(8)::a(3,4,16),b(3,4,16),k(3,8),weights(8),pi,err
 complex(8)::source(8,2,8),target(8,3,8),reference(8,3,8),callback_action(8,3,8)
 complex(8),allocatable::expanded_source(:,:,:),expanded_target(:,:,:)
 integer::s,d,j,g,o,ierr
 pi=acos(-1d0);a=0
 do s=1,16
 do d=1,3
 a(d,d,s)=1-2*ibits(s-1,d-1,1)
 enddo
 a(1,4,s)=0.5d0*ibits(s-1,3,1)
 enddo
 b=a;weights=1d0/8
 do j=1,8
 do d=1,3
 k(d,j)=pi/2*ibits(j-1,d-1,1)
 enddo
 do o=1,2
 do g=1,8
 source(g,o,j)=cmplx(sin(real(g+o*3+j*7,8)),cos(real(g*3+o*7+j,8)),8)
 enddo
 enddo
 do o=1,3
 do g=1,8
 target(g,o,j)=cmplx(cos(real(g+o*5+j*2,8)),sin(real(g*2+o+j*3,8)),8)
 enddo
 enddo
 enddo
 call symmetry_init(map,[2,2,2],[1d0,1d0,1d0],k,weights,a,b,2,ierr)
 if(ierr/=0)stop 1
 if(map%max_little/=16)stop 2
 call symmetry_expand(map,1,source,target,expanded_source,expanded_target,ierr)
 if(ierr/=0)stop 3
 call hse_kernel_init(kernel,2,2,1d0,map%full_k,.11d0,3,ierr,1,8)
 if(ierr/=0)stop 4
 call hse_kernel_apply_distributed(kernel,expanded_source,expanded_target,reference,[1],[8],0,transpose_tiles,ierr)
 if(ierr/=0)stop 5
 call hse_kernel_apply_distributed(kernel,expanded_target,expanded_target,callback_action,[1],[8],0, &
                                 transpose_tiles,ierr,fill_density)
 if(ierr/=0)stop 6
 err=maxval(abs(reference-callback_action))/maxval(abs(reference))
 if(.not.(err<1d-12))stop 7
 print *, 'symmetry exchange callback passed',err
 call hse_kernel_destroy(kernel)
contains
 subroutine transpose_tiles(send,recv,count)
 complex(8),intent(in)::send(:)
 complex(8),intent(out)::recv(:)
 integer,intent(in)::count
 if(count/=size(send))stop 8
 recv=send
 end subroutine
 subroutine fill_density(j,lo,rows,density)
 integer,intent(in)::j,lo,rows
 complex(8),intent(out)::density(:,:)
 complex(8)::u(8,2),ph
 real(8)::rk(3),Gvec(3),r(3),rin(3),angle
 integer::m,op,g,q,c(3),x,y,rep
 density=0;rep=map%owner(j)
 do m=1,map%multiplicity(j)
 op=map%operations(m,j);rk=matmul(b(:,1:3,op),k(:,rep));Gvec=rk-map%full_k(:,j)
 do g=1,8
 r=real([mod(g-1,2),mod((g-1)/2,2),(g-1)/4],8)
 rin=matmul(transpose(a(:,1:3,op)),r-2*a(:,4,op));c=modulo(nint(rin),2)
 q=1+c(1)+2*c(2)+4*c(3)
 angle=dot_product(Gvec,r)-dot_product(rk,2*a(:,4,op))
 ph=cmplx(cos(angle),sin(angle),8)*kernel%phase(g,j)
 u(g,:)=source(q,:,rep)*ph
 enddo
 do y=1,8
 do x=1,rows
 density(x,y)=density(x,y)+sum(u(lo+x-1,:)*conjg(u(y,:)))/map%multiplicity(j)
 enddo
 enddo
 enddo
 end subroutine
end program
