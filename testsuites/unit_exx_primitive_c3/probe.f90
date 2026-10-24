program primitive_c3
 use exx_symmetry
 implicit none
 type(s_exx_symmetry_map)::map,old
 real(8)::cell(3,3),inv(3,3),ops(3,4,3),id(3,3),p(3,3),k(3,4),f(3,4),wt(4),pi2
 real(8)::badcell(3,3),one(3,4,1),zero_k(3,1),one_w(1)
 real(8)::origin(3),atom(3,1),m(3),mr(3),expected(3),pos(3),rk(3),angle
 complex(8)::u(2,2),spinor(2),z(128,1),out(128,1),phase
 integer::i,j,s,g,r,err,c(3)
 id=0d0
 do i=1,3
 id(i,i)=1d0
 enddo
 cell=1d0-id;inv=0.5d0- id;pi2=2d0*acos(-1d0)
 p=0d0;p(2,1)=1d0;p(3,2)=1d0;p(1,3)=1d0
 origin=[0.25d0,0d0,0d0];ops=0d0
 ops(:,1:3,1)=id;ops(:,1:3,2)=p;ops(:,1:3,3)=matmul(p,p)
 do s=1,3
 ops(:,4,s)=origin-matmul(ops(:,1:3,s),origin)
 enddo
 f(:,1)=0d0;f(:,2)=-0.5d0;f(:,3)=[-.5d0,0d0,0d0];f(:,4)=[-.5d0,-.5d0,0d0]
 wt=[1d0,1d0,3d0,3d0]/8d0;k=pi2*matmul(transpose(inv),f)
 call symmetry_init(map,[4,4,4],[0.25d0,0.25d0,0.25d0],k,wt,ops,ops,2,err,cell=cell)
 if(err/=0) stop 'primitive C3 initialization rejected'
 if(maxval(abs(map%cell-cell))>1d-13.or.maxval(abs(map%inverse_cell-inv))>1d-13) stop 'cell inverse mismatch'
 atom(:,1)=matmul(cell,origin)
 call symmetry_validate_atoms(map,atom,[1],err)
 if(err/=0) stop 'primitive C3 atom rejected'
 call symmetry_spin_rotation(p,u,err)
 if(err/=0.or.maxval(abs(matmul(conjg(transpose(u)),u)-id(:2,:2)))>1d-13) stop 'C3 SU2 invalid'
 spinor=[cmplx(.7d0,.2d0,8),cmplx(-.1d0,.4d0,8)]
 m=mag(spinor);mr=mag(matmul(u,spinor))
 if(maxval(abs(mr-matmul(p,m)))>1d-13) stop 'C3 spin covariance failed'
 do g=1,64
 z(g,1)=spinor(1)*cmplx(cos(.13d0*g),sin(.13d0*g),8)
 z(64+g,1)=spinor(2)*cmplx(cos(.13d0*g),sin(.13d0*g),8)
 enddo
 do j=1,map%nfull
 r=map%owner(j)
 do i=1,map%multiplicity(j)
 s=map%operations(i,j)
 call symmetry_transform(map,j,i,z,out,err)
 if(err/=0) stop 'primitive C3 transform failed'
 call symmetry_spin_rotation(map%cart_rotation(:,:,s),u,err)
 rk=matmul(map%cart_rotation(:,:,s),map%rep_k(:,r))
 do g=1,64
 c=[mod(g-1,4),mod((g-1)/4,4),(g-1)/16]
 pos=matmul(cell,real(c,8)/4d0)
 angle=dot_product(rk-map%full_k(:,j),pos)-dot_product(rk,matmul(cell,map%translation(:,s)))
 phase=cmplx(cos(angle),sin(angle),8)
 spinor=phase*matmul(u,[z(map%grid(g,s),1),z(64+map%grid(g,s),1)])
 if(maxval(abs([out(g,1),out(64+g,1)]-spinor))>1d-13) stop 'primitive phase/rotation mismatch'
 enddo
 enddo
 enddo
 ! Existing signed diagonal operation retains its exact established SU2 sign.
 p=id;p(2,2)=-1d0;p(3,3)=-1d0
 call symmetry_spin_rotation(p,u,err)
 if(err/=0.or.abs(u(1,2)-cmplx(0d0,-1d0,8))>1d-14.or.abs(u(2,1)-u(1,2))>1d-14) stop 'old signed rotation changed'
 ! The old no-cell interface still constructs the diagonal Cartesian cell.
 one=0d0;one(:,1:3,1)=id;zero_k=0d0;one_w=1d0
 call symmetry_init(old,[2,3,4],[.7d0,.8d0,.9d0],zero_k,one_w,one,one,1,err)
 if(err/=0) stop 'legacy initialization changed'
 do i=1,3
 if(abs(old%cell(i,i)-old%n(i)*old%h(i))>1d-14) stop 'legacy diagonal cell changed'
 enddo
 ! Check inverse orientation for a nonsymmetric cell with identity symmetry.
 badcell=cell;badcell(1,1)=0.13d0;badcell(2,3)=1.17d0
 call symmetry_init(old,[2,2,2],[1d0,1d0,1d0],zero_k,one_w,one,one,1,err,cell=badcell)
 if(err/=0.or.maxval(abs(matmul(old%cell,old%inverse_cell)-id))>1d-13) stop 'general cell inverse mismatch'
 call symmetry_init(old,[4,4,4],[.25d0,.25d0,.25d0],k,wt,ops,ops,2,err,cell=badcell)
 if(err==0) stop 'nonisometric C3 accepted'
 print *,'primitive translated FCC C3 and legacy diagonal PASS'
contains
 function mag(v) result(m)
 complex(8),intent(in)::v(2)
 real(8)::m(3)
 m=[2d0*real(conjg(v(1))*v(2),8),2d0*aimag(conjg(v(1))*v(2)),abs(v(1))**2-abs(v(2))**2]
 end function
end program
