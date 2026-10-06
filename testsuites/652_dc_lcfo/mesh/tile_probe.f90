program tile_probe
 use lcfo_mesh_tile
 implicit none
 integer,parameter :: n(3)=[7,5,4],core(3)=[3,3,2],nb=2
 integer :: map(4,3),lo(3),m(3),first,last,i,j,k,g,b,x,y,z,d
 integer :: full_count(7,5,4),expected_count(7,5,4)
 integer,allocatable :: count(:)
 complex(8) :: basis(3,3,2,nb),coef(nb),full(7,5,4),expected(7,5,4)
 complex(8),allocatable :: tile(:)
 map=1;map(1:3,1)=[7,1,3];map(1:3,2)=[5,2,3];map(1:2,3)=[4,1]
 coef=[cmplx(.3d0,.7d0,8),cmplx(-.2d0,.4d0,8)]
 do b=1,nb;do k=1,2;do j=1,3;do i=1,3
 basis(i,j,k,b)=cmplx(.13d0*(i+2*j+b),.09d0*(k+b*i),8)
 enddo;enddo;enddo;enddo
 expected=0d0;expected_count=0
 do k=1,2;do j=1,3;do i=1,3
 x=map(i,1);y=map(j,2);z=map(k,3)
 expected_count(x,y,z)=expected_count(x,y,z)+1
 do b=1,nb
 expected(x,y,z)=expected(x,y,z)+basis(i,j,k,b)*coef(b)
 enddo
 enddo;enddo;enddo
 ! Unequal destination slabs; seven-point chunks cross local rows/planes.
 full=0d0;full_count=0
 do d=1,2
 lo=[1,1,1];m=n
 if(d==1)then
 m(2)=2
 else
 lo(2)=3;m(2)=3
 endif
 do first=1,product(m),7
 last=min(product(m),first+6)
 allocate(tile(last-first+1),count(last-first+1));tile=0d0;count=0
 call lcfo_tile_coverage(core,map,lo,m,first,count)
 call lcfo_tile_contract(map,lo,m,first,basis,coef,tile)
 ! Zero-band fragments must preserve all existing contributions.
 call lcfo_tile_contract(map,lo,m,first,basis(:,:,:,:0),coef(:0),tile)
 do g=first,last
 x=modulo(g-1,m(1))+lo(1)
 y=modulo((g-1)/m(1),m(2))+lo(2)
 z=(g-1)/(m(1)*m(2))+lo(3)
 full(x,y,z)=tile(g-first+1);full_count(x,y,z)=count(g-first+1)
 enddo
 deallocate(tile,count)
 enddo
 enddo
 if(maxval(abs(full-expected))>1d-14)error stop 'complex tile reconstruction'
 if(any(full_count/=expected_count))error stop 'tile coverage'
 ! Duplicate fragment coverage is retained for collective preflight rejection.
 allocate(count(product(n)));count=0
 call lcfo_tile_coverage(core,map,[1,1,1],n,1,count)
 call lcfo_tile_coverage(core,map,[1,1,1],n,1,count)
 if(maxval(count)/=2.or.minval(count)/=0)error stop 'missing/duplicate coverage'
 print *, 'PASS bounded tile reconstruction'
end program
