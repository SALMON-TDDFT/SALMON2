program stream_probe
 use iso_fortran_env, only: int64
 use lcfo_mesh_stream
 implicit none
 integer,parameter :: core(3)=[80,72,1],nb=35,no=3
 integer :: ub,uc,i,j,g,b,io,status,lo(3),m(3),first,count,ix,iy,idx
 integer :: mapping(80,3)
 integer(int64) :: bpos,cpos
 complex(8) :: basis(5760,nb),coef(nb,no)
 complex(8),allocatable :: tile(:),expected(:)
 real(8) :: wire(2)
 character(13) :: prefix='offset-marker'
 do j=1,3;do b=1,nb
 coef(b,j)=cmplx(.003d0*(b+j),.002d0*(b-j),8)
 enddo;enddo
 do b=1,nb;do g=1,5760
 basis(g,b)=cmplx(sin(.013d0*g*b),cos(.017d0*(g+b)),8)
 enddo;enddo
 open(newunit=ub,status='scratch',access='stream',form='unformatted')
 open(newunit=uc,status='scratch',access='stream',form='unformatted')
 write(ub)prefix;inquire(unit=ub,pos=bpos)
 write(uc)prefix;inquire(unit=uc,pos=cpos)
 do b=1,nb;do g=1,5760
 write(ub)real(basis(g,b),8),aimag(basis(g,b))
 enddo;enddo
 do j=1,no;do b=1,nb
 write(uc)real(coef(b,j),8),aimag(coef(b,j))
 enddo;enddo
 mapping=1
 do i=1,80
 mapping(i,1)=modulo(i+2,80)+1
 enddo
 do i=1,72
 mapping(i,2)=modulo(i+4,72)+1
 enddo
 do j=1,2
 lo=[1,1,1];m=core
 if(j==2)then
 lo=[3,2,1];m=[70,69,1]
 endif
 first=1;count=product(m)
 allocate(tile(count),expected(count));tile=0d0;expected=0d0
 io=2
 call lcfo_stream_contract(ub,uc,bpos,cpos,core,nb,io,mapping,lo,m,first,tile,status)
 if(status/=0)error stop 'stream read failed'
 do iy=1,72;do ix=1,80
 if(mapping(ix,1)<lo(1).or.mapping(ix,1)>=lo(1)+m(1))cycle
 if(mapping(iy,2)<lo(2).or.mapping(iy,2)>=lo(2)+m(2))cycle
 idx=1+mapping(ix,1)-lo(1)+m(1)*(mapping(iy,2)-lo(2))
 g=ix+80*(iy-1)
 do b=1,nb
 expected(idx)=expected(idx)+basis(g,b)*coef(b,io)
 enddo
 enddo;enddo
 if(maxval(abs(tile-expected))>1d-12)error stop 'stream contraction mismatch'
 ! Read a chunk that starts inside a row; zero bands preserve previous contributions.
 first=17;count=product(m)-31
 tile=0d0
 call lcfo_stream_contract(ub,uc,bpos,cpos,core,nb,io,mapping,lo,m,first,tile(:count),status)
 if(status/=0.or.maxval(abs(tile(:count)-expected(first:first+count-1)))>1d-12)error stop 'offset chunk'
 call lcfo_stream_contract(ub,uc,bpos,cpos,core,0,io,mapping,lo,m,first,tile(:count),status)
 if(status/=0.or.maxval(abs(tile(:count)-expected(first:first+count-1)))>1d-12)error stop 'zero bands'
 deallocate(tile,expected)
 enddo
 ! A contiguous requested domain larger than the removed65536-point cap.
 block
   integer,allocatable :: large_map(:,:)
   complex(8),allocatable :: large(:),actual(:)
   integer :: q
   allocate(large_map(70000,3),large(70000),actual(70000))
   large_map=1
   do q=1,70000
     large_map(q,1)=q
     large(q)=cmplx(sin(.01d0*q),cos(.02d0*q),8)
   enddo
   write(ub,pos=bpos)large
   write(uc,pos=cpos)1d0,0d0
   actual=0d0
   call lcfo_stream_contract(ub,uc,bpos,cpos,[70000,1,1],1,1,large_map, &
     [1,1,1],[70000,1,1],1,actual,status)
   if(status/=0.or.maxval(abs(actual-large))>1d-14)error stop 'domain-dependent read size'
 end block
 close(ub);close(uc)
 print *, 'PASS streamed complex fragment payload'
end program
