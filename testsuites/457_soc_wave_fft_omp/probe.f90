program probe
 use iso_c_binding
 use exx_spinor_wave_fft
 !$ use omp_lib,only:omp_get_max_threads
 implicit none
 integer,parameter::b=3,nc=5,nk=6
 complex(c_double_complex),target::rho(b,nc,nk,4)
 complex(c_double_complex)::initial(b,nc,nk,4),expected(b,nc,nk,4),v(nk)
 real(8)::kernel(0:3,0:5,0:0),pi,angle,error,t0,t1
 integer::point(3,20),shift(3,nk),first(3),k,j,c,r,g,a,d,i,pass,err,threads,want
 integer::mesh(3),n(3),index(3),kk(3),jj(3),rows,cols,lo,col
 type(wave_fft_plans)::plans
 pi=acos(-1d0);mesh=[2,3,1];n=[2,2,1];first=[1,2,3];lo=1;col=3
 do g=1,20
 point(:,g)=[mod(g-1,2),mod((g-1)/2,2),0]
 enddo
 do k=1,nk
 shift(:,k)=[mod(k-1,2),(k-1)/2,0]
 enddo
 do j=0,5;do i=0,3
 kernel(i,j,0)=.2d0+sin(real(i+2*j,8))
 enddo;enddo
 do c=1,4;do k=1,nk;do g=1,nc;do r=1,b
 initial(r,g,k,c)=cmplx(sin(real(r+3*g+2*k+c,8)),cos(real(2*r-g+k+3*c,8)),8)
 enddo;enddo;enddo;enddo
 rho=initial
 call wave_fft_create(plans,mesh,rho,err)
 if(err/=0)stop 1
 want=1
 !$ want=omp_get_max_threads()
 error=0
 do pass=1,5
 rows=b;cols=nc
 if(pass==2)then
 rows=2;cols=4
 else if(pass==3)then
 rows=0
 else if(pass==4)then
 cols=0
 else if(pass==5)then
 rows=0;cols=0
 endif
 expected=initial
 do c=1,4;do g=1,nc;do r=1,b
 v=0
 do k=1,nk
 kk=shift(:,k)
 do j=1,nk
 jj=shift(:,j);angle=-2*pi*sum(real(kk*jj,8)/mesh)
 v(k)=v(k)+initial(r,g,j,c)*exp(cmplx(0d0,angle,8))
 enddo
 if(r<=rows.and.g<=cols)then
 index=modulo(point(:,first(lo+r-1))-point(:,col+g-1)-shift(:,k),n*mesh)
 v(k)=v(k)*kernel(index(1),index(2),index(3))
 endif
 enddo
 do j=1,nk
 expected(r,g,j,c)=0
 do k=1,nk
 angle=2*pi*sum(real(shift(:,k)*shift(:,j),8)/mesh)
 expected(r,g,j,c)=expected(r,g,j,c)+v(k)*exp(cmplx(0d0,angle,8))
 enddo
 enddo
 enddo;enddo;enddo
 rho=initial
 call wave_fft_execute(plans,rho,n,point,shift,kernel,first,lo,col,rows,cols,threads)
 if(threads/=want)stop 2
 error=max(error,maxval(abs(rho-expected)))
 enddo
 if(error>1d-11)stop 3
 call cpu_time(t0)
 do pass=1,40
 rho=initial
 call wave_fft_execute(plans,rho,n,point,shift,kernel,first,lo,col,b,nc,threads)
 enddo
 call cpu_time(t1)
 call wave_fft_destroy(plans)
 call wave_fft_destroy(plans)
 print *, 'WAVE_FFT_OMP_PASS',threads,error,' CPU_SECONDS ',t1-t0,' RHO_BYTES ',size(rho)*16
end program
