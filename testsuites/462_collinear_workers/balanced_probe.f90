program balanced_probe
 use iso_fortran_env,only:int64
 use exx_worker_layout,only:choose_balanced_strips
 implicit none
 integer::n,ng,nk,nc,np,err,groups,chunk,i
 integer,allocatable::r(:),c(:),count(:),group(:),leader(:)
 integer(int64)::peak,active
 read(*,*)n,ng,nk,nc,np
 allocate(r(np),c(np),count(np),group(np),leader(np))
 call choose_balanced_strips(n,ng,nk,nc,np,r,c,count,group,leader,groups,chunk,peak,active,err)
 print *,err,groups,chunk,peak,active
 if(err/=0)stop
 do i=1,np
 print *,i-1,r(i),c(i),count(i),group(i),leader(group(i)+1)
 enddo
end program
