program worker_layout_probe
 use iso_fortran_env,only:int64
 use exx_worker_layout,only:choose_worker_layout
 implicit none
 integer::nfirst,ng,nk,nc,np,groups,divisor,rows,cols,b,chunk,ierr,ios
 integer(int64)::peak,active
 do
  read(*,*,iostat=ios)nfirst,ng,nk,nc,np,groups,divisor
  if(ios/=0)exit
  call choose_worker_layout(nfirst,ng,nk,nc,np,groups,divisor,rows,cols,b,chunk,peak,active,ierr)
  write(*,'(7(i0,1x))')rows,cols,b,chunk,peak,active,ierr
 enddo
end program worker_layout_probe
