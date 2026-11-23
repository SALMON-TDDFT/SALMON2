program band_util

  implicit none

  integer,parameter :: unit1=5, unit2=1
  integer :: nband, nkpt, nblk
  integer :: iseg,iband,ikpt,idat,jdat
  real(8),allocatable :: esp(:,:),kpt(:,:)
  integer :: idummy(2)
  real(8) :: rdummy(3),d,p
  character(80) :: cbuf

  open(unit1,file='band.dat',status='old')

  read(unit1,*) cbuf, nband
  read(unit1,*) cbuf, nkpt
  read(unit1,*) cbuf, nblk

  allocate( kpt(3,nkpt*nblk) ); kpt=0.0d0
  allocate( esp(nband,nkpt*nblk) ); esp=0.0d0 

  idat=0
  jdat=0
  do iseg=1,nblk
     do ikpt=1,nkpt
        idat=idat+1
        read(unit1,*) idummy(1),rdummy(1:3),kpt(1:3,idat)
     end do
     do ikpt=1,nkpt
        jdat=jdat+1
        do iband=1,nband
           read(unit1,*) idummy(1:2),esp(iband,jdat)
        end do
     end do
  end do
  close(unit1)

  p=0.0d0
  do idat=1,nkpt*nblk
     if( idat > 1 )then
        d=sqrt( sum( (kpt(1:3,idat)-kpt(1:3,idat-1))**2 ) )
        p=p+d
     end if
     do iband=1,nband
        write(1000+iband,*) p, esp(iband,idat)
     end do
  end do

end program band_util
