! Copyright 2026 SALMON developers. Apache-2.0.
! Explicit compiler-bound single-rank Gamma restart writer. No SALMON solver.
program gs_gamma_converter
  use iso_fortran_env,only:int8,int64,real64
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use salmonpy_gs_stream,only:gs_read_chunk
  implicit none
  integer :: nx,ny,nz,nb,u,ur,ui,uw,ios,err,b
  integer(int64) :: ng,bytes,i
  real(real64) :: electrons,total,v
  real(real64),allocatable :: occ(:)
  complex(real64),allocatable :: wave(:)
  integer(int8),allocatable :: raw(:)
  logical :: real_orbital
  open(newunit=u,file='request.txt',status='old',action='read',iostat=ios)
  if(ios/=0)stop 1
  read(u,*,iostat=ios)nx,ny,nz,nb
  if(ios/=0)stop 2
  read(u,*,iostat=ios)electrons
  if(ios/=0)stop 3
  close(u)
  if(min(nx,ny,nz,nb)<1)stop 4
  if(.not.ieee_is_finite(electrons))stop 5
  ng=int(nx,int64)
  if(ng>huge(ng)/ny)stop 6
  ng=ng*ny
  if(ng>huge(ng)/nz)stop 7
  ng=ng*nz
  if(ng>huge(ng)/16/nb)stop 8
  inquire(file='portable_real.bin',size=bytes,iostat=ios)
  if(ios/=0)stop 13
  if(bytes/=8*ng*nb)stop 14
  inquire(file='portable_imag.bin',size=bytes,iostat=ios)
  if(ios/=0)stop 15
  if(bytes/=8*ng*nb)stop 16
  allocate(occ(nb),wave(ng),raw(8*ng),stat=err)
  if(err/=0)stop 9
  call gs_read_chunk('portable_occupations.bin',occ,ios)
  if(ios/=0)stop 10
  total=0
  do b=1,nb
    if(occ(b)<0 .or. occ(b)>2)stop 11
    total=total+occ(b)
  end do
  if(abs(total-electrons)>1.d-12*max(1.d0,electrons))stop 12
  open(newunit=ur,file='portable_real.bin',status='old',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 17
  open(newunit=ui,file='portable_imag.bin',status='old',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 18
  open(newunit=uw,file='data_for_restart/wfn.bin',status='new',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 19
  do b=1,nb
    read(ur,iostat=ios)raw
    if(ios/=0)stop 20
    do i=1,ng
      v=transfer(raw(8*(i-1)+1:8*i),v)
      if(.not.ieee_is_finite(v))stop 21
      wave(i)=cmplx(v,0.d0,real64)
    end do
    read(ui,iostat=ios)raw
    if(ios/=0)stop 22
    do i=1,ng
      v=transfer(raw(8*(i-1)+1:8*i),v)
      if(.not.ieee_is_finite(v))stop 23
      wave(i)=cmplx(real(wave(i),real64),v,real64)
    end do
    write(uw,iostat=ios)wave
    if(ios/=0)stop 24
  end do
  close(ur);close(ui);close(uw,iostat=ios)
  if(ios/=0)stop 25
  real_orbital=.false.
  open(newunit=u,file='data_for_restart/info.bin',status='new',form='unformatted',iostat=ios)
  if(ios/=0)stop 26
  write(u,iostat=ios)1
  if(ios/=0)stop 27
  write(u,iostat=ios)nb
  if(ios/=0)stop 28
  write(u,iostat=ios)0
  if(ios/=0)stop 29
  write(u,iostat=ios)1
  if(ios/=0)stop 30
  write(u,iostat=ios)real_orbital
  if(ios/=0)stop 31
  close(u,iostat=ios)
  if(ios/=0)stop 32
  open(newunit=u,file='data_for_restart/occupation.bin',status='new',form='unformatted',iostat=ios)
  if(ios/=0)stop 33
  write(u,iostat=ios)occ
  if(ios/=0)stop 34
  close(u,iostat=ios)
  if(ios/=0)stop 35
  print *, 'SALMONPY_GS_GAMMA_CONVERTER_V1 PASS'
end program
