! Copyright 2026 SALMON developers. Apache-2.0.
! Compiler-paired Gamma native records→portable real64 streams. No solver.
program gs_gamma_exporter
  use iso_fortran_env,only:int8,int64,real64
  use, intrinsic :: ieee_arithmetic,only:ieee_is_finite
  use salmonpy_gs_stream,only:gs_write_chunk
  implicit none
  integer :: nx,ny,nz,nb,expected_iter,nk,mo,iter,nprocs,u,uw,ur,ui,ios,err,extra,expected_nk
  integer(int64) :: ng,bytes,i,factor,nc,idx,b,k
  real(real64) :: electrons,total,v,weight_sum
  real(real64),allocatable :: occ(:),real_wave(:),weights(:)
  complex(real64),allocatable :: wave(:)
  integer(int8),allocatable :: raw(:)
  logical :: real_orbital
  open(newunit=u,file='request.txt',status='old',action='read',iostat=ios)
  if(ios/=0)stop 1
  read(u,*,iostat=ios)nx,ny,nz,nb
  if(ios/=0)stop 2
  read(u,*,iostat=ios)electrons
  if(ios/=0)stop 3
  read(u,*,iostat=ios)expected_iter
  if(ios/=0)stop 4
  expected_nk=1
  read(u,*,iostat=ios)expected_nk
  if(ios>0)stop 43
  if(ios<0)expected_nk=1
  close(u)
  if(min(nx,ny,nz,nb,expected_nk)<1 .or. expected_iter<0)stop 5
  if(.not.ieee_is_finite(electrons))stop 6
  if(electrons<=0)stop 6
  ng=int(nx,int64)
  if(ng>huge(ng)/ny)stop 7
  ng=ng*ny
  if(ng>huge(ng)/nz)stop 8
  ng=ng*nz
  if(ng>huge(ng)/16/nb/expected_nk)stop 9
  open(newunit=u,file='native_input/info.bin',status='old',form='unformatted',iostat=ios)
  if(ios/=0)stop 10
  read(u,iostat=ios)nk
  if(ios/=0)stop 11
  read(u,iostat=ios)mo
  if(ios/=0)stop 12
  read(u,iostat=ios)iter
  if(ios/=0)stop 13
  read(u,iostat=ios)nprocs
  if(ios/=0)stop 14
  read(u,iostat=ios)real_orbital
  if(ios/=0)stop 15
  read(u,iostat=ios)extra
  if(ios>=0)stop 16
  close(u)
  if(nk/=expected_nk .or. mo/=nb .or. iter/=expected_iter .or. nprocs/=1)stop 17
  factor=16
  if(real_orbital)factor=8
  inquire(file='native_input/wfn.bin',size=bytes,iostat=ios)
  if(ios/=0)stop 18
  if(bytes/=factor*ng*nb*nk)stop 19
  nc=int(nb,int64)*nk
  allocate(occ(nc),weights(nk),raw(8*ng),stat=err)
  if(err/=0)stop 20
  if(real_orbital)then
    allocate(real_wave(ng),stat=err)
  else
    allocate(wave(ng),stat=err)
  end if
  if(err/=0)stop 20
  open(newunit=u,file='native_input/occupation.bin',status='old',form='unformatted',iostat=ios)
  if(ios/=0)stop 21
  read(u,iostat=ios)occ
  if(ios/=0)stop 22
  read(u,iostat=ios)extra
  if(ios>=0)stop 23
  close(u)
  weights=1.d0
  if(nk>1)then
    call read_weights()
  end if
  weight_sum=0
  do k=1,nk
    if(.not.ieee_is_finite(weights(k)))stop 45
    if(weights(k)<0)stop 46
    weight_sum=weight_sum+weights(k)
  end do
  if(.not.ieee_is_finite(weight_sum))stop 47
  if(abs(weight_sum-1.d0)>1.d-12)stop 48
  total=0
  do k=1,nk
  do b=1,nb
    idx=int(k-1,int64)*nb+b
    if(.not.ieee_is_finite(occ(idx)))stop 24
    if(occ(idx)<0 .or. occ(idx)>2)stop 25
    total=total+weights(k)*occ(idx)
  end do
  end do
  if(.not.ieee_is_finite(total))stop 49
  if(abs(total-electrons)>1.d-12*max(1.d0,electrons))stop 26
  call gs_write_chunk('portable_occupations.bin',occ,ios)
  if(ios/=0)stop 27
  open(newunit=uw,file='native_input/wfn.bin',status='old',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 28
  open(newunit=ur,file='portable_real.bin',status='new',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 29
  open(newunit=ui,file='portable_imag.bin',status='new',access='stream',form='unformatted',iostat=ios)
  if(ios/=0)stop 30
  do k=1,nk
  do b=1,nb
    if(real_orbital)then
      read(uw,iostat=ios)real_wave
      if(ios/=0)stop 31
    else
      read(uw,iostat=ios)wave
      if(ios/=0)stop 32
    end if
    do i=1,ng
      if(real_orbital)then
        v=real_wave(i)
      else
        v=real(wave(i),real64)
      end if
      if(.not.ieee_is_finite(v))stop 33
      raw(8*(i-1)+1:8*i)=transfer(v,raw(1:8))
    end do
    write(ur,iostat=ios)raw
    if(ios/=0)stop 34
    do i=1,ng
      v=0.d0
      if(.not.real_orbital)v=aimag(wave(i))
      if(.not.ieee_is_finite(v))stop 35
      raw(8*(i-1)+1:8*i)=transfer(v,raw(1:8))
    end do
    write(ui,iostat=ios)raw
    if(ios/=0)stop 36
  end do
  end do
  close(uw,iostat=ios)
  if(ios/=0)stop 37
  close(ur,iostat=ios)
  if(ios/=0)stop 38
  close(ui,iostat=ios)
  if(ios/=0)stop 39
  open(newunit=u,file='native_kind.txt',status='new',action='write',iostat=ios)
  if(ios/=0)stop 40
  if(real_orbital)then
    write(u,'(a)',iostat=ios)'real'
  else
    write(u,'(a)',iostat=ios)'complex'
  end if
  if(ios/=0)stop 41
  close(u,iostat=ios)
  if(ios/=0)stop 42
  print *, 'SALMONPY_GS_GAMMA_EXPORTER_V1 PASS'
contains
  subroutine read_weights()
    use salmonpy_gs_stream,only:gs_read_chunk
    call gs_read_chunk('portable_weights.bin',weights,ios)
    if(ios/=0)stop 44
  end subroutine
end program
