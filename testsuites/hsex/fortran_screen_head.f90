program test_head
  use hsex_screen_math, only: screen_metal_head,screen_channels,screen_transform,screen_finite
  implicit none
  integer::qm(3),mesh(3),na,u,status,nsample
  integer,allocatable :: samples(:)
  character(1024) :: input_path
  character(4) :: magic
  real(8)::rec(3,3),tolerance,error
  complex(8),allocatable::values(:,:,:),point(:,:),average(:,:),x(:,:)
  if(command_argument_count()==0)then
    open(newunit=u,file='head.bin',access='stream',form='unformatted',status='old')
    read(u)qm,mesh,na,rec,tolerance
    allocate(values(na,na,product(qm)),point(na,na),average(na,na))
    read(u)values
  else
    call get_command_argument(1,input_path)
    open(newunit=u,file=trim(input_path),access='stream',form='unformatted',status='old')
    read(u)magic,qm,mesh,na,nsample,rec,tolerance
    if(magic/='HXQ1'.or.nsample<1.or.nsample>product(qm).or.na<1)error stop 'invalid saved head samples'
    allocate(values(na,na,nsample),point(na,na),average(na,na),samples(nsample))
    read(u)samples,values
  endif
  close(u)
  call screen_metal_head(values,qm,rec,mesh,tolerance,point,average,status,samples)
  if(allocated(samples))write(*,'(a,i0)')'HSEX_HEAD_REPLAY_STATUS ',status
  open(newunit=u,file='head-output.bin',access='stream',form='unformatted',status='replace')
  write(u)status
  if(status==0)write(u)point,average
  close(u)
end program
