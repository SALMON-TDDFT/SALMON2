! Offline Gamma localization. Default cutoff zero removes only exactly empty states.
! A positive cutoff is a separate occupation approximation, reported explicitly.
program localize_snapshot
  use iso_fortran_env, only: int32
  use hse_wannier
  implicit none
  type(s_hse_wannier) :: op
  integer(int32) :: header(14)
  real(8) :: metadata(9),cutoff,tolerance,exchange,k(3,1),discarded_occupation
  real(8),allocatable :: occupation(:,:),selected(:,:)
  complex(8),allocatable :: u(:,:,:),phi(:,:,:),psi(:,:,:),q(:,:),active(:,:,:),action(:,:,:)
  integer,allocatable :: indices(:)
  integer :: iu,status,n(3),ng,no,na,i,maxiter
  character(2048) :: input,output,argument
  call get_command_argument(1,input);call get_command_argument(2,output)
  call get_command_argument(3,argument);read(argument,*)cutoff
  call get_command_argument(4,argument);read(argument,*)maxiter
  call get_command_argument(5,argument);read(argument,*)tolerance
  if(cutoff<0d0.or.maxiter<0.or.tolerance<=0d0)error stop 'invalid localization controls'
  open(newunit=iu,file=trim(input),status='old',access='stream',form='unformatted')
  read(iu)header
  if(header(1)/=16909060.or.header(2)/=1)error stop 'unsupported snapshot version/endian'
  if(any(header(6:8)/=1))error stop 'offline localization currently requires Gamma'
  n=header(3:5);ng=product(n);no=header(9)
  read(iu)metadata
  allocate(occupation(no,1),u(no,no,1),phi(ng,no,1),q(ng,no),psi(ng,no,1))
  read(iu)occupation,u,phi,q
  close(iu)
  na=count(occupation(:,1)>cutoff)
  if(na<1)error stop 'no retained source states'
  allocate(indices(na),active(ng,na,1),selected(na,1),action(ng,na,1))
  indices=pack([(i,i=1,no)],occupation(:,1)>cutoff)
  psi(:,:,1)=matmul(phi(:,:,1),conjg(transpose(u(:,:,1))))
  active(:,:,1)=psi(:,indices,1);selected(:,1)=occupation(indices,1)
  discarded_occupation=sum(occupation,mask=occupation<=cutoff)
  write(*,'(a,2i8,a,es24.15)')'source states total/active: ',no,na,' dropped occupation: ', &
    discarded_occupation
  deallocate(u,phi,psi,q)
  k=0d0
  call wannier_init(op,n,[1,1,1],metadata(1:3),k,metadata(4),status)
  if(status/=0)error stop 'init failed'
  call wannier_localize(op,active,maxiter,tolerance,status)
  if(status/=0)error stop 'localization setup failed'
  call wannier_set_source(op,active,selected,op%gauge,status)
  if(status/=0)error stop 'source failed'
  call wannier_apply(op,active,action,status)
  if(status/=0)error stop 'exchange failed'
  exchange=0d0
  do i=1,na
    exchange=exchange+.125d0*selected(i,1)*op%dv*real(sum(conjg(active(:,i,1))*action(:,i,1)),8)
  enddo
  ! The original SCF flag cannot certify a density altered by positive cutoff.
  call wannier_snapshot(op,selected,metadata(4),exchange,metadata(9),int(header(13)), &
    header(14)==1.and.discarded_occupation==0d0,trim(output),status)
  if(status/=0)error stop 'write failed'
  write(*,'(a,i8,a,2es24.15)')'localization status: ',op%localization_status,' spread/gradient: ',op%spread,op%gradient
  write(*,'(a,2es24.15)')'exchange original/relocalized Ha: ',metadata(8),exchange
  call wannier_destroy(op)
end program
