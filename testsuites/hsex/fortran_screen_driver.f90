program test_screen
  use hsex_screen_math
  implicit none
  integer::n(3),boxes,ng,nk,u,status
  real(8)::vol,rec(3,3),head,qhead
  real(8),allocatable::q(:,:)
  complex(8),allocatable::chi(:,:),mat(:,:)
  open(newunit=u,file='input.bin',access='stream',form='unformatted',status='old')
  read(u)n,boxes,ng,nk,vol,rec
  allocate(q(3,ng),chi(ng,ng),mat(boxes**3,boxes**3))
  read(u)q,chi
  close(u)
  call screen_assemble(n,boxes,q,rec,vol,nk,chi,.true.,mat,head,qhead,status)
  if(status/=0)then
    write(*,'(a,i0)')'SCREEN_REJECTED status=',status
    stop
  endif
  open(newunit=u,file='output.bin',access='stream',form='unformatted',status='replace')
  write(u)mat
  close(u)
end program
