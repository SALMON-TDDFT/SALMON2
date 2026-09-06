program test_automatic_distribution
  use automatic_distribution, only: select_numcpu_dft
  implicit none
  integer :: pk,po,pr(3),p,nk,no,ng(3),j,width,last
  logical :: found

  call expect(4,1,32,[29,28,10],.true.,.false.,1,2,[1,2,1])
  call expect(4,1,32,[17,16,16],.true.,.false.,1,2,[2,1,1])
  call expect(64,4,32,[16,16,16],.true.,.false.,4,16,[1,1,1])
  call expect(16,10,32,[16,16,16],.true.,.false.,8,2,[1,1,1])
  call expect(64,4,32,[32,32,32],.true.,.false.,4,2,[2,2,2])
  call expect(8,1,32,[64,16,16],.true.,.false.,1,2,[4,1,1])
  call expect(4,1,32,[64,64,64],.true.,.false.,1,1,[2,2,1])
  call expect(8,1,1,[16,16,16],.true.,.false.,1,1,[2,2,2])
  call expect(7,1,1,[28,8,8],.true.,.false.,1,1,[7,1,1])
  call expect(12,8,2,[16,16,16],.false.,.false.,6,2,[1,1,1])
  call expect(2,1,1,[15,16,15],.true.,.true.,1,1,[2,1,1])
  call select_numcpu_dft(7,1,1,[16,16,16],4,.true.,.false.,pk,po,pr,found)
  if(found) stop 1 ! Seven blocks cannot each supply a four-point halo.
  call select_numcpu_dft(17,4,4,[16,16,16],4,.false.,.false.,pk,po,pr,found)
  if(found) stop 2
  call select_numcpu_dft(0,1,1,[16,16,16],4,.true.,.false.,pk,po,pr,found)
  if(found) stop 3

  ! Verify invariants over uneven grids and rank counts, including primes.
  do p=1,96
    do nk=1,12
      do no=1,8
        ng=[19,28,35]
        call select_numcpu_dft(p,nk,no,ng,4,.true.,.false.,pk,po,pr,found)
        if(.not.found) cycle
        if(pk*po*product(pr)/=p.or.pk>nk.or.po>no.or.any(pr<1)) stop 4
        do j=1,3
          if(pr(j)==1) cycle
          width=(ng(j)-1)/pr(j)+1
          last=ng(j)-width*(pr(j)-1)
          if(last<4) stop 5
        end do
      end do
    end do
  end do
  print *, 'Automatic distribution tests passed'
contains
  subroutine expect(p,nk,no,ng,space,ffte,ek,eo,er)
    implicit none
    integer,intent(in) :: p,nk,no,ng(3),ek,eo,er(3)
    logical,intent(in) :: space,ffte
    call select_numcpu_dft(p,nk,no,ng,4,space,ffte,pk,po,pr,found)
    if(.not.found.or.pk/=ek.or.po/=eo.or.any(pr/=er)) then
      print *, 'Unexpected decomposition:',p,nk,no,ng,found,pk,po,pr
      print *, 'Expected:',ek,eo,er
      stop 6
    end if
  end subroutine expect
end program test_automatic_distribution
