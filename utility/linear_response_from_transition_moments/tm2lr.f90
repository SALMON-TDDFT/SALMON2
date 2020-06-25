! atomic unit
program main
  implicit none
  real(8),parameter :: pi=3.141592653589793d0
  complex(8),parameter :: zi=(0d0,1d0)
  !
  character(256) :: sysname
  integer :: num_kgrid(3),nstate,nelec
  real(8) :: al(3)
  integer :: mu,nu,nomega
  real(8) :: omega_max,delta
  namelist/input/sysname,num_kgrid,nstate,nelec,al,mu,nu,nomega,omega_max,delta
  !
  integer :: nk,nb
  real(8) :: n_e,V
  !
  character(256) :: filename
  real(8),allocatable :: eigen(:,:),occ(:,:)
  complex(8),allocatable :: p_matrix(:,:,:,:),q_matrix(:,:,:,:)
  integer :: ib1,ib2,ik,iw
  real(8) :: domega,w
  complex(8) :: sigma,eps,eps_inv

!===================================================================================================================================

 read(*,input)

  nk = num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
  nb = nstate
  V = al(1)*al(2)*al(3)* dble(nk) ! volume
  n_e = dble(nelec)/(al(1)*al(2)*al(3)) ! averaged electron number density
  
  domega = omega_max/dble(nomega)
  if(mu<1 .or. mu>3 .or. nu<1 .or. nu>3) then
    write(*,*) "error: mu & nu must be 1 or 2 or 3 (x,y,z)"
    stop
  end if

!===================================================================================================================================

  allocate(eigen(nb,nk),occ(nb,nk),p_matrix(nb,nb,3,nk),q_matrix(nb,nb,3,nk))
  call read_eigen_data
  call read_tm_data

  filename = trim(sysname) // '_sigma.data'
  open(777, file=filename, status='replace')
  write(777,*) "#1:omega[a.u.], 2:Re(sigma)[a.u.], 3:Im(sigma)[a.u.], 4:Re(epsilon), 5:Im(epsilon), 6:Re(1/epsilon), 7:Im(1/epsilon)"
  do iw=1,nomega
    w = dble(iw)*domega
    sigma = 0d0
    if(mu==nu) sigma = sigma + zi*n_e/w
    do ik=1,nk
      do ib1=1,nb
        if(occ(ib1,ik)==0d0) cycle
        do ib2=1,nb
          if(ib2==ib1) cycle
          sigma = sigma + (zi/(w*V))* occ(ib1,ik)*   &
          & ( q_matrix(ib1,ib2,mu,ik) * p_matrix(ib2,ib1,nu,ik) / ( w + eigen(ib1,ik) - eigen(ib2,ik) + zi*delta ) &
          & + q_matrix(ib2,ib1,mu,ik) * p_matrix(ib1,ib2,nu,ik) / (-w + eigen(ib1,ik) - eigen(ib2,ik) - zi*delta ) )
        end do
      end do
    end do
    eps = cmplx(1d0) + (4d0*pi*zi/w)* sigma
    eps_inv = cmplx(1d0)/eps
    write(777,'(7(1X,E23.15E3))') w,dble(sigma),imag(sigma),dble(eps),imag(eps),dble(eps_inv),imag(eps_inv)
  end do
  close(777)

!===================================================================================================================================

contains

  subroutine read_eigen_data
    implicit none
    integer,parameter :: fh=11
    character(256) :: dummy
    integer :: i, ik, ib, idummy, idummy2
    filename = trim(sysname) // '_eigen.data'
    open(fh,file=filename,status='old')
    do i=1, 3
      read(fh, '(a)') dummy !Skip
    end do
    do ik=1, nk
      read(fh, '(a)') dummy !Skip
      do ib=1, nb
        read(fh, *) idummy, eigen(ib, ik), occ(ib,ik)
      end do
    end do
    close(fh)
  end subroutine read_eigen_data

  subroutine read_tm_data
    implicit none
    integer,parameter :: fh=12
    character(256) :: dummy
    integer :: i, j, ik, ib, jb, idummy
    real(8) :: tmp(1:6)
    filename = trim(sysname) // '_tm.data'
    open(fh,file=filename,status='old')
    do i=1, 3
      read(fh, '(a)') dummy !Skip
    end do
    do ik=1, nk
      do ib=1, nb
        do jb=1, nb
          read(fh, *) idummy, idummy, idummy, tmp(1:6)
          p_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2)) ! <u_ib,k|p_1|u_jb,k>
          p_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4)) ! <u_ib,k|p_2|u_jb,k>
          p_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6)) ! <u_ib,k|p_3|u_jb,k>
        end do !jb
      end do !ib
    end do !ik
    read(fh, '(a)') dummy !Skip
    do ik=1,NK
      do ib=1,NB
        do i=1,3
          read(fh, '(a)') dummy !Skip
        enddo
      enddo
    enddo
    read(fh, '(a)') dummy !Skip
    do ik =1,NK
      do ib=1,NB
        do jb=1,NB
          read(fh, *) idummy, idummy, idummy,tmp(1:6)
          q_matrix(ib, jb, 1, ik) = p_matrix(ib, jb, 1, ik) - zi* dcmplx(tmp(1), tmp(2)) ! <u_ib,k|[r_1,V_nl]|u_jb,k>/i
          q_matrix(ib, jb, 2, ik) = p_matrix(ib, jb, 2, ik) - zi* dcmplx(tmp(3), tmp(4)) ! <u_ib,k|[r_2,V_nl]|u_jb,k>/i
          q_matrix(ib, jb, 3, ik) = p_matrix(ib, jb, 3, ik) - zi* dcmplx(tmp(5), tmp(6)) ! <u_ib,k|[r_3,V_nl]|u_jb,k>/i
        enddo
      enddo
    enddo
    close(fh)
  end subroutine read_tm_data

end program main
