program gen_plot

  implicit none
  integer :: i, unit,unit0,ix,iy,uo,ndata
  real(8) :: a, amin, amax
  character(12) :: fname

  unit0=4000

  amin= 1.d100
  amax=-1.d100

  write(*,*) "number of data = ?"
  read(*,*) ndata
! write(*,*) "the 0th unit = ?"
! read(*,*) unit0

  do i=1,ndata

     unit=unit0+i
     do
        read(unit,*,END=9) ix,iy,a
        amin = min( amin, a )
        amax = max( amax, a )
     end do
9    continue

     write(*,*) i, amin, amax

  end do

  uo=10
  open(uo,file="plot")
  write(uo,'("#set pm3d at b")')
  write(uo,'("#set view map")')
  write(uo,'("#unset surface")')
  write(uo,'("#set cbrange [",g14.5,":",g14.5,"]")') amin, amax
  write(uo,'("set zrange [",g14.5,":",g14.5,"]")') amin, amax
  write(uo,'(5x,a)') 
  do i=1,ndata
     write(fname,'("fort.",i4)') unit0+i
     fname = "'"//trim(fname)//"'"
     write(uo,'("splot ",a)') fname 
     write(uo,'("pause -1")')
  end do
  close(uo)

end program gen_plot
