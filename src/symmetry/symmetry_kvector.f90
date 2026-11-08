module sym_kvector

  use sym_sub, only: SymMatB, use_symmetry, DISPLAY 

  implicit none

  private
  public :: init_sym_kvector

contains

  subroutine init_sym_kvector( kvec_io, weight_io, nkvec, blatvec )
    use salmon_global, only: dm_unfold_option
    use parallelization, only: end_parallel
    implicit none
    real(8),intent(inout) :: kvec_io(:,:), weight_io(:)
    integer,intent(inout) :: nkvec
    real(8),optional,intent(in) :: blatvec(3,3)
    real(8),allocatable :: kvec(:,:), weight(:)
    real(8) :: blatinv(3,3), transformed(3), delta(3)
    integer :: isym, ik, n_unique, i, representative

    if ( .not.use_symmetry ) return

    if( dm_unfold_option == 'super' .or. dm_unfold_option == 'primitive' ) then
      if ( DISPLAY ) write(*,"(A)") 'symmetry option not supported for dm_unfold calculation'
      call end_parallel
      stop
    end if

    if ( DISPLAY ) write(*,'(a60)') repeat("-",36)//" init_sym_kvector(start)"

    allocate( kvec(3,nkvec) ); kvec=0.0d0
    allocate( weight(nkvec) )

    if ( present(blatvec) ) then
       call get_inverse_lattice( blatvec, blatinv )
       kvec=matmul( blatinv, kvec_io )
    else
       kvec=kvec_io
    end if
    weight = weight_io

    kvec_io=0.0d0
    weight_io=0.0d0
    n_unique=0

    do ik=1,nkvec
       representative=0
       do isym=1,size(SymMatB,3)
          transformed=matmul(SymMatB(:,1:3,isym),kvec(:,ik))
          do i=1,n_unique
             delta=transformed-kvec(:,i)
             ! Reciprocal-lattice translations leave a k point equivalent.
             delta=delta-anint(delta)
             if (sum(delta**2)<1.d-5) then
                representative=i
                exit
             end if
          end do
          if (representative/=0) exit
       end do
       if (representative==0) then
          n_unique=n_unique+1
          representative=n_unique
          kvec(:,representative)=kvec(:,ik)
          call shift_vec_into_bz(kvec(:,representative))
       end if
       weight_io(representative)=weight_io(representative)+weight(ik)
    end do
    nkvec=n_unique
    if (present(blatvec)) then
       kvec_io(:,:nkvec)=matmul(blatvec,kvec(:,:nkvec))
    else
       kvec_io(:,:nkvec)=kvec(:,:nkvec)
    end if

    if ( DISPLAY ) then
       do i=1,nkvec
          write(*,'(i4,f10.5,2x,3f10.5)') i, weight_io(i), kvec_io(:,i)
       end do
       write(*,*) "sum(weight)=",sum(weight_io)
    end if

    deallocate( kvec )
    deallocate( weight )

    if ( DISPLAY ) write(*,'(a60)') repeat("-",38)//" init_sym_kvector(end)"

  end subroutine init_sym_kvector

  subroutine get_inverse_lattice( a, ainv )
    implicit none
    real(8),intent(in)  :: a(3,3)
    real(8),intent(out) :: ainv(3,3)
    real(8) :: b(3,3), v
    b(1,1) = a(2,2)*a(3,3) - a(3,2)*a(2,3)
    b(2,1) = a(3,2)*a(1,3) - a(1,2)*a(3,3)
    b(3,1) = a(1,2)*a(2,3) - a(2,2)*a(1,3)
    b(1,2) = a(2,3)*a(3,1) - a(3,3)*a(2,1)
    b(2,2) = a(3,3)*a(1,1) - a(1,3)*a(3,1)
    b(3,2) = a(1,3)*a(2,1) - a(2,3)*a(1,1)
    b(1,3) = a(2,1)*a(3,2) - a(3,1)*a(2,2)
    b(2,3) = a(3,1)*a(1,2) - a(1,1)*a(3,2)
    b(3,3) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
    v=a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1) &
     +a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1) &
     -a(1,2)*a(2,1)*a(3,3)-a(1,1)*a(2,3)*a(3,2)
    ainv(:,:) = transpose( b(:,:) )/v
  end subroutine get_inverse_lattice

  subroutine shift_vec_into_bz( vec )
    implicit none
    real(8),intent(inout) :: vec(:)
    integer :: i
    do i=1,size(vec)
       do
          if ( vec(i) > 0.5d0 ) then
             vec(i) = vec(i) - 1.0d0
          else if ( vec(i) <= -0.5d0 ) then
             vec(i) = vec(i) + 1.0d0
          else
             exit
          end if
       end do
    end do
  end subroutine shift_vec_into_bz

end module sym_kvector
