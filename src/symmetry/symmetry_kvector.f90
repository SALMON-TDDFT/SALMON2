module sym_kvector

  use sym_sub, only: SymMatB, use_symmetry 

  implicit none

  private
  public :: init_sym_kvector

contains

  subroutine init_sym_kvector( kvec_io, weight, nkvec, blatvec )
    implicit none
    real(8),intent(inout) :: kvec_io(:,:), weight(:)
    integer,intent(inout) :: nkvec
    real(8),optional,intent(in) :: blatvec(3,3)
    real(8),allocatable :: kvec(:,:), k_list(:,:), Skvec(:,:)
    real(8) :: blatinv(3,3)
    integer :: isym, nsym, ik, n_k_list,i_list_0,i
    integer,allocatable :: i_list(:)

    if ( .not.use_symmetry ) return

    write(*,'(a60)') repeat("-",36)//" init_sym_kvector(start)"

    allocate( kvec(3,nkvec) ); kvec=0.0d0

    if ( present(blatvec) ) then
       call get_inverse_lattice( blatvec, blatinv )
       kvec=matmul( blatinv, kvec_io )
    else
       kvec=kvec_io
    end if

    nsym=size(SymMatB,3)

    allocate( Skvec(3,nsym)        ); Skvec=0.0d0
    allocate( k_list(3,nsym*nkvec) ); k_list=0.0d0
    allocate( i_list(nsym*nkvec)   ); i_list=0.0d0
    n_k_list=0

    do ik=1,nkvec
       do isym=1,nsym
          Skvec(:,isym)=matmul( SymMatB(:,1:3,isym), kvec(:,ik) )
          call shift_vec_into_bz( Skvec(:,isym) )
       end do ! isym
       call check_k_list( n_k_list, k_list, i_list, Skvec )
    end do ! ik

    kvec_io=0.0d0
    weight=0.0d0

    i_list_0=0
    i=0
    do ik=1,n_k_list
       if ( i_list_0 /= i_list(ik) ) then
          i_list_0=i_list(ik)
          i=i+1
          if ( present(blatvec) ) then
             kvec_io(:,i)=matmul( blatvec, k_list(:,i_list_0) )
          else
             kvec_io(:,i)=k_list(:,i_list_0)
          end if
          weight(i)=dble(count(i_list==i_list_0))/dble(nkvec)
       end if
    end do
    nkvec=i

    do i=1,nkvec
       write(*,'(i4,f10.5,2x,3f10.5)') i, weight(i), kvec_io(:,i)
    end do
    write(*,*) "sum(weight)=",sum(weight)

    deallocate( i_list )
    deallocate( k_list )
    deallocate( Skvec )
    deallocate( kvec )

    write(*,'(a60)') repeat("-",38)//" init_sym_kvector(end)"

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

  subroutine check_k_list( n_k_list, k_list, i_list, kvec )
    implicit none
    integer,intent(inout) :: n_k_list
    real(8),intent(inout) :: k_list(:,:)
    integer,intent(inout) :: i_list(:)
    real(8),intent(in) :: kvec(:,:)
    real(8) :: d
    integer :: i,isym,i_k_list

    do i=1,n_k_list
       d=sum((k_list(:,i)-kvec(:,1))**2)
       if ( d < 1.d-5 ) return
    end do

    n_k_list = n_k_list + 1
    i_k_list = n_k_list

    k_list(:,n_k_list) = kvec(:,1)

    i_list(n_k_list)   = i_k_list

    do isym=2,size(kvec,2)
       do i=1,n_k_list
          d=sum((k_list(:,i)-kvec(:,isym))**2)
          if ( d < 1.d-5 ) exit
       end do ! i
       if ( i > n_k_list ) then
          n_k_list = n_k_list + 1
          k_list(:,n_k_list) = kvec(:,isym)
          i_list(n_k_list) = i_k_list
       end if
    end do ! isym

  end subroutine check_k_list

end module sym_kvector
