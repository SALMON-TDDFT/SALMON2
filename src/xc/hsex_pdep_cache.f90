! A dielectric eigensubspace is only a guess for the current Hamiltonian.
module hsex_pdep_cache
  implicit none
  private
  public :: pdep_cache_read,pdep_cache_write
contains
  subroutine pdep_cache_read(path,n,g,cell,cutoff,nv,basis,comm,rank,status)
    use hsex_screen_math, only: screen_finite
    use communication, only: comm_bcast,comm_is_root
    implicit none
    character(*),intent(in) :: path
    integer,intent(in) :: n(3),g(:,:),nv,comm,rank
    real(8),intent(in) :: cell(3,3),cutoff
    complex(8),intent(out) :: basis(:,:)
    integer,intent(out) :: status
    integer :: u,ios,close_code,header(6),stored_g(3,size(g,2)),j
    real(8) :: stored_cell(3,3),stored_cutoff
    complex(8),allocatable :: gram(:,:)
    character(4) :: magic
    basis=0d0;status=0;ios=0
    if(comm_is_root(rank))then
      open(newunit=u,file=path,status='old',access='stream',form='unformatted',iostat=ios)
      if(ios==0)then
        read(u,iostat=ios)magic,header
        if(ios==0)then
          if(magic/='GPC1'.or.any(header/=[n,size(g,2),size(basis,2),nv]))ios=1
        endif
        if(ios==0)then
          read(u,iostat=ios)stored_g,stored_cell,stored_cutoff,basis
          if(ios==0)then
            if(any(stored_g/=g).or.any(stored_cell/=cell).or.stored_cutoff/=cutoff)ios=1
            if(.not.screen_finite(basis))ios=1
            if(maxval(abs(aimag(basis)))>1d-10)ios=1
            if(ios==0)then
              gram=matmul(conjg(transpose(basis)),basis)
              do j=1,size(basis,2)
                gram(j,j)=gram(j,j)-1d0
              enddo
              if(.not.screen_finite(gram))ios=1
              if(maxval(abs(gram))>1d-8)ios=1
            endif
          endif
        endif
        close(u,iostat=close_code)
        if(ios==0)ios=close_code
      endif
    endif
    call comm_bcast(ios,comm,0)
    status=ios
    if(status==0)call comm_bcast(basis,comm,0)
  end subroutine

  subroutine pdep_cache_write(path,n,g,cell,cutoff,nv,basis,comm,rank,status)
    use communication, only: comm_bcast,comm_is_root
    use hsex_screen_math, only: screen_finite
    use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
    implicit none
    character(*),intent(in) :: path
    integer,intent(in) :: n(3),g(:,:),nv,comm,rank
    real(8),intent(in) :: cell(3,3),cutoff
    complex(8),intent(in) :: basis(:,:)
    integer,intent(out) :: status
    integer :: u,ios,close_code
    logical :: exists
    interface
      function link_file(old,new) bind(C,name='link') result(code)
        import c_int,c_char
        implicit none
        character(c_char),intent(in) :: old(*),new(*)
        integer(c_int) :: code
      end function
      function unlink_file(path) bind(C,name='unlink') result(code)
        import c_int,c_char
        implicit none
        character(c_char),intent(in) :: path(*)
        integer(c_int) :: code
      end function
    end interface
    ios=0
    if(comm_is_root(rank))then
      if(.not.screen_finite(basis))ios=1
      inquire(file=path,exist=exists)
      if(exists)ios=1
      if(ios==0)then
        open(newunit=u,file=trim(path)//'.partial',status='new',access='stream',form='unformatted',iostat=ios)
        if(ios==0)then
          write(u,iostat=ios)'GPC1',n,size(g,2),size(basis,2),nv,g,cell,cutoff,basis
          close(u,iostat=close_code)
          if(ios==0)ios=close_code
        endif
        if(ios==0)ios=link_file(trim(path)//'.partial'//c_null_char,trim(path)//c_null_char)
        if(ios==0)ios=unlink_file(trim(path)//'.partial'//c_null_char)
      endif
    endif
    call comm_bcast(ios,comm,0);status=ios
  end subroutine
end module
