! Transpose k-owned density matrices to pair owners, as in exx_k_exchange.
! Each pair owner performs both k transforms using the collinear HSE FFT plans.
#include "config.h"
module hsex_k_fft
  use exx_collinear_fft, only: collinear_fft_plans
  implicit none
  private
  public :: hsex_fft_init,hsex_forward,hsex_backward
  ! Private plans and transpose buffers persist across HSEX Hamiltonian calls.
  ! They belong to the single active system initialized by hsex_fft_init.
  type(collinear_fft_plans) :: plans
  integer,allocatable :: starts(:),counts(:),slots(:)
  integer :: rows,total_rows,nk,np,my_first,my_count,max_k
  complex(8),allocatable :: send(:),recv(:),work(:)
contains
  subroutine hsex_fft_init(nrows,mesh,cell,k,info,first,count)
    use structures, only: s_parallel_info
    use communication, only: comm_summation
    use exx_collinear_fft, only: collinear_fft_create
    implicit none
    integer,intent(in)::nrows,mesh(3)
    real(8),intent(in)::cell(3,3),k(:,:)
    type(s_parallel_info),intent(in)::info
    integer,intent(out)::first,count
    integer::layout(2,info%isize_ko),all_layout(2,info%isize_ko),ik,loc(3),ierr
    integer,allocatable::hits(:)
    real(8)::u(3)
    nk=size(k,2);np=info%isize_ko;total_rows=nrows
    layout=0;layout(:,info%id_ko+1)=[info%ik_s,info%numk]
    call comm_summation(layout,all_layout,size(layout),info%icomm_ko)
    starts=all_layout(1,:);counts=all_layout(2,:);max_k=maxval(counts)
    rows=(nrows+np-1)/np;my_first=1+info%id_ko*rows;my_count=max(0,min(rows,nrows-my_first+1))
    first=my_first;count=my_count
    allocate(slots(nk),hits(nk),send(rows*max_k*np),recv(rows*max_k*np),work(rows*nk))
    hits=0
    do ik=1,nk
      u=matmul(transpose(cell),k(:,ik)-k(:,1))*mesh/(2*acos(-1d0))
      if(maxval(abs(u-nint(u)))>1d-8)error stop 'HSEX FFT requires uniform k mesh'
      loc=modulo(nint(u),mesh)
      slots(ik)=1+loc(1)+mesh(1)*(loc(2)+mesh(2)*loc(3))
      hits(slots(ik))=hits(slots(ik))+1
    enddo
    if(any(hits/=1))error stop 'incomplete HSEX k FFT mesh'
    call collinear_fft_create(plans,mesh,nk,rows,my_count,.true.,work,ierr)
    if(ierr/=0)error stop 'HSEX k FFT plan failed'
  end subroutine

  subroutine exchange_tiles(comm)
#ifdef USE_MPI
    use communication, only: comm_alltoall
#endif
    implicit none
    integer,intent(in)::comm
    if(np==1)then
      recv=send
    else
#ifdef USE_MPI
      call comm_alltoall(send,recv,comm,rows*max_k)
#else
      error stop 'HSEX transpose requires MPI'
#endif
    endif
  end subroutine

  subroutine hsex_forward(pk,offset,comm,density)
    use exx_collinear_fft, only: collinear_fft_execute,FFTW_FORWARD_DIRECTION
    implicit none
    complex(8),intent(in)::pk(:,:),offset(:)
    integer,intent(in)::comm
    complex(8),intent(out)::density(:,:)
    integer::p,j,l,lo,m,i,threads
    ! A tile belongs to one destination rank and local k column. Padding is
    ! initialized by its owner, without a separate serial full-buffer clear.
    !$omp parallel do collapse(2) schedule(static) if(size(send)>=65536) private(lo,m,i)
    do p=1,np
      do j=1,max_k
        lo=(p-1)*rows+1;m=min(rows,total_rows-lo+1)
        i=rows*((p-1)*max_k+j-1)
        send(i+1:i+rows)=0d0
        if(m>0.and.j<=size(pk,2))send(i+1:i+m)=pk(lo:lo+m-1,j)
      enddo
    enddo
    !$omp end parallel do
    call exchange_tiles(comm)
    ! Each row owns a contiguous k FFT. Walk k inside the row rather than
    ! striding between FFTs for every received k column.
    !$omp parallel do schedule(static) if(size(work)>=65536) private(p,j,i)
    do l=1,my_count
      do p=1,np
        do j=1,counts(p)
          i=rows*((p-1)*max_k+j-1)
          work(slots(starts(p)+j-1)+(l-1)*nk)=recv(i+l)
        enddo
      enddo
    enddo
    !$omp end parallel do
    call collinear_fft_execute(plans,work,FFTW_FORWARD_DIRECTION,threads)
    !$omp parallel do schedule(static) if(size(work)>=65536)
    do l=1,my_count
      density(l,:)=work((l-1)*nk+1:l*nk)*conjg(offset)/nk
    enddo
    !$omp end parallel do
  end subroutine

  subroutine hsex_backward(kernel,offset,comm,act)
    use exx_collinear_fft, only: collinear_fft_execute,FFTW_BACKWARD_DIRECTION
    implicit none
    complex(8),intent(in)::kernel(:,:),offset(:)
    integer,intent(in)::comm
    complex(8),intent(out)::act(:,:)
    integer::p,j,l,lo,m,i,threads
    !$omp parallel do schedule(static) if(size(work)>=65536)
    do l=1,my_count
      work((l-1)*nk+1:l*nk)=kernel(l,:)*offset
    enddo
    !$omp end parallel do
    call collinear_fft_execute(plans,work,FFTW_BACKWARD_DIRECTION,threads)
    !$omp parallel do collapse(2) schedule(static) if(size(send)>=65536) private(i,l)
    do p=1,np
      do j=1,max_k
        i=rows*((p-1)*max_k+j-1)
        send(i+1:i+rows)=0d0
        if(j>counts(p))cycle
        do l=1,my_count
          send(i+l)=work(slots(starts(p)+j-1)+(l-1)*nk)
        enddo
      enddo
    enddo
    !$omp end parallel do
    call exchange_tiles(comm)
    !$omp parallel do collapse(2) schedule(static) if(size(recv)>=65536) private(lo,m,i)
    do p=1,np
      do j=1,size(act,2)
        lo=(p-1)*rows+1;m=min(rows,total_rows-lo+1)
        if(m<=0)cycle
        i=rows*((p-1)*max_k+j-1)
        act(lo:lo+m-1,j)=recv(i+1:i+m)
      enddo
    enddo
    !$omp end parallel do
  end subroutine
end module
