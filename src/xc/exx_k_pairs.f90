! Hermitian real-space pair reuse on a full uniform k mesh. Only upper
! triangular density entries are transformed; targets need not be sources.
module exx_k_pairs
  use exx_collinear_fft,only:collinear_fft_create,collinear_fft_execute
  use iso_c_binding
  implicit none
  private
  include 'fftw3.f03'
  public :: apply_k_pairs
contains
  subroutine apply_k_pairs(op,source,target,action,starts,counts,rank,comm,transpose_tiles,ierr)
    use iso_fortran_env, only: int64
    use communication, only: comm_get_groupinfo,comm_summation,comm_bcast
    use exx_k_exchange, only: exx_k_kernel,execute_kernel_fft
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    complex(c_double_complex),intent(in) :: source(:,:,:),target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,count,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: count,comm
      end subroutine
    end interface
    complex(c_double_complex),parameter :: zero=(0d0,0d0),one=(1d0,0d0)
    complex(c_double_complex),allocatable :: s(:,:,:),t(:,:,:),reverse(:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:),rb(:,:,:),wk(:,:),kw(:,:),pack(:)
    complex(c_double_complex),pointer,contiguous :: flat_work(:)
    integer,allocatable :: rows(:,:),cols(:,:),active(:),inverse(:),slot(:),layout(:)
    integer :: np,id,ng,nk,batch,km,nmsg,nlocal,no,nt,bad,total_bad,p,j,i,a,r,c,ki,offset(3),dims(3)
    integer :: settings(13),expected(13),flags,tail,length,diagonal
    integer(int64) :: total_pairs,base,first,remaining,stride
    real(c_double) :: stamp
    external :: zgemm,zgeru
    ierr=1;action=zero
    op%seconds=0d0;stamp=0d0
    if(op%profile)stamp=walltime()
    call comm_get_groupinfo(comm,id,np)
    bad=0
    if(size(starts)/=np.or.size(counts)/=np.or.rank/=id)bad=1
    call comm_summation(bad,total_bad,comm)
    if(total_bad/=0)return
    ng=op%ng;nk=op%nk;no=size(source,2);nt=size(target,2)
    if(ng<1.or.nk<2.or.op%block<1.or.op%nq/=nk.or.allocated(op%accelerator))bad=1
    if(.not.c_associated(op%forward).or..not.c_associated(op%backward))bad=1
    if(.not.allocated(op%phase).or..not.allocated(op%work))bad=1
    if(any(counts<0).or.sum(int(counts,int64))/=int(nk,int64).or.starts(1)/=1)bad=1
    do p=2,np
      if(int(starts(p),int64)/=int(starts(p-1),int64)+counts(p-1))bad=1
    enddo
    nlocal=counts(id+1)
    if(min(no,nt)<1.or.size(source,1)/=ng.or.size(target,1)/=ng)bad=1
    if(size(source,3)/=nlocal.or.size(target,3)/=nlocal.or.any(shape(action)/=shape(target)))bad=1
    if(allocated(op%phase))then
      if(op%phase_start/=starts(id+1).or.size(op%phase,2)/=nlocal)bad=1
    endif
    if(.not.finite_complex(source).or..not.finite_complex(target))bad=1
    call comm_summation(bad,total_bad,comm)
    if(total_bad/=0)return
    settings=[op%n,op%mesh,ng,nk,op%nq,op%block,merge(1,0,op%contiguous_fft),no,nt]
    expected=settings
    call comm_bcast(expected,comm,0)
    if(any(settings/=expected))bad=1
    allocate(layout(2*np));layout=[starts,counts]
    call comm_bcast(layout,comm,0)
    if(any(layout/=[starts,counts]))bad=1
    km=maxval(counts)
    if(int(op%block,int64)*ng>int(huge(nmsg),int64)/np/km)bad=1
    call comm_summation(bad,total_bad,comm)
    if(total_bad/=0)return
    batch=op%block*ng;nmsg=batch*km
    flat_work(1:size(op%work))=>op%work
    total_pairs=int(ng,int64)*(int(ng,int64)+1)/2
    stride=int(np,int64)*batch
    ! Only one rank can own a partial final tile. Cache its exact FFT batch.
    remaining=modulo(total_pairs,stride)-int(id,int64)*batch
    tail=int(max(0_int64,min(int(batch,int64),remaining)))
    if(tail>0.and.tail<batch)then
      if(tail/=op%pair_tail_count.or..not.c_associated(op%pair_forward).or. &
         .not.c_associated(op%pair_backward).or..not.allocated(op%omp_pair_tail%forward))then
        if(c_associated(op%pair_forward))call fftw_destroy_plan(op%pair_forward)
        if(c_associated(op%pair_backward))call fftw_destroy_plan(op%pair_backward)
        dims=op%mesh(3:1:-1);flags=ior(FFTW_ESTIMATE,FFTW_UNALIGNED)
        if(op%contiguous_fft)then
          op%pair_forward=fftw_plan_many_dft(3,dims,tail,op%work,dims,1,nk, &
            op%work,dims,1,nk,FFTW_FORWARD,flags)
          op%pair_backward=fftw_plan_many_dft(3,dims,tail,op%work,dims,1,nk, &
            op%work,dims,1,nk,FFTW_BACKWARD,flags)
        else
          op%pair_forward=fftw_plan_many_dft(3,dims,tail,op%work,dims,batch,1, &
            op%work,dims,batch,1,FFTW_FORWARD,flags)
          op%pair_backward=fftw_plan_many_dft(3,dims,tail,op%work,dims,batch,1, &
            op%work,dims,batch,1,FFTW_BACKWARD,flags)
        endif
        call collinear_fft_create(op%omp_pair_tail,op%mesh,nk,batch,tail,op%contiguous_fft,flat_work,bad)
        if(bad==0)then
          op%pair_tail_count=tail
        else
          op%pair_tail_count=0
        endif
      endif
      if(.not.c_associated(op%pair_forward).or..not.c_associated(op%pair_backward))bad=1
    endif
    call comm_summation(bad,total_bad,comm)
    if(total_bad/=0)return
    allocate(send(nmsg*np),recv(nmsg*np),s(ng,no,nlocal),t(ng,nt,nlocal))
    allocate(reverse(batch))
    allocate(rows(batch,np),cols(batch,np),active(np),inverse(nk),slot(nk))
    sb(1:batch,1:km,1:np)=>send;rb(1:batch,1:km,1:np)=>recv
    wk(1:batch,1:nk)=>op%work;kw(1:nk,1:batch)=>op%work
    do ki=1,nk
      inverse(op%order(ki))=ki
    enddo
    do p=1,np
      do j=1,counts(p)
        slot(inverse(starts(p)+j-1))=j+(p-1)*km
      enddo
    enddo
    do j=1,nlocal
      do i=1,no
        s(:,i,j)=source(:,i,j)*op%phase(:,j)
      enddo
      do i=1,nt
        t(:,i,j)=target(:,i,j)*op%phase(:,j)
      enddo
    enddo
    if(op%profile)call mark(5)
    op%pair_fft_pairs=0_int64;op%threads_used=1
    do base=1_int64,total_pairs,stride
      send=zero
      do p=1,np
        first=base+int(p-1,int64)*batch
        active(p)=int(max(0_int64,min(int(batch,int64),total_pairs-first+1)))
        if(active(p)==0)cycle
        call unrank_pair(ng,first,r,c)
        do a=1,active(p)
          rows(a,p)=r;cols(a,p)=c
          c=c+1
          if(c>ng)then
            r=r+1;c=r
          endif
        enddo
        do j=1,nlocal
          pack(1:batch)=>send(((p-1)*km+j-1)*batch+1:((p-1)*km+j)*batch)
          a=1
          do while(a<=active(p))
            r=rows(a,p);c=cols(a,p);length=min(ng-c+1,active(p)-a+1)
            call zgemm('N','C',1,length,no,one,s(r,1,j),ng,s(c,1,j),ng,zero,pack(a),1)
            a=a+length
          enddo
        enddo
      enddo
      if(op%profile)call mark(1)
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark(2)
      do ki=1,nk
        p=(slot(ki)-1)/km+1;j=modulo(slot(ki)-1,km)+1
        if(op%contiguous_fft)then
          kw(ki,:)=rb(:,j,p)
        else
          wk(:,ki)=rb(:,j,p)
        endif
      enddo
      if(op%profile)call mark(5)
      if(active(id+1)>0)then
        if(active(id+1)<batch)then
          call collinear_fft_execute(op%omp_pair_tail,flat_work,FFTW_FORWARD,op%threads_used)
        else
          call execute_kernel_fft(op,op%forward)
        endif
        if(op%profile)call mark(3)
        !$omp parallel do private(r,c,ki,offset) schedule(static)
        do a=1,active(id+1)
          r=rows(a,id+1);c=cols(a,id+1)
          do ki=1,nk
            offset=modulo(op%point(:,r)-op%point(:,c)-op%shift(:,ki),op%n*op%mesh)
            if(op%contiguous_fft)then
              kw(ki,a)=kw(ki,a)*op%kernel(offset(1),offset(2),offset(3))
            else
              wk(a,ki)=wk(a,ki)*op%kernel(offset(1),offset(2),offset(3))
            endif
          enddo
        enddo
        !$omp end parallel do
        if(op%profile)call mark(4)
        if(active(id+1)<batch)then
          call collinear_fft_execute(op%omp_pair_tail,flat_work,FFTW_BACKWARD,op%threads_used)
        else
          call execute_kernel_fft(op,op%backward)
        endif
        if(op%profile)call mark(3)
        op%pair_fft_pairs=op%pair_fft_pairs+active(id+1)
      endif
      do ki=1,nk
        p=(slot(ki)-1)/km+1;j=modulo(slot(ki)-1,km)+1
        if(op%contiguous_fft)then
          sb(:,j,p)=kw(ki,:)
        else
          sb(:,j,p)=wk(:,ki)
        endif
      enddo
      if(op%profile)call mark(5)
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark(2)
      ! A contiguous upper-row segment contributes to both Hermitian rows.
      ! BLAS stays outside application OpenMP regions and owns its threading.
      do p=1,np
        do j=1,nlocal
          pack(1:batch)=>recv(((p-1)*km+j-1)*batch+1:((p-1)*km+j)*batch)
          a=1
          do while(a<=active(p))
            r=rows(a,p);c=cols(a,p);length=min(ng-c+1,active(p)-a+1)
            call zgemm('N','N',1,nt,length,-one/real(nk,c_double),pack(a),1, &
              t(c,1,j),ng,one,action(r,1,j),ng)
            diagonal=merge(1,0,r==c)
            if(length>diagonal)then
              reverse(1:length-diagonal)=conjg(pack(a+diagonal:a+length-1))
              call zgeru(length-diagonal,nt,-one/real(nk,c_double),reverse,1, &
                t(r,1,j),ng,action(c+diagonal,1,j),ng)
            endif
            a=a+length
          enddo
        enddo
      enddo
      if(op%profile)call mark(6)
    enddo
    do j=1,nlocal
      do i=1,nt
        action(:,i,j)=action(:,i,j)*conjg(op%phase(:,j))
      enddo
    enddo
    if(op%profile)call mark(5)
    ierr=0
  contains
    subroutine mark(stage)
      implicit none
      integer,intent(in) :: stage
      real(c_double) :: now
      now=walltime()
      op%seconds(stage)=op%seconds(stage)+now-stamp
      stamp=now
    end subroutine mark
  end subroutine apply_k_pairs

  real(c_double) function walltime()
    use iso_fortran_env, only: int64
    implicit none
    integer(int64) :: count,rate
    call system_clock(count,rate)
    walltime=real(count,c_double)/real(rate,c_double)
  end function walltime

  subroutine unrank_pair(ng,index,r,c)
    use iso_fortran_env, only: int64
    implicit none
    integer,intent(in) :: ng
    integer(int64),intent(in) :: index
    integer,intent(out) :: r,c
    integer(int64) :: lo,hi,mid,before
    lo=1;hi=ng
    do while(lo<hi)
      mid=(lo+hi)/2
      if(mid*(2_int64*ng-mid+1)/2>=index)then
        hi=mid
      else
        lo=mid+1
      endif
    enddo
    r=int(lo);before=(lo-1)*(2_int64*ng-lo+2)/2
    c=r+int(index-before-1)
  end subroutine unrank_pair

  logical function finite_complex(values) result(finite)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    complex(c_double_complex),intent(in) :: values(:,:,:)
    integer :: i,j,k
    finite=.false.
    do k=1,size(values,3)
      do j=1,size(values,2)
        do i=1,size(values,1)
          if(.not.ieee_is_finite(real(values(i,j,k))))return
          if(.not.ieee_is_finite(aimag(values(i,j,k))))return
        enddo
      enddo
    enddo
    finite=.true.
  end function finite_complex
end module exx_k_pairs
