module dg_hybrid_lcfo
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  private
  abstract interface
    subroutine dg_hybrid_lcfo_apply(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine dg_hybrid_lcfo_apply
  end interface
  public::assemble_dg_hybrid_lcfo_rows
contains
  subroutine assemble_dg_hybrid_lcfo_rows(comm,bases,row_ids,apply_h,apply_s,&
      hrows,srows,peak_elements,operator_fingerprint,ok,message)
    integer,intent(in)::comm
    type(s_dg_hybrid_fragment_basis),intent(in)::bases
    integer(int64),intent(in)::row_ids(:)
    procedure(dg_hybrid_lcfo_apply)::apply_h,apply_s
    complex(real64),allocatable,intent(out)::hrows(:,:),srows(:,:)
    integer(int64),intent(out)::peak_elements,operator_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,p,target,ierr,nlocal,global_count,global_point_count,local_bad,global_bad
    integer,allocatable::ownership(:)
    integer(int64)::local_count,global_count64,entry_hash,bits
    complex(real64),allocatable::vector(:,:),hvector(:,:),svector(:,:)
    complex(real64)::local_value,global_value
    logical::callback_ok

    ok=.false.;message='';peak_elements=0_int64;operator_fingerprint=0_int64
    if(.not.allocated(bases%global_ids).or..not.allocated(bases%buffer_point_ids).or.&
        .not.allocated(bases%buffer_values))then
      message='invalid distributed LCFO basis contract';return
    endif
    nlocal=size(bases%global_ids);local_count=int(nlocal,int64)
    call MPI_Allreduce(local_count,global_count64,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_count64<1_int64.or.global_count64>int(huge(global_count),int64))then
      message='invalid distributed LCFO basis extent';return
    endif
    global_count=int(global_count64);local_bad=0
    if(size(bases%buffer_values,1)/=size(bases%buffer_point_ids).or.&
        size(bases%buffer_values,2)/=nlocal.or.size(row_ids)/=nlocal)then
      local_bad=1
    elseif(any(row_ids/=bases%global_ids))then
      local_bad=1
    endif
    if(any(row_ids<1_int64).or.any(row_ids>global_count64))local_bad=1
    if(any(bases%buffer_point_ids<1_int64))local_bad=1
    if(.not.all(ieee_is_finite(real(bases%buffer_values))).or.&
        .not.all(ieee_is_finite(aimag(bases%buffer_values))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid distributed LCFO basis contract';return;endif
    global_point_count=0
    if(size(bases%buffer_point_ids)>0)global_point_count=int(maxval(bases%buffer_point_ids))
    call MPI_Allreduce(MPI_IN_PLACE,global_point_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_point_count<1)then;message='invalid LCFO physical point extent';return;endif
    allocate(ownership(global_count));ownership=0
    do i=1,nlocal;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then;message='LCFO basis rows are not owned exactly once';return;endif
    allocate(vector(size(bases%buffer_values,1),1),hvector(size(bases%buffer_values,1),1),&
      svector(size(bases%buffer_values,1),1),hrows(nlocal,global_count),srows(nlocal,global_count))
    hrows=(0d0,0d0);srows=(0d0,0d0)
    do j=1,global_count
      vector=(0d0,0d0)
      do target=1,global_point_count
        local_value=(0d0,0d0)
        do i=1,nlocal
          if(row_ids(i)/=int(j,int64))cycle
          do p=1,size(bases%buffer_point_ids)
            if(bases%buffer_point_ids(p)==int(target,int64))local_value=bases%buffer_values(p,i)
          enddo
        enddo
        call MPI_Allreduce(local_value,global_value,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='LCFO basis column redistribution failed';return;endif
        do p=1,size(bases%buffer_point_ids)
          if(bases%buffer_point_ids(p)==int(target,int64))vector(p,1)=global_value
        enddo
      enddo
      call apply_h(vector,hvector,callback_ok)
      call collective_callback_success(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='LCFO Hamiltonian application failed';return;endif
      call apply_s(vector,svector,callback_ok)
      call collective_callback_success(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='LCFO metric application failed';return;endif
      if(.not.all(ieee_is_finite(real(hvector))).or..not.all(ieee_is_finite(aimag(hvector))).or.&
          .not.all(ieee_is_finite(real(svector))).or..not.all(ieee_is_finite(aimag(svector))))then
        message='LCFO operator produced non-finite values';return
      endif
      do i=1,nlocal
        hrows(i,j)=sum(conjg(bases%buffer_values(:,i))*hvector(:,1))
        srows(i,j)=sum(conjg(bases%buffer_values(:,i))*svector(:,1))
      enddo
    enddo
    peak_elements=int(size(hrows)+size(srows),int64)
    entry_hash=0_int64
    do i=1,nlocal;do j=1,global_count
      bits=transfer(real(hrows(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,mod(7*int(row_ids(i))+11*j,63)))
      bits=transfer(aimag(hrows(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,mod(13*int(row_ids(i))+17*j,63)))
      bits=transfer(real(srows(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,mod(19*int(row_ids(i))+23*j,63)))
      bits=transfer(aimag(srows(i,j)),bits);entry_hash=ieor(entry_hash,ishftc(bits,mod(29*int(row_ids(i))+31*j,63)))
    enddo;enddo
    call MPI_Allreduce(entry_hash,operator_fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='LCFO operator fingerprint reduction failed';return;endif
    operator_fingerprint=ieor(operator_fingerprint,bases%provenance_fingerprint)
    if(operator_fingerprint==0_int64)operator_fingerprint=1709_int64
    ok=.true.;message=''
  contains
    subroutine collective_callback_success(local_ok,bad,status)
      logical,intent(in)::local_ok;integer,intent(out)::bad,status
      integer::local_value
      local_value=merge(0,1,local_ok)
      call MPI_Allreduce(local_value,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine collective_callback_success
  end subroutine assemble_dg_hybrid_lcfo_rows
end module dg_hybrid_lcfo
