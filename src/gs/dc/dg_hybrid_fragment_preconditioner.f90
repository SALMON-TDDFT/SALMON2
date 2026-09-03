module dg_hybrid_fragment_preconditioner
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,&
    ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  implicit none
  private
  type,public::s_dg_hybrid_preconditioner_key
    integer::fragment_id=0,basis_generation=0,operator_epoch=0
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64,&
      operator_fingerprint=0_int64,reference_fingerprint=0_int64
  end type
  type,public::s_dg_hybrid_fragment_preconditioner
    private
    logical::valid=.false.
    integer::global_count=0,comm_rank=-1,comm_size=0
    type(s_dg_hybrid_preconditioner_key)::key
    integer(int64)::fingerprint=0_int64
    real(real64)::tolerance=0d0
    integer(int64),allocatable::row_ids(:)
    complex(real64),allocatable::q_rows(:,:)
    real(real64),allocatable::h_diagonal(:),s_diagonal(:)
  end type
  public::prepare_dg_hybrid_fragment_preconditioner,apply_dg_hybrid_fragment_preconditioner
contains
  subroutine prepare_dg_hybrid_fragment_preconditioner(comm,global_count,row_ids,q_rows,h_rows,s_rows,&
      key,tolerance,cache,fingerprint,ok,message)
    ! F=B Q is an immutable physical reference. Q changes covariantly with B;
    ! it is NOT the complete-system null-compression map or an occupied subset.
    integer,intent(in)::comm,global_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::q_rows(:,:),h_rows(:,:),s_rows(:,:)
    type(s_dg_hybrid_preconditioner_key),intent(in)::key
    real(real64),intent(in)::tolerance
    type(s_dg_hybrid_fragment_preconditioner),intent(inout)::cache
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_fragment_preconditioner)::working
    complex(real64),allocatable::column_data(:,:),overlap(:),action(:)
    integer,allocatable::owners(:)
    complex(real64)::diagonal(2)
    real(real64)::operator_scale(2),error
    integer::n,nr,i,j,p,stat,ierr
    integer(int64)::local_hash,global_hash,row_hash,words(7)
    logical::valid,halting(3)
    ok=.false.;message='invalid fixed-frame preconditioner contract';fingerprint=0_int64
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      valid=agree_integer(comm,global_count)
      valid=agree_key(comm,key).and.valid
      valid=agree_bits(comm,transfer(tolerance,0_int64)).and.valid
      n=global_count;nr=size(row_ids)
      if(.not.collective_valid(comm,valid.and.n>0.and.ieee_is_finite(tolerance)))return
      if(tolerance<64d0*epsilon(1d0).or.tolerance>1d-2)return
      valid=all(shape(q_rows)==[nr,n]).and.all(shape(h_rows)==[nr,n]).and.all(shape(s_rows)==[nr,n])
      valid=valid.and.all(row_ids>=1_int64).and.all(row_ids<=int(n,int64))
      valid=valid.and.finite(q_rows).and.finite(h_rows).and.finite(s_rows)
      valid=valid.and.int(n,int64)<=int(huge(0),int64)/3_int64
      if(.not.collective_valid(comm,valid))return
      allocate(owners(n),column_data(n,3),overlap(n),action(nr),working%row_ids(nr),&
        working%q_rows(nr,n),working%h_diagonal(n),working%s_diagonal(n),stat=stat)
      if(.not.collective_valid(comm,stat==0))then;message='cannot allocate fixed-frame workspace';return;endif
      owners=0
      do i=1,nr;owners(int(row_ids(i)))=owners(int(row_ids(i)))+1;enddo
      call MPI_Allreduce(MPI_IN_PLACE,owners,n,MPI_INTEGER,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.any(owners/=1))then;message='preconditioner rows must be owned exactly once';return;endif
      operator_scale=0d0
      if(nr>0)operator_scale=[maxval(abs(h_rows)),maxval(abs(s_rows))]
      call MPI_Allreduce(MPI_IN_PLACE,operator_scale,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      ! Stream one Q column and one H/S row at a time: O(n) replicated workspace,
      ! no full matrix gathering and no Hamiltonian diagonalization.
      do j=1,n
        column_data=0d0
        do i=1,nr
          p=int(row_ids(i));column_data(p,1)=q_rows(i,j)
          if(p==j)then;column_data(:,2)=h_rows(i,:);column_data(:,3)=s_rows(i,:);endif
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,column_data,3*n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(ierr/=MPI_SUCCESS)then;message='reference column transfer failed';return;endif
        error=0d0
        if(nr>0)error=maxval(abs(h_rows(:,j)-conjg(column_data(int(row_ids),2))))/max(1d0,operator_scale(1))
        valid=error<=tolerance
        error=0d0
        if(nr>0)error=maxval(abs(s_rows(:,j)-conjg(column_data(int(row_ids),3))))/max(1d0,operator_scale(2))
        if(.not.collective_valid(comm,valid.and.error<=tolerance))then
          message='non-Hermitian fixed-frame operator rows';return
        endif
        overlap=matmul(conjg(transpose(q_rows)),column_data(int(row_ids),1))
        call MPI_Allreduce(MPI_IN_PLACE,overlap,n,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        overlap(j)=overlap(j)-1d0
        if(.not.collective_valid(comm,ierr==MPI_SUCCESS.and.finite_vector(overlap)))return
        if(maxval(abs(overlap))>tolerance)then;message='fixed reference map is not unitary';return;endif
        action=matmul(h_rows,column_data(:,1));diagonal(1)=sum(conjg(q_rows(:,j))*action)
        action=matmul(s_rows,column_data(:,1));diagonal(2)=sum(conjg(q_rows(:,j))*action)
        call MPI_Allreduce(MPI_IN_PLACE,diagonal,2,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
        if(.not.collective_valid(comm,ierr==MPI_SUCCESS.and.finite_vector(diagonal)))then
          message='nonfinite reference diagonal';return
        endif
        if(any(abs(aimag(diagonal))>tolerance*max(1d0,operator_scale)))then
          message='complex reference expectation for Hermitian operator';return
        endif
        working%h_diagonal(j)=real(diagonal(1),real64);working%s_diagonal(j)=real(diagonal(2),real64)
      enddo
      if(any(working%s_diagonal<=tolerance*maxval(abs(working%s_diagonal))))then
        message='nonpositive or unresolved reference metric norm';return
      endif
      call MPI_Comm_rank(comm,working%comm_rank,ierr);valid=ierr==MPI_SUCCESS
      call MPI_Comm_size(comm,working%comm_size,ierr)
      if(.not.collective_valid(comm,valid.and.ierr==MPI_SUCCESS))return
      local_hash=0_int64
      do i=1,nr
        row_hash=mix(701_int64,row_ids(i))
        ! Bind ownership and row order, not just the partition-independent
        ! mathematical matrix: caches from different layouts must not mix.
        row_hash=mix(row_hash,int(working%comm_rank,int64));row_hash=mix(row_hash,int(i,int64))
        do j=1,n
          row_hash=mix_complex(row_hash,q_rows(i,j))
          row_hash=mix_complex(row_hash,h_rows(i,j))
          row_hash=mix_complex(row_hash,s_rows(i,j))
        enddo
        local_hash=ieor(local_hash,row_hash)
      enddo
      call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      words=key_words(key)
      do j=1,size(words);global_hash=mix(global_hash,words(j));enddo
      global_hash=mix(global_hash,int(n,int64));global_hash=mix(global_hash,int(working%comm_size,int64))
      global_hash=mix(global_hash,transfer(tolerance,0_int64))
      if(global_hash==0_int64)global_hash=1_int64
      working%global_count=n;working%key=key;working%tolerance=tolerance
      working%fingerprint=global_hash;working%q_rows=q_rows;working%row_ids=row_ids
      ! Two-phase publication: failed rebuilds cannot destroy a previous cache.
      call move_alloc(working%row_ids,cache%row_ids);call move_alloc(working%q_rows,cache%q_rows)
      call move_alloc(working%h_diagonal,cache%h_diagonal);call move_alloc(working%s_diagonal,cache%s_diagonal)
      cache%global_count=n;cache%comm_rank=working%comm_rank;cache%comm_size=working%comm_size
      cache%key=key;cache%tolerance=tolerance;cache%fingerprint=global_hash;cache%valid=.true.
      fingerprint=global_hash;ok=.true.;message=''
    end subroutine
  end subroutine prepare_dg_hybrid_fragment_preconditioner

  subroutine apply_dg_hybrid_fragment_preconditioner(comm,row_ids,key,cache,shifts,residual,output,&
      fingerprint,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    type(s_dg_hybrid_preconditioner_key),intent(in)::key
    type(s_dg_hybrid_fragment_preconditioner),intent(in)::cache
    real(real64),intent(in)::shifts(:)
    complex(real64),intent(in)::residual(:,:)
    complex(real64),allocatable,intent(out)::output(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::reference_residual(:,:),candidate(:,:)
    real(real64),allocatable::denominator(:)
    real(real64)::scale,roundoff,floor
    integer::n,nr,ns,i,j,stat,ierr,rank,nproc
    logical::valid,halting(3)
    ok=.false.;message='invalid or stale fixed-frame preconditioner';fingerprint=0_int64
    call suspend_traps(halting);call execute();call restore_traps(halting)
  contains
    subroutine execute()
      valid=agree_key(comm,key)
      valid=agree_integer(comm,size(shifts)).and.valid
      if(.not.collective_valid(comm,valid.and.cache%valid))return
      valid=agree_bits(comm,cache%fingerprint)
      valid=valid.and.all(key_words(key)==key_words(cache%key))
      call MPI_Comm_rank(comm,rank,ierr);valid=valid.and.ierr==MPI_SUCCESS
      call MPI_Comm_size(comm,nproc,ierr)
      valid=valid.and.ierr==MPI_SUCCESS.and.rank==cache%comm_rank.and.nproc==cache%comm_size
      n=cache%global_count;nr=size(row_ids);ns=size(shifts)
      valid=valid.and.nr==size(cache%row_ids).and.all(shape(residual)==[nr,ns]).and.ns>0.and.ns<=n
      valid=valid.and.finite(residual).and.all(ieee_is_finite(shifts))
      if(.not.collective_valid(comm,valid))return
      valid=all(row_ids==cache%row_ids)
      if(.not.collective_valid(comm,valid))then;message='preconditioner local row layout changed';return;endif
      do j=1,ns
        if(.not.agree_bits(comm,transfer(shifts(j),0_int64)))then;message='rank-disagreeing Rayleigh shifts';return;endif
      enddo
      if(int(n,int64)*int(ns,int64)>int(huge(0),int64))return
      allocate(reference_residual(n,ns),candidate(nr,ns),denominator(n),stat=stat)
      if(.not.collective_valid(comm,stat==0))then;message='cannot allocate preconditioner application';return;endif
      reference_residual=matmul(conjg(transpose(cache%q_rows)),residual)
      call MPI_Allreduce(MPI_IN_PLACE,reference_residual,n*ns,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(.not.collective_valid(comm,ierr==MPI_SUCCESS.and.finite(reference_residual)))then
        message='nonfinite transformed residual';return
      endif
      do j=1,ns
        scale=max(1d0,maxval(abs(cache%h_diagonal))+abs(shifts(j))*maxval(abs(cache%s_diagonal)))
        roundoff=64d0*epsilon(1d0)*real(n,real64)*scale
        floor=max(cache%tolerance*scale,roundoff)
        denominator=cache%h_diagonal-shifts(j)*cache%s_diagonal
        valid=ieee_is_finite(scale).and.ieee_is_finite(floor).and.all(ieee_is_finite(denominator))
        if(.not.collective_valid(comm,valid))then;message='nonfinite preconditioner denominator';return;endif
        do i=1,n
          if(abs(denominator(i))<=roundoff)then
            denominator(i)=floor
          else
            denominator(i)=sign(max(abs(denominator(i)),floor),denominator(i))
          endif
        enddo
        reference_residual(:,j)=reference_residual(:,j)/denominator
      enddo
      candidate=matmul(cache%q_rows,reference_residual)
      if(.not.collective_valid(comm,finite(reference_residual).and.finite(candidate)))then
        message='nonfinite preconditioned residual';return
      endif
      fingerprint=cache%fingerprint
      do j=1,ns;fingerprint=mix(fingerprint,transfer(shifts(j),0_int64));enddo
      if(fingerprint==0_int64)fingerprint=1_int64
      call move_alloc(candidate,output);ok=.true.;message=''
    end subroutine
  end subroutine apply_dg_hybrid_fragment_preconditioner

  function key_words(key) result(words)
    type(s_dg_hybrid_preconditioner_key),intent(in)::key
    integer(int64)::words(7)
    words=[int(key%fragment_id,int64),int(key%basis_generation,int64),int(key%operator_epoch,int64),&
      key%basis_fingerprint,key%metric_fingerprint,key%operator_fingerprint,key%reference_fingerprint]
  end function
  logical function agree_key(comm,key)
    integer,intent(in)::comm
    type(s_dg_hybrid_preconditioner_key),intent(in)::key
    integer(int64)::words(7),lo(7),hi(7)
    integer::ierr
    words=key_words(key)
    call MPI_Allreduce(words,lo,7,MPI_INTEGER8,MPI_MIN,comm,ierr);agree_key=ierr==MPI_SUCCESS
    call MPI_Allreduce(words,hi,7,MPI_INTEGER8,MPI_MAX,comm,ierr)
    agree_key=agree_key.and.ierr==MPI_SUCCESS.and.all(lo==hi).and.all(lo(:3)>0_int64).and.all(lo(4:)/=0_int64)
  end function
  logical function agree_integer(comm,value)
    integer,intent(in)::comm,value
    agree_integer=agree_bits(comm,int(value,int64))
  end function
  logical function agree_bits(comm,value)
    integer,intent(in)::comm
    integer(int64),intent(in)::value
    integer(int64)::lo,hi
    integer::ierr
    call MPI_Allreduce(value,lo,1,MPI_INTEGER8,MPI_MIN,comm,ierr);agree_bits=ierr==MPI_SUCCESS
    call MPI_Allreduce(value,hi,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    agree_bits=agree_bits.and.ierr==MPI_SUCCESS.and.lo==hi
  end function
  logical function collective_valid(comm,valid)
    integer,intent(in)::comm
    logical,intent(in)::valid
    integer::bad,ierr
    call MPI_Allreduce(merge(0,1,valid),bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    collective_valid=ierr==MPI_SUCCESS.and.bad==0
  end function
  logical function finite(x)
    complex(real64),intent(in)::x(:,:)
    finite=all(ieee_is_finite(real(x,real64))).and.all(ieee_is_finite(aimag(x)))
  end function
  logical function finite_vector(x)
    complex(real64),intent(in)::x(:)
    finite_vector=all(ieee_is_finite(real(x,real64))).and.all(ieee_is_finite(aimag(x)))
  end function
  integer(int64) function mix(seed,word) result(hash)
    integer(int64),intent(in)::seed,word
    hash=ieor(ishftc(seed,17),word)
    hash=ieor(hash,ishftc(iand(hash,ishftc(hash,23)),9))
    hash=ieor(hash,ishftc(hash,31))
  end function
  integer(int64) function mix_complex(seed,value) result(hash)
    integer(int64),intent(in)::seed
    complex(real64),intent(in)::value
    hash=mix(seed,transfer(real(value,real64),0_int64));hash=mix(hash,transfer(aimag(value),0_int64))
  end function
  subroutine suspend_traps(halting)
    logical,intent(out)::halting(3)
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
  end subroutine
  subroutine restore_traps(halting)
    logical,intent(in)::halting(3)
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1))
    call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  end subroutine
end module dg_hybrid_fragment_preconditioner
