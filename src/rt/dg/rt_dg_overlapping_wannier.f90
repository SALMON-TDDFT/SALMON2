#include "config.h"
module rt_dg_overlapping_wannier
  use iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_invalid,ieee_divide_by_zero,&
    ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_overlapping_wannier_rt_state
    integer::nwann=0,basis_generation=0,geometry_generation=0,step=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64
    integer(int64)::hamiltonian_fingerprint=0_int64,observable_fingerprint=0_int64,fingerprint=0_int64
    character(64)::field_coupling_convention=''
    real(real64)::time=0d0
    complex(real64),allocatable::metric(:,:),hamiltonian0(:,:),position(:,:,:),velocity(:,:,:)
    complex(real64),allocatable::field_free_vectors(:,:)
    real(real64),allocatable::field_free_eigenvalues(:)
  end type
  public::initialize_dg_overlapping_wannier_rt,advance_dg_overlapping_wannier_rt
  public::write_dg_overlapping_wannier_rt_restart,read_dg_overlapping_wannier_rt_restart
contains
  subroutine initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,&
      basis_generation,geometry_generation,basis_fingerprint,operator_fingerprint,&
      hamiltonian_fingerprint,observable_fingerprint,field_coupling_convention,&
      expected_basis_generation,expected_geometry_generation,expected_basis_fingerprint,&
      expected_operator_fingerprint,coefficients,state,ok,message)
    integer,intent(in)::comm,row_ids(:),basis_generation,geometry_generation,&
      expected_basis_generation,expected_geometry_generation
    integer(int64),intent(in)::basis_fingerprint,operator_fingerprint,&
      hamiltonian_fingerprint,observable_fingerprint,&
      expected_basis_fingerprint,expected_operator_fingerprint
    character(*),intent(in)::field_coupling_convention
    complex(real64),intent(in)::srows(:,:),hrows(:,:),xrows(:,:,:),vrows(:,:,:),coefficients(:,:)
    type(s_dg_overlapping_wannier_rt_state),intent(out)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nwann,nlocal,i,j,axis,ierr,local_bad,global_bad,rank,nproc,r,nblock
    integer,allocatable::ownership(:),block_ids(:)
    integer::integer_metadata(4),integer_min(4),integer_max(4)
    integer(int64)::fingerprints(6),fingerprint_min(6),fingerprint_max(6)
    character(len(field_coupling_convention))::convention_reference
    complex(real64),allocatable::hcopy(:,:),scopy(:,:),block2(:,:),block3(:,:,:),coefficient_reference(:,:)
    real(real64)::defect,scale,coefficient_defect

    state=s_dg_overlapping_wannier_rt_state();ok=.false.;message=''
    call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
    nwann=size(srows,2);nlocal=size(row_ids)
    integer_metadata=[basis_generation,geometry_generation,expected_basis_generation,&
      expected_geometry_generation]
    fingerprints=[basis_fingerprint,operator_fingerprint,hamiltonian_fingerprint,&
      observable_fingerprint,expected_basis_fingerprint,expected_operator_fingerprint]
    call MPI_Allreduce(integer_metadata,integer_min,4,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(integer_metadata,integer_max,4,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(fingerprints,fingerprint_min,6,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(fingerprints,fingerprint_max,6,MPI_INTEGER8,MPI_MAX,comm,ierr)
    local_bad=0
    if(any(integer_min/=integer_max).or.any(fingerprint_min/=fingerprint_max))local_bad=1
    if(basis_generation/=expected_basis_generation.or.geometry_generation/=expected_geometry_generation.or.&
       basis_fingerprint/=expected_basis_fingerprint.or.operator_fingerprint/=expected_operator_fingerprint)&
      local_bad=1
    convention_reference=field_coupling_convention
    call MPI_Bcast(convention_reference,len(convention_reference),MPI_CHARACTER,0,comm,ierr)
    if(convention_reference/=field_coupling_convention.or.hamiltonian_fingerprint==0_int64.or.&
       observable_fingerprint==0_int64.or.&
       trim(field_coupling_convention)/='cell_wrapped_length_velocity')local_bad=1
    if(nwann<1.or.size(hrows,1)/=nlocal.or.size(hrows,2)/=nwann.or.&
       size(srows,1)/=nlocal.or.any(shape(xrows)/=[3,nlocal,nwann]).or.&
       any(shape(vrows)/=[3,nlocal,nwann]).or.size(coefficients,1)/=nwann.or.&
       size(coefficients,2)<1)local_bad=1
    if(any(row_ids<1).or.any(row_ids>nwann))local_bad=1
    do i=1,nlocal
      if(count(row_ids==row_ids(i))/=1)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid or stale overlapping-Wannier RT initialization';return;endif
    allocate(coefficient_reference(size(coefficients,1),size(coefficients,2)))
    coefficient_reference=coefficients
    call MPI_Bcast(coefficient_reference,size(coefficient_reference),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    coefficient_defect=maxval(abs(coefficients-coefficient_reference))
    call MPI_Allreduce(MPI_IN_PLACE,coefficient_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(coefficient_defect/=0d0)then
      message='rank-inconsistent overlapping-Wannier initial coefficients';return
    endif

    allocate(ownership(nwann));ownership=0
    do i=1,nlocal
      ownership(row_ids(i))=1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,nwann,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(rank==0)then
      allocate(state%metric(nwann,nwann),state%hamiltonian0(nwann,nwann),&
        state%position(3,nwann,nwann),state%velocity(3,nwann,nwann))
      state%metric=(0d0,0d0);state%hamiltonian0=(0d0,0d0)
      state%position=(0d0,0d0);state%velocity=(0d0,0d0)
    else
      allocate(state%metric(0,0),state%hamiltonian0(0,0),state%position(3,0,0),state%velocity(3,0,0))
    endif
    do r=0,nproc-1
      nblock=nlocal
      call MPI_Bcast(nblock,1,MPI_INTEGER,r,comm,ierr)
      allocate(block_ids(nblock),block2(nblock,nwann),block3(3,nblock,nwann))
      if(rank==r)block_ids=row_ids
      call MPI_Bcast(block_ids,nblock,MPI_INTEGER,r,comm,ierr)
      if(rank==r)block2=srows
      call MPI_Bcast(block2,nblock*nwann,MPI_DOUBLE_COMPLEX,r,comm,ierr)
      if(rank==0)then
        do i=1,nblock;state%metric(block_ids(i),:)=block2(i,:);enddo
      endif
      if(rank==r)block2=hrows
      call MPI_Bcast(block2,nblock*nwann,MPI_DOUBLE_COMPLEX,r,comm,ierr)
      if(rank==0)then
        do i=1,nblock;state%hamiltonian0(block_ids(i),:)=block2(i,:);enddo
      endif
      if(rank==r)block3=xrows
      call MPI_Bcast(block3,3*nblock*nwann,MPI_DOUBLE_COMPLEX,r,comm,ierr)
      if(rank==0)then
        do i=1,nblock;state%position(:,block_ids(i),:)=block3(:,i,:);enddo
      endif
      if(rank==r)block3=vrows
      call MPI_Bcast(block3,3*nblock*nwann,MPI_DOUBLE_COMPLEX,r,comm,ierr)
      if(rank==0)then
        do i=1,nblock;state%velocity(:,block_ids(i),:)=block3(:,i,:);enddo
      endif
      deallocate(block_ids,block2,block3)
    enddo
    local_bad=merge(0,1,all(ownership==1).and.finite_complex_2d(coefficients))
    if(rank==0)then
      if(.not.finite_complex_2d(state%metric).or..not.finite_complex_2d(state%hamiltonian0).or.&
         .not.finite_complex_3d(state%position).or..not.finite_complex_3d(state%velocity))local_bad=1
      scale=max(1d0,maxval(abs(state%metric)),maxval(abs(state%hamiltonian0)))
      defect=max(maxval(abs(state%metric-conjg(transpose(state%metric)))),&
        maxval(abs(state%hamiltonian0-conjg(transpose(state%hamiltonian0)))))
      do axis=1,3
        defect=max(defect,maxval(abs(state%position(axis,:,:)-&
          conjg(transpose(state%position(axis,:,:))))))
        defect=max(defect,maxval(abs(state%velocity(axis,:,:)-&
          conjg(transpose(state%velocity(axis,:,:))))))
      enddo
      if(defect>1d-11*scale)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid or non-Hermitian overlapping-Wannier RT matrices';return;endif

    if(rank==0)then
      hcopy=state%hamiltonian0;scopy=state%metric
      call generalized_hermitian_eigensystem(hcopy,scopy,state%field_free_eigenvalues,&
        state%field_free_vectors,ok,message)
    else
      allocate(state%field_free_eigenvalues(0),state%field_free_vectors(0,0));ok=.false.
    endif
    local_bad=merge(0,1,rank==0.and.ok)
    call MPI_Bcast(local_bad,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(message,len(message),MPI_CHARACTER,0,comm,ierr)
    ok=local_bad==0
    if(.not.ok)return
    state%nwann=nwann;state%basis_generation=basis_generation
    state%geometry_generation=geometry_generation
    state%basis_fingerprint=basis_fingerprint;state%operator_fingerprint=operator_fingerprint
    state%hamiltonian_fingerprint=hamiltonian_fingerprint
    state%observable_fingerprint=observable_fingerprint
    state%field_coupling_convention=field_coupling_convention
    state%fingerprint=ieor(ieor(ieor(basis_fingerprint,operator_fingerprint),&
      ieor(hamiltonian_fingerprint,observable_fingerprint)),&
      int(104729*nwann+1009*basis_generation+9176*geometry_generation,int64))
    state%step=0;state%time=0d0;ok=.true.
#else
    state=s_dg_overlapping_wannier_rt_state();ok=.false.
    message='overlapping-Wannier coefficient RT requires MPI'
#endif
  end subroutine

  subroutine advance_dg_overlapping_wannier_rt(comm,dt,electric_field,vector_potential,&
      coefficients,state,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::dt,electric_field(3),vector_potential(3)
    complex(real64),intent(inout)::coefficients(:,:)
    type(s_dg_overlapping_wannier_rt_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::axis,i,j,ierr,local_bad,global_bad,rank
    complex(real64),allocatable::hmid(:,:),scopy(:,:),vectors(:,:),metric_coefficients(:,:),modal(:,:)
    real(real64),allocatable::eigenvalues(:)
    real(real64)::contract(7),contract_min(7),contract_max(7)
    logical::field_free
    ok=.false.;message='';local_bad=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(state%nwann<1.or.dt<=0d0.or..not.ieee_is_finite(dt).or.&
       any(.not.ieee_is_finite(electric_field)).or.any(.not.ieee_is_finite(vector_potential)))local_bad=1
    if(size(coefficients,1)/=state%nwann.or.size(coefficients,2)<1)local_bad=1
    if(.not.allocated(state%metric).or..not.allocated(state%field_free_vectors))local_bad=1
    if(rank==0)then
      if(any(shape(state%metric)/=[state%nwann,state%nwann]).or.&
         any(shape(state%field_free_vectors)/=[state%nwann,state%nwann]))local_bad=1
    else
      if(size(state%metric)/=0.or.size(state%field_free_vectors)/=0)local_bad=1
    endif
    contract=[dt,electric_field,vector_potential]
    call MPI_Allreduce(contract,contract_min,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(contract,contract_max,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(any(contract_min/=contract_max))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid overlapping-Wannier RT step contract';return;endif
    if(rank==0)then
      field_free=all(electric_field==0d0).and.all(vector_potential==0d0)
      if(field_free)then
        vectors=state%field_free_vectors;eigenvalues=state%field_free_eigenvalues
      else
        hmid=state%hamiltonian0
        do axis=1,3
          hmid=hmid+electric_field(axis)*state%position(axis,:,:)+&
            vector_potential(axis)*state%velocity(axis,:,:)
        enddo
        scopy=state%metric
        call generalized_hermitian_eigensystem(hmid,scopy,eigenvalues,vectors,ok,message)
      endif
      if(field_free.or.ok)then
        metric_coefficients=matmul(state%metric,coefficients)
        modal=matmul(conjg(transpose(vectors)),metric_coefficients)
        do j=1,size(modal,2)
          do i=1,state%nwann
            modal(i,j)=exp(cmplx(0d0,-dt*eigenvalues(i),real64))*modal(i,j)
          enddo
        enddo
        coefficients=matmul(vectors,modal);local_bad=0
      else
        local_bad=1
      endif
    endif
    call MPI_Bcast(local_bad,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(message,len(message),MPI_CHARACTER,0,comm,ierr)
    if(local_bad/=0)return
    call MPI_Bcast(coefficients,size(coefficients),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    state%step=state%step+1;state%time=state%time+dt;ok=.true.
#else
    ok=.false.;message='overlapping-Wannier coefficient RT requires MPI'
#endif
  end subroutine

  subroutine write_dg_overlapping_wannier_rt_restart(comm,prefix,coefficients,state,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    complex(real64),intent(in)::coefficients(:,:)
    type(s_dg_overlapping_wannier_rt_state),intent(in)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    character(32),parameter::magic='SALMON_OW_COEFFICIENT_RT_V1'
    integer::rank,ierr,unit,ios,close_ios
    integer(int64)::stored_size,stored_digest,file_size
    character(512)::temporary
    call MPI_Comm_rank(comm,rank,ierr);ok=.false.;message='';ios=0
    if(size(coefficients,1)/=state%nwann.or.size(coefficients,2)<1.or.state%step<0.or.&
       .not.ieee_is_finite(state%time).or.state%time<0d0.or..not.finite_complex_2d(coefficients))ios=1
    call MPI_Allreduce(MPI_IN_PLACE,ios,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ios/=0)then;message='invalid overlapping-Wannier RT restart payload';return;endif
    stored_size=rt_restart_size(coefficients)
    stored_digest=rt_restart_digest(coefficients,state)
    if(rank==0)then
      temporary=trim(prefix)//'.temporary'
      open(newunit=unit,file=trim(temporary),status='replace',access='stream',form='unformatted',&
        action='write',iostat=ios)
      if(ios==0)then
        write(unit,iostat=ios)magic,stored_size,stored_digest,state%step,state%time,state%nwann,&
          size(coefficients,2),state%basis_generation,state%geometry_generation,state%basis_fingerprint,&
          state%operator_fingerprint,state%hamiltonian_fingerprint,state%observable_fingerprint,&
          state%fingerprint,state%field_coupling_convention,coefficients
        if(ios==0)flush(unit,iostat=ios)
        close_ios=0;close(unit,iostat=close_ios)
        if(ios==0)ios=close_ios
      endif
      if(ios==0)then
        inquire(file=trim(temporary),size=file_size)
        if(file_size/=stored_size)ios=1
      endif
      if(ios==0)call rename(trim(temporary),trim(prefix),ios)
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0)then;message='cannot atomically publish overlapping-Wannier RT restart';return;endif
    ok=.true.
#else
    ok=.false.;message='overlapping-Wannier RT restart requires MPI'
#endif
  end subroutine

  subroutine read_dg_overlapping_wannier_rt_restart(comm,prefix,coefficients,state,ok,message)
    integer,intent(in)::comm
    character(*),intent(in)::prefix
    complex(real64),intent(out)::coefficients(:,:)
    type(s_dg_overlapping_wannier_rt_state),intent(inout)::state
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    character(32),parameter::expected_magic='SALMON_OW_COEFFICIENT_RT_V1'
    character(32)::magic
    character(64)::convention
    integer::rank,ierr,unit,ios,step,nwann,nstate,basis_generation,geometry_generation
    integer(int64)::stored_size,stored_digest,file_size,basis_fingerprint,operator_fingerprint,&
      hamiltonian_fingerprint,observable_fingerprint,fingerprint
    real(real64)::time
    call MPI_Comm_rank(comm,rank,ierr);ok=.false.;message='';ios=0
    if(rank==0)then
      open(newunit=unit,file=trim(prefix),status='old',access='stream',form='unformatted',&
        action='read',iostat=ios)
      if(ios==0)then
        inquire(unit=unit,size=file_size)
        read(unit,iostat=ios)magic,stored_size,stored_digest,step,time,nwann,nstate,basis_generation,&
          geometry_generation,basis_fingerprint,operator_fingerprint,hamiltonian_fingerprint,&
          observable_fingerprint,fingerprint,convention,coefficients
        close(unit)
      endif
      if(ios==0)then
        if(magic/=expected_magic.or.file_size/=stored_size.or.stored_size/=rt_restart_size(coefficients).or.&
           nwann/=state%nwann.or.nstate/=size(coefficients,2).or.basis_generation/=state%basis_generation.or.&
           geometry_generation/=state%geometry_generation.or.basis_fingerprint/=state%basis_fingerprint.or.&
           operator_fingerprint/=state%operator_fingerprint.or.&
           hamiltonian_fingerprint/=state%hamiltonian_fingerprint.or.&
           observable_fingerprint/=state%observable_fingerprint.or.fingerprint/=state%fingerprint.or.&
           convention/=state%field_coupling_convention.or.step<0.or.time<0d0.or.&
           .not.ieee_is_finite(time))ios=1
      endif
      if(ios==0)then
        state%step=step;state%time=time
        if(stored_digest/=rt_restart_digest(coefficients,state))ios=1
      endif
    endif
    call MPI_Bcast(ios,1,MPI_INTEGER,0,comm,ierr)
    if(ios/=0)then;message='overlapping-Wannier RT restart is missing, corrupt, or stale';return;endif
    call MPI_Bcast(state%step,1,MPI_INTEGER,0,comm,ierr)
    call MPI_Bcast(state%time,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    call MPI_Bcast(coefficients,size(coefficients),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    ok=.true.
#else
    ok=.false.;message='overlapping-Wannier RT restart requires MPI'
#endif
  end subroutine

  integer(int64) function rt_restart_size(coefficients)
    complex(real64),intent(in)::coefficients(:,:)
    rt_restart_size=32_int64+64_int64+5_int64*int(storage_size(0)/8,int64)+&
      7_int64*int(storage_size(0_int64)/8,int64)+int(storage_size(0d0)/8,int64)+&
      int(size(coefficients),int64)*int(storage_size((0d0,0d0))/8,int64)
  end function

  integer(int64) function rt_restart_digest(coefficients,state)
    complex(real64),intent(in)::coefficients(:,:)
    type(s_dg_overlapping_wannier_rt_state),intent(in)::state
    integer::i,j
    integer(int64)::hash,bits
    hash=ieor(state%fingerprint,int(state%step,int64))
    bits=transfer(state%time,bits);hash=ieor(hash,ishftc(bits,7))
    do j=1,size(coefficients,2);do i=1,size(coefficients,1)
      bits=transfer(real(coefficients(i,j),real64),bits);hash=ieor(hash,ishftc(bits,modulo(i+3*j,63)))
      bits=transfer(aimag(coefficients(i,j)),bits);hash=ieor(hash,ishftc(bits,modulo(5*i+j,63)))
    enddo;enddo
    rt_restart_digest=hash
  end function

  subroutine generalized_hermitian_eigensystem(hamiltonian,metric,eigenvalues,vectors,ok,message)
    complex(real64),intent(inout)::hamiltonian(:,:),metric(:,:)
    real(real64),allocatable,intent(out)::eigenvalues(:)
    complex(real64),allocatable,intent(out)::vectors(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::work(:)
    real(real64),allocatable::rwork(:)
    complex(real64)::query(1)
    integer::n,lwork,info
    logical::invalid_halting,zero_halting
    external::zhegv
    n=size(hamiltonian,1);ok=.false.;message=''
    if(n<1.or.any(shape(hamiltonian)/=[n,n]).or.any(shape(metric)/=[n,n]))then
      message='invalid generalized eigensystem extent';return
    endif
    allocate(eigenvalues(n),rwork(max(1,3*n-2)))
    lwork=-1
    call ieee_get_halting_mode(ieee_invalid,invalid_halting)
    call ieee_get_halting_mode(ieee_divide_by_zero,zero_halting)
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call zhegv(1,'V','U',n,hamiltonian,n,metric,n,eigenvalues,query,lwork,rwork,info)
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_invalid,invalid_halting)
    call ieee_set_halting_mode(ieee_divide_by_zero,zero_halting)
    if(info/=0)then;message='generalized eigensystem workspace query failed';return;endif
    lwork=max(1,int(real(query(1),real64)));allocate(work(lwork))
    call ieee_set_halting_mode(ieee_invalid,.false.)
    call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call zhegv(1,'V','U',n,hamiltonian,n,metric,n,eigenvalues,work,lwork,rwork,info)
    call ieee_set_flag(ieee_invalid,.false.)
    call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_invalid,invalid_halting)
    call ieee_set_halting_mode(ieee_divide_by_zero,zero_halting)
    if(info/=0)then;message='generalized Hermitian eigensystem failed';return;endif
    if(any(.not.ieee_is_finite(eigenvalues)))then
      message='nonfinite generalized eigenvalues';return
    endif
    allocate(vectors(n,n));vectors=hamiltonian;ok=.true.
  end subroutine

  logical function finite_complex_2d(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex_2d=all(ieee_is_finite(real(values,real64))).and.&
      all(ieee_is_finite(aimag(values)))
  end function

  logical function finite_complex_3d(values)
    complex(real64),intent(in)::values(:,:,:)
    finite_complex_3d=all(ieee_is_finite(real(values,real64))).and.&
      all(ieee_is_finite(aimag(values)))
  end function
end module
