#include "config.h"
module dg_fragment_scdm_gauge
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_dg_fragment_scdm_gauge
contains
  subroutine build_dg_fragment_scdm_gauge(comm,fragment_id,basis_generation,grid_ids,values,weights,&
      tolerance,byte_limit,selected_grid_ids,a_matrix,unitarity_defect,projector_defect,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,fragment_id,basis_generation
    integer(int64),intent(in)::grid_ids(:),byte_limit
    complex(real64),intent(in)::values(:,:)
    real(real64),intent(in)::weights(:),tolerance
    integer(int64),allocatable,intent(out)::selected_grid_ids(:)
    complex(real64),allocatable,intent(out)::a_matrix(:,:)
    real(real64),intent(out)::unitarity_defect,projector_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,status,global_status,allocation_status,nband,nlocal,global_npoint
    integer::minimum_integer,maximum_integer,p,i,j,k,owner_token,owner_rank,svd_info,svd_lwork
    integer,allocatable::grid_counts(:)
    integer(int64)::minimum_int64,maximum_int64,local_npoint64,global_npoint64,&
      local_required,root_extra,global_required,&
      complex_bytes,real_bytes,integer_bytes,logical_bytes,bits,hash,local_candidate,selected_id
    real(real64)::minimum_real,maximum_real,local_maximum,global_maximum,tie_tolerance,&
      gram_defect,singular_scale
    real(real64),allocatable::residual_norms(:),singular_values(:),svd_rwork(:)
    complex(real64),allocatable::residual(:,:),selected_columns(:,:),gram(:,:),pivot_vector(:),&
      original_vector(:),q(:),svd_input(:,:),svd_u(:,:),svd_vt(:,:),svd_work(:)
    logical,allocatable::selected_local(:)
    logical::arithmetic_ok
    interface
      subroutine zgesvd(jobu,jobvt,m,n,a,lda,s,u,ldu,vt,ldvt,work,lwork,rwork,info)
        character,intent(in)::jobu,jobvt
        integer,intent(in)::m,n,lda,ldu,ldvt,lwork
        complex(8),intent(inout)::a(lda,*),work(*)
        real(8),intent(out)::s(*),rwork(*)
        complex(8),intent(out)::u(ldu,*),vt(ldvt,*)
        integer,intent(out)::info
      end subroutine zgesvd
    end interface

    ok=.false.;message='';unitarity_defect=huge(1d0);projector_defect=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;status=0
    call MPI_Comm_rank(comm,rank,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SCDM communicator rank query failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr)
    if(ierr/=MPI_SUCCESS)then;message='SCDM communicator size query failed';return;endif
    nband=size(values,1);nlocal=size(values,2)
    if(fragment_id<1.or.basis_generation<1.or.nband<1.or.nband>huge(0)/8.or.nlocal<0.or.&
        size(grid_ids)/=nlocal.or.size(weights)/=nlocal.or.byte_limit<=0_int64.or.&
        .not.ieee_is_finite(tolerance).or.tolerance<1d-15.or.tolerance>1d-2.or.&
        .not.all(ieee_is_finite(real(values))).or..not.all(ieee_is_finite(aimag(values))).or.&
        .not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or.any(grid_ids<1_int64))status=1
    if(nband>0)then
      if(nband>huge(0)/nband)status=1
    endif
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='invalid distributed fragment SCDM gauge contract';return
    endif
    call agree_integer(fragment_id,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)status=1
    call agree_integer(basis_generation,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)status=1
    call agree_integer(nband,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)status=1
    call MPI_Allreduce(byte_limit,minimum_int64,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(byte_limit,maximum_int64,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_int64/=maximum_int64)status=1
    call MPI_Allreduce(tolerance,minimum_real,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(tolerance,maximum_real,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_real/=maximum_real)status=1
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM controls or identity disagree across ranks';return
    endif
    local_npoint64=int(nlocal,int64)
    call MPI_Allreduce(local_npoint64,global_npoint64,1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    status=0
    if(ierr/=MPI_SUCCESS.or.global_npoint64<int(nband,int64).or.&
        global_npoint64>int(huge(0),int64))status=1
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM retained rank exceeds the distributed grid';return
    endif
    global_npoint=int(global_npoint64)

    complex_bytes=int(storage_size((0d0,0d0))/8,int64)
    real_bytes=int(storage_size(0d0)/8,int64)
    integer_bytes=int(storage_size(0)/8,int64)
    logical_bytes=int(storage_size(.false.)/8,int64)
    arithmetic_ok=.true.;local_required=0_int64
    call add_product(local_required,[int(nband,int64),int(nlocal,int64),complex_bytes],arithmetic_ok)
    call add_product(local_required,[2_int64,int(nband,int64),int(nband,int64),complex_bytes],arithmetic_ok)
    call add_product(local_required,[3_int64,int(nband,int64),complex_bytes],arithmetic_ok)
    call add_product(local_required,[int(nlocal,int64),real_bytes],arithmetic_ok)
    call add_product(local_required,[int(nlocal,int64),logical_bytes],arithmetic_ok)
    call add_product(local_required,[int(global_npoint,int64),integer_bytes],arithmetic_ok)
    call add_product(local_required,[int(nband,int64),8_int64],arithmetic_ok)
    root_extra=0_int64
    call add_product(root_extra,[5_int64,int(nband,int64),int(nband,int64),complex_bytes],arithmetic_ok)
    call add_product(root_extra,[8_int64,int(nband,int64),complex_bytes],arithmetic_ok)
    call add_product(root_extra,[6_int64,int(nband,int64),real_bytes],arithmetic_ok)
    if(rank==0)call checked_add(local_required,root_extra,arithmetic_ok)
    call MPI_Allreduce(local_required,global_required,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    status=merge(0,1,ierr==MPI_SUCCESS.and.arithmetic_ok.and.global_required<=byte_limit)
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM byte limit or workspace estimate rejected';return
    endif
    workspace_peak_bytes=global_required

    allocate(selected_grid_ids(nband),residual(nband,nlocal),residual_norms(nlocal),&
      selected_local(nlocal),selected_columns(nband,nband),gram(nband,nband),&
      pivot_vector(nband),original_vector(nband),q(nband),grid_counts(global_npoint),stat=allocation_status)
    call MPI_Allreduce(allocation_status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM distributed workspace allocation failed';return
    endif
    grid_counts=0
    do p=1,nlocal
      if(grid_ids(p)>int(global_npoint,int64))then
        status=1
      else
        grid_counts(int(grid_ids(p)))=grid_counts(int(grid_ids(p)))+1
      endif
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,grid_counts,global_npoint,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(grid_counts/=1))status=1
    call MPI_Allreduce(status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM grid IDs are not a unique complete layout';return
    endif

    gram=(0d0,0d0)
    do p=1,nlocal;do j=1,nband;do i=1,nband
      gram(i,j)=gram(i,j)+weights(p)*conjg(values(i,p))*values(j,p)
    enddo;enddo;enddo
    call MPI_Allreduce(MPI_IN_PLACE,gram,nband*nband,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment SCDM Gram reduction failed';return;endif
    do i=1,nband;gram(i,i)=gram(i,i)-(1d0,0d0);enddo
    gram_defect=maxval(abs(gram))
    call MPI_Allreduce(MPI_IN_PLACE,gram_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(gram_defect).or.gram_defect>10d0*tolerance)then
      message='fragment SCDM retained frame is not orthonormal or full rank';return
    endif

    do p=1,nlocal
      residual(:,p)=sqrt(weights(p))*conjg(values(:,p))
    enddo
    selected_local=.false.;selected_columns=(0d0,0d0)
    do k=1,nband
      local_maximum=-1d0
      do p=1,nlocal
        if(selected_local(p))cycle
        residual_norms(p)=max(0d0,real(dot_product(residual(:,p),residual(:,p)),real64))
        local_maximum=max(local_maximum,residual_norms(p))
      enddo
      call MPI_Allreduce(local_maximum,global_maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(global_maximum).or.global_maximum<=tolerance*tolerance)then
        message='fragment SCDM pivot sequence lost numerical rank';return
      endif
      tie_tolerance=64d0*epsilon(1d0)*max(1d0,global_maximum)
      local_candidate=huge(0_int64)
      do p=1,nlocal
        if(selected_local(p))cycle
        if(residual_norms(p)>=global_maximum-tie_tolerance)&
          local_candidate=min(local_candidate,grid_ids(p))
      enddo
      call MPI_Allreduce(local_candidate,selected_id,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.selected_id==huge(0_int64))then
        message='fragment SCDM pivot selection failed';return
      endif
      selected_grid_ids(k)=selected_id;owner_token=0
      pivot_vector=(0d0,0d0);original_vector=(0d0,0d0)
      do p=1,nlocal
        if(grid_ids(p)/=selected_id)cycle
        selected_local(p)=.true.;owner_token=rank+1
        pivot_vector=residual(:,p)
        original_vector=sqrt(weights(p))*conjg(values(:,p))
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,owner_token,1,MPI_INTEGER,MPI_SUM,comm,ierr)
      owner_rank=owner_token-1
      if(ierr/=MPI_SUCCESS.or.owner_rank<0.or.owner_rank>=nproc)then
        message='fragment SCDM pivot ownership failed';return
      endif
      call MPI_Bcast(pivot_vector,nband,MPI_DOUBLE_COMPLEX,owner_rank,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(original_vector,nband,MPI_DOUBLE_COMPLEX,owner_rank,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='fragment SCDM pivot broadcast failed';return;endif
      q=pivot_vector/sqrt(global_maximum);selected_columns(:,k)=original_vector
      do p=1,nlocal
        if(selected_local(p))cycle
        residual(:,p)=residual(:,p)-q*dot_product(q,residual(:,p))
      enddo
    enddo

    allocation_status=0
    if(rank==0)then
      svd_lwork=max(1,8*nband)
      allocate(a_matrix(nband,nband),svd_input(nband,nband),svd_u(nband,nband),&
        svd_vt(nband,nband),svd_work(svd_lwork),singular_values(nband),&
        svd_rwork(max(1,5*nband)),stat=allocation_status)
    else
      allocate(a_matrix(0,0),stat=allocation_status)
    endif
    call MPI_Allreduce(allocation_status,global_status,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_status/=0)then
      message='fragment SCDM polar workspace allocation failed';return
    endif
    status=0;unitarity_defect=0d0;projector_defect=0d0
    if(rank==0)then
      svd_input=selected_columns
      call zgesvd('A','A',nband,nband,svd_input,nband,singular_values,svd_u,nband,svd_vt,nband,&
        svd_work,svd_lwork,svd_rwork,svd_info)
      if(svd_info/=0.or..not.all(ieee_is_finite(singular_values)))then
        status=1
      else
        singular_scale=max(1d0,maxval(singular_values))
        if(minval(singular_values)<=tolerance*singular_scale)status=1
      endif
      if(status==0)then
        a_matrix=matmul(svd_u,svd_vt)
        gram=matmul(conjg(transpose(a_matrix)),a_matrix)
        do i=1,nband;gram(i,i)=gram(i,i)-(1d0,0d0);enddo
        unitarity_defect=maxval(abs(gram))
        if(.not.ieee_is_finite(unitarity_defect).or.unitarity_defect>10d0*tolerance)status=2
        gram=matmul(conjg(a_matrix),transpose(a_matrix))
        do i=1,nband;gram(i,i)=gram(i,i)-(1d0,0d0);enddo
        projector_defect=max(gram_defect,maxval(abs(gram)))
        if(.not.ieee_is_finite(projector_defect).or.projector_defect>10d0*tolerance)status=2
      endif
    endif
    call MPI_Bcast(status,1,MPI_INTEGER,0,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Bcast(unitarity_defect,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Bcast(projector_defect,1,MPI_DOUBLE_PRECISION,0,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.status/=0)then
      message='fragment SCDM selected overlap has no safe full-rank unitary polar gauge';return
    endif

    hash=int(z'6A09E667F3BCC909',int64)
    if(rank==0)then
      hash=ieor(ishftc(hash,7),int(fragment_id,int64))
      hash=ieor(ishftc(hash,7),int(basis_generation,int64))
      hash=ieor(ishftc(hash,7),int(nband,int64))
      do i=1,nband
        hash=ieor(ishftc(hash,11),selected_grid_ids(i))
      enddo
      do j=1,nband;do i=1,nband
        bits=transfer(real(a_matrix(i,j),real64),bits);hash=ieor(ishftc(hash,7),bits)
        bits=transfer(aimag(a_matrix(i,j)),bits);hash=ieor(ishftc(hash,11),bits)
      enddo;enddo
      if(hash==0_int64)hash=1_int64
    endif
    call MPI_Bcast(hash,1,MPI_INTEGER8,0,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment SCDM fingerprint broadcast failed';return;endif
    fingerprint=hash;ok=.true.;message=''
#else
    ok=.false.;message='fragment SCDM gauge requires MPI';unitarity_defect=huge(1d0)
    projector_defect=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine agree_integer(value,minimum_value,maximum_value,mpi_status)
      integer,intent(in)::value
      integer,intent(out)::minimum_value,maximum_value,mpi_status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,mpi_status)
      if(mpi_status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,mpi_status)
    end subroutine agree_integer

    subroutine checked_add(total,addition,valid)
      integer(int64),intent(inout)::total
      integer(int64),intent(in)::addition
      logical,intent(inout)::valid
      if(.not.valid)return
      if(addition<0_int64.or.total>huge(total)-addition)then;valid=.false.;return;endif
      total=total+addition
    end subroutine checked_add

    subroutine add_product(total,factors,valid)
      integer(int64),intent(inout)::total
      integer(int64),intent(in)::factors(:)
      logical,intent(inout)::valid
      integer::factor
      integer(int64)::product
      if(.not.valid)return
      product=1_int64
      do factor=1,size(factors)
        if(factors(factor)<0_int64)then;valid=.false.;return;endif
        if(factors(factor)>0_int64.and.product>huge(product)/factors(factor))then
          valid=.false.;return
        endif
        product=product*factors(factor)
      enddo
      call checked_add(total,product,valid)
    end subroutine add_product
#endif
  end subroutine build_dg_fragment_scdm_gauge
end module dg_fragment_scdm_gauge
