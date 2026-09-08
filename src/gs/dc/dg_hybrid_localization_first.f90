#include "config.h"
module dg_hybrid_localization_first
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_overlapping_wannier_construction,only:orthonormalize_dg_distributed_seed_space
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private

  type,public :: s_dg_hybrid_localization_receipt
    logical :: valid=.false.,symmetry_constrained=.false.,converged=.false.
    integer :: raw_rank=0,retained_rank=0,iterations=0
    real(real64) :: spread_min=0d0,spread_max=0d0
    real(real64) :: spread_mean=0d0,spread_total=0d0
    real(real64) :: transform_unitarity_defect=huge(0d0)
    integer(int64) :: seed_fingerprint=0_int64,transform_fingerprint=0_int64
  end type s_dg_hybrid_localization_receipt

  public::prepare_dg_hybrid_localization_first_seed
  public::build_dg_hybrid_localization_receipt
contains
  subroutine prepare_dg_hybrid_localization_first_seed(comm,point_ids,raw_seed_values,&
      weights,tolerance,prepared_seed,retained_rank,seed_fingerprint,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::point_ids(:)
    complex(real64),intent(in)::raw_seed_values(:,:)
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),allocatable,intent(out)::prepared_seed(:,:)
    integer,intent(out)::retained_rank
    integer(int64),intent(out)::seed_fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nlocal,nraw,ierr,local_bad,global_bad,minimum_integer,maximum_integer
    integer(int64)::tolerance_bits,minimum_bits,maximum_bits
    logical::fingerprint_ok

    ok=.false.;message='';retained_rank=0;seed_fingerprint=0_int64
    nlocal=size(raw_seed_values,2);nraw=size(raw_seed_values,1)
    local_bad=merge(0,1,nraw>0.and.nlocal>0.and.size(point_ids)==nlocal.and.&
      size(weights)==nlocal)
    if(local_bad==0)then
      if(.not.ieee_is_finite(tolerance))then
        local_bad=1
      else if(tolerance<=0d0.or.tolerance>=1d0)then
        local_bad=1
      endif
    endif
    if(local_bad==0)then
      if(any(point_ids<=0_int64).or..not.all(ieee_is_finite(weights)).or.&
          .not.all(ieee_is_finite(real(raw_seed_values))).or.&
          .not.all(ieee_is_finite(aimag(raw_seed_values))))then
        local_bad=1
      else if(any(weights<=0d0))then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid localization-first raw seed contract';return
    endif
    call agree_integer(comm,nraw,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='localization-first raw rank disagrees across ranks';return
    endif
    call agree_integer(comm,nlocal,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='localization-first local point count disagrees across ranks';return
    endif
    tolerance_bits=transfer(tolerance,tolerance_bits)
    call agree_int64(comm,tolerance_bits,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='localization-first metric tolerance disagrees across ranks';return
    endif

    call orthonormalize_dg_distributed_seed_space(comm,raw_seed_values,weights,tolerance,&
      prepared_seed,retained_rank,ok,message)
    local_bad=0
    if(.not.ok.or.retained_rank/=nraw)then
      local_bad=1
    else if(.not.allocated(prepared_seed))then
      local_bad=1
    else if(any(shape(prepared_seed)/=[nraw,nlocal]))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(prepared_seed))deallocate(prepared_seed)
      retained_rank=0;ok=.false.
      message='localization-first metric preparation lost raw seed rank';return
    endif
    call fingerprint_distributed_seed(comm,point_ids,raw_seed_values,weights,&
      seed_fingerprint,fingerprint_ok)
    if(.not.fingerprint_ok.or.seed_fingerprint==0_int64)then
      if(allocated(prepared_seed))deallocate(prepared_seed)
      retained_rank=0;seed_fingerprint=0_int64;ok=.false.
      message='localization-first seed fingerprint failed';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='localization-first seed preparation requires MPI'
    retained_rank=0;seed_fingerprint=0_int64
#endif
  end subroutine prepare_dg_hybrid_localization_first_seed

  subroutine build_dg_hybrid_localization_receipt(comm,raw_rank,retained_rank,&
      seed_fingerprint,transform,centers,spreads,converged,iterations,tolerance,&
      receipt,ok,message)
    integer,intent(in)::comm,raw_rank,retained_rank,iterations
    integer(int64),intent(in)::seed_fingerprint
    complex(real64),intent(in)::transform(:,:)
    real(real64),intent(in)::centers(:,:),spreads(:),tolerance
    logical,intent(in)::converged
    type(s_dg_hybrid_localization_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::gram(:,:)
    real(real64)::unitarity_defect,total,mean_value,minimum_value,maximum_value
    integer::i,ierr,local_bad,global_bad,minimum_integer,maximum_integer
    integer(int64)::minimum_hash,maximum_hash,tolerance_bits,transform_hash,payload_hash

    receipt=s_dg_hybrid_localization_receipt();ok=.false.;message=''
    local_bad=merge(0,1,raw_rank>0.and.retained_rank==raw_rank.and.converged.and.&
      iterations>=0.and.seed_fingerprint/=0_int64)
    if(local_bad==0)then
      if(.not.ieee_is_finite(tolerance))then
        local_bad=1
      else if(tolerance<=0d0.or.tolerance>=1d0)then
        local_bad=1
      endif
    endif
    if(local_bad==0)then
      if(any(shape(transform)/=[raw_rank,retained_rank]).or.&
          any(shape(centers)/=[3,retained_rank]).or.size(spreads)/=retained_rank)then
        local_bad=1
      else if(.not.all(ieee_is_finite(real(transform))).or.&
          .not.all(ieee_is_finite(aimag(transform))).or.&
          .not.all(ieee_is_finite(centers)).or..not.all(ieee_is_finite(spreads)))then
        local_bad=1
      else if(any(spreads<0d0))then
        local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid fixed-rank localization result';return
    endif
    call agree_integer(comm,raw_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='localization raw rank disagrees across ranks';return
    endif
    call agree_integer(comm,retained_rank,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='localization retained rank disagrees across ranks';return
    endif
    call agree_integer(comm,iterations,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='localization iteration receipt disagrees across ranks';return
    endif
    call agree_int64(comm,seed_fingerprint,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='localization seed fingerprint disagrees across ranks';return
    endif
    tolerance_bits=transfer(tolerance,tolerance_bits)
    call agree_int64(comm,tolerance_bits,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='localization receipt tolerance disagrees across ranks';return
    endif

    allocate(gram(retained_rank,retained_rank))
    gram=matmul(conjg(transpose(transform)),transform)
    do i=1,retained_rank;gram(i,i)=gram(i,i)-cmplx(1d0,0d0,real64);enddo
    unitarity_defect=maxval(abs(gram))
    total=sum(spreads);mean_value=total/real(retained_rank,real64)
    minimum_value=minval(spreads);maximum_value=maxval(spreads)
    local_bad=merge(0,1,ieee_is_finite(unitarity_defect).and.&
      unitarity_defect<=tolerance*real(max(1,retained_rank),real64).and.&
      ieee_is_finite(total).and.ieee_is_finite(mean_value))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid localization unitarity or spread aggregate';return
    endif

    transform_hash=fingerprint_transform(transform)
    payload_hash=fingerprint_localization_payload(centers,spreads)
    call agree_int64(comm,transform_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='localization transform disagrees across ranks';return
    endif
    call agree_int64(comm,payload_hash,minimum_hash,maximum_hash,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_hash/=maximum_hash)then
      message='localization center or spread payload disagrees across ranks';return
    endif

    receipt%symmetry_constrained=.false.;receipt%converged=.true.
    receipt%raw_rank=raw_rank;receipt%retained_rank=retained_rank
    receipt%iterations=iterations;receipt%spread_min=minimum_value
    receipt%spread_max=maximum_value;receipt%spread_mean=mean_value
    receipt%spread_total=total;receipt%transform_unitarity_defect=unitarity_defect
    receipt%seed_fingerprint=seed_fingerprint
    receipt%transform_fingerprint=transform_hash
    receipt%valid=.true.;ok=.true.
#else
    receipt=s_dg_hybrid_localization_receipt()
    ok=.false.;message='localization-first receipt requires MPI'
#endif
  end subroutine build_dg_hybrid_localization_receipt

#ifdef USE_MPI
  subroutine fingerprint_distributed_seed(comm,point_ids,seed_values,weights,&
      fingerprint,ok)
    integer,intent(in)::comm
    integer(int64),intent(in)::point_ids(:)
    complex(real64),intent(in)::seed_values(:,:)
    real(real64),intent(in)::weights(:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    integer::i,p,ierr,local_count,global_count
    integer(int64)::bits,entry,local_hash,global_hash
    local_hash=0_int64
    do p=1,size(point_ids)
      entry=ieor(int(z'243F6A8885A308D3',int64),ishftc(point_ids(p),7))
      bits=transfer(weights(p),bits);entry=ieor(ishftc(entry,11),bits)
      do i=1,size(seed_values,1)
        entry=ieor(ishftc(entry,9),int(i,int64))
        bits=transfer(real(seed_values(i,p),real64),bits)
        entry=ieor(ishftc(entry,13),bits)
        bits=transfer(aimag(seed_values(i,p)),bits)
        entry=ieor(ishftc(entry,17),bits)
      enddo
      local_hash=ieor(local_hash,entry)
    enddo
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_int64;ok=.false.;return;endif
    local_count=size(point_ids)
    call MPI_Allreduce(local_count,global_count,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;fingerprint=0_int64;ok=.false.;return;endif
    fingerprint=ieor(int(z'13198A2E03707344',int64),global_hash)
    fingerprint=ieor(ishftc(fingerprint,7),int(size(seed_values,1),int64))
    fingerprint=ieor(ishftc(fingerprint,11),int(global_count,int64))
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
  end subroutine fingerprint_distributed_seed

  integer(int64) function fingerprint_transform(transform) result(fingerprint)
    complex(real64),intent(in)::transform(:,:)
    integer::i,j
    integer(int64)::bits
    fingerprint=int(z'A4093822299F31D0',int64)
    fingerprint=ieor(ishftc(fingerprint,7),int(size(transform,1),int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(size(transform,2),int64))
    do j=1,size(transform,2)
      do i=1,size(transform,1)
        bits=transfer(real(transform(i,j),real64),bits)
        fingerprint=ieor(ishftc(fingerprint,11),bits)
        bits=transfer(aimag(transform(i,j)),bits)
        fingerprint=ieor(ishftc(fingerprint,13),bits)
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
  end function fingerprint_transform

  integer(int64) function fingerprint_localization_payload(centers,spreads) result(fingerprint)
    real(real64),intent(in)::centers(:,:),spreads(:)
    integer::i,j
    integer(int64)::bits
    fingerprint=int(z'082EFA98EC4E6C89',int64)
    do j=1,size(centers,2)
      do i=1,size(centers,1)
        bits=transfer(centers(i,j),bits)
        fingerprint=ieor(ishftc(fingerprint,9),bits)
      enddo
      bits=transfer(spreads(j),bits)
      fingerprint=ieor(ishftc(fingerprint,15),bits)
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
  end function fingerprint_localization_payload

  subroutine agree_integer(comm,value,minimum,maximum,ierr)
    integer,intent(in)::comm,value
    integer,intent(out)::minimum,maximum,ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer

  subroutine agree_int64(comm,value,minimum,maximum,ierr)
    integer,intent(in)::comm
    integer(int64),intent(in)::value
    integer(int64),intent(out)::minimum,maximum
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module dg_hybrid_localization_first
