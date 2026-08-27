#include "config.h"
module dg_hybrid_continuation_residuals
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_residuals
    real(real64)::r_h=huge(1d0),r_rho=huge(1d0),r_t=huge(1d0),r_s=huge(1d0)
  end type s_dg_hybrid_residuals
  type,public::s_dg_hybrid_metric_receipt
    logical::valid=.false.
    integer::global_count=0,numerical_rank=0
    integer(int64)::fingerprint=0_int64
  end type s_dg_hybrid_metric_receipt
  public::validate_dg_hybrid_occupied_rows,build_dg_hybrid_interface_observables,&
    evaluate_dg_hybrid_residuals,evaluate_dg_hybrid_projector_change,dg_hybrid_electron_count,&
    validate_cluster_occupations
contains
  subroutine validate_dg_hybrid_occupied_rows(comm,row_ids,coefficients,s_coefficients,occupations,receipt,ok,message)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::coefficients(:,:),s_coefficients(:,:)
    real(real64),intent(in)::occupations(:)
    type(s_dg_hybrid_metric_receipt),intent(in)::receipt
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::local_overlap(:,:),overlap(:,:)
    integer::i,m,local_bad,global_bad,ierr,local_rows,global_rows,min_count,max_count,min_rank,max_rank
    integer(int64)::min_hash,max_hash
    real(real64)::defect
    ok=.false.;message='';m=size(coefficients,2);local_bad=0
    if(.not.receipt%valid.or.receipt%global_count<1.or.receipt%numerical_rank/=receipt%global_count.or.&
        receipt%fingerprint==0_int64.or.size(row_ids)/=size(coefficients,1).or.&
        any(shape(s_coefficients)/=shape(coefficients)).or.m<1.or.size(occupations)/=m.or.&
        any(row_ids<1_int64).or.any(row_ids>int(receipt%global_count,int64)).or.&
        .not.finite_complex(coefficients).or..not.finite_complex(s_coefficients).or.&
        .not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))local_bad=1
    do i=1,size(row_ids);if(count(row_ids==row_ids(i))/=1)local_bad=1;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid distributed occupied-row algebra or metric receipt';return;endif
    local_rows=size(row_ids);call MPI_Allreduce(local_rows,global_rows,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(receipt%global_count,min_count,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(receipt%global_count,max_count,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(receipt%numerical_rank,min_rank,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(receipt%numerical_rank,max_rank,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(receipt%fingerprint,min_hash,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(receipt%fingerprint,max_hash,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_rows/=receipt%global_count.or.min_count/=max_count.or.min_rank/=max_rank.or.&
        min_hash/=max_hash)then;message='rank-disagreeing frozen metric or row distribution';return;endif
    allocate(local_overlap(m,m),overlap(m,m));local_overlap=matmul(conjg(transpose(coefficients)),s_coefficients)
    call MPI_Allreduce(local_overlap,overlap,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='occupied overlap reduction failed';return;endif
    do i=1,m;overlap(i,i)=overlap(i,i)-1d0;enddo
    defect=frobenius(overlap)
    if(defect>1d-10)then;message='occupied rows are not globally S-orthonormal';return;endif
    ok=.true.
#else
    ok=.false.;message='distributed occupied algebra requires MPI'
#endif
  end subroutine validate_dg_hybrid_occupied_rows

  subroutine evaluate_dg_hybrid_projector_change(comm,previous_coefficients,current_s_coefficients,receipt,residual,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::previous_coefficients(:,:),current_s_coefficients(:,:)
    type(s_dg_hybrid_metric_receipt),intent(in)::receipt
    real(real64),intent(out)::residual
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::local_overlap(:,:),overlap(:,:)
    integer::m,ierr,local_bad,global_bad
    real(real64)::captured,lost
    ok=.false.;message='';residual=huge(1d0);m=size(previous_coefficients,2)
    local_bad=merge(0,1,receipt%valid.and.receipt%fingerprint/=0_int64.and.m>0.and.&
      size(current_s_coefficients,1)==size(previous_coefficients,1).and.size(current_s_coefficients,2)==m.and.&
      finite_complex(previous_coefficients).and.finite_complex(current_s_coefficients))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid distributed projector comparison';return;endif
    allocate(local_overlap(m,m),overlap(m,m))
    local_overlap=matmul(conjg(transpose(previous_coefficients)),current_s_coefficients)
    call MPI_Allreduce(local_overlap,overlap,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projector overlap reduction failed';return;endif
    captured=sum(abs(overlap)**2);lost=max(0d0,real(m,real64)-captured)
    if(lost<=1d-12*real(m,real64))lost=0d0
    residual=sqrt(lost/real(m,real64));ok=.true.
#else
    ok=.false.;message='distributed projector comparison requires MPI';residual=huge(1d0)
#endif
  end subroutine evaluate_dg_hybrid_projector_change

  real(real64) function dg_hybrid_electron_count(occupations) result(count)
    real(real64),intent(in)::occupations(:)
    if(size(occupations)<1.or..not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))then
      count=huge(1d0)
    else;count=sum(occupations)
    endif
  end function dg_hybrid_electron_count

  subroutine validate_cluster_occupations(occupations,cluster_ids,ok,message)
    real(real64),intent(in)::occupations(:)
    integer,intent(in)::cluster_ids(:)
    logical,intent(out)::ok;character(*),intent(out)::message
    integer::i,j
    ok=.false.;message=''
    if(size(occupations)<1.or.size(cluster_ids)/=size(occupations).or.any(cluster_ids<1).or.&
        .not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))then
      message='invalid occupation cluster';return
    endif
    do i=1,size(occupations);do j=i+1,size(occupations)
      if(cluster_ids(i)==cluster_ids(j).and.abs(occupations(i)-occupations(j))>1d-12)then
        message='symmetry-incompatible occupations in a degenerate cluster';return
      endif
    enddo;enddo
    ok=.true.
  end subroutine validate_cluster_occupations

  subroutine build_dg_hybrid_interface_observables(value_coefficients,normal_coefficients,occupations,&
      value_density,normal_density,cross_density,ok,message)
    complex(real64),intent(in)::value_coefficients(:,:),normal_coefficients(:,:)
    real(real64),intent(in)::occupations(:)
    complex(real64),intent(out)::value_density(:,:),normal_density(:,:),cross_density(:,:)
    logical,intent(out)::ok;character(*),intent(out)::message
    complex(real64),allocatable::weighted_value(:,:),weighted_normal(:,:)
    integer::i,ntrace,nocc
    ok=.false.;message='';ntrace=size(value_coefficients,1);nocc=size(value_coefficients,2)
    if(ntrace<1.or.nocc<1.or.any(shape(normal_coefficients)/=shape(value_coefficients)).or.size(occupations)/=nocc.or.&
        any(shape(value_density)/=[ntrace,ntrace]).or.any(shape(normal_density)/=[ntrace,ntrace]).or.&
        any(shape(cross_density)/=[ntrace,ntrace]).or..not.finite_complex(value_coefficients).or.&
        .not.finite_complex(normal_coefficients).or..not.all(ieee_is_finite(occupations)).or.any(occupations<0d0))then
      message='invalid occupied interface coefficient traces';return
    endif
    allocate(weighted_value(ntrace,nocc),weighted_normal(ntrace,nocc));weighted_value=value_coefficients
    weighted_normal=normal_coefficients
    do i=1,nocc
      weighted_value(:,i)=occupations(i)*weighted_value(:,i)
      weighted_normal(:,i)=occupations(i)*weighted_normal(:,i)
    enddo
    value_density=matmul(weighted_value,conjg(transpose(value_coefficients)))
    normal_density=matmul(weighted_normal,conjg(transpose(normal_coefficients)))
    cross_density=matmul(weighted_value,conjg(transpose(normal_coefficients)));ok=.true.
  end subroutine build_dg_hybrid_interface_observables

  subroutine evaluate_dg_hybrid_residuals(comm,hc,sc_epsilon,coefficients,s_coefficients,rho_output,rho_input,&
      trace_output,trace_input,residuals,ok,message)
    integer,intent(in)::comm
    complex(real64),intent(in)::hc(:,:),sc_epsilon(:,:),coefficients(:,:),s_coefficients(:,:),trace_output(:,:),trace_input(:,:)
    real(real64),intent(in)::rho_output(:),rho_input(:)
    type(s_dg_hybrid_residuals),intent(out)::residuals
    logical,intent(out)::ok;character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::local_overlap(:,:),overlap(:,:)
    real(real64)::local_sums(6),global_sums(6),local_trace(2),global_trace(2)
    integer::i,m,ierr,local_bad,global_bad
    ok=.false.;message='';m=size(coefficients,2)
    local_bad=merge(0,1,all(shape(hc)==shape(sc_epsilon)).and.all(shape(hc)==shape(coefficients)).and.&
      all(shape(s_coefficients)==shape(coefficients)).and.size(rho_output)==size(rho_input).and.&
      all(shape(trace_output)==shape(trace_input)).and.finite_complex(hc).and.finite_complex(sc_epsilon).and.&
      finite_complex(coefficients).and.finite_complex(s_coefficients).and.finite_complex(trace_output).and.&
      finite_complex(trace_input).and.all(ieee_is_finite(rho_output)).and.all(ieee_is_finite(rho_input)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid distributed continuation residual input';return;endif
    local_sums=[sum(abs(hc-sc_epsilon)**2),sum(abs(hc)**2),sum(abs(sc_epsilon)**2),&
      sum((rho_output-rho_input)**2),sum(rho_input**2),0d0]
    call MPI_Allreduce(local_sums,global_sums,6,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    local_trace=[sum(abs(trace_output-trace_input)**2),sum(abs(trace_input)**2)]
    call MPI_Allreduce(local_trace,global_trace,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    residuals%r_h=sqrt(global_sums(1))/max(1d0,sqrt(global_sums(2)),sqrt(global_sums(3)))
    residuals%r_rho=sqrt(global_sums(4))/max(1d0,sqrt(global_sums(5)))
    residuals%r_t=sqrt(global_trace(1))/max(1d0,sqrt(global_trace(2)))
    allocate(local_overlap(m,m),overlap(m,m));local_overlap=matmul(conjg(transpose(coefficients)),s_coefficients)
    call MPI_Allreduce(local_overlap,overlap,m*m,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do i=1,m;overlap(i,i)=overlap(i,i)-1d0;enddo
    residuals%r_s=frobenius(overlap);ok=ierr==MPI_SUCCESS
    if(.not.ok)message='distributed residual reduction failed'
#else
    ok=.false.;message='distributed continuation residuals require MPI'
#endif
  end subroutine evaluate_dg_hybrid_residuals

  real(real64) function frobenius(values) result(norm)
    complex(real64),intent(in)::values(:,:);norm=sqrt(sum(abs(values)**2))
  end function frobenius
  logical function finite_complex(values) result(ok)
    complex(real64),intent(in)::values(:,:);ok=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex
end module dg_hybrid_continuation_residuals
