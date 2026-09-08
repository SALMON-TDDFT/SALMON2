#include "config.h"
module dg_hybrid_block_cg
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero,ieee_overflow
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine dg_hybrid_block_apply(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::ok
    end subroutine dg_hybrid_block_apply
  end interface
  public::solve_dg_hybrid_block_cg
contains
  subroutine solve_dg_hybrid_block_cg(comm,global_count,row_ids,initial_vectors,apply_h,apply_s,&
      outer_density_residual,final_tolerance,maximum_iterations,coefficients,eigenvalues,iterations,&
      maximum_residual,stop_reason,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_count,maximum_iterations
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::initial_vectors(:,:)
    procedure(dg_hybrid_block_apply)::apply_h,apply_s
    real(real64),intent(in)::outer_density_residual,final_tolerance
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::eigenvalues(:),maximum_residual
    integer,intent(out)::iterations
    character(*),intent(out)::stop_reason
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,nowned,nstate,local_bad,global_bad,minimum_integer,maximum_integer,allocation_status
    integer::i,j,k,kdim,residual_rank,stagnant_count,info,lwork
    integer,allocatable::ownership_count(:)
    integer(int64)::bits,minimum_bits,maximum_bits,complex_elements,real_elements,integer_elements,quantized,term
    complex(real64),allocatable::c(:,:),hc(:,:),sc(:,:),residual(:,:),trial(:,:),htrial(:,:),strial(:,:),&
      projected_h(:,:),projected_s(:,:),vectors(:,:),work(:),residual_metric(:,:),overlap(:,:)
    real(real64),allocatable::small_eigenvalues(:),rwork(:)
    real(real64)::adaptive_target,previous_residual,local_value,global_value,scale,quantization_scale
    logical::callback_ok,halt_invalid,halt_zero,halt_overflow
    ok=.false.;message='';stop_reason='invalid_contract';iterations=0;maximum_residual=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;eigenvalues=0d0
    nowned=size(row_ids);nstate=size(initial_vectors,2);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='block CG communicator failed';return;endif
    call agree_integer(global_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing block CG extent';return;endif
    call agree_integer(nstate,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing block CG state count';return;endif
    call agree_integer(maximum_iterations,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing block CG cap';return;endif
    bits=transfer(outer_density_residual,bits);call agree_int64(bits,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing outer density residual';return;endif
    bits=transfer(final_tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing block CG tolerance';return;endif
    if(global_count<1.or.nstate<1.or.nstate>global_count.or.maximum_iterations<1.or.maximum_iterations>256)local_bad=1
    if(nstate>huge(0)/6)local_bad=1
    term=int(nstate,int64)*int(nstate,int64)
    if(term>int(huge(0),int64)/4_int64)local_bad=1
    if(size(initial_vectors,1)/=nowned.or.size(eigenvalues)/=nstate)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_count),int64)))local_bad=1
    if(.not.ieee_is_finite(outer_density_residual).or.outer_density_residual<0d0.or.&
      .not.ieee_is_finite(final_tolerance).or.final_tolerance<1d-15.or.final_tolerance>1d-2)local_bad=1
    if(.not.finite_matrix(initial_vectors))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid adaptive block CG contract';return;endif
    complex_elements=0_int64
    if(int(nowned,int64)>0_int64.and.int(2*nstate,int64)>huge(term)/int(nowned,int64))local_bad=1
    if(local_bad==0)then
      term=int(nowned,int64)*int(2*nstate,int64)
      if(term>huge(term)/10_int64)then;local_bad=1;else;complex_elements=10_int64*term;endif
    endif
    if(local_bad==0)then
      term=int(nstate,int64)*int(nstate,int64)
      if(term>huge(term)/12_int64.or.complex_elements>huge(complex_elements)-12_int64*term)then
        local_bad=1
      else
        complex_elements=complex_elements+12_int64*term
      endif
    endif
    real_elements=4_int64*int(nstate,int64);integer_elements=int(global_count,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.real_elements>huge(workspace_peak_bytes)/8_int64.or.&
      integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0)then
      workspace_peak_bytes=16_int64*complex_elements+8_int64*real_elements
      if(workspace_peak_bytes>huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
      if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+4_int64*integer_elements
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='adaptive block CG workspace receipt overflow';return;endif
    allocate(ownership_count(global_count),c(nowned,nstate),hc(nowned,nstate),sc(nowned,nstate),&
      residual(nowned,nstate),residual_metric(nowned,nstate),trial(nowned,2*nstate),htrial(nowned,2*nstate),&
      strial(nowned,2*nstate),projected_h(2*nstate,2*nstate),projected_s(2*nstate,2*nstate),&
      vectors(2*nstate,2*nstate),overlap(nstate,nstate),small_eigenvalues(2*nstate),&
      rwork(max(1,6*nstate-2)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate adaptive block CG workspace';return;endif
    ownership_count=0;do i=1,nowned;ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;call cleanup();message='block CG rows are not owned exactly once';return;endif
    c=initial_vectors;call orthonormalize(c,nstate,callback_ok)
    call callback_consensus(callback_ok,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='initial block CG vectors are rank deficient';return;endif
    adaptive_target=max(final_tolerance,min(1d-3,0.1d0*outer_density_residual))
    previous_residual=huge(1d0);stagnant_count=0;stop_reason='iteration_cap'
    do iterations=1,maximum_iterations
      trial(:,1:nstate)=c;kdim=nstate;call rayleigh_ritz(trial,kdim,callback_ok)
      call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='block CG Rayleigh-Ritz failed';stop_reason='solver_failure';return;endif
      residual=hc
      do j=1,nstate;residual(:,j)=residual(:,j)-eigenvalues(j)*sc(:,j);enddo
      local_bad=merge(0,1,finite_matrix(residual).and.all(ieee_is_finite(eigenvalues)))
      call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
      call ieee_set_flag(ieee_overflow,.false.)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='nonfinite block CG residual';stop_reason='solver_failure';return;endif
      local_value=0d0;if(nowned>0)local_value=maxval(abs(residual))
      call disable_halting();call MPI_Allreduce(local_value,maximum_residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      call restore_halting()
      scale=max(1d0,maxval(abs(eigenvalues)))
      if(maximum_residual<=adaptive_target*scale)then;stop_reason='adaptive_target';ok=.true.;exit;endif
      if(iterations==maximum_iterations)exit
      if(iterations>1)then
        if(maximum_residual>=0.995d0*previous_residual)then;stagnant_count=stagnant_count+1;else;stagnant_count=0;endif
        if(maximum_residual>1.5d0*previous_residual)then;stop_reason='residual_growth';exit;endif
        if(stagnant_count>=3)then;stop_reason='residual_stagnation';exit;endif
      endif
      previous_residual=maximum_residual
      residual=-residual;call orthogonalize_residual(residual,residual_rank,callback_ok)
      call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;stop_reason='residual_rank_loss';exit;endif
      trial(:,1:nstate)=c;trial(:,nstate+1:nstate+residual_rank)=residual(:,1:residual_rank);kdim=nstate+residual_rank
      call rayleigh_ritz(trial,kdim,callback_ok);call callback_consensus(callback_ok,global_bad,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;stop_reason='solver_failure';exit;endif
    enddo
    if(.not.ok)then;call cleanup();message='adaptive block CG stopped before its target';return;endif
    allocate(coefficients(nowned,nstate),stat=allocation_status);local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate block CG output';return;endif
    coefficients=c;quantization_scale=1000d0*final_tolerance;fingerprint=2017_int64
    do j=1,nstate
      if(abs(eigenvalues(j))/quantization_scale>0.25d0*real(huge(0_int64),real64))then
        call cleanup();message='block CG fingerprint range is unsafe';return
      endif
      quantized=nint(eigenvalues(j)/quantization_scale,int64);fingerprint=ieor(fingerprint,ishftc(quantized,mod(13*j,63)))
    enddo
    call cleanup(.true.);message=''
  contains
    subroutine rayleigh_ritz(basis,dimension,step_ok)
      complex(real64),intent(in)::basis(:,:);integer,intent(in)::dimension;logical,intent(out)::step_ok
      complex(real64)::query(1);integer::lapack_info
      external::zhegv
      htrial(:,1:dimension)=(0d0,0d0);strial(:,1:dimension)=(0d0,0d0)
      call apply_h(basis(:,1:dimension),htrial(:,1:dimension),callback_ok);step_ok=callback_ok
      call apply_s(basis(:,1:dimension),strial(:,1:dimension),callback_ok);step_ok=step_ok.and.callback_ok
      if(.not.step_ok)return
      projected_h=(0d0,0d0);projected_s=(0d0,0d0)
      projected_h(1:dimension,1:dimension)=matmul(conjg(transpose(basis(:,1:dimension))),htrial(:,1:dimension))
      projected_s(1:dimension,1:dimension)=matmul(conjg(transpose(basis(:,1:dimension))),strial(:,1:dimension))
      call MPI_Allreduce(MPI_IN_PLACE,projected_h,4*nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      call MPI_Allreduce(MPI_IN_PLACE,projected_s,4*nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      vectors(1:dimension,1:dimension)=projected_h(1:dimension,1:dimension);lwork=-1
      call disable_halting();call zhegv(1,'V','U',dimension,vectors,2*nstate,projected_s,2*nstate,&
        small_eigenvalues,query,lwork,rwork,lapack_info)
      if(lapack_info==0)then
        lwork=max(1,int(real(query(1))));if(allocated(work))deallocate(work);allocate(work(lwork),stat=allocation_status)
        if(allocation_status==0)then
          vectors(1:dimension,1:dimension)=projected_h(1:dimension,1:dimension)
          call zhegv(1,'V','U',dimension,vectors,2*nstate,projected_s,2*nstate,small_eigenvalues,work,lwork,rwork,lapack_info)
        else;lapack_info=-999;endif
      endif
      call restore_halting();step_ok=lapack_info==0.and.all(ieee_is_finite(small_eigenvalues(1:dimension)))
      if(.not.step_ok)return
      c=matmul(basis(:,1:dimension),vectors(1:dimension,1:nstate))
      hc=matmul(htrial(:,1:dimension),vectors(1:dimension,1:nstate))
      sc=matmul(strial(:,1:dimension),vectors(1:dimension,1:nstate));eigenvalues=small_eigenvalues(1:nstate)
    end subroutine rayleigh_ritz
    subroutine orthonormalize(block,dimension,step_ok)
      complex(real64),intent(inout)::block(:,:);integer,intent(in)::dimension;logical,intent(out)::step_ok
      complex(real64)::query(1);integer::lapack_info,column
      external::zheev
      call apply_s(block(:,1:dimension),strial(:,1:dimension),callback_ok);step_ok=callback_ok;if(.not.step_ok)return
      projected_s=(0d0,0d0)
      projected_s(1:dimension,1:dimension)=matmul(conjg(transpose(block(:,1:dimension))),strial(:,1:dimension))
      call MPI_Allreduce(MPI_IN_PLACE,projected_s,4*nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      vectors(1:dimension,1:dimension)=projected_s(1:dimension,1:dimension);lwork=-1;call disable_halting()
      call zheev('V','U',dimension,vectors,2*nstate,small_eigenvalues,query,lwork,rwork,lapack_info)
      if(lapack_info==0)then
        lwork=max(1,int(real(query(1))));if(allocated(work))deallocate(work);allocate(work(lwork),stat=allocation_status)
        if(allocation_status==0)then
          vectors(1:dimension,1:dimension)=projected_s(1:dimension,1:dimension)
          call zheev('V','U',dimension,vectors,2*nstate,small_eigenvalues,work,lwork,rwork,lapack_info)
        else;lapack_info=-999;endif
      endif
      call restore_halting();step_ok=lapack_info==0
      if(.not.step_ok.or.minval(small_eigenvalues(1:dimension))<=final_tolerance)then;step_ok=.false.;return;endif
      projected_h(1:dimension,1:dimension)=vectors(1:dimension,1:dimension)
      do column=1,dimension
        projected_h(1:dimension,column)=projected_h(1:dimension,column)/sqrt(small_eigenvalues(column))
      enddo
      projected_h(1:dimension,1:dimension)=matmul(projected_h(1:dimension,1:dimension),&
        conjg(transpose(vectors(1:dimension,1:dimension))))
      block(:,1:dimension)=matmul(block(:,1:dimension),projected_h(1:dimension,1:dimension))
    end subroutine orthonormalize
    subroutine orthogonalize_residual(block,retained_rank,step_ok)
      complex(real64),intent(inout)::block(:,:);integer,intent(out)::retained_rank;logical,intent(out)::step_ok
      complex(real64)::query(1);integer::lapack_info,column,source
      real(real64)::rank_threshold
      external::zheev
      call apply_s(block,residual_metric,callback_ok);step_ok=callback_ok;if(.not.step_ok)return
      overlap=matmul(conjg(transpose(c)),residual_metric);call MPI_Allreduce(MPI_IN_PLACE,overlap,nstate*nstate,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);block=block-matmul(c,overlap)
      call apply_s(block,residual_metric,callback_ok);step_ok=callback_ok;if(.not.step_ok)return
      projected_s=(0d0,0d0);projected_s(1:nstate,1:nstate)=matmul(conjg(transpose(block)),residual_metric)
      call MPI_Allreduce(MPI_IN_PLACE,projected_s,4*nstate*nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      vectors(1:nstate,1:nstate)=projected_s(1:nstate,1:nstate);lwork=-1;call disable_halting()
      call zheev('V','U',nstate,vectors,2*nstate,small_eigenvalues,query,lwork,rwork,lapack_info)
      if(lapack_info==0)then
        lwork=max(1,int(real(query(1))));if(allocated(work))deallocate(work);allocate(work(lwork),stat=allocation_status)
        if(allocation_status==0)then
          vectors(1:nstate,1:nstate)=projected_s(1:nstate,1:nstate)
          call zheev('V','U',nstate,vectors,2*nstate,small_eigenvalues,work,lwork,rwork,lapack_info)
        else;lapack_info=-999;endif
      endif
      call restore_halting();step_ok=lapack_info==0;if(.not.step_ok)return
      rank_threshold=max(final_tolerance**2,100d0*epsilon(1d0)**2*&
        max(1d0,maxval(abs(small_eigenvalues(1:nstate)))))
      retained_rank=count(small_eigenvalues(1:nstate)>rank_threshold)
      if(retained_rank<1)then;step_ok=.false.;return;endif
      projected_h(1:nstate,1:retained_rank)=(0d0,0d0);column=0
      do source=1,nstate
        if(small_eigenvalues(source)<=rank_threshold)cycle
        column=column+1;projected_h(1:nstate,column)=vectors(1:nstate,source)/sqrt(small_eigenvalues(source))
      enddo
      residual_metric(:,1:retained_rank)=matmul(block,projected_h(1:nstate,1:retained_rank))
      block(:,1:retained_rank)=residual_metric(:,1:retained_rank)
    end subroutine orthogonalize_residual
    subroutine disable_halting()
      call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_get_halting_mode(ieee_overflow,halt_overflow);call ieee_set_halting_mode(ieee_invalid,.false.)
      call ieee_set_halting_mode(ieee_divide_by_zero,.false.);call ieee_set_halting_mode(ieee_overflow,.false.)
    end subroutine disable_halting
    subroutine restore_halting()
      call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.);call ieee_set_flag(ieee_overflow,.false.)
      call ieee_set_halting_mode(ieee_invalid,halt_invalid);call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero)
      call ieee_set_halting_mode(ieee_overflow,halt_overflow)
    end subroutine restore_halting
    subroutine agree_integer(value,minimum_value,maximum_value,status)
      integer,intent(in)::value;integer,intent(out)::minimum_value,maximum_value,status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine agree_integer
    subroutine agree_int64(value,minimum_value,maximum_value,status)
      integer(int64),intent(in)::value;integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
    end subroutine agree_int64
    subroutine callback_consensus(callback_result,bad,status)
      logical,intent(in)::callback_result;integer,intent(out)::bad,status;integer::local
      local=merge(0,1,callback_result);call MPI_Allreduce(local,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine callback_consensus
    logical function finite_matrix(values)
      complex(real64),intent(in)::values(:,:)
      finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
    end function finite_matrix
    subroutine cleanup(keep_output)
      logical,intent(in),optional::keep_output;logical::keep
      keep=.false.;if(present(keep_output))keep=keep_output
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(c))deallocate(c);if(allocated(hc))deallocate(hc);if(allocated(sc))deallocate(sc)
      if(allocated(residual))deallocate(residual);if(allocated(residual_metric))deallocate(residual_metric)
      if(allocated(trial))deallocate(trial);if(allocated(htrial))deallocate(htrial);if(allocated(strial))deallocate(strial)
      if(allocated(projected_h))deallocate(projected_h);if(allocated(projected_s))deallocate(projected_s)
      if(allocated(vectors))deallocate(vectors);if(allocated(overlap))deallocate(overlap)
      if(allocated(work))deallocate(work);if(allocated(rwork))deallocate(rwork)
      if(allocated(small_eigenvalues))deallocate(small_eigenvalues)
      if(.not.keep.and.allocated(coefficients))deallocate(coefficients)
    end subroutine cleanup
#else
    ok=.false.;message='MPI is required for adaptive hybrid block CG';iterations=0;stop_reason='no_mpi'
    maximum_residual=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64;eigenvalues=0d0
#endif
  end subroutine solve_dg_hybrid_block_cg
end module dg_hybrid_block_cg
