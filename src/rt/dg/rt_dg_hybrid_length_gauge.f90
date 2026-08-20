#include "config.h"
module rt_dg_hybrid_length_gauge
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_sparse_metric,only:s_dg_hybrid_sparse_metric
  use dg_hybrid_sparse_operators,only:s_dg_hybrid_sparse_operators
  use rt_dg_hybrid_metric_solver,only:solve_rt_dg_hybrid_metric
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,build_rt_dg_sparse_exchange,&
    exchange_rt_dg_sparse_values,clear_rt_dg_sparse_exchange
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::propagate_rt_dg_hybrid_length_gauge
contains
  subroutine propagate_rt_dg_hybrid_length_gauge(comm,metric,operators,coefficients_owned,electric_field,&
      time_step,tolerance,max_order,previous_polarization,polarization_periods,next_coefficients_owned,&
      metric_norm,energy,polarization,metric_iterations,workspace_peak_bytes,fingerprint,ok,message,&
      cached_metric_exchange,cached_operator_exchange)
    integer,intent(in)::comm,max_order
    type(s_dg_hybrid_sparse_metric),intent(in)::metric
    type(s_dg_hybrid_sparse_operators),intent(in)::operators
    complex(real64),intent(in)::coefficients_owned(:)
    real(real64),intent(in)::electric_field(3),time_step,tolerance,previous_polarization(3),polarization_periods(3)
    complex(real64),allocatable,intent(out)::next_coefficients_owned(:)
    real(real64),intent(out)::metric_norm,energy,polarization(3)
    integer,intent(out)::metric_iterations
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_rt_dg_sparse_exchange),optional,target,intent(inout)::cached_metric_exchange,cached_operator_exchange
#ifdef USE_MPI
    integer::n,nowned,order,i,k,row,ierr,local_bad,global_bad,allocation_status,inner_iterations,inner_cap,max_degree
    integer::minimum_integer,maximum_integer,reduction_pair(2)
    integer(int64)::bits,minimum_bits,maximum_bits,base_bytes,solver_bytes,solver_fingerprint,quantized,&
      local_hash,global_hash,row_hash
    complex(real64),allocatable::term(:),rhs(:,:),solved(:,:),operator_edge_values(:),metric_edge_values(:),&
      applied(:),metric_applied(:)
    complex(real64)::local_dot,global_dot
    real(real64)::term_norm,result_norm,local_real,inner_residual,scale,safe_scale,wrapped,quantization_limit,&
      operator_scale,hamiltonian_scale,position_scale,field_scale
    logical::converged,solver_ok,exchange_ok,owns_metric_exchange,owns_operator_exchange
    character(256)::solver_message,exchange_message
    type(s_rt_dg_sparse_exchange),target::local_metric_exchange,local_operator_exchange
    type(s_rt_dg_sparse_exchange),pointer::metric_exchange,operator_exchange
    ok=.false.;message='';metric_norm=huge(1d0);energy=huge(1d0);polarization=huge(1d0)
    metric_iterations=0;workspace_peak_bytes=0_int64;fingerprint=0_int64
    owns_metric_exchange=.not.present(cached_metric_exchange)
    owns_operator_exchange=.not.present(cached_operator_exchange)
    n=metric%global_count;nowned=size(metric%owned_row_ids);local_bad=0
    call agree_integer(n,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent length-gauge basis extent';return;endif
    call agree_integer(max_order,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent length-gauge expansion order';return;endif
    bits=transfer(time_step,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent length-gauge time step';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent length-gauge tolerance';return;endif
    call agree_int64(metric%fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.metric%fingerprint==0_int64)then
      message='inconsistent length-gauge metric provenance';return
    endif
    call agree_int64(operators%fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.operators%fingerprint==0_int64)then
      message='inconsistent length-gauge operator provenance';return
    endif
    do i=1,3
      bits=transfer(electric_field(i),bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent length-gauge electric field';return;endif
      bits=transfer(previous_polarization(i),bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent previous polarization';return;endif
      bits=transfer(polarization_periods(i),bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent polarization period';return;endif
    enddo
    if(.not.metric%valid.or..not.operators%valid.or.n<1.or.max_order<2)local_bad=1
    if(operators%global_count/=n.or.operators%metric_fingerprint/=metric%fingerprint)local_bad=1
    if(size(coefficients_owned)/=nowned.or.size(operators%owned_row_ids)/=nowned)local_bad=1
    if(nowned==huge(0))local_bad=1
    if(local_bad==0)then
      if(size(metric%row_offsets)/=nowned+1.or.size(operators%row_offsets)/=nowned+1)local_bad=1
    endif
    if(any(metric%owned_row_ids/=operators%owned_row_ids))local_bad=1
    if(local_bad==0)then
      do i=1,nowned
        row=int(metric%owned_row_ids(i))
        if(.not.metric%active_rows(row).and.abs(coefficients_owned(i))>tolerance)local_bad=1
      enddo
    endif
    if(size(metric%column_ids)/=size(metric%values))local_bad=1
    if(size(operators%column_ids)/=size(operators%hamiltonian_values).or.&
      size(operators%metric_values)/=size(operators%column_ids))local_bad=1
    if(any(shape(operators%position_values)/=[3,size(operators%column_ids)]))local_bad=1
    if(.not.finite_vector(coefficients_owned).or..not.finite_vector(metric%values).or.&
      .not.finite_vector(operators%hamiltonian_values).or..not.finite_matrix(operators%position_values))local_bad=1
    if(.not.all(ieee_is_finite(electric_field)).or..not.ieee_is_finite(time_step).or..not.ieee_is_finite(tolerance).or.&
      .not.all(ieee_is_finite(previous_polarization)).or..not.all(ieee_is_finite(polarization_periods)))local_bad=1
    if(time_step<=0d0.or.tolerance<1d-15.or.tolerance>1d-2.or.any(polarization_periods<=0d0))local_bad=1
    if(time_step>sqrt(huge(1d0)).or.any(abs(electric_field)>sqrt(huge(1d0))).or.&
      any(abs(previous_polarization)>sqrt(huge(1d0))).or.any(polarization_periods<sqrt(tiny(1d0))))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid generalized length-gauge contract';return;endif
    max_degree=0;hamiltonian_scale=0d0;position_scale=0d0
    if(metric%row_offsets(1)/=1.or.metric%row_offsets(nowned+1)/=size(metric%values)+1)local_bad=1
    if(operators%row_offsets(1)/=1.or.operators%row_offsets(nowned+1)/=size(operators%column_ids)+1)local_bad=1
    do i=1,nowned
      if(metric%row_offsets(i)<1.or.metric%row_offsets(i+1)<metric%row_offsets(i).or.&
        metric%row_offsets(i+1)>size(metric%values)+1)local_bad=1
      if(operators%row_offsets(i)<1.or.operators%row_offsets(i+1)<operators%row_offsets(i).or.&
        operators%row_offsets(i+1)>size(operators%column_ids)+1)local_bad=1
      max_degree=max(max_degree,operators%row_offsets(i+1)-operators%row_offsets(i))
    enddo
    if(any(metric%column_ids<1).or.any(metric%column_ids>n).or.any(operators%column_ids<1).or.&
      any(operators%column_ids>n))local_bad=1
    if(size(operators%column_ids)>0)then
      hamiltonian_scale=max(maxval(abs(real(operators%hamiltonian_values))),&
        maxval(abs(aimag(operators%hamiltonian_values))))
      position_scale=max(maxval(abs(real(operators%position_values))),&
        maxval(abs(aimag(operators%position_values))))
      if(hamiltonian_scale>huge(1d0)/2d0.or.position_scale>huge(1d0)/2d0)local_bad=1
    endif
    reduction_pair=[max_degree,local_bad]
    call MPI_Allreduce(MPI_IN_PLACE,reduction_pair,2,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='length-gauge degree reduction failed';return;endif
    max_degree=reduction_pair(1);global_bad=reduction_pair(2)
    if(global_bad/=0)then;message='invalid length-gauge sparse graph';return;endif
    field_scale=sum(abs(electric_field));operator_scale=2d0*hamiltonian_scale
    if(position_scale>0d0)then
      if(field_scale>(huge(1d0)-operator_scale)/(2d0*position_scale))then
        message='unsafe length-gauge field/operator scale';return
      endif
      operator_scale=operator_scale+2d0*field_scale*position_scale
    endif
    if(int(nowned,int64)>huge(base_bytes)/80_int64)then
      local_bad=1
    else
      base_bytes=80_int64*int(nowned,int64)
      if(int(size(metric%column_ids),int64)>huge(base_bytes)/16_int64.or.&
        int(size(operators%column_ids),int64)>huge(base_bytes)/16_int64)local_bad=1
      if(local_bad==0)then
        if(16_int64*int(size(metric%column_ids),int64)>huge(base_bytes)-base_bytes)local_bad=1
      endif
      if(local_bad==0)base_bytes=base_bytes+16_int64*int(size(metric%column_ids),int64)
      if(local_bad==0)then
        if(16_int64*int(size(operators%column_ids),int64)>huge(base_bytes)-base_bytes)local_bad=1
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='length-gauge workspace receipt overflow';return;endif
    base_bytes=base_bytes+16_int64*int(size(operators%column_ids),int64)
    allocate(next_coefficients_owned(nowned),term(nowned),rhs(nowned,1),operator_edge_values(size(operators%column_ids)),&
      metric_edge_values(size(metric%column_ids)),applied(nowned),metric_applied(nowned),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate length-gauge workspace';return;endif
    if(present(cached_metric_exchange))then
      metric_exchange=>cached_metric_exchange
      local_bad=merge(0,1,metric_exchange%valid.and.metric_exchange%catalog_fingerprint==metric%fingerprint)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='stale cached metric halo';return;endif
    else
      metric_exchange=>local_metric_exchange
      call build_rt_dg_sparse_exchange(comm,n,metric%fingerprint,metric%owned_row_ids,metric%column_ids,&
        metric_exchange,exchange_ok,exchange_message)
      if(.not.exchange_ok)then;call cleanup();message='metric halo setup failed: '//trim(exchange_message);return;endif
    endif
    if(present(cached_operator_exchange))then
      operator_exchange=>cached_operator_exchange
      local_bad=merge(0,1,operator_exchange%valid.and.operator_exchange%catalog_fingerprint==operators%fingerprint)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='stale cached operator halo';return;endif
    else
      operator_exchange=>local_operator_exchange
      call build_rt_dg_sparse_exchange(comm,n,operators%fingerprint,operators%owned_row_ids,operators%column_ids,&
        operator_exchange,exchange_ok,exchange_message)
      if(.not.exchange_ok)then;call cleanup();message='operator halo setup failed: '//trim(exchange_message);return;endif
    endif
    local_bad=0
    if(metric_exchange%workspace_peak_bytes>huge(base_bytes)-base_bytes)then
      local_bad=1
    else
      base_bytes=base_bytes+metric_exchange%workspace_peak_bytes
    endif
    if(local_bad==0)then
      if(operator_exchange%workspace_peak_bytes>huge(base_bytes)-base_bytes)then
        local_bad=1
      else
        base_bytes=base_bytes+operator_exchange%workspace_peak_bytes
      endif
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='length-gauge halo receipt overflow';return;endif
    local_real=0d0;if(nowned>0)local_real=maxval(abs(coefficients_owned))
    safe_scale=sqrt(huge(1d0))/(16d0*sqrt(real(n,real64)))
    safe_scale=safe_scale/real(max(1,max_degree),real64)/max(1d0,operator_scale)
    local_bad=merge(0,1,local_real<=safe_scale)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='unsafe length-gauge coefficient magnitude';return;endif
    next_coefficients_owned=coefficients_owned;term=coefficients_owned;converged=.false.;fingerprint=operators%fingerprint
    do i=1,nowned
      row=int(metric%owned_row_ids(i))
      if(.not.metric%active_rows(row))then
        next_coefficients_owned(i)=(0d0,0d0);term(i)=(0d0,0d0)
      endif
    enddo
    if(n>(huge(0)-50)/4)then;inner_cap=huge(0);else;inner_cap=max(50,4*n);endif
    do order=1,max_order
      local_real=0d0;if(nowned>0)local_real=maxval(abs(term))
      call MPI_Allreduce(local_real,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge coefficient scale failed';return;endif
      if(scale>safe_scale)then;call cleanup();message='unsafe length-gauge exponential term magnitude';return;endif
      call exchange_rt_dg_sparse_values(comm,operator_exchange,term,operator_edge_values,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge coefficient exchange failed';return;endif
      call apply_operator(operator_edge_values,electric_field,applied)
      rhs(:,1)=applied
      call solve_rt_dg_hybrid_metric(comm,metric,rhs,tolerance,inner_cap,solved,inner_iterations,inner_residual,&
        solver_bytes,solver_fingerprint,solver_ok,solver_message,metric_exchange)
      if(.not.solver_ok)then;call cleanup();message='length-gauge metric solve failed: '//trim(solver_message);return;endif
      metric_iterations=max(metric_iterations,inner_iterations)
      if(solver_bytes>huge(workspace_peak_bytes)-base_bytes)then
        deallocate(solved);call cleanup();message='length-gauge combined workspace receipt overflow';return
      endif
      workspace_peak_bytes=max(workspace_peak_bytes,base_bytes+solver_bytes)
      term=cmplx(0d0,-time_step/real(order,real64),real64)*solved(:,1);deallocate(solved)
      next_coefficients_owned=next_coefficients_owned+term
      if(.not.finite_vector(term).or..not.finite_vector(next_coefficients_owned))local_bad=1
      local_real=0d0;if(nowned>0)local_real=max(maxval(abs(term)),maxval(abs(next_coefficients_owned)))
      if(local_real>sqrt(huge(1d0))/(16d0*sqrt(real(n,real64))))local_bad=1
      local_real=sum(abs(term)**2);call MPI_Allreduce(local_real,term_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      local_real=sum(abs(next_coefficients_owned)**2);call MPI_Allreduce(local_real,result_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='nonfinite length-gauge exponential term';return;endif
      fingerprint=ieor(ishftc(fingerprint,7),solver_fingerprint)
      if(sqrt(term_norm)<=tolerance*max(1d0,sqrt(result_norm)))then;converged=.true.;exit;endif
    enddo
    if(.not.converged)then;call cleanup();message='length-gauge exponential expansion did not converge';return;endif
    local_real=0d0;if(nowned>0)local_real=maxval(abs(next_coefficients_owned))
    call MPI_Allreduce(local_real,scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge observable scale failed';return;endif
    if(scale>safe_scale)then;call cleanup();message='unsafe propagated observable magnitude';return;endif
    call exchange_rt_dg_sparse_values(comm,metric_exchange,next_coefficients_owned,metric_edge_values,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge metric exchange failed';return;endif
    call apply_metric(metric_edge_values,metric_applied)
    local_dot=sum(conjg(next_coefficients_owned)*metric_applied)
    call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    metric_norm=real(global_dot)
    if(ierr/=MPI_SUCCESS.or.metric_norm<=tiny(1d0).or.abs(aimag(global_dot))>100d0*tolerance*metric_norm)then
      call cleanup();message='invalid propagated metric norm';return
    endif
    call exchange_rt_dg_sparse_values(comm,operator_exchange,next_coefficients_owned,operator_edge_values,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge operator exchange failed';return;endif
    call apply_operator(operator_edge_values,[0d0,0d0,0d0],applied)
    local_dot=sum(conjg(next_coefficients_owned)*applied)
    call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr);energy=real(global_dot)/metric_norm
    if(ierr/=MPI_SUCCESS.or.abs(aimag(global_dot))>100d0*tolerance*max(1d0,abs(real(global_dot))))then
      call cleanup();message='invalid propagated energy';return
    endif
    do i=1,3
      call apply_position(operator_edge_values,i,applied)
      local_dot=sum(conjg(next_coefficients_owned)*applied)
      call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.abs(aimag(global_dot))>100d0*tolerance*max(1d0,abs(real(global_dot))))then
        call cleanup();message='invalid propagated polarization';return
      endif
      wrapped=real(global_dot)/metric_norm
      polarization(i)=wrapped+anint((previous_polarization(i)-wrapped)/polarization_periods(i))*polarization_periods(i)
    enddo
    quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*tolerance
    if(scale>quantization_limit)then;call cleanup();message='length-gauge fingerprint range is unsafe';return;endif
    local_hash=0_int64
    do i=1,nowned
      row_hash=ishftc(int(metric%owned_row_ids(i),int64),11)
      quantized=nint(real(next_coefficients_owned(i))/(100d0*tolerance),int64)
      row_hash=ieor(ishftc(row_hash,9),quantized)
      quantized=nint(aimag(next_coefficients_owned(i))/(100d0*tolerance),int64)
      row_hash=ieor(ishftc(row_hash,9),quantized);local_hash=ieor(local_hash,row_hash)
    enddo
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='length-gauge fingerprint reduction failed';return;endif
    fingerprint=ieor(ishftc(fingerprint,9),global_hash)
    if(fingerprint==0_int64)fingerprint=1_int64;ok=.true.
#else
    ok=.false.;message='hybrid length-gauge propagation requires MPI';metric_norm=huge(1d0);energy=huge(1d0)
    polarization=huge(1d0);metric_iterations=0;workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine apply_metric(values_by_edge,local_values)
      complex(real64),intent(in)::values_by_edge(:);complex(real64),intent(out)::local_values(:)
      local_values=(0d0,0d0)
      do k=1,nowned
        if(.not.metric%active_rows(int(metric%owned_row_ids(k))))cycle
        do row=metric%row_offsets(k),metric%row_offsets(k+1)-1
          if(.not.metric%active_rows(metric%column_ids(row)))cycle
          local_values(k)=local_values(k)+metric%values(row)*values_by_edge(row)
        enddo
      enddo
    end subroutine apply_metric
    subroutine apply_operator(values_by_edge,field,local_values)
      complex(real64),intent(in)::values_by_edge(:);real(real64),intent(in)::field(3)
      complex(real64),intent(out)::local_values(:)
      local_values=(0d0,0d0)
      do k=1,nowned
        if(.not.metric%active_rows(int(metric%owned_row_ids(k))))cycle
        do row=operators%row_offsets(k),operators%row_offsets(k+1)-1
          if(.not.metric%active_rows(operators%column_ids(row)))cycle
          local_values(k)=local_values(k)+(operators%hamiltonian_values(row)+&
            sum(field*operators%position_values(:,row)))*values_by_edge(row)
        enddo
      enddo
    end subroutine apply_operator
    subroutine apply_position(values_by_edge,component,local_values)
      complex(real64),intent(in)::values_by_edge(:);integer,intent(in)::component
      complex(real64),intent(out)::local_values(:)
      local_values=(0d0,0d0)
      do k=1,nowned
        if(.not.metric%active_rows(int(metric%owned_row_ids(k))))cycle
        do row=operators%row_offsets(k),operators%row_offsets(k+1)-1
          if(.not.metric%active_rows(operators%column_ids(row)))cycle
          local_values(k)=local_values(k)+operators%position_values(component,row)*values_by_edge(row)
        enddo
      enddo
    end subroutine apply_position
    subroutine cleanup()
      if(allocated(next_coefficients_owned))deallocate(next_coefficients_owned)
      if(allocated(term))deallocate(term)
      if(allocated(rhs))deallocate(rhs)
      if(allocated(solved))deallocate(solved)
      if(allocated(operator_edge_values))deallocate(operator_edge_values)
      if(allocated(metric_edge_values))deallocate(metric_edge_values)
      if(allocated(applied))deallocate(applied)
      if(allocated(metric_applied))deallocate(metric_applied)
      if(owns_metric_exchange)call clear_rt_dg_sparse_exchange(local_metric_exchange)
      if(owns_operator_exchange)call clear_rt_dg_sparse_exchange(local_operator_exchange)
    end subroutine cleanup
#endif
  end subroutine propagate_rt_dg_hybrid_length_gauge
  logical function finite_vector(values)
    complex(real64),intent(in)::values(:)
    finite_vector=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_vector
  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix
#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm;integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer
  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value;integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module rt_dg_hybrid_length_gauge
