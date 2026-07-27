#include "config.h"
module dg_overlapping_wannier_scf
  use iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_overlapping_wannier_solver,only:solve_dg_overlapping_wannier_coefficients
  use dg_overlapping_wannier_density,only:reconstruct_dg_overlapping_wannier_density
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_overlapping_wannier_scf_state
    real(8),allocatable::density(:),potential(:),eigenvalues(:),density_history(:,:)
    complex(8),allocatable::coefficients(:,:)
    integer::history_count=0,basis_generation=0,geometry_generation=0
    integer(int64)::basis_fingerprint=0_int64,operator_fingerprint=0_int64
    logical::accepted=.false.
  end type
  type,public::s_dg_overlapping_wannier_scf_result
    logical::converged=.false.
    integer::iterations=0,hamiltonian_rebuilds=0
    real(8)::density_residual=huge(1d0),unmixed_density_residual=huge(1d0)
    real(8)::coefficient_residual=huge(1d0),orthogonality_defect=huge(1d0)
    real(8)::integrated_charge=huge(1d0),trace_charge=huge(1d0)
    real(8)::symmetry_closure_residual=huge(1d0)
  end type
  abstract interface
    subroutine dg_overlapping_wannier_hamiltonian_builder(comm,density,potential,hrows,new_potential,&
        fingerprint,ok,message)
      import int64
      integer,intent(in)::comm
      real(8),intent(in)::density(:),potential(:)
      complex(8),intent(out)::hrows(:,:)
      real(8),intent(out)::new_potential(:)
      integer(int64),intent(out)::fingerprint
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine
    subroutine dg_overlapping_wannier_density_mixer(comm,mixing,current_density,raw_density,history,&
        history_count,mixed_density,new_history,new_history_count,ok,message)
      integer,intent(in)::comm,history_count
      real(8),intent(in)::mixing,current_density(:),raw_density(:),history(:,:)
      real(8),intent(out)::mixed_density(:),new_history(:,:)
      integer,intent(out)::new_history_count
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine
    subroutine dg_overlapping_wannier_transaction(action,ok,message)
      integer,intent(in)::action
      logical,intent(out)::ok
      character(*),intent(out)::message
    end subroutine
    subroutine dg_overlapping_wannier_transaction_commit()
    end subroutine
  end interface
  public::run_dg_overlapping_wannier_scf
#ifdef USE_MPI
  public::compute_dg_overlapping_wannier_scf_fingerprint
#endif
contains
  subroutine run_dg_overlapping_wannier_scf(comm,row_ids,srows,physical_ids,weights,values,tail_generation,&
      expected_generation,expected_geometry_generation,expected_basis_fingerprint,occupations,&
      symmetry_closure_residual,symmetry_tolerance,symmetry_fingerprint,expected_core_count,mixing,max_outer,max_inner,&
      density_tolerance,coefficient_tolerance,&
      build_hamiltonian,mix_density,transaction,commit_transaction,state,result,ok,message)
    integer,intent(in)::comm,expected_generation,expected_geometry_generation,max_outer,max_inner
    integer(int64),intent(in)::expected_basis_fingerprint
    integer(int64),intent(in)::symmetry_fingerprint
    integer(int64),intent(in)::row_ids(:),physical_ids(:),expected_core_count
    complex(8),intent(in)::srows(:,:),values(:,:)
    real(8),intent(in)::weights(:),occupations(:),symmetry_closure_residual,symmetry_tolerance,mixing,&
      density_tolerance,coefficient_tolerance
    integer,intent(in)::tail_generation(:,:)
    procedure(dg_overlapping_wannier_hamiltonian_builder)::build_hamiltonian
    procedure(dg_overlapping_wannier_density_mixer)::mix_density
    procedure(dg_overlapping_wannier_transaction)::transaction
    procedure(dg_overlapping_wannier_transaction_commit)::commit_transaction
    type(s_dg_overlapping_wannier_scf_state),intent(inout)::state
    type(s_dg_overlapping_wannier_scf_result),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_overlapping_wannier_scf_state)::candidate
    complex(8),allocatable::hrows(:,:),overlap(:,:),coeff(:,:),check_coeff(:,:)
    real(8),allocatable::raw_density(:),check_density(:),mixed_density(:),new_potential(:),check_potential(:),&
      eval(:),check_eval(:),new_history(:,:)
    real(8)::condition,charge,trace_charge,residual,orthogonality,check_residual,check_orthogonality
    integer::iteration,nwann,nstate,i,ierr,local_bad,global_bad,new_history_count
    integer::integer_contract(7),integer_min(7),integer_max(7)
    integer(int64)::fingerprint,check_fingerprint,computed_basis_fingerprint,core_min,core_max,basis_min,basis_max,&
      operator_min,operator_max,symmetry_min,symmetry_max
    real(8)::real_contract(5),real_min(5),real_max(5)
    integer,allocatable::row_owners(:)
    logical::step_ok
    character(256)::step_message
    result=s_dg_overlapping_wannier_scf_result();ok=.false.;message=''
    nwann=size(srows,2);nstate=size(occupations);local_bad=0
    integer_contract=[nwann,nstate,max_outer,max_inner,expected_generation,expected_geometry_generation,&
      state%history_count]
    real_contract=[mixing,density_tolerance,coefficient_tolerance,symmetry_closure_residual,symmetry_tolerance]
    call MPI_Allreduce(integer_contract,integer_min,7,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(integer_contract,integer_max,7,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(real_contract,real_min,5,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    call MPI_Allreduce(real_contract,real_max,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_core_count,core_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_core_count,core_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(expected_basis_fingerprint,basis_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(expected_basis_fingerprint,basis_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(state%operator_fingerprint,operator_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(state%operator_fingerprint,operator_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    call MPI_Allreduce(symmetry_fingerprint,symmetry_min,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(symmetry_fingerprint,symmetry_max,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(any(integer_min/=integer_max).or.any(real_min/=real_max).or.core_min/=core_max.or.&
       basis_min/=basis_max.or.operator_min/=operator_max.or.symmetry_min/=symmetry_max.or.&
       symmetry_fingerprint==0_int64)then
      message='rank-inconsistent overlapping-Wannier SCF contract';return
    endif
    if(nwann<1.or.nstate<1.or.nstate>nwann.or.max_outer<1.or.max_inner<1)local_bad=1
    if(int(nwann,int64)*int(nwann,int64)>100000000_int64)local_bad=1
    if(mixing<=0d0.or.mixing>1d0.or.density_tolerance<=0d0.or.coefficient_tolerance<=0d0)local_bad=1
    if(expected_core_count<=0_int64.or.expected_generation<=0.or.expected_geometry_generation<=0.or.&
       expected_basis_fingerprint==0_int64)local_bad=1
    if(symmetry_tolerance<=0d0.or.symmetry_closure_residual<0d0.or.&
       symmetry_closure_residual>symmetry_tolerance)local_bad=1
    if(size(srows,1)/=size(row_ids).or.size(values,1)/=nwann.or.size(values,2)/=size(physical_ids).or.&
       size(weights)/=size(physical_ids).or.any(shape(tail_generation)/=shape(values)))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nwann,int64)))local_bad=1
    if(.not.finite_complex(srows).or..not.finite_complex(values).or..not.all(ieee_is_finite(weights)).or.&
       .not.all(ieee_is_finite(occupations)).or..not.all(ieee_is_finite(real_contract)))local_bad=1
    if(.not.allocated(state%density).or..not.allocated(state%potential).or.&
       .not.allocated(state%coefficients).or..not.allocated(state%eigenvalues).or.&
       .not.allocated(state%density_history))local_bad=1
    if(local_bad==0)then
      if(size(state%density)/=size(physical_ids).or.size(state%potential)/=size(physical_ids))local_bad=1
      if(any(shape(state%coefficients)/=[nwann,nstate]).or.size(state%eigenvalues)/=nstate)local_bad=1
      if(size(state%density_history,1)/=size(physical_ids).or.size(state%density_history,2)<1)local_bad=1
      if(state%history_count<1.or.state%history_count>size(state%density_history,2))local_bad=1
      if(state%basis_generation/=expected_generation.or.state%geometry_generation<=0.or.&
         state%basis_fingerprint==0_int64.or.state%operator_fingerprint==0_int64)local_bad=1
      if(state%geometry_generation/=expected_geometry_generation.or.&
         state%basis_fingerprint/=expected_basis_fingerprint)local_bad=1
      if(.not.all(ieee_is_finite(state%density)).or..not.all(ieee_is_finite(state%potential)).or.&
         .not.all(ieee_is_finite(state%eigenvalues)).or..not.all(ieee_is_finite(state%density_history)).or.&
         .not.finite_complex(state%coefficients))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then;message='invalid overlapping-Wannier SCF transaction';return;endif
    call compute_dg_overlapping_wannier_scf_fingerprint(comm,row_ids,srows,physical_ids,weights,values,&
      tail_generation,computed_basis_fingerprint,symmetry_fingerprint)
    if(computed_basis_fingerprint/=expected_basis_fingerprint)then
      message='overlapping-Wannier S/basis content fingerprint mismatch';return
    endif
    allocate(row_owners(nwann));row_owners=0
    do i=1,size(row_ids);row_owners(int(row_ids(i)))=row_owners(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_owners,nwann,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(any(row_owners/=1))then;message='SCF coefficient rows do not form a unique global partition';return;endif
    call transaction(0,step_ok,step_message)
    call collective_step_status(comm,step_ok,step_message,ok,message)
    if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
    candidate=state;result%symmetry_closure_residual=symmetry_closure_residual
    allocate(hrows(size(row_ids),nwann),overlap(nwann,nwann),coeff(nwann,nstate),check_coeff(nwann,nstate))
    allocate(raw_density(size(physical_ids)),check_density(size(physical_ids)),mixed_density(size(physical_ids)))
    allocate(new_potential(size(state%potential)),check_potential(size(state%potential)))
    allocate(eval(nstate),check_eval(nstate),new_history(size(state%density_history,1),&
      size(state%density_history,2)));overlap=(0d0,0d0)
    do i=1,size(row_ids);overlap(int(row_ids(i)),:)=srows(i,:);enddo
    call MPI_Allreduce(MPI_IN_PLACE,overlap,nwann*nwann,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    do iteration=1,max_outer
      call build_hamiltonian(comm,candidate%density,candidate%potential,hrows,new_potential,fingerprint,step_ok,step_message)
      result%hamiltonian_rebuilds=result%hamiltonian_rebuilds+1
      call collective_step_status(comm,step_ok,step_message,ok,message)
      if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
      call validate_builder_output(comm,hrows,new_potential,fingerprint,step_ok,step_message)
      if(.not.step_ok)then;message=step_message;call rollback_transaction(comm,transaction);return;endif
      call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,nstate,max_inner,&
        coefficient_tolerance,coefficient_tolerance,coeff,eval,residual,orthogonality,condition,step_ok,step_message)
      call collective_step_status(comm,step_ok,step_message,ok,message)
      if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
      call reconstruct_dg_overlapping_wannier_density(comm,physical_ids,weights,values,tail_generation,&
        expected_generation,coeff,occupations,expected_core_count,overlap,density_tolerance,raw_density,&
        charge,trace_charge,step_ok,step_message)
      call collective_step_status(comm,step_ok,step_message,ok,message)
      if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
      result%density_residual=relative_density(comm,raw_density,candidate%density)
      result%coefficient_residual=residual;result%orthogonality_defect=orthogonality
      result%integrated_charge=charge;result%trace_charge=trace_charge;result%iterations=iteration
      if(result%density_residual<density_tolerance.and.residual<coefficient_tolerance.and.&
          orthogonality<10d0*coefficient_tolerance)exit
      call mix_density(comm,mixing,candidate%density,raw_density,candidate%density_history,&
        candidate%history_count,mixed_density,new_history,new_history_count,step_ok,step_message)
      call collective_step_status(comm,step_ok,step_message,ok,message)
      if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
      call validate_mixer_output(comm,mixed_density,new_history,new_history_count,step_ok,step_message)
      if(.not.step_ok)then;message=step_message;call rollback_transaction(comm,transaction);return;endif
      candidate%density=mixed_density
      candidate%potential=new_potential;candidate%coefficients=coeff;candidate%eigenvalues=eval
      candidate%operator_fingerprint=fingerprint
      candidate%density_history=new_history;candidate%history_count=new_history_count
    enddo
    if(iteration>max_outer)then
      message='overlapping-Wannier SCF did not converge; no fallback'
      call rollback_transaction(comm,transaction);return
    endif
    call build_hamiltonian(comm,raw_density,new_potential,hrows,check_potential,check_fingerprint,step_ok,step_message)
    result%hamiltonian_rebuilds=result%hamiltonian_rebuilds+1
    call collective_step_status(comm,step_ok,step_message,ok,message)
    if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
    call validate_builder_output(comm,hrows,check_potential,check_fingerprint,step_ok,step_message)
    if(.not.step_ok)then;message=step_message;call rollback_transaction(comm,transaction);return;endif
    call solve_dg_overlapping_wannier_coefficients(comm,row_ids,hrows,srows,nstate,max_inner,&
      coefficient_tolerance,coefficient_tolerance,check_coeff,check_eval,check_residual,check_orthogonality,&
      condition,step_ok,step_message)
    call collective_step_status(comm,step_ok,step_message,ok,message)
    if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
    call reconstruct_dg_overlapping_wannier_density(comm,physical_ids,weights,values,tail_generation,&
      expected_generation,check_coeff,occupations,expected_core_count,overlap,density_tolerance,check_density,&
      charge,trace_charge,step_ok,step_message)
    call collective_step_status(comm,step_ok,step_message,ok,message)
    if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
    result%unmixed_density_residual=relative_density(comm,check_density,raw_density)
    if(result%unmixed_density_residual>=density_tolerance.or.check_residual>=coefficient_tolerance.or.&
       check_orthogonality>=10d0*coefficient_tolerance)then
      message='overlapping-Wannier unmixed fixed-point gate failed; no fallback'
      call rollback_transaction(comm,transaction);return
    endif
    candidate%density=check_density;candidate%potential=check_potential
    candidate%coefficients=check_coeff;candidate%eigenvalues=check_eval
    candidate%operator_fingerprint=check_fingerprint;candidate%accepted=.true.
    call transaction(1,step_ok,step_message)
    call collective_step_status(comm,step_ok,step_message,ok,message)
    if(.not.ok)then;call rollback_transaction(comm,transaction);return;endif
    call commit_transaction()
    state=candidate;result%converged=.true.;result%coefficient_residual=check_residual
    result%orthogonality_defect=check_orthogonality;result%integrated_charge=charge;result%trace_charge=trace_charge
    ok=.true.;message=''
#else
    result=s_dg_overlapping_wannier_scf_result();ok=.false.
    message='overlapping-Wannier SCF requires MPI'
#endif
  end subroutine

#ifdef USE_MPI
  subroutine compute_dg_overlapping_wannier_scf_fingerprint(comm,row_ids,srows,physical_ids,weights,&
      values,tail_generation,fingerprint,symmetry_fingerprint)
    integer,intent(in)::comm
    integer(int64),intent(in)::row_ids(:),physical_ids(:)
    complex(8),intent(in)::srows(:,:),values(:,:)
    real(8),intent(in)::weights(:)
    integer,intent(in)::tail_generation(:,:)
    integer(int64),intent(out)::fingerprint
    integer(int64),intent(in),optional::symmetry_fingerprint
    integer(int64)::local_hash,bits
    integer::i,j,ierr,shift
    local_hash=0_int64
    do i=1,size(row_ids)
      shift=int(modulo(row_ids(i),63_int64))
      local_hash=ieor(local_hash,ishftc(row_ids(i),shift))
      do j=1,size(srows,2)
        bits=transfer(real(srows(i,j),8),bits);shift=mod(shift+7*j,63)
        local_hash=ieor(local_hash,ishftc(bits,shift))
        bits=transfer(aimag(srows(i,j)),bits);local_hash=ieor(local_hash,ishftc(bits,mod(shift+19,63)))
      enddo
    enddo
    do i=1,size(physical_ids)
      shift=int(modulo(physical_ids(i),63_int64))
      local_hash=ieor(local_hash,ishftc(physical_ids(i),shift))
      bits=transfer(weights(i),bits);local_hash=ieor(local_hash,ishftc(bits,mod(shift+31,63)))
      do j=1,size(values,1)
        bits=transfer(real(values(j,i),8),bits);local_hash=ieor(local_hash,ishftc(bits,mod(shift+11*j,63)))
        bits=transfer(aimag(values(j,i)),bits);local_hash=ieor(local_hash,ishftc(bits,mod(shift+17*j,63)))
        local_hash=ieor(local_hash,ishftc(int(tail_generation(j,i),int64),mod(shift+29*j,63)))
      enddo
    enddo
    call MPI_Allreduce(local_hash,fingerprint,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(present(symmetry_fingerprint))fingerprint=ieor(fingerprint,symmetry_fingerprint)
    fingerprint=ieor(fingerprint,int(z'6A09E667F3BCC909',int64))
    if(fingerprint==0_int64)fingerprint=1_int64
  end subroutine

  subroutine rollback_transaction(comm,transaction)
    integer,intent(in)::comm
    procedure(dg_overlapping_wannier_transaction)::transaction
    logical::rollback_ok
    character(256)::rollback_message
    integer::bad,global_bad,ierr
    call transaction(-1,rollback_ok,rollback_message)
    bad=merge(0,1,rollback_ok)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)error stop 'overlapping-Wannier external rollback failed collectively'
  end subroutine

  subroutine validate_builder_output(comm,hrows,potential,fingerprint,ok,message)
    integer,intent(in)::comm
    complex(8),intent(in)::hrows(:,:)
    real(8),intent(in)::potential(:)
    integer(int64),intent(in)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::bad,global_bad,ierr
    integer(int64)::minimum_fingerprint,maximum_fingerprint
    bad=0
    if(.not.finite_complex(hrows).or..not.all(ieee_is_finite(potential)).or.fingerprint==0_int64)bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call MPI_Allreduce(fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    ok=global_bad==0.and.minimum_fingerprint==maximum_fingerprint
    if(ok)then;message='';else;message='invalid or rank-inconsistent staged DC Hamiltonian';endif
  end subroutine

  subroutine validate_mixer_output(comm,density,history,history_count,ok,message)
    integer,intent(in)::comm,history_count
    real(8),intent(in)::density(:),history(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::bad,global_bad,count_min,count_max,ierr
    bad=0
    if(history_count<1.or.history_count>size(history,2))bad=1
    if(.not.all(ieee_is_finite(density)).or..not.all(ieee_is_finite(history)))bad=1
    if(any(density<0d0))bad=1
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(history_count,count_min,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(history_count,count_max,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=global_bad==0.and.count_min==count_max
    if(ok)then;message='';else;message='invalid or rank-inconsistent seeded DC mixer output';endif
  end subroutine

  real(8) function relative_density(comm,a,b)
    integer,intent(in)::comm
    real(8),intent(in)::a(:),b(:)
    real(8)::local_values(2),global_values(2),local_scale(2),global_scale(2),scale_ratio,factor
    integer::ierr
    local_scale=0d0
    if(size(a)>0)local_scale=[maxval(abs(a-b)),maxval(abs(a))]
    call MPI_Allreduce(local_scale,global_scale,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(global_scale(1)==0d0)then;relative_density=0d0;return;endif
    if(global_scale(2)==0d0)then;relative_density=huge(1d0);return;endif
    local_values=[sum(((a-b)/global_scale(1))**2),sum((a/global_scale(2))**2)]
    call MPI_Allreduce(local_values,global_values,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    factor=sqrt(global_values(1)/global_values(2))
    if(global_scale(2)<1d0)then
      if(global_scale(1)>huge(1d0)*global_scale(2))then
        relative_density=huge(1d0);return
      endif
    endif
    scale_ratio=global_scale(1)/global_scale(2)
    if(factor>1d0)then
      if(scale_ratio>huge(1d0)/factor)then
        relative_density=huge(1d0)
      else
        relative_density=scale_ratio*factor
      endif
    else
      relative_density=scale_ratio*factor
    endif
  end function

  logical function finite_complex(values)
    complex(8),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function

  subroutine collective_step_status(comm,local_ok,local_message,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::local_ok
    character(*),intent(in)::local_message
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::bad,global_bad,ierr
    bad=merge(0,1,local_ok)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=global_bad==0
    if(ok)then
      message=''
    else if(len_trim(local_message)>0)then
      message=local_message
    else
      message='collective overlapping-Wannier SCF step failed'
    endif
  end subroutine
#endif
end module
