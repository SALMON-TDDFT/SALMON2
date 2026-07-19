module dg_wpw_lcfo_potential_map
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private
  type,public::s_dg_wpw_operator_candidate
    complex(8),allocatable::ww(:),wp(:),pp(:)
    integer(int64)::fingerprint=0_int64
    logical::valid=.false.
  end type
  type,public::s_dg_wpw_grid_route
    integer::comm=MPI_COMM_NULL,nlocal=0
    integer,allocatable::send_counts(:),recv_counts(:),send_displs(:),recv_displs(:)
    integer,allocatable::send_original(:),recv_ids(:)
    logical::valid=.false.
  end type
  type,public::s_dg_wpw_map_candidate
    real(8),allocatable::rho_in(:),rho_raw(:),rho_candidate(:),potential_candidate(:)
    real(8)::energy_candidate=0d0
    integer(int64)::operator_candidate=0_int64
    complex(8),allocatable::operator_ww(:),operator_wp(:),operator_pp(:)
    integer::publication_epoch=0,potential_epoch=0,energy_epoch=0,operator_epoch=0
    logical::density_mix_only=.true.,valid=.false.
  end type
  type,public::s_dg_wpw_potential_map_state
    integer::comm=MPI_COMM_NULL,publication_epoch=0,iteration=0
    real(8),allocatable::rho(:),potential(:),mixer_history(:),residual_history(:)
    real(8)::energy=0d0
    integer(int64)::operator_fingerprint=0_int64
    complex(8),allocatable::operator_ww(:),operator_wp(:),operator_pp(:)
  end type
  public::initialize_dg_wpw_potential_map_state,run_dg_wpw_lcfo_potential_map
  public::rollback_dg_wpw_lcfo_potential_map,publish_dg_wpw_lcfo_potential_map
  public::build_dg_wpw_raw_density,redistribute_dg_wpw_core_density
  public::build_dg_wpw_core_density
  public::prepare_dg_wpw_core_density_route,return_dg_wpw_core_potential
  public::evaluate_and_run_dg_wpw_lcfo_potential_map
  abstract interface
    subroutine dg_wpw_potential_energy_evaluator(rho,epoch,potential,energy,info)
      real(8),intent(in)::rho(:);integer,intent(in)::epoch
      real(8),intent(out)::potential(:),energy;integer,intent(out)::info
    end subroutine
    subroutine dg_wpw_operator_rebuilder(potential,epoch,candidate,info)
      import s_dg_wpw_operator_candidate
      real(8),intent(in)::potential(:);integer,intent(in)::epoch
      type(s_dg_wpw_operator_candidate),intent(out)::candidate;integer,intent(out)::info
    end subroutine
  end interface
contains
  subroutine evaluate_and_run_dg_wpw_lcfo_potential_map(state,input_epoch,rho_raw,mixing,evaluator,&
      rebuilder,info)
    type(s_dg_wpw_potential_map_state),intent(inout)::state
    integer,intent(in)::input_epoch
    real(8),intent(in)::rho_raw(:),mixing
    procedure(dg_wpw_potential_energy_evaluator)::evaluator
    procedure(dg_wpw_operator_rebuilder)::rebuilder
    integer,intent(out)::info
    real(8),allocatable::rho_candidate(:),potential(:)
    real(8)::energy
    type(s_dg_wpw_operator_candidate)::operator_candidate
    integer::stage_info,local_bad,global_bad,ierr
    info=1
    local_bad=merge(0,1,input_epoch==state%publication_epoch.and.size(rho_raw)==size(state%rho).and.&
      ieee_is_finite(mixing).and.mixing>=0d0.and.mixing<=1d0.and.all(ieee_is_finite(rho_raw)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    rho_candidate=(1d0-mixing)*state%rho+mixing*rho_raw;allocate(potential(size(state%potential)))
    call evaluator(rho_candidate,input_epoch,potential,energy,stage_info)
    local_bad=merge(0,1,stage_info==0.and.all(ieee_is_finite(potential)).and.ieee_is_finite(energy))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call rebuilder(potential,input_epoch,operator_candidate,stage_info)
    local_bad=0
    if(stage_info/=0.or..not.operator_candidate%valid.or.operator_candidate%fingerprint==0_int64.or.&
       .not.allocated(operator_candidate%ww).or..not.allocated(operator_candidate%wp).or.&
       .not.allocated(operator_candidate%pp))local_bad=1
    if(local_bad==0)then
      if(.not.finite_complex(operator_candidate%ww).or..not.finite_complex(operator_candidate%wp).or.&
         .not.finite_complex(operator_candidate%pp))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call run_dg_wpw_lcfo_potential_map(state,input_epoch,rho_candidate,1d0,potential,energy,&
      operator_candidate%fingerprint,0,info,potential_epoch=input_epoch,energy_epoch=input_epoch,&
      operator_epoch=input_epoch,operator_ww=operator_candidate%ww,operator_wp=operator_candidate%wp,&
      operator_pp=operator_candidate%pp)
  end subroutine
  subroutine build_dg_wpw_raw_density(windows,coefficients,occupations,weights,rho_raw,charge,info)
    complex(8),intent(in)::windows(:,:),coefficients(:,:)
    real(8),intent(in)::occupations(:),weights(:)
    real(8),intent(out)::rho_raw(:),charge
    integer,intent(out)::info
    complex(8),allocatable::psi(:,:)
    integer::i
    info=1;rho_raw=0d0;charge=0d0
    if(size(windows,1)/=size(weights).or.size(rho_raw)/=size(weights).or.&
       size(windows,2)/=size(coefficients,1).or.size(coefficients,2)/=size(occupations))return
    if(any(weights<0d0).or.any(occupations<0d0).or..not.all(ieee_is_finite(weights)).or.&
       .not.all(ieee_is_finite(occupations)).or..not.finite_complex_2(windows).or.&
       .not.finite_complex_2(coefficients))return
    psi=matmul(windows,coefficients)
    do i=1,size(occupations)
      rho_raw=rho_raw+occupations(i)*real(conjg(psi(:,i))*psi(:,i),8)
    enddo
    charge=sum(weights*rho_raw)
    if(.not.all(ieee_is_finite(rho_raw)).or..not.ieee_is_finite(charge))return
    info=0
  end subroutine
  subroutine build_dg_wpw_core_density(w_points,p_points,cw,cp,occupations,weights,rho_raw,charge,info)
    complex(8),intent(in)::w_points(:,:),p_points(:,:),cw(:,:),cp(:,:)
    real(8),intent(in)::occupations(:),weights(:)
    real(8),intent(out)::rho_raw(:),charge
    integer,intent(out)::info
    complex(8),allocatable::psi(:,:)
    integer::iorb
    info=1;rho_raw=0d0;charge=0d0
    if(size(w_points,2)/=size(p_points,2).or.size(rho_raw)/=size(weights).or.&
       size(w_points,2)/=size(weights).or.size(w_points,1)/=size(cw,1).or.&
       size(p_points,1)/=size(cp,1).or.size(cw,2)/=size(occupations).or.&
       size(cp,2)/=size(occupations).or.any(weights<0d0).or.any(occupations<0d0).or.&
       .not.finite_complex(reshape(w_points,[size(w_points)] )).or.&
       .not.finite_complex(reshape(p_points,[size(p_points)] )).or.&
       .not.finite_complex(reshape(cw,[size(cw)] )).or.&
       .not.finite_complex(reshape(cp,[size(cp)] )).or.&
       .not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(occupations)))return
    psi=matmul(transpose(w_points),cw)+matmul(transpose(p_points),cp)
    do iorb=1,size(occupations)
      rho_raw=rho_raw+occupations(iorb)*abs(psi(:,iorb))**2
    enddo
    charge=sum(weights*rho_raw)
    if(.not.all(ieee_is_finite(rho_raw)).or..not.ieee_is_finite(charge))return
    info=0
  end subroutine

  subroutine redistribute_dg_wpw_core_density(comm,grid_ids,rho_local,weights,nglobal,owned_ids,rho_owned,&
      charge_global,info)
    integer,intent(in)::comm,grid_ids(:),nglobal,owned_ids(:)
    real(8),intent(in)::rho_local(:),weights(:)
    real(8),intent(out)::rho_owned(:),charge_global
    integer,intent(out)::info
    type(s_dg_wpw_grid_route)::route
    call prepare_dg_wpw_core_density_route(comm,grid_ids,rho_local,weights,nglobal,owned_ids,route,&
      rho_owned,charge_global,info)
  end subroutine

  subroutine prepare_dg_wpw_core_density_route(comm,grid_ids,rho_local,weights,nglobal,owned_ids,route,&
      rho_owned,charge_global,info,destinations)
    integer,intent(in)::comm,grid_ids(:),nglobal,owned_ids(:)
    real(8),intent(in)::rho_local(:),weights(:)
    type(s_dg_wpw_grid_route),intent(out)::route
    real(8),intent(out)::rho_owned(:),charge_global
    integer,intent(out)::info
    integer,intent(in),optional::destinations(:)
    integer::nproc,rank,nlocal,ierr,i,p,total,local_bad,global_bad,pos,expected
    integer,allocatable::send_counts(:),recv_counts(:),send_displs(:),recv_displs(:),cursor(:),send_ids(:),recv_ids(:)
    integer,allocatable::send_original(:)
    real(8),allocatable::send_rho(:),recv_rho(:),send_weight(:),recv_weight(:)
    integer,allocatable::multiplicity(:)
    real(8)::charge_local
    info=1;rho_owned=0d0;charge_global=0d0;nlocal=size(grid_ids)
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    expected=(nglobal+nproc-1-rank)/nproc
    local_bad=merge(0,1,nlocal==size(rho_local).and.nlocal==size(weights).and.nglobal>0.and.&
      size(rho_owned)==size(owned_ids).and.&
      all(grid_ids>=1).and.all(grid_ids<=nglobal).and.all(owned_ids>=1).and.all(owned_ids<=nglobal).and.&
      all(weights>=0d0).and.all(ieee_is_finite(rho_local)).and.all(ieee_is_finite(weights)))
    if(present(destinations))then
      if(size(destinations)/=nlocal)then
        local_bad=1
      else if(any(destinations<0).or.any(destinations>=nproc))then
        local_bad=1
      endif
    else if(size(owned_ids)/=expected.or.any(owned_ids/=[(rank+1+(i-1)*nproc,i=1,expected)]))then
      local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(send_counts(nproc),recv_counts(nproc),send_displs(nproc),recv_displs(nproc),cursor(nproc))
    send_counts=0
    do i=1,nlocal
      if(present(destinations))then;p=destinations(i)+1;else;p=modulo(grid_ids(i)-1,nproc)+1;endif
      send_counts(p)=send_counts(p)+1
    enddo
    call MPI_Alltoall(send_counts,1,MPI_INTEGER,recv_counts,1,MPI_INTEGER,comm,ierr);if(ierr/=MPI_SUCCESS)return
    send_displs(1)=0;recv_displs(1)=0
    do p=2,nproc
      send_displs(p)=send_displs(p-1)+send_counts(p-1);recv_displs(p)=recv_displs(p-1)+recv_counts(p-1)
    enddo
    total=sum(recv_counts);cursor=send_displs
    allocate(send_ids(nlocal),send_original(nlocal),send_rho(nlocal),send_weight(nlocal),&
      recv_ids(total),recv_rho(total),recv_weight(total))
    do i=1,nlocal
      if(present(destinations))then;p=destinations(i)+1;else;p=modulo(grid_ids(i)-1,nproc)+1;endif
      pos=cursor(p)+1;cursor(p)=pos
      send_ids(pos)=grid_ids(i);send_original(pos)=i;send_rho(pos)=rho_local(i);send_weight(pos)=weights(i)
    enddo
    call MPI_Alltoallv(send_ids,send_counts,send_displs,MPI_INTEGER,recv_ids,recv_counts,recv_displs,MPI_INTEGER,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(send_rho,send_counts,send_displs,MPI_DOUBLE_PRECISION,&
      recv_rho,recv_counts,recv_displs,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Alltoallv(send_weight,send_counts,send_displs,MPI_DOUBLE_PRECISION,&
      recv_weight,recv_counts,recv_displs,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    allocate(multiplicity(size(owned_ids)));multiplicity=0
    do i=1,total
      pos=find_integer(owned_ids,recv_ids(i))
      if(pos==0)return
      multiplicity(pos)=multiplicity(pos)+1;rho_owned(pos)=recv_rho(i)
    enddo
    if(any(multiplicity/=1))return
    charge_local=sum(recv_weight*recv_rho)
    call MPI_Allreduce(charge_local,charge_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    if(.not.ieee_is_finite(charge_global))return
    route%comm=comm;route%nlocal=nlocal;route%send_counts=send_counts;route%recv_counts=recv_counts
    route%send_displs=send_displs;route%recv_displs=recv_displs
    route%send_original=send_original;route%recv_ids=recv_ids;route%valid=.true.
    info=0
  end subroutine

  subroutine return_dg_wpw_core_potential(route,owned_ids,potential_owned,potential_core,info)
    type(s_dg_wpw_grid_route),intent(in)::route
    integer,intent(in)::owned_ids(:)
    real(8),intent(in)::potential_owned(:)
    real(8),intent(out)::potential_core(:)
    integer,intent(out)::info
    real(8),allocatable::reply_send(:),reply_recv(:)
    integer::i,pos,ierr,local_bad,global_bad
    info=1;potential_core=0d0;local_bad=0
    if(.not.route%valid.or.size(owned_ids)/=size(potential_owned).or.&
       size(potential_core)/=route%nlocal.or..not.all(ieee_is_finite(potential_owned)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,route%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(reply_send(size(route%recv_ids)),reply_recv(route%nlocal));reply_send=0d0;reply_recv=0d0
    do i=1,size(route%recv_ids)
      pos=find_integer(owned_ids,route%recv_ids(i));if(pos==0)then;local_bad=1;exit;endif
      reply_send(i)=potential_owned(pos)
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,route%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Alltoallv(reply_send,route%recv_counts,route%recv_displs,MPI_DOUBLE_PRECISION,&
      reply_recv,route%send_counts,route%send_displs,MPI_DOUBLE_PRECISION,route%comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    do i=1,route%nlocal;potential_core(route%send_original(i))=reply_recv(i);enddo
    if(.not.all(ieee_is_finite(potential_core)))return
    info=0
  end subroutine

  integer function find_integer(values,target)result(position)
    integer,intent(in)::values(:),target
    integer::i
    position=0
    do i=1,size(values);if(values(i)==target)then;position=i;return;endif;enddo
  end function

  logical function finite_complex_2(a)result(ok)
    complex(8),intent(in)::a(:,:)
    ok=all(ieee_is_finite(real(a,8))).and.all(ieee_is_finite(aimag(a)))
  end function
  subroutine initialize_dg_wpw_potential_map_state(state,comm,epoch,rho,potential,energy,operator_fingerprint,&
      mixer_history,residual_history,iteration,info)
    type(s_dg_wpw_potential_map_state),intent(out)::state
    integer,intent(in)::comm,epoch,iteration
    real(8),intent(in)::rho(:),potential(:),energy,mixer_history(:),residual_history(:)
    integer(int64),intent(in)::operator_fingerprint
    integer,intent(out)::info
    info=1
    if(comm==MPI_COMM_NULL.or.epoch<0.or.iteration<0.or.size(rho)/=size(potential))return
    if(.not.all(ieee_is_finite(rho)).or..not.all(ieee_is_finite(potential)).or.&
       .not.ieee_is_finite(energy).or..not.all(ieee_is_finite(mixer_history)).or.&
       .not.all(ieee_is_finite(residual_history)))return
    state%comm=comm;state%publication_epoch=epoch;state%iteration=iteration
    state%rho=rho;state%potential=potential;state%energy=energy
    state%operator_fingerprint=operator_fingerprint
    allocate(state%operator_ww(0),state%operator_wp(0),state%operator_pp(0))
    state%mixer_history=mixer_history;state%residual_history=residual_history
    info=0
  end subroutine

  subroutine run_dg_wpw_lcfo_potential_map(state,input_epoch,rho_raw,mixing,potential_candidate,&
      energy_candidate,operator_candidate,failure_stage,info,potential_epoch,energy_epoch,operator_epoch,&
      operator_ww,operator_wp,operator_pp)
    type(s_dg_wpw_potential_map_state),intent(inout)::state
    integer,intent(in)::input_epoch,failure_stage
    real(8),intent(in)::rho_raw(:),mixing,potential_candidate(:),energy_candidate
    integer(int64),intent(in)::operator_candidate
    integer,intent(out)::info
    integer,intent(in),optional::potential_epoch,energy_epoch,operator_epoch
    complex(8),intent(in),optional::operator_ww(:),operator_wp(:),operator_pp(:)
    type(s_dg_wpw_map_candidate)::candidate
    integer::local_bad,global_bad,ierr,vepoch,eepoch,oepoch
    info=1
    vepoch=input_epoch;eepoch=input_epoch;oepoch=input_epoch
    if(present(potential_epoch))vepoch=potential_epoch
    if(present(energy_epoch))eepoch=energy_epoch
    if(present(operator_epoch))oepoch=operator_epoch
    local_bad=0
    if(state%comm==MPI_COMM_NULL.or.input_epoch/=state%publication_epoch.or.&
       size(rho_raw)/=size(state%rho).or.size(potential_candidate)/=size(state%potential).or.&
       .not.ieee_is_finite(mixing).or.mixing<0d0.or.mixing>1d0.or.&
       .not.all(ieee_is_finite(rho_raw)).or..not.all(ieee_is_finite(potential_candidate)).or.&
       .not.ieee_is_finite(energy_candidate).or.failure_stage/=0.or.vepoch/=input_epoch.or.&
       eepoch/=input_epoch.or.oepoch/=input_epoch)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,state%comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call rollback_dg_wpw_lcfo_potential_map(candidate)
      return
    endif
    candidate%rho_in=state%rho;candidate%rho_raw=rho_raw
    candidate%rho_candidate=(1d0-mixing)*state%rho+mixing*rho_raw
    candidate%potential_candidate=potential_candidate
    candidate%energy_candidate=energy_candidate;candidate%operator_candidate=operator_candidate
    if(present(operator_ww))candidate%operator_ww=operator_ww
    if(present(operator_wp))candidate%operator_wp=operator_wp
    if(present(operator_pp))candidate%operator_pp=operator_pp
    if(.not.allocated(candidate%operator_ww))allocate(candidate%operator_ww(0))
    if(.not.allocated(candidate%operator_wp))allocate(candidate%operator_wp(0))
    if(.not.allocated(candidate%operator_pp))allocate(candidate%operator_pp(0))
    candidate%potential_epoch=vepoch;candidate%energy_epoch=eepoch;candidate%operator_epoch=oepoch
    candidate%publication_epoch=input_epoch+1;candidate%valid=.true.
    call publish_dg_wpw_lcfo_potential_map(state,candidate,info)
  end subroutine

  subroutine rollback_dg_wpw_lcfo_potential_map(candidate)
    type(s_dg_wpw_map_candidate),intent(inout)::candidate
    candidate%valid=.false.
  end subroutine

  subroutine publish_dg_wpw_lcfo_potential_map(state,candidate,info)
    type(s_dg_wpw_potential_map_state),intent(inout)::state
    type(s_dg_wpw_map_candidate),intent(in)::candidate
    integer,intent(out)::info
    info=1
    if(.not.candidate%valid.or.candidate%publication_epoch/=state%publication_epoch+1.or.&
       candidate%potential_epoch/=state%publication_epoch.or.&
       candidate%energy_epoch/=state%publication_epoch.or.candidate%operator_epoch/=state%publication_epoch)return
    state%rho=candidate%rho_candidate;state%potential=candidate%potential_candidate
    state%energy=candidate%energy_candidate;state%operator_fingerprint=candidate%operator_candidate
    state%operator_ww=candidate%operator_ww;state%operator_wp=candidate%operator_wp
    state%operator_pp=candidate%operator_pp
    state%publication_epoch=candidate%publication_epoch;state%iteration=state%iteration+1
    state%mixer_history=state%rho;state%residual_history=candidate%rho_raw-candidate%rho_in
    info=0
  end subroutine
  logical function finite_complex(values)result(ok)
    complex(8),intent(in)::values(:)
    ok=all(ieee_is_finite(real(values,8))).and.all(ieee_is_finite(aimag(values)))
  end function
end module
