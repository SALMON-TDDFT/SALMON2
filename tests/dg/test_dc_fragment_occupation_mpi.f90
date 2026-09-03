#include "config.h"
program test_dc_fragment_occupation_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use occupation_kernel,only:solve_weighted_state_occupations
  use dc_fragment_occupation,only:determine_dc_fragment_occupations,assess_dc_fragment_occupation_capacity
  implicit none
  integer::comm,rank,nproc,ierr,fragment
  integer(int64)::fingerprint,minimum_fingerprint,maximum_fingerprint
  real(real64),allocatable::energies(:,:),core_norms(:,:),occupations(:,:),flat_occupations(:)
  real(real64)::chemical_potential,electron_count,direct_mu,direct_count
  real(real64),parameter::tolerance=1d-8
  real(real64)::flat_energies(4),flat_weights(4)
  logical,allocatable::representative_mask(:)
  logical::ok,capacity_sufficient,needs_extension(2),can_extend(2)
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  allocate(energies(2,2),core_norms(2,2),representative_mask(2))
  call test_unordered_spectrum()

  call distribute_two_fragment_case(&
    reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
    reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(ok,'gapped zero-temperature fragment occupation failed: '//trim(message))
  call require(all(abs(occupations-reshape([2d0,0d0,2d0,0d0],[2,2]))<1d-14),&
    'gapped zero-temperature fragment occupations changed')
  call require(abs(electron_count-2d0)<tolerance,&
    'gapped zero-temperature weighted electron count changed')
  call require(chemical_potential>=-0.5d0.and.chemical_potential<0.5d0,&
    'gapped fragments did not share one chemical potential')
  fingerprint=1469598103934665603_int64
  call hash_real(fingerprint,chemical_potential);call hash_real(fingerprint,electron_count)
  do fragment=1,size(occupations,2)
    call hash_vector(fingerprint,occupations(:,fragment))
  enddo
  call MPI_Allreduce(fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
  call MPI_Allreduce(fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  call require(minimum_fingerprint==maximum_fingerprint,&
    'fragment occupation fingerprint is rank-dependent')

  call distribute_two_fragment_case(&
    reshape([-1d0,3.5d0,-0.5d0,4d0],[2,2]),&
    reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0.1d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(ok,'finite-temperature fragment occupation failed: '//trim(message))
  call require(all(occupations>=0d0).and.all(occupations<=2d0),&
    'finite-temperature occupations exceed spin bounds')
  call require(occupations(2,1)>0d0.and.occupations(2,2)>0d0.and.&
    maxval(occupations(2,:))<=tolerance,&
    'explicit finite-temperature guard-state tail was discarded')
  call require(abs(sum(occupations*reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))-2d0)<tolerance,&
    'finite-temperature core-weighted electron count changed')
  call hash_real(fingerprint,chemical_potential);call hash_real(fingerprint,electron_count)
  do fragment=1,size(occupations,2)
    call hash_vector(fingerprint,occupations(:,fragment))
  enddo
  flat_energies=[-1d0,3.5d0,-0.5d0,4d0]
  flat_weights=[0.75d0,0.25d0,0.25d0,0.75d0]
  call solve_weighted_state_occupations(flat_energies,flat_weights,2d0,0.1d0,2d0,&
    flat_occupations,direct_mu,direct_count,ok,message)
  call require(ok,'weighted-state occupation kernel failed: '//trim(message))
  call require(maxval(abs(reshape(occupations,[4])-flat_occupations))<1d-14.and.&
    abs(chemical_potential-direct_mu)<1d-14.and.abs(electron_count-direct_count)<1d-14,&
    'distributed adapter diverged from the authoritative weighted kernel')

  call distribute_two_fragment_case(&
    reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
    reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0.1d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(.not.ok.and.index(message,'tail')>0,&
    'finite-temperature solve accepted an occupied spectrum boundary')
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0.1d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message,needs_extension)
  call require(ok.and.all(needs_extension),'finite-temperature tail did not request per-fragment extension')

  deallocate(energies,core_norms);allocate(energies(3,2),core_norms(3,2))
  call distribute_two_fragment_case(&
    reshape([-1d0,3.28d0,3.28d0,-0.5d0,3.28d0,3.28d0],[3,2]),&
    reshape([0.75d0,0.75d0,0.75d0,0.25d0,0.75d0,0.75d0],[3,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0.1d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(.not.ok.and.index(message,'tail')>0,&
    'finite-temperature solve ignored a degenerate boundary-tail sum')
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0.1d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message,needs_extension)
  call require(ok.and.all(needs_extension),'degenerate terminal shell charge was not summed')

  deallocate(energies,core_norms);allocate(energies(2,2),core_norms(2,2))
  call distribute_two_fragment_case(&
    reshape([-1d0,0d0,0d0,1d0],[2,2]),&
    reshape([0.5d0,0.75d0,0.25d0,0.5d0],[2,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,3d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(ok,'Fermi-edge degeneracy solve failed: '//trim(message))
  call require(all(abs(occupations-reshape([2d0,2d0,2d0,0d0],[2,2]))<1d-14).and.&
    chemical_potential==0d0.and.abs(electron_count-3d0)<tolerance,&
    'filled Fermi-degenerate shell convention changed')
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(ok,'partial Fermi-edge degeneracy solve failed: '//trim(message))
  call require(all(abs(occupations-reshape([2d0,1d0,1d0,0d0],[2,2]))<1d-14).and.&
    chemical_potential==0d0.and.abs(electron_count-2d0)<tolerance,&
    'partial Fermi-degenerate shell was not occupied uniformly')
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message,needs_extension)
  call require(ok.and.needs_extension(1).and..not.needs_extension(2),&
    'zero-temperature incomplete Fermi terminal shell was not extended')
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message,needs_extension,&
    terminal_shell_complete=[.true.,.true.])
  call require(ok.and..not.any(needs_extension),'explicitly complete zero-temperature shell extended unnecessarily')
  flat_energies=[-1d0,-0.5d0,0d0,1d0]
  flat_weights=[0.7d0,0.2d0,0d0,0d0]
  call solve_weighted_state_occupations(flat_energies,flat_weights,1d0,0d0,2d0,&
    flat_occupations,direct_mu,direct_count,ok,message)
  call require(ok.and.abs(flat_occupations(1)-1d0/0.7d0)<1d-14.and.&
    all(abs(flat_occupations(2:))<1d-14).and.direct_mu==-1d0.and.&
    abs(direct_count-1d0)<tolerance,&
    'zero-temperature fractional core norm could not carry the electron target')
  flat_energies=[0d0,-1d0,1d0,1d-14]
  flat_weights=[0.75d0,0.5d0,0.5d0,0.25d0]
  call solve_weighted_state_occupations(flat_energies,flat_weights,2d0,0d0,2d0,&
    flat_occupations,direct_mu,direct_count,ok,message)
  call require(ok.and.all(abs(flat_occupations-[1d0,2d0,0d0,1d0])<1d-14).and.&
    direct_mu==0d0.and.abs(direct_count-2d0)<tolerance,&
    'unsorted weighted states changed the uniform Fermi-shell ensemble')

  deallocate(energies,core_norms);allocate(energies(3,2),core_norms(3,2))
  call distribute_two_fragment_case(&
    reshape([-1d0,1.5d0,1.5d0,-0.5d0,1d0,2d0],[3,2]),&
    reshape([0.75d0,0.25d0,0d0,0.25d0,0.5d0,0.25d0],[3,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(ok,'zero-weight state padding failed: '//trim(message))
  call require(all(abs(occupations(:,1)-[2d0,0d0,0d0])<1d-14).and.&
    all(abs(occupations(:,2)-[2d0,0d0,0d0])<1d-14),&
    'different fragment state counts changed physical occupations')

  deallocate(energies,core_norms);allocate(energies(1,2),core_norms(1,2))
  call distribute_two_fragment_case(reshape([-1d0,-0.5d0],[1,2]),&
    reshape([0.75d0,0.5d0],[1,2]))
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2.6d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(.not.ok.and.index(message,'capacity')>0,&
    'collective fragment solve accepted insufficient spectrum capacity')
  can_extend=[.true.,.true.]
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2.6d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(ok.and..not.capacity_sufficient.and.all(needs_extension),&
    'capacity preflight must extend every nonexhausted fragment before occupation solve')
  can_extend=[.false.,.true.]
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2.6d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(ok.and..not.needs_extension(1).and.needs_extension(2),'exhausted fragment marked for capacity extension')
  can_extend=.false.
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2.6d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(.not.ok.and.index(message,'insufficient')>0,'exhausted capacity must fail collectively')
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(ok.and.capacity_sufficient.and..not.any(needs_extension),'sufficient capacity triggered growth')
  if(representative_mask(1))core_norms(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(.not.ok,'nonfinite capacity input accepted')
  call distribute_two_fragment_case(reshape([-1d0,-0.5d0],[1,2]),reshape([0.75d0,0.5d0],[1,2]))
  if(representative_mask(1))core_norms(1,1)=-0.75d0
  call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
    2d0,2d0,tolerance,capacity_sufficient,needs_extension,ok,message)
  call require(.not.ok,'negative capacity input accepted')
  call distribute_two_fragment_case(reshape([-1d0,-0.5d0],[1,2]),reshape([0.75d0,0.5d0],[1,2]))
  if(nproc>1)then
    call assess_dc_fragment_occupation_capacity(comm,core_norms,representative_mask,can_extend,&
      merge(1d0,2d0,rank==0),2d0,tolerance,capacity_sufficient,needs_extension,ok,message)
    call require(.not.ok,'rank-disagreeing capacity controls accepted')
  endif

  deallocate(energies,core_norms);allocate(energies(2,2),core_norms(2,2))
  call distribute_two_fragment_case(&
    reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
    reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
  representative_mask(2)=.false.
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(.not.ok.and.index(message,'representative')>0,&
    'collective fragment solve accepted a missing representative')

  if(nproc>1)then
    call distribute_two_fragment_case(&
      reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
      reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
    representative_mask(1)=.true.
    call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
      0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
    call require(.not.ok.and.index(message,'representative')>0,&
      'collective fragment solve accepted duplicate representatives')

    call distribute_two_fragment_case(&
      reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
      reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
    call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
      merge(0d0,0.1d0,rank==0),2d0,2d0,tolerance,chemical_potential,&
      occupations,electron_count,ok,message)
    call require(.not.ok.and.index(message,'rank')>0,&
      'collective fragment solve accepted rank-dependent temperature')
  endif

  call distribute_two_fragment_case(&
    reshape([-1d0,0.5d0,-0.5d0,1d0],[2,2]),&
    reshape([0.75d0,0.25d0,0.25d0,0.75d0],[2,2]))
  if(representative_mask(1))core_norms(1,1)=-0.75d0
  call determine_dc_fragment_occupations(comm,energies,core_norms,representative_mask,&
    0d0,2d0,2d0,tolerance,chemical_potential,occupations,electron_count,ok,message)
  call require(.not.ok,'collective fragment solve accepted a negative core norm')

  if(rank==0)then
    write(*,'(a,i0,a,i0)')'DC_FRAGMENT_OCCUPATION ranks=',nproc,' fingerprint=',fingerprint
    write(*,'(a,i0,a)')'PASS DC fragment occupation on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine test_unordered_spectrum()
    real(real64)::e(4,2),w(4,2),saved_e(4,2),saved_w(4,2),expected(4,2)
    real(real64),allocatable::answer(:,:),direct(:)
    real(real64)::reference_mu,reference_electrons
    logical::representatives(2),tail(2)
    integer::f
    e=0d0;w=0d0
    representatives=[rank==0,rank==mod(1,nproc)]
    do f=1,2
      if(.not.representatives(f))cycle
      ! Unequal weights expose energy-only sorting; zero-weight padding is not a state.
      e(:,f)=[4d0,-1d0,0d0,0d0]
      w(:,f)=[0.25d0,0.5d0,0.25d0,0d0]
    enddo
    saved_e=e;saved_w=w
    call determine_dc_fragment_occupations(comm,e,w,representatives,0d0,2d0,2d0,tolerance,&
      chemical_potential,answer,electron_count,ok,message)
    call require(.not.ok,'default sorted-spectrum contract changed')
    call determine_dc_fragment_occupations(comm,e,w,representatives,0d0,2d0,2d0,tolerance,&
      chemical_potential,answer,electron_count,ok,message,tail,allow_unordered=.true.)
    call require(ok,'unordered measured spectrum rejected: '//trim(message))
    expected=0d0;expected(2,:)=2d0
    call require(maxval(abs(answer-expected))<1d-14.and..not.any(tail),&
      'occupations were not returned in coefficient-column order')
    call require(all(e==saved_e).and.all(w==saved_w),'occupation packing mutated caller data')
    call require(abs(electron_count-2d0)<tolerance,'unordered spectrum changed electron count')
    call determine_dc_fragment_occupations(comm,e,w,representatives,0.1d0,2d0,2d0,tolerance,&
      chemical_potential,answer,electron_count,ok,message,tail,allow_unordered=.true.)
    call require(ok.and..not.any(tail),'unordered finite-temperature tail failed: '//trim(message))
    call solve_weighted_state_occupations([4d0,-1d0,0d0,0d0,4d0,-1d0,0d0,0d0],&
      [0.25d0,0.5d0,0.25d0,0d0,0.25d0,0.5d0,0.25d0,0d0],2d0,0.1d0,2d0,&
      direct,reference_mu,reference_electrons,ok,message)
    call require(ok.and.maxval(abs(answer-reshape(direct,[4,2])))<1d-12.and.&
      abs(chemical_potential-reference_mu)<1d-12,'unordered adapter changed weighted-kernel physics')
    do f=1,2
      if(representatives(f))then
        e(:,f)=[0d0,-1d0,0d0,0d0]
        w(:,f)=[0.25d0,0.5d0,0.25d0,0d0]
      endif
    enddo
    call determine_dc_fragment_occupations(comm,e,w,representatives,0d0,2d0,3d0,tolerance,&
      chemical_potential,answer,electron_count,ok,message,tail,allow_unordered=.true.)
    ! Zero-weight padding shares the shell's formal occupation but carries no charge.
    expected=1d0;expected(2,:)=2d0
    call require(ok.and.maxval(abs(answer-expected))<1d-14.and.all(tail),&
      'unordered degenerate terminal shell lost its occupations or extension mask')
    if(nproc>1)then
      call determine_dc_fragment_occupations(comm,e,w,representatives,0d0,2d0,3d0,tolerance,&
        chemical_potential,answer,electron_count,ok,message,tail,allow_unordered=rank==0)
      call require(.not.ok,'rank-disagreeing unordered-spectrum policy accepted')
    endif
    if(representatives(1))e(1,1)=ieee_value(0d0,ieee_quiet_nan)
    call determine_dc_fragment_occupations(comm,e,w,representatives,0d0,2d0,3d0,tolerance,&
      chemical_potential,answer,electron_count,ok,message,tail,allow_unordered=.true.)
    call require(.not.ok.and..not.allocated(answer),'unordered option bypassed finite-spectrum validation')
  end subroutine
  subroutine distribute_two_fragment_case(global_energies,global_core_norms)
    real(real64),intent(in)::global_energies(:,:),global_core_norms(:,:)
    energies=100d0+real(rank,real64)
    core_norms=0.125d0*real(rank+1,real64)
    representative_mask=.false.
    do fragment=1,size(global_energies,2)
      if(rank==mod(fragment-1,nproc))then
        representative_mask(fragment)=.true.
        energies(:,fragment)=global_energies(:,fragment)
        core_norms(:,fragment)=global_core_norms(:,fragment)
      endif
    enddo
  end subroutine distribute_two_fragment_case

  subroutine hash_real(hash,value)
    integer(int64),intent(inout)::hash
    real(real64),intent(in)::value
    hash=ieor(ishftc(hash,7),transfer(value,0_int64))
  end subroutine hash_real

  subroutine hash_vector(hash,values)
    integer(int64),intent(inout)::hash
    real(real64),intent(in)::values(:)
    integer::i
    do i=1,size(values);call hash_real(hash,values(i));enddo
  end subroutine hash_vector

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      call MPI_Abort(comm,1,ierr)
    endif
  end subroutine require
end program test_dc_fragment_occupation_mpi
