#include "config.h"
program test_rt_dg_overlapping_wannier_mpi
  use mpi
  use iso_fortran_env,only:int64,real64
  use ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use rt_dg_overlapping_wannier,only:s_dg_overlapping_wannier_rt_state,&
    initialize_dg_overlapping_wannier_rt,advance_dg_overlapping_wannier_rt,&
    write_dg_overlapping_wannier_rt_restart,read_dg_overlapping_wannier_rt_restart,&
    evaluate_dg_overlapping_wannier_observables,&
    write_dg_overlapping_wannier_rt_observable_sample
  implicit none
  integer::comm,rank,nproc,ierr,i,nlocal,unit
  integer,allocatable::row_ids(:)
  complex(real64),allocatable::srows(:,:),hrows(:,:),xrows(:,:,:),vrows(:,:,:)
  complex(real64)::initial(2,1),continuous(2,1),split(2,1),zero_ref(2,1),field_ref(2,1),&
    velocity_ref(2,1),phase_coeff(2,1),&
    rotated_coeff(2,1),metric(2,2),hamiltonian(2,2),rotation(2,2),rotated_metric(2,2),&
    rotated_hamiltonian(2,2),full_position(3,2,2),full_velocity(3,2,2),&
    rotated_position(3,2,2),rotated_velocity(3,2,2),phase
  complex(real64),allocatable::rotated_srows(:,:),rotated_hrows(:,:),rotated_xrows(:,:,:),rotated_vrows(:,:,:)
  real(real64)::norm0,norm1,energy0,energy1
  type(s_dg_overlapping_wannier_rt_state)::state,split_state,restart_state,zero_state,field_state,&
    velocity_state,phase_state
  type(s_dg_overlapping_wannier_rt_state)::rotation_state
  logical::ok
  character(256)::message
  character(512)::restart_prefix
  character(512)::observable_prefix

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_environment_variable('OW_RT_RESTART_PREFIX',restart_prefix)
  call get_environment_variable('OW_RT_OBSERVABLE_PREFIX',observable_prefix)
  nlocal=count([(mod(i-1,nproc)==rank,i=1,2)])
  allocate(row_ids(nlocal),srows(nlocal,2),hrows(nlocal,2),xrows(3,nlocal,2),vrows(3,nlocal,2))
  nlocal=0
  do i=1,2
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;row_ids(nlocal)=i
  enddo
  metric=reshape([cmplx(2d0,0d0,8),cmplx(0.2d0,0d0,8),&
    cmplx(0.2d0,0d0,8),cmplx(1.5d0,0d0,8)],[2,2])
  hamiltonian=reshape([cmplx(0.3d0,0d0,8),cmplx(0.05d0,0d0,8),&
    cmplx(0.05d0,0d0,8),cmplx(0.8d0,0d0,8)],[2,2])
  do i=1,size(row_ids)
    srows(i,:)=metric(row_ids(i),:)
    hrows(i,:)=hamiltonian(row_ids(i),:)
    xrows(:,i,:)=(0d0,0d0);vrows(:,i,:)=(0d0,0d0)
    xrows(1,i,row_ids(i))=cmplx(real(row_ids(i)-1,real64),0d0,real64)
    vrows(1,i,3-row_ids(i))=cmplx(0d0,merge(0.1d0,-0.1d0,row_ids(i)==1),real64)
  enddo
  initial(:,1)=[cmplx(0.6d0,0.2d0,8),cmplx(-0.3d0,0.1d0,8)]
  initial=initial/sqrt(real(dot_product(initial(:,1),matmul(metric,initial(:,1))),real64))
  continuous=initial;split=initial
  norm0=real(dot_product(initial(:,1),matmul(metric,initial(:,1))),real64)
  energy0=real(dot_product(initial(:,1),matmul(hamiltonian,initial(:,1))),real64)

  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,continuous,state,ok,message)
  call require(ok,trim(message))
  call require((rank==0.and.size(state%metric)==4).or.(rank/=0.and.size(state%metric)==0),&
    'dense eigensystem storage is root-owned')
  call test_observable_contracts(state,metric,full_position,full_velocity,trim(observable_prefix))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0d0,0d0,0d0],&
    continuous,state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0d0,0d0,0d0],&
    continuous,state,ok,message)
  call require(ok,trim(message))
  norm1=real(dot_product(continuous(:,1),matmul(metric,continuous(:,1))),real64)
  energy1=real(dot_product(continuous(:,1),matmul(hamiltonian,continuous(:,1))),real64)
  call require(abs(norm1-norm0)<1d-12,'S norm conservation')
  call require(abs(energy1-energy0)<1d-12,'field-free energy conservation')

  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,split_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0d0,0d0,0d0],&
    split,split_state,ok,message)
  call require(ok,trim(message))
  call write_dg_overlapping_wannier_rt_restart(comm,trim(restart_prefix),split,split_state,ok,message)
  call require(ok,trim(message))
  split=initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,restart_state,ok,message)
  call require(ok,trim(message))
  call read_dg_overlapping_wannier_rt_restart(comm,trim(restart_prefix),split,restart_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0d0,0d0,0d0],&
    split,restart_state,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(split-continuous))<1d-13,'restart-split coefficient identity')
  call require(restart_state%step==state%step.and.abs(restart_state%time-state%time)<1d-15,&
    'restart-split state identity')
  restart_state%hamiltonian_fingerprint=ieor(restart_state%hamiltonian_fingerprint,1_int64)
  call read_dg_overlapping_wannier_rt_restart(comm,trim(restart_prefix),split,restart_state,ok,message)
  call require(.not.ok,'stale RT restart matrix provenance rejection')
  restart_state%hamiltonian_fingerprint=ieor(restart_state%hamiltonian_fingerprint,1_int64)

  zero_ref=initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,zero_ref,zero_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0d0,0d0,0d0],&
    zero_ref,zero_state,ok,message)
  call require(ok,trim(message))
  field_ref=initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,field_ref,field_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0.2d0,0d0,0d0],[0d0,0d0,0d0],&
    field_ref,field_state,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(field_ref-zero_ref))>1d-6,'position field coupling')
  call require(abs(real(dot_product(field_ref(:,1),matmul(metric,field_ref(:,1))),real64)-norm0)<1d-12,&
    'field-coupled S norm conservation')
  velocity_ref=initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,velocity_ref,velocity_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0d0,0d0,0d0],[0.1d0,0d0,0d0],&
    velocity_ref,velocity_state,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(velocity_ref-zero_ref))>1d-6,'velocity field coupling')

  rotation=reshape([cmplx(cos(0.31d0),0d0,8),cmplx(-sin(0.31d0),0d0,8),&
    cmplx(sin(0.31d0),0d0,8),cmplx(cos(0.31d0),0d0,8)],[2,2])
  full_position=(0d0,0d0);full_velocity=(0d0,0d0)
  full_position(1,2,2)=1d0
  full_velocity(1,1,2)=cmplx(0d0,0.1d0,8);full_velocity(1,2,1)=cmplx(0d0,-0.1d0,8)
  rotated_metric=matmul(conjg(transpose(rotation)),matmul(metric,rotation))
  rotated_hamiltonian=matmul(conjg(transpose(rotation)),matmul(hamiltonian,rotation))
  do i=1,3
    rotated_position(i,:,:)=matmul(conjg(transpose(rotation)),matmul(full_position(i,:,:),rotation))
    rotated_velocity(i,:,:)=matmul(conjg(transpose(rotation)),matmul(full_velocity(i,:,:),rotation))
  enddo
  allocate(rotated_srows(size(row_ids),2),rotated_hrows(size(row_ids),2),&
    rotated_xrows(3,size(row_ids),2),rotated_vrows(3,size(row_ids),2))
  do i=1,size(row_ids)
    rotated_srows(i,:)=rotated_metric(row_ids(i),:)
    rotated_hrows(i,:)=rotated_hamiltonian(row_ids(i),:)
    rotated_xrows(:,i,:)=rotated_position(:,row_ids(i),:)
    rotated_vrows(:,i,:)=rotated_velocity(:,row_ids(i),:)
  enddo
  rotated_coeff=matmul(conjg(transpose(rotation)),initial)
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,rotated_srows,rotated_hrows,&
    rotated_xrows,rotated_vrows,17,23,7101_int64,8101_int64,9102_int64,10102_int64,&
    'cell_wrapped_length_velocity',17,23,7101_int64,8101_int64,rotated_coeff,rotation_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0.2d0,0d0,0d0],[0d0,0d0,0d0],&
    rotated_coeff,rotation_state,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(matmul(rotation,rotated_coeff)-field_ref))<1d-12,&
    'retained-basis unitary rotation covariance')

  phase=exp(cmplx(0d0,0.37d0,real64));phase_coeff=phase*initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,phase_coeff,phase_state,ok,message)
  call require(ok,trim(message))
  call advance_dg_overlapping_wannier_rt(comm,0.05d0,[0.2d0,0d0,0d0],[0d0,0d0,0d0],&
    phase_coeff,phase_state,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(phase_coeff-phase*field_ref))<1d-12,'coefficient phase covariance')

  split=initial
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,18,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,split_state,ok,message)
  call require(.not.ok,'stale basis generation rejection')
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8102_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,split_state,ok,message)
  call require(.not.ok,'stale operator fingerprint rejection')
  do i=1,size(row_ids)
    if(row_ids(i)==1)hrows(i,2)=hrows(i,2)+0.2d0
  enddo
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,split_state,ok,message)
  call require(.not.ok,'non-Hermitian Hamiltonian rejection')
  do i=1,size(row_ids)
    if(row_ids(i)==1)hrows(i,2)=hrows(i,2)-0.2d0
    if(row_ids(i)==2)row_ids(i)=1
  enddo
  call initialize_dg_overlapping_wannier_rt(comm,row_ids,srows,hrows,xrows,vrows,17,23,&
    7101_int64,8101_int64,9101_int64,10101_int64,'cell_wrapped_length_velocity',&
    17,23,7101_int64,8101_int64,split,split_state,ok,message)
  call require(.not.ok,'incomplete row ownership rejection')

  if(rank==0)then
    open(newunit=unit,file=trim(restart_prefix),status='old',access='stream',form='unformatted',action='write')
    write(unit,pos=1)'CORRUPT'
    close(unit)
  endif
  call MPI_Barrier(comm,ierr)
  call read_dg_overlapping_wannier_rt_restart(comm,trim(restart_prefix),split,restart_state,ok,message)
  call require(.not.ok,'corrupt RT restart rejection')

  if(rank==0)write(*,'(a,i0,a,z16.16)')'PASS overlapping-Wannier coefficient RT on ',&
    nproc,' ranks fingerprint=',state%fingerprint
  call MPI_Finalize(ierr)
contains
  subroutine test_observable_contracts(observable_state,metric_reference,position_reference,&
      velocity_reference,output_prefix)
    type(s_dg_overlapping_wannier_rt_state),intent(in)::observable_state
    complex(real64),intent(in)::metric_reference(2,2)
    complex(real64),intent(out)::position_reference(3,2,2),velocity_reference(3,2,2)
    character(*),intent(in)::output_prefix
    type(s_dg_overlapping_wannier_rt_state)::sample_state,invalid_state
    complex(real64)::coefficients(2,2),bad_coefficients(2,2),operator_product(2)
    real(real64)::occupations(2),bad_occupations(2),polarization(3),current(3),&
      expected_polarization(3),expected_current(3),field(3),sample_p(3),sample_j(3)
    integer::axis,band,corrupt_unit

    position_reference=(0d0,0d0);velocity_reference=(0d0,0d0)
    position_reference(1,2,2)=1d0
    velocity_reference(1,1,2)=cmplx(0d0,0.1d0,real64)
    velocity_reference(1,2,1)=cmplx(0d0,-0.1d0,real64)
    coefficients(:,1)=[cmplx(0.6d0,0.2d0,real64),cmplx(-0.3d0,0.1d0,real64)]
    coefficients(:,1)=coefficients(:,1)/sqrt(real(dot_product(coefficients(:,1),&
      matmul(metric_reference,coefficients(:,1))),real64))
    coefficients(:,2)=[cmplx(0.1d0,-0.4d0,real64),cmplx(0.5d0,0.2d0,real64)]
    coefficients(:,2)=coefficients(:,2)/sqrt(real(dot_product(coefficients(:,2),&
      matmul(metric_reference,coefficients(:,2))),real64))
    occupations=[2d0,0.5d0]
    expected_polarization=0d0;expected_current=0d0
    do axis=1,3
      do band=1,2
        operator_product=matmul(position_reference(axis,:,:),coefficients(:,band))
        expected_polarization(axis)=expected_polarization(axis)-occupations(band)*&
          real(dot_product(coefficients(:,band),operator_product),real64)/4d0
        operator_product=matmul(velocity_reference(axis,:,:),coefficients(:,band))
        expected_current(axis)=expected_current(axis)-occupations(band)*&
          real(dot_product(coefficients(:,band),operator_product),real64)/4d0
      enddo
    enddo
    call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations,4d0,&
      observable_state,polarization,current,ok,message)
    call require(ok,trim(message))
    call require(maxval(abs(polarization-expected_polarization))<1d-14,&
      'occupation-weighted polarization divided by volume')
    call require(maxval(abs(current-expected_current))<1d-14,&
      'occupation-weighted current divided by volume')

    if(nproc>1)then
      bad_occupations=occupations
      if(rank==0)bad_occupations(2)=0.25d0
      call evaluate_dg_overlapping_wannier_observables(comm,coefficients,bad_occupations,4d0,&
        observable_state,polarization,current,ok,message)
      call require(.not.ok,'collectively inconsistent occupations rejection')
    endif
    call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations(:1),4d0,&
      observable_state,polarization,current,ok,message)
    call require(.not.ok,'wrong occupation count rejection')
    call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations,0d0,&
      observable_state,polarization,current,ok,message)
    call require(.not.ok,'nonpositive volume rejection')
    call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations,&
      ieee_value(4d0,ieee_quiet_nan),observable_state,polarization,current,ok,message)
    call require(.not.ok,'nonfinite volume rejection')
    bad_coefficients=coefficients
    bad_coefficients(1,1)=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,real64)
    call evaluate_dg_overlapping_wannier_observables(comm,bad_coefficients,occupations,4d0,&
      observable_state,polarization,current,ok,message)
    call require(.not.ok,'nonfinite coefficients rejection')
    call evaluate_dg_overlapping_wannier_observables(comm,coefficients,occupations,4d0,&
      invalid_state,polarization,current,ok,message)
    call require(.not.ok,'uninitialized observable state rejection')

    sample_state=observable_state;field=[0.01d0,0d0,0d0]
    sample_p=expected_polarization;sample_j=expected_current
    sample_state%step=0;sample_state%time=0d0
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-one.dat',&
      field,sample_p,sample_j,4d0,sample_state,.false.,ok,message)
    call require(ok,trim(message))
    sample_state%step=1;sample_state%time=0.05d0;sample_p(1)=sample_p(1)+0.001d0
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-one.dat',&
      field,sample_p,sample_j,4d0,sample_state,.false.,ok,message)
    call require(ok,trim(message))
    sample_state%step=2;sample_state%time=0.10d0;sample_p(1)=sample_p(1)+0.001d0
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-one.dat',&
      field,sample_p,sample_j,4d0,sample_state,.false.,ok,message)
    call require(ok,trim(message))

    sample_state%step=0;sample_state%time=0d0;sample_p=expected_polarization
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,4d0,sample_state,.false.,ok,message)
    call require(ok,trim(message))
    sample_state%step=1;sample_state%time=0.05d0;sample_p(1)=sample_p(1)+0.001d0
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,4d0,sample_state,.false.,ok,message)
    call require(ok,trim(message))
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,4d0,sample_state,.true.,ok,message)
    call require(ok,trim(message))
    invalid_state=sample_state
    invalid_state%observable_fingerprint=ieor(invalid_state%observable_fingerprint,1_int64)
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,4d0,invalid_state,.true.,ok,message)
    call require(.not.ok,'stale observable provenance rejection')
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,5d0,sample_state,.true.,ok,message)
    call require(.not.ok,'mismatched observable volume rejection')
    sample_state%step=2;sample_state%time=0.10d0;sample_p(1)=sample_p(1)+0.001d0
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-split.dat',&
      field,sample_p,sample_j,4d0,sample_state,.true.,ok,message)
    call require(ok,trim(message))

    if(rank==0)then
      open(newunit=corrupt_unit,file=trim(output_prefix)//'-corrupt.dat',status='replace')
      write(corrupt_unit,'(a)')'truncated observable evidence'
      close(corrupt_unit)
    endif
    call MPI_Barrier(comm,ierr)
    call write_dg_overlapping_wannier_rt_observable_sample(comm,trim(output_prefix)//'-corrupt.dat',&
      field,sample_p,sample_j,4d0,sample_state,.true.,ok,message)
    call require(.not.ok,'truncated observable file rejection')
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_failure/=0)error stop label
  end subroutine
end program
