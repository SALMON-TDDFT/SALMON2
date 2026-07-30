#include "config.h"
program test_dg_overlapping_wannier_checkpoint_mpi
  use mpi
  use iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_positive_inf
  use dg_overlapping_wannier_checkpoint
  implicit none
  integer::comm,rank,nproc,ierr,nlocal,i,unit,ntail,j,nrow
  integer(int64)::original_publication_id
  character(512)::prefix,prefix2,prefix_bad,shard
  type(s_dg_overlapping_wannier_checkpoint)::a,b
  logical::ok,reusable,exists
  character(256)::message
  real(8),parameter::current_gates(6)=[1d-9,1d-9,1d-9,1d-9,10d0,1d-9]
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  call get_environment_variable('OW_CHECKPOINT_PREFIX',prefix)
  write(prefix(len_trim(prefix)+1:),'(a,i0)')'-',nproc
  prefix2=trim(prefix)//'-roundtrip'
  prefix_bad=trim(prefix)//'-bad'
  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  a%basis_generation=7;a%geometry_generation=3
  a%basis_fingerprint=700_int64;a%operator_fingerprint=900_int64
  a%density_residual=1d-12;a%coefficient_residual=2d-12;a%charge_error=3d-13;a%accepted=.true.
  a%unmixed_density_residual=4d-12;a%orthogonality_defect=5d-12;a%metric_condition=2d0
  a%density_tolerance=1d-9;a%coefficient_tolerance=1d-9;a%orthogonality_tolerance=1d-9
  a%charge_tolerance=1d-9;a%condition_limit=10d0
  a%symmetry_closure_residual=1d-12;a%symmetry_tolerance=1d-9
  allocate(a%center_owner(2),a%core_physical_ids(nlocal),a%coefficients(2,1),&
    a%occupations(1),a%density(nlocal))
  a%center_owner=[0,mod(1,nproc)]
  nrow=count(a%center_owner==rank)
  allocate(a%overlap_row_ids(nrow),a%overlap(nrow,2));a%overlap=(0d0,0d0)
  j=0
  do i=1,2
    if(a%center_owner(i)/=rank)cycle
    j=j+1;a%overlap_row_ids(j)=i;a%overlap(j,i)=1d0
  enddo
  a%coefficients=reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)],[2,1]);a%occupations=1d0
  ntail=count(a%center_owner==rank)
  allocate(a%tail_center(ntail),a%tail_generation(ntail),a%tail_offsets(ntail+1),&
    a%tail_physical_ids(2*ntail))
  j=0;a%tail_offsets(1)=1
  do i=1,2
    if(a%center_owner(i)/=rank)cycle
    j=j+1;a%tail_center(j)=i;a%tail_generation(j)=7;a%tail_offsets(j+1)=2*j+1
    a%tail_physical_ids(2*j-1:2*j)=[1_int64,2_int64]
  enddo
  nlocal=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;a%core_physical_ids(nlocal)=i;a%density(nlocal)=0.25d0*i
  enddo
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix),a,ok,message)
  call require(ok,trim(message))
  original_publication_id=a%publication_id
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(ok.and.reusable,trim(message))
  call require(all(b%center_owner==a%center_owner).and.all(b%tail_physical_ids==a%tail_physical_ids).and.&
    all(b%tail_center==a%tail_center).and.all(b%core_physical_ids==a%core_physical_ids).and.&
    all(b%overlap_row_ids==a%overlap_row_ids),&
    'checkpoint ownership round trip')
  call require(all(b%overlap==a%overlap).and.all(b%coefficients==a%coefficients).and.all(b%density==a%density),&
    'checkpoint payload round trip')
  do i=1,size(a%overlap_row_ids)
    if(a%overlap_row_ids(i)==2_int64)a%overlap_row_ids(i)=1_int64
  enddo
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix_bad),a,ok,message)
  call require(.not.ok,'duplicate/missing overlap row rejected')
  do i=1,size(a%overlap_row_ids)
    if(a%center_owner(2)==rank.and.a%overlap_row_ids(i)==1_int64.and.i>1)&
      a%overlap_row_ids(i)=2_int64
    if(a%center_owner(2)==rank.and.nproc>1.and.a%overlap_row_ids(i)==1_int64)&
      a%overlap_row_ids(i)=2_int64
  enddo
  a%accepted=.false.;a%charge_error=2d-8
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix),a,ok,message)
  call require(.not.ok,'unaccepted checkpoint publication rejected')
  call require(index(message,'code=65')>0,'checkpoint rejection identifies failed invariant')
  a%accepted=.true.;a%charge_error=3d-13
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(ok.and.reusable,'failed publication preserves prior atomic checkpoint')
  a%density=a%density+0.125d0
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix),a,ok,message,failure_injection_rank=0)
  call require(.not.ok,'injected post-shard failure rejected before manifest commit')
  call require(a%publication_id==original_publication_id,'failed publication restores publication id')
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(ok.and.reusable.and.all(b%density/=a%density),&
    'same-provenance failed publication preserves prior checkpoint')
  a%density=a%density-0.125d0
  if(rank==max(0,nproc-1))a%operator_fingerprint=901_int64
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix_bad),a,ok,message)
  if(nproc>1)call require(.not.ok,'rank-inconsistent manifest payload rejected')
  a%operator_fingerprint=900_int64
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix2),b,ok,message)
  call require(ok,'write/read/write equivalence')
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),8,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(ok.and..not.reusable,'stale generation rejection')
  if(rank==0)then
    write(shard,'(a,".g7-",z16.16,".t",z16.16,".rank",i8.8)')trim(prefix),900_int64,&
      original_publication_id,max(0,nproc-1)
    open(newunit=unit,file=trim(shard),status='replace',access='stream',form='unformatted');write(unit)'partial';close(unit)
  endif
  call MPI_Barrier(comm,ierr)
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(.not.ok.and..not.reusable,'partial shard rejection')
  if(rank==0)then
    open(newunit=unit,file=trim(prefix)//'.manifest',status='replace',access='stream',form='unformatted')
    write(unit)'NORMAL_DC_CHECKPOINT';close(unit)
  endif
  call MPI_Barrier(comm,ierr)
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,current_gates,&
    b,reusable,ok,message)
  call require(.not.ok.and..not.reusable,'normal/direct checkpoint rejection')
  if(rank==0)then
    ! Recreate an accepted route checkpoint, then require current (tighter) policy gates.
    call remove_test_file(trim(prefix)//'.manifest')
  endif
  call MPI_Barrier(comm,ierr)
  call write_dg_overlapping_wannier_checkpoint(comm,trim(prefix),a,ok,message)
  call require(ok,'checkpoint rewrite for current-policy gate')
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,&
    [1d-13,1d-9,1d-9,1d-9,10d0,1d-9],b,reusable,ok,message)
  call require(ok.and..not.reusable,'current tighter acceptance policy rejects checkpoint reuse')
  call read_dg_overlapping_wannier_checkpoint(comm,trim(prefix),7,3,700_int64,900_int64,&
    merge([ieee_value(0d0,ieee_positive_inf),1d-9,1d-9,1d-9,10d0,1d-9],current_gates,rank==0),&
    b,reusable,ok,message)
  call require(.not.ok.and..not.reusable,'rank-skewed non-finite current policy rejected collectively')
  if(rank==0)then
    inquire(file=trim(prefix)//'.manifest.temporary',exist=exists)
    call require(.not.exists,'no shared temporary manifest artifact')
  endif
  if(rank==0)write(*,'(a,i0,a)')'PASS overlapping-Wannier checkpoint on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    if(.not.condition)then;write(0,'(a)')trim(text);error stop 1;endif
  end subroutine
  subroutine remove_test_file(filename)
    character(*),intent(in)::filename
    integer::delete_unit,delete_ios
    open(newunit=delete_unit,file=filename,status='old',iostat=delete_ios)
    if(delete_ios==0)close(delete_unit,status='delete')
  end subroutine
end program
