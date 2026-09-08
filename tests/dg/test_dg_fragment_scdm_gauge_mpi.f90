#include "config.h"
program test_dg_fragment_scdm_gauge_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_fragment_scdm_gauge,only:build_dg_fragment_scdm_gauge
  implicit none
  integer,parameter::nglobal=8,nband=2
  integer::rank,nproc,ierr,nlocal,p,j,global_id
  integer(int64),allocatable::grid_ids(:),saved_grid_ids(:),selected_ids(:),mixed_selected_ids(:)
  real(real64),allocatable::weights(:)
  complex(real64),allocatable::values(:,:),saved_values(:,:),mixed_values(:,:),a_matrix(:,:),&
    mixed_a_matrix(:,:),replicated_a(:,:),global_values(:,:),global_rotated(:,:),rotated_values(:,:)
  complex(real64)::mixing(nband,nband),phase
  complex(real64)::raw_selected(nband,nband),positive_overlap(nband,nband)
  real(real64)::unitarity_defect,projector_defect,mixed_unitarity,mixed_projector
  integer(int64)::workspace_peak,fingerprint,mixed_workspace,mixed_fingerprint
  logical::ok
  character(512)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(mod(nglobal,nproc)==0,'fixture requires a divisor of eight ranks')
  nlocal=nglobal/nproc
  allocate(grid_ids(nlocal),weights(nlocal),values(nband,nlocal))
  do p=1,nlocal
    grid_ids(p)=int(rank*nlocal+nlocal-p+1,int64)
    global_id=int(grid_ids(p));weights(p)=1d0;values(:,p)=(0d0,0d0)
    if(global_id<=4)then
      values(1,p)=0.5d0
      phase=exp(cmplx(0d0,0.5d0*acos(-1d0)*real(global_id-1,real64),real64))
      values(2,p)=0.5d0*phase
    endif
  enddo
  saved_values=values;saved_grid_ids=grid_ids
  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
    huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,workspace_peak,&
    fingerprint,ok,message)
  call require(ok,'valid SCDM gauge: '//trim(message))
  call require(all(selected_ids==[1_int64,3_int64]),'canonical global-ID pivot order')
  call require(unitarity_defect<1d-11.and.projector_defect<1d-11,'unitary/projector receipts')
  call require(workspace_peak>0_int64.and.workspace_peak<huge(0_int64),'bounded workspace receipt')
  allocate(replicated_a(nband,nband));replicated_a=(0d0,0d0)
  if(rank==0)then
    call require(all(shape(a_matrix)==[nband,nband]),'coordinator owns square SCDM gauge')
    replicated_a=a_matrix
  else
    call require(size(a_matrix)==0,'noncoordinator owns no dense SCDM gauge')
  endif
  call MPI_Bcast(replicated_a,size(replicated_a),MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
  rotated_values=matmul(transpose(replicated_a),values)
  allocate(global_values(nband,nglobal),global_rotated(nband,nglobal))
  global_values=(0d0,0d0);global_rotated=(0d0,0d0)
  do p=1,nlocal
    global_values(:,int(grid_ids(p)))=values(:,p)
    global_rotated(:,int(grid_ids(p)))=rotated_values(:,p)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,global_values,size(global_values),MPI_DOUBLE_COMPLEX,MPI_SUM,&
    MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(MPI_IN_PLACE,global_rotated,size(global_rotated),MPI_DOUBLE_COMPLEX,MPI_SUM,&
    MPI_COMM_WORLD,ierr)
  if(rank==0)then
    do j=1,nband
      raw_selected(:,j)=conjg(global_values(:,int(selected_ids(j))))
    enddo
    positive_overlap=matmul(conjg(transpose(a_matrix)),raw_selected)
    call require(maxval(abs(positive_overlap-conjg(transpose(positive_overlap))))<1d-11,&
      'SCDM polar gauge uses the production A=<retained|trial> conjugation convention')
    call require(all(real([(positive_overlap(j,j),j=1,nband)],real64)>0d0),&
      'SCDM polar overlap is positive in the production orientation')
  endif
  call require(maxval(abs(matmul(conjg(transpose(global_values)),global_values)-&
    matmul(conjg(transpose(global_rotated)),global_rotated)))<1d-11,&
    'SCDM unitary rotation preserves the retained projector')

  mixing=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64),&
    cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)],[nband,nband])/sqrt(2d0)
  mixed_values=matmul(mixing,values)
  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,mixed_values,weights,1d-12,&
    huge(0_int64),mixed_selected_ids,mixed_a_matrix,mixed_unitarity,mixed_projector,&
    mixed_workspace,mixed_fingerprint,ok,message)
  call require(ok.and.all(mixed_selected_ids==selected_ids),&
    'SCDM exact and roundoff-near ties are stable under a retained-subspace gauge rotation')
  call require(mixed_unitarity<1d-11.and.mixed_projector<1d-11,'rotated-input SCDM receipts')

  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
    workspace_peak-1_int64,selected_ids,a_matrix,unitarity_defect,projector_defect,&
    mixed_workspace,mixed_fingerprint,ok,message)
  call require(.not.ok.and.index(message,'byte')>0,'insufficient workspace cap fails closed')
  values=saved_values;values(2,:)=values(1,:)
  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
    huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,mixed_workspace,&
    mixed_fingerprint,ok,message)
  call require(.not.ok.and.(index(message,'rank')>0.or.index(message,'orthonormal')>0),&
    'rank-deficient retained frame fails closed')
  values=saved_values
  if(rank==0)values(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
    huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,mixed_workspace,&
    mixed_fingerprint,ok,message)
  call require(.not.ok.and.index(message,'invalid')>0,'nonfinite retained values fail collectively')
  values=saved_values
  if(nproc>1)then
    call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11+rank,7,grid_ids,values,weights,1d-12,&
      huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,mixed_workspace,&
      mixed_fingerprint,ok,message)
    call require(.not.ok.and.index(message,'disagree')>0,'rank-disagreeing fragment identity fails closed')
    if(rank==nproc-1)grid_ids(1)=1_int64
    call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
      huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,mixed_workspace,&
      mixed_fingerprint,ok,message)
    call require(.not.ok.and.index(message,'grid')>0,'duplicate global grid ID fails collectively')
  endif
  grid_ids=saved_grid_ids
  if(rank==nproc-1)grid_ids(1)=int(nglobal+1,int64)
  call build_dg_fragment_scdm_gauge(MPI_COMM_WORLD,11,7,grid_ids,values,weights,1d-12,&
    huge(0_int64),selected_ids,a_matrix,unitarity_defect,projector_defect,mixed_workspace,&
    mixed_fingerprint,ok,message)
  call require(.not.ok.and.index(message,'grid')>0,'out-of-range global grid ID fails collectively')
  if(rank==0)write(*,'(a,i0,a,i0)')'SCDM_GAUGE ranks=',nproc,' fingerprint=',fingerprint
  if(rank==0)write(*,'(a,i0,a)')'PASS fragment SCDM gauge on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    if(.not.condition)then
      write(*,'(3a)')'FAIL: ',trim(label),''
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_dg_fragment_scdm_gauge_mpi
