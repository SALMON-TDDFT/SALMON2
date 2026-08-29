subroutine build_dg_hybrid_retained_basis_representation(comm_arg,global_count_arg,core_ids_arg,&
    core_weights_arg,pencil_maps_arg,fragment_basis_arg,row_ids_arg,s_rows_arg,tolerance_arg,&
    representation_arg,callback_ok,callback_message)
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_overlapping_wannier_construction,only:assemble_dg_distributed_basis_symmetry_overlap
  implicit none
  integer,intent(in)::comm_arg,global_count_arg
  integer(int64),intent(in)::core_ids_arg(:),pencil_maps_arg(:,:),row_ids_arg(:)
  real(real64),intent(in)::core_weights_arg(:),tolerance_arg
  complex(real64),intent(in)::s_rows_arg(:,:)
  type(s_dg_hybrid_fragment_basis),intent(in)::fragment_basis_arg
  complex(real64),allocatable,intent(out)::representation_arg(:,:,:)
  logical,intent(out)::callback_ok
  character(*),intent(out)::callback_message
  complex(real64),allocatable::local_basis(:,:),overlap(:,:,:),local_metric(:,:),metric(:,:),metric_work(:,:)
  integer,allocatable::ownership(:),pivot(:)
  integer::nbasis,ncore,noperation,i,j,operation,point,position,info,ierr,local_bad,global_bad
  real(real64)::local_defect,global_defect,scale

  callback_ok=.false.;callback_message=''
  nbasis=size(s_rows_arg,2);ncore=size(core_ids_arg);noperation=size(pencil_maps_arg,2)
  local_bad=merge(0,1,global_count_arg>0.and.nbasis>0.and.ncore>0.and.noperation>0.and.&
    size(core_weights_arg)==ncore.and.size(pencil_maps_arg,1)==ncore.and.&
    size(s_rows_arg,1)==size(row_ids_arg).and.tolerance_arg>0d0.and.&
    all(core_ids_arg>0_int64).and.all(core_ids_arg<=int(global_count_arg,int64)).and.&
    all(core_weights_arg>=0d0).and.all(ieee_is_finite(core_weights_arg)).and.&
    all(ieee_is_finite(real(s_rows_arg))).and.all(ieee_is_finite(aimag(s_rows_arg))))
  call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr)
  if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
    callback_message='invalid retained-basis symmetry representation contract';return
  endif
  allocate(local_basis(nbasis,ncore),local_metric(nbasis,nbasis),metric(nbasis,nbasis),&
    metric_work(nbasis,nbasis),ownership(nbasis),pivot(nbasis))
  local_basis=(0d0,0d0);local_metric=(0d0,0d0);ownership=0
  do j=1,size(fragment_basis_arg%global_ids)
    i=int(fragment_basis_arg%global_ids(j))
    if(i<1.or.i>nbasis)then;local_bad=1;cycle;endif
    do point=1,ncore
      position=findloc(fragment_basis_arg%buffer_point_ids,core_ids_arg(point),dim=1)
      if(position>0)local_basis(i,point)=fragment_basis_arg%buffer_values(position,j)
    enddo
  enddo
  do i=1,size(row_ids_arg)
    j=int(row_ids_arg(i))
    if(j<1.or.j>nbasis)then;local_bad=1;cycle;endif
    ownership(j)=ownership(j)+1;local_metric(j,:)=s_rows_arg(i,:)
  enddo
  call MPI_Allreduce(MPI_IN_PLACE,ownership,nbasis,MPI_INTEGER,MPI_SUM,comm_arg,ierr)
  if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_metric,metric,nbasis*nbasis,MPI_DOUBLE_COMPLEX,MPI_SUM,comm_arg,ierr)
  if(ierr/=MPI_SUCCESS.or.any(ownership/=1).or.local_bad/=0)then
    callback_message='retained-basis metric rows are incomplete';return
  endif
  call assemble_dg_distributed_basis_symmetry_overlap(comm_arg,local_basis,core_weights_arg,&
    pencil_maps_arg,overlap,callback_ok,callback_message)
  if(.not.callback_ok)return
  allocate(representation_arg(nbasis,nbasis,noperation))
  do operation=1,noperation
    metric_work=metric;representation_arg(:,:,operation)=overlap(:,:,operation)
    call zgesv(nbasis,nbasis,metric_work,nbasis,pivot,representation_arg(:,:,operation),nbasis,info)
    local_bad=merge(0,1,info==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm_arg,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      callback_message='retained-basis metric solve failed';return
    endif
  enddo
  local_defect=0d0;scale=max(1d0,sqrt(sum(abs(metric)**2)))
  do operation=1,noperation
    local_defect=max(local_defect,sqrt(sum(abs(matmul(metric,representation_arg(:,:,operation))-&
      overlap(:,:,operation))**2))/scale)
    local_defect=max(local_defect,sqrt(sum(abs(matmul(conjg(transpose(representation_arg(:,:,operation))),&
      matmul(metric,representation_arg(:,:,operation)))-metric)**2))/scale)
  enddo
  call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm_arg,ierr)
  callback_ok=ierr==MPI_SUCCESS.and.ieee_is_finite(global_defect).and.global_defect<=tolerance_arg
  if(callback_ok)then;callback_message='';else;callback_message='retained WF+PW basis is not symmetry closed';endif
end subroutine build_dg_hybrid_retained_basis_representation
