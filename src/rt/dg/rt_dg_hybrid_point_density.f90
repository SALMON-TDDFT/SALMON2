#include "config.h"
module rt_dg_hybrid_point_density
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use rt_dg_hybrid_sparse_exchange,only:s_rt_dg_sparse_exchange,exchange_rt_dg_sparse_matrix
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::reconstruct_rt_dg_point_csr_density
contains
  subroutine reconstruct_rt_dg_point_csr_density(comm,plan,point_offsets,support_slots,support_values,&
      coefficients_owned,occupations,density,workspace_peak_bytes,payload_collective_count,ok,message)
    integer,intent(in)::comm,point_offsets(:),support_slots(:)
    type(s_rt_dg_sparse_exchange),intent(in)::plan
    complex(real64),intent(in)::support_values(:),coefficients_owned(:,:)
    real(real64),intent(in)::occupations(:)
    real(real64),intent(out)::density(:)
    integer(int64),intent(out)::workspace_peak_bytes
    integer,intent(out)::payload_collective_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::coefficients_by_halo(:,:),orbital_values(:)
    integer::p,edge,nocc,ierr,local_bad,global_bad,allocation_status
    integer(int64)::exchange_workspace,local_workspace
    ok=.false.;message='';workspace_peak_bytes=0_int64;payload_collective_count=0
    nocc=size(coefficients_owned,2);local_bad=0
    if(.not.plan%valid.or.nocc<1.or.size(occupations)/=nocc.or.size(point_offsets)/=size(density)+1.or.&
      size(support_slots)/=size(support_values))local_bad=1
    if(local_bad==0)then
      if(point_offsets(1)/=1.or.point_offsets(size(point_offsets))/=size(support_slots)+1.or.&
        any(point_offsets<1).or.any(point_offsets>size(support_slots)+1))local_bad=1
      if(size(point_offsets)>1)then
        if(any(point_offsets(2:)<point_offsets(:size(point_offsets)-1)))local_bad=1
      endif
      if(any(support_slots<1).or.any(support_slots>size(plan%value_slots)))local_bad=1
    endif
    if(.not.all(ieee_is_finite(real(support_values))).or..not.all(ieee_is_finite(aimag(support_values))).or.&
      .not.all(ieee_is_finite(real(coefficients_owned))).or.&
      .not.all(ieee_is_finite(aimag(coefficients_owned))).or..not.all(ieee_is_finite(occupations)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid point-CSR Hybrid density contract';return;endif
    allocate(coefficients_by_halo(size(plan%value_slots),nocc),orbital_values(nocc),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='cannot allocate point-CSR Hybrid density workspace';return;endif
    call exchange_rt_dg_sparse_matrix(comm,plan,coefficients_owned,coefficients_by_halo,exchange_workspace,&
      payload_collective_count,ok,message)
    if(.not.ok)then;message='point-CSR Hybrid coefficient halo failed: '//trim(message);return;endif
    do p=1,size(density)
      orbital_values=(0d0,0d0)
      do edge=point_offsets(p),point_offsets(p+1)-1
        orbital_values=orbital_values+support_values(edge)*coefficients_by_halo(support_slots(edge),:)
      enddo
      density(p)=sum(occupations*abs(orbital_values)**2)
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(density)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;ok=.false.;message='nonfinite point-CSR Hybrid density';return;endif
    local_workspace=16_int64*(int(size(coefficients_by_halo),int64)+int(size(orbital_values),int64))
    workspace_peak_bytes=exchange_workspace+local_workspace
    ok=.true.;message=''
#else
    ok=.false.;message='point-CSR Hybrid density requires MPI'
    workspace_peak_bytes=0_int64;payload_collective_count=0
#endif
  end subroutine reconstruct_rt_dg_point_csr_density
end module rt_dg_hybrid_point_density
