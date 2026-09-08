#include "config.h"
program test_dg_hybrid_continuation_residuals_mpi
  use mpi, only: MPI_Allreduce, MPI_Comm_rank, MPI_Comm_size, MPI_COMM_WORLD, MPI_Finalize, MPI_Init, MPI_INTEGER, MPI_MAX, &
    MPI_SUCCESS
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_continuation_residuals,only:s_dg_hybrid_residuals,s_dg_hybrid_metric_receipt,&
    validate_dg_hybrid_occupied_rows,build_dg_hybrid_interface_observables,&
    evaluate_dg_hybrid_residuals,evaluate_dg_hybrid_projector_change,dg_hybrid_electron_count,&
    validate_cluster_occupations
  implicit none
  integer::icomm,id_rank,nproc,ierr,i,nlocal
  integer(int64),allocatable::row_ids(:)
  complex(real64)::s(3,3),c(3,2),rotated(3,2),u(2,2),h(3,3),hc_all(3,2),sc_all(3,2),s_c(3,2),rotated_s(3,2)
  complex(real64),allocatable::c_local(:,:),rotated_local(:,:),sc_local(:,:),rotated_sc_local(:,:),hc_local(:,:)
  complex(real64)::values(2,2),normals(2,2),tv(2,2),tn(2,2),tc(2,2),tv2(2,2),tn2(2,2),tc2(2,2)
  real(real64),allocatable::rho_input(:),rho_output(:)
  real(real64)::occupations(2),eigenvalues(2),projector_change,electron_count
  type(s_dg_hybrid_metric_receipt)::metric_receipt
  type(s_dg_hybrid_residuals)::residuals
  logical::ok
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  row_ids=pack([1_int64,2_int64,3_int64],[(mod(i-1,nproc)==id_rank,i=1,3)]);nlocal=size(row_ids)
  s=(0d0,0d0);s(1,1)=2d0;s(2,2)=3d0;s(3,3)=4d0
  c=(0d0,0d0);c(1,1)=1d0/sqrt(2d0);c(2,2)=1d0/sqrt(3d0)
  u=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64),cmplx(0d0,1d0,real64),&
    cmplx(1d0,0d0,real64)],[2,2])/sqrt(2d0)
  rotated=matmul(c,u);occupations=[2d0,2d0]
  allocate(c_local(nlocal,2),rotated_local(nlocal,2),sc_local(nlocal,2),rotated_sc_local(nlocal,2))
  c_local=c(int(row_ids),:);rotated_local=rotated(int(row_ids),:)
  s_c=matmul(s,c);sc_all=s_c;sc_local=sc_all(int(row_ids),:)
  rotated_s=matmul(s,rotated);rotated_sc_local=rotated_s(int(row_ids),:)
  metric_receipt%valid=.true.;metric_receipt%global_count=3;metric_receipt%numerical_rank=3
  metric_receipt%fingerprint=771_int64
  call validate_dg_hybrid_occupied_rows(icomm,row_ids,c_local,sc_local,occupations,metric_receipt,ok,message)
  call require(ok,trim(message))
  call evaluate_dg_hybrid_projector_change(icomm,c_local,rotated_sc_local,metric_receipt,&
    projector_change,ok,message)
  call require(ok.and.projector_change<1d-13,'principal-angle occupied-subspace change is gauge dependent')
  electron_count=dg_hybrid_electron_count(occupations)
  call require(abs(electron_count-sum(occupations))<1d-13,'electron count was not computed from occupations')
  call validate_cluster_occupations([2d0,1d0],[1,1],ok,message)
  call require(.not.ok,'unequal occupations inside one degenerate cluster were accepted')
  values=reshape([(cmplx(0.1d0*i,0.03d0*i,real64),i=1,4)],[2,2])
  normals=reshape([(cmplx(-0.04d0*i,0.02d0*i,real64),i=1,4)],[2,2])
  call build_dg_hybrid_interface_observables(values,normals,occupations,tv,tn,tc,ok,message)
  call require(ok,trim(message))
  call build_dg_hybrid_interface_observables(matmul(values,u),matmul(normals,u),occupations,tv2,tn2,tc2,ok,message)
  call require(ok.and.maxval(abs(tv-tv2))+maxval(abs(tn-tn2))+maxval(abs(tc-tc2))<1d-13,&
    'gauge-invariant interface observables changed under occupied rotation')

  eigenvalues=[-0.7d0,-0.2d0];h=(0d0,0d0);h(1,1)=-1.4d0;h(2,2)=-0.6d0;h(3,3)=1d0
  hc_all=matmul(h,c);sc_all(:,1)=sc_all(:,1)*eigenvalues(1);sc_all(:,2)=sc_all(:,2)*eigenvalues(2)
  allocate(hc_local(nlocal,2),rho_input(nlocal),rho_output(nlocal))
  hc_local=hc_all(int(row_ids),:);sc_local=sc_all(int(row_ids),:)
  rho_input=[(1d0/(real(row_ids(i),real64)),i=1,nlocal)];rho_output=rho_input
  call evaluate_dg_hybrid_residuals(icomm,hc_local,sc_local,c_local,s_c(int(row_ids),:),&
    rho_output,rho_input,tv,tv,residuals,ok,message)
  call require(ok.and.residuals%r_h<1d-14.and.residuals%r_rho==0d0.and.residuals%r_t==0d0.and.&
    residuals%r_s<1d-14,'distributed exact continuation residual fixture did not vanish')
  metric_receipt%fingerprint=0_int64
  call validate_dg_hybrid_occupied_rows(icomm,row_ids,c_local,s_c(int(row_ids),:),occupations,&
    metric_receipt,ok,message)
  call require(.not.ok,'invalid frozen-metric receipt was accepted')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS hybrid continuation residuals on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;if(id_rank==0)write(0,'(a)')trim(label);error stop 1;endif
  end subroutine require
end program test_dg_hybrid_continuation_residuals_mpi
