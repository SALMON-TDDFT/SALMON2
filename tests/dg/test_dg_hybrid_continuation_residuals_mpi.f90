#include "config.h"
program test_dg_hybrid_continuation_residuals_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:real64
  use dg_hybrid_continuation_residuals,only:s_dg_hybrid_residuals,build_dg_hybrid_occupied_algebra,&
    build_dg_hybrid_interface_observables,evaluate_dg_hybrid_residuals,validate_cluster_occupations,&
    evaluate_dg_hybrid_projector_change,dg_hybrid_electron_count
  implicit none
  integer::comm,rank,nproc,ierr,i
  complex(real64)::s(3,3),c(3,2),rotated(3,2),u(2,2),q(3,3),q2(3,3),gamma(3,3),gamma2(3,3)
  complex(real64)::indefinite_metric(3,3)
  complex(real64)::values(2,3),normals(2,3),tv(2,2),tn(2,2),tc(2,2),tv2(2,2),tn2(2,2),tc2(2,2)
  complex(real64)::h(3,3),hc(3,2),sc_eps(3,2),trace_input(2,2),trace_output(2,2)
  real(real64)::occupations(2),eigenvalues(2),rho_input(3),rho_output(3)
  real(real64)::projector_change,electron_count
  type(s_dg_hybrid_residuals)::residuals
  logical::ok
  character(256)::message

  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  s=(0d0,0d0);s(1,1)=2d0;s(2,2)=3d0;s(3,3)=4d0
  c=(0d0,0d0);c(1,1)=1d0/sqrt(2d0);c(2,2)=1d0/sqrt(3d0)
  u=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64),&
    cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)],[2,2])/sqrt(2d0)
  rotated=matmul(c,u);occupations=[2d0,2d0]
  call build_dg_hybrid_occupied_algebra(c,s,occupations,q,gamma,ok,message);call require(ok,trim(message))
  call build_dg_hybrid_occupied_algebra(rotated,s,occupations,q2,gamma2,ok,message);call require(ok,trim(message))
  call require(maxval(abs(c-rotated))>1d-2,'gauge fixture did not change raw coefficients')
  call require(maxval(abs(q-q2))<1d-13.and.maxval(abs(gamma-gamma2))<1d-13,&
    'occupied algebra changed under an equally occupied unitary rotation')
  call evaluate_dg_hybrid_projector_change(q,q2,s,projector_change,ok,message)
  call require(ok.and.projector_change<1d-13,'S-metric occupied-projector change is gauge dependent')
  electron_count=dg_hybrid_electron_count(gamma,s)
  call require(abs(electron_count-sum(occupations))<1d-13,'electron count was not computed from Gamma and S')
  values=reshape([(cmplx(0.1d0*i,0.03d0*i,real64),i=1,6)],[2,3])
  normals=reshape([(cmplx(-0.04d0*i,0.02d0*i,real64),i=1,6)],[2,3])
  call build_dg_hybrid_interface_observables(values,normals,gamma,tv,tn,tc,ok,message)
  call require(ok,trim(message))
  call build_dg_hybrid_interface_observables(values,normals,gamma2,tv2,tn2,tc2,ok,message)
  call require(ok.and.maxval(abs(tv-tv2))+maxval(abs(tn-tn2))+maxval(abs(tc-tc2))<1d-13,&
    'gauge-invariant interface observables changed under occupied rotation')
  call validate_cluster_occupations([2d0,1d0],[1,1],ok,message)
  call require(.not.ok,'unequal occupations inside one degenerate cluster were accepted')

  eigenvalues=[-0.7d0,-0.2d0];h=(0d0,0d0);h(1,1)=-1.4d0;h(2,2)=-0.6d0;h(3,3)=1d0
  hc=matmul(h,c);sc_eps=matmul(s,c);sc_eps(:,1)=sc_eps(:,1)*eigenvalues(1);sc_eps(:,2)=sc_eps(:,2)*eigenvalues(2)
  rho_input=[1d0,0.5d0,0.25d0];rho_output=rho_input
  trace_input=tv;trace_output=tv
  call evaluate_dg_hybrid_residuals(hc,sc_eps,c,s,rho_output,rho_input,trace_output,trace_input,&
    residuals,ok,message)
  call require(ok.and.residuals%r_h<1d-14.and.residuals%r_rho==0d0.and.residuals%r_t==0d0.and.&
    residuals%r_s<1d-14,'exact continuation residual fixture did not vanish')
  rho_output(1)=rho_output(1)+0.1d0;trace_output(1,1)=trace_output(1,1)+(0.2d0,0d0)
  call evaluate_dg_hybrid_residuals(hc,sc_eps,c,s,rho_output,rho_input,trace_output,trace_input,&
    residuals,ok,message)
  call require(ok.and.residuals%r_rho>0d0.and.residuals%r_t>0d0.and.residuals%r_h<1d-14,&
    'density and interface residual channels are not independent')
  indefinite_metric=s;indefinite_metric(3,3)=-4d0
  call build_dg_hybrid_occupied_algebra(c,indefinite_metric,occupations,q,gamma,ok,message)
  call require(.not.ok,'indefinite DG metric was accepted outside the occupied subspace')
  if(rank==0)write(*,'(a,i0,a)')'PASS hybrid continuation residuals on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_continuation_residuals_mpi
