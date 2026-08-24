#include "config.h"
program test_dg_spatial_mlwf_seed_mpi
  use mpi
  use iso_fortran_env, only: real64, int64
  use dg_spatial_mlwf_seed, only: spatial_seed_options, spatial_seed_diagnostics,spatial_seed_sparse_matrix, &
    make_spatial_seed_sparse_matrix,construct_spatial_mlwf_seed,apply_s_projector,certified_sparse_metric_solve
  implicit none
  integer,parameter::n=4,no=2,ne=2,ng=2
  complex(real64)::h(n,n),s(n,n),yocc(n,no),yempty(n,ne),bg(n,n,ng)
  complex(real64)::docc(no,no,ng),dempty(ne,ne,ng),a(n,n),tg(n,n,ng)
  complex(real64),allocatable::cocc(:,:),cempty(:,:),pocc(:,:),pempty(:,:)
  complex(real64),allocatable::metric_solution(:)
  complex(real64)::metric_rhs(n)
  type(spatial_seed_options)::options
  type(spatial_seed_diagnostics)::diagnostics
  type(spatial_seed_sparse_matrix)::sparse_h,sparse_s
  logical::ok
  character(256)::message
  integer::i,ierr,rank,nproc
  integer(int64)::fingerprint,minimum_fingerprint,maximum_fingerprint
  complex(real64)::saved_h(n,n),saved_a(n,n),saved_yocc(n,no),saved_tg(n,n,ng)

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  h=(0d0,0d0);s=(0d0,0d0);a=(0d0,0d0);bg=(0d0,0d0);tg=(0d0,0d0)
  do i=1,n
    h(i,i)=cmplx(real(i-1,real64),0d0,real64)
    s(i,i)=cmplx(1d0+0.1d0*real(i-1,real64),0d0,real64)
    a(i,i)=sqrt(s(i,i))
    bg(i,i,1)=1d0;tg(i,i,1)=1d0
  end do
  bg(2,1,2)=sqrt(real(s(1,1)/s(2,2),real64));bg(1,2,2)=sqrt(real(s(2,2)/s(1,1),real64))
  bg(4,3,2)=sqrt(real(s(3,3)/s(4,4),real64));bg(3,4,2)=sqrt(real(s(4,4)/s(3,3),real64))
  tg(2,1,2)=1d0;tg(1,2,2)=1d0;tg(4,3,2)=1d0;tg(3,4,2)=1d0
  docc=(0d0,0d0);dempty=(0d0,0d0)
  docc(1,1,1)=1d0;docc(2,2,1)=1d0;dempty(1,1,1)=1d0;dempty(2,2,1)=1d0
  docc(2,1,2)=1d0;docc(1,2,2)=1d0;dempty(2,1,2)=1d0;dempty(1,2,2)=1d0
  yocc=(0d0,0d0);yempty=(0d0,0d0)
  yocc(1,1)=1d0;yocc(2,2)=0.7d0;yocc(1,2)=0.3d0
  yempty(3,1)=1d0;yempty(4,2)=0.8d0;yempty(3,2)=0.2d0
  options%lower_bound=-0.1d0;options%upper_bound=3.1d0
  options%occupied_upper=1.1d0;options%empty_lower=1.4d0
  options%empty_upper=3.1d0;options%filter_degree=8
  options%tolerance=1d-9;options%rank_tolerance=1d-10
  options%condition_limit=10d0;options%fill_limit=n*n
  call make_spatial_seed_sparse_matrix(h,1d-14,sparse_h,ok,message);call require(ok,'sparse H fixture')
  call make_spatial_seed_sparse_matrix(s,1d-14,sparse_s,ok,message);call require(ok,'sparse S fixture')
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,sparse_h,sparse_s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(ok,'valid generalized sparse seed: '//trim(message))
  call require(size(cocc,2)==no.and.size(cempty,2)==ne,'complete occupied and empty ranks')
  call require(diagnostics%metric_defect<options%tolerance,'metric purification residual')
  call require(diagnostics%back_transform_defect<options%tolerance,'C=XQ back-transformation')
  call require(diagnostics%dense_seed_workspace_bytes==0_int64,'sparse route forms no dense H/S/X/Hbar')
  call require(diagnostics%sparse_nnz_peak<=options%fill_limit,'sparse intermediate fill is bounded')
  call require(diagnostics%sparse_workspace_peak_bytes<=512_int64*n,&
    'sparse seed workspace grows linearly for the replicated diagonal fixture')
  metric_rhs=[(1d0,0d0),(2d0,0d0),(-1d0,0.5d0),(0.3d0,-0.2d0)]
  call certified_sparse_metric_solve(sparse_s,metric_rhs,options%tolerance,40,metric_solution,i,&
    diagnostics%metric_solve_residual,ok,message)
  call require(ok.and.i<=4,'certified sparse S metric solve is bounded')
  call require(maxval(abs(matmul(s,metric_solution)-metric_rhs))<options%tolerance,&
    'certified sparse S metric solve residual')
  call require(diagnostics%occupied_projector_defect<options%tolerance,'occupied projector idempotency')
  call require(diagnostics%empty_projector_defect<options%tolerance,'empty projector idempotency')
  call require(diagnostics%mutual_projector_defect<options%tolerance,'occupied-empty S orthogonality')
  call require(diagnostics%map_metric_defect<options%tolerance,'stitching metric pullback')
  call require(diagnostics%map_linearity_defect<options%tolerance,'stitching linearity')
  call require(diagnostics%map_equivariance_defect<options%tolerance,'stitching equivariance')
  call require(diagnostics%intertwining_defect<options%tolerance,'Reynolds-projected intertwiner')
  call apply_s_projector(cocc,s,pocc,ok,message);call require(ok,'occupied S projector')
  call apply_s_projector(cempty,s,pempty,ok,message);call require(ok,'empty S projector')
  call require(maxval(abs(matmul(pocc,pocc)-pocc))<1d-9,'P_C idempotency')
  call require(maxval(abs(matmul(conjg(transpose(pocc)),s)-matmul(s,pocc)))<1d-9,'P_C S self-adjointness')
  call require(maxval(abs(matmul(pocc,pempty)))<1d-9,'P_occ P_empty')
  call MPI_Allreduce(fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierr)
  call require(minimum_fingerprint==maximum_fingerprint,'MPI-layout invariant fingerprint')

  saved_h=h;saved_a=a;saved_yocc=yocc;saved_tg=tg
  h=(0d0,0d0);h(1,1)=1d0;h(2,2)=1.1d0;h(3,3)=2.4d0;h(4,4)=2.6d0
  h(1,2)=1.1d0;h(2,1)=1.1d0
  options%lower_bound=0d0;options%upper_bound=3.1d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'spectral')>0,&
    'verified spectral bound detects an off-diagonal lower-bound violation')
  h=saved_h;options%lower_bound=-0.1d0

  options%lower_bound=0.2d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'spectral')>0,'invalid spectral bound fails closed')
  options%lower_bound=-0.1d0;options%empty_lower=1.1d0;options%occupied_upper=1.1d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'gap')>0,'closed gap fails closed')
  options%empty_lower=1.4d0;s(4,4)=-1d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'positive definite')>0,'indefinite metric fails closed')
  s(4,4)=1.3d0;a(1,1)=1.1d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'map metric')>0,'map metric failure fails closed')
  a=saved_a;tg=saved_tg;tg(1,2,2)=-1d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'equivariance')>0,'map equivariance failure fails closed')
  tg=saved_tg;yocc=saved_yocc;yocc(:,2)=yocc(:,1)
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'rank')>0,'complete-irrep rank loss fails closed')
  yocc=saved_yocc;options%condition_limit=1.05d0
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'condition')>0,'ill-conditioned metric fails closed')
  options%condition_limit=10d0;options%fill_limit=1
  call construct_spatial_mlwf_seed(MPI_COMM_WORLD,h,s,yocc,yempty,bg,docc,dempty,a,tg,&
    options,cocc,cempty,diagnostics,fingerprint,ok,message)
  call require(.not.ok.and.index(message,'fill')>0,'fill growth fails closed')
  if(rank==0)write(*,'(a,i0,a)')'PASS spatial MLWF sparse seed on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    if(.not.condition)then
      write(*,'(3a)')'FAIL: ',trim(label),'';call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    end if
  end subroutine
end program
