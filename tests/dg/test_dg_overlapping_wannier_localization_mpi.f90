#include "config.h"
program test_dg_overlapping_wannier_localization_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_overlapping_wannier_localization,only:evaluate_dg_periodic_localization,&
    optimize_dg_wannier_pair,build_dg_overlapping_pair_graph
  use dg_overlapping_wannier_localization,only:localize_dg_overlapping_wannier_basis
  implicit none
  complex(8)::values(2,4),phases(3,4),shifted_phases(3,4),moment(3,2),shifted_moment(3,2)
  complex(8)::gradients(3,2,4),rotation(2,2),identity(2,2)
  complex(8)::graph_values(4,4)
  complex(8)::sweep_values(4,4),sweep_gradients(3,4,4),sweep_representation(4,4,2)
  complex(8)::sweep_identity(4,4)
  complex(8),allocatable::sweep_transform(:,:)
  integer,allocatable::pair_first(:),pair_second(:)
  integer::sweep_product(2,2),sweep_iterations
  real(8),allocatable::pair_support(:)
  real(8)::initial_spread,final_spread,maximum_pair_gradient
  logical::converged
  real(8)::weights(4),norm(2),shifted_norm(2),spread,shifted_spread,nan_value,&
    before,after,pair_gradient,density_before(4),gradient_norm_before
  logical::ok,accepted
  character(256)::message
  integer::ierr,rank,nproc,point

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  weights=1d0;values=(0d0,0d0);phases=(1d0,0d0)
  phases(1,:)=[(1d0,0d0),(0d0,1d0),(-1d0,0d0),(0d0,-1d0)]
  values(1,1)=1d0
  values(2,1)=1d0/sqrt(2d0);values(2,3)=1d0/sqrt(2d0)
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(ok,'valid periodic localization payload')
  call require(maxval(abs(norm-1d0))<1d-14,'Wannier norms')
  call require(abs(moment(1,1)-1d0)<1d-14,'delta-localized phase moment')
  call require(abs(moment(1,2))<1d-14,'opposite-site delocalized phase moment')
  call require(abs(spread-1d0)<1d-14,'bounded periodic spread distinguishes localization')

  do point=1,4
    shifted_phases(1,point)=phases(1,point)*exp(cmplx(0d0,0.37d0,8))
    shifted_phases(2,point)=phases(2,point)*exp(cmplx(0d0,-0.21d0,8))
    shifted_phases(3,point)=phases(3,point)*exp(cmplx(0d0,0.13d0,8))
  end do
  call evaluate_dg_periodic_localization(values,weights,shifted_phases,shifted_norm,&
    shifted_moment,shifted_spread,ok,message)
  call require(ok.and.abs(shifted_spread-spread)<1d-14,'periodic spread is origin invariant')

  values=(0d0,0d0);gradients=(0d0,0d0)
  values(1,1)=cos(0.3d0);values(1,3)=sin(0.3d0)
  values(2,1)=-sin(0.3d0);values(2,3)=cos(0.3d0)
  gradients(1,1,1)=cos(0.3d0);gradients(1,1,3)=sin(0.3d0)
  gradients(1,2,1)=-sin(0.3d0);gradients(1,2,3)=cos(0.3d0)
  density_before=sum(abs(values)**2,dim=1)
  gradient_norm_before=sum(abs(gradients)**2)
  call optimize_dg_wannier_pair(values,gradients,weights,phases,1,2,1d-13,&
    rotation,before,after,pair_gradient,accepted,ok,message)
  call require(ok.and.accepted,'mixed pair accepts localization rotation')
  call require(after<before-1d-10,'pair rotation strictly lowers spread')
  call require(pair_gradient>1d-10,'mixed pair has nonzero localization gradient')
  identity=matmul(conjg(transpose(rotation)),rotation)
  identity(1,1)=identity(1,1)-1d0;identity(2,2)=identity(2,2)-1d0
  call require(maxval(abs(identity))<1d-12,'accepted pair rotation is unitary')
  call require(maxval(abs(sum(abs(values)**2,dim=1)-density_before))<1d-12,&
    'pair rotation preserves pointwise density')
  call require(abs(sum(abs(gradients)**2)-gradient_norm_before)<1d-12,&
    'pair rotation transforms all gradients unitarily')
  call optimize_dg_wannier_pair(values,gradients,weights,phases,1,2,1d-13,&
    rotation,before,after,pair_gradient,accepted,ok,message)
  call require(ok.and..not.accepted.and.abs(after-before)<1d-13,&
    'stationary localized pair remains unchanged')

  graph_values=(0d0,0d0)
  graph_values(1,1)=1d0;graph_values(2,1)=0.8d0;graph_values(2,2)=0.6d0
  graph_values(3,3)=1d0;graph_values(4,3)=0.8d0;graph_values(4,4)=0.6d0
  call build_dg_overlapping_pair_graph(MPI_COMM_WORLD,graph_values,weights,0.2d0,&
    pair_first,pair_second,pair_support,ok,message)
  call require(ok.and.size(pair_first)==2,'sparse overlapping-pair graph size')
  call require(all(pair_first==[1,3]).and.all(pair_second==[2,4]),&
    'pair graph retains only shared buffered support')
  call require(all(pair_support>0.8d0).and.all(pair_support<=1d0),&
    'pair support is rank-independent and normalized')

  sweep_values=(0d0,0d0);sweep_gradients=(0d0,0d0)
  sweep_values(1,1)=cos(0.3d0);sweep_values(1,3)=sin(0.3d0)
  sweep_values(3,1)=-sin(0.3d0);sweep_values(3,3)=cos(0.3d0)
  sweep_values(2,2)=cos(0.3d0);sweep_values(2,4)=sin(0.3d0)
  sweep_values(4,2)=-sin(0.3d0);sweep_values(4,4)=cos(0.3d0)
  sweep_gradients(1,:,:)=sweep_values
  sweep_representation=(0d0,0d0)
  do point=1,4;sweep_representation(point,point,1)=1d0;end do
  sweep_representation(2,1,2)=1d0;sweep_representation(1,2,2)=1d0
  sweep_representation(4,3,2)=1d0;sweep_representation(3,4,2)=1d0
  sweep_product=reshape([1,2,2,1],[2,2])
  call localize_dg_overlapping_wannier_basis(MPI_COMM_WORLD,sweep_values,sweep_gradients,&
    weights,phases,sweep_representation,sweep_product,0.1d0,1d-16,1d-7,1d-12,32,&
    initial_spread,final_spread,maximum_pair_gradient,sweep_iterations,converged,&
    sweep_transform,ok,message)
  if(rank==0.and..not.ok)write(*,'(2a,3(a,es12.4),a,i0)')'SWEEP-DIAGNOSTIC ',trim(message),&
    ' initial=',initial_spread,' final=',final_spread,' gradient=',maximum_pair_gradient,&
    ' iterations=',sweep_iterations
  call require(ok.and.converged,'symmetry-constrained localization sweep converges')
  call require(final_spread<initial_spread-1d-8,'localization sweep lowers total spread')
  sweep_identity=matmul(conjg(transpose(sweep_transform)),sweep_transform)
  do point=1,4;sweep_identity(point,point)=sweep_identity(point,point)-1d0;end do
  call require(maxval(abs(sweep_identity))<1d-11,'sweep transform is unitary')
  call require(maxval(abs(matmul(sweep_transform,sweep_representation(:,:,2))-&
    matmul(sweep_representation(:,:,2),sweep_transform)))<1d-11,&
    'accepted sweep transform commutes with exact symmetry')
  call require(maximum_pair_gradient<1d-7,'published localization gradient is converged')

  sweep_values=(0d0,0d0);sweep_gradients=(0d0,0d0)
  sweep_values(1,1)=cos(0.3d0);sweep_values(1,3)=sin(0.3d0)
  sweep_values(3,1)=-sin(0.3d0);sweep_values(3,3)=cos(0.3d0)
  sweep_values(2,2)=cos(0.3d0);sweep_values(2,4)=sin(0.3d0)
  sweep_values(4,2)=-sin(0.3d0);sweep_values(4,4)=cos(0.3d0)
  sweep_gradients(1,:,:)=sweep_values
  call localize_dg_overlapping_wannier_basis(MPI_COMM_WORLD,sweep_values,sweep_gradients,&
    weights,phases,sweep_representation,sweep_product,0.1d0,1d-16,1d-14,1d-12,1,&
    initial_spread,final_spread,maximum_pair_gradient,sweep_iterations,converged,&
    sweep_transform,ok,message)
  call require(.not.ok.and..not.converged.and.index(message,'converge')>0,&
    'nonconverged localization sweep rejects publication')

  sweep_values=(0d0,0d0);sweep_gradients=(0d0,0d0)
  sweep_values(1,1)=cos(0.3d0);sweep_values(1,3)=sin(0.3d0)
  sweep_values(3,1)=-sin(0.3d0);sweep_values(3,3)=cos(0.3d0)
  sweep_values(2,2)=cos(0.3d0);sweep_values(2,4)=sin(0.3d0)
  sweep_values(4,2)=-sin(0.3d0);sweep_values(4,4)=cos(0.3d0)
  sweep_gradients(1,:,:)=sweep_values
  call localize_dg_overlapping_wannier_basis(MPI_COMM_WORLD,sweep_values,sweep_gradients,&
    weights,phases,sweep_representation(:,:,1:1),reshape([1],[1,1]),&
    0.1d0,1d-16,1d-7,1d-12,32,initial_spread,final_spread,maximum_pair_gradient,&
    sweep_iterations,converged,sweep_transform,ok,message)
  call require(ok.and.converged.and.final_spread<initial_spread,&
    'identity-only symmetry permits independent local localization')

  nan_value=ieee_value(0d0,ieee_quiet_nan);weights(2)=nan_value
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'finite')>0,'nonfinite weights rejected')
  weights=1d0;values(1,:)=0d0
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'norm')>0,'zero-norm Wannier rejected')
  values=(0d0,0d0);values(1,1)=1d0;values(2,2)=1d0;phases(1,1)=2d0
  call evaluate_dg_periodic_localization(values,weights,phases,norm,moment,spread,ok,message)
  call require(.not.ok.and.index(message,'phase')>0,'non-unit periodic phase rejected')

  if(rank==0)write(*,'(a,i0,a)')'PASS buffer-local periodic Wannier spread on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_failure,global_failure
    local_failure=merge(0,1,condition)
    call MPI_Allreduce(local_failure,global_failure,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_failure/=0)error stop label
  end subroutine require
end program test_dg_overlapping_wannier_localization_mpi
