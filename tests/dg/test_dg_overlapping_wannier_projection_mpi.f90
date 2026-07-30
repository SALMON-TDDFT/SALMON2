#include "config.h"
program test_dg_overlapping_wannier_projection_mpi
  use mpi
  use dg_overlapping_wannier_projection,only:t_dg_projection_channel,&
    build_dg_complete_sp_manifest,dg_complete_sp_target_count,validate_dg_complete_sp_manifest,&
    evaluate_dg_periodic_sp_projectors,dg_periodic_grid_point_owned,&
    select_dg_sp_projector_ordinals,select_dg_sp_atomic_orbital_ordinals
  implicit none
  type(t_dg_projection_channel),allocatable::channels(:)
  integer::ierr,rank,target
  real(8)::lattice(3,3),inverse(3,3),points(3,4),translated(3,4),atoms(3,1)
  real(8)::radial_grid(4,2,1),radial_value(4,2,1)
  integer::radial_count(2,1)
  integer::nproj(0:1,2),inorm(0:3,2),ordinals(2,2)
  logical::atomic_orbital_fallback(2,2)
  real(8),allocatable::projector(:,:),translated_projector(:,:)
  logical::ok
  character(256)::message

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)

  call build_dg_complete_sp_manifest([7,11],channels,ok,message)
  call require(ok,trim(message))
  call require(size(channels)==8,'two complete s+p shells')
  call require(all(channels%atom==[7,7,7,7,11,11,11,11]),'atom-major ownership')
  call require(all(channels%l==[0,1,1,1,0,1,1,1]),'complete s+p angular ordering')
  call require(all(channels%m==[1,1,2,3,1,1,2,3]),'complete real-harmonic ordering')
  call validate_dg_complete_sp_manifest(channels,[7,11],ok,message)
  call require(ok,trim(message))
  channels(4)%m=2
  call validate_dg_complete_sp_manifest(channels,[7,11],ok,message)
  call require(.not.ok,'missing p member rejected')
  channels(4)%m=3
  call validate_dg_complete_sp_manifest(channels(:7),[7,11],ok,message)
  call require(.not.ok,'truncated shell rejected')
  call build_dg_complete_sp_manifest([7,7],channels,ok,message)
  call require(.not.ok,'duplicate core-owned atom rejected')

  call dg_complete_sp_target_count(2,target,ok,message)
  call require(ok.and.target==8,'complete s+p target count')
  call dg_complete_sp_target_count(0,target,ok,message)
  call require(.not.ok,'empty core atom set rejected')
  call require(dg_periodic_grid_point_owned([3,0,0],[0,0,0],[4,4,4],[8,8,8]),&
    'last zero-based cell remains in first fragment')
  call require(.not.dg_periodic_grid_point_owned([4,0,0],[0,0,0],[4,4,4],[8,8,8]),&
    'half-open fragment boundary belongs to next fragment')
  call require(dg_periodic_grid_point_owned([7,0,0],[4,0,0],[4,4,4],[8,8,8]),&
    'translated fragment owns its last cell')
  nproj(:,1)=[2,1];nproj(:,2)=[1,3]
  call select_dg_sp_atomic_orbital_ordinals([1,1],nproj,ordinals,ok,message)
  call require(ok.and.all(ordinals(:,1)==[0,2]).and.all(ordinals(:,2)==[0,1]),&
    'pseudo-atomic s+p selector honors flattened multi-projector channel offsets')
  nproj=2;inorm=0
  inorm(:,1)=[0,1,0,1];inorm(:,2)=[1,0,1,0]
  call select_dg_sp_projector_ordinals([1,1],nproj,inorm,ordinals,ok,message,&
    atomic_orbital_fallback)
  call require(ok.and.all(ordinals(:,1)==[1,3]).and.all(ordinals(:,2)==[0,2]),&
    'mixed-species selector skips zero-norm first projectors')
  inorm(2:3,1)=0
  call select_dg_sp_projector_ordinals([1,1],nproj,inorm,ordinals,ok,message)
  call require(.not.ok,'local-reference channel requires an explicit atomic-orbital fallback')
  call select_dg_sp_projector_ordinals([1,1],nproj,inorm,ordinals,ok,message,&
    atomic_orbital_fallback)
  call require(ok.and.ordinals(2,1)==2.and.atomic_orbital_fallback(2,1),&
    'local-reference p channel uses pseudopotential atomic-orbital radial data')

  lattice=0d0;inverse=0d0
  lattice(1,1)=4d0;lattice(2,2)=5d0;lattice(3,3)=6d0
  inverse(1,1)=0.25d0;inverse(2,2)=0.2d0;inverse(3,3)=1d0/6d0
  atoms(:,1)=[3.8d0,0.2d0,5.7d0]
  radial_grid(:,1,1)=[0d0,0.5d0,1d0,2d0]
  radial_grid(:,2,1)=[0d0,0.25d0,0.75d0,2d0]
  radial_value(:,1,1)=[1d0,0.8d0,0.3d0,0d0]
  radial_value(:,2,1)=[0d0,2d0,1d0,0d0]
  radial_count=4
  points(:,1)=[0.1d0,0.2d0,5.7d0]
  points(:,2)=[3.8d0,0.5d0,5.7d0]
  points(:,3)=[3.8d0,0.2d0,0.1d0]
  points(:,4)=atoms(:,1)
  call build_dg_complete_sp_manifest([1],channels,ok,message)
  call evaluate_dg_periodic_sp_projectors(lattice,inverse,points,atoms,radial_grid,&
    radial_value,radial_count,918273645_8,channels,projector,ok,message)
  call require(ok,trim(message))
  call require(all(shape(projector)==[4,4]),'periodic s+p projector shape')
  call require(abs(projector(3,1)-1.9d0)<1d-13,&
    'mixed local-reference and nonlocal channels retain distinct radial grids')
  call require(abs(projector(2,1))<1d-14.and.abs(projector(3,1))>1d-6,&
    'Wannier90 real p(z,x,y) ordering')
  call require(maxval(abs(projector(2:4,4)))<1d-14,'p projectors vanish at atom')
  translated=points
  translated(1,:)=translated(1,:)+lattice(1,1)
  call evaluate_dg_periodic_sp_projectors(lattice,inverse,translated,atoms,radial_grid,&
    radial_value,radial_count,918273645_8,channels,translated_projector,ok,message)
  call require(ok.and.maxval(abs(projector-translated_projector))<1d-13,&
    'projectors are covariant under lattice translation')
  call evaluate_dg_periodic_sp_projectors(lattice,inverse,points,atoms,radial_grid,&
    radial_value,radial_count,0_8,channels,translated_projector,ok,message)
  call require(.not.ok,'missing pseudopotential provenance rejected')
  call evaluate_dg_periodic_sp_projectors(lattice,inverse,points,atoms,radial_grid,&
    radial_value,radial_count,918273645_8,channels(:3),translated_projector,ok,message)
  call require(.not.ok,'incomplete shell rejected at projector boundary')

  if(rank==0)write(*,'(a)')'PASS complete-s+p projection manifest'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,text)
    logical,intent(in)::condition
    character(*),intent(in)::text
    if(.not.condition)then
      write(0,'(a,i0,2a)')'rank ',rank,': ',trim(text)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine
end program
