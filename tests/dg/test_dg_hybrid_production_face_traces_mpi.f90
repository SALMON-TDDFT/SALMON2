#include "config.h"
program test_dg_hybrid_production_face_traces_mpi
  use mpi,only:MPI_Allreduce,MPI_Comm_rank,MPI_Comm_size,MPI_COMM_WORLD,MPI_Finalize,MPI_Init,&
    MPI_DOUBLE_PRECISION,MPI_INTEGER,MPI_MAX,MPI_SUCCESS,MPI_SUM
  use mpi,only:MPI_IN_PLACE
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_production_face_traces,only:s_dg_hybrid_production_face_trace,&
    build_dg_hybrid_production_face_trace,assemble_dg_hybrid_production_face,&
    validate_dg_hybrid_production_face_collection,materialize_dg_hybrid_production_face_collection,&
    assemble_dg_hybrid_production_interface_rows,assemble_dg_hybrid_production_interface_component_rows,&
    reconstruct_dg_hybrid_production_interface_actions,freeze_dg_hybrid_basis_directory,&
    materialize_dg_hybrid_production_interior,&
    reconstruct_dg_hybrid_production_interface_state
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_sipg_operator,only:s_dg_hybrid_sipg_face_operator
  use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows
  implicit none
  integer::icomm,id_rank,nproc,ierr,request_count,response_count,global_request_count,global_response_count,&
    expected_request_count
  integer,allocatable::basis_owner(:),basis_fragment(:)
  integer(int64)::point_ids(2)
  real(real64)::weights(2),normal(3)
  real(real64)::local_interface_norm,global_interface_norm,volume_diagnostics(4)
  real(real64),allocatable::interior_weights(:),interior_potential(:)
  complex(real64)::minus_values(2,2),plus_values(2,1),minus_outward(2,2),plus_outward(2,1)
  complex(real64),allocatable::empty_values(:,:),empty_derivatives(:,:)
  complex(real64),allocatable::interface_rows(:,:),component_rows(:,:,:),component_actions(:,:,:)
  complex(real64),allocatable::interior_values(:,:),interior_gradients(:,:,:)
  complex(real64),allocatable::interior_kinetic_action(:,:)
  complex(real64),allocatable::production_kinetic(:,:),production_local(:,:)
  complex(real64),allocatable::occupied_coefficients(:,:),rotated_coefficients(:,:),interface_state(:,:),&
    rotated_interface_state(:,:)
  integer,allocatable::empty_ids(:)
  type(s_dg_hybrid_production_face_trace)::trace,bad_trace
  type(s_dg_hybrid_production_face_trace),allocatable::production_faces(:)
  type(s_dg_hybrid_fragment_basis)::fragment_bases(2)
  type(s_dg_hybrid_sipg_face_operator)::face
  integer::origins(3,2),sizes(3,2),grid_size(3),p,x,y
  integer(int64)::all_point_ids(8)
  integer(int64),allocatable::owned_row_ids(:)
  integer(int64),allocatable::interior_ids(:)
  integer,allocatable::interior_fragment(:)
  complex(real64)::analytic_values(8,2)
  logical::ok,participant_checks_ok,volume_blocks_ok,interface_state_ok
  complex(real64)::occupied_rotation(2,2)
  character(256)::message

  call MPI_Init(ierr);icomm=MPI_COMM_WORLD
  call MPI_Comm_rank(icomm,id_rank,ierr);call MPI_Comm_size(icomm,nproc,ierr)
  point_ids=[101_int64,109_int64];weights=[0.7d0,1.1d0];normal=[1d0,0d0,0d0]
  minus_values=reshape([cmplx(1d0,0.1d0,real64),cmplx(0.8d0,-0.2d0,real64),&
    cmplx(-0.3d0,0.2d0,real64),cmplx(0.4d0,0.1d0,real64)],[2,2])
  plus_values(:,1)=[cmplx(0.6d0,-0.1d0,real64),cmplx(0.5d0,0.3d0,real64)]
  minus_outward=0.25d0*minus_values;plus_outward(:,1)=-0.4d0*plus_values(:,1)

  call build_dg_hybrid_production_face_trace(icomm,17,1,2,[1,0,0],normal,0.8d0,point_ids,point_ids,&
    weights,[1,2],[3],minus_values,minus_outward,plus_values,plus_outward,trace,ok,message)
  call require(ok,trim(message))
  call require(trace%frozen.and.trace%fingerprint/=0_int64,'production face was not frozen')
  call require(trace%minus_fragment==1.and.trace%plus_fragment==2,'canonical fragment orientation changed')
  call require(all(trace%periodic_shift==[1,0,0]),'periodic image shift was not retained')
  call require(maxval(abs(trace%derivative_minus-minus_outward))<1d-14,&
    'minus derivative is not expressed in the canonical normal')
  call require(maxval(abs(trace%derivative_plus+plus_outward))<1d-14,&
    'plus outward derivative was not converted to the canonical normal')
  call assemble_dg_hybrid_production_face(icomm,trace,6d0,face,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs(face%total(1:2,3)))>1d-12,'cross-fragment SIPG block is zero')
  call require(maxval(abs(face%total-conjg(transpose(face%total))))<1d-13,'production SIPG block is not Hermitian')
  call validate_dg_hybrid_production_face_collection(icomm,[trace,trace],ok,message)
  call require(.not.ok,'duplicate physical face was accepted')

  call build_dg_hybrid_production_face_trace(icomm,19,1,2,[0,0,0],normal,0.8d0,point_ids,[102_int64,110_int64],&
    weights,[1,2],[3],minus_values,minus_outward,plus_values,plus_outward,bad_trace,ok,message)
  call require(ok,'paired cell-centered face points were rejected')
  call build_dg_hybrid_production_face_trace(icomm,21,1,2,[0,0,0],normal,0.8d0,point_ids,[101_int64],&
    weights,[1,2],[3],minus_values,minus_outward,plus_values,plus_outward,bad_trace,ok,message)
  call require(.not.ok,'nonconformable face point arrays were accepted')
  allocate(empty_ids(0),empty_values(2,0),empty_derivatives(2,0))
  call build_dg_hybrid_production_face_trace(icomm,23,1,2,[0,0,0],normal,0.8d0,point_ids,point_ids,&
    weights,empty_ids,[1],empty_values,empty_derivatives,plus_values,plus_outward,bad_trace,ok,message)
  call require(.not.ok,'one-sided production face basis was accepted')

  grid_size=[4,2,1];origins=reshape([0,0,0,1,0,0],[3,2]);sizes=reshape([1,2,1,3,2,1],[3,2])
  p=0
  do y=0,1;do x=0,3
    p=p+1;all_point_ids(p)=int(1+x+4*y,int64)
    analytic_values(p,1)=cmplx(sin(0.5d0*acos(-1d0)*real(x,real64))+0.2d0*y,0.1d0*x,real64)
    analytic_values(p,2)=cmplx(cos(0.5d0*acos(-1d0)*real(x,real64))-0.1d0*y,-0.2d0*x,real64)
  enddo;enddo
  do p=1,2
    fragment_bases(p)%fragment_id=merge(p,0,mod(p-1,nproc)==id_rank);fragment_bases(p)%generation=1
    allocate(fragment_bases(p)%global_ids(merge(merge(2,1,p==1),0,mod(p-1,nproc)==id_rank)),&
      fragment_bases(p)%sector(merge(merge(2,1,p==1),0,mod(p-1,nproc)==id_rank)),&
      fragment_bases(p)%buffer_point_ids(8),&
      fragment_bases(p)%buffer_values(8,merge(merge(2,1,p==1),0,mod(p-1,nproc)==id_rank)))
    if(mod(p-1,nproc)==id_rank)then
      if(p==1)then
        fragment_bases(p)%global_ids=[1_int64,3_int64];fragment_bases(p)%sector=[1,1]
        fragment_bases(p)%buffer_values(:,2)=2d0*analytic_values(:,p)
      else
        fragment_bases(p)%global_ids=[2_int64];fragment_bases(p)%sector=[1]
      endif
      fragment_bases(p)%buffer_values(:,1)=analytic_values(:,p)
    endif
    fragment_bases(p)%buffer_point_ids=all_point_ids
    fragment_bases(p)%provenance_fingerprint=int(100+p,int64)
  enddo
  call freeze_dg_hybrid_basis_directory(icomm,fragment_bases,[1,2,3],basis_owner,basis_fragment,ok,message)
  call require(ok,trim(message))
  call require(all(basis_owner==[mod(0,nproc),mod(1,nproc),mod(0,nproc)]).and.&
    all(basis_fragment==[1,2,1]),&
    'frozen basis owner or fragment directory is incorrect')
  allocate(interior_ids(count([(mod(p-1,nproc)==id_rank,p=1,8)])),&
    interior_fragment(count([(mod(p-1,nproc)==id_rank,p=1,8)])))
  interior_ids=pack(all_point_ids,[(mod(p-1,nproc)==id_rank,p=1,8)])
  do p=1,size(interior_ids)
    interior_fragment(p)=merge(1,2,modulo(int(interior_ids(p)-1_int64),4)==0)
  enddo
  call materialize_dg_hybrid_production_interior(icomm,grid_size,reshape([0.5d0,0.25d0,0.125d0],[1,3]),&
    0.75d0,reshape([0.5d0,0.25d0,0.125d0],[1,3]),fragment_bases,basis_owner,basis_fragment,&
    [1,2,3],interior_ids,interior_fragment,interior_values,interior_gradients,interior_kinetic_action,&
    ok,message,request_count,response_count)
  call require(ok,trim(message))
  call MPI_Allreduce(request_count,global_request_count,1,MPI_INTEGER,MPI_SUM,icomm,ierr)
  call MPI_Allreduce(response_count,global_response_count,1,MPI_INTEGER,MPI_SUM,icomm,ierr)
  expected_request_count=0
  if(any(interior_fragment==1).and.id_rank/=basis_owner(1))expected_request_count=expected_request_count+1
  if(any(interior_fragment==2).and.id_rank/=basis_owner(2))expected_request_count=expected_request_count+1
  call MPI_Allreduce(MPI_IN_PLACE,expected_request_count,1,MPI_INTEGER,MPI_SUM,icomm,ierr)
  call require(global_request_count==expected_request_count.and.global_response_count==expected_request_count,&
    'grouped materialization did not exchange each destination point-ID set exactly once')
  do p=1,size(interior_ids)
    x=modulo(int(interior_ids(p)-1_int64),4);y=int((interior_ids(p)-1_int64)/4_int64)
    call require(abs(interior_values(interior_fragment(p),p)-analytic_values(int(interior_ids(p)),&
      interior_fragment(p)))<1d-13,'owned interior basis value is incorrect')
    call require(abs(interior_values(3-interior_fragment(p),p))<1d-14,&
      'foreign-fragment basis leaked into broken volume')
    call require(abs(interior_gradients(1,interior_fragment(p),p)-0.5d0*(&
      analytic_values(1+modulo(x+1,4)+4*y,interior_fragment(p))-&
      analytic_values(1+modulo(x-1,4)+4*y,interior_fragment(p))))<1d-13,&
      'owned interior basis gradient is incorrect')
    call require(abs(interior_kinetic_action(interior_fragment(p),p)-(&
      0.75d0*analytic_values(int(interior_ids(p)),interior_fragment(p))-0.5d0*(&
      0.5d0*(analytic_values(1+modulo(x+1,4)+4*y,interior_fragment(p))+&
      analytic_values(1+modulo(x-1,4)+4*y,interior_fragment(p)))+&
      0.5d0*analytic_values(1+x+4*modulo(y+1,2),interior_fragment(p))+&
      0.25d0*analytic_values(int(interior_ids(p)),interior_fragment(p)))))<1d-13,&
      'owned interior strong kinetic action is incorrect')
    call require(abs(interior_kinetic_action(3-interior_fragment(p),p))<1d-14,&
      'foreign-fragment kinetic action leaked into broken volume')
  enddo
  call require(maxval(abs(interior_values(3,:)-2d0*interior_values(1,:)))<1d-13,&
    'multiple columns in one fragment were not returned in one grouped response')
  allocate(interior_weights(size(interior_ids)),interior_potential(size(interior_ids)))
  interior_weights=0.25d0;interior_potential=1d0
  allocate(owned_row_ids(count([(mod(p-1,nproc)==id_rank,p=1,2)])))
  owned_row_ids=pack([1_int64,2_int64],[(mod(p-1,nproc)==id_rank,p=1,2)])
  call assemble_dg_hybrid_broken_volume_rows(icomm,2,owned_row_ids,basis_fragment(1:2),interior_ids,&
    interior_fragment,interior_weights,interior_values(1:2,:),interior_gradients(:,1:2,:),interior_potential,&
    production_kinetic,production_local,volume_diagnostics,ok,message)
  call require(ok,trim(message))
  volume_blocks_ok=.true.
  do p=1,size(owned_row_ids)
    volume_blocks_ok=volume_blocks_ok.and.abs(production_kinetic(p,3-int(owned_row_ids(p))))<1d-14.and.&
      abs(production_local(p,3-int(owned_row_ids(p))))<1d-14
  enddo
  call require(volume_blocks_ok,'production interior materialization created a volume cross-fragment block')
  if(size(fragment_bases(1)%global_ids)>0)then
    fragment_bases(1)%global_ids=fragment_bases(1)%global_ids(1:1)
    fragment_bases(1)%sector=fragment_bases(1)%sector(1:1)
    fragment_bases(1)%buffer_values=fragment_bases(1)%buffer_values(:,1:1)
  endif
  call freeze_dg_hybrid_basis_directory(icomm,fragment_bases,[1,2],basis_owner,basis_fragment,ok,message)
  call require(ok,trim(message))
  call materialize_dg_hybrid_production_face_collection(icomm,origins,sizes,grid_size,[1d0,2d0,3d0],&
    reshape([0.5d0,0.25d0,0.125d0],[1,3]),fragment_bases,basis_owner,basis_fragment,[1,2],&
    production_faces,ok,message)
  call require(ok,trim(message))
  call require(size(production_faces)==2,'production topology did not group internal and periodic interfaces')
  participant_checks_ok=.true.
  if(any(basis_owner==id_rank))then
    participant_checks_ok=all([(size(production_faces(p)%weights)==2,p=1,size(production_faces))]).and.&
      any([(abs(production_faces(p)%periodic_shift(1))==1,p=1,size(production_faces))]).and.&
      all([(production_faces(p)%frozen,p=1,size(production_faces))]).and.&
      abs(production_faces(1)%derivative_minus(1,1)-&
      0.5d0*(analytic_values(2,1)-analytic_values(4,1)))<1d-13
    call assemble_dg_hybrid_production_face(icomm,production_faces(1),6d0,face,ok,message)
    participant_checks_ok=participant_checks_ok.and.ok.and.abs(face%total(1,2))>1d-12
  endif
  call require(participant_checks_ok,'face participant materialization or SIPG coupling failed')
  allocate(occupied_coefficients(size(owned_row_ids),2),rotated_coefficients(size(owned_row_ids),2))
  occupied_coefficients=(0d0,0d0)
  do p=1,size(owned_row_ids);occupied_coefficients(p,int(owned_row_ids(p)))=1d0;enddo
  occupied_rotation=reshape([cmplx(1d0,0d0,real64),cmplx(0d0,1d0,real64),&
    cmplx(0d0,1d0,real64),cmplx(1d0,0d0,real64)],[2,2])/sqrt(2d0)
  rotated_coefficients=matmul(occupied_coefficients,occupied_rotation)
  call reconstruct_dg_hybrid_production_interface_state(icomm,2,owned_row_ids,occupied_coefficients,&
    [1d0,1d0],production_faces,interface_state,ok,message)
  call require(ok,trim(message))
  call reconstruct_dg_hybrid_production_interface_state(icomm,2,owned_row_ids,rotated_coefficients,&
    [1d0,1d0],production_faces,rotated_interface_state,ok,message)
  call require(ok,trim(message))
  interface_state_ok=all(shape(interface_state)==shape(rotated_interface_state))
  if(size(interface_state)>0)interface_state_ok=interface_state_ok.and.&
    maxval(abs(interface_state-rotated_interface_state))<1d-12
  call require(interface_state_ok,'occupied gauge rotation changed the production interface state')
  call assemble_dg_hybrid_production_interface_rows(icomm,2,owned_row_ids,production_faces,6d0,&
    interface_rows,ok,message)
  call require(ok,trim(message))
  call require(size(interface_rows,1)==size(owned_row_ids).and.size(interface_rows,2)==2,&
    'production interface rows have the wrong distributed shape')
  local_interface_norm=sum(abs(interface_rows)**2)
  call MPI_Allreduce(local_interface_norm,global_interface_norm,1,MPI_DOUBLE_PRECISION,MPI_SUM,icomm,ierr)
  call require(ierr==MPI_SUCCESS.and.global_interface_norm>1d-24,&
    'production interface row assembly lost all SIPG coupling')
  call assemble_dg_hybrid_production_interface_component_rows(icomm,2,owned_row_ids,production_faces,6d0,&
    component_rows,ok,message)
  call require(ok.and.(size(owned_row_ids)==0.or.maxval(abs(interface_rows-sum(component_rows,dim=3)))<1d-13),&
    'SIPG component rows do not reproduce the frozen total interface action')
  call reconstruct_dg_hybrid_production_interface_actions(icomm,2,owned_row_ids,occupied_coefficients,&
    production_faces,6d0,component_actions,ok,message)
  call require(ok.and.all(shape(component_actions)==[size(owned_row_ids),2,3]),&
    'reconstructed SIPG component actions have the wrong distributed shape')
  do p=1,3
    call require(size(owned_row_ids)==0.or.maxval(abs(component_actions(:,:,p)-component_rows(:,:,p)))<1d-13,&
      'reconstructed occupied SIPG action disagrees with its frozen component rows')
  enddo
  origins(1,2)=2;sizes(1,2)=3
  call materialize_dg_hybrid_production_face_collection(icomm,origins,sizes,grid_size,[1d0,2d0,3d0],&
    reshape([0.5d0,0.25d0,0.125d0],[1,3]),fragment_bases,basis_owner,basis_fragment,[1,2],&
    production_faces,ok,message)
  call require(.not.ok,'fragment box extending outside the basic cell was accepted')
  if(id_rank==0)write(*,'(a,i0,a)')'PASS production face traces on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,icomm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(id_rank==0)write(0,'(a)')trim(label)
      error stop 1
    endif
  end subroutine require
end program test_dg_hybrid_production_face_traces_mpi
