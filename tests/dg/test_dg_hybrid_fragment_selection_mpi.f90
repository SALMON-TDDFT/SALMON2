program test_fragment_selection
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_fragment_selection
  use dg_hybrid_fragment_wannier
  use dg_hybrid_fragment_wannier_test_stubs
  implicit none
  integer::rank,np,ierr,f,j,a,nx,variant,owner,expected(4),nexpected,setup_saved,run_saved
  integer::raw_shape(3),core_shape(3),total_shape(3),mapping(8,3),bad_mapping(8,3)
  real(real64)::lattice(3,3),reciprocal(3,3),total_lattice(3,3),origin(3),raw_origin(3)
  real(real64),allocatable::lower(:,:),extent(:,:)
  real(real64)::fractional(3,8),weights(8),energies(2),occupations(2),atoms(3,1)
  real(real64)::points(3,8),bad_points(3,8),local_lower(3),local_extent(3),probe(3)
  complex(real64)::seeds(2,8),buffer(1,8),projector(1,8)
  complex(real64),allocatable::selected_values(:,:)
  integer(int64)::ids(8),first_fp
  type(s_dg_hybrid_fragment_wannier_cache)::cache,bad_cache,snapshot,empty_cache
  complex(real64)::values(8,8),rotated(8,8),physical_projector(8,8),rotated_projector(8,8)
  integer::permutation(8),k
  integer::dims(3),xyz(3)
  type(s_dg_hybrid_core_selection)::selection,again
  logical::ok
  character(256)::message
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  f=np-rank;nx=np*(np+1)/2
  allocate(lower(3,np),extent(3,np))
  origin=[3d0,-5d0,7d0];lower=spread(origin,2,np);extent=8d0
  do j=1,np
    lower(1,j)=origin(1)+real(j*(j-1)/2,real64);extent(1,j)=real(j,real64)
  enddo
  raw_origin=lower(:,f);raw_shape=[8,1,1];core_shape=[f,1,1];total_shape=[nx,1,1]
  lattice=0d0;reciprocal=0d0;total_lattice=0d0
  do a=1,3
    lattice(a,a)=8d0;reciprocal(a,a)=2d0*acos(-1d0)/8d0;total_lattice(a,a)=8d0
  enddo
  total_lattice(1,1)=real(nx,real64)
  mapping=1
  do j=1,8;mapping(j,1)=1+modulo(f*(f-1)/2+j-1,nx);enddo
  seeds=0d0;seeds(1,1)=1d0;seeds(2,2)=1d0;buffer=0d0;buffer(1,8)=1d0
  projector=0d0;projector(1,7)=1d0;fractional=0d0;weights=1d0;atoms=0d0
  energies=[-2d0,-1d0];occupations=[2d0,1d0]
  do j=1,8;ids(j)=j;fractional(1,j)=real(j-1,real64)/8d0;enddo
  call reset_w90_stub_state;expected_fragment_id=f
  call build_dg_hybrid_fragment_wannier(MPI_COMM_WORLD,MPI_COMM_SELF,f,1,'selection-cache',&
    ids,weights,seeds,energies,occupations,buffer,projector,1d-12,lattice,reciprocal,&
    ['H '],atoms,fractional,20,1d-10,10000000_int64,cache,ok,message)
  call require(ok,'raw cache build: '//trim(message));snapshot=cache
  setup_saved=setup_calls;run_saved=run_calls
  call invoke(cache,mapping)
  call require(ok,'selection from real cache: '//trim(message))
  nexpected=0
  ! Stub final x centers are 0.1,0.2,0.3,0.875 in the eight-unit raw cell.
  do j=1,4
    probe=raw_origin;probe(1)=probe(1)+8d0*cache%centers_fractional(1,j)
    probe(1)=origin(1)+modulo(probe(1)-origin(1),real(nx,real64))
    if(probe(1)>=lower(1,f).and.probe(1)<lower(1,f)+extent(1,f))then
      nexpected=nexpected+1;expected(nexpected)=j
    endif
  enddo
  call require(selection%valid.and.selection%raw_count==4.and.selection%selected_count==nexpected,&
    'wrong selection counts')
  call require(all(selection%raw_column_ids==int(expected(:nexpected),int64)), 'wrong raw selected IDs')
  call require(selection%fingerprint/=0_int64.and.selection%raw_cache_fingerprint==&
    cache%receipt%replicated_payload_fingerprint,'missing selection provenance')
  first_fp=selection%fingerprint
  call invoke(cache,mapping)
  call require(ok.and.selection%fingerprint==first_fp,'selection not deterministic')
  call require(all(cache%wannier_values==snapshot%wannier_values).and.&
    all(cache%centers_fractional==snapshot%centers_fractional),'selection mutated raw cache')
  call require(setup_calls==setup_saved.and.run_calls==run_saved,'selection reran W90')

  call export_dg_hybrid_selected_wannier(MPI_COMM_WORLD,f,cache,selection,selected_values,ok,message)
  call require(ok,'selected value export: '//trim(message))
  call require(all(shape(selected_values)==[nexpected,8]),'selected export dropped buffer rows')
  call require(all(selected_values==cache%wannier_values(expected(:nexpected),:)),&
    'selected export changed raw values/order')
  do variant=1,7
    again=selection;bad_cache=cache
    if(rank==0)then
      select case(variant)
      case(1);again%fingerprint=ieor(again%fingerprint,1_int64)
      case(2);again%basis_generation=again%basis_generation+1
      case(3);again%raw_column_ids(1)=0_int64
      case(4);again%center_owner(1)=0
      case(5);bad_cache%wannier_values(1,1)=bad_cache%wannier_values(1,1)+0.1d0
      case(6);again%geometry_fingerprint=ieor(again%geometry_fingerprint,1_int64)
      case(7);deallocate(again%raw_column_ids)
      end select
    endif
    call export_dg_hybrid_selected_wannier(MPI_COMM_WORLD,f,bad_cache,again,selected_values,ok,message)
    call require(.not.ok.and..not.allocated(selected_values),'corrupt selection published WF values')
  enddo
  call require(all(cache%wannier_values==snapshot%wannier_values),'export mutated raw cache')
  call require(setup_calls==setup_saved.and.run_calls==run_saved,'export reran W90')

  do variant=1,4
    bad_cache=cache;bad_mapping=mapping
    if(rank==0)then
      select case(variant)
      case(1);bad_cache%centers_fractional(1,1)=ieee_value(0d0,ieee_quiet_nan)
      case(2);bad_cache%centers_fractional(1,1)=0.9d0
      case(3);bad_mapping(1,1)=0
      case(4);bad_mapping(2,1)=1+modulo(bad_mapping(2,1),nx)
      end select
    endif
    if(variant==4.and.nx==1)cycle
    call invoke(bad_cache,bad_mapping)
    call require(.not.ok.and..not.selection%valid.and..not.allocated(selection%raw_column_ids),&
      'invalid input published selection')
  enddo

  ! Exercise the same ownership kernel directly with face, edge and corner
  ! probes; this does not fabricate a trusted raw-cache receipt.
  points=0d0
  points(:,1)=origin
  points(:,2)=origin+[real(nx,real64),8d0,8d0]
  points(:,3)=origin+[real(nx,real64)-epsilon(1d0),-epsilon(1d0),0d0]
  points(:,4)=lower(:,np)
  points(:,5)=lower(:,np)+[0d0,8d0,8d0]
  points(:,6)=lower(:,np)+[0.25d0,0d0,0d0]
  points(:,7)=points(:,6)
  points(:,8)=origin+[0.25d0,4d0,4d0]
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,points,selection,ok,message)
  call require(ok,'ownership probes: '//trim(message))
  call require(all(selection%center_owner==[1,1,1,np,np,np,np,1]),'half-open periodic ownership incorrect')
  values=0d0;physical_projector=0d0;rotated_projector=0d0
  permutation=[8,3,5,1,7,2,6,4]
  do j=1,8
    do k=1,8;values(k,j)=cmplx(sin(real(j*k,real64)),cos(real(j+2*k,real64)),real64);enddo
    if(selection%center_owner(j)==f)then
      do k=1,8;physical_projector(:,k)=physical_projector(:,k)+values(:,j)*conjg(values(k,j));enddo
    endif
  enddo
  do j=1,8;rotated(:,j)=values(:,permutation(j))*cmplx(cos(real(j,real64)),sin(real(j,real64)),real64);enddo
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,points(:,permutation),again,ok,message)
  call require(ok.and.all(again%center_owner==selection%center_owner(permutation)),&
    'column permutation changed physical center ownership')
  do j=1,8
    if(again%center_owner(j)/=f)cycle
    do k=1,8;rotated_projector(:,k)=rotated_projector(:,k)+rotated(:,j)*conjg(rotated(k,j));enddo
  enddo
  call require(maxval(abs(rotated_projector-physical_projector))<1d-12,&
    'phase/permutation changed selected physical density kernel')

  ! A second valid construction epoch puts all centers at the same buffer point.
  ! The selector must return an empty, valid catalog for fragment one (np>1).
  override_test_centers=.true.
  call build_dg_hybrid_fragment_wannier(MPI_COMM_WORLD,MPI_COMM_SELF,f,2,'selection-cache',&
    ids,weights,seeds,energies,occupations,buffer,projector,1d-12,lattice,reciprocal,&
    ['H '],atoms,fractional,20,1d-10,10000000_int64,empty_cache,ok,message)
  call require(ok,'empty-selection raw cache construction: '//trim(message))
  override_test_centers=.false.;setup_saved=setup_calls;run_saved=run_calls
  call invoke(empty_cache,mapping)
  call require(ok.and.selection%valid,'empty or coincident-center selection failed')
  if(np>1.and.f==1)then
    call require(selection%selected_count==0.and.size(selection%raw_column_ids)==0,&
      'outside centers admitted to core because their WFs have core tails')
  else
    call require(selection%selected_count==0.or.selection%selected_count==4,&
      'coincident centers were split or deduplicated')
  endif
  call require(setup_calls==setup_saved.and.run_calls==run_saved,'empty selection reran W90')
  call export_dg_hybrid_selected_wannier(MPI_COMM_WORLD,f,empty_cache,selection,selected_values,ok,message)
  if(np>1)then
    call require(.not.ok.and..not.allocated(selected_values).and.index(message,'PW-only')>0,&
      'empty selection silently enabled PW-only path')
  else
    call require(ok,'all-retained single fragment export failed')
  endif
  bad_points=points
  if(rank==0)bad_points(1,1)=ieee_value(0d0,ieee_quiet_nan)
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,bad_points,again,ok,message)
  call require(.not.ok.and..not.again%valid,'NaN center accepted')
  if(np>1)then
    local_lower=lower(:,1);lower(:,1)=lower(:,2)
    call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
      total_shape,points,again,ok,message)
    call require(.not.ok,'overlapping/missing core partition accepted');lower(:,1)=local_lower
    if(rank==0)lower(1,1)=lower(1,1)+1d0
    call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
      total_shape,points,again,ok,message)
    call require(.not.ok,'rank-disagreeing geometry accepted')
  endif
  ! A genuinely 3-D partition: 2x2x2 on eight ranks, reduced for smaller runs.
  dims=[min(2,np),max(1,min(2,np/2)),max(1,np/4)]
  total_shape=4*dims;total_lattice=0d0
  do a=1,3;total_lattice(a,a)=real(total_shape(a),real64);enddo
  origin=[3d0,-5d0,7d0];extent=4d0
  do j=1,np
    xyz=[modulo(j-1,dims(1)),modulo((j-1)/dims(1),dims(2)),(j-1)/(dims(1)*dims(2))]
    lower(:,j)=origin+4d0*real(xyz,real64);points(:,j)=lower(:,j)+0.25d0
  enddo
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,points(:,:np),again,ok,message)
  call require(ok.and.all(again%center_owner==[(j,j=1,np)]),'3-D core interior ownership incorrect')
  points(:,1)=origin+4d0;points(:,2)=origin+real(total_shape,real64)
  points(:,3)=origin+4d0-4d0*epsilon(1d0)
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,points(:,:3),again,ok,message)
  call require(ok.and.all(again%center_owner==[np,1,np]),'3-D corner wrapping or snapping incorrect')

  ! Finite but unresolved geometry must return a collective error, not overflow.
  total_lattice=0d0
  do a=1,3;total_lattice(a,a)=1d-154;enddo
  origin=1d154;lower=-1d154;extent=1d-154;total_shape=1;points=0d0
  call classify_dg_hybrid_core_centers(MPI_COMM_WORLD,total_lattice,origin,lower,extent,&
    total_shape,points,again,ok,message)
  call require(.not.ok.and..not.again%valid,'unresolved finite geometry accepted')
  if(rank==0)write(*,'(a,i0,a)')'PASS core-center selection on ',np,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine invoke(input_cache,input_mapping)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::input_cache
    integer,intent(in)::input_mapping(:,:)
    call select_dg_hybrid_core_wannier(MPI_COMM_WORLD,f,input_cache,lattice,raw_origin,&
      total_lattice,origin,lower,extent,raw_shape,core_shape,total_shape,input_mapping,&
      selection,ok,message)
  end subroutine
  subroutine require(condition,description)
    logical,intent(in)::condition
    character(*),intent(in)::description
    integer::bad,global_bad
    bad=merge(0,1,condition)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(global_bad/=0)then
      if(.not.condition)write(*,'(a)')trim(description)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine
end program
