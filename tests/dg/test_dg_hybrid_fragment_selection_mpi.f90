program test_fragment_selection
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_hybrid_fragment_selection
  use dg_hybrid_fragment_wannier
  use dg_hybrid_fragment_wannier_test_stubs
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_projected_fragment_pipeline,only:build_dg_hybrid_projected_local_fragment_basis,&
    project_dg_hybrid_core_seeds,s_dg_hybrid_core_projection_report,&
    s_dg_hybrid_support_samples,check_dg_hybrid_seed_support
  use dg_hybrid_broken_volume,only:assemble_dg_hybrid_broken_volume_rows
  use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,&
    initialize_dg_hybrid_fragment_density_checked
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
  type(s_dg_hybrid_selected_catalog)::unequal_catalog
  type(s_dg_hybrid_dc_reference)::dc_oracle
  integer,allocatable::selected_counts(:)
  complex(real64)::density_metric(4,4)
  logical::density_callback_failure=.false.
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
  call prepare_dg_hybrid_selected_catalog(MPI_COMM_WORLD,f,cache,selection,unequal_catalog,ok,message)
  call require(ok.and.unequal_catalog%valid,'unequal selected catalog: '//trim(message))
  allocate(selected_counts(np))
  call MPI_Allgather(nexpected,1,MPI_INTEGER,selected_counts,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
  call require(size(unequal_catalog%wannier_owner)==sum(selected_counts),'unequal global count mismatch')
  call require(count(unequal_catalog%wannier_owner==f)==nexpected,'unequal local count mismatch')
  call require(all(unequal_catalog%local_active_ids==&
    int([(sum(selected_counts(np-f+2:))+j,j=1,nexpected)],int64)),&
    'unequal active ID offsets follow ranks rather than fragments')

  call export_dg_hybrid_selected_wannier(MPI_COMM_WORLD,f,cache,selection,selected_values,ok,message)
  call require(ok,'selected value export: '//trim(message))
  call require(all(shape(selected_values)==[nexpected,8]),'selected export dropped buffer rows')
  call require(all(selected_values==cache%wannier_values(expected(:nexpected),:)),&
    'selected export changed raw values/order')
  call export_dg_hybrid_dc_reference(MPI_COMM_WORLD,f,cache,selection,dc_oracle,ok,message)
  call require(ok.and.dc_oracle%valid,'bound raw DC reference export: '//trim(message))
  call require(all(dc_oracle%core_row_slots==[(j,j=1,f)]),'raw core row slots disagree with DC mapping')
  call require(all(dc_oracle%physical_grid_ids==int(mapping(:,1),int64)),&
    'bound DC reference lost periodic physical grid mapping')
  call require(maxval(abs(dc_oracle%buffer_orbitals-transpose(seeds)))<1d-11,&
    'raw DC oracle was sliced by selected WF IDs')
  call require(maxval(abs(dc_oracle%core_orbitals-transpose(seeds(:,:f))))<1d-11,&
    'raw DC core oracle used the wrong grid rows')
  call require(all(dc_oracle%energies==energies).and.all(dc_oracle%occupations==occupations).and.&
    dc_oracle%selection_fingerprint==selection%fingerprint,'raw DC metadata binding lost')
  do variant=1,12
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
      case(8);again%core_row_slots(1)=8
      case(9);again%physical_grid_ids(1)=0_int64
      case(10);deallocate(again%core_row_slots)
      case(11);bad_cache%dc_seed_coefficients_in_wannier(1,1)=&
        bad_cache%dc_seed_coefficients_in_wannier(1,1)+0.1d0
      case(12);bad_cache%physical_dc_seed_occupations(1)=bad_cache%physical_dc_seed_occupations(1)+0.1d0
      end select
    endif
    call export_dg_hybrid_selected_wannier(MPI_COMM_WORLD,f,bad_cache,again,selected_values,ok,message)
    call require(.not.ok.and..not.allocated(selected_values),'corrupt selection published WF values')
    call export_dg_hybrid_dc_reference(MPI_COMM_WORLD,f,bad_cache,again,dc_oracle,ok,message)
    call require(.not.ok.and..not.dc_oracle%valid.and..not.allocated(dc_oracle%buffer_orbitals),&
      'corrupt selection published raw DC reference')
  enddo
  call require(all(cache%wannier_values==snapshot%wannier_values),'export mutated raw cache')
  call require(setup_calls==setup_saved.and.run_calls==run_saved,'export reran W90')
  call test_permuted_reference()

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
  if(np>=2)call test_selected_pw_connection()
  if(rank==0)write(*,'(a,i0,a)')'PASS core-center selection on ',np,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine test_permuted_reference()
    type(s_dg_hybrid_fragment_wannier_cache)::reordered
    type(s_dg_hybrid_core_selection)::receipt
    type(s_dg_hybrid_dc_reference)::oracle
    integer::order(8),p,setup_before,run_before
    logical::passed
    character(256)::why
    order=[8,1,7,2,6,3,5,4]
    call build_dg_hybrid_fragment_wannier(MPI_COMM_WORLD,MPI_COMM_SELF,f,4,'reference-permuted-cache',&
      ids(order),weights(order),seeds(:,order),energies,occupations,buffer(:,order),projector(:,order),&
      1d-12,lattice,reciprocal,['H '],atoms,fractional(:,order),20,1d-10,10000000_int64,reordered,passed,why)
    call require(passed,'permuted raw DC reference construction: '//trim(why))
    call select_dg_hybrid_core_wannier(MPI_COMM_WORLD,f,reordered,lattice,raw_origin,total_lattice,origin,&
      lower,extent,raw_shape,core_shape,total_shape,mapping,receipt,passed,why)
    call require(passed,'permuted raw DC selection: '//trim(why))
    setup_before=setup_calls;run_before=run_calls
    call export_dg_hybrid_dc_reference(MPI_COMM_WORLD,f,reordered,receipt,oracle,passed,why)
    call require(passed,'permuted raw DC export: '//trim(why))
    call require(all(oracle%core_row_slots==pack([(p,p=1,8)],order<=f)),&
      'reference assumed core rows precede buffer rows')
    call require(all(oracle%physical_grid_ids==int(mapping(order,1),int64)),&
      'reference discarded raw row permutation')
    call require(maxval(abs(oracle%buffer_orbitals-transpose(seeds(:,order))))<1d-11.and.&
      maxval(abs(oracle%core_orbitals-transpose(seeds(:,pack(order,order<=f)))))<1d-11,&
      'reference orbitals do not follow raw core row slots')
    call require(setup_calls==setup_before.and.run_calls==run_before,'permuted reference reran W90')
  end subroutine
  subroutine test_selected_pw_connection()
    type(s_dg_hybrid_fragment_wannier_cache)::raw
    type(s_dg_hybrid_dc_reference)::bound_reference
    type(s_dg_hybrid_core_selection)::chosen,corrupt
    type(s_dg_hybrid_selected_catalog)::active,invalid
    type(s_dg_hybrid_basis_catalog)::packets
    type(s_dg_hybrid_fragment_basis)::selected_basis,unselected_basis
    type(s_dg_hybrid_core_projection_report)::report
    type(s_dg_hybrid_fragment_subspace_state)::initial,saved
    integer,allocatable::chosen_seeds(:)
    real(real64)::reference_density(4),density_errors(2)
    real(real64)::density_limit
    real(real64)::cell(3,3),tot(3,3),recip(3,3),boxes(3,np),widths(3,np),start(3)
    real(real64)::coords(3,8),windows(np,8),g(3,1),core_weights(4),quad(8),frac(3,8)
    complex(real64)::dc(2,8),aux(1,8),proj(1,8),target(4)
    complex(real64)::reference(4,1),dependent(4,4),expected_coeff(4,1),mixed_basis(4,4)
    complex(real64)::dc_reference(4,2)
    complex(real64)::full_mixed_basis(8,4)
    complex(real64)::volume_values(4*np,4),gradients(3,4*np,4),rhs(4,1)
    complex(real64),allocatable::metric_rows(:,:),kinetic_rows(:,:)
    real(real64)::volume_diagnostics(4)
    complex(real64),allocatable::coeff(:,:)
    real(real64)::metric_weights(4),limits(4)
    integer::map(8,3),p,t,owner0,raw_owner(4*np),run0,setup0,local_n,local_m
    integer(int64)::physical(8),core(4),mem,fp,raw_union_fp
    logical::passed
    character(256)::why
    if(np==1)return
    cell=0d0;recip=0d0;tot=0d0;boxes=0d0;widths=8d0;g=0d0;quad=1d0;core_weights=1d0
    do t=1,3;cell(t,t)=8d0;recip(t,t)=2d0*acos(-1d0)/8d0;tot(t,t)=8d0;enddo
    tot(1,1)=4d0*np;widths(1,:)=4d0
    do p=1,np;boxes(1,p)=4d0*(p-1);enddo
    start=boxes(:,f);map=1;coords=0d0;frac=0d0;windows=0d0
    do p=1,8
      map(p,1)=1+modulo(4*(f-1)+p-1,4*np)
      physical(p)=int(map(p,1),int64);coords(1,p)=real(map(p,1)-1,real64)
      frac(1,p)=real(p-1,real64)/8d0
      owner0=1+(map(p,1)-1)/4;windows(owner0,p)=1d0
    enddo
    core=physical(:4);dc=0d0;aux=0d0;proj=0d0
    dc(1,1)=1d0;dc(2,2)=1d0;aux(1,3)=1d0
    ! The excluded raw WF has a finite core tail and an extended-domain part.
    proj(1,4)=sqrt(0.5d0);proj(1,7)=sqrt(0.5d0)
    call build_dg_hybrid_fragment_wannier(MPI_COMM_WORLD,MPI_COMM_SELF,f,3,'selected-pw-cache',&
      ids,quad,dc,energies,occupations,aux,proj,1d-12,cell,recip,&
      ['H '],atoms,frac,20,1d-10,10000000_int64,raw,passed,why)
    call require(passed,'PW fixture raw build: '//trim(why));run0=run_calls;setup0=setup_calls
    call select_dg_hybrid_core_wannier(MPI_COMM_WORLD,f,raw,cell,start,tot,[0d0,0d0,0d0],&
      boxes,widths,[8,1,1],[4,1,1],[4*np,1,1],map,chosen,passed,why)
    call require(passed.and.chosen%selected_count==3,'PW fixture center selection: '//trim(why))
    call prepare_dg_hybrid_selected_catalog(MPI_COMM_WORLD,f,raw,chosen,active,passed,why)
    call require(passed.and.active%valid,'selected active catalog: '//trim(why))
    call require(all(active%local_active_ids==int([(3*(f-1)+p,p=1,3)],int64)),&
      'active IDs must be compact fragment-major, independent of MPI rank order')
    call require(all(active%raw_column_ids==int([(modulo(p-1,3)+1,p=1,3*np)],int64)),&
      'raw provenance IDs were replaced by global active IDs')
    call require(all(active%wannier_owner==[(1+(p-1)/3,p=1,3*np)]),'wrong selected WF inventory')
    call require(all(active%local_values==raw%wannier_values(:3,:)),'selected catalog changed buffer values')
    allocate(packets%packets(np));packets%valid=.true.
    packets%packet_fingerprint=991_int64;packets%catalog_fingerprint=997_int64
    do p=1,np
      packets%packets(p)%fragment_id=p;packets%packets(p)%owner_rank=np-p;packets%packets(p)%star_id=1
      allocate(packets%packets(p)%g_indices(1),source=[1])
    enddo
    call build_dg_hybrid_projected_local_fragment_basis(MPI_COMM_WORLD,4*np,np,f,core,core_weights,&
      coords(:,:4),windows(:,:4),physical,active%local_values,coords,windows,packets,g,&
      active%wannier_owner,2,1d-12,active%fingerprint,selected_basis,mem,fp,passed,why,basis_generation=3)
    call require(passed,'selected WF to production PW pipeline: '//trim(why))
    call require(all(selected_basis%global_ids(:3)==active%local_active_ids),'pipeline lost active IDs')
    call require(maxval(abs(selected_basis%buffer_values(:,:3)-transpose(active%local_values)))<1d-13,&
      'pipeline clipped selected buffer values')
    target=0d0;target(4)=1d0
    call require(maxval(abs(selected_basis%buffer_values(:4,4)-target))<1d-11,&
      'selected-only projection removed the needed core PW component')
    raw_owner=[(1+(p-1)/4,p=1,4*np)]
    call MPI_Allreduce(raw%receipt%distributed_wannier_fingerprint,raw_union_fp,1,&
      MPI_INTEGER8,MPI_BXOR,MPI_COMM_WORLD,ierr)
    if(raw_union_fp==0_int64)raw_union_fp=1_int64
    call build_dg_hybrid_projected_local_fragment_basis(MPI_COMM_WORLD,4*np,np,f,core,core_weights,&
      coords(:,:4),windows(:,:4),physical,raw%wannier_values,coords,windows,packets,g,&
      raw_owner,2,1d-12,raw_union_fp,&
      unselected_basis,mem,fp,passed,why,basis_generation=3)
    call require(passed,'unselected comparison pipeline: '//trim(why))
    call require(maxval(abs(unselected_basis%buffer_values(:4,5)))<1d-11,&
      'raw-union comparison did not expose lost PW component')
    ! Independent oracle: the excluded column's core restriction is proportional
    ! to e4. Slicing its old raw coordinate vector gives zero, whereas the PW
    ! column represents it exactly. This is a projection probe, not a DC run.
    reference(:,1)=raw%wannier_values(4,:4)
    call require(sum(abs(reference)**2)>0.1d0,'excluded WF lacks the required finite core tail')
    metric_weights=[0.5d0,1.25d0,2d0,0.75d0];limits=[1d-12,1d-10,1d-10,1d-10]
    mixed_basis=selected_basis%buffer_values(:4,:)
    ! A nonorthogonal coordinate change makes using an identity S fail.
    mixed_basis(:,4)=2d0*mixed_basis(:,4)+0.3d0*mixed_basis(:,1)
    expected_coeff=0d0;expected_coeff(4,1)=reference(4,1)/2d0
    expected_coeff(1,1)=-0.3d0*expected_coeff(4,1)
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis,metric_weights,reference,&
      [1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(passed.and.report%metric_rank==4.and.report%measured,'core seed projection: '//trim(why))
    call require(maxval(abs(coeff-expected_coeff))<1d-11,'core solve differs from independent coefficient oracle')
    volume_values=0d0;gradients=0d0
    volume_values(selected_basis%global_ids,:)=transpose(mixed_basis)
    call assemble_dg_hybrid_broken_volume_rows(MPI_COMM_WORLD,4*np,selected_basis%global_ids,&
      [active%wannier_owner,[(p,p=1,np)]],core,[f,f,f,f],metric_weights,volume_values,gradients,&
      [1d0,1d0,1d0,1d0],kinetic_rows,metric_rows,volume_diagnostics,passed,why)
    call require(passed,'actual broken-volume core metric: '//trim(why))
    rhs=matmul(conjg(transpose(mixed_basis)),reference*spread(metric_weights,2,1))
    call require(maxval(abs(matmul(metric_rows(:,selected_basis%global_ids),coeff)-rhs))<1d-11,&
      'core projection normal equations disagree with broken-volume metric')
    call require(report%orbital_residual<1d-11.and.report%density_defect<1d-11.and.&
      report%electron_defect<1d-11,'exact core projection changed density or electron count')
    call export_dg_hybrid_dc_reference(MPI_COMM_WORLD,f,raw,chosen,bound_reference,passed,why)
    call require(passed,'production raw DC reference adapter: '//trim(why))
    call require(all(bound_reference%physical_grid_ids==physical).and.&
      bound_reference%selection_fingerprint==chosen%fingerprint,'reference/selected grid binding differs')
    dc_reference=bound_reference%core_orbitals
    call require(maxval(abs(dc_reference-transpose(dc(:,:4))))<1d-11,'raw DC reconstruction oracle failed')
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis,metric_weights,dc_reference,&
      bound_reference%occupations,limits,3,0d0,coeff,report,passed,why)
    call require(passed,'immutable raw DC seed projection: '//trim(why))
    call require(maxval(abs(matmul(mixed_basis,coeff)-dc_reference))<1d-11,'raw DC states were not reproduced')
    full_mixed_basis=selected_basis%buffer_values
    full_mixed_basis(:,4)=2d0*full_mixed_basis(:,4)+0.3d0*full_mixed_basis(:,1)
    call test_support_checks(full_mixed_basis,bound_reference%buffer_orbitals,coeff,raw%wannier_values(4,:))
    density_metric=matmul(conjg(transpose(mixed_basis)),mixed_basis)
    reference_density=0d0
    do p=1,size(bound_reference%occupations)
      reference_density=reference_density+bound_reference%occupations(p)*abs(bound_reference%core_orbitals(:,p))**2
    enddo
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,bound_reference%basis_generation,&
      11_int64,13_int64,coeff,bound_reference%energies,bound_reference%occupations,0,1d-8,1d-10,1d-10,&
      mixed_basis,core_weights,reference_density,1d-10,1d-10,apply_density_metric,initial,&
      chosen_seeds,density_errors,passed,why)
    call require(passed.and.maxval(density_errors)<1d-10,'bound raw DC reference to initializer: '//trim(why))
    ! The half-norm negative below starts without a previously published state.
    deallocate(initial%vectors,initial%directions)
    call require(all(raw%physical_dc_seed_energies==energies).and.&
      all(raw%physical_dc_seed_occupations==occupations),'projection changed raw physical seed metadata')
    local_n=3+modulo(rank,2);local_m=1+modulo(rank,2)
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis(:,:local_n),metric_weights,&
      dc_reference(:,:local_m),raw%physical_dc_seed_occupations(:local_m),limits,3,0d0,&
      coeff,report,passed,why)
    call require(passed.and.report%metric_rank==local_n,'unequal local basis/seed counts failed: '//trim(why))
    ! Insufficient selected-only span: positive S is not enough for admission.
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis(:,:3),metric_weights,reference,&
      [1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(.not.passed.and..not.allocated(coeff).and.report%measured.and.&
      report%orbital_residual>0.9d0.and.index(why,'insufficient core span')>0,&
      'insufficient span was not separately diagnosed')
    dependent=mixed_basis
    if(rank==0)dependent(:,4)=dependent(:,1)
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,dependent,metric_weights,reference,&
      [1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(.not.passed.and..not.allocated(coeff).and.index(why,'core metric')>0,&
      'dependent core metric silently compressed or published coefficients')
    dependent=mixed_basis
    if(rank==0)dependent(1,1)=cmplx(huge(1d0)/2d0,0d0,real64)
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,dependent,metric_weights,reference,&
      [1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(.not.passed.and..not.allocated(coeff),'overflowing finite metric input escaped admission')
    if(rank==0)limits(2)=2d-10
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis,metric_weights,reference,&
      [1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(.not.passed.and.index(why,'controls differ')>0,'rank-disagreeing projection controls accepted')
    limits=[1d-12,1d-10,1d-10,1d-10]
    density_metric=matmul(conjg(transpose(mixed_basis)),mixed_basis)
    ! The exact projected seed has core norm 1/2: no PW-span error exists.
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis,core_weights,reference,&
      [2d0],limits,3,0d0,coeff,report,passed,why)
    call require(passed.and.report%orbital_residual<1d-11,'half-norm seed is not exactly representable')
    reference_density=2d0*abs(reference(:,1))**2
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0],[2d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and.index(why,'post-initializer density mismatch')>0.and.&
      .not.allocated(initial%vectors).and..not.allocated(chosen_seeds),&
      'half-norm initializer published density-changing state or mislabeled span failure')
    call require(abs(density_errors(2)-1d0)<1d-10,'half-norm density did not double as expected')
    ! Unit-core-norm state is admissible; energy/occupation follow returned IDs.
    coeff=coeff*sqrt(2d0);reference_density=reference_density*2d0
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0],[2d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(passed.and.all(chosen_seeds==[1]).and.maxval(density_errors)<1d-10,&
      'density-preserving initialization failed: '//trim(why))
    saved=initial
    if(rank==0)reference_density=reference_density*0.5d0
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0],[2d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and..not.allocated(chosen_seeds).and.&
      all(initial%vectors==saved%vectors).and.all(initial%directions==saved%directions),&
      'one-rank density failure changed an existing solver state')
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,mixed_basis,core_weights,dc_reference,&
      [0.5d0,1.5d0],limits,3,0d0,coeff,report,passed,why)
    call require(passed,'fractionally occupied seed projection failed')
    reference_density=0.5d0*abs(dc_reference(:,1))**2+1.5d0*abs(dc_reference(:,2))**2
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(passed.and.all(chosen_seeds==[2,1]).and.maxval(density_errors)<1d-10,&
      'energy sorting lost the fractional occupation mapping')
    saved=initial
    if(rank==0)reference_density=cshift(reference_density,1)
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and.density_errors(2)<1d-10.and.all(initial%vectors==saved%vectors),&
      'density redistribution was accepted solely because electron count matched')
    reference_density=0.5d0*abs(dc_reference(:,1))**2+1.5d0*abs(dc_reference(:,2))**2
    density_limit=1d-10;if(rank==0)density_limit=2d-10
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      density_limit,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and.index(why,'controls differ')>0,'different density controls accepted')
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],merge(1,0,rank==0),1d-8,1d-10,1d-10,&
      mixed_basis,core_weights,reference_density,1d-10,1d-10,apply_density_metric,initial,&
      chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and.index(why,'controls differ')>0,'different guard counts accepted')
    if(rank==0)then
      call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
        coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
        1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why,energy_cutoff=-1d0)
    else
      call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
        coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
        1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    endif
    call require(.not.passed.and.index(why,'controls differ')>0,'different cutoff presence accepted')
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why,&
      energy_cutoff=merge(-1d0,0d0,rank==0))
    call require(.not.passed.and.index(why,'controls differ')>0,'different energy cutoffs accepted')
    density_callback_failure=rank==0
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and..not.allocated(chosen_seeds).and.all(initial%vectors==saved%vectors),&
      'one-rank metric callback failure published state')
    density_callback_failure=.false.
    if(rank==0)reference_density(1)=ieee_value(0d0,ieee_quiet_nan)
    call initialize_dg_hybrid_fragment_density_checked(MPI_COMM_WORLD,f,3,11_int64,13_int64,&
      coeff,[-1d0,-2d0],[0.5d0,1.5d0],0,1d-8,1d-10,1d-10,mixed_basis,core_weights,reference_density,&
      1d-10,1d-10,apply_density_metric,initial,chosen_seeds,density_errors,passed,why)
    call require(.not.passed.and.all(initial%vectors==saved%vectors),'nonfinite reference density accepted')
    corrupt=chosen;if(rank==0)corrupt%fingerprint=ieor(corrupt%fingerprint,1_int64)
    call prepare_dg_hybrid_selected_catalog(MPI_COMM_WORLD,f,raw,corrupt,invalid,passed,why)
    call require(.not.passed.and..not.invalid%valid.and..not.allocated(invalid%local_values),&
      'invalid selection published active catalog')
    call require(run_calls==run0.and.setup_calls==setup0,'PW catalog reran W90')
  end subroutine
  subroutine test_support_checks(full_basis,raw_seeds,dc_coeff,excluded)
    complex(real64),intent(in)::full_basis(8,4),raw_seeds(:,:),dc_coeff(:,:),excluded(8)
    type(s_dg_hybrid_support_samples)::samples(3),changed(3)
    type(s_dg_hybrid_core_projection_report)::core_report
    complex(real64)::boundary(3,8),derivative(2,8),projector_map(1,8)
    complex(real64),allocatable::probe_coeff(:,:)
    real(real64)::defects(3),tolerances(3),expected_errors(3)
    integer::channel,counts(3)
    logical::passed,measured
    character(256)::why
    boundary=0d0;boundary(1,1)=1d0;boundary(2,4)=1d0;boundary(3,7)=1d0
    ! Explicit small difference stencils and a normalized projector functional;
    ! these test admission arithmetic, not production SIPG assembly (Task C5).
    derivative=0d0;derivative(1,1)=-1d0;derivative(1,2)=1d0
    derivative(2,5)=-0.5d0;derivative(2,7)=0.5d0
    projector_map=0d0;projector_map(1,[1,4,7])=1d0/sqrt(3d0)
    samples(1)%basis=matmul(boundary,full_basis);samples(1)%reference=matmul(boundary,raw_seeds)
    samples(2)%basis=matmul(derivative,full_basis);samples(2)%reference=matmul(derivative,raw_seeds)
    samples(3)%basis=matmul(projector_map,full_basis);samples(3)%reference=matmul(projector_map,raw_seeds)
    counts=[3,2,1];tolerances=1d-10
    do channel=1,3;allocate(samples(channel)%weights(counts(channel)),source=1d0);enddo
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,samples,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(passed.and.measured.and.maxval(defects)<1d-11,'exact DC support reconstruction: '//trim(why))
    do channel=1,3
      changed=samples
      if(rank==0)changed(channel)%reference(1,1)=changed(channel)%reference(1,1)+0.25d0
      call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
        defects,measured,passed,why)
      call require(.not.passed.and.measured.and.index(why,'required support mismatch')>0,&
        'single-rank support defect was not rejected collectively')
      if(rank==0)then
        call require(abs(defects(channel)-0.25d0)<1d-10.and.count(defects>1d-10)==1,&
          'support channels were mixed or hidden')
      else
        call require(maxval(defects)<1d-10,'support report lost local error information')
      endif
    enddo
    changed=samples
    if(rank==0)deallocate(changed(2)%basis)
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and..not.measured,'missing derivative evidence accepted')
    counts(1)=4
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,samples,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and..not.measured,'missing required sample count accepted')
    counts=[3,2,1]
    changed=samples
    if(rank==0)then
      changed(1)%reference(1,1)=changed(1)%reference(1,1)+0.25d0
      changed(1)%weights(1)=4d0
    endif
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and.measured,'weighted support defect accepted')
    if(rank==0)then
      call require(abs(defects(1)-0.5d0)<1d-10,'support norm ignored quadrature weights')
    else
      call require(defects(1)<1d-10,'weighted support report changed another fragment')
    endif
    changed=samples
    if(rank==0)changed(3)%weights(1)=ieee_value(0d0,ieee_quiet_nan)
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and..not.measured,'nonfinite support weight accepted')
    changed=samples
    if(rank==0)changed(1)%reference(1,1)=cmplx(huge(1d0)/2d0,0d0,real64)
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and..not.measured,'overflowing finite support error escaped admission')
    if(rank==0)tolerances(2)=2d-10
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,samples,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(.not.passed.and.index(why,'controls differ')>0,'different support tolerances accepted')
    tolerances=1d-10
    changed=samples
    if(modulo(rank,2)==0)then
      counts(3)=0
      deallocate(changed(3)%basis,changed(3)%reference,changed(3)%weights)
      allocate(changed(3)%basis(0,4),changed(3)%reference(0,size(dc_coeff,2)),changed(3)%weights(0))
    endif
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,dc_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    call require(passed.and.measured.and.maxval(defects)<1d-10,'declared zero projector inventory failed')
    counts=[3,2,1]
    ! Exactly representable core restriction, but the excluded buffer tail is
    ! needed by the declared boundary/stencil/projector inventory.
    call project_dg_hybrid_core_seeds(MPI_COMM_WORLD,full_basis(:4,:),[1d0,1d0,1d0,1d0],&
      reshape(excluded(:4),[4,1]),[1d0],[1d-12,1d-10,1d-10,1d-10],3,0d0,&
      probe_coeff,core_report,passed,why)
    call require(passed,'support-tail probe failed before support check')
    changed=samples
    changed(1)%reference=matmul(boundary,reshape(excluded,[8,1]))
    changed(2)%reference=matmul(derivative,reshape(excluded,[8,1]))
    changed(3)%reference=matmul(projector_map,reshape(excluded,[8,1]))
    call check_dg_hybrid_seed_support(MPI_COMM_WORLD,probe_coeff,changed,counts,tolerances,3,0d0,&
      defects,measured,passed,why)
    expected_errors=[sqrt(0.5d0),sqrt(0.5d0)/2d0,sqrt(0.5d0/3d0)]
    call require(.not.passed.and.measured.and.maxval(abs(defects-expected_errors))<1d-10,&
      'core-exact state concealed required buffer-tail loss')
  end subroutine
  subroutine apply_density_metric(input,output,valid)
    complex(real64),intent(in)::input(:,:)
    complex(real64),intent(out)::output(:,:)
    logical,intent(out)::valid
    output=matmul(density_metric,input);valid=.not.density_callback_failure
  end subroutine
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
