#include "config.h"
module dg_hybrid_fragment_selection
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_wannier,only:s_dg_hybrid_fragment_wannier_cache,&
    export_dg_hybrid_fragment_coordinates,map_dg_hybrid_fragment_dc_grid
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::s_dg_hybrid_core_selection,select_dg_hybrid_core_wannier,classify_dg_hybrid_core_centers
  integer,parameter::center_convention=1
  real(real64),parameter::roundoff_factor=64d0*epsilon(1d0)
  type s_dg_hybrid_core_selection
    logical::valid=.false.
    integer::fragment_id=0,basis_generation=0,raw_count=0,selected_count=0,convention=center_convention
    integer(int64)::raw_cache_fingerprint=0_int64,geometry_fingerprint=0_int64,fingerprint=0_int64
    integer(int64),allocatable::raw_column_ids(:)
    integer,allocatable::center_owner(:)
  end type
contains
  ! Geometry-only entry: classifies caller-local centers; it does not certify a WF cache.
  subroutine classify_dg_hybrid_core_centers(comm,lattice,origin,lower,extent,grid,centers,result,ok,message)
    integer,intent(in)::comm,grid(3)
    real(real64),intent(in)::lattice(3,3),origin(3),lower(:,:),extent(:,:),centers(:,:)
    type(s_dg_hybrid_core_selection),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_core_selection)::work
    integer::np,ierr,a,b,f,g,j,status,n,match
    integer,allocatable::starts(:,:),widths(:,:)
    real(real64),allocatable::metadata(:),reference(:)
    real(real64)::length(3),h(3),lo(3),wide(3),x(3),distance(3),tol(3),diagonal(3,3)
    integer(int64)::volume,total_volume,box_volume
    logical::valid,overlap
    call MPI_Comm_size(comm,np,ierr)
    valid=ierr==MPI_SUCCESS.and.all(grid>0).and.all(ieee_is_finite(lattice)).and.&
      all(ieee_is_finite(origin)).and.all(ieee_is_finite(lower)).and.all(ieee_is_finite(extent)).and.&
      all(ieee_is_finite(centers)).and.all(shape(lower)==[3,np]).and.all(shape(extent)==[3,np]).and.&
      size(centers,1)==3
    call gate(comm,valid,'invalid center geometry extents or finite values',ok,message);if(.not.ok)return
    length=[(lattice(a,a),a=1,3)];diagonal=0d0
    do a=1,3;diagonal(a,a)=length(a);enddo
    valid=all(length>tiny(1d0)).and.all(length<sqrt(huge(1d0))).and.&
      maxval(abs(lattice-diagonal))==0d0.and.all(abs(origin)<sqrt(huge(1d0))).and.&
      all(abs(lower)<sqrt(huge(1d0))).and.all(abs(extent)<sqrt(huge(1d0))).and.&
      all(abs(centers)<sqrt(huge(1d0)))
    ! Keep integer coordinate products compatible with the existing DC grid mapper.
    valid=valid.and.product(real(grid,real64))<=real(huge(0),real64)
    call gate(comm,valid,'unsupported non-axis-aligned or out-of-range center geometry',ok,message)
    if(.not.ok)return
    allocate(metadata(15+6*np),reference(15+6*np),starts(3,np),widths(3,np),stat=status)
    call gate(comm,status==0,'center geometry allocation failed',ok,message);if(.not.ok)return
    metadata=[reshape(lattice,[9]),origin,real(grid,real64),reshape(lower,[3*np]),reshape(extent,[3*np])]
    reference=metadata
    call MPI_Bcast(reference,size(reference),MPI_DOUBLE_PRECISION,0,comm,ierr)
    call gate(comm,ierr==MPI_SUCCESS.and.all(metadata==reference),&
      'center geometry differs between ranks',ok,message);if(.not.ok)return
    h=length/real(grid,real64)
    ! Reject extreme scale ratios before converting positions to integer grid indices.
    valid=all(h>=1d0/sqrt(huge(1d0)))
    call gate(comm,valid,'unresolved center grid spacing',ok,message);if(.not.ok)return
    tol=roundoff_factor*max(1d0,real(grid,real64),abs(origin)/h)
    valid=all(tol<0.25d0)
    call gate(comm,valid,'origin precision cannot resolve the core grid',ok,message);if(.not.ok)return
    volume=product(int(grid,int64));total_volume=0_int64
    do f=1,np
      lo=(lower(:,f)-origin)/h;wide=extent(:,f)/h
      if(any(lo < -tol).or.any(lo>=real(grid,real64)+tol).or.&
          any(wide<1d0-tol).or.any(wide>real(grid,real64)+tol))then
        valid=.false.;cycle
      endif
      starts(:,f)=modulo(nint(lo),grid);widths(:,f)=nint(wide)
      valid=valid.and.all(abs(lo-anint(lo))<=tol).and.all(abs(wide-anint(wide))<=tol)
      box_volume=product(int(widths(:,f),int64))
      if(box_volume>volume-total_volume)then
        valid=.false.
      else
        total_volume=total_volume+box_volume
      endif
    enddo
    valid=valid.and.total_volume==volume
    call gate(comm,valid,'core boxes are not a grid-aligned complete partition',ok,message);if(.not.ok)return
    do f=1,np;do g=f+1,np
      overlap=.true.
      do a=1,3
        overlap=overlap.and.(modulo(starts(a,g)-starts(a,f),grid(a))<widths(a,f).or.&
          modulo(starts(a,f)-starts(a,g),grid(a))<widths(a,g))
      enddo
      if(overlap)valid=.false.
    enddo;enddo
    call gate(comm,valid,'core boxes overlap',ok,message);if(.not.ok)return
    n=size(centers,2);allocate(work%center_owner(n),stat=status)
    call gate(comm,status==0,'center ownership allocation failed',ok,message);if(.not.ok)return
    work%center_owner=0
    do j=1,n
      ! Wrap physical coordinates first; every box uses the same snapped grid coordinate.
      x=modulo(centers(:,j)-origin,length)/h
      where(abs(x-anint(x))<=tol)x=anint(x)
      x=modulo(x,real(grid,real64));match=0
      do f=1,np
        distance=modulo(x-real(starts(:,f),real64),real(grid,real64))
        if(all(distance<real(widths(:,f),real64)))then
          match=match+1;work%center_owner(j)=f
        endif
      enddo
      valid=valid.and.match==1
    enddo
    call gate(comm,valid,'center has no unique core owner',ok,message);if(.not.ok)return
    work%geometry_fingerprint=hash_reals(metadata,917_int64)
    work%geometry_fingerprint=mix(work%geometry_fingerprint,int(center_convention,int64))
    work%geometry_fingerprint=mix(work%geometry_fingerprint,transfer(roundoff_factor,0_int64))
    work%raw_count=n;work%valid=.true.;result=work;ok=.true.;message=''
#else
    ok=.false.;message='center ownership requires MPI'
#endif
  end subroutine

  subroutine select_dg_hybrid_core_wannier(comm,fragment_id,cache,fragment_lattice,raw_origin,&
      total_lattice,total_origin,core_lower,core_extent,raw_grid,core_grid,total_grid,dc_indices,result,ok,message)
    integer,intent(in)::comm,fragment_id,raw_grid(3),core_grid(3),total_grid(3),dc_indices(:,:)
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::cache
    real(real64),intent(in)::fragment_lattice(3,3),raw_origin(3),total_lattice(3,3),total_origin(3),&
      core_lower(:,:),core_extent(:,:)
    type(s_dg_hybrid_core_selection),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    type(s_dg_hybrid_core_selection)::work
    integer::np,rank,ierr,status,a,j,base(3),slot
    integer,allocatable::fragments(:)
    complex(real64),allocatable::q(:,:),seed(:,:)
    real(real64),allocatable::centers(:,:)
    integer(int64),allocatable::physical_ids(:)
    logical,allocatable::core_mask(:)
    real(real64)::h(3),length(3),expected_lattice(3,3),offset(3),tol(3)
    integer(int64)::fingerprint
    logical::valid,cache_ok
    character(256)::cache_message
    call MPI_Comm_size(comm,np,ierr);call MPI_Comm_rank(comm,rank,ierr)
    valid=ierr==MPI_SUCCESS.and.fragment_id>=1.and.fragment_id<=np.and.cache%valid.and.&
      allocated(cache%local_grid_ids).and.allocated(cache%centers_fractional).and.&
      all(ieee_is_finite(fragment_lattice)).and.all(ieee_is_finite(raw_origin)).and.&
      all(raw_grid>0).and.all(core_grid>0).and.all(core_grid<=raw_grid).and.all(total_grid>0).and.&
      size(dc_indices,1)>=maxval(raw_grid).and.size(dc_indices,2)==3
    call gate(comm,valid,'invalid single-owner WF selection input',ok,message);if(.not.ok)return
    allocate(fragments(np),stat=status)
    call gate(comm,status==0,'fragment ownership allocation failed',ok,message);if(.not.ok)return
    call MPI_Allgather(fragment_id,1,MPI_INTEGER,fragments,1,MPI_INTEGER,comm,ierr)
    valid=ierr==MPI_SUCCESS
    do j=1,np;valid=valid.and.count(fragments==j)==1;enddo
    call gate(comm,valid,'WF selection requires one rank per fragment',ok,message);if(.not.ok)return
    ! Reuse the authoritative cache validator; these temporary maps are not new seed projections.
    call export_dg_hybrid_fragment_coordinates(MPI_COMM_SELF,fragment_id,cache%receipt%basis_generation,&
      cache%receipt%seed_fingerprint,cache%receipt%basis_fingerprint,cache%local_grid_ids,&
      cache%local_row_layout_fingerprint,cache,q,seed,cache_ok,cache_message)
    call gate(comm,cache_ok,'WF selection rejected invalid raw cache integrity',ok,message);if(.not.ok)return
    deallocate(q,seed)
    call map_dg_hybrid_fragment_dc_grid(MPI_COMM_SELF,raw_grid,core_grid,total_grid,dc_indices,&
      cache%local_grid_ids,physical_ids,core_mask,cache_ok,cache_message)
    call gate(comm,cache_ok,'WF selection rejected DC grid mapping',ok,message);if(.not.ok)return
    ! Classify a zero-column probe first to certify global geometry before arithmetic/indexing.
    allocate(centers(3,0),stat=status)
    call gate(comm,status==0,'center workspace allocation failed',ok,message);if(.not.ok)return
    call classify_dg_hybrid_core_centers(comm,total_lattice,total_origin,core_lower,core_extent,&
      total_grid,centers,work,ok,message)
    if(.not.ok)return
    length=[(total_lattice(a,a),a=1,3)];h=length/real(total_grid,real64)
    tol=roundoff_factor*max(1d0,real(total_grid,real64),abs(total_origin)/h)
    expected_lattice=0d0
    do a=1,3;expected_lattice(a,a)=h(a)*real(raw_grid(a),real64);enddo
    valid=all(abs(raw_origin)<sqrt(huge(1d0))).and.all(abs(fragment_lattice)<sqrt(huge(1d0)))
    call gate(comm,valid,'out-of-range fragment center geometry',ok,message);if(.not.ok)return
    valid=maxval(abs(fragment_lattice-expected_lattice))<=roundoff_factor*max(1d0,maxval(abs(expected_lattice)))
    offset=modulo(raw_origin-total_origin,length)/h
    valid=valid.and.all(abs(offset-anint(offset))<=tol)
    base=modulo(nint(offset),total_grid)
    valid=valid.and.all(abs(modulo(core_lower(:,fragment_id)-total_origin,length)/h-&
      real(base,real64))<=tol).and.all(abs(core_extent(:,fragment_id)/h-real(core_grid,real64))<=tol)
    do a=1,3;do j=1,raw_grid(a)
      valid=valid.and.dc_indices(j,a)==1+int(modulo(int(base(a),int64)+int(j-1,int64),int(total_grid(a),int64)))
    enddo;enddo
    call gate(comm,valid,'fragment lattice/origin disagrees with raw DC core-first mapping',ok,message)
    if(.not.ok)return
    deallocate(centers);allocate(centers(3,cache%receipt%retained_rank),stat=status)
    call gate(comm,status==0,'WF center workspace allocation failed',ok,message);if(.not.ok)return
    centers=spread(raw_origin,2,size(centers,2))+matmul(fragment_lattice,cache%centers_fractional)
    call classify_dg_hybrid_core_centers(comm,total_lattice,total_origin,core_lower,core_extent,&
      total_grid,centers,work,ok,message)
    if(.not.ok)return
    work%fragment_id=fragment_id;work%basis_generation=cache%receipt%basis_generation
    work%raw_cache_fingerprint=cache%receipt%replicated_payload_fingerprint
    work%selected_count=count(work%center_owner==fragment_id)
    allocate(work%raw_column_ids(work%selected_count),stat=status)
    call gate(comm,status==0,'selected column allocation failed',ok,message);if(.not.ok)return
    slot=0
    do j=1,work%raw_count
      if(work%center_owner(j)/=fragment_id)cycle
      slot=slot+1;work%raw_column_ids(slot)=int(j,int64)
    enddo
    fingerprint=hash_reals([reshape(fragment_lattice,[9]),raw_origin],work%geometry_fingerprint)
    do a=1,3
      fingerprint=mix(fingerprint,int(raw_grid(a),int64));fingerprint=mix(fingerprint,int(core_grid(a),int64))
      do j=1,raw_grid(a);fingerprint=mix(fingerprint,int(dc_indices(j,a),int64));enddo
    enddo
    work%geometry_fingerprint=fingerprint
    fingerprint=mix(fingerprint,work%raw_cache_fingerprint)
    fingerprint=mix(fingerprint,cache%receipt%distributed_wannier_fingerprint)
    fingerprint=mix(fingerprint,int(fragment_id,int64));fingerprint=mix(fingerprint,int(rank,int64))
    do j=1,np;fingerprint=mix(fingerprint,int(fragments(j),int64));enddo
    do j=1,work%selected_count;fingerprint=mix(fingerprint,work%raw_column_ids(j));enddo
    work%fingerprint=fingerprint;result=work;ok=.true.;message=''
#else
    ok=.false.;message='core-centered WF selection requires MPI'
#endif
  end subroutine
#ifdef USE_MPI
  subroutine gate(comm,valid,description,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::valid
    character(*),intent(in)::description
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::bad,global_bad,ierr
    bad=merge(0,1,valid)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.global_bad==0;message=''
    if(.not.ok)message=description
  end subroutine
  integer(int64) function mix(hash,value)result(next)
    integer(int64),intent(in)::hash,value
    next=ieor(ishftc(hash,7),value)
    if(next==0_int64)next=1_int64
  end function
  integer(int64) function hash_reals(values,initial)result(hash)
    real(real64),intent(in)::values(:)
    integer(int64),intent(in)::initial
    integer::j
    hash=mix(initial,int(size(values),int64))
    do j=1,size(values);hash=mix(hash,transfer(values(j),0_int64));enddo
  end function
#endif
end module
