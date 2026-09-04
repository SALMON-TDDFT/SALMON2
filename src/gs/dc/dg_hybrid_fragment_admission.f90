module dg_hybrid_fragment_admission
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,&
    ieee_set_flag,ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_fragment_wannier,only:s_dg_hybrid_fragment_wannier_cache
  use dg_hybrid_fragment_selection,only:s_dg_hybrid_core_selection,s_dg_hybrid_selected_catalog,&
    s_dg_hybrid_dc_reference,prepare_dg_hybrid_selected_catalog,export_dg_hybrid_dc_reference,&
    export_dg_hybrid_selected_frame
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  use dg_hybrid_projected_fragment_pipeline,only:s_dg_hybrid_projection_factorization_receipt,&
    s_dg_hybrid_core_projection_report,s_dg_hybrid_support_samples,validate_dg_hybrid_projected_basis,&
    project_dg_hybrid_core_seeds,check_dg_hybrid_seed_support,dg_hybrid_core_quadrature_binding
  use dg_hybrid_fragment_subspace,only:s_dg_hybrid_fragment_subspace_state,&
    initialize_dg_hybrid_fragment_density_checked
  implicit none
  private
  type,public::s_dg_hybrid_support_operator
    private
    logical::valid=.false.
    integer::fragment_id=0,generation=0,channel=0
    integer(int64)::fingerprint=0_int64
    integer(int64),allocatable::sample_ids(:),point_ids(:)
    integer,allocatable::offsets(:)
    complex(real64),allocatable::coefficients(:)
    real(real64),allocatable::weights(:)
  end type
  type,public::s_dg_hybrid_admission_report
    logical::valid=.false.,support_measured=.false.
    type(s_dg_hybrid_core_projection_report)::core
    real(real64)::support_defects(3)=huge(1d0),density_defects(2)=huge(1d0)
    real(real64)::final_support_defects(3)=huge(1d0)
    integer(int64)::basis_fingerprint=0_int64,metric_fingerprint=0_int64,selection_fingerprint=0_int64
  end type
  public::prepare_dg_hybrid_support_operator,admit_dg_hybrid_selected_fragment
  public::export_dg_hybrid_selected_basis_frame
contains
  ! Bind the projected reference to the actual selected WF+PW payload. This
  ! exports coordinates only; C3 metric/support/state admission is still
  ! required before production CG, and H remains the operator producer's job.
  subroutine export_dg_hybrid_selected_basis_frame(comm,fragment_id,cache,selection,basis,receipt,&
      frame,fingerprint,ok,message)
    integer,intent(in)::comm,fragment_id
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::cache
    type(s_dg_hybrid_core_selection),intent(in)::selection
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    type(s_dg_hybrid_projection_factorization_receipt),intent(in)::receipt
    complex(real64),allocatable,intent(out)::frame(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_selected_catalog)::catalog
    complex(real64),allocatable::wf_frame(:,:),work(:,:)
    integer::n,nw,nraw,npw,nref,status,j
    integer(int64)::wf_fp,hash,nref64
    logical::valid
    fingerprint=0_int64
    call prepare_dg_hybrid_selected_catalog(comm,fragment_id,cache,selection,catalog,ok,message)
    if(.not.ok)return
    call validate_dg_hybrid_projected_basis(comm,basis,receipt,catalog%fingerprint,ok,message)
    if(.not.ok)return
    n=size(basis%global_ids);nw=selection%selected_count;nraw=selection%raw_count
    valid=basis%fragment_id==fragment_id.and.basis%generation==selection%basis_generation.and.n>=nw.and.&
      size(basis%buffer_point_ids)==size(selection%physical_grid_ids)
    call gate(comm,valid,'selected reference basis extent/generation mismatch',ok,message);if(.not.ok)return
    valid=all(basis%global_ids(:nw)==catalog%local_active_ids).and.all(basis%sector(:nw)==1).and.&
      all(basis%sector(nw+1:)==2).and.all(basis%buffer_point_ids==selection%physical_grid_ids).and.&
      all(basis%buffer_values(:,:nw)==transpose(catalog%local_values))
    call gate(comm,valid,'selected reference basis payload mismatch',ok,message);if(.not.ok)return
    call export_dg_hybrid_selected_frame(comm,fragment_id,cache,selection,wf_frame,wf_fp,ok,message)
    if(.not.ok)return
    npw=n-nw;nref64=int(nraw,int64)+int(npw,int64)
    call gate(comm,nref64<=int(huge(0),int64),'selected reference extent overflow',ok,message);if(.not.ok)return
    nref=int(nref64)
    allocate(work(n,nref),stat=status)
    call gate(comm,status==0,'selected WF+PW reference allocation failed',ok,message);if(.not.ok)return
    work=0d0;work(:nw,:nraw)=wf_frame
    do j=1,npw;work(nw+j,nraw+j)=1d0;enddo
    hash=mix(1913_int64,wf_fp);hash=mix(hash,receipt%payload_fingerprint)
    hash=mix(hash,int(n,int64));hash=mix(hash,int(nref,int64));hash=mix(hash,int(npw,int64))
    if(hash==0_int64)hash=1_int64
    call move_alloc(work,frame);fingerprint=hash
  end subroutine

  ! Freeze sparse linear functionals on physical grid IDs. required_ids must
  ! come from the operator's independent required inventory, not from filtering
  ! available samples. C5 connects the actual face/stencil/projector producers.
  subroutine prepare_dg_hybrid_support_operator(comm,fragment_id,generation,channel,required_ids,&
      sample_ids,offsets,point_ids,coefficients,weights,operator,fingerprint,ok,message)
    integer,intent(in)::comm,fragment_id,generation,channel,offsets(:)
    integer(int64),intent(in)::required_ids(:),sample_ids(:),point_ids(:)
    complex(real64),intent(in)::coefficients(:)
    real(real64),intent(in)::weights(:)
    type(s_dg_hybrid_support_operator),intent(out)::operator
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(s_dg_hybrid_support_operator)::work
    integer::np,ierr,n,nz,j,k,status
    integer(int64)::hash
    logical::valid
    fingerprint=0_int64
    call MPI_Comm_size(comm,np,ierr);n=size(sample_ids);nz=size(point_ids)
    valid=ierr==MPI_SUCCESS.and.fragment_id>=1.and.fragment_id<=np.and.generation>0.and.&
      channel>=1.and.channel<=3.and.n<huge(0).and.nz<huge(0).and.size(required_ids)==n.and.&
      size(weights)==n.and.size(coefficients)==nz.and.size(offsets,kind=int64)==int(n,int64)+1_int64.and.&
      all(sample_ids>0_int64).and.all(point_ids>0_int64).and.all(ieee_is_finite(weights)).and.&
      all(weights>0d0).and.all(ieee_is_finite(real(coefficients))).and.all(ieee_is_finite(aimag(coefficients)))
    call gate(comm,valid,'invalid support operator manifest',ok,message);if(.not.ok)return
    valid=all(sample_ids==required_ids).and.offsets(1)==1.and.offsets(n+1)==nz+1.and.&
      all(offsets>=1).and.all(offsets<=nz+1)
    do j=1,n
      valid=valid.and.offsets(j+1)>offsets(j)
      do k=1,j-1;valid=valid.and.sample_ids(k)/=sample_ids(j);enddo
    enddo
    call gate(comm,valid,'support operator IDs/order or sparse rows disagree with required inventory',ok,message)
    if(.not.ok)return
    allocate(work%sample_ids(n),work%point_ids(nz),work%offsets(n+1),work%coefficients(nz),work%weights(n),stat=status)
    call gate(comm,status==0,'support operator allocation failed',ok,message);if(.not.ok)return
    work%sample_ids=sample_ids;work%point_ids=point_ids;work%offsets=offsets
    work%coefficients=coefficients;work%weights=weights
    hash=mix(1901_int64,int(fragment_id,int64));hash=mix(hash,int(generation,int64))
    hash=mix(hash,int(channel,int64));hash=mix(hash,int(n,int64));hash=mix(hash,int(nz,int64))
    do j=1,n;hash=mix(hash,sample_ids(j));hash=mix(hash,transfer(weights(j),0_int64));enddo
    do j=1,n+1;hash=mix(hash,int(offsets(j),int64));enddo
    do j=1,nz
      hash=mix(hash,point_ids(j));hash=mix(hash,transfer(real(coefficients(j),real64),0_int64))
      hash=mix(hash,transfer(aimag(coefficients(j)),0_int64))
    enddo
    work%fragment_id=fragment_id;work%generation=generation;work%channel=channel
    work%fingerprint=hash;work%valid=.true.;operator=work;fingerprint=hash
  end subroutine

  ! Atomic, single-owner admission. References and support values are rebuilt
  ! here, so callers cannot substitute precomputed densities or sample arrays.
  ! The accepted operator fingerprints identify the current frozen inventory.
  subroutine admit_dg_hybrid_selected_fragment(comm,fragment_id,cache,selection,basis,receipt,operators,&
      expected_operator_fingerprints,weights,core_limits,support_limits,pw_cutoff,guard_count,&
      occupation_tolerance,energy_tolerance,orthogonality_tolerance,state,selected_seeds,report,ok,message,&
      energy_cutoff)
    integer,intent(in)::comm,fragment_id,guard_count
    type(s_dg_hybrid_fragment_wannier_cache),intent(in)::cache
    type(s_dg_hybrid_core_selection),intent(in)::selection
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    type(s_dg_hybrid_projection_factorization_receipt),intent(in)::receipt
    type(s_dg_hybrid_support_operator),intent(in)::operators(3)
    integer(int64),intent(in)::expected_operator_fingerprints(3)
    real(real64),intent(in)::weights(:),core_limits(4),support_limits(3),pw_cutoff,&
      occupation_tolerance,energy_tolerance,orthogonality_tolerance
    real(real64),optional,intent(in)::energy_cutoff
    type(s_dg_hybrid_fragment_subspace_state),intent(inout)::state
    integer,allocatable,intent(out)::selected_seeds(:)
    type(s_dg_hybrid_admission_report),intent(out)::report
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::metric(:,:)
    logical::halting(3)
    call ieee_get_halting_mode(ieee_invalid,halting(1))
    call ieee_get_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_get_halting_mode(ieee_overflow,halting(3))
    call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    call execute()
    call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
    call ieee_set_flag(ieee_overflow,.false.)
    call ieee_set_halting_mode(ieee_invalid,halting(1));call ieee_set_halting_mode(ieee_divide_by_zero,halting(2))
    call ieee_set_halting_mode(ieee_overflow,halting(3))
  contains
    subroutine execute()
      type(s_dg_hybrid_selected_catalog)::catalog
      type(s_dg_hybrid_dc_reference)::reference
      type(s_dg_hybrid_support_samples)::samples(3)
      type(s_dg_hybrid_fragment_subspace_state)::candidate
      integer,allocatable::chosen(:)
      complex(real64),allocatable::core_basis(:,:),coefficients(:,:)
      real(real64),allocatable::density(:)
      integer(int64),allocatable::sorted_ids(:)
      integer,allocatable::point_slots(:)
      integer::nw,n,npoint,ncore,nseed,status,c,j,k,row,slot,counts(3)
      integer(int64)::metric_fp
      logical::valid,final_measured
      character(512)::support_message
      ok=.false.;message=''
      call prepare_dg_hybrid_selected_catalog(comm,fragment_id,cache,selection,catalog,ok,message)
      if(.not.ok)return
      call export_dg_hybrid_dc_reference(comm,fragment_id,cache,selection,reference,ok,message)
      if(.not.ok)return
      call validate_dg_hybrid_projected_basis(comm,basis,receipt,catalog%fingerprint,ok,message)
      if(.not.ok)return
      nw=selection%selected_count;n=size(basis%global_ids);npoint=size(reference%physical_grid_ids)
      ncore=size(reference%core_row_slots);nseed=size(reference%occupations)
      valid=basis%fragment_id==fragment_id.and.basis%generation==selection%basis_generation.and.&
        size(basis%buffer_point_ids)==npoint.and.n>=nw.and.size(weights)==ncore.and.&
        all(ieee_is_finite(weights)).and.all(weights>0d0)
      call gate(comm,valid,'selected admission basis/grid extent mismatch',ok,message);if(.not.ok)return
      valid=receipt%core_quadrature_fingerprint==dg_hybrid_core_quadrature_binding(&
        reference%physical_grid_ids(reference%core_row_slots),weights)
      call gate(comm,valid,'selected admission core quadrature mismatch',ok,message);if(.not.ok)return
      valid=all(basis%buffer_point_ids==reference%physical_grid_ids).and.&
        all(basis%global_ids(:nw)==catalog%local_active_ids).and.all(basis%sector(:nw)==1).and.&
        all(basis%sector(nw+1:)==2).and.&
        all(basis%buffer_values(:,:nw)==transpose(catalog%local_values))
      call gate(comm,valid,'selected admission basis/raw-reference binding mismatch',ok,message);if(.not.ok)return
      do c=1,3
        valid=valid.and.operators(c)%valid.and.operators(c)%fragment_id==fragment_id.and.&
          operators(c)%generation==selection%basis_generation.and.operators(c)%channel==c.and.&
          expected_operator_fingerprints(c)/=0_int64.and.&
          operators(c)%fingerprint==expected_operator_fingerprints(c)
      enddo
      call gate(comm,valid,'selected admission support operator context mismatch',ok,message);if(.not.ok)return
      allocate(core_basis(ncore,n),metric(n,n),density(ncore),sorted_ids(npoint),point_slots(npoint),stat=status)
      call gate(comm,status==0,'selected admission allocation failed',ok,message);if(.not.ok)return
      sorted_ids=reference%physical_grid_ids;point_slots=[(j,j=1,npoint)]
      call sort_points(sorted_ids,point_slots,1,npoint)
      do j=2,npoint;valid=valid.and.sorted_ids(j)>sorted_ids(j-1);enddo
      call gate(comm,valid,'ambiguous physical support in selected admission',ok,message);if(.not.ok)return
      core_basis=basis%buffer_values(reference%core_row_slots,:)
      call project_dg_hybrid_core_seeds(comm,core_basis,weights,reference%core_orbitals,reference%occupations,&
        core_limits,nw,pw_cutoff,coefficients,report%core,ok,message)
      if(.not.ok)return
      do c=1,3
        counts(c)=size(operators(c)%sample_ids)
        allocate(samples(c)%basis(counts(c),n),samples(c)%reference(counts(c),nseed),&
          samples(c)%weights(counts(c)),stat=status)
        call gate(comm,status==0,'selected support evaluation allocation failed',ok,message);if(.not.ok)return
        samples(c)%basis=0d0;samples(c)%reference=0d0;samples(c)%weights=operators(c)%weights
        do row=1,counts(c)
          do k=operators(c)%offsets(row),operators(c)%offsets(row+1)-1
            j=lookup_point(sorted_ids,operators(c)%point_ids(k))
            if(j==0)then;valid=.false.;cycle;endif
            slot=point_slots(j)
            samples(c)%basis(row,:)=samples(c)%basis(row,:)+&
              operators(c)%coefficients(k)*basis%buffer_values(slot,:)
            samples(c)%reference(row,:)=samples(c)%reference(row,:)+&
              operators(c)%coefficients(k)*reference%buffer_orbitals(slot,:)
          enddo
        enddo
      enddo
      call gate(comm,valid,'required operator point missing from raw/selected buffer',ok,message);if(.not.ok)return
      call check_dg_hybrid_seed_support(comm,coefficients,samples,counts,support_limits,nw,pw_cutoff,&
        report%support_defects,report%support_measured,ok,message)
      if(.not.ok)return
      metric=matmul(conjg(transpose(core_basis)),core_basis*spread(weights,2,n))
      metric_fp=mix(receipt%payload_fingerprint,selection%fingerprint)
      do j=1,ncore
        metric_fp=mix(metric_fp,int(reference%core_row_slots(j),int64))
        metric_fp=mix(metric_fp,transfer(weights(j),0_int64))
      enddo
      density=0d0
      do j=1,nseed;density=density+reference%occupations(j)*abs(reference%core_orbitals(:,j))**2;enddo
      report%basis_fingerprint=receipt%payload_fingerprint;report%metric_fingerprint=metric_fp
      report%selection_fingerprint=selection%fingerprint
      call initialize_dg_hybrid_fragment_density_checked(comm,fragment_id,selection%basis_generation,&
        receipt%payload_fingerprint,metric_fp,coefficients,reference%energies,reference%occupations,guard_count,&
        occupation_tolerance,energy_tolerance,orthogonality_tolerance,core_basis,weights,density,&
        core_limits(3),core_limits(4),apply_metric,candidate,chosen,report%density_defects,ok,message,energy_cutoff)
      if(.not.ok)return
      do c=1,3
        samples(c)%reference=samples(c)%reference(:,chosen)
      enddo
      call check_dg_hybrid_seed_support(comm,candidate%vectors,samples,counts,support_limits,nw,pw_cutoff,&
        report%final_support_defects,final_measured,ok,support_message)
      if(.not.ok)then
        message='post-initialization support: '//trim(support_message)
        return
      endif
      state=candidate
      call move_alloc(chosen,selected_seeds)
      report%valid=ok
    end subroutine
    subroutine apply_metric(input,output,valid)
      complex(real64),intent(in)::input(:,:)
      complex(real64),intent(out)::output(:,:)
      logical,intent(out)::valid
      output=matmul(metric,input)
      valid=all(ieee_is_finite(real(output))).and.all(ieee_is_finite(aimag(output)))
    end subroutine
  end subroutine

  recursive subroutine sort_points(ids,slots,left,right)
    integer(int64),intent(inout)::ids(:)
    integer,intent(inout)::slots(:)
    integer,intent(in)::left,right
    integer::i,j,tmp
    integer(int64)::pivot,word
    if(left>=right)return
    i=left;j=right;pivot=ids(left+(right-left)/2)
    do
      do while(ids(i)<pivot);i=i+1;enddo
      do while(ids(j)>pivot);j=j-1;enddo
      if(i>j)exit
      word=ids(i);ids(i)=ids(j);ids(j)=word;tmp=slots(i);slots(i)=slots(j);slots(j)=tmp
      i=i+1;j=j-1
      if(i>j)exit
    enddo
    if(left<j)call sort_points(ids,slots,left,j)
    if(i<right)call sort_points(ids,slots,i,right)
  end subroutine
  integer function lookup_point(ids,target)result(slot)
    integer(int64),intent(in)::ids(:),target
    integer::lo,hi,mid
    lo=1;hi=size(ids);slot=0
    do while(lo<=hi)
      mid=lo+(hi-lo)/2
      if(ids(mid)==target)then;slot=mid;return
      else if(ids(mid)<target)then;lo=mid+1
      else;hi=mid-1
      endif
    enddo
  end function
  subroutine gate(comm,valid,description,ok,message)
    integer,intent(in)::comm
    logical,intent(in)::valid
    character(*),intent(in)::description
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::bad,total,ierr
    bad=merge(0,1,valid)
    call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.total==0;message=''
    if(.not.ok)message=description
  end subroutine
  integer(int64) function mix(hash,word)result(next)
    integer(int64),intent(in)::hash,word
    next=ieor(ishftc(hash,11),word)
    next=ieor(ishftc(next,7),int(z'9E3779B97F4A7C15',int64))
    if(next==0_int64)next=1_int64
  end function
end module dg_hybrid_fragment_admission
