module dg_overlapping_wannier_projection
  use,intrinsic::iso_fortran_env,only:int64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  implicit none
  private

  type,public::t_dg_projection_channel
    integer::atom=0
    integer::l=-1
    integer::m=0
  end type

  public::build_dg_complete_sp_manifest,dg_complete_sp_target_count,&
    validate_dg_complete_sp_manifest
  public::build_dg_complete_shell_manifest,dg_complete_shell_target_count,&
    validate_dg_complete_shell_manifest
  public::evaluate_dg_periodic_sp_projectors
  public::dg_real_harmonic_value
  public::dg_periodic_grid_point_owned,select_dg_sp_projector_ordinals,&
    select_dg_sp_atomic_orbital_ordinals
contains
  pure logical function dg_periodic_grid_point_owned(point,origin,extent,period) result(owned)
    integer,intent(in)::point(3),origin(3),extent(3),period(3)
    owned=all(period>0).and.all(extent>0).and.all(extent<=period)
    if(owned)owned=all(modulo(point-origin,period)<extent)
  end function

  subroutine select_dg_sp_atomic_orbital_ordinals(species_lmax,nproj,ordinals,ok,message)
    integer,intent(in)::species_lmax(:),nproj(0:,:)
    integer,intent(out)::ordinals(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::species

    ok=.false.;message='';ordinals=-1
    if(size(ordinals,1)/=2.or.size(ordinals,2)/=size(species_lmax).or.&
        size(nproj,2)<size(species_lmax).or.ubound(nproj,1)<1)then
      message='invalid s+p pseudo-atomic orbital selector contract';return
    endif
    do species=1,size(species_lmax)
      if(species_lmax(species)<1.or.any(nproj(0:1,species)<1))then
        message='pseudopotential lacks a complete s+p atomic-orbital shell';return
      endif
      ordinals(1,species)=0
      ordinals(2,species)=nproj(0,species)
    enddo
    ok=.true.
  end subroutine

  subroutine select_dg_sp_projector_ordinals(species_lmax,nproj,inorm,ordinals,ok,message,&
      atomic_orbital_fallback)
    integer,intent(in)::species_lmax(:),nproj(0:,:),inorm(0:,:)
    integer,intent(out)::ordinals(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    logical,intent(out),optional::atomic_orbital_fallback(:,:)
    integer::species,ll,first,last,ordinal
    ok=.false.;message='';ordinals=-1
    if(present(atomic_orbital_fallback))then
      if(any(shape(atomic_orbital_fallback)/=shape(ordinals)))then
        message='invalid atomic-orbital fallback selector shape';return
      endif
      atomic_orbital_fallback=.false.
    endif
    if(size(ordinals,1)/=2.or.size(ordinals,2)/=size(species_lmax).or.&
        size(nproj,2)<size(species_lmax).or.size(inorm,2)<size(species_lmax).or.&
        ubound(nproj,1)<1)then
      message='invalid s+p pseudopotential selector contract';return
    endif
    do species=1,size(species_lmax)
      if(species_lmax(species)<1.or.any(nproj(0:1,species)<1))then
        message='pseudopotential lacks a complete s+p projector shell';return
      endif
      do ll=0,1
        first=0;if(ll>0)first=sum(nproj(0:ll-1,species))
        last=first+nproj(ll,species)-1
        if(last>ubound(inorm,1))then
          message='pseudopotential s+p projector ordinal is out of range';return
        endif
        do ordinal=first,last
          if(inorm(ordinal,species)/=0)then;ordinals(ll+1,species)=ordinal;exit;endif
        enddo
        if(ordinals(ll+1,species)<0)then
          if(.not.present(atomic_orbital_fallback))then
            message='pseudopotential local-reference s+p channel requires atomic-orbital data';return
          endif
          ordinals(ll+1,species)=first
          atomic_orbital_fallback(ll+1,species)=.true.
        endif
      enddo
    enddo
    ok=.true.
  end subroutine

  subroutine dg_complete_shell_target_count(core_atom_count,max_l,target_count,ok,message)
    integer,intent(in)::core_atom_count,max_l
    integer,intent(out)::target_count
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::channels_per_atom

    target_count=0;ok=.false.;message=''
    if(core_atom_count<=0)then
      message='complete shell manifest requires a positive core-owned atom count'
      return
    endif
    select case(max_l)
    case(1);channels_per_atom=4
    case(2);channels_per_atom=9
    case default
      message='complete shell manifest supports only s+p or s+p+d'
      return
    end select
    if(core_atom_count>huge(target_count)/channels_per_atom)then
      message='complete shell manifest target count overflow'
      return
    endif
    target_count=channels_per_atom*core_atom_count
    ok=.true.
  end subroutine

  subroutine dg_complete_sp_target_count(core_atom_count,target_count,ok,message)
    integer,intent(in)::core_atom_count
    integer,intent(out)::target_count
    logical,intent(out)::ok
    character(*),intent(out)::message

    call dg_complete_shell_target_count(core_atom_count,1,target_count,ok,message)
  end subroutine

  subroutine build_dg_complete_shell_manifest(core_atom_ids,max_l,channels,ok,message)
    integer,intent(in)::core_atom_ids(:)
    integer,intent(in)::max_l
    type(t_dg_projection_channel),allocatable,intent(out)::channels(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::target_count,atom_index,channel,m,other

    call dg_complete_shell_target_count(size(core_atom_ids),max_l,target_count,ok,message)
    if(.not.ok)return
    if(any(core_atom_ids<=0))then
      ok=.false.;message='complete shell manifest contains an invalid atom ID';return
    endif
    do atom_index=1,size(core_atom_ids)
      do other=atom_index+1,size(core_atom_ids)
        if(core_atom_ids(atom_index)==core_atom_ids(other))then
          ok=.false.;message='complete shell manifest contains a duplicate atom ID';return
        endif
      enddo
    enddo
    allocate(channels(target_count));channel=0
    do atom_index=1,size(core_atom_ids)
      channel=channel+1
      channels(channel)=t_dg_projection_channel(core_atom_ids(atom_index),0,1)
      do m=1,3
        channel=channel+1
        channels(channel)=t_dg_projection_channel(core_atom_ids(atom_index),1,m)
      enddo
      if(max_l==2)then
        do m=1,5
          channel=channel+1
          channels(channel)=t_dg_projection_channel(core_atom_ids(atom_index),2,m)
        enddo
      endif
    enddo
    ok=.true.;message=''
  end subroutine

  subroutine build_dg_complete_sp_manifest(core_atom_ids,channels,ok,message)
    integer,intent(in)::core_atom_ids(:)
    type(t_dg_projection_channel),allocatable,intent(out)::channels(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    call build_dg_complete_shell_manifest(core_atom_ids,1,channels,ok,message)
  end subroutine

  subroutine validate_dg_complete_shell_manifest(channels,core_atom_ids,max_l,ok,message)
    type(t_dg_projection_channel),intent(in)::channels(:)
    integer,intent(in)::core_atom_ids(:),max_l
    logical,intent(out)::ok
    character(*),intent(out)::message
    type(t_dg_projection_channel),allocatable::expected(:)
    integer::i

    call build_dg_complete_shell_manifest(core_atom_ids,max_l,expected,ok,message)
    if(.not.ok)return
    if(size(channels)/=size(expected))then
      ok=.false.;message='complete shell manifest has a truncated or extra shell';return
    endif
    do i=1,size(expected)
      if(channels(i)%atom/=expected(i)%atom.or.channels(i)%l/=expected(i)%l.or.&
          channels(i)%m/=expected(i)%m)then
        ok=.false.;message='complete shell manifest ordering or member is invalid';return
      endif
    enddo
    ok=.true.;message=''
  end subroutine

  subroutine validate_dg_complete_sp_manifest(channels,core_atom_ids,ok,message)
    type(t_dg_projection_channel),intent(in)::channels(:)
    integer,intent(in)::core_atom_ids(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    call validate_dg_complete_shell_manifest(channels,core_atom_ids,1,ok,message)
  end subroutine

  subroutine evaluate_dg_periodic_sp_projectors(lattice,lattice_inverse,positions,&
      atom_positions,radial_grid,radial_projector,radial_count,pseudopotential_fingerprint,&
      channels,values,ok,message)
    real(8),intent(in)::lattice(3,3),lattice_inverse(3,3),positions(:,:),&
      atom_positions(:,:),radial_grid(:,:,:),radial_projector(:,:,:)
    integer,intent(in)::radial_count(:,:)
    integer(int64),intent(in)::pseudopotential_fingerprint
    type(t_dg_projection_channel),intent(in)::channels(:)
    real(8),allocatable,intent(out)::values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8)::fractional_delta(3),cartesian_delta(3),direction(3),radius,radial_value
    integer,allocatable::core_atom_ids(:)
    integer::channel,point,atom,shell_count
    logical::image_ok

    ok=.false.;message=''
    if(size(positions,1)/=3.or.size(atom_positions,1)/=3.or.size(channels)<1.or.&
        size(radial_grid,2)/=2.or.size(radial_grid,3)/=size(atom_positions,2).or.&
        size(radial_projector,1)/=size(radial_grid,1).or.&
        size(radial_projector,2)<2.or.size(radial_projector,3)/=size(atom_positions,2).or.&
        any(shape(radial_count)/=[2,size(atom_positions,2)]).or.&
        pseudopotential_fingerprint==0_int64.or..not.all(ieee_is_finite(lattice)).or.&
        .not.all(ieee_is_finite(lattice_inverse)).or.&
        .not.all(ieee_is_finite(positions)).or..not.all(ieee_is_finite(atom_positions)).or.&
        .not.all(ieee_is_finite(radial_grid)).or..not.all(ieee_is_finite(radial_projector)))then
      message='invalid periodic s+p projector contract';return
    endif
    if(maxval(abs(matmul(lattice_inverse,lattice)-identity3()))>1d-10)then
      message='periodic s+p projector lattice inverse mismatch';return
    endif
    if(any(channels%atom<1).or.any(channels%atom>size(atom_positions,2)).or.&
        any(channels%l<0).or.any(channels%l>1))then
      message='periodic s+p projector channel is out of range';return
    endif
    shell_count=count(channels%l==0)
    if(shell_count<1)then;message='periodic s+p projector has no complete shell';return;endif
    allocate(core_atom_ids(shell_count));core_atom_ids=pack(channels%atom,channels%l==0)
    call validate_dg_complete_sp_manifest(channels,core_atom_ids,ok,message)
    if(.not.ok)return
    do atom=1,size(atom_positions,2)
      if(any(radial_count(:,atom)<2).or.any(radial_count(:,atom)>size(radial_grid,1)))then
        message='invalid pseudopotential radial count';ok=.false.;return
      endif
      do channel=1,2
        if(any(radial_grid(2:radial_count(channel,atom),channel,atom)<=&
            radial_grid(1:radial_count(channel,atom)-1,channel,atom)))then
          message='pseudopotential radial grid is not strictly increasing';ok=.false.;return
        endif
      enddo
    enddo
    allocate(values(size(channels),size(positions,2)));values=0d0
    do channel=1,size(channels)
      atom=channels(channel)%atom
      do point=1,size(positions,2)
        fractional_delta=matmul(lattice_inverse,positions(:,point)-atom_positions(:,atom))
        call closest_periodic_cartesian(lattice,lattice_inverse,fractional_delta,&
          cartesian_delta,image_ok)
        if(.not.image_ok)then
          message='periodic s+p projector image search failed';return
        endif
        radius=sqrt(sum(cartesian_delta**2))
        radial_value=interpolate_radial(radius,radial_grid(:,channels(channel)%l+1,atom),&
          radial_projector(:,channels(channel)%l+1,atom),&
          radial_count(channels(channel)%l+1,atom))
        if(radius>tiny(1d0))then
          direction=cartesian_delta/radius
        else
          direction=0d0
        endif
        values(channel,point)=radial_value*&
          dg_real_harmonic_value(channels(channel)%l,channels(channel)%m,direction)
      enddo
    enddo
    ok=.true.;message=''
  end subroutine

  real(8) function interpolate_radial(radius,grid,radial,count) result(value)
    real(8),intent(in)::radius,grid(:),radial(:)
    integer,intent(in)::count
    integer::i
    value=0d0
    if(radius<0d0.or.radius>grid(count))return
    if(radius<=grid(1))then;value=radial(1);return;endif
    do i=2,count
      if(radius<=grid(i))then
        value=radial(i-1)+(radial(i)-radial(i-1))*&
          (radius-grid(i-1))/(grid(i)-grid(i-1))
        return
      endif
    enddo
  end function

  pure real(8) function dg_real_harmonic_value(l,m,r) result(value)
    integer,intent(in)::l,m
    real(8),intent(in)::r(3)
    real(8),parameter::sqrt2=sqrt(2d0),inv_sqrt2=1d0/sqrt(2d0),&
      inv_sqrt6=1d0/sqrt(6d0)
    value=0d0
    select case(l)
    case(0)
      if(m==1)value=1d0
    case(1)
      select case(m)
      case(1);value=r(3)
      case(2);value=r(1)
      case(3);value=r(2)
      end select
    case(2)
      select case(m)
      case(1);value=inv_sqrt6*(2d0*r(3)*r(3)-r(1)*r(1)-r(2)*r(2))
      case(2);value=sqrt2*r(3)*r(1)
      case(3);value=sqrt2*r(2)*r(3)
      case(4);value=inv_sqrt2*(r(1)*r(1)-r(2)*r(2))
      case(5);value=sqrt2*r(1)*r(2)
      end select
    end select
  end function

  subroutine closest_periodic_cartesian(lattice,lattice_inverse,fractional_delta,cartesian,ok)
    real(8),intent(in)::lattice(3,3),lattice_inverse(3,3),fractional_delta(3)
    real(8),intent(out)::cartesian(3)
    logical,intent(out)::ok
    real(8)::candidate(3),candidate_cartesian(3),best2,bound
    integer::lower(3),upper(3),n1,n2,n3

    candidate=fractional_delta-dnint(fractional_delta)
    cartesian=matmul(lattice,candidate);best2=sum(cartesian**2)
    bound=sqrt(sum(lattice_inverse**2))*sqrt(best2)+1d-12
    lower=floor(fractional_delta-bound)-1
    upper=ceiling(fractional_delta+bound)+1
    ok=.false.
    if(any(upper-lower>200))return
    do n3=lower(3),upper(3);do n2=lower(2),upper(2);do n1=lower(1),upper(1)
      candidate=fractional_delta-real([n1,n2,n3],8)
      candidate_cartesian=matmul(lattice,candidate)
      if(sum(candidate_cartesian**2)<best2)then
        best2=sum(candidate_cartesian**2);cartesian=candidate_cartesian
      endif
    enddo;enddo;enddo
    ok=all(ieee_is_finite(cartesian))
  end subroutine

  pure function identity3() result(identity)
    real(8)::identity(3,3)
    integer::i
    identity=0d0
    do i=1,3;identity(i,i)=1d0;enddo
  end function
end module
