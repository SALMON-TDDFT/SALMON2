module lcfo_wannier_sawf_seed
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: write_sawf_local_eig_amn_mmn
  public :: build_sawf_local_seed_matrices
  public :: select_sawf_local_complete_shells

contains

  subroutine select_sawf_local_complete_shells(channel_atom,expected_per_atom,inside_atom, &
      selected_channel,ok,message)
    integer,intent(in)::channel_atom(:),expected_per_atom(:)
    logical,intent(in)::inside_atom(:)
    integer,allocatable,intent(out)::selected_channel(:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::channel_count(:)
    integer::channel,nselected

    ok=.false.;message=''
    if(size(channel_atom)<=0.or.size(expected_per_atom)<=0.or. &
        size(inside_atom)/=size(expected_per_atom).or.any(expected_per_atom<=0).or. &
        any(channel_atom<1).or.any(channel_atom>size(expected_per_atom)))then
      message='SAWF local projection-shell dimensions are invalid';return
    end if
    do channel=2,size(channel_atom)
      if(channel_atom(channel)<channel_atom(channel-1))then
        message='SAWF projection channels are not atom-major';return
      end if
    end do
    allocate(channel_count(size(expected_per_atom)));channel_count=0
    do channel=1,size(channel_atom)
      channel_count(channel_atom(channel))=channel_count(channel_atom(channel))+1
    end do
    if(any(channel_count/=expected_per_atom))then
      message='SAWF projection channels do not contain complete atomic shells';return
    end if
    nselected=count(inside_atom(channel_atom));allocate(selected_channel(nselected));nselected=0
    do channel=1,size(channel_atom)
      if(.not.inside_atom(channel_atom(channel)))cycle
      nselected=nselected+1;selected_channel(nselected)=channel
    end do
    if(nselected<=0)then
      message='SAWF local core and buffer contain no complete projection shell';return
    end if
    ok=.true.
  end subroutine select_sawf_local_complete_shells

  subroutine build_sawf_local_seed_matrices(states,projections,phase_factor,weight,amn,mmn,ok,message)
    complex(8),intent(in)::states(:,:),projections(:,:),phase_factor(:,:)
    real(8),intent(in)::weight
    complex(8),intent(out)::amn(:,:),mmn(:,:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(8),allocatable::normalized_projection(:,:),phased_states(:,:)
    real(8)::norm2
    integer::projection,neighbor

    ok=.false.;message='';amn=(0d0,0d0);mmn=(0d0,0d0)
    if(size(states,1)<=0.or.size(states,2)<=0.or.size(projections,1)/=size(states,1).or. &
        size(phase_factor,1)/=size(states,1).or.size(phase_factor,2)<=0.or. &
        size(amn,1)/=size(states,2).or.size(amn,2)/=size(projections,2).or. &
        size(mmn,1)/=size(states,2).or.size(mmn,2)/=size(states,2).or. &
        size(mmn,3)/=size(phase_factor,2).or. &
        .not.ieee_is_finite(weight).or.weight<=0d0)then
      message='SAWF local seed matrix dimensions or integration weight are invalid';return
    end if
    if(.not.all(ieee_is_finite(real(states))).or..not.all(ieee_is_finite(aimag(states))).or. &
        .not.all(ieee_is_finite(real(projections))).or..not.all(ieee_is_finite(aimag(projections))).or. &
        .not.all(ieee_is_finite(real(phase_factor))).or..not.all(ieee_is_finite(aimag(phase_factor))))then
      message='SAWF local seed matrix input contains a non-finite value';return
    end if
    allocate(normalized_projection(size(projections,1),size(projections,2)), &
      phased_states(size(states,1),size(states,2)))
    do projection=1,size(projections,2)
      norm2=weight*sum(abs(projections(:,projection))**2)
      if(norm2<=1d-300)then
        message='SAWF local projection has zero norm';return
      end if
      normalized_projection(:,projection)=projections(:,projection)/sqrt(norm2)
    end do
    amn=weight*matmul(conjg(transpose(states)),normalized_projection)
    do neighbor=1,size(phase_factor,2)
      phased_states=spread(phase_factor(:,neighbor),2,size(states,2))*states
      mmn(:,:,neighbor)=weight*matmul(conjg(transpose(states)),phased_states)
    end do
    ok=.true.
  end subroutine build_sawf_local_seed_matrices

  subroutine write_sawf_local_eig_amn_mmn(directory,seedname,energy_ev,amn,mmn,neighbor_gvec,ok,message)
    character(*),intent(in)::directory,seedname
    real(8),intent(in)::energy_ev(:)
    complex(8),intent(in)::amn(:,:),mmn(:,:,:)
    integer,intent(in)::neighbor_gvec(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::unit,ios,iband,iproj,ineighbor,jband
    character(1024)::filename

    ok=.false.;message=''
    if(len_trim(directory)==0.or.len_trim(seedname)==0.or.size(energy_ev)<=0.or. &
        size(amn,1)/=size(energy_ev).or.size(amn,2)<=0.or. &
        size(mmn,1)/=size(energy_ev).or.size(mmn,2)/=size(energy_ev).or. &
        size(neighbor_gvec,1)/=3.or.size(neighbor_gvec,2)/=size(mmn,3).or.size(mmn,3)<=0)then
      message='SAWF local seed dimensions are inconsistent';return
    end if
    if(.not.all(ieee_is_finite(energy_ev)).or. &
        .not.all(ieee_is_finite(real(amn))).or..not.all(ieee_is_finite(aimag(amn))).or. &
        .not.all(ieee_is_finite(real(mmn))).or..not.all(ieee_is_finite(aimag(mmn))))then
      message='SAWF local seed contains a non-finite value';return
    end if

    filename=trim(directory)//'/'//trim(seedname)//'.eig'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .eig open failed';return;end if
    do iband=1,size(energy_ev);write(unit,'(2i8,1x,es23.15)',iostat=ios)iband,1,energy_ev(iband);if(ios/=0)exit;end do
    close(unit);if(ios/=0)then;message='SAWF local .eig write failed';return;end if

    filename=trim(directory)//'/'//trim(seedname)//'.amn'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .amn open failed';return;end if
    write(unit,'(a)',iostat=ios)'SALMON local SAWF projections'
    if(ios==0)write(unit,'(3i10)',iostat=ios)size(amn,1),1,size(amn,2)
    do iproj=1,size(amn,2);do iband=1,size(amn,1)
      if(ios==0)write(unit,'(3i8,2(1x,es23.15))',iostat=ios)iband,iproj,1,amn(iband,iproj)
    end do;end do
    close(unit);if(ios/=0)then;message='SAWF local .amn write failed';return;end if

    filename=trim(directory)//'/'//trim(seedname)//'.mmn'
    open(newunit=unit,file=trim(filename),status='replace',action='write',iostat=ios)
    if(ios/=0)then;message='SAWF local .mmn open failed';return;end if
    write(unit,'(a)',iostat=ios)'SALMON local SAWF overlaps'
    if(ios==0)write(unit,'(3i10)',iostat=ios)size(mmn,1),1,size(mmn,3)
    do ineighbor=1,size(mmn,3)
      if(ios==0)write(unit,'(5i8)',iostat=ios)1,1,neighbor_gvec(:,ineighbor)
      do jband=1,size(mmn,2);do iband=1,size(mmn,1)
        if(ios==0)write(unit,'(2(1x,es23.15))',iostat=ios)mmn(iband,jband,ineighbor)
      end do;end do
    end do
    close(unit);if(ios/=0)then;message='SAWF local .mmn write failed';return;end if
    ok=.true.
  end subroutine write_sawf_local_eig_amn_mmn

end module lcfo_wannier_sawf_seed
