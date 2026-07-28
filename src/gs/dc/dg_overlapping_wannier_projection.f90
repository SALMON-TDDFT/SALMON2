module dg_overlapping_wannier_projection
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
contains
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
end module
