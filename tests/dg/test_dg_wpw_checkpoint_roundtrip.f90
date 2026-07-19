program test_dg_wpw_checkpoint_roundtrip
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  use dg_wpw_checkpoint,only:s_dg_wpw_checkpoint_state,write_dg_wpw_checkpoint,read_dg_wpw_checkpoint,&
    write_dg_wpw_checkpoint_manifest,read_dg_wpw_checkpoint_manifest,dg_wpw_checkpoint_checksum
  implicit none
  type(s_dg_wpw_checkpoint_state)::state,loaded,bad
  integer::info,unit
  integer(8),allocatable::checksums(:)
  character(256)::path

  path='wpw.chk'
  call make_state(state)
  call write_dg_wpw_checkpoint(path,state,info)
  if(info/=0)error stop 1
  call read_dg_wpw_checkpoint(path,loaded,info,expected_fingerprint=state%layout_fingerprint)
  if(info/=0.or..not.same_state(state,loaded))error stop 2

  call read_dg_wpw_checkpoint(path,bad,info,expected_fingerprint=state%layout_fingerprint+1_8)
  if(info==0)error stop 3

  bad=state;bad%eigenvalues(1)=ieee_value(0d0,ieee_quiet_nan)
  call write_dg_wpw_checkpoint('nonfinite.chk',bad,info)
  if(info==0)error stop 4

  open(newunit=unit,file='truncated.chk',access='stream',form='unformatted',status='replace')
  write(unit)'DG_WPW_CHECKPT'
  close(unit)
  call read_dg_wpw_checkpoint('truncated.chk',bad,info)
  if(info==0)error stop 5

  open(newunit=unit,file=path,access='stream',form='unformatted',status='old',position='append')
  write(unit)123456789
  close(unit)
  call read_dg_wpw_checkpoint(path,bad,info)
  if(info==0)error stop 6

  call write_dg_wpw_checkpoint('missing/parent/wpw.chk',state,info)
  if(info==0)error stop 7
  call write_dg_wpw_checkpoint_manifest('wpw.manifest',[dg_wpw_checkpoint_checksum(state),123_8],info)
  if(info/=0)error stop 8
  call read_dg_wpw_checkpoint_manifest('wpw.manifest',checksums,info)
  if(info/=0.or.any(checksums/=[dg_wpw_checkpoint_checksum(state),123_8]))error stop 9
  open(newunit=unit,file='bad.manifest',access='stream',form='unformatted',status='replace')
  write(unit)'DG_WPW_MANIFEST1',1,2,[1_8]
  close(unit)
  call read_dg_wpw_checkpoint_manifest('bad.manifest',checksums,info)
  if(info==0)error stop 10
  print '(a)','PASS versioned bounded WPW checkpoint roundtrip and fail-closed validation'
contains
  subroutine make_state(s)
    type(s_dg_wpw_checkpoint_state),intent(out)::s
    integer::i
    s%schema_version=2;s%operator_epoch=7;s%layout_fingerprint=987654321_8
    s%ownership_kind='fragment_root_v1';s%metric_convention='orthonormal_ww'
    s%operator_convention='windowed_kg_sipg_v1';s%n_occ=1
    s%peer_ranks=[integer::]
    s%owned_w_ids=[1,2];s%owned_p_ids=[5,6];s%support_w_ids=[1,2,3];s%support_p_ids=[5,6,7]
    s%eigenvalues=[-0.5d0,0.25d0];s%occupations=[2d0,0d0]
    allocate(s%coeff_w(2,2),s%coeff_p(2,2));s%coeff_w=reshape([(cmplx(i,-i,8),i=1,4)],[2,2])
    s%coeff_p=0.5d0*s%coeff_w;s%potential=[0.1d0,0.2d0,0.3d0]
    s%ww_r=[1,2];s%ww_c=[1,2];s%wp_w=[1,2];s%wp_p=[5,6];s%pp_r=[5,6];s%pp_c=[5,6]
    s%ww_h=[(1d0,0d0),(2d0,0d0)];s%ww_s=[(1d0,0d0),(1d0,0d0)]
    s%wp_h=[(0.1d0,0.2d0),(0.3d0,0.4d0)];s%wp_s=[(0.01d0,0d0),(0.02d0,0d0)]
    s%pp_h=[(3d0,0d0),(4d0,0d0)];s%pp_s=[(1d0,0d0),(1d0,0d0)]
    s%ww_z=[(0.1d0,0d0),(0.2d0,0d0)];s%wp_z=[(0.3d0,0d0),(0.4d0,0d0)]
    s%pp_z=[(0.5d0,0d0),(0.6d0,0d0)]
  end subroutine
  logical function same_state(a,b)
    type(s_dg_wpw_checkpoint_state),intent(in)::a,b
    same_state=a%schema_version==b%schema_version.and.a%operator_epoch==b%operator_epoch.and.&
      a%layout_fingerprint==b%layout_fingerprint.and.a%n_occ==b%n_occ.and.&
      trim(a%ownership_kind)==trim(b%ownership_kind).and.all(a%owned_w_ids==b%owned_w_ids).and.&
      all(a%peer_ranks==b%peer_ranks).and.&
      all(a%owned_p_ids==b%owned_p_ids).and.all(a%support_w_ids==b%support_w_ids).and.&
      all(a%support_p_ids==b%support_p_ids).and.maxval(abs(a%coeff_w-b%coeff_w))==0d0.and.&
      maxval(abs(a%coeff_p-b%coeff_p))==0d0.and.all(a%eigenvalues==b%eigenvalues).and.&
      all(a%occupations==b%occupations).and.all(a%potential==b%potential).and.&
      all(a%ww_r==b%ww_r).and.all(a%ww_c==b%ww_c).and.all(a%wp_w==b%wp_w).and.&
      all(a%wp_p==b%wp_p).and.all(a%pp_r==b%pp_r).and.all(a%pp_c==b%pp_c).and.&
      all(a%ww_h==b%ww_h).and.all(a%ww_s==b%ww_s).and.all(a%wp_h==b%wp_h).and.&
      all(a%wp_s==b%wp_s).and.all(a%pp_h==b%pp_h).and.all(a%pp_s==b%pp_s).and.&
      all(a%ww_z==b%ww_z).and.all(a%wp_z==b%wp_z).and.all(a%pp_z==b%pp_z)
  end function
end program
