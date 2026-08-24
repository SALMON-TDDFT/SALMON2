#include "config.h"
module dg_hybrid_production_pw_basis
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
  use dg_hybrid_reciprocal_catalog,only:build_dg_hybrid_reciprocal_catalog
  use dg_hybrid_window_distribution,only:prepare_dg_hybrid_window_distribution
  use dg_hybrid_windowed_pw_basis,only:build_dg_hybrid_windowed_pw_basis
  implicit none
  private
  public::build_dg_hybrid_production_pw_basis
contains
  subroutine build_dg_hybrid_production_pw_basis(comm,global_point_count,fragment_count,&
      fragment_ids,box_ids,box_windows,core_ids,core_fragment_ids,coordinates,row_action,&
      reciprocal_lattice,reciprocal_rotation,cutoff,tile_width,tolerance,windows,g_vectors,&
      catalog,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,fragment_count,fragment_ids(:),core_fragment_ids(:)
    integer(int64),intent(in)::box_ids(:),core_ids(:)
    real(real64),intent(in)::box_windows(:,:),coordinates(:,:),reciprocal_lattice(3,3),&
      reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,intent(in)::row_action(:,:),tile_width
    real(real64),allocatable,intent(out)::windows(:,:),g_vectors(:,:)
    type(s_dg_hybrid_basis_catalog),intent(out)::catalog
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer,allocatable::fragment_action(:,:),g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable::raw_windows(:,:)
    integer(int64)::window_workspace,window_fingerprint,reciprocal_fingerprint,basis_workspace,basis_fingerprint
    logical::stage_ok
    character(256)::stage_message
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    catalog%valid=.false.
    call prepare_dg_hybrid_window_distribution(comm,global_point_count,fragment_count,fragment_ids,&
      box_ids,box_windows,core_ids,core_fragment_ids,row_action,raw_windows,fragment_action,&
      window_workspace,window_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,cutoff,tolerance,&
      g_integer,g_vectors,g_action,g_star,g_conjugate,reciprocal_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);return;endif
    call build_dg_hybrid_windowed_pw_basis(comm,global_point_count,core_ids,coordinates,raw_windows,&
      fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,&
      tile_width,tolerance,windows,catalog,basis_workspace,basis_fingerprint,stage_ok,stage_message)
    if(.not.stage_ok)then;message=trim(stage_message);catalog%valid=.false.;return;endif
    catalog%window_fingerprint=window_fingerprint
    catalog%packet_fingerprint=ieor(reciprocal_fingerprint,basis_fingerprint)
    catalog%catalog_fingerprint=ieor(ieor(ishftc(window_fingerprint,11),&
      ishftc(reciprocal_fingerprint,23)),basis_fingerprint)
    if(catalog%catalog_fingerprint==0_int64)catalog%catalog_fingerprint=1_int64
    fingerprint=catalog%catalog_fingerprint
    workspace_peak_bytes=max(window_workspace,basis_workspace)
    deallocate(raw_windows,fragment_action,g_integer,g_action,g_star,g_conjugate)
    ok=.true.
  end subroutine build_dg_hybrid_production_pw_basis
end module dg_hybrid_production_pw_basis
