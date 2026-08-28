module dg_hybrid_windowed_pw_types
  implicit none
  private

  type,public::s_dg_hybrid_pw_packet
    integer::fragment_id=0
    integer::star_id=0
    integer::owner_rank=-1
    integer,allocatable::g_indices(:)
  end type s_dg_hybrid_pw_packet

  type,public::s_dg_hybrid_basis_catalog
    logical::valid=.false.
    integer(8)::wannier_fingerprint=0_8
    integer(8)::window_fingerprint=0_8
    integer(8)::packet_fingerprint=0_8
    integer(8)::catalog_fingerprint=0_8
    integer,allocatable::accepted_wannier_blocks(:)
    integer,allocatable::rejected_wannier_blocks(:)
    type(s_dg_hybrid_pw_packet),allocatable::packets(:)
  end type s_dg_hybrid_basis_catalog

  type,public::s_dg_hybrid_production_selection
    logical::analysis_complete=.false.
    logical::identity_only=.false.
    integer::operation_count=0
    integer(8)::analysis_fingerprint=0_8
    integer(8)::window_fingerprint=0_8
    integer(8)::packet_fingerprint=0_8
    integer,allocatable::fragment_action(:,:)
    integer,allocatable::row_action(:,:)
    integer,allocatable::reciprocal_action(:,:)
    integer,allocatable::packet_ids(:)
    integer,allocatable::packet_action(:,:)
    real(8),allocatable::reciprocal_rotation(:,:,:)
    type(s_dg_hybrid_pw_packet),allocatable::packets(:)
  end type s_dg_hybrid_production_selection
end module dg_hybrid_windowed_pw_types
