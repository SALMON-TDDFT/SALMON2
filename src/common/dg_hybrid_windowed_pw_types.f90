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
end module dg_hybrid_windowed_pw_types
