module dg_hybrid_sparse_operators
  use,intrinsic::iso_fortran_env,only:int64,real64
  implicit none
  private
  type,public::s_dg_hybrid_sparse_operators
    logical::valid=.false.
    integer::global_count=0
    integer(int64)::selection_fingerprint=0_int64,window_fingerprint=0_int64,&
      packet_fingerprint=0_int64,complement_fingerprint=0_int64,metric_fingerprint=0_int64,&
      fingerprint=0_int64,persistent_bytes=0_int64,transient_peak_bytes=0_int64
    integer(int64),allocatable::owned_row_ids(:)
    integer,allocatable::row_offsets(:),column_ids(:)
    complex(real64),allocatable::metric_values(:),hamiltonian_values(:),position_values(:,:)
  end type s_dg_hybrid_sparse_operators
end module dg_hybrid_sparse_operators
