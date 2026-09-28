! Input and basis fixture; production communication/localization/transport are linked unchanged.
module salmon_global
 implicit none
 integer :: exx_mlwf_maxiter=100,hse_lcfo_u_interval=1
 character(1) :: yn_hse_wannier='y',yn_hse_lcfo_seed_distributed='n'
 real(8) :: exx_mlwf_tolerance=1d-7,hse_lcfo_wf_radius=0d0
end module
module lcfo_rt_basis
 implicit none
 logical :: lcfo_direct_wf=.false.
 complex(8),allocatable :: lcfo_basis(:,:)
 integer,allocatable :: lcfo_counts(:),lcfo_offsets(:),lcfo_origins(:,:)
 integer :: lcfo_grid(3)=[8,2,2],lcfo_core(3)=[8,2,2],lcfo_buffer(3)=0,lcfo_rank=0,lcfo_comm=0
 real(8) :: lcfo_h(3)=1d0,lcfo_dv=1d0
end module
