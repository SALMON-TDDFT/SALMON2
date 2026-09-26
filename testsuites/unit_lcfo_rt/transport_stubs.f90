! Single-rank communication fixture; production localization/transport is linked unchanged.
module salmon_global
 integer :: hse_mlwf_maxiter=100
 real(8) :: hse_mlwf_tolerance=1d-7
end module
module lcfo_rt_basis
 implicit none
 complex(8),allocatable :: lcfo_basis(:,:)
 integer,allocatable :: lcfo_counts(:),lcfo_offsets(:),lcfo_origins(:,:)
 integer :: lcfo_grid(3)=[8,2,2],lcfo_core(3)=[8,2,2],lcfo_buffer(3)=0,lcfo_rank=0,lcfo_comm=0
 real(8) :: lcfo_h(3)=1d0,lcfo_dv=1d0
end module
module communication
 implicit none
 interface comm_bcast
  module procedure b_i,b_r,b_l,b_c3
 end interface
 interface comm_summation
  module procedure s_r1,s_c1,s_c2,s_c4,s_i
 end interface
contains
 subroutine b_i(a,comm,root)
 integer,intent(inout)::a
 integer::comm,root
 end subroutine
 subroutine b_r(a,comm,root)
 real(8),intent(inout)::a
 integer::comm,root
 end subroutine
 subroutine b_l(a,comm,root)
 logical,intent(inout)::a
 integer::comm,root
 end subroutine
 subroutine b_c3(a,comm,root)
 complex(8),intent(inout)::a(:,:,:)
 integer::comm,root
 end subroutine
 subroutine s_r1(a,b,n,comm)
 real(8),intent(in)::a(:)
 real(8),intent(out)::b(:)
 integer::n,comm
 b=a
 end subroutine
 subroutine s_c1(a,b,n,comm)
 complex(8),intent(in)::a(:)
 complex(8),intent(out)::b(:)
 integer::n,comm
 b=a
 end subroutine
 subroutine s_c2(a,b,n,comm)
 complex(8),intent(in)::a(:,:)
 complex(8),intent(out)::b(:,:)
 integer::n,comm
 b=a
 end subroutine
 subroutine s_c4(a,b,n,comm)
 complex(8),intent(in)::a(:,:,:,:)
 complex(8),intent(out)::b(:,:,:,:)
 integer::n,comm
 b=a
 end subroutine
 subroutine s_i(a,b,comm)
 integer,intent(in)::a
 integer,intent(out)::b
 integer::comm
 b=a
 end subroutine
end module
