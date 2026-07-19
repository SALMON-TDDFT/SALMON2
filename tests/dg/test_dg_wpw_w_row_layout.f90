program test_dg_wpw_w_row_layout
  use dg_wpw_w_row_layout,only:build_dg_wpw_w_row_layout
  implicit none
  integer::n_basis(8),index_basis(2,8),info,k
  integer,allocatable::owned(:),support(:)
  n_basis=2
  do k=1,8;index_basis(:,k)=[2*k-1,2*k];enddo
  call build_dg_wpw_w_row_layout(1,[1,2,3,4,5,6,7,8],n_basis,index_basis,owned,support,info)
  if(info/=0.or.any(owned/=[1,2]))error stop 1
  if(size(support)/=16.or.any(support/=[(k,k=1,16)]))error stop 2
  call build_dg_wpw_w_row_layout(1,[1,2,3,4],n_basis,index_basis,owned,support,info)
  if(info/=0.or.size(support)/=8.or.any(support/=[(k,k=1,8)]))error stop 3
  call build_dg_wpw_w_row_layout(1,[1,3,3],n_basis,index_basis,owned,support,info)
  if(info==0.or.allocated(owned).or.allocated(support))error stop 4
  print '(a)','PASS face-edge-corner WPW W-row layout'
end program
