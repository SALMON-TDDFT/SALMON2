! Total-grid rVV10 potential with native spatial halo derivatives.
module rvv10_native
  use structures, only: s_rgrid,s_parallel_info,s_sendrecv_grid,s_stencil,s_dft_system,s_scalar
  use sendrecv_grid, only: update_overlap_real8
  use stencil_sub, only: calc_gradient_field
  use rvv10_distributed, only: rvv10_evaluate_distributed
  implicit none
  private
  public :: rvv10_native_periodic
contains
  subroutine rvv10_native_periodic(n,mg,info,srg,system,stencil,rho,b,c,nq,energy,potential,used,status)
    integer,intent(in) :: n(3),nq
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_dft_system),intent(in) :: system
    type(s_stencil),intent(in) :: stencil
    type(s_scalar),intent(in) :: rho
    real(8),intent(in) :: b,c
    real(8),intent(out) :: energy(:,:,:),potential(:,:,:)
    logical,intent(out) :: used
    integer,intent(out) :: status
    real(8),allocatable :: halo(:,:,:),grad(:,:,:,:),fluxgrad(:,:,:,:),r(:),s(:),e(:),v(:),w(:)
    integer :: lo(3),hi(3),ng,d
    lo=mg%is;hi=mg%ie;ng=product(mg%num)
    allocate(halo(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3)))
    allocate(grad(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
      fluxgrad(3,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    allocate(r(ng),s(ng),e(ng),v(ng),w(ng))
    halo=0d0;halo(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=rho%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    if(info%if_divide_rspace)call update_overlap_real8(srg,mg,halo)
    call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,halo,grad)
    r=reshape(rho%f(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),[ng]);s=reshape(sum(grad**2,dim=1),[ng])
    call rvv10_evaluate_distributed(n,lo,mg%num, &
      [info%isize_x,info%isize_y,info%isize_z],[info%id_x,info%id_y,info%id_z], &
      [info%icomm_x,info%icomm_y,info%icomm_z],info%icomm_r,system%hgs,r,s,b,c,nq,e,v,w,used,status)
    if(.not.used.or.status/=0)return
    energy=reshape(e,mg%num);potential=reshape(v,mg%num)
    do d=1,3
      halo=0d0
      halo(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=2*reshape(w,mg%num)*grad(d,:,:,:)
      if(info%if_divide_rspace)call update_overlap_real8(srg,mg,halo)
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,halo,fluxgrad)
      potential=potential-fluxgrad(d,:,:,:)
    enddo
  end subroutine
end module
