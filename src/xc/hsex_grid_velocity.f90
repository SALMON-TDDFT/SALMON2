! Analytic k derivative of the real-space kinetic and nonlocal pseudopotential action.
module hsex_grid_velocity
  implicit none
  private
  public :: grid_velocity
contains
  subroutine grid_velocity(sys,mg,stencil,ppg,ik,vector,velocity)
    use structures, only: s_dft_system,s_rgrid,s_stencil,s_pp_grid
    use stencil_sub, only: calc_gradient_psi
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_stencil),intent(in) :: stencil
    type(s_pp_grid),intent(in) :: ppg
    integer,intent(in) :: ik
    complex(8),intent(in) :: vector(:)
    complex(8),intent(out) :: velocity(:,:)
    complex(8),allocatable :: grid(:,:,:,:),gradient(:,:,:,:),action(:,:,:,:,:)
    complex(8) :: overlap,position_overlap(3),beta
    integer :: n(3),nc,component,axis,channel,atom,j,ix,iy,iz,nchannel
    real(8) :: inverse_norm
    complex(8),parameter :: zi=(0d0,1d0)
    n=mg%ie-mg%is+1
    nc=size(vector)/product(n)
    allocate(grid(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3),nc))
    allocate(gradient(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3)),action(n(1),n(2),n(3),nc,3))
    grid=0d0
    grid(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:)=reshape(vector,[n,nc])
    do component=1,nc
      call calc_gradient_psi(grid(:,:,:,component),gradient,mg%is_array,mg%ie_array,mg%is,mg%ie, &
        mg%idx,mg%idy,mg%idz,stencil%coef_nab,sys%rmatrix_B)
      do axis=1,3
        action(:,:,:,component,axis)=-zi*gradient(axis,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
          mg%is(3):mg%ie(3))+(sys%vec_k(axis,ik)+sys%vec_Ac(axis))* &
          grid(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),component)
      enddo
    enddo
    nchannel=ppg%nlma
    if(nc==2)nchannel=size(ppg%ia_tbl_so)
    do channel=1,nchannel
      if(nc==2)then
        atom=ppg%ia_tbl_so(channel);inverse_norm=ppg%rinv_uvu_so(channel)
      else
        atom=ppg%ia_tbl(channel);inverse_norm=ppg%rinv_uvu(channel)
      endif
      overlap=0d0;position_overlap=0d0
      do component=1,nc
        do j=1,ppg%mps(atom)
          ix=ppg%jxyz(1,j,atom);iy=ppg%jxyz(2,j,atom);iz=ppg%jxyz(3,j,atom)
          if(nc==2)then
            beta=ppg%zekr_uv_so(j,channel,ik,component,1)
          else
            beta=ppg%zekr_uv(j,channel,ik)
          endif
          overlap=overlap+conjg(beta)*grid(ix,iy,iz,component)
          position_overlap=position_overlap+ppg%rxyz(:,j,atom)*conjg(beta)*grid(ix,iy,iz,component)
        enddo
      enddo
      do component=1,nc
        do j=1,ppg%mps(atom)
          ix=ppg%jxyz(1,j,atom);iy=ppg%jxyz(2,j,atom);iz=ppg%jxyz(3,j,atom)
          if(nc==2)then
            beta=ppg%zekr_uv_so(j,channel,ik,component,1)
          else
            beta=ppg%zekr_uv(j,channel,ik)
          endif
          do axis=1,3
            action(ix-mg%is(1)+1,iy-mg%is(2)+1,iz-mg%is(3)+1,component,axis)= &
              action(ix-mg%is(1)+1,iy-mg%is(2)+1,iz-mg%is(3)+1,component,axis)+ &
              zi*beta*inverse_norm*(position_overlap(axis)-ppg%rxyz(axis,j,atom)*overlap)
          enddo
        enddo
      enddo
    enddo
    velocity=reshape(action,[size(vector),3])
  end subroutine grid_velocity
end module hsex_grid_velocity
