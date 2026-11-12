program test_grid_velocity
  use structures, only: s_dft_system,s_rgrid,s_stencil,s_pp_grid
  use hsex_grid_velocity, only: grid_velocity
  use stencil_sub, only: calc_gradient_psi
  implicit none
  type(s_dft_system) :: sys
  type(s_rgrid) :: mg
  type(s_stencil) :: stencil
  type(s_pp_grid) :: ppg
  complex(8) :: grid(4,4,4,2),gradient(3,4,4,4,2),projectors(8,2),uv(8,2)
  complex(8),allocatable :: vector(:),velocity(:,:),plus(:),minus(:)
  real(8) :: kv(3),delta,maxerror,error
  integer :: nc,i,j,c,axis,ix,iy,iz
  mg%is=1;mg%ie=4;mg%is_array=1;mg%ie_array=4
  allocate(mg%idx(-3:8),mg%idy(-3:8),mg%idz(-3:8),sys%vec_k(3,1))
  do i=-3,8
    mg%idx(i)=1+modulo(i-1,4);mg%idy(i)=mg%idx(i);mg%idz(i)=mg%idx(i)
  enddo
  sys%vec_k(:,1)=[.17d0,-.21d0,.11d0];sys%vec_Ac=[.03d0,.01d0,-.02d0]
  sys%rmatrix_B=reshape([1d0,.2d0,0d0,.1d0,1.1d0,.2d0,0d0,.1d0,.9d0],[3,3])
  stencil%coef_nab=0d0;stencil%coef_nab(1,:)=.5d0
  ppg%nlma=1
  allocate(ppg%ia_tbl(1),ppg%ia_tbl_so(1),ppg%rinv_uvu(1),ppg%rinv_uvu_so(1),ppg%mps(1))
  allocate(ppg%jxyz(3,8,1),ppg%rxyz(3,8,1),ppg%zekr_uv(8,1,1),ppg%zekr_uv_so(8,1,1,2,1))
  ppg%ia_tbl=1;ppg%ia_tbl_so=1;ppg%rinv_uvu=.7d0;ppg%rinv_uvu_so=.7d0;ppg%mps=8
  do j=1,8
    ppg%jxyz(:,j,1)=1+[mod(j,4),mod(j/2,4),mod(j/3,4)]
    ppg%rxyz(:,j,1)=[.31d0*j,-.17d0*j,.13d0*j]
    do c=1,2
      uv(j,c)=cmplx(cos(dble(j+c)),sin(dble(2*j+c)),8)/4d0
    enddo
  enddo
  maxerror=0d0;delta=1d-5
  do nc=1,2
    do c=1,nc
      do iz=1,4
        do iy=1,4
          do ix=1,4
            grid(ix,iy,iz,c)=cmplx(sin(dble(ix+iy+2*iz+c)),cos(dble(ix+2*iy+iz+c)),8)/8d0
          enddo
        enddo
      enddo
      call calc_gradient_psi(grid(:,:,:,c),gradient(:,:,:,:,c),mg%is_array,mg%ie_array,mg%is,mg%ie, &
        mg%idx,mg%idy,mg%idz,stencil%coef_nab,sys%rmatrix_B)
    enddo
    allocate(vector(64*nc),velocity(64*nc,3),plus(64*nc),minus(64*nc))
    vector=reshape(grid(:,:,:,:nc),[64*nc])
    call phase(sys%vec_k(:,1)+sys%vec_Ac)
    call grid_velocity(sys,mg,stencil,ppg,1,vector,velocity)
    do axis=1,3
      kv=sys%vec_k(:,1)+sys%vec_Ac;kv(axis)=kv(axis)+delta
      call model(kv,plus)
      kv(axis)=kv(axis)-2*delta
      call model(kv,minus)
      error=maxval(abs((plus-minus)/(2*delta)-velocity(:,axis)))
      maxerror=max(maxerror,error)
    enddo
    deallocate(vector,velocity,plus,minus)
  enddo
  if(maxerror>1d-8)error stop 'grid velocity differs from finite-difference Hamiltonian'
  print *, 'GRID_VELOCITY_SCALAR_SPINOR_PASS',maxerror
contains
  subroutine phase(k)
    implicit none
    real(8),intent(in) :: k(3)
    integer :: j
    do j=1,8
      projectors(j,:)=uv(j,:)*exp(cmplx(0d0,-dot_product(k,ppg%rxyz(:,j,1)),8))
    enddo
    ppg%zekr_uv(:,1,1)=projectors(:,1)
    ppg%zekr_uv_so(:,1,1,:,1)=projectors
  end subroutine
  subroutine model(k,result)
    implicit none
    real(8),intent(in) :: k(3)
    complex(8),intent(out) :: result(:)
    complex(8) :: out(4,4,4,nc),overlap
    integer :: c,j,axis,ix,iy,iz
    call phase(k)
    out=.5d0*sum(k*k)*grid(:,:,:,:nc)
    do c=1,nc
      do axis=1,3
        out(:,:,:,c)=out(:,:,:,c)-cmplx(0d0,k(axis),8)*gradient(axis,:,:,:,c)
      enddo
    enddo
    overlap=0d0
    do c=1,nc
      do j=1,8
        ix=ppg%jxyz(1,j,1);iy=ppg%jxyz(2,j,1);iz=ppg%jxyz(3,j,1)
        overlap=overlap+conjg(projectors(j,c))*grid(ix,iy,iz,c)*.7d0
      enddo
    enddo
    do c=1,nc
      do j=1,8
        ix=ppg%jxyz(1,j,1);iy=ppg%jxyz(2,j,1);iz=ppg%jxyz(3,j,1)
        out(ix,iy,iz,c)=out(ix,iy,iz,c)+projectors(j,c)*overlap
      enddo
    enddo
    result=reshape(out,[64*nc])
  end subroutine
end program
