! Static limit helpers. An occupied projector need not consist of eigenvectors,
! but it must be invariant under H before the insulating optical solver is used.
module hsex_frozen_limit
  implicit none
  private
  public :: frozen_optical_tensor,frozen_eta_limit,frozen_static_part
  abstract interface
    subroutine h_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine velocity_action(x,y,status)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine frozen_static_part(raw,mode,tol,kept,defect,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: raw(:,:)
    character(*),intent(in) :: mode
    real(8),intent(in) :: tol
    complex(8),intent(out) :: kept(:,:)
    real(8),intent(out) :: defect
    integer,intent(out) :: status
    real(8) :: norm
    status=1;kept=0d0;defect=huge(1d0)
    if(size(raw,1)<1.or.size(raw,1)/=size(raw,2).or.any(shape(kept)/=shape(raw)))return
    if(.not.screen_finite(raw).or..not.screen_finite([tol]))return
    if(tol<=0d0.or.tol>=1d0)return
    if(mode/='strict'.and.mode/='hermitian')return
    norm=max(1d-12,sqrt(sum(abs(raw)**2)))
    defect=sqrt(sum(abs(raw-conjg(transpose(raw)))**2))/norm
    if(.not.screen_finite([norm,defect]))return
    if(mode=='strict'.and.defect>tol)then
      status=2;return
    endif
    ! Explicit static approximation: preserve Re[x^H chi x] for every x.
    ! This is not the full causal response or a proof of nonequilibrium W.
    kept=.5d0*(raw+conjg(transpose(raw)))
    if(.not.screen_finite(kept))return
    status=0
  end subroutine
  subroutine frozen_optical_tensor(apply,velocity,occupied,vertices,tol,chi,rate,status)
    use hsex_screen_math, only: screen_finite
    use hsex_frozen_response, only: frozen_background_rate
    use hsex_longwave, only: longwave_band
    implicit none
    procedure(h_action) :: apply
    procedure(velocity_action) :: velocity
    complex(8),intent(in) :: occupied(:,:),vertices(:,:)
    real(8),intent(in) :: tol
    complex(8),intent(out) :: chi(:,:)
    real(8),intent(out) :: rate
    integer,intent(out) :: status
    complex(8),allocatable :: hp(:,:),matrix(:,:),psi(:,:),v(:,:),body(:,:),part(:,:),work(:)
    real(8),allocatable :: energy(:),rwork(:)
    integer :: n,nv,nb,a,code
    status=1;chi=0d0;rate=0d0;n=size(occupied,1);nv=size(occupied,2);nb=size(vertices,2)
    if(nv<1.or.n<nv.or.size(vertices,1)/=n.or.any(shape(chi)/=[nb+3,nb+3]))return
    if(.not.screen_finite(occupied).or..not.screen_finite(vertices).or..not.screen_finite([tol]))return
    if(tol<=0d0)return
    allocate(hp(n,nv),matrix(nv,nv),psi(n,nv),v(n,3),body(n,nb),part(nb+3,nb+3), &
      energy(nv),work(max(1,4*nv)),rwork(max(1,3*nv-2)),stat=code)
    if(code/=0)return
    status=2
    do a=1,nv
      call apply(occupied(:,a),hp(:,a),code)
      if(code/=0)return
      if(.not.screen_finite(hp(:,a)))return
    enddo
    call frozen_background_rate(occupied,hp,rate,code)
    if(code/=0)return
    status=3
    if(rate>tol)return
    matrix=matmul(conjg(transpose(occupied)),hp)
    if(maxval(abs(matrix-conjg(transpose(matrix))))>tol)return
    matrix=.5d0*(matrix+conjg(transpose(matrix)))
    call zheev('V','U',nv,matrix,nv,energy,work,size(work),rwork,code)
    if(code/=0)return
    psi=matmul(occupied,matrix)
    status=4
    do a=1,nv
      call velocity(psi(:,a),v,code)
      if(code/=0)return
      if(.not.screen_finite(v))return
      body=vertices*spread(psi(:,a),2,nb)
      call longwave_band(apply,psi,energy(a),v,body,part,code)
      if(code/=0)return
      chi=chi+part
    enddo
    if(.not.screen_finite(chi))return
    status=0
  end subroutine

  subroutine frozen_eta_limit(samples,tol,limit,error,status)
    use hsex_screen_math, only: screen_finite
    implicit none
    complex(8),intent(in) :: samples(:,:,:)
    real(8),intent(in) :: tol
    complex(8),intent(out) :: limit(:,:)
    real(8),intent(out) :: error
    integer,intent(out) :: status
    real(8) :: norm
    status=1;limit=0d0;error=huge(1d0)
    if(size(samples,3)/=3.or.any(shape(limit)/=shape(samples(:,:,1))))return
    if(.not.screen_finite([tol]).or.tol<=0d0.or.tol>=1d0)return
    if(.not.screen_finite(samples))return
    ! Samples are eta, eta/2, eta/4 at the SAME nonzero q. Include linear
    ! eta terms as well as eta**2; neither analyticity in eta**2 nor a gap
    ! is assumed here. This is an acceptance estimate, not a proof of a limit.
    limit=(samples(:,:,1)-6d0*samples(:,:,2)+8d0*samples(:,:,3))/3d0
    if(.not.screen_finite(limit))return
    norm=max(1d-12,sqrt(sum(abs(limit)**2)))
    if(.not.screen_finite([norm]))return
    error=max(sqrt(sum(abs(limit-samples(:,:,3))**2)), &
      sqrt(sum(abs(2d0*samples(:,:,3)-3d0*samples(:,:,2)+samples(:,:,1))**2)))/norm
    if(.not.screen_finite([error]))return
    status=2
    if(error>tol)return
    status=0
  end subroutine
end module
