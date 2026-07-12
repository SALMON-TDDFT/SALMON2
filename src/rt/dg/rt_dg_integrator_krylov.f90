  subroutine time_evolution_krylov(dg_frag,system,mg,ppg,rt,itt,dt)
    use structures
    use communication, only: comm_summation,comm_is_root
    use eigen_subdiag_sub, only: eigen_zheev
    implicit none
    type(s_dg_fragment_rt),intent(inout) :: dg_frag
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_pp_grid),intent(in) :: ppg
    type(s_rt),intent(inout) :: rt
    integer,intent(in) :: itt
    real(8),intent(in) :: dt

    integer,parameter :: krylov_dim_max=64
    real(8),parameter :: breakdown_tol=1.0d-11
    integer :: nlocal,nstate_prop,istate,j,i,it0,it1,kdim_current
    integer :: used_dim_max,used_dim_sum
    real(8) :: Ac_mid(3),phase_c,phase_s,final_beta_max
    complex(8),allocatable :: initial(:,:),q(:,:,:),w(:,:),deriv(:,:,:)
    complex(8),allocatable :: hessenberg(:,:,:),hsmall(:,:),hvec(:,:),coeff(:),tmp(:)
    complex(8),allocatable :: dots(:)
    real(8),allocatable :: heval(:),norms(:),beta(:)
    integer,allocatable :: kdim(:)
    logical,allocatable :: active(:)

    if (dg_frag%nspin /= 1) stop 'DG Krylov prototype currently requires nspin=1'
    if (dg_frag%use_plane_wave_basis .or. allocated(dg_frag%coef_pw)) &
      stop 'DG Krylov prototype supports the pure fragment basis only'
    nlocal=size(dg_frag%coef,1)
    nstate_prop=dg_frag%nstate_tot
    if (allocated(dg_frag%nocc_spin)) nstate_prop=min(nstate_prop,max(1,dg_frag%nocc_spin(1)))
    it0=max(lbound(rt%Ac_tot,2),itt-1)
    it1=min(ubound(rt%Ac_tot,2),itt)
    Ac_mid=0.5d0*(rt%Ac_tot(:,it0)+rt%Ac_tot(:,it1))

    allocate(initial(nlocal,nstate_prop),q(nlocal,krylov_dim_max+1,nstate_prop))
    allocate(w(nlocal,nstate_prop),deriv(nlocal,nstate_prop,1))
    allocate(hessenberg(krylov_dim_max+1,krylov_dim_max,nstate_prop))
    allocate(hsmall(krylov_dim_max,krylov_dim_max),hvec(krylov_dim_max,krylov_dim_max))
    allocate(coeff(krylov_dim_max),tmp(krylov_dim_max),dots(nstate_prop))
    allocate(heval(krylov_dim_max),norms(nstate_prop),beta(nstate_prop),kdim(nstate_prop))
    allocate(active(nstate_prop))

    initial=dg_frag%coef(:,1:nstate_prop,1)
    call global_dots(initial,initial,dots)
    do istate=1,nstate_prop
      norms(istate)=sqrt(max(0.0d0,real(dots(istate),8)))
      if (norms(istate) <= 1.0d-14) stop 'DG Krylov encountered a zero-norm occupied state'
    end do
    q=(0.0d0,0.0d0)
    hessenberg=(0.0d0,0.0d0)
    do istate=1,nstate_prop
      q(:,1,istate)=initial(:,istate)/norms(istate)
    end do
    active=.true.
    kdim=krylov_dim_max
    beta=0.0d0

    do j=1,krylov_dim_max
      do istate=1,nstate_prop
        if (active(istate)) then
          dg_frag%coef(:,istate,1)=q(:,j,istate)
        else
          dg_frag%coef(:,istate,1)=(0.0d0,0.0d0)
        end if
      end do
      deriv=(0.0d0,0.0d0)
      call calculate_time_derivative(dg_frag,system,mg,ppg,Ac_mid,deriv,1,nstate_prop)
      w=deriv(:,:,1)

      do i=1,j
        call global_dots(q(:,i,1:nstate_prop),w,dots)
        do istate=1,nstate_prop
          if (.not. active(istate)) cycle
          hessenberg(i,j,istate)=dots(istate)
          w(:,istate)=w(:,istate)-dots(istate)*q(:,i,istate)
        end do
      end do
      do i=1,j
        call global_dots(q(:,i,1:nstate_prop),w,dots)
        do istate=1,nstate_prop
          if (.not. active(istate)) cycle
          hessenberg(i,j,istate)=hessenberg(i,j,istate)+dots(istate)
          w(:,istate)=w(:,istate)-dots(istate)*q(:,i,istate)
        end do
      end do
      call global_dots(w,w,dots)
      do istate=1,nstate_prop
        if (.not. active(istate)) cycle
        beta(istate)=sqrt(max(0.0d0,real(dots(istate),8)))
        if (j < krylov_dim_max) hessenberg(j+1,j,istate)=cmplx(beta(istate),0.0d0,kind=8)
        if (beta(istate) <= breakdown_tol) then
          kdim(istate)=j
          active(istate)=.false.
        else if (j < krylov_dim_max) then
          q(:,j+1,istate)=w(:,istate)/beta(istate)
        end if
      end do
      if (.not. any(active)) exit
    end do

    do istate=1,nstate_prop
      kdim_current=kdim(istate)
      hsmall(1:kdim_current,1:kdim_current)=cmplx(0.0d0,1.0d0,kind=8)* &
        hessenberg(1:kdim_current,1:kdim_current,istate)
      do j=1,kdim_current
        hsmall(j,j)=cmplx(real(hsmall(j,j),8),0.0d0,kind=8)
        do i=j+1,kdim_current
          hsmall(i,j)=0.5d0*(hsmall(i,j)+conjg(hsmall(j,i)))
          hsmall(j,i)=conjg(hsmall(i,j))
        end do
      end do
      call eigen_zheev(hsmall(1:kdim_current,1:kdim_current),heval(1:kdim_current), &
                       hvec(1:kdim_current,1:kdim_current))
      tmp(1:kdim_current)=conjg(hvec(1,1:kdim_current))
      do i=1,kdim_current
        phase_c=cos(heval(i)*dt)
        phase_s=sin(heval(i)*dt)
        tmp(i)=cmplx(phase_c,-phase_s,kind=8)*tmp(i)
      end do
      coeff(1:kdim_current)=matmul(hvec(1:kdim_current,1:kdim_current),tmp(1:kdim_current))
      dg_frag%coef(:,istate,1)=norms(istate)* &
        matmul(q(:,1:kdim_current,istate),coeff(1:kdim_current))
    end do
    used_dim_max=maxval(kdim)
    used_dim_sum=sum(kdim)
    final_beta_max=maxval(beta)
    if (comm_is_root(dg_frag%id) .and. itt == 1) then
      write(*,'(1x,a,i0,a,1pe13.5,a,i0,a,0pf8.3,a,1pe13.5)') '[DG-KRYLOV] max_dim=', &
        krylov_dim_max,' breakdown_tol=',breakdown_tol,' used_dim_max=',used_dim_max, &
        ' used_dim_avg=',dble(used_dim_sum)/dble(max(1,nstate_prop)), &
        ' final_beta_max=',final_beta_max
      flush(6)
    end if
    deallocate(initial,q,w,deriv,hessenberg,hsmall,hvec,coeff,tmp,dots)
    deallocate(heval,norms,beta,kdim,active)

  contains
    subroutine global_dots(x,y,values)
      complex(8),intent(in) :: x(:,:),y(:,:)
      complex(8),intent(out) :: values(:)
      real(8) :: local_parts(2*size(values)),global_parts(2*size(values))
      complex(8) :: local_value
      integer :: ilocal,istate_local,gid
      local_parts=0.0d0
      do istate_local=1,min(size(values),size(x,2),size(y,2))
        local_value=(0.0d0,0.0d0)
        do ilocal=1,min(size(x,1),size(y,1))
          gid=ilocal
          if (allocated(dg_frag%local_coef_global_ids)) then
            if (ilocal > size(dg_frag%local_coef_global_ids,1)) cycle
            gid=dg_frag%local_coef_global_ids(ilocal,1)
          end if
          if (allocated(dg_frag%coef_owner)) then
            if (gid < 1 .or. gid > size(dg_frag%coef_owner,1)) cycle
            if (dg_frag%coef_owner(gid,1) /= dg_frag%id) cycle
          end if
          local_value=local_value+conjg(x(ilocal,istate_local))*y(ilocal,istate_local)
        end do
        local_parts(2*istate_local-1)=real(local_value,8)
        local_parts(2*istate_local)=aimag(local_value)
      end do
      call comm_summation(local_parts,global_parts,size(local_parts),dg_frag%icomm)
      do istate_local=1,size(values)
        values(istate_local)=cmplx(global_parts(2*istate_local-1),global_parts(2*istate_local),kind=8)
      end do
    end subroutine global_dots
  end subroutine time_evolution_krylov
