! Discrete variational noncollinear BR89-CS generalized Kohn-Sham bridge.
! Density, tau and paramagnetic current are Pauli traces (0,x,y,z).
! tau follows the paper's convention WITHOUT a factor 1/2.
module nc_mgga23
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_stencil,s_sendrecv_grid,s_pp_nlcc,s_scalar
  use nc_mgga23_kernel, only: nc_mgga23_eval
  use nc_mgga23_core, only: core23_fields,core23_energy_deda
  use stencil_sub, only: calc_gradient_psi,calc_gradient_field,calc_laplacian_field
  use sendrecv_grid, only: update_overlap_real8,update_overlap_complex8
  use communication, only: comm_summation
  use math_constants, only: zi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: nc_mgga23_refresh,nc_mgga23_add_action
contains
  function traces(a,b) result(v)
    implicit none
    complex(8),intent(in) :: a(2),b(2)
    complex(8) :: v(4)
    v(1)=conjg(a(1))*b(1)+conjg(a(2))*b(2)
    v(2)=conjg(a(1))*b(2)+conjg(a(2))*b(1)
    v(3)=-zi*conjg(a(1))*b(2)+zi*conjg(a(2))*b(1)
    v(4)=conjg(a(1))*b(1)-conjg(a(2))*b(2)
  end function
  function act(v,p) result(q)
    implicit none
    real(8),intent(in) :: v(:)
    complex(8),intent(in) :: p(2)
    complex(8) :: q(2)
    q(1)=(v(1)+v(4))*p(1)+(v(2)-zi*v(3))*p(2)
    q(2)=(v(2)+zi*v(3))*p(1)+(v(1)-v(4))*p(2)
  end function

  subroutine nc_mgga23_refresh(system,mg,info,spsi,stencil,srg_scalar,srg,ppn,E_xc,eexc)
    implicit none
    type(s_dft_system),intent(inout) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(inout) :: spsi
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg_scalar,srg
    type(s_pp_nlcc),intent(in) :: ppn
    real(8),intent(out) :: E_xc
    type(s_scalar),optional,intent(inout) :: eexc
    real(8),allocatable :: fields(:,:,:,:),total(:,:,:,:),valence(:,:,:,:),grad(:,:,:,:,:),lap(:,:,:,:)
    real(8),allocatable :: work(:,:,:),gwork(:,:,:,:),lwork(:,:,:),dg(:,:,:,:,:),dl(:,:,:,:)
    complex(8),allocatable :: gp(:,:,:,:,:)
    complex(8) :: p(2),dp(2)
    real(8) :: n(4),g(4,3),l(4),t(4),j(4,3),dn(4),dgrad(4,3),dlap(4),dt(4),dj(4,3)
    real(8) :: occ,k(3),e,esum,eglob,band_local,band_global
    real(8) :: cn(4),ct(4),cj(4,3),deda(3),deda_local(3),deda_global(3)
    integer :: x,y,z,a,c,sp,io,ik,ierr,bad,global_bad
    if(allocated(spsi%rwf).or.system%nspin/=2) error stop 'nc_mgga23 requires complex two-component spinors'
    if(info%im_s/=1.or.info%im_e/=1) error stop 'nc_mgga23 requires one replica'
    ! Only V2 metadata certifies joint auxiliary-orbital rho/tau/current provenance.
    ! Historical tau_nlcc may be a von Weizsaecker estimate and is never used here.
    bad=0;system%nc_mgga23_payload%core_deda=0d0
    if(ppn%core23_complete) then
      if(.not.allocated(ppn%rho_nlcc).or..not.allocated(ppn%core23_tau).or. &
         .not.allocated(ppn%core23_current)) bad=1
    endif
    if(allocated(ppn%rho_nlcc)) then
      if(any(lbound(ppn%rho_nlcc)/=mg%is).or.any(ubound(ppn%rho_nlcc)/=mg%ie)) bad=1
    endif
    if(allocated(ppn%core23_tau)) then
      if(any(lbound(ppn%core23_tau)/=mg%is).or.any(ubound(ppn%core23_tau)/=mg%ie)) bad=1
    endif
    if(allocated(ppn%core23_current)) then
      if(any(lbound(ppn%core23_current)/=[1,1,mg%is]).or. &
         any(ubound(ppn%core23_current)/=[4,3,mg%ie])) bad=1
    endif
    if(allocated(ppn%tau_nlcc)) then
      if(any(lbound(ppn%tau_nlcc)/=mg%is).or.any(ubound(ppn%tau_nlcc)/=mg%ie)) bad=1
    endif
    call comm_summation(bad,global_bad,info%icomm_rko)
    if(global_bad/=0) error stop 'nc_mgga23: incomplete V2 core tau/current metadata or grid layout'
    do z=mg%is(3),mg%ie(3)
    do y=mg%is(2),mg%ie(2)
    do x=mg%is(1),mg%ie(1)
      if(ppn%core23_complete) then
        call core23_fields(ppn%rho_nlcc(x,y,z),ppn%core23_tau(x,y,z), &
          ppn%core23_current(:,:,x,y,z),system%vec_Ac,cn,ct,cj,ierr)
        if(ierr/=0) bad=1
      else
        if(allocated(ppn%rho_nlcc)) then
          if(.not.ieee_is_finite(ppn%rho_nlcc(x,y,z))) bad=1
          if(ppn%rho_nlcc(x,y,z)/=0d0) bad=1
        endif
        if(allocated(ppn%tau_nlcc)) then
          if(.not.ieee_is_finite(ppn%tau_nlcc(x,y,z))) bad=1
          if(ppn%tau_nlcc(x,y,z)/=0d0) bad=1
        endif
        if(allocated(ppn%core23_tau)) then
          if(.not.ieee_is_finite(ppn%core23_tau(x,y,z))) bad=1
          if(ppn%core23_tau(x,y,z)/=0d0) bad=1
        endif
        if(allocated(ppn%core23_current)) then
          do c=1,3
          do a=1,4
            if(.not.ieee_is_finite(ppn%core23_current(a,c,x,y,z))) bad=1
            if(ppn%core23_current(a,c,x,y,z)/=0d0) bad=1
          enddo
          enddo
        endif
      endif
    enddo
    enddo
    enddo
    call comm_summation(bad,global_bad,info%icomm_rko)
    if(global_bad/=0) error stop 'nc_mgga23: consistent core tau/current model requires valid complete V2 metadata'
    allocate(fields(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),20))
    allocate(total,mold=fields); fields=0d0
    allocate(gp(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3),2))
    if(info%if_divide_rspace) call update_overlap_complex8(srg,mg,spsi%zwf)
    do ik=info%ik_s,info%ik_e
      k=system%vec_k(:,ik)
      do io=info%io_s,info%io_e
        occ=system%rocc(io,ik,1)*system%wtk(ik)
        do sp=1,2
          call calc_gradient_psi(spsi%zwf(:,:,:,sp,io,ik,1),gp(:,:,:,:,sp), &
            mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
        end do
        do z=mg%is(3),mg%ie(3)
        do y=mg%is(2),mg%ie(2)
        do x=mg%is(1),mg%ie(1)
          p=spsi%zwf(x,y,z,:,io,ik,1)
          fields(x,y,z,1:4)=fields(x,y,z,1:4)+occ*real(traces(p,p),8)
          do c=1,3
            dp=gp(c,x,y,z,:)+zi*k(c)*p
            fields(x,y,z,5:8)=fields(x,y,z,5:8)+occ*real(traces(dp,dp),8)
            fields(x,y,z,8+4*(c-1)+1:8+4*c)= &
              fields(x,y,z,8+4*(c-1)+1:8+4*c)+occ*aimag(traces(p,dp))
          end do
        end do
        end do
        end do
      end do
    end do
    call comm_summation(fields,total,size(fields),info%icomm_ko)
    allocate(valence,source=total)
    if(ppn%core23_complete) then
      do z=mg%is(3),mg%ie(3)
      do y=mg%is(2),mg%ie(2)
      do x=mg%is(1),mg%ie(1)
        call core23_fields(ppn%rho_nlcc(x,y,z),ppn%core23_tau(x,y,z), &
          ppn%core23_current(:,:,x,y,z),system%vec_Ac,cn,ct,cj,ierr)
        total(x,y,z,1:4)=total(x,y,z,1:4)+cn
        total(x,y,z,5:8)=total(x,y,z,5:8)+ct
        do c=1,3
          total(x,y,z,8+4*(c-1)+1:8+4*c)=total(x,y,z,8+4*(c-1)+1:8+4*c)+cj(:,c)
        enddo
      enddo
      enddo
      enddo
    endif
    allocate(work(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3)))
    allocate(gwork(3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(lwork(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(grad(4,3,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(lap(4,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(dg,mold=grad); allocate(dl,mold=lap)
    do a=1,4
      work=0d0
      work(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))=total(:,:,:,a)
      if(info%if_divide_rspace) call update_overlap_real8(srg_scalar,mg,work)
      call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,work,gwork)
      call calc_laplacian_field(mg,stencil,work,lwork)
      do c=1,3
        grad(a,c,:,:,:)=gwork(c,:,:,:)
      end do
      lap(a,:,:,:)=lwork
    end do
    if(allocated(system%nc_mgga23_payload%local)) then
      deallocate(system%nc_mgga23_payload%local,system%nc_mgga23_payload%tau,system%nc_mgga23_payload%current)
    end if
    allocate(system%nc_mgga23_payload%local(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),4))
    allocate(system%nc_mgga23_payload%tau(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),4))
    allocate(system%nc_mgga23_payload%current(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),4,3))
    esum=0d0; bad=0;deda_local=0d0
    do z=mg%is(3),mg%ie(3)
    do y=mg%is(2),mg%ie(2)
    do x=mg%is(1),mg%ie(1)
      n=total(x,y,z,1:4);t=total(x,y,z,5:8)
      g=grad(:,:,x,y,z);l=lap(:,x,y,z)
      do c=1,3
        j(:,c)=total(x,y,z,8+4*(c-1)+1:8+4*c)
      end do
      call nc_mgga23_eval(n,g,l,t,j,e,dn,dgrad,dlap,dt,dj,ierr)
      if(ierr/=0) bad=1
      if(ppn%core23_complete) then
        cn=0d0;cn(1)=ppn%rho_nlcc(x,y,z)
        call core23_energy_deda(cn,ppn%core23_current(:,:,x,y,z),system%vec_Ac,dt,dj,deda)
        deda_local=deda_local+deda
      endif
      system%nc_mgga23_payload%local(x,y,z,:)=dn
      system%nc_mgga23_payload%tau(x,y,z,:)=dt
      system%nc_mgga23_payload%current(x,y,z,:,:)=dj
      dg(:,:,x,y,z)=dgrad;dl(:,x,y,z)=dlap
      esum=esum+e
      if(present(eexc)) eexc%f(x,y,z)=e
    end do
    end do
    end do
    call comm_summation(bad,global_bad,info%icomm_rko)
    if(global_bad/=0) error stop 'nc_mgga23 kernel rejected density/invariants'
    ! Euler derivative of E[n,grad n,lap n]: dn - div(dgrad) + lap(dlap).
    do a=1,4
      do c=1,3
        work=0d0
        work(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))=dg(a,c,:,:,:)
        if(info%if_divide_rspace) call update_overlap_real8(srg_scalar,mg,work)
        call calc_gradient_field(mg,stencil%coef_nab,system%rmatrix_B,work,gwork)
        system%nc_mgga23_payload%local(:,:,:,a)=system%nc_mgga23_payload%local(:,:,:,a)-gwork(c,:,:,:)
      end do
      work=0d0
      work(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))=dl(a,:,:,:)
      if(info%if_divide_rspace) call update_overlap_real8(srg_scalar,mg,work)
      call calc_laplacian_field(mg,stencil,work,lwork)
      system%nc_mgga23_payload%local(:,:,:,a)=system%nc_mgga23_payload%local(:,:,:,a)+lwork
    end do
    ! Occupied-state expectation of the complete variational XC operator.
    band_local=0d0
    do z=mg%is(3),mg%ie(3)
    do y=mg%is(2),mg%ie(2)
    do x=mg%is(1),mg%ie(1)
      n=valence(x,y,z,1:4)
      band_local=band_local+sum(system%nc_mgga23_payload%local(x,y,z,:)*n) &
        +sum(system%nc_mgga23_payload%tau(x,y,z,:)*valence(x,y,z,5:8))
      do c=1,3
        band_local=band_local+sum(system%nc_mgga23_payload%current(x,y,z,:,c) &
          *valence(x,y,z,8+4*(c-1)+1:8+4*c))
      end do
    end do
    end do
    end do
    call comm_summation(band_local,band_global,info%icomm_r)
    system%nc_mgga23_payload%band_expectation=band_global*system%hvol
    call comm_summation(deda_local,deda_global,3,info%icomm_r)
    system%nc_mgga23_payload%core_deda=deda_global*system%hvol
    call comm_summation(esum,eglob,info%icomm_r)
    E_xc=eglob*system%hvol
    system%nc_mgga23_payload%active=.true.
    system%xc_payload%use_tau_operator=.false.
  end subroutine

  subroutine nc_mgga23_add_action(system,mg,info,tpsi,htpsi,stencil,srg)
    implicit none
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: tpsi
    type(s_orbital),intent(inout) :: htpsi
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    complex(8),allocatable :: gp(:,:,:,:,:),flux(:,:,:,:,:,:,:,:),df(:,:,:,:)
    complex(8) :: p(2),dp(2),q(2),f(2,3)
    real(8) :: k(3)
    integer :: x,y,z,c,d,sp,io,ik,im
    if(.not.system%nc_mgga23_payload%active) return
    allocate(gp(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3),2))
    allocate(df(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3)))
    allocate(flux(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
      mg%is_array(3):mg%ie_array(3),2,info%io_s:info%io_e,info%ik_s:info%ik_e, &
      info%im_s:info%im_e,3)); flux=0d0
    do im=info%im_s,info%im_e
    do ik=info%ik_s,info%ik_e
      k=system%vec_k(:,ik)
      do io=info%io_s,info%io_e
        do sp=1,2
          call calc_gradient_psi(tpsi%zwf(:,:,:,sp,io,ik,im),gp(:,:,:,:,sp), &
            mg%is_array,mg%ie_array,mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
        end do
        do z=mg%is(3),mg%ie(3)
        do y=mg%is(2),mg%ie(2)
        do x=mg%is(1),mg%ie(1)
          p=tpsi%zwf(x,y,z,:,io,ik,im)
          q=act(system%nc_mgga23_payload%local(x,y,z,:),p)
          do c=1,3
            dp=gp(c,x,y,z,:)+zi*k(c)*p
            ! H_j = -i/2 [J D + D J]. Flux is T D psi + i J psi/2.
            f(:,c)=act(system%nc_mgga23_payload%tau(x,y,z,:),dp) &
              +0.5d0*zi*act(system%nc_mgga23_payload%current(x,y,z,:,c),p)
            q=q-0.5d0*zi*act(system%nc_mgga23_payload%current(x,y,z,:,c),dp)-zi*k(c)*f(:,c)
          end do
          do d=1,3
            flux(x,y,z,:,io,ik,im,d)=matmul(f,system%rmatrix_B(d,:))
          end do
          htpsi%zwf(x,y,z,:,io,ik,im)=htpsi%zwf(x,y,z,:,io,ik,im)+q
        end do
        end do
        end do
      end do
    end do
    end do
    do d=1,3
      if(info%if_divide_rspace) call update_overlap_complex8(srg,mg,flux(:,:,:,:,:,:,:,d))
      do im=info%im_s,info%im_e
      do ik=info%ik_s,info%ik_e
      do io=info%io_s,info%io_e
      do sp=1,2
        ! Flux is already transformed into grid coordinates: use identity here.
        call calc_gradient_psi(flux(:,:,:,sp,io,ik,im,d),df,mg%is_array,mg%ie_array, &
          mg%is,mg%ie,mg%idx,mg%idy,mg%idz,stencil%coef_nab,identity())
        htpsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),sp,io,ik,im)= &
          htpsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),sp,io,ik,im) &
          -df(d,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))
      end do
      end do
      end do
      end do
    end do
  end subroutine
  function identity() result(b)
    implicit none
    real(8) :: b(3,3)
    integer :: i
    b=0d0
    do i=1,3
      b(i,i)=1d0
    end do
  end function
end module
