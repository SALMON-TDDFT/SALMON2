#include "config.h"
! Matched optical Sternheimer diagnostic of an unexcited scalar RT reference.
module hsex_sternheimer_optical
  implicit none
  private
  public :: sternheimer_optical_probe
contains
  subroutine sternheimer_optical_probe(lg,mg,system,info,stencil,xc_func,srg,srg_scalar,pp,ppg,ppn, &
      states,vl,fg,poisson)
    use structures
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi,update_kvector_nonlocalpt
    use density_matrix, only: calc_density,calc_current
    use hartree_sub, only: hartree
    use salmon_xc, only: exchange_correlation
    use communication, only: comm_get_max,comm_summation,comm_bcast,comm_is_root
    use hsex_screen_math, only: screen_finite
    use hsex_grid_velocity, only: grid_velocity,grid_mass
    use hsex_coherent_optical, only: stationary_current_reference
    use hsex_equilibrium_q, only: equilibrium_q_response
    use hsex_sternheimer, only: sternheimer_solve
    use hsex_impulse_resolvent, only: impulse_model,impulse_build,impulse_resolve,impulse_trajectory
    use salmon_global, only: hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,hsex_coherent_axis,nenergy,de,nt,dt,num_kgrid
    use salmon_global, only: yn_fix_func,hsex_sternheimer_kernel
    implicit none
    type(s_rgrid),intent(in) :: lg,mg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_xc_functional),intent(in) :: xc_func
    type(s_sendrecv_grid),intent(inout) :: srg,srg_scalar
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(in) :: ppg
    type(s_pp_nlcc),intent(in) :: ppn
    type(s_orbital),intent(in) :: states
    type(s_scalar),intent(in) :: vl(:)
    type(s_reciprocal_grid),intent(in) :: fg
    type(s_poisson),intent(inout) :: poisson
    type(s_dft_system) :: ss,qsys
    type(s_parallel_info) :: one_info
    type(s_pp_grid) :: pp0,qpp
    type(s_orbital) :: trial,acted
    type(s_scalar) :: rs(1),rho,vh,vxc(1)
    type(impulse_model) :: model
    complex(8),allocatable :: psi(:,:,:),h0(:,:,:),vpsi(:,:,:),hp(:,:,:),source(:),work(:)
    complex(8),allocatable :: shifts(:),transfer(:),sigma(:),epsilon(:),v(:,:),mass(:)
    complex(8),allocatable :: plus(:),minus(:),vp(:),vm(:)
    complex(8),allocatable :: density_shifts(:)
    complex(8),allocatable :: density_chi(:,:)
    complex(8),allocatable :: occupied_rotation(:,:,:),ritz_work(:)
    real(8),allocatable :: occupied_energy(:,:),ritz_real(:)
    real(8),allocatable :: residual(:),linear_current(:)
    real(8),parameter :: difference=1d-4
    real(8) :: volume,local_contact,contact,rate,global_rate,unused,current(3,1,1),current_pm(2)
    real(8) :: current_error,derivative_error,global_error,scale,analytic_pm(2),finite_pm(2)
    real(8) :: rate_result(1)
    real(8) :: reference_local(2),reference_total(2),reference(2)
    integer :: n(3),ng,nv,nk,nlocal,ns,ik,k,a,code,status,mode,unit,j,side,offset,reference_k,nprobe,generator_calls
    real(8) :: step_error,step_scale
    real(8) :: tangent_difference
    real(8) :: density_residual
    integer :: q_source,q_target,qpoint(3),qmode,direction
    logical :: self_consistent
    character(16) :: label
    n=mg%ie-mg%is+1;ng=product(n);nv=system%no;nk=info%numk;nlocal=ng*nv*nk;ns=nenergy+1
    ss=system;ss%vec_Ac=0d0;pp0=ppg
    call update_kvector_nonlocalpt(info%ik_s,info%ik_e,ss,pp0)
    if(nv<1.or.system%nspin/=1.or.info%isize_r/=1.or.info%isize_o/=1) &
      error stop 'optical Sternheimer needs scalar k-only MPI'
    if(maxval(abs(system%wtk-1d0/system%nk))>1d-12.or.maxval(abs(system%rocc(:,:,1)-2d0))>1d-12) &
      error stop 'optical Sternheimer needs uniform full k mesh and occupied-only reference'
    volume=ng*system%hvol
    allocate(psi(ng,nv,nk),h0(nv,nv,nk),vpsi(ng,nv,nk),hp(ng,nv,nk),source(nlocal),work(nlocal), &
      plus(ng),minus(ng),vp(ng),vm(ng),v(ng,3),mass(ng),shifts(ns),transfer(ns),sigma(ns),epsilon(ns),residual(ns))
    call allocate_orbital_complex(1,mg,info,trial);call allocate_orbital_complex(1,mg,info,acted)
    call allocate_scalar(mg,rs(1));call allocate_scalar(mg,rho);call allocate_scalar(mg,vh)
    call allocate_scalar(mg,vxc(1))
    psi=reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
      1,:,info%ik_s:info%ik_e,1),[ng,nv,nk])*sqrt(system%hvol)
    do j=1,ns
      shifts(j)=cmplx(hsex_w_eta,-(j-1)*de,8)
    enddo
    call fixed_action(reshape(psi,[nlocal]),work)
    hp=reshape(work,[ng,nv,nk]);rate=0d0;local_contact=0d0
    do k=1,nk
      ik=info%ik_s+k-1;h0(:,:,k)=matmul(conjg(transpose(psi(:,:,k))),hp(:,:,k))
      rate=max(rate,sqrt(sum(abs(hp(:,:,k)-matmul(psi(:,:,k),h0(:,:,k)))**2)))
      do a=1,nv
        call grid_velocity(ss,mg,stencil,pp0,ik,psi(:,a,k),v)
        vpsi(:,a,k)=v(:,hsex_coherent_axis)
        call grid_mass(ss,mg,pp0,ik,hsex_coherent_axis,psi(:,a,k),mass,code)
        if(code/=0)error stop 'optical Sternheimer mass action failed'
        local_contact=local_contact+2d0*real(dot_product(psi(:,a,k),mass),8)/(volume*system%nk)
      enddo
    enddo
    call comm_get_max([rate],rate_result,1,info%icomm_ko);global_rate=rate_result(1)
    call comm_summation(local_contact,contact,info%icomm_ko)
    if(comm_is_root(info%id_ko))write(*,'(a,2es20.10)')'STERNHEIMER_OPTICAL_REFERENCE rate/contact ',global_rate,contact
    if(global_rate>1d-6)error stop 'optical Sternheimer reference not stationary'
    ! Validate the native current at fixed orbitals against its A derivative.
    do side=1,2
      ss%vec_Ac(hsex_coherent_axis)=merge(-difference,difference,side==1)
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,ss,pp0)
      call unpack(reshape(psi,[nlocal]))
      call calc_current(ss,mg,stencil,info,srg,trial,pp0,current)
      current_pm(side)=current(hsex_coherent_axis,1,1)
    enddo
    ss%vec_Ac=0d0;call update_kvector_nonlocalpt(info%ik_s,info%ik_e,ss,pp0)
    current_error=abs((current_pm(2)-current_pm(1))/(2d0*difference)-contact)/max(abs(contact),1d-12)
    if(current_error>1d-6)error stop 'optical Sternheimer contact current check failed'
    reference_local=0d0
    do k=1,nk
      reference_k=info%ik_s+k-1
      call stationary_current_reference(reference_action,reference_velocity,reference_mass,psi(:,:,k), &
        2d0,volume*system%nk,1d-7,reference(1),reference(2),status)
      if(status/=0)error stop 'optical Sternheimer static density reference failed'
      reference_local=reference_local+reference
    enddo
    call comm_summation(reference_local,reference_total,2,info%icomm_ko)
    allocate(occupied_rotation(nv,nv,nk),occupied_energy(nv,nk),ritz_work(4*nv),ritz_real(3*nv))
    do k=1,nk
      occupied_rotation(:,:,k)=h0(:,:,k)
      call zheev('V','L',nv,occupied_rotation(:,:,k),nv,occupied_energy(:,k),ritz_work,size(ritz_work),ritz_real,status)
      if(status/=0)error stop 'occupied rotation for Sternheimer preconditioning failed'
    enddo
    if(comm_is_root(info%id_ko))write(*,'(a,2es20.10)')'STERNHEIMER_DENSITY_OPTICAL_REFERENCE D_raw/eps_IP ', &
      reference_total(1),1d0+4d0*acos(-1d0)*reference_total(2)
    if(hsex_sternheimer_kernel=='reference')then
      if(comm_is_root(info%id_ko))then
        open(newunit=unit,file='sternheimer-reference.dat',status='new',iostat=status)
        if(status/=0)error stop 'cannot write stationary current reference'
        write(unit,'(a)')'# D_raw eps_IP_static contact reference_rate_Ha'
        write(unit,'(4es24.15)')reference_total(1),1d0+4d0*acos(-1d0)*reference_total(2),contact,global_rate
        close(unit)
      endif
      call deallocate_orbital(trial);call deallocate_orbital(acted)
      return
    endif
    tangent_difference=difference
    if(yn_fix_func/='y'.and.(hsex_sternheimer_kernel=='both'.or. &
        hsex_sternheimer_kernel=='selfconsistent'))call select_tangent_step()
    if(hsex_sternheimer_kernel=='both'.or.hsex_sternheimer_kernel=='density')then
    nprobe=min(201,ns);allocate(density_shifts(nprobe),density_chi(nprobe,3))
    do j=1,nprobe
      density_shifts(j)=cmplx((j-1)*de*(ns-1)/(nprobe-1),hsex_w_eta,8)
    enddo
    qsys=ss;qpp=pp0;q_source=-1;q_target=-1
    do qmode=0,1
      qpoint=[qmode,0,0]
      call equilibrium_q_response(n,num_kgrid,qpoint,ss%primitive_a,ss%vec_k(:,:ss%nk),psi,h0, &
        info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,min(1024,hsex_frozen_rank), &
        hsex_frozen_tol,density_shifts,sector_action,density_chi,density_residual,status)
      if(status/=0)error stop 'equilibrium density chi0 diagnostic failed'
      if(comm_is_root(info%id_ko))then
        label=merge('q0              ','qfinite         ',qmode==0)
        open(newunit=unit,file='sternheimer-density-'//trim(label)//'.dat',status='new',iostat=status)
        if(status/=0)error stop 'cannot write equilibrium density chi0'
        write(unit,'(a)')'# trial omega_Ha eta_Ha Re_chi0 Im_chi0 true_residual; q0 G=100,010,001; finite G=0,100,010'
        do a=1,3
          do j=1,nprobe
            write(unit,'(i6,5es24.15)')a,real(density_shifts(j),8),aimag(density_shifts(j)), &
              real(density_chi(j,a),8),aimag(density_chi(j,a)),density_residual
          enddo
        enddo
        close(unit)
      endif
    enddo
    endif
    if(hsex_sternheimer_kernel=='density')then
      call deallocate_orbital(trial);call deallocate_orbital(acted)
      return
    endif
    do mode=2,1,-1
      if(mode==1.and.yn_fix_func=='y')cycle
      if(mode==1.and.hsex_sternheimer_kernel=='independent')cycle
      if(mode==2.and.hsex_sternheimer_kernel=='selfconsistent')cycle
      generator_calls=0;step_error=0d0
      self_consistent=mode==1;label=merge('selfconsistent  ','independent     ',self_consistent)
      source=reshape(vpsi,[nlocal])
      if(self_consistent)then
        call potential(reshape(psi,[nlocal]),-difference,plus)
        call potential(reshape(psi,[nlocal]),difference,minus)
        do k=1,nk
          do a=1,nv
            offset=ng*((k-1)*nv+a-1)
            source(offset+1:offset+ng)=source(offset+1:offset+ng)+(minus-plus)/(2d0*difference)*psi(:,a,k)
          enddo
        enddo
      endif
      source=cmplx(0d0,-1d0,8)*source;call project(source)
      ! The impulse direction is current-null by time reversal. Check its
      ! quadrature too, and scale both errors by the nonzero current derivative.
      do direction=1,2
        work=source
        if(direction==2)work=cmplx(0d0,1d0,8)*source
        call observation(work,analytic_pm(direction),code)
        if(code/=0)error stop 'optical Sternheimer current observation failed'
        do side=1,2
          call unpack(reshape(psi,[nlocal])+merge(-difference,difference,side==1)*work)
          call calc_current(ss,mg,stencil,info,srg,trial,pp0,current)
          current_pm(side)=current(hsex_coherent_axis,1,1)
        enddo
        finite_pm(direction)=(current_pm(2)-current_pm(1))/(2d0*difference)
      enddo
      if(.not.screen_finite(analytic_pm).or..not.screen_finite(finite_pm)) &
        error stop 'nonfinite optical Sternheimer current derivative'
      derivative_error=maxval(abs(finite_pm-analytic_pm))/max(maxval(abs(analytic_pm)),1d-8)
      if(comm_is_root(info%id_ko))write(*,'(a,5es20.10)')'STERNHEIMER_CURRENT_DERIVATIVE_CHECK ', &
        analytic_pm,finite_pm,derivative_error
      if(derivative_error>1d-5)error stop 'optical Sternheimer orbital current check failed'
      call impulse_build(generator,observation,source,shifts,hsex_frozen_rank,hsex_frozen_tol, &
        info%icomm_ko,info%id_ko,model,status,precondition=static_precondition)
      if(status/=0)error stop 'optical Sternheimer tangent solve not converged'
      if(comm_is_root(info%id_ko))then
        call impulse_resolve(model,shifts,transfer,residual,status)
        if(status/=0.or.maxval(residual)>hsex_frozen_tol)error stop 'optical Sternheimer true residual'
        sigma=(contact+transfer)/shifts
        do j=1,ns
          epsilon(j)=1d0+4d0*acos(-1d0)*sigma(j)/shifts(j)
        enddo
        open(newunit=unit,file='sternheimer-'//trim(label)//'.dat',status='new',iostat=status)
        if(status/=0)error stop 'cannot open optical Sternheimer output'
        write(unit,'(a)')'# omega_Ha eta_Ha Re_sigma Im_sigma Re_eps Im_eps full_tangent_residual'
        do j=1,ns
          write(unit,'(7es24.15)')(j-1)*de,hsex_w_eta,real(sigma(j),8),aimag(sigma(j)), &
            real(epsilon(j),8),aimag(epsilon(j)),residual(j)
        enddo
        close(unit)
        open(newunit=unit,file='sternheimer-'//trim(label)//'-operator.bin',status='new', &
          access='stream',form='unformatted',iostat=status)
        if(status/=0)error stop 'cannot write Sternheimer operator record'
        write(unit)'HIO2',size(model%g),model%source_norm,model%next_norm,contact,model%h,model%g, &
          model%validation_residual,merge(1,0,allocated(model%residual_factor))
        if(allocated(model%residual_factor))write(unit)model%residual_factor
        close(unit)
        if(mode==2)then
          open(newunit=unit,file='sternheimer-density-optical-chi0.dat',status='new',iostat=status)
          if(status/=0)error stop 'cannot write density optical chi0'
          write(unit,'(a)')'# omega eta Re_eps_IP_density Im_eps_IP_density; velocity/dipole density head, no contact term'
          do j=1,ns
            epsilon(j)=1d0+4d0*acos(-1d0)*(contact+transfer(j)-reference_total(1))/shifts(j)**2
            write(unit,'(4es24.15)')(j-1)*de,hsex_w_eta,real(epsilon(j),8),aimag(epsilon(j))
          enddo
          close(unit)
        endif
        allocate(linear_current(0:nt))
        call impulse_trajectory(model,contact,dt,linear_current,status)
        if(status==0)then
          open(newunit=unit,file='sternheimer-'//trim(label)//'-time.dat',status='new',iostat=status)
          if(status/=0)error stop 'cannot open Sternheimer window companion'
          write(unit,'(a)')'# time_Ha_inverse J_per_A; finite-time companion of Sternheimer operator, not physical RT'
          do j=0,nt
            write(unit,'(2es24.15)')j*dt,linear_current(j)
          enddo
          close(unit)
          write(*,'(a,a)')'STERNHEIMER_WINDOW_PASS ',trim(label)
        else
          write(*,'(a,a,i8)')'STERNHEIMER_WINDOW_REJECT ',trim(label),status
        endif
        deallocate(linear_current)
        write(*,'(a,a,4es20.10)')'STERNHEIMER_OPTICAL_PASS ',trim(label), &
          maxval(residual),current_error,derivative_error,global_rate
        write(*,'(a,a,es20.10)')'STERNHEIMER_XC_DIRECTIONAL_STEP_ERROR ',trim(label),step_error
      endif
    enddo
    call deallocate_orbital(trial);call deallocate_orbital(acted)
  contains
    subroutine static_precondition(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      complex(8) :: block(ng,nv),rotated(ng,nv),solved(ng,nv)
      integer :: kk,band,lo,hi
      action=0d0;ierr=0
      do kk=1,nk
        reference_k=info%ik_s+kk-1;lo=(kk-1)*ng*nv+1;hi=kk*ng*nv
        block=reshape(vector(lo:hi),[ng,nv]);rotated=matmul(block,occupied_rotation(:,:,kk))
        do band=1,nv
          call sternheimer_solve(reference_action,psi(:,:,kk),occupied_energy(band,kk), &
            rotated(:,band),solved(:,band),ierr)
          if(ierr/=0)return
        enddo
        block=matmul(solved,conjg(transpose(occupied_rotation(:,:,kk))))
        action(lo:hi)=cmplx(0d0,-1d0,8)*reshape(block,[ng*nv])
      enddo
      if(.not.screen_finite(action))ierr=1
    end subroutine

    subroutine sector_action(source,target,vector,action,ierr)
      implicit none
      integer,intent(in) :: source,target
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      if(source/=q_source.or.target/=q_target)then
        qsys%vec_k(:,source)=ss%vec_k(:,target)
        call update_kvector_nonlocalpt(source,source,qsys,qpp)
        q_source=source;q_target=target
      endif
      one_info=info;one_info%ik_s=source;one_info%ik_e=source;one_info%numk=1
      one_info%io_s=1;one_info%io_e=1;one_info%numo=1
      trial%zwf=0d0
      trial%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,source,1)= &
        reshape(vector/sqrt(system%hvol),n)
      trial%update_zwf_overlap=.false.
      call hpsi(trial,acted,one_info,mg,vl,qsys,stencil,srg,qpp)
      action=reshape(acted%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
        1,1,source,1),[ng])*sqrt(system%hvol)
      ierr=0
      if(.not.screen_finite(action))ierr=1
    end subroutine

    subroutine single_action(index,vector,action,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      one_info=info;one_info%ik_s=index;one_info%ik_e=index;one_info%numk=1
      one_info%io_s=1;one_info%io_e=1;one_info%numo=1
      trial%zwf=0d0
      trial%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,index,1)= &
        reshape(vector/sqrt(system%hvol),n)
      trial%update_zwf_overlap=.false.
      call hpsi(trial,acted,one_info,mg,vl,ss,stencil,srg,pp0)
      action=reshape(acted%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
        1,1,index,1),[ng])*sqrt(system%hvol)
      ierr=0
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine reference_action(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      call single_action(reference_k,vector,action,ierr)
    end subroutine
    subroutine reference_velocity(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      call grid_velocity(ss,mg,stencil,pp0,reference_k,vector,v)
      action=v(:,hsex_coherent_axis);ierr=0
    end subroutine
    subroutine reference_mass(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      call grid_mass(ss,mg,pp0,reference_k,hsex_coherent_axis,vector,action,ierr)
    end subroutine
    subroutine unpack(vector)
      implicit none
      complex(8),intent(in) :: vector(:)
      trial%zwf=0d0
      trial%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
        1,:,info%ik_s:info%ik_e,1)=reshape(vector/sqrt(system%hvol),[n,nv,nk])
      trial%update_zwf_overlap=.false.
    end subroutine
    subroutine fixed_action(vector,action)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      call unpack(vector)
      call hpsi(trial,acted,info,mg,vl,ss,stencil,srg,pp0)
      action=reshape(acted%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
        1,:,info%ik_s:info%ik_e,1),[nlocal])*sqrt(system%hvol)
    end subroutine
    subroutine project(vector)
      implicit none
      complex(8),intent(inout) :: vector(:)
      complex(8) :: block(ng,nv)
      integer :: kk,lo,hi
      do kk=1,nk
        lo=(kk-1)*ng*nv+1;hi=kk*ng*nv;block=reshape(vector(lo:hi),[ng,nv])
        block=block-matmul(psi(:,:,kk),matmul(conjg(transpose(psi(:,:,kk))),block))
        vector(lo:hi)=reshape(block,[ng*nv])
      enddo
    end subroutine
    subroutine potential(vector,ac,values)
      implicit none
      complex(8),intent(in) :: vector(:)
      real(8),intent(in) :: ac
      complex(8),intent(out) :: values(:)
      ss%vec_Ac=0d0;ss%vec_Ac(hsex_coherent_axis)=ac
      call unpack(vector)
      call calc_density(ss,rs,trial,info,mg);rho%f=rs(1)%f
      call hartree(lg,mg,info,ss,fg,poisson,srg_scalar,stencil,rho,vh)
      call exchange_correlation(ss,xc_func,mg,srg_scalar,srg,rs,pp,ppn,info,trial,stencil,vxc,unused)
      values=cmplx(reshape(vh%f+vxc(1)%f,[ng]),0d0,8)
      ss%vec_Ac=0d0
      if(.not.screen_finite(values))error stop 'nonfinite native optical response potential'
    end subroutine
    subroutine select_tangent_step()
      implicit none
      real(8),parameter :: steps(8)=[1d-3,5d-4,1d-4,5d-5,1d-5,5d-6,1d-6,5d-7]
      real(8) :: local_norm,global_norm,error,scale,selected,best_error
      complex(8) :: previous(ng),derivative(ng)
      integer :: i
      source=cmplx(0d0,-1d0,8)*reshape(vpsi,[nlocal]);call project(source)
      local_norm=sum(abs(source)**2);call comm_summation(local_norm,global_norm,info%icomm_ko)
      source=source/sqrt(global_norm);selected=0d0;best_error=huge(1d0)
      do i=1,size(steps)
        call potential(reshape(psi,[nlocal])+steps(i)*source,0d0,plus)
        call potential(reshape(psi,[nlocal])-steps(i)*source,0d0,minus)
        derivative=(plus-minus)/(2d0*steps(i))
        if(mod(i,2)==0)then
          scale=max(1d-8,sqrt(sum(abs(derivative)**2)))
          error=sqrt(sum(abs(derivative-previous)**2))/scale
          if(comm_is_root(info%id_ko))write(*,'(a,3es20.10)')'STERNHEIMER_XC_STEP_SCAN delta/error/scale ', &
            steps(i-1),error,scale
          if(error<1d-4.and.error<best_error)then
            selected=steps(i-1);best_error=error
          endif
        endif
        previous=derivative
      enddo
      if(selected==0d0)error stop 'native TB-mBJ derivative has no converged step'
      tangent_difference=selected
      if(comm_is_root(info%id_ko))write(*,'(a,es20.10)')'STERNHEIMER_XC_SELECTED_STEP ',tangent_difference
    end subroutine
    subroutine generator(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      complex(8) :: block(ng,nv)
      complex(8) :: error_action(nlocal),kernel_action(nlocal),base_action(nlocal)
      real(8) :: local_norms(3),global_norms(3)
      integer :: kk,lo,hi
      call fixed_action(vector,action)
      if(self_consistent)then
        call potential(reshape(psi,[nlocal])+tangent_difference*vector,0d0,vp)
        call potential(reshape(psi,[nlocal])-tangent_difference*vector,0d0,vm)
        if(generator_calls<2)then
          call potential(reshape(psi,[nlocal])+.5d0*tangent_difference*vector,0d0,plus)
          call potential(reshape(psi,[nlocal])-.5d0*tangent_difference*vector,0d0,minus)
          do kk=1,nk
            lo=(kk-1)*ng*nv+1;hi=kk*ng*nv
            kernel_action(lo:hi)=reshape(spread((vp-vm)/(2d0*tangent_difference),2,nv)*psi(:,:,kk),[ng*nv])
            error_action(lo:hi)=kernel_action(lo:hi)- &
              reshape(spread((plus-minus)/tangent_difference,2,nv)*psi(:,:,kk),[ng*nv])
            block=reshape(vector(lo:hi),[ng,nv])
            base_action(lo:hi)=action(lo:hi)-reshape(matmul(block,h0(:,:,kk)),[ng*nv])
          enddo
          call project(error_action);call project(kernel_action);call project(base_action)
          local_norms=[sum(abs(error_action)**2),sum(abs(kernel_action)**2),sum(abs(base_action)**2)]
          call comm_summation(local_norms,global_norms,3,info%icomm_ko)
          ! Compare operator actions, using a componentwise backward-error
          ! scale. A density-null tangent must not fail solely because its
          ! exact induced potential is zero while the native evaluation rounds.
          step_scale=max(1d-8,sqrt(global_norms(2)),sqrt(global_norms(3)))
          step_error=max(step_error,sqrt(global_norms(1))/step_scale)
          if(comm_is_root(info%id_ko))write(*,'(a,3es20.10)')'STERNHEIMER_XC_DIRECTIONAL_CHECK relative/abs/scale ', &
            step_error,sqrt(global_norms(1)),step_scale
          if(step_error>1d-4)error stop 'TB-mBJ directional derivative step not converged'
        endif
      endif
      generator_calls=generator_calls+1
      do kk=1,nk
        lo=(kk-1)*ng*nv+1;hi=kk*ng*nv
        block=reshape(vector(lo:hi),[ng,nv])
        block=reshape(action(lo:hi),[ng,nv])-matmul(block,h0(:,:,kk))
        if(self_consistent)block=block+spread((vp-vm)/(2d0*tangent_difference),2,nv)*psi(:,:,kk)
        action(lo:hi)=reshape(cmplx(0d0,-1d0,8)*block,[ng*nv])
      enddo
      call project(action);ierr=0
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine observation(vector,value,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      real(8),intent(out) :: value
      integer,intent(out) :: ierr
      real(8) :: local
      local=4d0*real(dot_product(reshape(vpsi,[nlocal]),vector),8)/(volume*system%nk)
      call comm_summation(local,value,info%icomm_ko);ierr=0
      if(.not.screen_finite([value]))ierr=1
    end subroutine
  end subroutine
end module
