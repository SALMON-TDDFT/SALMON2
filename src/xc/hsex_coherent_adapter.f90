#include "config.h"
module hsex_coherent_adapter
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
  implicit none
  private
  public :: coherent_snapshot
  abstract interface
    subroutine additional_action(sys,mg,info,x,hx)
      import s_dft_system,s_rgrid,s_parallel_info,s_orbital
      implicit none
      type(s_dft_system),intent(in) :: sys
      type(s_rgrid),intent(in) :: mg
      type(s_parallel_info),intent(in) :: info
      type(s_orbital),intent(in) :: x
      type(s_orbital),intent(inout) :: hx
    end subroutine
    subroutine shifted_action(index,x,shift,y,status)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: x(:)
      real(8),intent(in) :: shift
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine additional_derivative(index,x,y,status)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine coherent_snapshot(sys,mg,info,states,vl,stencil,srg,ppg,extra,extra_velocity,extra_mass,extra_shift)
    use structures, only: s_scalar,s_stencil,s_pp_grid,allocate_orbital_complex,deallocate_orbital
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi,projection_action,update_kvector_nonlocalpt
    use communication, only: comm_is_root,comm_get_max,comm_summation,comm_bcast
    use hsex_screen_math, only: screen_finite,screen_shift
    use hsex_frozen_response, only: frozen_background_rate
    use hsex_coherent_optical, only: coherent_optical,stationary_current_reference
    use hsex_grid_velocity, only: grid_velocity,grid_mass
    use salmon_global, only: nelec,hsex_w_eta,hsex_frozen_rank,hsex_frozen_tol,hsex_coherent_omega,hsex_coherent_axis
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    type(s_scalar),intent(in) :: vl(:)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    procedure(additional_action) :: extra
    procedure(additional_derivative) :: extra_velocity,extra_mass
    procedure(shifted_action) :: extra_shift
    integer,parameter :: nz=4,nb=4
    type(s_orbital) :: x,hx
    type(s_parallel_info) :: target_info
    type(s_dft_system) :: shifted_sys(2)
    type(s_pp_grid) :: shifted_pp(2)
    complex(8),allocatable :: psi(:,:),hp(:,:),vertices(:,:),seed(:,:),plus(:),minus(:),dv(:),dm(:)
    complex(8) :: z(nz),sigma(nz),density(nb,nz),local(nz),total(nz),dlocal(nb,nz),dtotal(nb,nz),epsilon(nz)
    real(8) :: residual(2,nz),maximum(2*nz+4),global_max(2*nz+4),r(3),rate,seconds,shift
    real(8) :: reference_local(2),reference_total(2),reference_value(2)
    real(8) :: translation_local(12),translation_max(12)
    real(8),parameter :: delta=1d-4
    integer :: n(3),ng,nv,ik,band,g,j,code,local_actions,actions,all_actions(2),work_actions(2),unit,ios,side,nreference
    integer(8) :: started,finished,clock_rate
    logical :: exists
    call system_clock(started,clock_rate)
    n=mg%ie-mg%is+1;ng=product(n);nv=nelec/2;code=0
    if(nv<1.or.sys%no/=nv.or.size(states%zwf,4)/=1.or.associated(projection_action))code=1
    if(.not.screen_finite(sys%wtk).or..not.screen_finite(sys%rocc(:,:,1)))then
      code=1
    else
      if(maxval(abs(sys%wtk-1d0/sys%nk))>1d-12.or.maxval(abs(sys%rocc(:,:,1)-2d0))>1d-12)code=1
    endif
    if(comm_is_root(info%id_ko))then
      inquire(file='hsex-coherent-optical.dat',exist=exists)
      if(exists)code=1
    endif
    call comm_get_max(code,info%icomm_ko)
    if(code/=0)error stop 'invalid coherent optical snapshot or existing output'
    allocate(psi(ng,nv),hp(ng,nv),vertices(ng,nb),seed(ng,6*nv),plus(ng),minus(ng),dv(ng),dm(ng),stat=code)
    call comm_get_max(code,info%icomm_ko)
    if(code/=0)error stop 'cannot allocate coherent optical snapshot'
    vertices(:,1)=1d0
    do g=0,ng-1
      r=dble([mod(g,n(1)),mod(g/n(1),n(2)),g/(n(1)*n(2))])/n
      vertices(g+1,2:)=exp(cmplx(0d0,2d0*acos(-1d0),8)*r)
    enddo
    z=[cmplx(0d0,hsex_w_eta,8),cmplx(0d0,hsex_w_eta/2d0,8), &
      cmplx(hsex_coherent_omega,hsex_w_eta,8),cmplx(hsex_coherent_omega,hsex_w_eta/2d0,8)]
    call allocate_orbital_complex(sys%nspin,mg,info,x);call allocate_orbital_complex(sys%nspin,mg,info,hx)
    do side=1,2
      shifted_sys(side)=sys;shifted_pp(side)=ppg
      shifted_sys(side)%vec_Ac(hsex_coherent_axis)=sys%vec_Ac(hsex_coherent_axis)+merge(-delta,delta,side==1)
      call update_kvector_nonlocalpt(info%ik_s,info%ik_e,shifted_sys(side),shifted_pp(side))
    enddo
    local=0d0;dlocal=0d0;maximum=0d0;local_actions=0;reference_local=0d0;nreference=0
    translation_local=0d0
    do ik=info%ik_s,info%ik_e
      psi=reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3), &
        1,:,ik,1),[ng,nv])*sqrt(sys%hvol)
      do band=1,nv
        call apply(psi(:,band),hp(:,band),ios)
        if(ios/=0)exit
        call shifted(psi(:,band),2,plus,ios)
        if(ios/=0)exit
        call shifted(psi(:,band),1,minus,ios)
        if(ios/=0)exit
        call velocity(psi(:,band),dv,ios)
        if(ios/=0)exit
        call mass(psi(:,band),dm,ios)
        if(ios/=0)exit
        maximum(2*nz+3)=max(maximum(2*nz+3),sqrt(sum(abs(dv-(plus-minus)/(2d0*delta))**2))/ &
          max(1d-12,sqrt(sum(abs(dv)**2))))
        maximum(2*nz+4)=max(maximum(2*nz+4),sqrt(sum(abs(dm-(plus-2d0*hp(:,band)+minus)/delta**2)**2))/ &
          max(1d-12,sqrt(sum(abs(dm)**2))))
        call translation_probe(psi(:,band),ios)
        if(ios/=0)exit
      enddo
      if(ios/=0)then
        code=1;exit
      endif
      call frozen_background_rate(psi,hp,rate,ios)
      if(ios/=0)then
        code=1;exit
      endif
      maximum(2*nz+2)=max(maximum(2*nz+2),rate)
      if(rate<=1d-8)then
        call stationary_current_reference(apply,velocity,mass,psi,2d0,sys%hvol*ng*sys%nk,1d-8, &
          reference_value(1),reference_value(2),ios)
        if(ios==0)then
          reference_local=reference_local+reference_value;nreference=nreference+1
        else
          write(*,'(a,2i8)')'HSEX_COHERENT_DC_REFERENCE_SKIP k/status ',ik,ios
        endif
      endif
      do j=1,3
        seed(:,(2*j-2)*nv+1:(2*j-1)*nv)=psi*spread(vertices(:,j+1),2,nv)
        seed(:,(2*j-1)*nv+1:2*j*nv)=psi*spread(conjg(vertices(:,j+1)),2,nv)
      enddo
      call coherent_optical(apply,velocity,mass,psi,[(2d0,j=1,nv)],vertices,z,hsex_frozen_rank, &
        hsex_frozen_tol,sys%hvol*ng*sys%nk,sigma,density,residual,actions,ios,seed)
      if(ios/=0)then
        write(*,'(a,3i8)')'HSEX_COHERENT_REJECT rank/k/status ',info%id_ko,ik,ios
        code=1;exit
      endif
      local=local+sigma;dlocal=dlocal+density/(sys%hvol*ng*sys%nk);local_actions=local_actions+actions
      maximum(:nz)=max(maximum(:nz),residual(1,:))
      maximum(nz+1:2*nz)=max(maximum(nz+1:2*nz),residual(2,:))
      maximum(2*nz+1)=max(maximum(2*nz+1),maxval(abs(density(1,:))))
    enddo
    call deallocate_orbital(x);call deallocate_orbital(hx)
    call comm_get_max(code,info%icomm_ko)
    if(code/=0)error stop 'coherent optical diagnostic failed; W unchanged'
    call comm_summation(local,total,nz,info%icomm_ko)
    call comm_summation(dlocal,dtotal,nb*nz,info%icomm_ko)
    call comm_get_max(maximum,global_max,size(maximum),info%icomm_ko)
    call comm_get_max(translation_local,translation_max,size(translation_local),info%icomm_ko)
    work_actions=[local_actions,nreference];call comm_summation(work_actions,all_actions,2,info%icomm_ko)
    call comm_summation(reference_local,reference_total,2,info%icomm_ko)
    epsilon=1d0+cmplx(0d0,4d0*acos(-1d0),8)*total/z
    if(.not.screen_finite(total).or..not.screen_finite(dtotal).or..not.screen_finite(epsilon)) &
      error stop 'nonfinite coherent optical sum'
    if(global_max(2*nz+1)>max(1d-8,10d0*hsex_frozen_tol)) &
      error stop 'coherent optical uniform charge not conserved'
    if(global_max(2*nz+3)>1d-6.or.global_max(2*nz+4)>1d-4) &
      error stop 'coherent optical operator derivative finite-difference check failed'
    if(maxval(translation_max(:3))>1d-10) &
      error stop 'coherent optical integer translation control failed'
    ios=0
    if(comm_is_root(info%id_ko))then
      open(newunit=unit,file='hsex-coherent-optical.dat',status='new',iostat=ios)
      if(ios==0)then
        write(unit,'(a)')'# IP electric-impulse optical response; frozen H, evolving coherent background; W unchanged'
        write(unit,'(a)')'# omega eta Re_sigma Im_sigma Re_epsilon_IP Im_epsilon_IP residual_Rrho residual_RvRrho'
        do j=1,nz
          write(unit,'(8es24.15)',iostat=ios)real(z(j),8),aimag(z(j)),real(total(j),8),aimag(total(j)), &
            real(epsilon(j),8),aimag(epsilon(j)),global_max(j),global_max(nz+j)
          if(ios/=0)exit
        enddo
        if(ios==0)then
          write(unit,'(a)')'# density/E impulse: G index, frequency index, real, imaginary (G=0,100,010,001)'
          do g=1,nb
            do j=1,nz
              write(unit,'(2i6,2es24.15)',iostat=ios)g,j,real(dtotal(g,j),8),aimag(dtotal(g,j))
              if(ios/=0)exit
            enddo
            if(ios/=0)exit
          enddo
        endif
        close(unit,iostat=code)
        if(ios==0)ios=code
      endif
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot write coherent optical diagnostic'
    call system_clock(finished);seconds=dble(finished-started)/clock_rate
    if(comm_is_root(info%id_ko))then
      write(*,'(a,2i8,4es20.10)')'HSEX_COHERENT_OPTICAL_PASS axis/H_actions/residual/charge/rate/seconds ', &
        hsex_coherent_axis,all_actions(1),maxval(global_max(:2*nz)),global_max(2*nz+1:2*nz+2),seconds
      write(*,'(a,2es20.10)')'HSEX_COHERENT_DERIVATIVE_ERROR velocity/mass ',global_max(2*nz+3:2*nz+4)
      do j=1,2
        write(*,'(a,i2,6es20.10)')'HSEX_GRID_LOCAL_TRANSLATION integer_or_half/relative_defect/energy_Ha ', &
          j,translation_max(3*j-2:3*j),translation_max(6+3*j-2:6+3*j)
      enddo
      if(all_actions(2)==sys%nk)then
        write(*,'(a,2es20.10)')'HSEX_COHERENT_STATIONARY_REFERENCE D_raw/epsilon_IP_static ', &
          reference_total(1),1d0+4d0*acos(-1d0)*reference_total(2)
      else
        write(*,'(a,i8)')'HSEX_COHERENT_STATIONARY_REFERENCE_NOT_AVAILABLE k_count ',all_actions(2)
      endif
    endif
  contains
    ! Translate the actual frozen local field and occupied periodic orbitals.
    ! This isolates collocation aliasing; it is not a resampled nonlocal-PP
    ! translation test and is not used to subtract any physical Drude response.
    subroutine translation_probe(vector,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      integer,intent(out) :: ierr
      complex(8) :: potential(ng),translated_potential(ng),translated_vector(ng),translated_product(ng)
      real(8) :: displacement(3),norm,energy
      integer :: axis,kind,index
      potential=cmplx(reshape(vl(1)%f(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
        mg%is(3):mg%ie(3)),[ng]),0d0,8)
      norm=max(1d-14,sqrt(sum(abs(potential*vector)**2)))
      energy=real(dot_product(vector,potential*vector),8)
      do kind=1,2
        do axis=1,3
          displacement=0d0;displacement(axis)=merge(1d0,.5d0,kind==1)/n(axis)
          call screen_shift(n,displacement,potential,translated_potential,ierr,real_field=.true.)
          if(ierr/=0)return
          call screen_shift(n,displacement,vector,translated_vector,ierr)
          if(ierr/=0)return
          call screen_shift(n,displacement,potential*vector,translated_product,ierr)
          if(ierr/=0)return
          index=axis+3*(kind-1)
          translation_local(index)=max(translation_local(index), &
            sqrt(sum(abs(translated_potential*translated_vector-translated_product)**2))/norm)
          translation_local(index+6)=max(translation_local(index+6), &
            abs(real(dot_product(translated_vector,translated_potential*translated_vector),8)-energy))
        enddo
      enddo
    end subroutine

    subroutine shifted(vector,side,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      integer,intent(in) :: side
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      complex(8) :: added(size(vector))
      target_info=info;target_info%ik_s=ik;target_info%ik_e=ik;target_info%numk=1
      target_info%io_s=1;target_info%io_e=1;target_info%numo=1
      x%zwf=0d0
      x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,ik,1)=reshape(vector,n)
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,target_info,mg,vl,shifted_sys(side),stencil,srg,shifted_pp(side))
      action=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,ik,1),[ng])
      call extra_shift(ik,vector,merge(-delta,delta,side==1),added,ierr)
      if(ierr/=0)return
      action=action+added
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine apply(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      target_info=info;target_info%ik_s=ik;target_info%ik_e=ik;target_info%numk=1
      target_info%io_s=1;target_info%io_e=1;target_info%numo=1
      x%zwf=0d0
      x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,ik,1)=reshape(vector,n)
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,target_info,mg,vl,sys,stencil,srg,ppg)
      call extra(sys,mg,target_info,x,hx)
      action=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,ik,1),[ng])
      ierr=0
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine velocity(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      complex(8) :: v(size(vector),3),added(size(vector))
      call grid_velocity(sys,mg,stencil,ppg,ik,vector,v)
      call extra_velocity(ik,vector,added,ierr)
      if(ierr/=0)return
      action=v(:,hsex_coherent_axis)+added
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine mass(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      complex(8) :: added(size(vector))
      call grid_mass(sys,mg,ppg,ik,hsex_coherent_axis,vector,action,ierr)
      if(ierr/=0)return
      call extra_mass(ik,vector,added,ierr)
      if(ierr/=0)return
      action=action+added
      if(.not.screen_finite(action))ierr=1
    end subroutine
  end subroutine
end module
