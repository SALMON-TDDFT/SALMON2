#include "config.h"
module hsex_frozen_adapter
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
  use hsex_frozen_full, only: s_frozen_window_cache
  implicit none
  private
  public :: frozen_snapshot_probe
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
    subroutine transfer_frame(index,owner,worker,slot,status)
      implicit none
      integer,intent(in) :: index,owner,worker,slot
      integer,intent(out) :: status
    end subroutine
    subroutine additional_velocity(index,vector,action,status)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine frozen_snapshot_probe(sys,mg,info,states,vl,stencil,srg,ppg,extra,boxes,first,nlocal,result, &
    extra_velocity,prepare_extra,apply_extra,window_cache,initialize_window)
    use structures, only: s_scalar,s_stencil,s_pp_grid,allocate_orbital_complex,deallocate_orbital
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi,projection_action,update_kvector_nonlocalpt
    use hsex_frozen_q, only: frozen_q_probe
    use hsex_frozen_full, only: frozen_full_w
    use communication, only: comm_get_max,comm_is_root
    use hsex_screen_math, only: screen_finite
    use hsex_grid_velocity, only: grid_velocity
    use update_kvector_plusU_sub, only: PLUS_U_ON
    use hsex_frozen_kinetic, only: frozen_kinetic_setup,frozen_kinetic_apply
    use salmon_global, only: nelec,num_kgrid,hsex_w_qgrid,hsex_w_eta,hsex_frozen_rank,hsex_frozen_tol, &
      hsex_frozen_qpoint,yn_hsex_frozen_reference,hsex_w_cutoff,hsex_w_eta_tol,hsex_w_head_tol, &
      yn_hsex_frozen_optical,hsex_frozen_gamma,hsex_frozen_parallel,hsex_frozen_nproc_q,hsex_frozen_precondition, &
      hsex_frozen_closure,hsex_frozen_window,hsex_frozen_window_time
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    type(s_scalar),intent(in) :: vl(:)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    procedure(additional_action),optional :: extra
    procedure(additional_velocity),optional :: extra_velocity
    procedure(transfer_frame),optional :: prepare_extra
    procedure(additional_velocity),optional :: apply_extra
    type(s_frozen_window_cache),optional,intent(inout) :: window_cache
    logical,optional,intent(in) :: initialize_window
    integer,optional,intent(in) :: boxes,first,nlocal
    real(8),allocatable,optional,intent(out) :: result(:,:)
    type(s_orbital) :: x,hx
    type(s_parallel_info) :: target_info,scratch_info
    type(s_dft_system) :: scratch_sys(2)
    type(s_pp_grid) :: scratch_pp(2)
    real(8),allocatable :: kinetic(:,:)
    integer :: frame_k(2),slot
    complex(8),allocatable :: occupied(:,:,:)
    integer :: n(3),ng,nv,ik,code,filter
    logical :: ref_init
    filter=0;ref_init=.false.
    if(hsex_frozen_window=='gaussian')filter=1
    if(hsex_frozen_window=='exponential')filter=2
    if(present(initialize_window))ref_init=initialize_window
    if(filter/=0.and..not.present(window_cache))error stop 'missing frozen window reference cache'
    n=mg%ie-mg%is+1;ng=product(n);nv=nelec/2
    code=0
    if(nelec<2.or.mod(nelec,2)/=0)code=1
    if(associated(projection_action).or.size(states%zwf,4)/=1.or.sys%no<nv)code=1
    if(.not.screen_finite(sys%wtk).or..not.screen_finite(sys%rocc(:,:,1)))then
      code=1
    else
      if(maxval(abs(sys%wtk-1d0/sys%nk))>1d-12)code=1
      if(nv>0.and.sys%no>=nv)then
        if(maxval(abs(sys%rocc(:nv,:,1)-2d0))>1d-12)code=1
        if(size(sys%rocc,1)>nv)then
          if(maxval(abs(sys%rocc(nv+1:,:,1)))>1d-12)code=1
        endif
      endif
    endif
    call comm_get_max(code,info%icomm_ko)
    if(code/=0)error stop 'frozen probe requires scalar orbitals and no active projection action'
    if(comm_is_root(info%id_ko))write(*,'(a,2i8,3es20.10)') &
      'HSEX_FROZEN_SNAPSHOT occupied/available/A ',nv,sys%no,sys%vec_Ac
    allocate(occupied(ng,nv,info%ik_s:info%ik_e))
    do ik=info%ik_s,info%ik_e
      occupied(:,:,ik)=reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
        mg%is(3):mg%ie(3),1,:nv,ik,1),[ng,nv])*sqrt(sys%hvol)
    enddo
    scratch_info=info
    if(hsex_frozen_parallel=='pairs')then
      if(PLUS_U_ON)error stop 'frozen pair MPI +U frame transfer is not validated'
      if(info%isize_r/=1.or.info%isize_o/=1)error stop 'frozen pair MPI requires k-only decomposition'
      if(present(extra).and.(.not.present(prepare_extra).or..not.present(apply_extra))) &
        error stop 'missing frozen HSEX frame callbacks'
      scratch_info%ik_s=1;scratch_info%ik_e=1;scratch_info%numk=1
      scratch_info%io_s=1;scratch_info%io_e=1;scratch_info%numo=1
      allocate(kinetic(ng,2));frame_k=0
      do slot=1,2
        scratch_sys(slot)=sys;scratch_pp(slot)=ppg
        scratch_sys(slot)%no=1
        if(allocated(scratch_pp(slot)%zekr_uv))deallocate(scratch_pp(slot)%zekr_uv)
      enddo
    endif
    call allocate_orbital_complex(sys%nspin,mg,scratch_info,x)
    call allocate_orbital_complex(sys%nspin,mg,scratch_info,hx)
    if(present(result))then
      if(.not.present(boxes).or..not.present(first).or..not.present(nlocal))error stop 'missing frozen W layout'
      allocate(result(nlocal,sys%nk))
      if(hsex_frozen_parallel=='pairs')then
        if(hsex_frozen_precondition=='kinetic')then
          call frozen_full_w(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a,sys%vec_k(:,:sys%nk), &
            occupied,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,first,nlocal, &
            hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,hsex_w_eta_tol,apply,result,code, &
            density_gamma=hsex_frozen_gamma=='q2',gamma_tol=hsex_w_head_tol, &
            prepare=prepare,nproc_q=hsex_frozen_nproc_q,precondition=precondition,closure=hsex_frozen_closure, &
            window_kind=filter,window_time=hsex_frozen_window_time,window_cache=window_cache,initialize_window=ref_init)
        else
          call frozen_full_w(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a,sys%vec_k(:,:sys%nk), &
            occupied,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,first,nlocal, &
            hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,hsex_w_eta_tol,apply,result,code, &
            density_gamma=hsex_frozen_gamma=='q2',gamma_tol=hsex_w_head_tol, &
            prepare=prepare,nproc_q=hsex_frozen_nproc_q,closure=hsex_frozen_closure, &
            window_kind=filter,window_time=hsex_frozen_window_time,window_cache=window_cache,initialize_window=ref_init)
        endif
      else if(yn_hsex_frozen_optical=='y')then
        if(present(extra).neqv.present(extra_velocity))error stop 'frozen optical H and velocity must match'
        call frozen_full_w(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a,sys%vec_k(:,:sys%nk), &
          occupied,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,first,nlocal, &
          hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,hsex_w_eta_tol,apply,result,code,velocity)
      else
        call frozen_full_w(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a,sys%vec_k(:,:sys%nk), &
          occupied,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,first,nlocal, &
          hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,hsex_w_eta_tol,apply,result,code, &
          density_gamma=hsex_frozen_gamma=='q2',gamma_tol=hsex_w_head_tol,closure=hsex_frozen_closure, &
            window_kind=filter,window_time=hsex_frozen_window_time,window_cache=window_cache,initialize_window=ref_init)
      endif
    else
      call frozen_q_probe(n,num_kgrid,hsex_w_qgrid,hsex_frozen_qpoint,sys%primitive_a,sys%vec_k(:,:sys%nk), &
        occupied,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko, &
        hsex_frozen_rank,hsex_frozen_tol,hsex_w_eta,yn_hsex_frozen_reference=='y',apply,code)
    endif
    call deallocate_orbital(x);call deallocate_orbital(hx)
    if(code/=0)then
      if(comm_is_root(info%id_ko))write(*,'(a,i8)')'HSEX_FROZEN_REJECT status=',code
      error stop 'frozen response failed; existing W preserved'
    endif
  contains
    subroutine prepare(source,target,worker,source_owner,target_owner,ierr)
      implicit none
      integer,intent(in) :: source,target,worker,source_owner,target_owner
      integer,intent(out) :: ierr
      integer :: side,index,owner
      ierr=0
      do side=1,2
        index=merge(source,target,side==1);owner=merge(source_owner,target_owner,side==1)
        if(present(prepare_extra))call prepare_extra(index,owner,worker,side,ierr)
        call comm_get_max(ierr,info%icomm_ko)
        if(ierr/=0)return
        if(info%id_ko==worker)then
          frame_k(side)=index;scratch_sys(side)%vec_k(:,1)=sys%vec_k(:,index)
          call update_kvector_nonlocalpt(1,1,scratch_sys(side),scratch_pp(side))
          if(hsex_frozen_precondition=='kinetic')then
            call frozen_kinetic_setup(n,sys%primitive_a,sys%vec_k(:,index)+sys%vec_Ac,kinetic(:,side),ierr)
          endif
        endif
        call comm_get_max(ierr,info%icomm_ko)
        if(ierr/=0)return
      enddo
    end subroutine
    subroutine precondition(index,vector,action,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      integer :: side
      ierr=1;side=0
      if(frame_k(1)==index)side=1
      if(frame_k(2)==index)side=2
      if(side==0)return
      call frozen_kinetic_apply(n,kinetic(:,side),vector,action,ierr)
    end subroutine
    subroutine velocity(index,vector,action,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:,:)
      integer,intent(out) :: ierr
      complex(8) :: added(size(vector),3)
      ierr=1
      if(index<info%ik_s.or.index>info%ik_e)return
      call grid_velocity(sys,mg,stencil,ppg,index,vector,action)
      if(present(extra_velocity))then
        call extra_velocity(index,vector,added,ierr)
        if(ierr/=0)return
        action=action+added
      endif
      ierr=0
      if(.not.screen_finite(action))ierr=1
    end subroutine
    subroutine apply(index,vector,action,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      integer :: side
      complex(8) :: added(size(vector),1)
      ierr=1
      if(hsex_frozen_parallel=='pairs')then
        side=0
        if(frame_k(1)==index)side=1
        if(frame_k(2)==index)side=2
        if(side==0)return
        x%zwf=0d0
        x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,1,1)=reshape(vector,n)
        x%update_zwf_overlap=.false.
        call hpsi(x,hx,scratch_info,mg,vl,scratch_sys(side),stencil,srg,scratch_pp(side))
        action=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,1,1),[ng])
        if(present(apply_extra))then
          call apply_extra(side,vector,added,ierr)
          if(ierr/=0)return
          action=action+added(:,1)
        endif
        ierr=0
        if(.not.screen_finite(action))ierr=1
        return
      endif
      if(index<info%ik_s.or.index>info%ik_e)return
      target_info=info
      target_info%ik_s=index;target_info%ik_e=index;target_info%numk=1
      target_info%io_s=1;target_info%io_e=1;target_info%numo=1
      x%zwf=0d0
      x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,index,1)=reshape(vector,n)
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,target_info,mg,vl,sys,stencil,srg,ppg)
      if(present(extra))call extra(sys,mg,target_info,x,hx)
      action=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,index,1),[ng])
      ierr=0
    end subroutine
  end subroutine
end module
