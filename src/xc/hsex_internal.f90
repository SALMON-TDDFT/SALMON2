#include "config.h"
! SALMON adapter: keep the Hamiltonian and projection manifold in memory.
module hsex_internal
  implicit none
  private
  public :: hsex_compute_internal,hsex_save_internal
  abstract interface
    subroutine extra_velocity_action(ik,vector,action,status)
      implicit none
      integer,intent(in) :: ik
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine hsex_compute_internal(sys,mg,info,states,x,hx,eigen,pop,nv,boxes,eta,vl,stencil,srg,ppg, &
                                   first,nlocal,result,extra_velocity,solver_override)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar,s_stencil,s_pp_grid
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi
    use communication, only: comm_get_max,comm_is_root
    use salmon_global, only: num_kgrid,hsex_w_qgrid,hsex_w_cutoff,hsex_w_head_tol,hsex_w_eta_tol,temperature,yn_hsex_head_only
    use salmon_global, only: hsex_w_solver
    use hsex_screen_driver, only: screen_compute
    use hsex_frozen_adapter, only: frozen_snapshot_probe
    use salmon_global, only: hsex_frozen_rank
    use hsex_pdep_q, only: pdep_q_test
    use hsex_reduced_q, only: reduced_q_test
    use salmon_global, only: hsex_reduced_target,hsex_reduced_qcount,hsex_reduced_maxrank,hsex_reduced_tol,hsex_reduced_batch, &
      yn_hsex_reduced_symmetry
    use salmon_global, only: hsex_pdep_rank,hsex_pdep_qpoint,hsex_pdep_reference
    use hsex_screen_math, only: screen_finite
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    type(s_orbital),intent(inout) :: x,hx
    real(8),intent(in) :: eigen(:,:),pop(:,:),eta
    integer,intent(in) :: nv,boxes,first,nlocal
    type(s_scalar),intent(in) :: vl(:)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    procedure(extra_velocity_action),optional :: extra_velocity
    character(*),optional,intent(in) :: solver_override
    character(16) :: solver
    real(8),allocatable,intent(out) :: result(:,:)
    complex(8),allocatable :: screen_states(:,:,:)
    type(s_parallel_info) :: target_info
    integer :: n(3),ng,nreal,ncomponent,ik,status,threads,astat,nkeep,band
    real(8) :: dv
    solver=hsex_w_solver
    if(hsex_frozen_rank>0)then
      if(present(extra_velocity).or.present(solver_override))error stop 'invalid frozen GS diagnostic dispatch'
      call frozen_snapshot_probe(sys,mg,info,states,vl,stencil,srg,ppg)
      allocate(result(0,0))
      return
    endif
    if(present(solver_override))then
      if(solver_override/='hybrid'.or.hsex_w_solver/='reduced'.or.present(extra_velocity)) &
        error stop 'invalid same-GS screening reference override'
      solver=solver_override
    endif
    if(all(num_kgrid==1).and.(solver=='pdep'.or.solver=='sternheimer'))then
      if(present(extra_velocity))error stop 'Gamma PDEP currently supports equilibrium screening only'
      call hsex_compute_gamma(sys,mg,info,states,eigen,pop,nv,boxes,vl,stencil,srg,ppg,first,nlocal,result)
      return
    endif
    n=mg%ie-mg%is+1
    nreal=product(n)
    ncomponent=size(states%zwf,4)
    ng=nreal*ncomponent
    dv=abs(dot_product(sys%primitive_a(:,1), &
       [sys%primitive_a(2,2)*sys%primitive_a(3,3)-sys%primitive_a(3,2)*sys%primitive_a(2,3), &
        sys%primitive_a(3,2)*sys%primitive_a(1,3)-sys%primitive_a(1,2)*sys%primitive_a(3,3), &
        sys%primitive_a(1,2)*sys%primitive_a(2,3)-sys%primitive_a(2,2)*sys%primitive_a(1,3)]))/nreal
    nkeep=sys%no
    if(ncomponent==2)then
      if(.not.screen_finite(pop).or..not.screen_finite(eigen))error stop 'nonfinite SOC screening reference'
      if(any(pop<0d0).or.any(pop>1d0))error stop 'invalid SOC screening occupation'
      ! Retain the occupied manifold and the first empty band. Higher empty
      ! states are handled by the projected Sternheimer solve.
      do band=sys%no,1,-1
        if(any(pop(band,:)>1d-12))exit
      enddo
      nkeep=min(sys%no,max(nv,band)+1)
      if(comm_is_root(info%id_ko))write(*,*)'HSEX_STERNHEIMER_BANDS retained/available',nkeep,sys%no
    endif
    allocate(screen_states(ng,nkeep,info%ik_s:info%ik_e),stat=astat)
    status=merge(1,0,astat/=0)
    call comm_get_max(status,info%icomm_ko)
    if(status/=0)error stop 'internal screening workspace allocation failed'
    do ik=info%ik_s,info%ik_e
      screen_states(:,:,ik)= &
        reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,:nkeep,ik,1),[ng,nkeep])*sqrt(dv)
    enddo
    if(hsex_reduced_target>0)then
      if(ncomponent/=1.or.extra_velocity_present())error stop 'fixed-k diagnostic requires scalar equilibrium Hamiltonian'
      if(maxval(abs(sys%wtk-1d0/sys%nk))>1d-12)error stop 'fixed-k diagnostic requires full uniform weights'
      if(maxval(abs(pop(:nv,:)-2d0))>1d-12.or.maxval(abs(pop(nv+1:,:)))>1d-12) &
        error stop 'fixed-k diagnostic requires integer insulating populations'
      call reduced_q_test(n,num_kgrid,hsex_w_qgrid,hsex_w_cutoff,sys%primitive_a,sys%vec_k(:,:sys%nk), &
        eigen,screen_states,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,hsex_reduced_target, &
        hsex_reduced_qcount,hsex_reduced_maxrank,hsex_reduced_tol,apply_grid,status, &
        yn_hsex_reduced_symmetry=='y',sys%Rion,sys%kion)
      call comm_get_max(status,info%icomm_ko)
      if(status/=0)error stop 'fixed-k response compression diagnostic rejected'
      allocate(result(0,0))
      return
    endif
    if(hsex_pdep_rank/=0)then
      if(ncomponent/=1.or.extra_velocity_present())error stop 'PDEP diagnostic requires scalar equilibrium Hamiltonian'
      if(maxval(abs(sys%wtk-1d0/sys%nk))>1d-12)error stop 'PDEP diagnostic requires full uniform k weights'
      if(maxval(abs(pop(:nv,:)-2d0))>1d-12.or.maxval(abs(pop(nv+1:,:)))>1d-12) &
        error stop 'PDEP diagnostic requires integer insulating populations'
      call pdep_q_test(n,num_kgrid,hsex_w_qgrid,hsex_pdep_qpoint,hsex_w_cutoff,sys%primitive_a, &
        sys%vec_k(:,:sys%nk),eigen,screen_states,nv,info%ik_s,info%ik_e,info%icomm_ko,info%id_ko, &
        hsex_pdep_rank,hsex_pdep_reference,apply_grid,status)
      call comm_get_max(status,info%icomm_ko)
      if(status/=0)error stop 'finite-q PDEP diagnostic rejected'
      allocate(result(0,0))
      return
    endif
    allocate(result(nlocal,sys%nk))
    threads=1
!$  threads=omp_get_max_threads()
    if(solver=='reduced')then
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
        sys%vec_k(:,:sys%nk),eigen(:nkeep,:),pop(:nkeep,:),screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
        info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
        apply_grid,temperature,distributed_states=.true.,apply_velocity=velocity_total, &
        reduced_rank=hsex_reduced_maxrank,reduced_batch=hsex_reduced_batch,reduced_tolerance=hsex_reduced_tol, &
        reduced_symmetry=yn_hsex_reduced_symmetry=='y',atoms=sys%Rion,species=sys%kion)
    else if(ncomponent==2.or.solver=='sternheimer')then
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
      sys%vec_k(:,:sys%nk),eigen(:nkeep,:),pop(:nkeep,:),screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
      info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
      apply_grid,temperature,yn_hsex_head_only=='y',distributed_states=.true.,apply_velocity=velocity_total)
    else if(solver=='hybrid')then
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
      sys%vec_k(:,:sys%nk),eigen,pop,screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
      info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
      distributed_states=.true.,apply_velocity=velocity_total,apply_optical_h=apply_grid)
    else
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
      sys%vec_k(:,:sys%nk),eigen,pop,screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
      info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
      distributed_states=.true.,apply_velocity=velocity_total)
    endif
    call comm_get_max(status,info%icomm_ko)
    if(status/=0)then
      if(comm_is_root(info%id_ko))write(*,*)'Internal screening failure code:',status
      error stop 'internal screening rejected W; previous W preserved'
    endif
  contains
    logical function extra_velocity_present()
      implicit none
      extra_velocity_present=present(extra_velocity)
    end function
    subroutine velocity_total(index,vector,action,ierr)
      use hsex_grid_velocity, only: grid_velocity
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:,:)
      integer,intent(out) :: ierr
      complex(8) :: extra(size(vector),3)
      call grid_velocity(sys,mg,stencil,ppg,index,vector,action)
      ierr=0
      if(present(extra_velocity))then
        call extra_velocity(index,vector,extra,ierr)
        if(ierr/=0)return
        action=action+extra
      endif
    end subroutine
    subroutine apply_grid(index,vector,action,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr

      ierr=1
      if(index<info%ik_s.or.index>info%ik_e)return
      target_info=info
      target_info%ik_s=index
      target_info%ik_e=index
      target_info%numk=1
      target_info%io_s=1
      target_info%io_e=1
      target_info%numo=1
      x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,1,index,1)= &
        reshape(vector,[n,ncomponent])
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,target_info,mg,vl,sys,stencil,srg,ppg)
      action=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),:,1,index,1),[ng])
      ierr=0
    end subroutine apply_grid

    subroutine load_h(index,h,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),intent(out) :: h(:,:)
      integer,intent(out) :: ierr
      integer :: g,j,col,count,component
      ierr=0
      if(index<info%ik_s.or.index>info%ik_e)then
        ierr=1
        return
      endif
      ! No spatial/orbital decomposition is allowed on this path. The grid
      ! Hamiltonian and the already staged HSEX action are k-local.
      target_info=info
      target_info%ik_s=index;target_info%ik_e=index;target_info%numk=1
      do g=1,ng,sys%no
        x%zwf=0d0
        count=min(sys%no,ng-g+1)
        do j=1,count
          col=mod(g+j-2,nreal)
          component=(g+j-2)/nreal+1
          x%zwf(mg%is(1)+mod(col,n(1)),mg%is(2)+mod(col/n(1),n(2)), &
            mg%is(3)+col/(n(1)*n(2)),component,j,index,1)=1d0
        enddo
        x%update_zwf_overlap=.false.
        call hpsi(x,hx,target_info,mg,vl,sys,stencil,srg,ppg)
        do j=1,count
          h(:,g+j-1)=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
            mg%is(3):mg%ie(3),:,j,index,1),[ng])
        enddo
      enddo
    end subroutine load_h
  end subroutine hsex_compute_internal

  subroutine hsex_compute_gamma(sys,mg,info,states,eigen,pop,nv,boxes,vl,stencil,srg,ppg,first,nlocal,result)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar,s_stencil,s_pp_grid, &
      allocate_orbital_complex,deallocate_orbital
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi
    use hsex_grid_velocity, only: grid_velocity
    use hsex_gamma_pdep, only: gamma_pdep_screen
    use hsex_screen_math, only: screen_finite
    use communication, only: comm_summation,comm_get_max,comm_is_root
    use salmon_global, only: hsex_w_solver,hsex_pdep_rank,hsex_w_cutoff,hsex_pdep_tol,hsex_pdep_maxiter, &
      hsex_pdep_basis_read,hsex_pdep_basis_write
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    real(8),intent(in) :: eigen(:,:),pop(:,:)
    integer,intent(in) :: nv,boxes,first,nlocal
    type(s_scalar),intent(in) :: vl(:)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    real(8),allocatable,intent(out) :: result(:,:)
    type(s_parallel_info) :: work_info
    type(s_orbital) :: work_x,work_hx
    complex(8),allocatable :: local(:,:),occupied(:,:),kernel(:,:)
    integer,allocatable :: hits(:),all_hits(:)
    real(8) :: macro(3,3)
    integer :: n(3),ng,band,a,b,local_code,requested,first_band,last_band,i
    n=mg%ie-mg%is+1;ng=product(n);local_code=0
    if(sys%nk/=1.or.info%npk/=1.or.info%if_divide_rspace.or.sys%nspin/=1)local_code=1
    if(maxval(abs(sys%vec_k(:,1)))>1d-12.or.maxval(abs(sys%vec_Ac))>1d-12)local_code=1
    if(size(eigen,1)<nv+1.or..not.screen_finite(eigen).or..not.screen_finite(pop))local_code=1
    call comm_get_max(local_code,info%icomm_ko)
    if(local_code/=0)error stop 'Gamma screening requires full-grid insulating scalar Gamma GS with a gap probe'
    if(maxval(abs(pop(:nv,1)-2d0))>1d-12.or.maxval(abs(pop(nv+1:,1)))>1d-12) &
      error stop 'Gamma screening requires integer populations'
    if(eigen(nv+1,1)-eigen(nv,1)<=1d-8)error stop 'Gamma screening requires a positive GS gap'
    allocate(local(ng,nv),occupied(ng,nv),hits(nv),all_hits(nv),stat=local_code)
    call comm_get_max(local_code,info%icomm_ko)
    if(local_code/=0)error stop 'Gamma occupied-projector workspace allocation failed'
    local=0d0;hits=0
    do band=max(1,info%io_s),min(nv,info%io_e)
      if(allocated(states%rwf))then
        local(:,band)=cmplx(reshape(states%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
          mg%is(3):mg%ie(3),1,band,1,1),[ng]),0d0,8)*sqrt(sys%hvol)
      else
        local(:,band)=reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
          mg%is(3):mg%ie(3),1,band,1,1),[ng])*sqrt(sys%hvol)
      endif
      hits(band)=1
    enddo
    call comm_summation(hits,all_hits,nv,info%icomm_ko)
    if(any(all_hits/=1))error stop 'Gamma occupied band ownership is not unique'
    call comm_summation(local,occupied,size(occupied),info%icomm_ko);deallocate(local)
    work_info=info;work_info%io_s=1;work_info%io_e=1;work_info%numo=1;work_info%if_divide_orbit=.false.
    call allocate_orbital_complex(1,mg,work_info,work_x)
    call allocate_orbital_complex(1,mg,work_info,work_hx)
    allocate(kernel(boxes**3,boxes**3),result(nlocal,1))
    first_band=1+nv*info%id_ko/info%isize_ko;last_band=nv*(info%id_ko+1)/info%isize_ko
    requested=-1
    if(hsex_w_solver=='pdep')requested=hsex_pdep_rank
    call gamma_pdep_screen(n,boxes,hsex_w_cutoff,sys%primitive_a,occupied,eigen(:nv,1),first_band,last_band, &
      info%icomm_ko,info%id_ko,requested,hsex_pdep_tol,hsex_pdep_maxiter,apply,velocity,macro,kernel,local_code, &
      hsex_pdep_basis_read,hsex_pdep_basis_write)
    call deallocate_orbital(work_x);call deallocate_orbital(work_hx)
    call comm_get_max(local_code,info%icomm_ko)
    if(local_code/=0)error stop 'Gamma PDEP screening rejected; previous W preserved'
    do i=1,nlocal
      a=1+mod(first+i-2,boxes**3);b=1+(first+i-2)/(boxes**3)
      result(i,1)=real(kernel(a,b),8)
    enddo
    if(comm_is_root(info%id_ko))write(*,'(a)')'HSEX_GAMMA_BAND_MPI_PASS occupied-only response'
  contains
    subroutine apply(x,y,code)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: code
      work_x%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,1,1)=reshape(x,n)
      work_x%update_zwf_overlap=.false.
      call hpsi(work_x,work_hx,work_info,mg,vl,sys,stencil,srg,ppg)
      y=reshape(work_hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,1,1,1),[ng])
      code=0
    end subroutine
    subroutine velocity(x,y,code)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: code
      call grid_velocity(sys,mg,stencil,ppg,1,x,y);code=0
    end subroutine
  end subroutine hsex_compute_gamma

  subroutine hsex_save_internal(path,magic,step,sys,mg,info,boxes,nv,vref,eref,eigen,pop,first,values)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info
    use communication, only: comm_summation,comm_bcast,comm_is_root
    use salmon_global, only: num_kgrid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
    implicit none
    character(*),intent(in) :: path
    character(4),intent(in) :: magic
    integer,intent(in) :: step,boxes,nv,first
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    real(8),intent(in) :: vref(:),eref,eigen(:,:),pop(:,:),values(:,:)
    real(8),allocatable :: partial(:),block(:),partner(:),r(:,:)
    integer :: u,ios,ik,j,ip,na,l(3),minus,last
    logical :: exists
    interface
      function link_file(old,new) bind(C,name='link') result(code)
        import c_int,c_char
        implicit none
        character(c_char),intent(in) :: old(*),new(*)
        integer(c_int) :: code
      end function link_file
      function unlink_file(pathname) bind(C,name='unlink') result(code)
        import c_int,c_char
        implicit none
        character(c_char),intent(in) :: pathname(*)
        integer(c_int) :: code
      end function unlink_file
    end interface
    na=boxes**3;last=first+size(values,1)-1
    allocate(partial(na*na),block(na*na),partner(na*na),r(3,sys%nk))
    do ik=1,sys%nk
      l=[mod(ik-1,num_kgrid(1)),mod((ik-1)/num_kgrid(1),num_kgrid(2)), &
         (ik-1)/(num_kgrid(1)*num_kgrid(2))]
      where(2*l>=num_kgrid)l=l-num_kgrid
      r(:,ik)=matmul(sys%primitive_a,real(l,8))
    enddo
    ! Validate every block collectively before creating a usable checkpoint.
    do ik=1,sys%nk
      partial=0d0;partial(first:last)=values(:,ik)
      call comm_summation(partial,block,size(block),info%icomm_ko)
      l=modulo(-[mod(ik-1,num_kgrid(1)),mod((ik-1)/num_kgrid(1),num_kgrid(2)), &
        (ik-1)/(num_kgrid(1)*num_kgrid(2))],num_kgrid)
      minus=1+l(1)+num_kgrid(1)*(l(2)+num_kgrid(2)*l(3))
      partial=0d0;partial(first:last)=values(:,minus)
      call comm_summation(partial,partner,size(partner),info%icomm_ko)
      do j=1,na*na
        ip=1+(j-1)/na+na*mod(j-1,na)
        if(.not.ieee_is_finite(block(j)))error stop 'nonfinite internal W checkpoint'
        if(.not.(abs(block(j)-partner(ip))<=1d-8))error stop 'internal W reciprocity failure'
      enddo
    enddo
    ios=0
    if(comm_is_root(info%id_ko))then
      inquire(file=path,exist=exists)
      if(exists)then
        ios=1
      else
        open(newunit=u,file=trim(path)//'.partial',status='new',access='stream',form='unformatted',iostat=ios)
        if(ios==0)write(u,iostat=ios)magic,mg%ie-mg%is+1,num_kgrid,boxes,nv,sys%nion, &
          sys%primitive_a,sys%vec_k(:,:sys%nk),sys%Rion,vref,eref,r
      endif
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot create internal W checkpoint'
    do ik=1,sys%nk
      partial=0d0;partial(first:last)=values(:,ik)
      call comm_summation(partial,block,size(block),info%icomm_ko,0)
      if(comm_is_root(info%id_ko))write(u,iostat=ios)block
      call comm_bcast(ios,info%icomm_ko,0)
      if(ios/=0)error stop 'cannot write internal W checkpoint'
    enddo
    if(comm_is_root(info%id_ko))then
      if(magic=='HXP1')write(u,iostat=ios)step,sys%no,eigen,pop
      if(ios==0)close(u,iostat=ios)
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot finish internal W checkpoint'
    if(comm_is_root(info%id_ko))then
      ios=link_file(trim(path)//'.partial'//c_null_char,trim(path)//c_null_char)
      if(ios==0)ios=unlink_file(trim(path)//'.partial'//c_null_char)
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot publish internal W checkpoint'
  end subroutine hsex_save_internal
end module hsex_internal
