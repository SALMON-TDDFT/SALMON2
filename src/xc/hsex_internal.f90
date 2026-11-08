#include "config.h"
! SALMON adapter: keep the Hamiltonian and projection manifold in memory.
module hsex_internal
  implicit none
  private
  public :: hsex_compute_internal,hsex_save_internal
contains
  subroutine hsex_compute_internal(sys,mg,info,states,x,hx,eigen,pop,nv,boxes,eta,vl,stencil,srg,ppg, &
                                   first,nlocal,result)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar,s_stencil,s_pp_grid
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi
    use communication, only: comm_get_max,comm_is_root
    use salmon_global, only: num_kgrid,hsex_w_qgrid,hsex_w_cutoff,hsex_w_head_tol,hsex_w_eta_tol,temperature,yn_hsex_head_only
    use hsex_screen_driver, only: screen_compute
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
    real(8),allocatable,intent(out) :: result(:,:)
    complex(8),allocatable :: screen_states(:,:,:)
    type(s_parallel_info) :: target_info
    integer :: n(3),ng,nreal,ncomponent,ik,status,threads,astat,nkeep,band
    real(8) :: dv
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
    allocate(result(nlocal,sys%nk))
    threads=1
!$  threads=omp_get_max_threads()
    if(ncomponent==2)then
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
      sys%vec_k(:,:sys%nk),eigen(:nkeep,:),pop(:nkeep,:),screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
      info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
      apply_spinor,temperature,yn_hsex_head_only=='y',distributed_states=.true.)
    else
      call screen_compute(n,num_kgrid,hsex_w_qgrid,boxes,hsex_w_cutoff,sys%primitive_a, &
      sys%vec_k(:,:sys%nk),eigen,pop,screen_states,nv,eta,hsex_w_head_tol,hsex_w_eta_tol,threads, &
      info%ik_s,info%ik_e,info%icomm_ko,info%id_ko,info%isize_ko,load_h,first,nlocal,result,status, &
      distributed_states=.true.)
    endif
    call comm_get_max(status,info%icomm_ko)
    if(status/=0)then
      if(comm_is_root(info%id_ko))write(*,*)'Internal screening failure code:',status
      error stop 'internal screening rejected W; previous W preserved'
    endif
  contains
    subroutine apply_spinor(index,vector,action,ierr)
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
    end subroutine apply_spinor

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
