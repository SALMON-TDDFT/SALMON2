#include "config.h"
! Static microscopic W, projected onto disjoint normalized real-space boxes.
! Only the exchange action is projected. SALMON propagates full grid orbitals.
module kli_hsex
  implicit none
  private
  public :: hsex_init,hsex_stage,hsex_add,hsex_current,hsex_xc_energy,hsex_induced_energy
  public :: hsex_snapshot,hsex_replace_w,hsex_current_decomposed,hsex_update_internal,hsex_checkpoint,hsex_restore
  ! Deliberately shared Hamiltonian framework state for the active SALMON system.
  ! Initialization, propagation, current evaluation and checkpoint hooks share
  ! this private state; callers access it only through the public procedures.
  real(8),allocatable :: base_vref(:)
  real(8) :: base_eref
  integer :: pair_start
  real(8),allocatable :: translation(:,:),w(:,:),position(:,:),local_xyz(:,:),vec_k(:,:),current(:)
  complex(8),allocatable :: reference(:,:),incoming(:,:),sigma(:,:,:),projector(:,:),offset(:),exchange_kernel(:,:)
  integer,allocatable :: box(:)
  integer :: na,ng,nk
  real(8) :: xc0,exchange_energy=0d0,projector_ac(3),endpoint_ac(3)
  logical :: have_projector=.false.,have_endpoint=.false.
  interface
    subroutine zgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
      implicit none
      character(1),intent(in) :: transa,transb
      integer,intent(in) :: m,n,k,lda,ldb,ldc
      complex(8),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
      complex(8),intent(inout) :: c(ldc,*)
    end subroutine
  end interface
contains
  integer function hsex_k_threads(count)
!$  use omp_lib, only: omp_get_max_threads
    implicit none
    integer,intent(in) :: count
    hsex_k_threads=1
!$  hsex_k_threads=max(1,min(count,omp_get_max_threads()))
  end function hsex_k_threads

  logical function hsex_parallel_blas()
#ifdef HAVE_EXX_OPENBLAS_THREADS
    use iso_c_binding, only: c_int
#endif
    implicit none
#ifdef HAVE_EXX_OPENBLAS_THREADS
    interface
      function get_parallel() bind(C,name='openblas_get_parallel') result(mode)
        import c_int
        implicit none
        integer(c_int) :: mode
      end function
    end interface
    ! A sequential OpenBLAS build may lack concurrent-call locking.
    hsex_parallel_blas=get_parallel()/=0
#ifdef HAVE_HSEX_APPLE_GNU
    ! The observed shutdown failure is inside this OpenMP OpenBLAS backend.
    ! Keep calls on the parent thread; MPI and independent FFT OMP remain active.
    if(get_parallel()==2)hsex_parallel_blas=.false.
#endif
#else
    hsex_parallel_blas=.true.
#endif
  end function hsex_parallel_blas

  subroutine hsex_init(sys,mg,info,orb,exc,exchange,tbxc,restarting)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
    use hsex_k_fft, only: hsex_fft_init
    use kli_scf, only: kli_check_reference
    use salmon_global, only: num_kgrid,cval
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    real(8),intent(in)::exc
    logical,intent(in)::exchange
    type(s_scalar),optional,intent(in)::tbxc
    logical,optional,intent(in)::restarting
    logical :: restart_state
    integer::u,dims(9),b,n(3),g,j,index3(3),l(3),ik,i,first_pair,nlocal
    integer(8) :: data_start
    character(4)::magic
    real(8)::cell(3,3),k(3,sys%nk),ions(3,sys%nion),e
    real(8),allocatable::v(:),centre(:,:)
    xc0=exc
    restart_state=.false.
    if(present(restarting))restart_state=restarting
    if(.not.exchange)return
    n=mg%ie-mg%is+1;ng=product(n);nk=sys%nk
    open(newunit=u,file='hsex-w.bin',access='stream',form='unformatted',status='old',action='read')
    read(u)magic,dims
    b=dims(7)
    if(magic/=merge('HXT1','HXR1',present(tbxc)).or.any(dims(1:3)/=n).or.any(dims(4:6)/=num_kgrid).or.product(num_kgrid)/=nk.or. &
       dims(8)/=sys%no.or.dims(9)/=sys%nion.or.b<1.or.b>minval(n))error stop 'incompatible HSEX screening grid'
    if(any(mod(n,b)/=0))error stop 'HSEX boxes must divide real-space grid'
    na=b**3
    allocate(v(ng),translation(3,nk),local_xyz(ng,3),position(3,na),centre(3,na),box(ng))
    read(u)cell,k,ions,v,e,translation
    inquire(unit=u,pos=data_start)
    if(any(cell/=sys%primitive_a).or.any(k/=sys%vec_k(:,:nk)).or.any(ions/=sys%Rion))error stop 'HSEX crystal mismatch'
    if(present(tbxc))then
      if(.not.restart_state)then
        if(e/=cval.or..not.(maxval(abs(v-reshape(tbxc%f,[ng])))<1d-5))error stop 'W mismatch with TB-mBJ GS potential'
      endif
    else
      call kli_check_reference(v,e)
    endif
    base_vref=v;base_eref=e
    call hsex_fft_init(na*na,dims(4:6),cell,k,info,first_pair,nlocal)
    do ik=1,nk
      l=[mod(ik-1,dims(4)),mod((ik-1)/dims(4),dims(5)),(ik-1)/(dims(4)*dims(5))]
      where(2*l>=dims(4:6))l=l-dims(4:6)
      if(maxval(abs(translation(:,ik)-matmul(cell,real(l,8))))>1d-9)error stop 'invalid HSEX translation order'
    enddo
    pair_start=first_pair
    allocate(w(nlocal,nk))
    do ik=1,nk
      if(nlocal>0)read(u,pos=data_start+8_8*(int(na,8)*na*(ik-1)+first_pair-1))w(:,ik)
      do i=1,nlocal
        if(.not.ieee_is_finite(w(i,ik)))error stop 'nonfinite W'
      enddo
    enddo
    close(u)
    centre=0
    do g=1,ng
      index3=[mod(g-1,n(1)),mod((g-1)/n(1),n(2)),(g-1)/(n(1)*n(2))]
      l=index3/(n/b);j=1+l(1)+b*(l(2)+b*l(3))
      box(g)=j;local_xyz(g,:)=matmul(cell,real(index3,8)/n)
      centre(:,j)=centre(:,j)+local_xyz(g,:)
    enddo
    position=centre/(ng/na)
    do g=1,ng
      local_xyz(g,:)=local_xyz(g,:)-position(:,box(g))
    enddo
    allocate(reference(nlocal,nk),incoming(nlocal,nk),offset(nk), &
             sigma(na,na,info%ik_s:info%ik_e),projector(ng,info%ik_s:info%ik_e),current(3))
    vec_k=k
    reference=0
    call hsex_stage(sys,mg,info,orb,[0d0,0d0,0d0],3)
  end subroutine

  subroutine project_grid(x,f,p)
    implicit none
    complex(8),intent(in)::x(:,:),f(:)
    complex(8),intent(out)::p(:,:)
    integer::j,g
    p=0
    do j=1,size(x,2)
      do g=1,ng
        p(box(g),j)=p(box(g),j)+f(g)*x(g,j)
      enddo
    enddo
  end subroutine

  subroutine hsex_stage(sys,mg,info,orb,ac,stage)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
    use hsex_k_fft, only: hsex_forward,hsex_backward
    use communication, only: comm_summation
    use exx_blas_threads, only: exx_blas_thread_control
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    real(8),intent(in)::ac(3)
    integer,intent(in)::stage ! 0 incoming, 1 predictor, 2 endpoint, 3 reference
    integer::ik,ib,l,axis,g,first,last,blas_parent(3),blas_worker(3)
    complex(8)::p(na,sys%no),hp(na,sys%no),d(na,na),z
    complex(8),allocatable::x(:,:),pk(:,:),density(:,:),act(:,:),projected(:,:,:)
    real(8)::kac(3),vol,reduced(4),local(4)
    logical::update_phase
    if(stage==0.and.have_endpoint)then
      if(all(ac==endpoint_ac))return ! endpoint stored its density and operator for the next predictor
    endif
    update_phase=.not.have_projector
    if(have_projector)update_phase=any(ac/=projector_ac)
    first=info%ik_s;last=info%ik_e
    allocate(x(ng,sys%no),pk(na*na,first:last),projected(na,sys%no,first:last))
    if(update_phase)offset=exp(cmplx(0d0,matmul(transpose(translation),vec_k(:,1)+ac),8))
    blas_parent=0
#ifdef HAVE_EXX_OPENBLAS_THREADS
    call exx_blas_thread_control(1,blas_parent)
#endif
    !$omp parallel if(hsex_parallel_blas()) num_threads(hsex_k_threads(info%numk)) &
    !$omp& default(shared) private(ik,kac,g,x,p,d,blas_worker)
    blas_worker=0
    call exx_blas_thread_control(1,blas_worker)
    !$omp do schedule(static)
    do ik=first,last
      if(update_phase)then
        kac=vec_k(:,ik)+ac
        do g=1,ng
          projector(g,ik)=sqrt(sys%hvol/(ng/na))*exp(cmplx(0d0,dot_product(local_xyz(g,:)+position(:,box(g)),kac),8))
        enddo
      endif
      x=reshape(orb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,:,ik,1),[ng,sys%no])
      call project_grid(x,projector(:,ik),p)
      projected(:,:,ik)=p
      call zgemm('N','C',na,na,sys%no,(1d0,0d0),p,na,p,na,(0d0,0d0),d,na)
      pk(:,ik)=reshape(d,[na*na])
    enddo
    !$omp end do
    call exx_blas_thread_control(0,blas_worker)
    !$omp end parallel
    call exx_blas_thread_control(0,blas_parent)
    projector_ac=ac
    have_projector=.true.
    allocate(density(size(w,1),nk),act(na*na,first:last))
    call hsex_forward(pk,offset,info%icomm_ko,density)
    if(stage==3)reference=density
    density=density-reference
    if(stage==0.or.stage==2)incoming=density
    if(stage==2)then
      endpoint_ac=ac
      have_endpoint=.true.
    endif
    if(stage==1)density=.5d0*(incoming+density)
    ! Spin degeneracy 2 cancels 1/2 in the energy. Action uses one spin projector.
    exchange_kernel=-w*density
    call hsex_backward(exchange_kernel,offset,info%icomm_ko,act)
    do ik=first,last
      d=reshape(act(:,ik),[na,na])
      sigma(:,:,ik)=.5d0*(d+conjg(transpose(d)))
    enddo
    if(stage/=2)return
    vol=ng*sys%hvol
    local=0d0
    local(4)=sum(real(conjg(density)*exchange_kernel,8))
    ! Differentiate the lattice phase in R space. The k sum is precisely
    ! conjg(density+reference), already available on this pair owner.
    do l=1,nk
      z=sum(exchange_kernel(:,l)*conjg(density(:,l)+reference(:,l)))
      local(1:3)=local(1:3)-2*aimag(z)*translation(:,l)/vol
    enddo
    blas_parent=0
#ifdef HAVE_EXX_OPENBLAS_THREADS
    call exx_blas_thread_control(1,blas_parent)
#endif
    !$omp parallel if(hsex_parallel_blas()) num_threads(hsex_k_threads(info%numk)) &
    !$omp& default(shared) private(ik,hp,x,g,z,blas_worker) reduction(+:local)
    blas_worker=0
    call exx_blas_thread_control(1,blas_worker)
    !$omp do schedule(static)
    do ik=first,last
      call zgemm('N','N',na,sys%no,na,(1d0,0d0),sigma(:,:,ik),na, &
        projected(:,:,ik),na,(0d0,0d0),hp,na)
      x=reshape(orb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,:,ik,1),[ng,sys%no])
      do g=1,ng
        z=sum(conjg(x(g,:)*projector(g,ik))*hp(box(g),:))
        local(1:3)=local(1:3)+4*aimag(z)*(local_xyz(g,:)+position(:,box(g)))/(nk*vol)
      enddo
    enddo
    !$omp end do
    call exx_blas_thread_control(0,blas_worker)
    !$omp end parallel
    call exx_blas_thread_control(0,blas_parent)
    call comm_summation(local,reduced,4,info%icomm_ko)
    current=reduced(1:3);exchange_energy=reduced(4)
  end subroutine

  subroutine hsex_add(sys,mg,info,orb,hpsi)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
    use exx_blas_threads, only: exx_blas_thread_control
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    type(s_orbital),intent(inout)::hpsi
    complex(8)::x(ng,sys%no),p(na,sys%no),hp(na,sys%no)
    integer::ik,ib,n(3),blas_parent(3),blas_worker(3)
    n=mg%ie-mg%is+1
    blas_parent=0
#ifdef HAVE_EXX_OPENBLAS_THREADS
    call exx_blas_thread_control(1,blas_parent)
#endif
    !$omp parallel if(hsex_parallel_blas()) num_threads(hsex_k_threads(info%numk)) &
    !$omp& default(shared) private(ik,ib,x,p,hp,blas_worker)
    blas_worker=0
    call exx_blas_thread_control(1,blas_worker)
    !$omp do schedule(static)
    do ik=info%ik_s,info%ik_e
      x=reshape(orb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,:,ik,1),[ng,sys%no])
      call project_grid(x,projector(:,ik),p)
      call zgemm('N','N',na,sys%no,na,(1d0,0d0),sigma(:,:,ik),na,p,na,(0d0,0d0),hp,na)
      do ib=1,sys%no
        associate(h=>hpsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,ib,ik,1))
          h=h+reshape(conjg(projector(:,ik))*hp(box,ib)/sys%hvol,n)
        end associate
      enddo
    enddo
    !$omp end do
    call exx_blas_thread_control(0,blas_worker)
    !$omp end parallel
    call exx_blas_thread_control(0,blas_parent)
  end subroutine

  subroutine hsex_snapshot(step,sys,mg,info,states,x,hx,eigen,weighted_occ,vl,stencil,srg,ppg)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar,s_stencil,s_pp_grid
    use hamiltonian, only: hpsi
    use sendrecv_grid, only: s_sendrecv_grid
    implicit none
    integer,intent(in)::step
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::states
    type(s_orbital),intent(inout)::x,hx
    real(8),intent(in)::eigen(:,:),weighted_occ(:,:)
    type(s_scalar),intent(in)::vl(:)
    type(s_stencil),intent(in)::stencil
    type(s_sendrecv_grid),intent(inout)::srg
    type(s_pp_grid),intent(in)::ppg
    complex(8),allocatable :: column(:)
    integer::ik,unit(info%ik_s:info%ik_e),n(3),g,j,col,count
    character(64)::name
    n=mg%ie-mg%is+1
    do ik=info%ik_s,info%ik_e
      write(name,'("hsex-pump-",i6.6,".bin")')ik
      open(newunit=unit(ik),file=name,status='new',access='stream',form='unformatted')
      write(unit(ik))'HSP2',n,sys%nk,ik,sys%no,step,sys%primitive_a,sys%vec_k(:,ik),sys%vec_Ac, &
        eigen(:,ik),weighted_occ(:,ik)/sys%wtk(ik), &
        states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,:,ik,1)
    enddo
    allocate(column(ng))
    do g=1,ng,sys%no
      x%zwf=0
      count=min(sys%no,ng-g+1)
      do j=1,count
        col=g+j-2
        x%zwf(mg%is(1)+mod(col,n(1)),mg%is(2)+mod(col/n(1),n(2)),mg%is(3)+col/(n(1)*n(2)),1,j,:,1)=1
      enddo
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,info,mg,vl,sys,stencil,srg,ppg)
      do ik=info%ik_s,info%ik_e
        do j=1,count
          column=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,j,ik,1),[ng])
          write(unit(ik))column(g+j-1:ng)
        enddo
      enddo
    enddo
    do ik=info%ik_s,info%ik_e
      close(unit(ik))
    enddo
  end subroutine

  subroutine hsex_update_internal(step,sys,mg,info,states,x,hx,eigen,weighted_occ,vl,stencil,srg,ppg)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar,s_stencil,s_pp_grid
    use sendrecv_grid, only: s_sendrecv_grid
    use hsex_internal, only: hsex_compute_internal,hsex_save_internal
    use salmon_global, only: hsex_w_eta,nelec
    use communication, only: comm_bcast,comm_is_root
    implicit none
    integer,intent(in) :: step
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    type(s_orbital),intent(inout) :: x,hx
    real(8),intent(in) :: eigen(:,:),weighted_occ(:,:)
    type(s_scalar),intent(in) :: vl(:)
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_pp_grid),intent(in) :: ppg
    real(8),allocatable :: next_w(:,:),pop(:,:),vref(:)
    real(8) :: cell(3,3),k(3,nk),ions(3,sys%nion),eref
    integer :: ik,u,ios,dims(9),boxes
    character(4) :: magic
    character(64) :: path
    allocate(pop(sys%no,nk),vref(ng))
    do ik=1,nk
      pop(:,ik)=weighted_occ(:,ik)/sys%wtk(ik)
    enddo
    boxes=nint(real(na,8)**(1d0/3d0))
    call hsex_compute_internal(sys,mg,info,states,x,hx,eigen,pop,nelec/2,boxes,hsex_w_eta, &
      vl,stencil,srg,ppg,pair_start,size(w,1),next_w)
    ios=0
    if(comm_is_root(info%id_ko))then
      open(newunit=u,file='hsex-w.bin',status='old',access='stream',form='unformatted',iostat=ios)
      if(ios==0)then
        read(u,iostat=ios)magic,dims,cell,k,ions,vref,eref
        close(u)
      endif
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot read W reference metadata'
    write(path,'("hsex-updated-",i9.9,".bin")')step
    call hsex_save_internal(trim(path),'HXP1',step,sys,mg,info,boxes,nelec/2,vref,eref, &
      eigen,pop,pair_start,next_w)
    call move_alloc(next_w,w)
    have_endpoint=.false.
    if(comm_is_root(info%id_ko))write(*,*)'HSEX_INTERNAL_W_UPDATE_PASS step=',step,' projection states=',sys%no
  end subroutine hsex_update_internal

  subroutine hsex_replace_w(path,step,sys,mg,info,eigen,weighted_occ)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info
    use salmon_global, only: num_kgrid
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in)::path
    integer,intent(in)::step
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    real(8),intent(in)::eigen(:,:),weighted_occ(:,:)
    integer::u,dims(9),ik,i,tag,nb
    character(4)::magic
    real(8)::cell(3,3),k(3,nk),ions(3,sys%nion),v(ng),e,r(3,nk),column(na*na)
    real(8)::values(sys%no,nk)
    open(newunit=u,file=path,status='old',access='stream',form='unformatted')
    read(u)magic,dims
    if(magic/='HXP1'.or.any(dims(1:3)/=mg%ie-mg%is+1).or.any(dims(4:6)/=num_kgrid).or. &
       dims(7)**3/=na.or.dims(9)/=sys%nion)error stop 'incompatible pump W'
    read(u)cell,k,ions,v,e,r
    if(any(cell/=sys%primitive_a).or.any(k/=vec_k).or.any(ions/=sys%Rion).or.any(r/=translation))error stop 'pump W mesh'
    do ik=1,nk
      read(u)column
      do i=1,size(column)
        if(.not.ieee_is_finite(column(i)))error stop 'nonfinite pump W'
      enddo
      w(:,ik)=column(pair_start:pair_start+size(w,1)-1)
    enddo
    read(u)tag,nb
    if(tag/=step.or.nb/=sys%no)error stop 'pump snapshot mismatch'
    read(u)values
    if(.not.(maxval(abs(values-eigen))<1d-7))error stop 'pump W eigenvalues mismatch'
    read(u)values
    do ik=1,nk
      if(.not.(maxval(abs(values(:,ik)-weighted_occ(:,ik)/sys%wtk(ik)))<1d-6))error stop 'pump W occupations mismatch'
    enddo
    close(u)
    have_endpoint=.false.
  end subroutine

  subroutine hsex_checkpoint(dir,step,sys,info)
    use structures, only: s_dft_system,s_parallel_info
    use communication, only: comm_summation,comm_bcast,comm_is_root
    implicit none
    character(*),intent(in) :: dir
    integer,intent(in) :: step
    type(s_dft_system),intent(in) :: sys
    type(s_parallel_info),intent(in) :: info
    real(8),allocatable :: piece(:),whole(:)
    complex(8),allocatable :: cpiece(:),cwhole(:)
    integer :: u,ios,ik,last
    if(.not.allocated(w))return
    allocate(piece(na*na),whole(na*na),cpiece(na*na),cwhole(na*na))
    last=pair_start+size(w,1)-1;ios=0
    if(comm_is_root(info%id_ko))then
      open(newunit=u,file=trim(dir)//'hsex-active.bin',status='replace',access='stream',form='unformatted',iostat=ios)
      if(ios==0)write(u,iostat=ios)'HXC2',step,na,nk,sys%primitive_a,vec_k,translation,xc0,base_eref,base_vref
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot create active W checkpoint'
    do ik=1,nk
      piece=0d0;piece(pair_start:last)=w(:,ik)
      cpiece=0d0;cpiece(pair_start:last)=reference(:,ik)
      call comm_summation(piece,whole,na*na,info%icomm_ko,0)
      call comm_summation(cpiece,cwhole,na*na,info%icomm_ko,0)
      if(comm_is_root(info%id_ko))write(u,iostat=ios)whole,cwhole
      call comm_bcast(ios,info%icomm_ko,0)
      if(ios/=0)error stop 'cannot write active W checkpoint'
    enddo
    if(comm_is_root(info%id_ko))close(u,iostat=ios)
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'cannot finish active W checkpoint'
  end subroutine hsex_checkpoint

  subroutine hsex_restore(dir,step,sys,mg,info,orb)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
    use communication, only: comm_bcast,comm_is_root
    use salmon_global, only: yn_hsex_update_w,hsex_screen_step,hsex_probe_w
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in) :: dir
    integer,intent(in) :: step
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: orb
    real(8),allocatable :: whole(:),new_w(:,:)
    complex(8),allocatable :: cwhole(:),new_ref(:,:)
    real(8) :: cell(3,3),k(3,nk),r(3,nk),saved_xc,saved_e,saved_v(ng)
    integer :: u,ios,ik,j,last,dims(3)
    character(4) :: magic
    logical :: exists
    if(.not.allocated(w))return
    ios=0;exists=.false.
    if(comm_is_root(info%id_ko))inquire(file=trim(dir)//'hsex-active.bin',exist=exists)
    call comm_bcast(exists,info%icomm_ko,0)
    if(.not.exists)then
      error stop 'HSEX restart missing active W and reference-density state'
    endif
    allocate(whole(na*na),cwhole(na*na),new_w(size(w,1),nk),new_ref(size(w,1),nk))
    if(comm_is_root(info%id_ko))then
      open(newunit=u,file=trim(dir)//'hsex-active.bin',status='old',access='stream',form='unformatted',iostat=ios)
      if(ios==0)then
        read(u,iostat=ios)magic,dims,cell,k,r,saved_xc,saved_e,saved_v
        if(ios==0)then
          if(magic/='HXC2'.or.any(dims/=[step,na,nk]))ios=1
          if(saved_e/=base_eref.or.any(saved_v/=base_vref))ios=1
          if(any(cell/=sys%primitive_a).or.any(k/=vec_k).or.any(r/=translation))ios=1
          if(.not.ieee_is_finite(saved_xc))ios=1
        endif
      endif
    endif
    call comm_bcast(ios,info%icomm_ko,0)
    if(ios/=0)error stop 'invalid active W checkpoint metadata'
    call comm_bcast(saved_xc,info%icomm_ko,0)
    last=pair_start+size(w,1)-1
    do ik=1,nk
      if(comm_is_root(info%id_ko))read(u,iostat=ios)whole,cwhole
      call comm_bcast(ios,info%icomm_ko,0)
      if(ios/=0)error stop 'truncated active W checkpoint'
      call comm_bcast(whole,info%icomm_ko,0);call comm_bcast(cwhole,info%icomm_ko,0)
      do j=1,na*na
        if(.not.ieee_is_finite(whole(j)).or..not.ieee_is_finite(real(cwhole(j),8)).or. &
           .not.ieee_is_finite(aimag(cwhole(j))))error stop 'nonfinite active W checkpoint'
      enddo
      new_w(:,ik)=whole(pair_start:last);new_ref(:,ik)=cwhole(pair_start:last)
    enddo
    if(comm_is_root(info%id_ko))close(u)
    call move_alloc(new_w,w);call move_alloc(new_ref,reference);xc0=saved_xc
    have_endpoint=.false.
    call hsex_stage(sys,mg,info,orb,sys%vec_Ac,2)
    if(comm_is_root(info%id_ko))write(*,*)'HSEX_ACTIVE_W_RESTORED step=',step
  end subroutine hsex_restore

  subroutine hsex_current(j)
    implicit none
    real(8),intent(inout)::j(3)
    j=j+current
  end subroutine

  subroutine hsex_current_decomposed(sys,mg,info,states,curr_decomp)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
    use communication, only: comm_summation
    use hsex_k_fft, only: hsex_backward
    use hsex_velocity, only: hsex_velocity_diagonal
    use exx_blas_threads, only: exx_blas_thread_control
    implicit none
    type(s_dft_system),intent(in) :: sys
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_orbital),intent(in) :: states
    real(8),intent(inout) :: curr_decomp(3,sys%nspin,sys%no,sys%nk)
    complex(8),allocatable :: kernel(:,:),act(:,:),dk(:,:,:,:)
    complex(8) :: x(ng,sys%no),p(na,sys%no),rp(na,sys%no,3),phase(ng)
    real(8) :: velocity(3,sys%no),local(3,sys%no,sys%nk),total(3,sys%no,sys%nk)
    integer :: ik,axis,l,first,last,blas_parent(3),blas_worker(3)
    if(.not.have_projector)error stop 'HSEX projected current before operator initialization'
    first=info%ik_s;last=info%ik_e
    allocate(kernel(size(w,1),nk),act(na*na,first:last),dk(na,na,3,first:last))
    ! Retain the original (pre-Hermitianization) R-space kernel. Recovering it
    ! from sigma loses finite-A Nyquist terms on even electronic meshes.
    do axis=1,3
      do l=1,nk
        kernel(:,l)=exchange_kernel(:,l)*cmplx(0d0,translation(axis,l),8)
      enddo
      call hsex_backward(kernel,offset,info%icomm_ko,act)
      do ik=first,last
        dk(:,:,axis,ik)=reshape(act(:,ik),[na,na])
      enddo
    enddo
    local=0d0
    blas_parent=0
#ifdef HAVE_EXX_OPENBLAS_THREADS
    call exx_blas_thread_control(1,blas_parent)
#endif
    !$omp parallel if(hsex_parallel_blas()) num_threads(hsex_k_threads(info%numk)) default(shared) &
    !$omp& private(ik,axis,x,p,rp,phase,velocity,blas_worker)
    blas_worker=0
    call exx_blas_thread_control(1,blas_worker)
    !$omp do schedule(static)
    do ik=first,last
      x=reshape(states%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,:,ik,1),[ng,sys%no])
      call project_grid(x,projector(:,ik),p)
      do axis=1,3
        phase=projector(:,ik)*(local_xyz(:,axis)+position(axis,box))
        call project_grid(x,phase,rp(:,:,axis))
      enddo
      call hsex_velocity_diagonal(p,rp,sigma(:,:,ik),dk(:,:,:,ik),velocity)
      local(:,:,ik)=velocity/(ng*sys%hvol)
    enddo
    !$omp end do
    call exx_blas_thread_control(0,blas_worker)
    !$omp end parallel
    call exx_blas_thread_control(0,blas_parent)
    call comm_summation(local,total,size(local),info%icomm_ko)
    curr_decomp(:,1,:,:)=curr_decomp(:,1,:,:)+total
  end subroutine hsex_current_decomposed

  real(8) function hsex_induced_energy()
    implicit none
    hsex_induced_energy=exchange_energy
  end function

  function hsex_xc_energy(v,rho,rho0,dv) result(e)
    use structures, only: s_scalar
    implicit none
    type(s_scalar),intent(in)::v,rho,rho0
    real(8),intent(in)::dv
    real(8)::e
    e=xc0+dv*sum(v%f*(rho%f-rho0%f))+exchange_energy
  end function
end module
