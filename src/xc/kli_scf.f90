module kli_scf
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use structures
  use kli_periodic
  implicit none
  private
  public :: apply_kli_scf,kli_write,kli_read,apply_frozen,kli_diagnostics
  real(8),allocatable,save::cached_v(:,:),cached_e(:)
  type(kli_frame),save :: frame(2)
contains
  subroutine apply_kli_scf(sys,mg,info,orb,vxc,exc,horb,predictor)
    use salmon_global, only: num_kgrid,kli_local_tol,kli_solve_tol,kli_hartree_tol,kli_maxiter
    use communication,only:comm_summation
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    type(s_scalar),intent(inout)::vxc(:)
    real(8),intent(inout)::exc
    complex(8),allocatable::part(:,:,:),psi(:,:,:)
    type(s_orbital),optional,intent(in)::horb
    logical,optional,intent(in)::predictor
    complex(8),allocatable::prev(:,:,:),gu(:,:,:),hp(:,:,:),hsum(:,:,:)
    integer,allocatable::oldcut(:,:)
    logical::trial
    real(8),allocatable::v(:)
    real(8)::energy,locerr,edge,occ
    integer::s,n(3),nb,ng,ik,io,ix,iy,iz,g,iter,ffts,istat
    n=mg%ie-mg%is+1;ng=product(n)
    if(info%isize_r/=1.or.sys%if_real_orbital)error stop 'KLI requires a full complex primitive grid on each k rank'
    if(any(abs(sys%wtk-1d0/sys%nk)>1d-12))error stop 'KLI requires equal weights'
    trial=.false.;if(present(predictor))trial=predictor
    occ=real(3-sys%nspin,8)
    allocate(v(ng))
    if(.not.allocated(cached_v))allocate(cached_v(ng,sys%nspin),cached_e(sys%nspin))
    do s=1,sys%nspin
      nb=count(sys%rocc(:,1,s)==occ)
      if(nb==0)error stop 'KLI requires a nonempty integer occupied spin space'
      if(any(sys%rocc(1:nb,:,s)/=occ).or.any(sys%rocc(nb+1:,:,s)/=0d0))error stop 'KLI requires integer valence occupations'
      allocate(part(ng,nb,sys%nk),psi(ng,nb,sys%nk));part=0
      do ik=info%ik_s,info%ik_e;do io=info%io_s,min(nb,info%io_e)
        g=0
        do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
          g=g+1;part(g,io,ik)=orb%zwf(ix,iy,iz,s,io,ik,1)
        enddo;enddo;enddo
      enddo;enddo
      call comm_summation(part,psi,size(part),info%icomm_ko)
      if(present(horb))then
        if(trial)oldcut=frame(s)%image_shift
        prev=frame(s)%w%previous;gu=frame(s)%w%gauge
        allocate(hp(ng,nb,sys%nk),hsum(ng,nb,sys%nk));hp=0
        do ik=info%ik_s,info%ik_e;do io=info%io_s,min(nb,info%io_e)
          hp(:,io,ik)=reshape(horb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),s,io,ik,1),[ng])
        enddo;enddo
        call comm_summation(hp,hsum,size(hp),info%icomm_ko)
        do ik=1,sys%nk
          call canonical(psi(:,:,ik),hsum(:,:,ik),sys%hvol,matmul(prev(:,:,ik),conjg(transpose(gu(:,:,ik)))))
        enddo
        deallocate(hp,hsum)
      endif
      call evaluate_kli(frame(s),psi,n,num_kgrid,sys%primitive_a,sys%vec_k,sys%hvol,kli_local_tol,kli_solve_tol, &
                        kli_hartree_tol,kli_maxiter,v,energy,locerr,edge,iter,ffts,istat)
      if(info%id_ko==0)write(*,'(a,4i8,3es18.9)')'KLI spin/status/iter/FFT E/loc/boundary: ',s,istat,iter,ffts,energy,locerr,edge
      if(info%id_ko==0)write(*,'(a,7es16.7)')'KLI orth/rho/Slater/row/compat/solve/box: ',frame(s)%audit
      if(istat/=0)error stop 'KLI rejected gauge, density, localization, or linear solve'
      g=0
      do iz=mg%is(3),mg%ie(3);do iy=mg%is(2),mg%ie(2);do ix=mg%is(1),mg%ie(1)
        g=g+1;vxc(s)%f(ix,iy,iz)=vxc(s)%f(ix,iy,iz)+v(g)
      enddo;enddo;enddo
      if(trial)then
        frame(s)%image_shift=oldcut
        frame(s)%w%previous=prev;frame(s)%w%gauge=gu
      else
        cached_v(:,s)=v;cached_e(s)=occ*energy
      endif
      exc=exc+occ*energy
      deallocate(part,psi)
    enddo
  end subroutine
  subroutine canonical(psi,hpsi,dv,ref)
    use exx_wannier_gauge, only:gauge_transport
    complex(8),intent(inout)::psi(:,:)
    complex(8),intent(in)::hpsi(:,:),ref(:,:)
    real(8),intent(in)::dv
    complex(8),allocatable::h(:,:),work(:),align(:,:,:)
    real(8),allocatable::ev(:),rw(:)
    real(8)::sigma
    integer::n,l,r,nb,ng,s
    ng=size(psi,1);n=size(psi,2)
    allocate(h(n,n),work(4*n),ev(n),rw(3*n))
    h=matmul(conjg(transpose(psi)),hpsi)*dv;h=.5d0*(h+conjg(transpose(h)))
    call zheev('V','U',n,h,n,ev,work,size(work),rw,s)
    if(s/=0)error stop 'KLI projected H diagonalization failed'
    psi=matmul(psi,h);l=1
    do while(l<=n)
      r=l
      do while(r<n)
        if(ev(r+1)-ev(l)>1d-7)exit
        r=r+1
      enddo
      nb=r-l+1
      if(nb>1)then
        allocate(align(nb,nb,1))
        call gauge_transport(reshape(psi(:,l:r),[ng,nb,1]),reshape(ref(:,l:r),[ng,nb,1]),dv,align,sigma,s)
        if(s/=0)error stop 'KLI degenerate-block tracking lost rank'
        psi(:,l:r)=matmul(psi(:,l:r),align(:,:,1));deallocate(align)
      endif
      l=r+1
    enddo
  end subroutine
  subroutine kli_diagnostics(sys,mg,info,orb,rho,rho0)
    use communication,only:comm_get_max
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    type(s_scalar),intent(in)::rho(:),rho0(:)
    real(8),allocatable,save::v0(:,:)
    real(8)::orth,dens,orthmax(1)
    complex(8),allocatable::a(:,:),over(:,:)
    integer::k,s,nb,i,ng
    logical::first
    first=.not.allocated(v0)
    if(first)v0=cached_v
    if(info%isize_o/=1)error stop "KLI RT requires an unsplit occupied band space"
    orth=0;dens=0;ng=product(mg%ie-mg%is+1)
    do s=1,sys%nspin
      dens=max(dens,maxval(abs(rho(s)%f-rho0(s)%f)))
      nb=count(sys%rocc(:,1,s)>0)
      do k=info%ik_s,info%ik_e
        a=reshape(orb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),s,1:nb,k,1),[ng,nb])
        if(first)then
          over=sys%hvol*matmul(conjg(transpose(a)),frame(s)%w%previous(:,:,k))
          over=matmul(over,conjg(transpose(over)))
          do i=1,nb
            over(i,i)=over(i,i)-1d0
          enddo
          if(.not.(maxval(abs(over))<1d-7))error stop 'KLI checkpoint belongs to a different occupied projector'
        endif
        over=sys%hvol*matmul(conjg(transpose(a)),a)
        do i=1,nb
          over(i,i)=over(i,i)-1d0
        enddo
        orth=max(orth,maxval(abs(over)))
      enddo
    enddo
    call comm_get_max([orth],orthmax,1,info%icomm_ko);orth=orthmax(1)
    if(info%id_ko==0)write(*,'(a,3es20.10)')'KLI RT density/orthogonality/potential: ',dens,orth,maxval(abs(cached_v-v0))
  end subroutine
  subroutine kli_write(dir,sys,mg)
    character(*),intent(in)::dir
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    integer::u,s
    open(newunit=u,file=trim(dir)//'kli-v1.bin',access='stream',form='unformatted',status='replace')
    write(u)'KLI1',mg%ie-mg%is+1,sys%nspin,sys%nk,sys%no,sys%nion,sys%primitive_a,sys%vec_k,sys%rocc,sys%Rion
    write(u)cached_v,cached_e
    do s=1,sys%nspin
      write(u)frame(s)%w%previous,frame(s)%w%gauge,frame(s)%image_shift
    enddo
    close(u)
  end subroutine

  subroutine kli_read(dir,sys,mg)
    character(*),intent(in)::dir
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    character(4)::magic
    integer::u,dims(7),n(3),nb,s,ng,i
    real(8)::cell(3,3),k(3,sys%nk),ions(3,sys%nion)
    real(8),allocatable::occ(:,:,:)
    n=mg%ie-mg%is+1;ng=product(n)
    open(newunit=u,file=trim(dir)//'kli-v1.bin',access='stream',form='unformatted',status='old',action='read')
    read(u)magic,dims
    if(magic/='KLI1'.or.any(dims(1:5)/=[n,sys%nspin,sys%nk]).or.dims(6)<sys%no.or.dims(7)/=sys%nion) &
      error stop 'KLI checkpoint grid mismatch'
    allocate(occ(dims(6),sys%nk,sys%nspin))
    read(u)cell,k,occ,ions
    if(any(cell/=sys%primitive_a).or.any(k/=sys%vec_k).or.any(occ(:sys%no,:,:)/=sys%rocc).or. &
       any(occ(sys%no+1:,:,:)/=0d0).or.any(ions/=sys%Rion)) &
      error stop 'KLI checkpoint crystal or occupation mismatch'
    allocate(cached_v(ng,sys%nspin),cached_e(sys%nspin))
    read(u)cached_v,cached_e
    do s=1,sys%nspin
      if(.not.ieee_is_finite(cached_e(s)))error stop 'KLI checkpoint contains nonfinite energy'
      do i=1,ng
        if(.not.ieee_is_finite(cached_v(i,s)))error stop 'KLI checkpoint contains nonfinite potential'
      enddo
      nb=count(occ(:,1,s)>0)
      allocate(frame(s)%w%previous(ng,nb,sys%nk),frame(s)%w%gauge(nb,nb,sys%nk),frame(s)%image_shift(3,nb))
      read(u)frame(s)%w%previous,frame(s)%w%gauge,frame(s)%image_shift
    enddo
    close(u)
  end subroutine
  subroutine apply_frozen(vxc,exc)
    type(s_scalar),intent(inout)::vxc(:)
    real(8),intent(inout)::exc
    integer::s
    do s=1,size(vxc)
      vxc(s)%f=vxc(s)%f+reshape(cached_v(:,s),shape(vxc(s)%f))
    enddo
    exc=exc+sum(cached_e)
  end subroutine
end module
