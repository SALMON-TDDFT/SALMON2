module kli_scf
  use exx_symmetry, only: s_exx_symmetry_map
  use kli_periodic, only: kli_frame
  implicit none
  private
  public :: kli_export_hamiltonian,kli_check_reference,kli_get_reference
  public :: apply_kli_scf,kli_write,kli_read,apply_frozen,kli_diagnostics
  ! Shared between SCF updates, checkpoint I/O and RT predictor/corrector calls.
  real(8),allocatable,save::cached_v(:,:),cached_e(:)
  type(s_exx_symmetry_map),save :: star
  type(kli_frame),save :: frame(2)
contains
  ! Grid-column export for the direct Sternheimer reference solver, one file per k.
  ! Probe hpsi itself, including finite differences and the actual pseudopotential.
  subroutine kli_export_hamiltonian(sys,mg,info,vl,stencil,srg,ppg,x,hx,tbxc)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_scalar,s_stencil,s_pp_grid,s_orbital
    use sendrecv_grid, only: s_sendrecv_grid
    use hamiltonian, only: hpsi
    use salmon_global, only: cval
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_scalar),intent(in)::vl(:)
    type(s_stencil),intent(in)::stencil
    type(s_sendrecv_grid),intent(inout)::srg
    type(s_pp_grid),intent(in)::ppg
    type(s_orbital),intent(inout)::x,hx
    complex(8),allocatable :: column(:)
    integer::u(info%ik_s:info%ik_e),k,g,j,col,ix,iy,iz,n(3),ng
    type(s_scalar),optional,intent(in)::tbxc
    character(48)::name
    if(info%isize_r/=1.or.info%isize_o/=1.or.sys%nspin/=1.or.any(sys%vec_Ac/=0d0)) &
      error stop 'Sternheimer export needs field-free closed-shell k-only layout'
    n=mg%ie-mg%is+1;ng=product(n)
    do k=info%ik_s,info%ik_e
      write(name,'("kli-hamiltonian-",i6.6,".bin")')k
      open(newunit=u(k),file=trim(name),status='new',form='unformatted',access='stream',action='write')
      if(present(tbxc))then
        write(u(k))'TSH2',n,k,sys%vec_k(:,k),sys%primitive_a,cval,tbxc%f
      else
        write(u(k))'KSH2',n,k,sys%vec_k(:,k),sys%primitive_a,cached_e(1),cached_v(:,1)
      endif
    enddo
    allocate(column(ng))
    do g=1,ng,sys%no
      x%zwf=0
      do j=1,min(sys%no,ng-g+1)
        col=g+j-2
        ix=mg%is(1)+mod(col,n(1));iy=mg%is(2)+mod(col/n(1),n(2));iz=mg%is(3)+col/(n(1)*n(2))
        x%zwf(ix,iy,iz,1,j,:,1)=1
      enddo
      x%update_zwf_overlap=.false.
      call hpsi(x,hx,info,mg,vl,sys,stencil,srg,ppg)
      do k=info%ik_s,info%ik_e
        do j=1,min(sys%no,ng-g+1)
          column=reshape(hx%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,j,k,1),[ng])
          write(u(k))column(g+j-1:ng)
        enddo
      enddo
    enddo
    do k=info%ik_s,info%ik_e
      close(u(k))
    enddo
  end subroutine

  subroutine kli_get_reference(v,e)
    implicit none
    real(8),intent(out) :: v(:),e
    if(.not.allocated(cached_v).or..not.allocated(cached_e))error stop 'KLI reference unavailable'
    if(size(cached_v,1)/=size(v))error stop 'KLI reference grid mismatch'
    v=cached_v(:,1);e=cached_e(1)
  end subroutine kli_get_reference

  subroutine kli_check_reference(v,e)
    implicit none
    real(8),intent(in)::v(:),e
    if(any(v/=cached_v(:,1)).or.e/=cached_e(1))error stop 'W does not match static KLI restart'
  end subroutine

  subroutine prepare_star(sys,mg)
    use structures, only: s_dft_system,s_rgrid
    use sym_sub, only: use_symmetry,SymMatA,SymMatB,symmetry_validate_field
    use exx_symmetry, only: symmetry_init,symmetry_validate_atoms
    use salmon_global, only: num_kgrid,theory,epdir_re1,epdir_im1
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    integer::status
    if(.not.use_symmetry.or.star%nfull/=0)return
    call symmetry_init(star,mg%num,sys%hgs,sys%vec_k(:,:sys%nk),sys%wtk(:sys%nk), &
                       SymMatA,SymMatB,num_kgrid,status,sys%primitive_a)
    if(status/=0)error stop 'KLI: incompatible crystal symmetry, grid or k mesh'
    call symmetry_validate_atoms(star,sys%Rion,sys%kion,status)
    if(status/=0)error stop 'KLI: symmetry does not preserve ions'
    if(theory/='dft')then
      call symmetry_validate_field(epdir_re1)
      call symmetry_validate_field(epdir_im1)
    endif
  end subroutine

  subroutine expand_star(orbitals)
    use sym_sub, only: use_symmetry
    use exx_symmetry, only: symmetry_transform
    implicit none
    complex(8),allocatable,intent(inout)::orbitals(:,:,:)
    complex(8),allocatable::full(:,:,:)
    integer::k,status
    if(.not.use_symmetry)return
    allocate(full(size(orbitals,1),size(orbitals,2),star%nfull))
    do k=1,star%nfull
      call symmetry_transform(star,k,1,orbitals(:,:,star%owner(k)),full(:,:,k),status)
      if(status/=0)error stop 'KLI: failed Bloch star expansion'
    enddo
    call move_alloc(full,orbitals)
  end subroutine

  subroutine apply_kli_scf(sys,mg,info,orb,vxc,exc,horb,predictor)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
    use kli_periodic, only: evaluate_kli
    use sym_sub, only: use_symmetry
    use ieee_arithmetic, only: ieee_is_finite
    use salmon_global, only: num_kgrid,kli_local_tol,kli_solve_tol,kli_hartree_tol,kli_maxiter,yn_kli_rt_box_check
    use communication, only: comm_is_root,comm_summation
    implicit none
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
    real(8),allocatable::v(:),full_k(:,:)
    real(8)::energy,locerr,edge,occ
    integer::s,n(3),nb,ng,ik,io,ix,iy,iz,g,iter,ffts,istat
    n=mg%ie-mg%is+1;ng=product(n)
    if(info%isize_r/=1.or.sys%if_real_orbital)error stop 'KLI requires a full complex primitive grid on each k rank'
    call prepare_star(sys,mg)
    if(.not.use_symmetry.and.any(abs(sys%wtk(:sys%nk)-1d0/sys%nk)>1d-12))error stop 'KLI requires equal weights'
    if(use_symmetry)then
      full_k=star%full_k
    else
      full_k=sys%vec_k(:,:sys%nk)
    endif
    trial=.false.;if(present(predictor))trial=predictor
    occ=real(3-sys%nspin,8)
    allocate(v(ng))
    if(use_symmetry.and.comm_is_root(info%id_ko).and..not.allocated(cached_v)) &
      write(*,*) "KLI representative/full k points:",sys%nk,star%nfull
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
      call expand_star(psi)
      if(present(horb))then
        if(trial)oldcut=frame(s)%image_shift
        prev=frame(s)%w%previous;gu=frame(s)%w%gauge
        allocate(hp(ng,nb,sys%nk),hsum(ng,nb,sys%nk));hp=0
        do ik=info%ik_s,info%ik_e;do io=info%io_s,min(nb,info%io_e)
          hp(:,io,ik)=reshape(horb%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),s,io,ik,1),[ng])
        enddo;enddo
        call comm_summation(hp,hsum,size(hp),info%icomm_ko)
        call expand_star(hsum)
        do ik=1,size(psi,3)
          call canonical(psi(:,:,ik),hsum(:,:,ik),sys%hvol,matmul(prev(:,:,ik),conjg(transpose(gu(:,:,ik)))))
        enddo
        deallocate(hp,hsum)
      endif
      call evaluate_kli(frame(s),psi,n,num_kgrid,sys%primitive_a,full_k,sys%hvol,kli_local_tol,kli_solve_tol, &
                        kli_hartree_tol,kli_maxiter,v,energy,locerr,edge,iter,ffts,istat,info%icomm_ko,info%id_ko,info%isize_ko, &
                        check_box=.not.present(horb).or.yn_kli_rt_box_check=='y')
      if(comm_is_root(info%id_ko))write(*,'(a,4i8,3es18.9)') &
        'KLI spin/status/iter/FFT E/loc/boundary: ',s,istat,iter,ffts,energy,locerr,edge
      if(comm_is_root(info%id_ko))write(*,'(a,7es16.7)')'KLI orth/rho/Slater/row/compat/solve/box: ',frame(s)%audit
      if(istat/=0)error stop 'KLI rejected gauge, density, localization, or linear solve'
      if(present(horb).and.comm_is_root(info%id_ko))then
        write(*,'(a,1x,a,1x,i0,3es23.14)') 'KLI RT localization',merge('predictor','corrector',trial), &
          iter,frame(s)%initial_pz_residual,locerr,frame(s)%w%min_singular
      endif
      if(use_symmetry)then
        ! The local Hamiltonian used on a representative must commute with
        ! the retained operations. Average numerical gauge noise only.
        block
          real(8)::average(ng),defect
          integer::op
          average=0d0
          do op=1,star%nsym
            average=average+v(star%grid(:,op))
          enddo
          average=average/star%nsym
          defect=maxval(abs(v-average))
          if(comm_is_root(info%id_ko))write(*,'(a,es14.6)') 'KLI potential symmetry defect [Ha]: ',defect
          if(.not.ieee_is_finite(defect).or.defect>kli_local_tol) &
            error stop 'KLI localized solution breaks requested symmetry; use full k mesh'
          v=average
        end block
      endif
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
    implicit none
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
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital,s_scalar
    use sym_sub, only: use_symmetry
    use communication, only: comm_is_root,comm_get_max
    implicit none
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    type(s_parallel_info),intent(in)::info
    type(s_orbital),intent(in)::orb
    type(s_scalar),intent(in)::rho(:),rho0(:)
    real(8),allocatable,save::v0(:,:)
    real(8)::orth,dens,orthmax(1)
    complex(8),allocatable::a(:,:),over(:,:)
    integer::k,s,nb,i,ng,kfull
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
          kfull=k
          if(use_symmetry)kfull=star%first(k)
          over=sys%hvol*matmul(conjg(transpose(a)),frame(s)%w%previous(:,:,kfull))
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
    if(comm_is_root(info%id_ko))write(*,'(a,3es20.10)')'KLI RT density/orthogonality/potential: ',dens,orth,maxval(abs(cached_v-v0))
  end subroutine
  subroutine kli_write(dir,sys,mg)
    use structures, only: s_dft_system,s_rgrid
    use sym_sub, only: use_symmetry,SymMatA
    implicit none
    character(*),intent(in)::dir
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    integer::u,s
    open(newunit=u,file=trim(dir)//'kli-v1.bin',access='stream',form='unformatted',status='replace')
    write(u)merge('KLS1','KLI1',use_symmetry),mg%ie-mg%is+1,sys%nspin,sys%nk,sys%no,sys%nion, &
      sys%primitive_a,sys%vec_k(:,:sys%nk),sys%rocc,sys%Rion
    if(use_symmetry)write(u)star%nfull,star%nsym,star%full_k,SymMatA
    write(u)cached_v,cached_e
    do s=1,sys%nspin
      write(u)frame(s)%w%previous,frame(s)%w%gauge,frame(s)%image_shift
    enddo
    close(u)
  end subroutine

  subroutine kli_read(dir,sys,mg)
    use structures, only: s_dft_system,s_rgrid
    use sym_sub, only: use_symmetry,SymMatA
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    character(*),intent(in)::dir
    type(s_dft_system),intent(in)::sys
    type(s_rgrid),intent(in)::mg
    character(4)::magic
    integer::u,dims(7),n(3),nb,s,ng,i,nfull,nsym
    real(8)::cell(3,3),k(3,sys%nk),ions(3,sys%nion)
    real(8),allocatable::occ(:,:,:),saved_k(:,:),ops(:,:,:)
    n=mg%ie-mg%is+1;ng=product(n)
    open(newunit=u,file=trim(dir)//'kli-v1.bin',access='stream',form='unformatted',status='old',action='read')
    read(u)magic,dims
    if(magic/=merge('KLS1','KLI1',use_symmetry).or.any(dims(1:5)/=[n,sys%nspin,sys%nk]).or.dims(6)<sys%no.or.dims(7)/=sys%nion) &
      error stop 'KLI checkpoint grid mismatch'
    allocate(occ(dims(6),sys%nk,sys%nspin))
    read(u)cell,k,occ,ions
    if(any(cell/=sys%primitive_a).or.any(k/=sys%vec_k(:,:sys%nk)).or.any(occ(:sys%no,:,:)/=sys%rocc).or. &
       any(occ(sys%no+1:,:,:)/=0d0).or.any(ions/=sys%Rion)) &
      error stop 'KLI checkpoint crystal or occupation mismatch'
    nfull=sys%nk
    if(use_symmetry)then
      call prepare_star(sys,mg)
      read(u)nfull,nsym
      if(nfull/=star%nfull.or.nsym/=star%nsym)error stop 'KLI restart symmetry size mismatch'
      allocate(saved_k(3,nfull),ops(3,4,nsym))
      read(u)saved_k,ops
      if(any(saved_k/=star%full_k).or.any(ops/=SymMatA))error stop 'KLI restart symmetry changed'
    endif
    allocate(cached_v(ng,sys%nspin),cached_e(sys%nspin))
    read(u)cached_v,cached_e
    do s=1,sys%nspin
      if(.not.ieee_is_finite(cached_e(s)))error stop 'KLI checkpoint contains nonfinite energy'
      do i=1,ng
        if(.not.ieee_is_finite(cached_v(i,s)))error stop 'KLI checkpoint contains nonfinite potential'
      enddo
      nb=count(occ(:,1,s)>0)
      allocate(frame(s)%w%previous(ng,nb,nfull),frame(s)%w%gauge(nb,nb,nfull),frame(s)%image_shift(3,nb))
      read(u)frame(s)%w%previous,frame(s)%w%gauge,frame(s)%image_shift
    enddo
    close(u)
  end subroutine
  subroutine apply_frozen(vxc,exc)
    use structures, only: s_scalar
    implicit none
    type(s_scalar),intent(inout)::vxc(:)
    real(8),intent(inout)::exc
    integer::s
    do s=1,size(vxc)
      vxc(s)%f=vxc(s)%f+reshape(cached_v(:,s),shape(vxc(s)%f))
    enddo
    exc=exc+sum(cached_e)
  end subroutine
end module
