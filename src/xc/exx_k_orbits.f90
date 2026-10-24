! CPU space-group reuse of full-q exchange FFTs, in bounded first-coordinate tiles.
! This experimental path uses representative target k points and the existing
! occupied-projector callback. MPI remains outside threaded local loops.
module exx_k_orbits
  use iso_c_binding
  implicit none
  private
  include 'fftw3.f03'
  public :: apply_k_orbits,apply_spinor_k_orbits
contains
  subroutine apply_k_orbits(op,map,target,action,starts,counts,rep_start,rank,comm,transpose_tiles,fill_density,ierr)
    use communication, only: comm_summation
    use exx_k_exchange, only: exx_k_kernel,execute_kernel_fft
    use exx_symmetry, only: s_exx_symmetry_map
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    type(s_exx_symmetry_map),intent(in) :: map
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rep_start,rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,count,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: count,comm
      end subroutine
      subroutine fill_density(j,points,density)
        import c_double_complex
        implicit none
        integer,intent(in) :: j,points(:)
        complex(c_double_complex),intent(out) :: density(:,:)
      end subroutine
      subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        import c_double_complex
        implicit none
        character,intent(in) :: ta,tb
        integer,intent(in) :: m,n,k,lda,ldb,ldc
        complex(c_double_complex),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
        complex(c_double_complex),intent(inout) :: c(ldc,*)
      end subroutine
    end interface
    integer :: ng,nk,b,np,km,nmsg,nfirst,nsym,nlocal,nt,nrep,ns(3),i,j,g,s,p,r,base,nactive,lo,rows,ki
    integer :: local_error,total_error
    integer :: a,dest,source_g,first_id,global_k,source_k,local_k,op_id,ip,offset(3),delta(3),c(3),idx(3)
    integer,allocatable :: forward(:,:),wrap(:,:,:),first(:),owner(:),operation(:),kperm(:,:),lookup(:)
    integer,allocatable :: inverse(:),slot(:),points(:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:)
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),fs(:,:,:),fr(:,:,:),kw(:,:,:)
    complex(c_double_complex),allocatable :: density(:,:),t(:,:,:),row(:,:),contracted(:,:),gauge(:,:),phase(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    real(c_double) :: pi,scaled(3),v(3),length(3),angle,k0(3)
    ierr=1;action=zero
    ng=op%ng;nk=op%nk;b=op%block;np=size(starts);nsym=map%nsym
    nlocal=counts(rank+1);nt=size(target,2);nrep=size(target,3);ns=op%n*op%mesh
    local_error=0
    if(allocated(op%accelerator).or.op%nq/=nk)local_error=1
    if(size(target,1)/=ng.or.sum(counts)/=nk.or.any(shape(target)/=shape(action)))local_error=1
    if(any(op%n/=map%n).or.any(op%mesh/=op%mesh(1)))local_error=1
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    allocate(forward(ng,nsym),wrap(3,ng,nsym),first(ng),owner(ng),operation(ng),gauge(ng,nsym))
    length=map%n*map%h;pi=acos(-1d0);k0=map%full_k(:,1)
    do s=1,nsym
      do g=1,ng
        v=matmul(map%rotation(:,:,s),real(op%point(:,g),c_double)/op%n)+map%translation(:,s)
        c=nint(v*op%n);wrap(:,g,s)=floor(real(c,c_double)/op%n);c=modulo(c,op%n)
        forward(g,s)=1+c(1)+op%n(1)*(c(2)+op%n(2)*c(3))
      enddo
      do g=1,ng
        source_g=map%grid(g,s)
        angle=sum(k0*map%h*real(op%point(:,source_g)-op%point(:,g),c_double))
        gauge(g,s)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    owner=0;operation=0;nfirst=0
    do g=1,ng
      if(owner(g)/=0)cycle
      nfirst=nfirst+1;first(nfirst)=g
      do s=1,nsym
        dest=forward(g,s)
        if(owner(dest)/=0)cycle
        owner(dest)=nfirst;operation(dest)=s
      enddo
    enddo
    ! Ownership is complete; do not overlap this map with exchange scratch.
    deallocate(forward)
    allocate(kperm(nk,nsym),lookup(nk))
    lookup=0
    do i=1,nk
      scaled=(map%full_k(:,i)-k0)*length*op%mesh/(2*pi)
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));lookup(j)=i
    enddo
    do s=1,nsym;do i=1,nk
      v=matmul(transpose(map%rotation(:,:,s)),map%full_k(:,i))
      scaled=(v-k0)*length*op%mesh/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)return
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));kperm(i,s)=lookup(j)
      if(kperm(i,s)<1)return
    enddo;enddo
    deallocate(lookup)
    ! A star remains on one k owner, including all inverse-rotated targets.
    local_error=0
    do j=1,nrep
      global_k=map%first(rep_start+j-1)
      do s=1,nsym
        local_k=kperm(global_k,s)-starts(rank+1)+1
        if(local_k<1.or.local_k>nlocal)local_error=1
      enddo
    enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    km=maxval(counts)
    allocate(inverse(nk),slot(nk))
    do i=1,nk;inverse(op%order(i))=i;enddo
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
    deallocate(inverse)
    nmsg=b*ng*km
    allocate(send(nmsg*np),recv(nmsg*np),density(np*b,ng),points(np*b),t(ng,nt,nrep),row(1,ng),phase(27,nrep))
    allocate(contracted(1,nt))
    sb(1:b,1:ng,1:km,1:np)=>send;rb(1:b,1:ng,1:km,1:np)=>recv
    fs(1:b,1:ng,1:km*np)=>send;fr(1:b,1:ng,1:km*np)=>recv
    do j=1,nrep
      global_k=map%first(rep_start+j-1);local_k=global_k-starts(rank+1)+1
      do i=1,nt;t(:,i,j)=target(:,i,j)*op%phase(:,local_k);enddo
      do i=1,27
        delta=[mod(i-1,3)-1,mod((i-1)/3,3)-1,(i-1)/9-1]
        angle=sum(map%full_k(:,global_k)*length*delta)
        phase(i,j)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    if(rank==0)write(*,'(a,2i12)')'EXX_FFT_ORBITS first points/full: ',nfirst,ng
    do base=1,nfirst,np*b
      nactive=min(np*b,nfirst-base+1);points(:nactive)=first(base:base+nactive-1);send=zero
      do j=1,nlocal
        call fill_density(j,points(:nactive),density(:nactive,:))
        do p=0,np-1
          lo=p*b+1;rows=min(b,nactive-lo+1)
          if(rows>0)sb(:rows,:,j,p+1)=density(lo:lo+rows-1,:)
        enddo
      enddo
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%contiguous_fft)then
        kw(1:nk,1:b,1:ng)=>op%work
        do ki=1,nk;do g=1,ng;do r=1,b
          kw(ki,r,g)=fr(r,g,slot(ki))
        enddo;enddo;enddo
      else
        do ki=1,nk;op%work(:,:,ki)=fr(:,:,slot(ki));enddo
      endif
      call execute_kernel_fft(op,op%forward)
      rows=max(0,min(b,nfirst-base-rank*b+1))
      !$omp parallel do collapse(2) private(r,a,offset) schedule(static)
      do ki=1,nk;do g=1,ng;do r=1,rows
        a=first(base+rank*b+r-1)
        offset=modulo(op%point(:,a)-op%point(:,g)-op%shift(:,ki),ns)
        if(op%contiguous_fft)then
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        else
          op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
        endif
      enddo;enddo;enddo
      !$omp end parallel do
      call execute_kernel_fft(op,op%backward)
      if(op%contiguous_fft)then
        do ki=1,nk;do g=1,ng;do r=1,b
          fs(r,g,slot(ki))=kw(ki,r,g)
        enddo;enddo;enddo
      else
        do ki=1,nk;fs(:,:,slot(ki))=op%work(:,:,ki);enddo
      endif
      call transpose_tiles(send,recv,nmsg,comm)
      ! Reconstruct only representative-target rows; no expanded target/action.
      do dest=1,ng
        first_id=owner(dest)
        if(first_id<base.or.first_id>=base+nactive)cycle
        p=(first_id-base)/b+1;r=mod(first_id-base,b)+1
        a=first(first_id);op_id=operation(dest)
        do j=1,nrep
          global_k=map%first(rep_start+j-1);source_k=kperm(global_k,op_id)
          local_k=source_k-starts(rank+1)+1
          do g=1,ng
            source_g=map%grid(g,op_id)
            delta=wrap(:,source_g,op_id)-wrap(:,a,op_id)
            if(any(abs(delta)>1))error stop 'EXX orbit wrap exceeds one cell'
            ip=1+delta(1)+1+3*(delta(2)+1+3*(delta(3)+1))
            row(1,g)=rb(r,source_g,local_k,p)*phase(ip,j)*gauge(dest,op_id)*conjg(gauge(g,op_id))
          enddo
          call zgemm('N','N',1,nt,ng,-one/real(nk,c_double),row,1,t(1,1,j),ng,zero,contracted,1)
          action(dest,:,j)=contracted(1,:)
          local_k=global_k-starts(rank+1)+1
          action(dest,:,j)=action(dest,:,j)*conjg(op%phase(dest,local_k))
        enddo
      enddo
    enddo
    ierr=0
  end subroutine
  subroutine apply_spinor_k_orbits(op,map,target,action,starts,counts,rep_start,rank,comm,transpose_tiles,fill_density,ierr)
    use communication, only: comm_summation,comm_bcast
    use,intrinsic :: ieee_arithmetic,only:ieee_is_finite
    use iso_fortran_env,only:int64
    use exx_k_exchange, only: exx_k_kernel
    use exx_symmetry, only: s_exx_symmetry_map,symmetry_spin_rotation
    implicit none
    type(exx_k_kernel),target,intent(inout) :: op
    type(s_exx_symmetry_map),intent(in) :: map
    complex(c_double_complex),intent(in) :: target(:,:,:)
    complex(c_double_complex),intent(out) :: action(:,:,:)
    integer,intent(in) :: starts(:),counts(:),rep_start,rank,comm
    integer,intent(out) :: ierr
    interface
      subroutine transpose_tiles(send,recv,count,comm)
        import c_double_complex
        implicit none
        complex(c_double_complex),intent(in) :: send(:)
        complex(c_double_complex),intent(out) :: recv(:)
        integer,intent(in) :: count,comm
      end subroutine
      subroutine fill_density(j,points,density,spin_left,spin_right,columns)
        import c_double_complex
        implicit none
        integer,intent(in) :: j,points(:),spin_left,spin_right,columns(:)
        complex(c_double_complex),intent(out) :: density(:,:)
      end subroutine
      subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
        import c_double_complex
        implicit none
        character,intent(in) :: ta,tb
        integer,intent(in) :: m,n,k,lda,ldb,ldc
        complex(c_double_complex),intent(in) :: alpha,beta,a(lda,*),b(ldb,*)
        complex(c_double_complex),intent(inout) :: c(ldc,*)
      end subroutine
    end interface
    integer :: ng,nk,b,np,km,nmsg,nfirst,nsym,nlocal,nt,nrep,ns(3),i,j,g,s,p,r,base,nactive,lo,rows,ki
    integer :: local_error,total_error,component,spin_left,spin_right,spin_out,spin_target
    integer :: settings(13),expected(13),divisor,chunk,colbase,ncols,stat,nchars
    integer,allocatable :: columns(:)
    character(32) :: setting
    type(c_ptr) :: forward_fft,backward_fft
    integer :: a,dest,source_g,first_id,global_k,source_k,local_k,op_id,ip,offset(3),delta(3),c(3),idx(3)
    integer,allocatable :: forward(:,:),wrap(:,:,:),first(:),owner(:),operation(:),kperm(:,:),lookup(:)
    integer,allocatable :: inverse(:),slot(:),points(:)
    complex(c_double_complex),allocatable,target :: send(:),recv(:),flat(:)
    complex(c_double_complex),allocatable :: spin_rotation(:,:,:)
    complex(c_double_complex) :: coefficient
    complex(c_double_complex),pointer :: sb(:,:,:,:),rb(:,:,:,:),fs(:,:,:),fr(:,:,:),kw(:,:,:)
    complex(c_double_complex),allocatable :: density(:,:),t(:,:),row(:,:),contracted(:,:),gauge(:,:),phase(:,:)
    complex(c_double_complex),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    real(c_double) :: pi,scaled(3),v(3),angle,k0(3),stamp
    ierr=1;action=zero;op%seconds=0d0;op%threads_used=1;op%pair_fft_pairs=0_int64
    stamp=clock_now();forward_fft=c_null_ptr;backward_fft=c_null_ptr
    ng=op%ng;nk=op%nk;b=op%block;np=size(starts);nsym=map%nsym
    nt=size(target,2);nrep=size(target,3);ns=op%n*op%mesh
    local_error=0
    if(size(counts)/=np.or.rank<0.or.rank>=np.or.np<1)local_error=1
    if(min(ng,nk,b,nsym)<1.or.any(counts<0).or.any(starts<1))local_error=1
    if(allocated(op%accelerator).or.op%nq/=nk)local_error=1
    if(.not.allocated(op%phase))local_error=1
    if(size(target,1)/=2*ng.or.sum(counts)/=nk.or.any(shape(target)/=shape(action)))local_error=1
    if(any(op%n/=map%n))local_error=1
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    nlocal=counts(rank+1)
    if(starts(1)/=1.or.op%phase_start/=starts(rank+1))local_error=1
    do p=2,np
      if(starts(p)/=starts(p-1)+counts(p-1))local_error=1
    enddo
    if(size(op%phase,2)/=nlocal.or.rep_start<1.or.rep_start+nrep>map%nrep+1)local_error=1
    if(maxval(counts)>0)then
      if(int(b,int64)*ng>int(huge(nmsg),int64)/np/maxval(counts))local_error=1
    endif
    divisor=4;setting=' '
    call get_environment_variable('SALMON_EXX_ORBIT_TILE_DIVISOR',setting,length=nchars,status=stat)
    if(stat/=1)then
      if(stat/=0.or.nchars<1.or.nchars>len(setting))then
        local_error=1
      else if(verify(setting(:nchars),'0123456789')/=0)then
        local_error=1
      else
        read(setting(:nchars),*,iostat=stat)divisor
        if(stat/=0.or.divisor<1)local_error=1
      endif
    endif
    settings=[op%n,op%mesh,ng,nk,b,nsym,nt,map%nrep,divisor];expected=settings
    call comm_bcast(expected,comm,0)
    if(any(settings/=expected))local_error=1
    do j=1,nrep;do i=1,nt;do g=1,2*ng
      if(.not.ieee_is_finite(real(target(g,i,j))).or..not.ieee_is_finite(aimag(target(g,i,j))))local_error=1
    enddo;enddo;enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    allocate(forward(ng,nsym),wrap(3,ng,nsym),first(ng),owner(ng),operation(ng),gauge(ng,nsym))
    allocate(spin_rotation(2,2,nsym))
    do s=1,nsym
      call symmetry_spin_rotation(map%cart_rotation(:,:,s),spin_rotation(:,:,s),local_error)
      if(local_error/=0)return
    enddo
    pi=acos(-1d0);k0=map%full_k(:,1)
    do s=1,nsym
      do g=1,ng
        v=matmul(map%rotation(:,:,s),real(op%point(:,g),c_double)/op%n)+map%translation(:,s)
        c=nint(v*op%n);wrap(:,g,s)=floor(real(c,c_double)/op%n);c=modulo(c,op%n)
        forward(g,s)=1+c(1)+op%n(1)*(c(2)+op%n(2)*c(3))
      enddo
      do g=1,ng
        source_g=map%grid(g,s)
        angle=dot_product(k0,matmul(map%cell,real(op%point(:,source_g)-op%point(:,g),c_double)/op%n))
        gauge(g,s)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    owner=0;operation=0;nfirst=0
    do g=1,ng
      if(owner(g)/=0)cycle
      nfirst=nfirst+1;first(nfirst)=g
      do s=1,nsym
        dest=forward(g,s)
        if(owner(dest)/=0)cycle
        owner(dest)=nfirst;operation(dest)=s
      enddo
    enddo
    ! Ownership is complete; do not overlap this map with exchange scratch.
    ! Retain the inverse grid map to stream second-coordinate tiles.
    allocate(kperm(nk,nsym),lookup(nk))
    lookup=0
    do i=1,nk
      scaled=matmul(transpose(map%cell),map%full_k(:,i)-k0)*op%mesh/(2*pi)
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));lookup(j)=i
    enddo
    do s=1,nsym;do i=1,nk
      v=matmul(transpose(map%cart_rotation(:,:,s)),map%full_k(:,i))
      scaled=matmul(transpose(map%cell),v-k0)*op%mesh/(2*pi)
      if(maxval(abs(scaled-anint(scaled)))>1d-8)return
      idx=modulo(nint(scaled),op%mesh)
      j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));kperm(i,s)=lookup(j)
      if(kperm(i,s)<1)return
    enddo;enddo
    deallocate(lookup)
    ! A star remains on one k owner, including all inverse-rotated targets.
    local_error=0
    do j=1,nrep
      global_k=map%first(rep_start+j-1)
      do s=1,nsym
        local_k=kperm(global_k,s)-starts(rank+1)+1
        if(local_k<1.or.local_k>nlocal)local_error=1
      enddo
    enddo
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)return
    km=maxval(counts)
    allocate(inverse(nk),slot(nk))
    do i=1,nk;inverse(op%order(i))=i;enddo
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
    deallocate(inverse)
    chunk=1+(ng-1)/divisor;nmsg=b*chunk*km
    allocate(send(nmsg*np),recv(nmsg*np),flat(b*chunk*nk),density(np*b,chunk),points(np*b), &
      t(chunk,nt),row(1,chunk),phase(27,nrep))
    allocate(contracted(1,nt))
    sb(1:b,1:chunk,1:km,1:np)=>send;rb(1:b,1:chunk,1:km,1:np)=>recv
    fs(1:b,1:chunk,1:km*np)=>send;fr(1:b,1:chunk,1:km*np)=>recv
    do j=1,nrep
      global_k=map%first(rep_start+j-1);local_k=global_k-starts(rank+1)+1
      do i=1,27
        delta=[mod(i-1,3)-1,mod((i-1)/3,3)-1,(i-1)/9-1]
        angle=dot_product(map%full_k(:,global_k),matmul(map%cell,real(delta,c_double)))
        phase(i,j)=exp(cmplx(0d0,angle,c_double))
      enddo
    enddo
    if(op%contiguous_fft)then
      forward_fft=fftw_plan_many_dft(3,op%mesh(3:1:-1),b*chunk,flat,op%mesh(3:1:-1),1,nk, &
        flat,op%mesh(3:1:-1),1,nk,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      backward_fft=fftw_plan_many_dft(3,op%mesh(3:1:-1),b*chunk,flat,op%mesh(3:1:-1),1,nk, &
        flat,op%mesh(3:1:-1),1,nk,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    else
      forward_fft=fftw_plan_many_dft(3,op%mesh(3:1:-1),b*chunk,flat,op%mesh(3:1:-1),b*chunk,1, &
        flat,op%mesh(3:1:-1),b*chunk,1,FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      backward_fft=fftw_plan_many_dft(3,op%mesh(3:1:-1),b*chunk,flat,op%mesh(3:1:-1),b*chunk,1, &
        flat,op%mesh(3:1:-1),b*chunk,1,FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    endif
    local_error=merge(0,1,c_associated(forward_fft).and.c_associated(backward_fft))
    call comm_summation(local_error,total_error,comm)
    if(total_error/=0)goto 900
    if(rank==0)write(*,'(a,2i12)')'EXX_SPINOR_ORBIT_TILE divisor/columns: ',divisor,chunk
    if(rank==0)write(*,'(a,2i12)')'EXX_SPINOR_FFT_ORBITS first points/full: ',nfirst,ng
    if(op%profile)call mark(5)
    do base=1,nfirst,np*b
      nactive=min(np*b,nfirst-base+1);points(:nactive)=first(base:base+nactive-1)
      do colbase=1,ng,chunk
      ncols=min(chunk,ng-colbase+1);columns=[(g,g=colbase,colbase+ncols-1)]
      do component=1,4
      spin_left=1+modulo(component-1,2);spin_right=1+(component-1)/2
      send=zero
      do j=1,nlocal
        call fill_density(j,points(:nactive),density(:nactive,:ncols),spin_left,spin_right,columns)
        do p=0,np-1
          lo=p*b+1;rows=min(b,nactive-lo+1)
          if(rows>0)sb(:rows,:ncols,j,p+1)=density(lo:lo+rows-1,:ncols)
        enddo
      enddo
      do g=1,size(send)
        if(.not.ieee_is_finite(real(send(g))).or..not.ieee_is_finite(aimag(send(g))))local_error=1
      enddo
      call comm_summation(local_error,total_error,comm)
      if(total_error/=0)goto 900
      if(op%profile)call mark(1)
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark(2)
      if(op%contiguous_fft)then
        kw(1:nk,1:b,1:chunk)=>flat
        do ki=1,nk;do g=1,chunk;do r=1,b
          kw(ki,r,g)=fr(r,g,slot(ki))
        enddo;enddo;enddo
      else
        do ki=1,nk
          flat((ki-1)*b*chunk+1:ki*b*chunk)=reshape(fr(:,:,slot(ki)),[b*chunk])
        enddo
      endif
      if(op%profile)call mark(5)
      call fftw_execute_dft(forward_fft,flat,flat)
      if(op%profile)call mark(3)
      rows=max(0,min(b,nfirst-base-rank*b+1))
      do ki=1,nk;do g=1,ncols;do r=1,rows
        a=first(base+rank*b+r-1)
        offset=modulo(op%point(:,a)-op%point(:,colbase+g-1)-op%shift(:,ki),ns)
        if(op%contiguous_fft)then
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        else
          flat((ki-1)*b*chunk+(g-1)*b+r)=flat((ki-1)*b*chunk+(g-1)*b+r)* &
            op%kernel(offset(1),offset(2),offset(3))
        endif
      enddo;enddo;enddo
      if(op%profile)call mark(4)
      call fftw_execute_dft(backward_fft,flat,flat)
      if(op%profile)call mark(3)
      op%pair_fft_pairs=op%pair_fft_pairs+int(rows,int64)*ncols
      if(op%contiguous_fft)then
        do ki=1,nk;do g=1,chunk;do r=1,b
          fs(r,g,slot(ki))=kw(ki,r,g)
        enddo;enddo;enddo
      else
        do ki=1,nk
          fs(:,:,slot(ki))=reshape(flat((ki-1)*b*chunk+1:ki*b*chunk),[b,chunk])
        enddo
      endif
      if(op%profile)call mark(5)
      call transpose_tiles(send,recv,nmsg,comm)
      if(op%profile)call mark(2)
      ! Linearity permits one spin block at a time: U rho U-dagger is
      ! contracted componentwise, avoiding four retained communication buffers.
      do dest=1,ng
        first_id=owner(dest)
        if(first_id<base.or.first_id>=base+nactive)cycle
        p=(first_id-base)/b+1;r=mod(first_id-base,b)+1
        a=first(first_id);op_id=operation(dest)
        do j=1,nrep
          global_k=map%first(rep_start+j-1);source_k=kperm(global_k,op_id)
          local_k=source_k-starts(rank+1)+1
          do g=1,ncols
            source_g=colbase+g-1
            i=forward(source_g,op_id)
            delta=wrap(:,source_g,op_id)-wrap(:,a,op_id)
            if(any(abs(delta)>1))error stop 'EXX spinor orbit wrap exceeds one cell'
            ip=1+delta(1)+1+3*(delta(2)+1+3*(delta(3)+1))
            row(1,g)=rb(r,g,local_k,p)*phase(ip,j)*gauge(dest,op_id)*conjg(gauge(i,op_id))
          enddo
          local_k=global_k-starts(rank+1)+1
          do spin_out=1,2;do spin_target=1,2
            coefficient=spin_rotation(spin_out,spin_left,op_id)*conjg(spin_rotation(spin_target,spin_right,op_id))
            if(coefficient==zero)cycle
            do g=1,ncols
              i=forward(colbase+g-1,op_id)
              t(g,:)=target(i+(spin_target-1)*ng,:,j)*op%phase(i,local_k)
            enddo
            call zgemm('N','N',1,nt,ncols,-coefficient/real(nk,c_double),row,1, &
              t,chunk,zero,contracted,1)
            action(dest+(spin_out-1)*ng,:,j)=action(dest+(spin_out-1)*ng,:,j)+ &
              contracted(1,:)*conjg(op%phase(dest,local_k))
          enddo;enddo
        enddo
      enddo
      if(op%profile)call mark(6)
      enddo
      enddo
    enddo
    ierr=0
900 continue
    if(c_associated(forward_fft))call fftw_destroy_plan(forward_fft)
    if(c_associated(backward_fft))call fftw_destroy_plan(backward_fft)
  contains
    real(c_double) function clock_now()
      integer(int64) :: tick,rate
      call system_clock(tick,rate)
      clock_now=real(tick,c_double)/real(rate,c_double)
    end function
    subroutine mark(stage)
      integer,intent(in) :: stage
      real(c_double) :: now
      now=clock_now();op%seconds(stage)=op%seconds(stage)+now-stamp;stamp=now
    end subroutine
  end subroutine
end module
