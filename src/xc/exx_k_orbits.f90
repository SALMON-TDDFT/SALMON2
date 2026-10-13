! CPU space-group reuse of full-q exchange FFTs, in bounded first-coordinate tiles.
! This experimental path uses representative target k points and the existing
! occupied-projector callback. MPI remains outside threaded local loops.
module exx_k_orbits
  use communication, only: comm_summation
  use iso_c_binding
  use exx_k_exchange, only: exx_k_kernel
  use exx_symmetry, only: s_exx_symmetry_map
  implicit none
  private
  include 'fftw3.f03'
  public :: apply_k_orbits
contains
  subroutine apply_k_orbits(op,map,target,action,starts,counts,rep_start,rank,comm,transpose_tiles,fill_density,ierr)
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
    allocate(contracted(1,size(target,2)))
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
    allocate(kperm(nk,nsym),lookup(nk),inverse(nk),slot(nk))
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
    km=maxval(counts);nmsg=b*ng*km
    allocate(send(nmsg*np),recv(nmsg*np),density(np*b,ng),points(np*b),t(ng,nt,nrep),row(1,ng),phase(27,nrep))
    sb(1:b,1:ng,1:km,1:np)=>send;rb(1:b,1:ng,1:km,1:np)=>recv
    fs(1:b,1:ng,1:km*np)=>send;fr(1:b,1:ng,1:km*np)=>recv
    do i=1,nk;inverse(op%order(i))=i;enddo
    do p=1,np;do j=1,counts(p)
      slot(inverse(starts(p)+j-1))=j+(p-1)*km
    enddo;enddo
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
      call fftw_execute_dft(op%forward,op%work,op%work)
      rows=max(0,min(b,nfirst-base-rank*b+1))
      do ki=1,nk;do g=1,ng;do r=1,rows
        a=first(base+rank*b+r-1)
        offset=modulo(op%point(:,a)-op%point(:,g)-op%shift(:,ki),ns)
        if(op%contiguous_fft)then
          kw(ki,r,g)=kw(ki,r,g)*op%kernel(offset(1),offset(2),offset(3))
        else
          op%work(r,g,ki)=op%work(r,g,ki)*op%kernel(offset(1),offset(2),offset(3))
        endif
      enddo;enddo;enddo
      call fftw_execute_dft(op%backward,op%work,op%work)
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
end module
