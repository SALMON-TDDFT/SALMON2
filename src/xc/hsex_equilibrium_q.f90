! Occupied-only, k-distributed finite-frequency chi0 of a stationary insulator.
module hsex_equilibrium_q
  implicit none
  private
  public :: equilibrium_q_response
  abstract interface
    subroutine sector_action(source,target,x,y,status)
      implicit none
      integer,intent(in) :: source,target
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine equilibrium_q_response(n,mesh,qp,cell,k,occupied,hocc,ks,ke,icomm,rank,np,capacity,tol,z, &
      apply,chi,residual,status)
    use communication, only: comm_summation,comm_bcast,comm_get_max,comm_is_root
    use hsex_screen_math, only: screen_geometry,screen_finite
    use hsex_frozen_response, only: frozen_solve,s_frozen_solution,s_frozen_stats
    implicit none
    integer,intent(in) :: n(3),mesh(3),qp(3),ks,ke,icomm,rank,np,capacity
    real(8),intent(in) :: cell(3,3),k(:,:),tol
    complex(8),intent(in) :: occupied(:,:,:),hocc(:,:,:),z(:)
    procedure(sector_action) :: apply
    complex(8),intent(out) :: chi(:,:)
    real(8),intent(out) :: residual
    integer,intent(out) :: status
    type(s_frozen_solution) :: solution
    type(s_frozen_stats) :: statistics
    complex(8),allocatable :: psi(:,:,:),ghost(:,:,:,:),buffer(:,:),matrix(:,:),work(:),rhs(:,:),vertices(:,:)
    complex(8),allocatable :: source_c(:,:),scalar(:,:),shifts(:),local(:,:),eigenvector(:)
    real(8),allocatable :: energy(:,:),e(:),rw(:),frac(:,:)
    integer,allocatable :: owners(:),lookup(:),order(:,:),target_k(:,:)
    integer :: sizes(2,np),ranges(2,np),ng,nv,nlocal,nk,nz,code,i,j,a,b,s,side,sign,target,source,trial
    integer :: coord(3),fold(3),points(3,3),active_source,active_target,unit,actions
    real(8) :: rec(3,3),volume,position(3),q(3),source_energy,worst(1),global_worst(1)
    complex(8) :: value
    status=1;chi=0d0;residual=huge(1d0);ng=product(n);nk=product(mesh);nv=size(occupied,2)
    nlocal=max(0,ke-ks+1);nz=size(z);code=0;actions=0
    if(nv<1.or.ng<1.or.nk<1.or.nz<1.or.capacity<1.or.tol<=0d0)code=1
    if(any(shape(occupied)/=[ng,nv,nlocal]).or.any(shape(hocc)/=[nv,nv,nlocal]))code=1
    if(any(shape(k)/=[3,nk]).or.any(shape(chi)/=[nz,3]))code=1
    if(.not.screen_finite(occupied).or..not.screen_finite(hocc).or..not.screen_finite(k))code=1
    if(.not.screen_finite(z).or.any(aimag(z)<=0d0).or..not.screen_finite([tol]))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    call screen_geometry(cell,rec,volume,code)
    if(code/=0)return
    allocate(psi(ng,nv,nlocal),ghost(ng,nv,nlocal,2),buffer(ng,nv),matrix(nv,nv),work(max(1,4*nv)), &
      rhs(ng,1),vertices(ng,3),source_c(1,1),scalar(1,1),shifts(nz),local(nz,3), &
      eigenvector(ng),energy(nv,nlocal),e(nv),rw(max(1,3*nv)),frac(3,nk),owners(nk),lookup(nk), &
      order(3,nk),target_k(nlocal,2))
    do j=1,nlocal
      matrix=hocc(:,:,j)
      if(maxval(abs(matrix-conjg(transpose(matrix))))>1d-9)code=1
      call zheev('V','L',nv,matrix,nv,e,work,size(work),rw,i)
      if(i/=0)code=1
      psi(:,:,j)=matmul(occupied(:,:,j),matrix);energy(:,j)=e
    enddo
    ranges=0;ranges(:,rank+1)=[ks,ke];call comm_summation(ranges,sizes,size(ranges),icomm)
    owners=-1
    do j=1,np
      do i=sizes(1,j),sizes(2,j)
        if(i<1.or.i>nk)then
          code=1;cycle
        endif
        if(owners(i)/=-1)code=1
        owners(i)=j-1
      enddo
    enddo
    if(any(owners<0))code=1
    frac=matmul(transpose(cell),k)/(2d0*acos(-1d0));lookup=0
    do i=1,nk
      position=(frac(:,i)-frac(:,1))*mesh;order(:,i)=nint(position)
      if(maxval(abs(position-order(:,i)))>1d-8)code=1
      coord=modulo(order(:,i),mesh);a=1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3))
      if(lookup(a)/=0)code=1
      lookup(a)=i
    enddo
    call comm_get_max(code,icomm)
    if(code/=0)return
    do side=1,2
      sign=merge(1,-1,side==1)
      do j=1,nlocal
        coord=modulo(order(:,ks+j-1)+sign*qp,mesh)
        target_k(j,side)=lookup(1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3)))
      enddo
    enddo
    ! Stream only the occupied projectors needed for this rank's +/-q pairs.
    do i=1,nk
      if(rank==owners(i))buffer=psi(:,:,i-ks+1)
      call comm_bcast(buffer,icomm,owners(i))
      do side=1,2
        do j=1,nlocal
          if(target_k(j,side)==i)ghost(:,:,j,side)=buffer
        enddo
      enddo
    enddo
    points=0
    if(all(qp==0))then
      points(1,1)=1;points(2,2)=1;points(3,3)=1
    else
      points(1,2)=1;points(2,3)=1
    endif
    local=0d0;worst=0d0;scalar=1d0;source_c=1d0
    pairs: do side=1,2
      sign=merge(1,-1,side==1);shifts=z
      if(side==2)shifts=-conjg(z)
      q=sign*dble(qp)/mesh
      do j=1,nlocal
        source=ks+j-1;target=target_k(j,side);active_source=source;active_target=target
        fold=nint(frac(:,source)+q-frac(:,target))
        do i=0,ng-1
          position=dble([mod(i,n(1)),mod(i/n(1),n(2)),i/(n(1)*n(2))])/n
          do trial=1,3
            vertices(i+1,trial)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*(sign*points(:,trial)+fold)),8))
          enddo
        enddo
        do b=1,nv
          source_energy=energy(b,j)
          ! Independent vertices share H but not a combined Krylov basis.
          ! Keeping one RHS avoids tripling the factorization dimension.
          do trial=1,3
            rhs(:,1)=vertices(:,trial)*psi(:,b,j)
            rhs=rhs-matmul(ghost(:,:,j,side),matmul(conjg(transpose(ghost(:,:,j,side))),rhs))
            call frozen_solve(projected_action,energy_action,rhs,source_c,scalar,shifts,capacity,tol, &
              solution,statistics,i)
            if(i/=0)then
              write(*,'(a,6i7)')'EQUILIBRIUM_CHI0_REJECT rank/source/side/band/trial/status ', &
                rank,source,side,b,trial,i
              code=1;exit pairs
            endif
            if(comm_is_root(rank))write(*,'(a,5i7,es18.8)') &
              'EQUILIBRIUM_CHI0_BAND source/side/band/trial/space/residual ', &
              source,side,b,trial,statistics%left_size,maxval(solution%residual)
            actions=actions+statistics%h_actions;worst(1)=max(worst(1),maxval(solution%residual))
            do s=1,nz
              eigenvector=matmul(solution%left,solution%coefficient(:,1,s))*conjg(solution%right(1,1))
              value=dot_product(rhs(:,1),eigenvector)
              if(side==2)value=conjg(value)
              local(s,trial)=local(s,trial)+2d0*value/(volume*nk)
            enddo
          enddo
        enddo
      enddo
    enddo pairs
    call comm_get_max(code,icomm)
    if(code/=0)return
    call comm_summation(local,chi,size(local),icomm)
    call comm_get_max(worst,global_worst,1,icomm);residual=global_worst(1)
    if(.not.screen_finite(chi).or.residual>tol)return
    if(comm_is_root(rank))write(*,'(a,3i4,es20.10)')'EQUILIBRIUM_CHI0_PASS qp/residual ',qp,residual
    status=0
  contains
    subroutine projected_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      complex(8) :: projected(ng)
      projected=x-matmul(ghost(:,:,j,side),matmul(conjg(transpose(ghost(:,:,j,side))),x))
      call apply(active_source,active_target,projected,y,ierr)
      if(ierr/=0)return
      y=y-matmul(ghost(:,:,j,side),matmul(conjg(transpose(ghost(:,:,j,side))),y))
    end subroutine
    subroutine energy_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      y=source_energy*x;ierr=0
    end subroutine
  end subroutine
end module
