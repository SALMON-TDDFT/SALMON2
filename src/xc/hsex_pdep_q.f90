! Opt-in finite-q equilibrium diagnostic. Production W is never replaced here.
module hsex_pdep_q
  implicit none
  private
  public :: pdep_q_test
  abstract interface
    subroutine grid_action(k,x,y,status)
      implicit none
      integer,intent(in) :: k
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine pdep_q_test(n,mesh,qmesh,qp,cutoff,cell,k,levels,states,nv,ks,ke,comm,rank,prank,path,apply,status)
    use hsex_sternheimer, only: sternheimer_band,sternheimer_solve,s_sternheimer_stats
    use hsex_pdep, only: pdep_basis,pdep_inverse
    use hsex_screen_math, only: screen_geometry,screen_finite,screen_inverse
    use communication, only: comm_summation,comm_get_max,comm_bcast,comm_is_root,comm_sync_all
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    integer,intent(in) :: n(3),mesh(3),qmesh(3),qp(3),nv,ks,ke,comm,rank,prank
    real(8),intent(in) :: cutoff,cell(3,3),k(:,:),levels(:,:)
    complex(8),intent(in) :: states(:,:,:)
    character(*),intent(in) :: path
    procedure(grid_action) :: apply
    integer,intent(out) :: status
    interface
      function hsex_peak_rss_mib() bind(C) result(rss)
        import c_double
        real(c_double) :: rss
      end function
    end interface
    integer,parameter :: nt=6
    type(s_sternheimer_stats) :: statistics
    complex(8),allocatable :: all_states(:,:,:),local_states(:,:,:),base(:,:,:),vertex(:,:),fold(:),solution(:)
    complex(8),allocatable :: response(:,:),chi(:,:),total(:,:),basis(:,:),trials(:,:),actual(:,:),reference(:,:)
    complex(8),allocatable :: v(:),y(:),density(:),saved_trials(:,:),gram(:,:)
    real(8),allocatable :: frac(:,:),b(:),eigen(:),ritz(:)
    integer,allocatable :: g(:,:,:),lookup(:),order(:,:)
    integer :: nk,ng,nc,imax(3),a,c,d,ig,j,source,orientation,band,shift(3),coord(3),idx,code,target,calls
    integer :: u,ios,close_code,header(14),saved_header(14),pair,ia,ib,info,hole_band
    integer,allocatable :: saved_g(:,:)
    real(8) :: rec(3,3),vol,q(3,2),r(3),saved_cell(3,3),saved_cutoff,error(nt),rss(2),counts(3),count_sum(3)
    real(8) :: timing(4),timing_max(4),worst(2),worst_max(2),rss_max(2),phases(4),phase_max(4),ksum_seconds
    integer(8) :: started,prepared,built,finished,rate,collective_start,collective_end
    complex(8) :: elements(2),expected(2)
    character(4) :: magic
    call comm_sync_all(comm);call system_clock(started,rate)
    status=1;nk=product(mesh);ng=product(n);calls=0;statistics=s_sternheimer_stats();ksum_seconds=0d0
    code=0
    if(any(qmesh<1))code=1
    if(nv<1.or.size(states,1)/=ng.or.size(states,2)<nv+1.or.size(states,3)/=max(0,ke-ks+1))code=1
    if(.not.screen_finite(states).or..not.screen_finite(levels))code=1
    if(prank==0.or.prank< -1)code=1
    call comm_get_max(code,comm)
    if(code/=0)return
    if(any(mod(mesh,qmesh)/=0).or.all(modulo(qp,qmesh)==0))return
    call screen_geometry(cell,rec,vol,code)
    if(code/=0)return
    q(:,1)=dble(modulo(qp,qmesh))/qmesh;q(:,2)=dble(modulo(-qp,qmesh))/qmesh
    imax=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell**2,dim=1))/(2d0*acos(-1d0))+1d0)
    nc=0
    do a=-imax(1),imax(1)
      do c=-imax(2),imax(2)
        do d=-imax(3),imax(3)
          if(sum(matmul(rec,dble([a,c,d])+q(:,1))**2)/2d0<=cutoff)nc=nc+1
        enddo
      enddo
    enddo
    if(nc<1.or.prank>nc)return
    allocate(g(3,nc,2),b(nc));ig=0
    do a=-imax(1),imax(1)
      do c=-imax(2),imax(2)
        do d=-imax(3),imax(3)
          if(sum(matmul(rec,dble([a,c,d])+q(:,1))**2)/2d0>cutoff)cycle
          ig=ig+1;g(:,ig,1)=[a,c,d]
          g(:,ig,2)=-g(:,ig,1)-nint(q(:,1)+q(:,2))
          b(ig)=sqrt(4d0*acos(-1d0)/sum(matmul(rec,g(:,ig,1)+q(:,1))**2))
        enddo
      enddo
    enddo
    allocate(frac(3,nk),lookup(nk),order(3,nk))
    frac=matmul(transpose(cell),k)/(2d0*acos(-1d0));lookup=0
    do j=1,nk
      r=(frac(:,j)-frac(:,1))*mesh;order(:,j)=nint(r)
      if(maxval(abs(r-order(:,j)))>1d-8)return
      coord=modulo(order(:,j),mesh);idx=1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3))
      if(lookup(idx)/=0)return
      lookup(idx)=j
    enddo
    allocate(local_states(ng,nv+1,nk),all_states(ng,nv+1,nk));local_states=0d0
    if(ke>=ks)local_states(:,:,ks:ke)=states(:,:nv+1,:)
    call comm_summation(local_states,all_states,size(all_states),comm)
    deallocate(local_states)
    if(minval(levels(nv+1,:))-maxval(levels(nv,:))<=1d-8)return
    allocate(gram(nv,nv))
    do j=1,nk
      gram=matmul(conjg(transpose(all_states(:,:nv,j))),all_states(:,:nv,j))
      do a=1,nv
        gram(a,a)=gram(a,a)-1d0
      enddo
      if(.not.screen_finite(gram))return
      if(maxval(abs(gram))>1d-8)return
    enddo
    deallocate(gram)
    allocate(base(ng,nc,2),vertex(ng,nc),fold(ng),solution(ng),v(nc),y(nc),density(ng))
    do orientation=1,2
      do ig=1,nc
        do a=0,ng-1
          r=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
          base(a+1,ig,orientation)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(r*g(:,ig,orientation)),8))
        enddo
      enddo
      do j=1,nk
        call source_vertex(j,orientation,source,shift)
        do ig=1,nc
          if(any(2*abs(g(:,ig,orientation)+shift)>=n))return
        enddo
      enddo
    enddo
    ! Validate occupied eigenpairs with the actual real-space Hamiltonian.
    code=0
    do target=ks,ke
      do band=1,nv
        call apply(target,all_states(:,band,target),solution,info)
        if(info/=0)code=1
        if(.not.screen_finite(solution))code=1
        if(norm2(abs(solution-levels(band,target)*all_states(:,band,target)))>1d-7)code=1
      enddo
    enddo
    call comm_get_max(code,comm)
    if(code/=0)return
    allocate(trials(nc,nt),actual(nc,nt),reference(nc,nt));trials=0d0
    do a=1,nc
      trials(a,1)=cmplx(sin(dble(a)),cos(2d0*a),8)
      trials(a,2)=cmplx(cos(3d0*a),sin(2d0*a),8)
    enddo
    do a=1,2
      trials(:,a)=trials(:,a)/norm2(abs(trials(:,a)))
    enddo
    do pair=1,2
      j=1+(pair-1)*(nk/2)
      call make_vertex(j,1,nv+1)
      trials(:,2*pair+1)=matmul(conjg(transpose(vertex)),all_states(:,nv+1,j))
      hole_band=max(1,nv-pair+1)
      call make_vertex(j,1,hole_band)
      trials(:,2*pair+2)=matmul(conjg(transpose(vertex)),all_states(:,hole_band,j))
    enddo
    header=[n,mesh,qmesh,qp,nc,nt]
    rss(1)=hsex_peak_rss_mib()
    code=merge(1,0,rss(1)<0d0);call comm_get_max(code,comm)
    if(code/=0)return
    call comm_sync_all(comm);call system_clock(prepared)
    if(prank<0)then
      allocate(chi(nc,nc),response(nc,nc),total(nc,nc));chi=0d0
      code=0
      reference_sides: do orientation=1,2
        do target=ks,ke
          call source_vertex(target,orientation,source,shift)
          do band=1,nv
            call make_vertex(target,orientation,band)
            call sternheimer_band(at_target,all_states(:,:nv,target),levels(:nv,target), &
              spread(1d0,1,nv),levels(band,source),1d0,1d0,vertex,response,code,statistics)
            if(code/=0)exit reference_sides
            if(orientation==1)then
              chi=chi+2d0*response
            else
              chi=chi+2d0*conjg(response)
            endif
          enddo
        enddo
      enddo reference_sides
      call system_clock(collective_start)
      call comm_get_max(code,comm)
      if(code/=0)return
      call comm_summation(chi,total,size(total),comm)
      call system_clock(collective_end)
      ksum_seconds=ksum_seconds+dble(collective_end-collective_start)/rate
      do a=1,nc
        total(:,a)=2d0/(nk*vol)*b*total(:,a)*b(a)
        total(a,a)=total(a,a)+1d0
      enddo
      call system_clock(built)
      call screen_inverse(total,chi,code)
      call comm_get_max(code,comm)
      if(code/=0)return
      do a=1,nt
        actual(:,a)=b*matmul(chi,b*trials(:,a))
      enddo
      reference=actual
      ios=0
      if(comm_is_root(rank))then
        open(newunit=u,file=path,status='new',access='stream',form='unformatted',iostat=ios)
        if(ios==0)then
          write(u,iostat=ios) 'PDQ1',header,g(:,:,1),cell,cutoff,trials,reference
          close(u,iostat=close_code)
          if(ios==0)ios=close_code
        endif
      endif
      call comm_bcast(ios,comm,0)
      if(ios/=0)return
      calls=nc
    else
      allocate(saved_g(3,nc),saved_trials(nc,nt),basis(nc,prank),eigen(prank),ritz(prank))
      ios=0
      if(comm_is_root(rank))then
        open(newunit=u,file=path,status='old',access='stream',form='unformatted',iostat=ios)
        if(ios==0)then
          read(u,iostat=ios)magic,saved_header,saved_g,saved_cell,saved_cutoff,saved_trials,reference
          close(u,iostat=close_code)
          if(ios==0)ios=close_code
          if(ios==0)then
            if(magic/='PDQ1'.or.any(saved_header/=header).or.any(saved_g/=g(:,:,1)))ios=1
            if(any(saved_cell/=cell).or.saved_cutoff/=cutoff)ios=1
            if(.not.screen_finite(reference).or..not.screen_finite(saved_trials))ios=1
            if(maxval(abs(saved_trials-trials))>1d-10)ios=1
          endif
        endif
      endif
      call comm_bcast(ios,comm,0)
      if(ios/=0)return
      call comm_bcast(reference,comm,0)
      call pdep_basis(polarization,nc,prank,basis,eigen,ritz,calls,code)
      if(code/=0)return
      call system_clock(built)
      do a=1,nt
        call pdep_inverse(basis,eigen,b*trials(:,a),y,code)
        call comm_get_max(code,comm)
        if(code/=0)return
        actual(:,a)=b*y
      enddo
    endif
    call comm_sync_all(comm);call system_clock(finished)
    timing=dble([prepared-started,built-prepared,finished-built,finished-started])/rate
    call comm_get_max(timing,timing_max,4,comm)
    phases=[statistics%solve_seconds,ksum_seconds, &
      max(0d0,timing(2)-statistics%solve_seconds-ksum_seconds),timing(3)]
    call comm_get_max(phases,phase_max,4,comm)
    counts=dble([statistics%rhs,statistics%iterations,statistics%h_actions])
    call comm_summation(counts,count_sum,3,comm)
    rss(2)=hsex_peak_rss_mib();call comm_get_max(rss,rss_max,2,comm)
    worst=[statistics%max_relative_residual,0d0]
    if(prank>0)worst(2)=maxval(ritz)
    call comm_get_max(worst,worst_max,2,comm)
    do a=1,nt
      error(a)=norm2(abs(actual(:,a)-reference(:,a)))/max(norm2(abs(reference(:,a))),1d-14)
    enddo
    do pair=1,2
      ia=2*pair+1;ib=ia+1
      elements(pair)=dot_product(trials(:,ia),actual(:,ib))/vol
      expected(pair)=dot_product(trials(:,ia),reference(:,ib))/vol
    enddo
    if(comm_is_root(rank))then
      write(*,'(a,4i8)')'PDEP_DIMENSIONS grid/k/occupied/G ',ng,nk,nv,nc
      write(*,'(a,3f14.8)')'PDEP_Q_FRACTIONAL ',q(:,1)
      write(*,'(a,4es18.8)')'PDEP_Q_CARTESIAN_AU xyz/norm ',matmul(rec,q(:,1)),norm2(matmul(rec,q(:,1)))
      write(*,'(a)')'PDEP_CONVENTIONS full uniform k; spin=2; +/-q explicit; finite-q only; Gamma unchanged'
      write(*,'(a,i0,a,i0,a)')'PDEP_RANK ',prank,' polarization_actions=',calls,' explicit_empty_sum=0'
      write(*,'(a,4es18.8)')'PDEP_SECONDS prep/build/inverse+checks/total ',timing_max
      write(*,'(a,4es18.8)')'PDEP_PHASE_SECONDS CG/k_reduce/vertices+assembly+Ritz/inverse ',phase_max
      write(*,'(a,3f20.0)')'PDEP_SOLVES rhs/iterations/H_actions ',count_sum
      write(*,'(a,2es18.8)')'PDEP_RESIDUAL solve/ritz ',worst_max
      if(prank>0)write(*,'(a,l1)')'PDEP_RITZ_CONVERGED ',worst_max(2)<=1d-7
      write(*,'(a,2f16.5)')'PDEP_PEAK_RSS_MIB before/after max_rank ',rss_max
      write(*,'(a,6es18.8)')'PDEP_W_ACTION_REL_ERRORS ',error
      do pair=1,2
        write(*,'(a,i0,4es20.10)')'PDEP_BSE_DIRECT_HA pair=',pair,elements(pair),expected(pair)
      enddo
      write(*,'(a)')'PDEP_Q_TEST_PASS'
    endif
    status=0
  contains
    subroutine source_vertex(destination,side,source_index,wrap)
      implicit none
      integer,intent(in) :: destination,side
      integer,intent(out) :: source_index,wrap(3)
      integer :: ix(3),slot
      ix=modulo(order(:,destination)-nint(q(:,side)*mesh),mesh)
      slot=1+ix(1)+mesh(1)*(ix(2)+mesh(2)*ix(3));source_index=lookup(slot)
      wrap=nint(frac(:,source_index)+q(:,side)-frac(:,destination))
    end subroutine
    subroutine make_vertex(destination,side,source_band)
      implicit none
      integer,intent(in) :: destination,side,source_band
      integer :: s,wrap(3),site,component
      real(8) :: position(3)
      call source_vertex(destination,side,s,wrap)
      do site=0,ng-1
        position=dble([mod(site,n(1)),mod(site/n(1),n(2)),site/(n(1)*n(2))])/n
        fold(site+1)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*wrap),8))
      enddo
      density=fold*all_states(:,source_band,s)
      do component=1,nc
        vertex(:,component)=base(:,component,side)*density
      enddo
    end subroutine
    subroutine at_target(x,y,code)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: code
      call apply(target,x,y,code)
    end subroutine
    subroutine polarization(x,y,code)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: code
      complex(8) :: local(nc),global(nc),weights(nc),action(nc),rhs(ng)
      integer :: side,source_k,wrap(3),source_band
      local=0d0;code=0
      action_sides: do side=1,2
        weights=b*x
        if(side==2)weights=conjg(weights)
        do target=ks,ke
          call source_vertex(target,side,source_k,wrap)
          do source_band=1,nv
            call make_vertex(target,side,source_band)
            rhs=matmul(vertex,weights)
            call sternheimer_solve(at_target,all_states(:,:nv,target),levels(source_band,source_k), &
              rhs,solution,code,statistics)
            if(code/=0)exit action_sides
            action=matmul(conjg(transpose(vertex)),solution)
            if(side==2)action=conjg(action)
            local=local+action
          enddo
        enddo
      enddo action_sides
      call system_clock(collective_start)
      call comm_get_max(code,comm)
      if(code/=0)then
        y=0d0
        return
      endif
      call comm_summation(local,global,nc,comm)
      call system_clock(collective_end)
      ksum_seconds=ksum_seconds+dble(collective_end-collective_start)/rate
      y=2d0/(nk*vol)*b*global;code=0
    end subroutine
  end subroutine
end module
