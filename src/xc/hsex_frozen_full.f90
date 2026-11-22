! Experimental coherent occupied-orbital, adiabatic static snapshot RPA.
module hsex_frozen_full
  implicit none
  private
  public :: frozen_full_w,s_frozen_window_cache,frozen_window_reset,frozen_window_anchor
  type :: q_response
    integer,allocatable :: g(:,:),reverse(:)
    complex(8),allocatable :: chi(:,:,:)
    real(8),allocatable :: error(:)
    integer :: minus
  end type
  type :: s_frozen_window_cache
    private
    logical :: ready=.false.,anchored=.false.,density_limit=.false.
    integer :: n(3)=0,mesh(3)=0,boxes=0,nv=0,first=0,nlocal=0,kind=0
    real(8) :: cutoff=0d0,eta=0d0,time=0d0,cell(3,3)=0d0,closure_tol=0d0,density_tol=0d0
    real(8),allocatable :: k(:,:),reference_w(:,:),anchor_w(:,:)
    type(q_response),allocatable :: blocks(:)
  end type
  abstract interface
    subroutine grid_action(k,x,y,status)
      implicit none
      integer,intent(in) :: k
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: status
    end subroutine
    subroutine prepare_pair(source,target,worker,source_owner,target_owner,status)
      implicit none
      integer,intent(in) :: source,target,worker,source_owner,target_owner
      integer,intent(out) :: status
    end subroutine
    subroutine grid_velocity(k,x,y,status)
      implicit none
      integer,intent(in) :: k
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine frozen_window_reset(cache)
    implicit none
    type(s_frozen_window_cache),intent(inout) :: cache
    cache=s_frozen_window_cache()
  end subroutine
  subroutine frozen_window_anchor(cache,w0,tol,icomm,status)
    use communication, only: comm_get_max,comm_summation
    use hsex_screen_math, only: screen_finite
    implicit none
    type(s_frozen_window_cache),intent(inout) :: cache
    real(8),intent(in) :: w0(:,:),tol
    integer,intent(in) :: icomm
    integer,intent(out) :: status
    real(8) :: local(2),total(2),error
    integer :: code
    status=1;code=0
    if(.not.cache%ready.or..not.allocated(cache%reference_w))code=1
    if(.not.screen_finite(w0).or..not.screen_finite([tol]))code=1
    if(tol<=0d0.or.tol>=1d0)code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(any(shape(w0)/=shape(cache%reference_w)))code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    local=[sum((w0-cache%reference_w)**2),sum(cache%reference_w**2)]
    call comm_summation(local,total,2,icomm)
    if(.not.screen_finite(total))return
    error=sqrt(total(1)/max(tiny(1d0),total(2)))
    if(.not.screen_finite([error]))return
    if(error>tol)return
    cache%anchor_w=w0;cache%anchored=.true.;status=0
  end subroutine
  subroutine frozen_full_w(n,mesh,qmesh,boxes,cutoff,cell,k,states,nv,ks,ke,icomm,rank,np, &
    first,nlocal,capacity,tol,eta,closure_tol,apply,w,status,velocity,density_gamma,gamma_tol,prepare,nproc_q, &
    precondition,closure,window_kind,window_time,window_cache,initialize_window)
    use communication, only: comm_bcast,comm_summation,comm_get_max,comm_is_root
    use hsex_screen_math, only: screen_geometry,screen_finite,screen_sites,screen_assemble,screen_transform,screen_channels
    use hsex_frozen_response, only: frozen_sources,frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats,frozen_background_rate
    use hsex_frozen_limit, only: frozen_optical_tensor,frozen_eta_limit,frozen_static_part
    use hsex_density_gamma, only: density_gamma_stencil,density_gamma_replace
    use hsex_longwave, only: longwave_cell
    use, intrinsic :: iso_fortran_env, only: output_unit
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    integer,intent(in) :: n(3),mesh(3),qmesh(3),boxes,nv,ks,ke,icomm,rank,np,first,nlocal,capacity
    real(8),intent(in) :: cutoff,cell(3,3),k(:,:),tol,eta,closure_tol
    complex(8),intent(in) :: states(:,:,:)
    procedure(grid_action) :: apply
    procedure(grid_velocity),optional :: velocity
    procedure(prepare_pair),optional :: prepare
    procedure(grid_action),optional :: precondition
    integer,optional,intent(in) :: nproc_q
    character(*),optional,intent(in) :: closure
    integer,optional,intent(in) :: window_kind
    real(8),optional,intent(in) :: window_time
    type(s_frozen_window_cache),optional,intent(inout) :: window_cache
    logical,optional,intent(in) :: initialize_window
    integer :: filter,static_samples
    real(8) :: area,delta_norm,chi_error
    logical :: window_active,window_init
    integer,allocatable :: filter_kind(:)
    real(8),allocatable :: filter_time(:)
    character(16) :: static_closure
    logical :: project_static
    interface
      function hsex_peak_rss_mib() bind(C) result(rss)
        import c_double
        implicit none
        real(c_double) :: rss
      end function
    end interface
    integer :: total_precondition
    real(8) :: peak_rss(1),max_rss(1)
    integer :: qworkers
    logical,optional,intent(in) :: density_gamma
    real(8),optional,intent(in) :: gamma_tol
    real(8) :: density_tolerance
    logical :: density_limit
    integer,allocatable :: density_samples(:)
    real(8),intent(out) :: w(nlocal,product(mesh))
    integer,intent(out) :: status
    type(q_response),allocatable :: blocks(:)
    type(s_frozen_solution) :: solution
    type(s_frozen_stats) :: stats
    integer :: ng,nk,nq,na,code,iq,jq,mg,a,b,c,ig,jg,j,source,target,idx,unit_shift(3),coord(3),limits(3)
    integer :: command(2),request_code,owner,shift(3),ios,total_actions,basis(2),ia,ib,neta,nb
    integer,allocatable :: owners(:),ranges(:,:),local_ranges(:,:),lookup(:),order(:,:),sites(:,:)
    complex(8),allocatable :: buffer(:),answer(:),orbitals(:,:),ps(:,:),pt(:,:),vertices(:,:),u(:,:),v(:,:),coeff(:,:,:)
    complex(8),allocatable :: flat(:),raw(:,:),box_value(:,:,:),coarse(:,:),local_matrix(:,:),current(:,:),partner(:,:)
    complex(8),allocatable :: gram(:,:)
    complex(8),allocatable :: velocity_answer(:,:),optical(:,:),part(:,:),optical_vertices(:,:),point(:,:),origin(:,:)
    complex(8),allocatable :: extrapolated(:,:),eta_chi(:,:,:),static_w(:,:)
    integer,allocatable :: optical_g(:,:),optical_reverse(:)
    real(8),allocatable :: frac(:,:),qcoords(:,:),heads(:),qheads(:)
    real(8) :: rec(3,3),vol,q(3),position(3),worst,sym_error,anti_error,error,einv,qmin,qc,seconds
    real(8) :: macro(3,3),background_rate
    complex(8),allocatable :: shifts(:)
    integer(8) :: started,finished,rate
    filter=0;area=0d0;window_init=.false.
    if(present(window_kind))filter=window_kind
    if(present(window_time))area=window_time
    if(present(initialize_window))window_init=initialize_window
    window_active=filter/=0
    static_closure='strict'
    if(present(closure))static_closure=closure
    project_static=static_closure=='hermitian'
    qworkers=1
    if(present(nproc_q))qworkers=nproc_q
    density_tolerance=closure_tol
    if(present(gamma_tol))density_tolerance=gamma_tol
    density_limit=.false.
    if(present(density_gamma))density_limit=density_gamma
    status=10;w=0d0;ng=product(n);nk=product(mesh);nq=product(qmesh);na=boxes**3;code=0
    if(any(n<1).or.any(mesh<1).or.any(qmesh<1).or.boxes<1.or.capacity<1.or.nv<1)code=1
    if(.not.screen_finite([cutoff,tol,eta,closure_tol]).or..not.screen_finite(cell))code=1
    if(cutoff<=0d0.or.tol<=0d0.or.eta<=0d0.or.closure_tol<=0d0.or.closure_tol>=1d0)code=1
    if(size(k,1)/=3.or.size(k,2)/=nk.or.size(states,1)/=ng.or.size(states,2)/=nv)code=1
    if(size(states,3)/=max(0,ke-ks+1).or..not.screen_finite(states).or..not.screen_finite(k))code=1
    if(first<1.or.nlocal<0.or.rank<0.or.rank>=np)code=1
    if(nlocal>0.and.first+nlocal-1>na*na)code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(window_active)then
      call validate_window(code)
      call comm_get_max(code,icomm)
      if(code/=0)return
    else
      if(window_init)return
    endif
    if(static_closure/='strict'.and..not.project_static)return
    if(project_static.and.present(velocity))return
    if(qworkers<1.or.qworkers>np)return
    if(mod(np,qworkers)/=0)return
    if(present(prepare).and.present(velocity))return
    if(present(precondition).and..not.present(prepare))return
    if(any(qmesh/=mesh).or.nq<2.or.any(mod(n,boxes)/=0))return
    do a=1,3
      j=mesh(a)
      do while(j>1.and.mod(j,2)==0)
        j=j/2
      enddo
      if(j/=1)return
    enddo
    call screen_geometry(cell,rec,vol,code)
    if(code/=0)return
    if(density_limit)then
      status=11
      if(present(velocity))return
      if(.not.screen_finite([density_tolerance]))return
      if(density_tolerance<=0d0.or.density_tolerance>=1d0)return
      call density_gamma_stencil(qmesh,rec,density_samples,code)
      call comm_get_max(code,icomm)
      if(code/=0)return
      deallocate(density_samples)
    endif
    allocate(owners(nk),ranges(2,np),local_ranges(2,np),lookup(nk),order(3,nk),sites(3,nq),frac(3,nk), &
      buffer(ng),answer(ng),orbitals(ng,nv),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    local_ranges=0;local_ranges(:,rank+1)=[ks,ke]
    call comm_summation(local_ranges,ranges,size(ranges),icomm)
    owners=-1
    do j=1,np
      do a=ranges(1,j),ranges(2,j)
        if(a<1.or.a>nk)then
          code=1;cycle
        endif
        if(owners(a)/=-1)code=1
        owners(a)=j-1
      enddo
    enddo
    if(any(owners<0))code=1
    do j=1,size(states,3)
      gram=matmul(conjg(transpose(states(:,:,j))),states(:,:,j))
      do a=1,nv
        gram(a,a)=gram(a,a)-1d0
      enddo
      if(maxval(abs(gram))>1d-8)code=1
    enddo
    frac=matmul(transpose(cell),k)/(2d0*acos(-1d0));lookup=0
    do j=1,nk
      position=(frac(:,j)-frac(:,1))*mesh;order(:,j)=nint(position)
      if(maxval(abs(position-order(:,j)))>1d-8)code=1
      coord=modulo(order(:,j),mesh);idx=1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3))
      if(lookup(idx)/=0)code=1
      lookup(idx)=j
    enddo
    call comm_get_max(code,icomm)
    if(code/=0)return
    call screen_sites(qmesh,sites)
    neta=2
    if(present(velocity).or.density_limit.or.project_static)neta=3
    if(window_active)then
      neta=1
      if(window_init)neta=4
    endif
    static_samples=neta
    if(window_init)static_samples=3
    allocate(shifts(neta),filter_kind(neta),filter_time(neta),velocity_answer(ng,3),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    filter_kind=0;filter_time=0d0
    do c=1,neta
      shifts(c)=cmplx(0d0,eta/2d0**(c-1),8)
    enddo
    if(window_active)then
      filter_kind(neta)=filter;filter_time(neta)=area;shifts(neta)=cmplx(0d0,1d0/area,8)
      if(window_init.and.comm_is_root(rank))then
        allocate(window_cache%blocks(nq),stat=code)
      endif
      call comm_get_max(code,icomm)
      if(code/=0)return
      if(comm_is_root(rank))write(*,'(a,i3,es20.10,l2)') &
        'HSEX_FROZEN_WINDOW kind/area_au/reference ',filter,area,window_init
    endif
    call system_clock(started,rate)
    if(comm_is_root(rank))write(*,'(a,a)')'HSEX_FROZEN_CLOSURE ',trim(static_closure)
    status=0;worst=0d0;total_actions=0;total_precondition=0;basis=0
    if(present(prepare))then
      call parallel_response()
    else
    if(comm_is_root(rank))then
      allocate(blocks(nq),ps(ng,nv),pt(ng,nv),stat=ios)
      if(ios/=0)status=20
      limits=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell**2,dim=1))/(2d0*acos(-1d0))+1d0)
      qloop: do iq=1,nq
        if(status/=0)exit qloop
        if(density_limit.and.iq==1)cycle qloop
        q=dble(sites(:,iq))/qmesh;mg=0
        do a=-limits(1),limits(1)
          do b=-limits(2),limits(2)
            do c=-limits(3),limits(3)
              if(sum(matmul(rec,[a,b,c]+q)**2)/2d0<=cutoff)mg=mg+1
            enddo
          enddo
        enddo
        if(mg<1)then
          status=20;exit qloop
        endif
        allocate(blocks(iq)%g(3,mg),blocks(iq)%chi(mg,mg,neta),blocks(iq)%reverse(mg), &
          vertices(ng,mg),flat(neta*mg),stat=ios)
        if(ios/=0)then
          status=20;exit qloop
        endif
        j=0;blocks(iq)%chi=0d0
        do a=-limits(1),limits(1)
          do b=-limits(2),limits(2)
            do c=-limits(3),limits(3)
              if(sum(matmul(rec,[a,b,c]+q)**2)/2d0>cutoff)cycle
              j=j+1;blocks(iq)%g(:,j)=[a,b,c]
            enddo
          enddo
        enddo
        if(iq==1.and.present(velocity))then
          nb=mg-1
          allocate(optical_g(3,nb),optical_reverse(nb+3),optical(nb+3,nb+3),part(nb+3,nb+3), &
            optical_vertices(ng,nb),point(na,na),origin(na,na),stat=ios)
          if(ios/=0)then
            status=20;exit qloop
          endif
          j=0
          do ig=1,mg
            if(all(blocks(iq)%g(:,ig)==0))cycle
            j=j+1;optical_g(:,j)=blocks(iq)%g(:,ig)
          enddo
          optical=0d0
          do ig=1,nb
            do a=0,ng-1
              position=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
              optical_vertices(a+1,ig)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*optical_g(:,ig)),8))
            enddo
          enddo
          do source=1,nk
            command=[2,source];call request();ps=orbitals
            call frozen_optical_tensor(right_action,velocity_action,ps,optical_vertices, &
              max(1d-8,tol),part,background_rate,code)
            write(*,'(a,2i6,es18.8)')'HSEX_FROZEN_OPTICAL_STATIONARITY k/status/rate ',source,code,background_rate
            if(code/=0)then
              status=50;exit qloop
            endif
            optical=optical+part
          enddo
          optical_reverse(:3)=[1,2,3]
          do ig=1,nb
            optical_reverse(ig+3)=0
            do jg=1,nb
              if(all(optical_g(:,jg)==-optical_g(:,ig)))optical_reverse(ig+3)=jg+3
            enddo
          enddo
          if(any(optical_reverse==0))then
            status=51;exit qloop
          endif
          do b=1,nb+3
            do a=1,nb+3
              part(a,b)=optical(a,b)+merge(-1d0,1d0,(a<=3).neqv.(b<=3))* &
                conjg(optical(optical_reverse(a),optical_reverse(b)))
            enddo
          enddo
          call longwave_cell(n,boxes,optical_g,rec,vol,nk,part,macro,point,origin,code)
          if(code/=0)then
            status=52;exit qloop
          endif
          write(*,'(a,3es18.8)')'HSEX_FROZEN_EPS_LONGWAVE diagonal ',(macro(a,a),a=1,3)
        endif
        coord=modulo(-sites(:,iq),qmesh)
        blocks(iq)%minus=1+coord(1)+qmesh(1)*(coord(2)+qmesh(2)*coord(3))
        do source=1,nk
          coord=modulo(order(:,source)+sites(:,iq),mesh)
          target=lookup(1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3)))
          shift=nint(frac(:,source)+q-frac(:,target))
          command=[2,source];call request();ps=orbitals
          command=[2,target];call request();pt=orbitals
          do ig=1,mg
            unit_shift=blocks(iq)%g(:,ig)+shift
            if(any(2*abs(unit_shift)>=n))then
              status=21;exit qloop
            endif
            do a=0,ng-1
              position=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
              vertices(a+1,ig)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*unit_shift),8))
            enddo
          enddo
          call frozen_sources(ps,pt,2d0*[(1d0,a=1,nv)],2d0*[(1d0,a=1,nv)],vertices,u,v,coeff,code)
          if(code/=0)then
            status=22;exit qloop
          endif
          call frozen_solve(left_action,right_action,u,coeff(:,:,1),v,shifts,capacity,tol,solution,stats,code, &
            extra_sources=coeff(:,:,2:),residual_bound=.true.,window_kind=filter_kind,window_time=filter_time)
          if(code/=0)then
            write(*,'(a,4i8)')'HSEX_FROZEN_FULL_SOLVE_REJECT q/k/code/cap ',iq,source,code,capacity
            status=23;exit qloop
          endif
          worst=max(worst,maxval(solution%residual));total_actions=total_actions+stats%h_actions
          basis=max(basis,[stats%left_size,stats%right_size])
          do ig=1,mg
            call frozen_trace(solution,vertices(:,ig),flat,code)
            if(code/=0)then
              status=24;exit qloop
            endif
            do j=1,mg
              blocks(iq)%chi(ig,j,:)=blocks(iq)%chi(ig,j,:)+flat(neta*(j-1)+1:neta*j)
            enddo
          enddo
          write(*,'(a,4i7,es16.7)')'HSEX_FROZEN_FULL_PAIR q/k/G/actions/residual ',iq,source,mg,stats%h_actions,worst
          flush(output_unit)
        enddo
        deallocate(vertices,flat)
      enddo qloop
      command=0;call comm_bcast(command,icomm,0)
    else
      do
        call comm_bcast(command,icomm,0)
        if(command(1)==0)exit
        call service()
      enddo
    endif
    endif
    call comm_bcast(status,icomm,0)
    if(status/=0)return
    allocate(coarse(nlocal,nq),box_value(na,na,neta),heads(nq),qheads(nq),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)then
      status=30;return
    endif
    heads=0d0;qheads=0d0;coarse=0d0
    do iq=1,nq
      if(density_limit.and.iq==1)cycle
      if(comm_is_root(rank))then
        jq=blocks(iq)%minus;mg=size(blocks(iq)%g,2)
        allocate(raw(mg,mg),qcoords(3,mg),eta_chi(mg,mg,neta),stat=ios)
        if(ios/=0)status=30
        if(status==0)then
          do ig=1,mg
            coord=-blocks(iq)%g(:,ig)-(sites(:,iq)+sites(:,jq))/qmesh
            blocks(iq)%reverse(ig)=0
            do jg=1,size(blocks(jq)%g,2)
              if(all(coord==blocks(jq)%g(:,jg)))blocks(iq)%reverse(ig)=jg
            enddo
          enddo
          if(any(blocks(iq)%reverse==0))status=31
        endif
        if(status==0)then
          qcoords=blocks(iq)%g+spread(dble(sites(:,iq))/qmesh,2,mg)
          if(window_active)then
            if(window_init)then
              window_cache%blocks(iq)%g=blocks(iq)%g
              allocate(window_cache%blocks(iq)%chi(mg,mg,2),window_cache%blocks(iq)%error(2),stat=ios)
              if(ios/=0)status=30
            else
              if(.not.allocated(window_cache%blocks(iq)%g))status=30
              if(status==0)then
                if(any(shape(window_cache%blocks(iq)%g)/=shape(blocks(iq)%g)))status=30
              endif
              if(status==0)then
                if(any(window_cache%blocks(iq)%g/=blocks(iq)%g))status=30
              endif
            endif
          endif
          if(status/=0)then
            deallocate(raw,qcoords,eta_chi)
          else
          do c=1,neta
            do b=1,mg
              do a=1,mg
                raw(a,b)=conjg(blocks(jq)%chi(blocks(iq)%reverse(a),blocks(iq)%reverse(b),c))
              enddo
            enddo
            sym_error=sqrt(sum(abs(raw-blocks(iq)%chi(:,:,c))**2))/ &
              max(1d-12,sqrt(sum(abs(blocks(iq)%chi(:,:,c))**2)))
            if(.not.screen_finite([sym_error]))then
              status=32;exit
            endif
            raw=.5d0*(raw+blocks(iq)%chi(:,:,c))
            call frozen_static_part(raw,static_closure,closure_tol,eta_chi(:,:,c),anti_error,code)
            if(window_active.and.filter_kind(c)/=0)then
              write(*,'(a,2i6,2es18.8)')'HSEX_FROZEN_WINDOW_DEFECT q/channel/reality/antihermitian ', &
                iq,c,sym_error,anti_error
            else
              write(*,'(a,2i6,2es18.8)')'HSEX_FROZEN_STATIC_DEFECT q/eta/reality/antihermitian ',iq,c,sym_error,anti_error
            endif
            if(sym_error>closure_tol.or.code/=0)then
              status=32;exit
            endif
            if(project_static)then
              if(window_active.and.filter_kind(c)/=0)then
                write(*,'(a,2i6,es18.8)')'HSEX_FROZEN_WINDOW_PROJECTION q/channel/discarded_relative ',iq,c,.5d0*anti_error
              else
                write(*,'(a,2i6,es18.8)')'HSEX_FROZEN_STATIC_PROJECTION q/eta/discarded_relative ',iq,c,.5d0*anti_error
              endif
            endif
            raw=eta_chi(:,:,c)
            if(window_active)then
              if(window_init.and.c==neta)then
                window_cache%blocks(iq)%chi(:,:,2)=raw
                window_cache%blocks(iq)%error(2)=.5d0*(blocks(iq)%error(c)+blocks(jq)%error(c))
                cycle
              else if(.not.window_init)then
                delta_norm=sqrt(sum(abs(raw-window_cache%blocks(iq)%chi(:,:,2))**2))
                write(*,'(a,i6,es20.10)')'HSEX_FROZEN_WINDOW_DELTA_CHI q/norm ',iq,delta_norm
                raw=window_cache%blocks(iq)%chi(:,:,1)+(raw-window_cache%blocks(iq)%chi(:,:,2))
                chi_error=.5d0*(blocks(iq)%error(c)+blocks(jq)%error(c))+sum(window_cache%blocks(iq)%error)
                if(.not.screen_finite([chi_error]))then
                  status=32;exit
                endif
                write(*,'(a,i6,es20.10)')'HSEX_FROZEN_WINDOW_CHI_ERROR_BOUND q/absolute_Frobenius ',iq,chi_error
              endif
            endif
            ! Legacy assembler uses +2 chi/(Nk*volume); this chi includes spin and both transitions.
            call screen_assemble(n,boxes,qcoords,rec,vol,nk,-.5d0*raw,.true.,box_value(:,:,c), &
              heads(iq),qheads(iq),code)
            if(code/=0)then
              status=33;exit
            endif
          enddo
          if(status==0)then
            if(window_active.and..not.window_init)then
              ! Width is a physical finite response window, not an eta->0 regulator.
              error=0d0
            else if(present(velocity).or.density_limit.or.project_static)then
              allocate(extrapolated(mg,mg),static_w(na,na))
              call frozen_eta_limit(eta_chi(:,:,:static_samples),closure_tol,extrapolated,error,code)
              write(*,'(a,i6,es18.8)')'HSEX_FROZEN_ETA_CHI_LIMIT q/error ',iq,error
              if(code/=0)status=34
              if(status==0)then
                if(window_init)then
                  window_cache%blocks(iq)%chi(:,:,1)=extrapolated
                  window_cache%blocks(iq)%error(1)=sum([1d0,6d0,8d0]* &
                    .5d0*(blocks(iq)%error(:3)+blocks(jq)%error(:3)))/3d0
                endif
                ! Extrapolate chi before inversion, retaining the existing
                ! positive dielectric guard on the eta=0 candidate as well.
                call screen_assemble(n,boxes,qcoords,rec,vol,nk,-.5d0*extrapolated,.true.,static_w, &
                  heads(iq),qheads(iq),code)
                if(code/=0)then
                  status=33
                else
                  call screen_channels(box_value(:,:,static_samples:static_samples),reshape(static_w,[na,na,1]),error,.false.)
                  box_value(:,:,neta)=static_w
                endif
              endif
              deallocate(extrapolated,static_w)
            else
              call screen_channels(box_value(:,:,1:1),box_value(:,:,2:2),error,.false.)
            endif
            if(.not.window_active.or.window_init)then
              write(*,'(a,i6,es18.8)')'HSEX_FROZEN_ETA_ERROR q/error ',iq,error
              if(error>closure_tol)status=34
            endif
          endif
          endif
        endif
        if(allocated(raw))deallocate(raw,qcoords,eta_chi)
      endif
      call comm_bcast(status,icomm,0)
      if(status/=0)return
      call comm_bcast(box_value,icomm,0)
      do j=1,nlocal
        ia=mod(first+j-2,na)+1;ib=(first+j-2)/na+1;coarse(j,iq)=box_value(ia,ib,neta)
      enddo
    enddo
    if(density_limit)then
      call density_gamma_replace(coarse,first,na,mesh,rec,vol,density_tolerance,icomm,rank,status)
      if(status/=0)then
        status=36;return
      endif
    else
    call comm_bcast(heads,icomm,0);call comm_bcast(qheads,icomm,0)
    if(.not.any(qheads>0d0))then
      status=35;return
    endif
    qmin=minval(qheads,mask=qheads>0d0);einv=0d0;j=0
    do iq=1,nq
      if(abs(qheads(iq)-qmin)>1d-8+1d-5*qmin)cycle
      einv=einv+heads(iq);j=j+1
    enddo
    einv=einv/j;qc=(6d0*acos(-1d0)**2/(vol*nk))**(1d0/3d0)
    if(present(velocity))then
      if(.not.comm_is_root(rank))allocate(origin(na,na))
      call comm_bcast(origin,icomm,0)
      do j=1,nlocal
        ia=mod(first+j-2,na)+1;ib=(first+j-2)/na+1;coarse(j,1)=origin(ia,ib)
      enddo
    endif
    endif
    call screen_transform(coarse,mesh)
    if(.not.present(velocity).and..not.density_limit)coarse=coarse+12d0*acos(-1d0)*einv/(qc*qc*vol*nk)
    if(window_active.and..not.window_init)then
      coarse=coarse+cmplx(window_cache%anchor_w-window_cache%reference_w,0d0,8)
    endif
    code=0
    if(.not.screen_finite(coarse))code=1
    if(size(coarse)>0)then
      if(maxval(abs(aimag(coarse)))>1d-7)code=1
    endif
    allocate(local_matrix(na,na),current(na,na),partner(na,na),stat=ios)
    code=max(code,merge(0,1,ios==0));call comm_get_max(code,icomm)
    if(code/=0)then
      status=40;return
    endif
    do iq=1,nq
      local_matrix=0d0
      do j=1,nlocal
        ia=mod(first+j-2,na)+1;ib=(first+j-2)/na+1;local_matrix(ia,ib)=coarse(j,iq)
      enddo
      call comm_summation(local_matrix,current,na*na,icomm)
      coord=modulo(-sites(:,iq),mesh);jq=1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3))
      local_matrix=0d0
      do j=1,nlocal
        ia=mod(first+j-2,na)+1;ib=(first+j-2)/na+1;local_matrix(ia,ib)=coarse(j,jq)
      enddo
      call comm_summation(local_matrix,partner,na*na,icomm)
      if(maxval(abs(current-transpose(partner)))>1d-8)code=1
    enddo
    call comm_get_max(code,icomm)
    if(code/=0)then
      status=41;return
    endif
    w=real(coarse,8);status=0
    if(window_init)then
      window_cache%n=n;window_cache%mesh=mesh;window_cache%boxes=boxes;window_cache%nv=nv
      window_cache%first=first;window_cache%nlocal=nlocal;window_cache%kind=filter;window_cache%time=area
      window_cache%cutoff=cutoff;window_cache%eta=eta;window_cache%cell=cell;window_cache%k=k
      window_cache%density_limit=density_limit;window_cache%reference_w=w
      window_cache%closure_tol=closure_tol;window_cache%density_tol=density_tolerance
      window_cache%ready=.true.;window_cache%anchored=.false.
    endif
    call system_clock(finished);seconds=dble(finished-started)/rate
    peak_rss(1)=hsex_peak_rss_mib()
    call comm_get_max(peak_rss,max_rss,1,icomm)
    if(comm_is_root(rank))then
      write(*,'(a,es18.8)')'HSEX_FROZEN_PEAK_RSS_MIB rankmax ',max_rss(1)
      write(*,'(a,i12)')'HSEX_FROZEN_PRECONDITION_ACTIONS ',total_precondition
      if(density_limit)then
        write(*,'(a)')'HSEX_FROZEN_DENSITY_LIMIT_PASS finite-q density eta limit and q2 Gamma cell'
      else if(present(velocity))then
        write(*,'(a)')'HSEX_FROZEN_LIMIT_PASS stationary optical head/wings; finite-q eta extrapolated'
      else
        write(*,'(a,es18.8)')'HSEX_FROZEN_HEAD_MIN_SHELL eps_M ',1d0/einv
      endif
      if(window_active)then
        write(*,'(a,4i8,2es18.8)')'HSEX_FROZEN_WINDOW_DONE nq/bases/actions/error_bound/seconds ', &
          nq,basis,total_actions,worst,seconds
      else
        write(*,'(a,4i8,2es18.8)')'HSEX_FROZEN_FULL_DONE nq/bases/actions/residual/seconds ', &
          nq,basis,total_actions,worst,seconds
      endif
    endif
  contains
    subroutine validate_window(code)
      implicit none
      integer,intent(out) :: code
      code=1
      if(filter<1.or.filter>2.or..not.screen_finite([area]))return
      if(area<=0d0.or..not.project_static.or..not.present(prepare).or.present(velocity))return
      if(.not.present(window_cache))return
      if(window_init)then
        if(window_cache%ready.or.allocated(window_cache%blocks))return
      else
        if(.not.window_cache%ready.or..not.window_cache%anchored)return
        if(any(window_cache%n/=n).or.any(window_cache%mesh/=mesh).or.window_cache%boxes/=boxes)return
        if(window_cache%nv/=nv.or.window_cache%first/=first.or.window_cache%nlocal/=nlocal)return
        if(window_cache%kind/=filter.or.window_cache%time/=area.or.window_cache%eta/=eta)return
        if(window_cache%cutoff/=cutoff.or.any(window_cache%cell/=cell))return
        if(window_cache%closure_tol/=closure_tol.or.window_cache%density_tol/=density_tolerance)return
        if(window_cache%density_limit.neqv.density_limit)return
        if(.not.allocated(window_cache%k))return
        if(any(shape(window_cache%k)/=shape(k)))return
        if(any(window_cache%k/=k))return
      endif
      code=0
    end subroutine
    subroutine left_precondition(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      call precondition(target,x,y,ierr)
    end subroutine
    subroutine right_precondition(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      call precondition(source,x,y,ierr)
    end subroutine
    subroutine parallel_response()
      implicit none
      integer :: qb,kb,worker,kworkers,qi,sx,tx,myq,mysource,mytarget,mm,co(3),uu(3),ii,jj,kk,gg,ss
      integer :: local_stats(4),sent_stats(4)
      real(8) :: local_worst,sent_worst,reference_rate
      real(8),allocatable :: pair_error(:),sent_error(:),response_error(:)
      complex(8),allocatable :: payload(:),pair_chi(:,:,:),href(:,:)
      kworkers=np/qworkers
      if(window_active)then
        allocate(pair_error(neta),sent_error(neta),stat=code)
        call comm_get_max(code,icomm)
        if(code/=0)then
          status=20;return
        endif
      endif
      allocate(blocks(nq),ps(ng,nv),pt(ng,nv),stat=code)
      call comm_get_max(code,icomm)
      if(code/=0)then
        status=20;return
      endif
      limits=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell**2,dim=1))/(2d0*acos(-1d0))+1d0)
      do qi=1,nq
        if(density_limit.and.qi==1)cycle
        q=dble(sites(:,qi))/qmesh;mm=0
        do ii=-limits(1),limits(1)
          do jj=-limits(2),limits(2)
            do kk=-limits(3),limits(3)
              if(sum(matmul(rec,[ii,jj,kk]+q)**2)/2d0<=cutoff)mm=mm+1
            enddo
          enddo
        enddo
        if(mm<1)then
          status=20;return
        endif
        allocate(blocks(qi)%g(3,mm),blocks(qi)%reverse(mm),stat=code)
        call comm_get_max(code,icomm)
        if(code/=0)then
          status=20;return
        endif
        co=modulo(-sites(:,qi),qmesh)
        blocks(qi)%minus=1+co(1)+qmesh(1)*(co(2)+qmesh(2)*co(3))
        gg=0
        do ii=-limits(1),limits(1)
          do jj=-limits(2),limits(2)
            do kk=-limits(3),limits(3)
              if(sum(matmul(rec,[ii,jj,kk]+q)**2)/2d0>cutoff)cycle
              gg=gg+1;blocks(qi)%g(:,gg)=[ii,jj,kk]
            enddo
          enddo
        enddo
        if(comm_is_root(rank))then
          allocate(blocks(qi)%chi(mm,mm,neta),stat=code)
          if(code==0)then
            blocks(qi)%chi=0d0
            if(window_active)then
              allocate(blocks(qi)%error(neta),stat=code)
              if(code==0)blocks(qi)%error=0d0
            endif
          endif
        endif
        call comm_get_max(code,icomm)
        if(code/=0)then
          status=20;return
        endif
      enddo
      do qb=1,nq,qworkers
        do kb=1,nk,kworkers
          myq=0;mysource=0;mytarget=0
          do worker=0,np-1
            qi=qb+worker/kworkers;sx=kb+mod(worker,kworkers)
            if(qi>nq.or.sx>nk)cycle
            if(density_limit.and.qi==1)cycle
            co=modulo(order(:,sx)+sites(:,qi),mesh)
            tx=lookup(1+co(1)+mesh(1)*(co(2)+mesh(2)*co(3)))
            call prepare(sx,tx,worker,owners(sx),owners(tx),code)
            call comm_get_max(code,icomm)
            if(code/=0)then
              status=25;return
            endif
            if(rank==owners(sx))orbitals=states(:,:,sx-ks+1)
            call comm_bcast(orbitals,icomm,owners(sx))
            if(rank==worker)ps=orbitals
            if(rank==owners(tx))orbitals=states(:,:,tx-ks+1)
            call comm_bcast(orbitals,icomm,owners(tx))
            if(rank==worker)then
              pt=orbitals;myq=qi;mysource=sx;mytarget=tx
            endif
          enddo
          local_stats=0;local_worst=0d0;code=0
          if(myq>0)then
            source=mysource;target=mytarget;mm=size(blocks(myq)%g,2)
            if(window_init)then
              allocate(href(ng,nv),stat=code)
              if(code==0)then
                do ii=1,nv
                  call right_action(ps(:,ii),href(:,ii),ss)
                  if(ss/=0)then
                    code=1;exit
                  endif
                enddo
                if(code==0)call frozen_background_rate(ps,href,reference_rate,code)
                if(code==0)then
                  if(reference_rate>max(1d-8,min(1d-6,tol)))code=1
                endif
                deallocate(href)
              endif
            endif
            if(code==0)allocate(vertices(ng,mm),flat(neta*mm),pair_chi(mm,mm,neta),stat=code)
            if(code==0)then
              q=dble(sites(:,myq))/qmesh
              uu=nint(frac(:,source)+q-frac(:,target))
              do gg=1,mm
                co=blocks(myq)%g(:,gg)+uu
                if(any(2*abs(co)>=n))code=21
                do ii=0,ng-1
                  position=dble([mod(ii,n(1)),mod(ii/n(1),n(2)),ii/(n(1)*n(2))])/n
                  vertices(ii+1,gg)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*co),8))
                enddo
              enddo
              if(code==0)call frozen_sources(ps,pt,2d0*[(1d0,ii=1,nv)],2d0*[(1d0,ii=1,nv)],vertices,u,v,coeff,code)
              if(code==0)then
                if(present(precondition))then
                  call frozen_solve(left_action,right_action,u,coeff(:,:,1),v,shifts,capacity,tol,solution,stats,code, &
                    extra_sources=coeff(:,:,2:),residual_bound=.true., &
                    precondition_left=left_precondition,precondition_right=right_precondition,new_vectors_only=.true., &
                    window_kind=filter_kind,window_time=filter_time)
                else
                  call frozen_solve(left_action,right_action,u,coeff(:,:,1),v,shifts,capacity,tol,solution,stats,code, &
                    extra_sources=coeff(:,:,2:),residual_bound=.true.,new_vectors_only=.true., &
                    window_kind=filter_kind,window_time=filter_time)
                endif
              endif
              if(code==0)then
                local_stats=[stats%left_size,stats%right_size,stats%h_actions,stats%precondition_actions]
                if(window_init)local_stats(3)=local_stats(3)+nv
                local_worst=maxval(solution%residual)
                if(window_active)then
                  allocate(response_error(mm),stat=code)
                  if(code==0)then
                    do c=1,neta
                      do jj=1,mm
                        ss=(jj-1)*neta+c
                        response_error(jj)=solution%residual(ss)*solution%error_scale(ss)
                      enddo
                      if(filter_kind(c)==0)response_error=response_error/aimag(shifts(c))
                      pair_error(c)=sqrt(dble(mm)*ng)*sqrt(sum(response_error**2))
                    enddo
                    if(.not.screen_finite(pair_error))code=1
                    deallocate(response_error)
                  endif
                endif
                do gg=1,mm
                  call frozen_trace(solution,vertices(:,gg),flat,ss)
                  code=max(code,ss)
                  do jj=1,mm
                    pair_chi(gg,jj,:)=flat(neta*(jj-1)+1:neta*jj)
                  enddo
                enddo
              endif
            endif
          endif
          call comm_get_max(code,icomm)
          if(code/=0)then
            status=23;return
          endif
          ! Only converged response blocks are transferred; no MPI in Krylov actions.
          do worker=0,np-1
            qi=qb+worker/kworkers;sx=kb+mod(worker,kworkers)
            if(qi>nq.or.sx>nk)cycle
            if(density_limit.and.qi==1)cycle
            mm=size(blocks(qi)%g,2);allocate(payload(mm*mm*neta),stat=code)
            call comm_get_max(code,icomm)
            if(code/=0)then
              status=20;return
            endif
            if(rank==worker)payload=reshape(pair_chi,[mm*mm*neta])
            call comm_bcast(payload,icomm,worker)
            if(window_active)then
              if(rank==worker)sent_error=pair_error
              call comm_bcast(sent_error,icomm,worker)
            endif
            if(rank==worker)then
              sent_stats=local_stats;sent_worst=local_worst
            endif
            call comm_bcast(sent_stats,icomm,worker)
            call comm_bcast(sent_worst,icomm,worker)
            if(comm_is_root(rank))then
              blocks(qi)%chi=blocks(qi)%chi+reshape(payload,[mm,mm,neta])
              if(window_active)blocks(qi)%error=blocks(qi)%error+sent_error
              basis=max(basis,sent_stats(:2));total_actions=total_actions+sent_stats(3)
              total_precondition=total_precondition+sent_stats(4)
              worst=max(worst,sent_worst)
              if(window_active)then
                write(*,'(a,5i7,es16.7)')'HSEX_FROZEN_WINDOW_PAIR q/k/rank/G/actions/error_bound ', &
                  qi,sx,worker,mm,sent_stats(3),sent_worst
              else
                write(*,'(a,5i7,es16.7)')'HSEX_FROZEN_MPI_PAIR q/k/rank/G/actions/residual ', &
                  qi,sx,worker,mm,sent_stats(3),sent_worst
              endif
              flush(output_unit)
            endif
            deallocate(payload)
          enddo
          if(allocated(vertices))deallocate(vertices,flat,pair_chi)
        enddo
      enddo
    end subroutine
    subroutine request()
      implicit none
      call comm_bcast(command,icomm,0)
      call service()
    end subroutine
    subroutine service()
      implicit none
      owner=owners(command(2));request_code=0
      if(command(1)==1)then
        call comm_bcast(buffer,icomm,0)
        if(rank==owner)then
          call apply(command(2),buffer,answer,request_code)
          if(request_code==0)then
            if(.not.screen_finite(answer))request_code=1
          endif
          if(request_code/=0)answer=0d0
        endif
        call comm_bcast(request_code,icomm,owner);call comm_bcast(answer,icomm,owner)
      else if(command(1)==3)then
        call comm_bcast(buffer,icomm,0)
        if(rank==owner)then
          request_code=1
          if(present(velocity))call velocity(command(2),buffer,velocity_answer,request_code)
          if(request_code==0)then
            if(.not.screen_finite(velocity_answer))request_code=1
          endif
          if(request_code/=0)velocity_answer=0d0
        endif
        call comm_bcast(request_code,icomm,owner);call comm_bcast(velocity_answer,icomm,owner)
      else
        if(rank==owner)orbitals=states(:,:,command(2)-ks+1)
        call comm_bcast(orbitals,icomm,owner)
      endif
    end subroutine
    subroutine left_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      if(present(prepare))then
        call apply(target,x,y,ierr)
      else
        command=[1,target];buffer=x;call request();y=answer;ierr=request_code
      endif
    end subroutine
    subroutine right_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      if(present(prepare))then
        call apply(source,x,y,ierr)
      else
        command=[1,source];buffer=x;call request();y=answer;ierr=request_code
      endif
    end subroutine
    subroutine velocity_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: ierr
      command=[3,source];buffer=x;call request();y=velocity_answer;ierr=request_code
    end subroutine
  end subroutine
end module
