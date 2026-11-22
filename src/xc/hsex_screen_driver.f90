! Screening orchestration. Input orbitals have Euclidean norm one (sqrt(dv)*psi).
! Dense operators are loaded only for destination k points owned by this rank.
module hsex_screen_driver
  implicit none
  private
  public :: screen_compute
  type :: q_block
    integer,allocatable :: g(:,:),reverse(:)
    complex(8),allocatable :: chi(:,:,:),box(:,:,:)
    integer :: minus,owner
  end type
  abstract interface
    subroutine h_loader(j,h,status)
      implicit none
      integer,intent(in)::j
      complex(8),intent(out)::h(:,:)
      integer,intent(out)::status
    end subroutine
    subroutine h_action(j,x,y,status)
      implicit none
      integer,intent(in)::j
      complex(8),intent(in)::x(:)
      complex(8),intent(out)::y(:)
      integer,intent(out)::status
    end subroutine
    subroutine velocity_action(j,x,y,status)
      implicit none
      integer,intent(in) :: j
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:,:)
      integer,intent(out) :: status
    end subroutine
  end interface
contains
  subroutine screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states,nv,eta,head_tol,eta_tol, &
      nthreads,ik_s,ik_e,icomm,rank,nproc,load_h,first_pair,nlocal,w,status,apply_h,temperature,head_only, &
      distributed_states,apply_velocity,apply_optical_h,reduced_rank,reduced_batch,reduced_tolerance, &
      reduced_symmetry,atoms,species)
    use, intrinsic :: iso_fortran_env, only: error_unit,output_unit
    use communication, only: comm_summation,comm_bcast,comm_get_max,comm_is_root,comm_alltoallv_complex, &
      comm_thread_funneled
    use hsex_screen_math, only: screen_sites,screen_geometry,screen_assemble,screen_population,screen_channels, &
      screen_metal_head,screen_forms,screen_transform,screen_finite
    use hsex_q_response, only: hsex_q_shared,hsex_q_workspace,hsex_q_prepare,hsex_q_release, &
      hsex_q_worker_prepare,hsex_q_worker_release,hsex_q_one
    use hsex_sternheimer, only: sternheimer_band
    use hsex_longwave, only: longwave_band,longwave_cell
    use hsex_response_space, only: response_space_solve,s_response_space_stats
    use hsex_response_symmetry, only: response_q_stars,s_response_stars
    implicit none
    integer,intent(in)::n(3),mesh(3),qmesh(3),boxes,nv,nthreads,ik_s,ik_e,icomm,rank,nproc,first_pair,nlocal
    real(8),intent(in)::cutoff,cell(3,3),k(:,:),levels(:,:),pop(:,:),eta,head_tol,eta_tol
    complex(8),intent(in)::states(:,:,:)
    procedure(h_loader)::load_h
    procedure(h_action),optional::apply_h
    procedure(velocity_action),optional :: apply_velocity
    procedure(h_action),optional :: apply_optical_h
    real(8),intent(in),optional::temperature
    integer,optional,intent(in) :: reduced_rank,reduced_batch,species(:)
    real(8),optional,intent(in) :: reduced_tolerance,atoms(:,:)
    logical,optional,intent(in) :: reduced_symmetry
    logical :: compressed,compressed_symmetry
    complex(8),allocatable :: occupied_all(:,:,:),occupied_part(:,:,:)
    integer :: compressed_batch,compressed_capacity
    real(8) :: compressed_tolerance
    complex(8),allocatable::vertex(:,:),band_response(:,:)
    logical,intent(in),optional::head_only,distributed_states
    logical,allocatable :: active_q(:)
    logical::matrix_free,only_head,distributed
    integer :: orbital_offset,local_k,k_slots,k_slot,peer,batch_threads,q_regions,q_batches
    integer,allocatable :: ranges(:,:),local_range(:,:),k_owner(:)
    integer::band,project_bands
    real(8),intent(out)::w(nlocal,product(mesh))
    integer,intent(out)::status
    type(q_block),allocatable::blocks(:)
    integer,allocatable::sites(:,:),qsites(:,:),order(:,:),lookup(:),gall(:,:),counts(:),istat(:),samples(:)
    integer::stage
    integer::nk,nq,ng,nreal,ncomponent,component,nb,na,nc,igmax(3),allg,i,j,l,a,b,c,iq,jq,ig,jg,owner,local_error,err,team
    integer::idx(3),shift(3),source,beg,last,batch,mq,mg,nu,slot,nsel,ia,ib,dx,dy,dz,lo(3),vtx(3)
    real(8)::rec(3,3),vol,dv,q(3),q2,head,qhead,einv,qmin,error,rmax,radius,qc,weight,fraction(3),point3(3)
    real(8) :: spin_weight,solver_temperature
    real(8),allocatable::frac(:,:),d(:),e(:),lev(:,:),occ(:,:),etas(:),residual(:),heads(:),qheads(:),radii(:)
    complex(8),allocatable::h(:,:),original_h(:,:),adj(:,:),u(:,:),tau(:),work(:),hu(:,:),huq(:,:),gram(:,:)
    complex(8),allocatable::psi(:,:,:),phase(:,:,:),response(:,:,:,:),reduced(:,:,:),qvec(:,:),matrix(:,:)
    complex(8),allocatable::coarse(:,:),values(:,:),tail(:,:),origin(:,:),point(:,:),head_samples(:,:,:,:)
    complex(8),allocatable::buf(:,:,:),coarse_origin(:,:),coarse_point(:,:),current(:,:),partner(:,:)
    real(8),allocatable::Qcoords(:,:)
    logical::pump,insulating
    logical :: optical
    integer,allocatable :: optical_g(:,:),optical_reverse(:)
    complex(8),allocatable :: optical_chi(:,:),optical_total(:,:),optical_h(:,:)
    real(8) :: macro(3,3)
    integer(8)::clock_start,clock_response,clock_assembly,clock_end,clock_rate,clock_phase
    real(8)::timing(3),timing_max(3),max_residual
    real(8) :: phase_seconds(11),phase_max(11)
    call system_clock(clock_start,clock_rate)
    phase_seconds=0d0;q_regions=0;q_batches=0
    stage=10
    status=1;w=0d0;nk=product(mesh);nq=product(qmesh);nreal=product(n);nb=size(levels,1);na=boxes**3
    ng=size(states,1)
    local_error=0
    compressed=present(reduced_rank);compressed_symmetry=.false.
    compressed_batch=4;compressed_capacity=384;compressed_tolerance=1d-6
    if(present(reduced_rank))compressed_capacity=reduced_rank
    if(present(reduced_batch))compressed_batch=reduced_batch
    if(present(reduced_tolerance))compressed_tolerance=reduced_tolerance
    if(present(reduced_symmetry))compressed_symmetry=reduced_symmetry
    if(compressed)then
      if(compressed_batch<1.or.compressed_capacity<1)local_error=1
      if(.not.screen_finite([compressed_tolerance]))local_error=1
      if(compressed_tolerance<=0d0.or.compressed_tolerance>=1d0)local_error=1
      if(.not.present(apply_h))local_error=1
      if(compressed_symmetry.and.(.not.present(atoms).or..not.present(species)))local_error=1
    endif
    distributed=.false.
    if(present(distributed_states))distributed=distributed_states
    local_k=max(0,ik_e-ik_s+1)
    orbital_offset=0
    if(distributed)orbital_offset=ik_s-1
    if(any(n<1).or.any(mesh<1).or.any(qmesh<1).or.boxes<1.or.nthreads<1) local_error=1
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    if(any(mod(mesh,qmesh)/=0).or.nq<2.or.any(mod(n,boxes)/=0).or.cutoff<=0d0) local_error=1
    ncomponent=ng/nreal
    if(ng/=ncomponent*nreal.or.ncomponent<1.or.ncomponent>2) local_error=1
    if(size(states,2)/=nb.or.size(states,3)/=merge(local_k,nk,distributed).or. &
       any(shape(pop)/=[nb,nk]).or.any(shape(k)/=[3,nk])) local_error=1
    if(size(levels,2)/=nk.or.nb<1.or.nv<1.or.nv>nb)local_error=1
    if(rank<0.or.rank>=nproc.or.nproc<1.or.ik_s<1.or.ik_e>nk)local_error=1
    if(.not.screen_finite(k))local_error=1
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    if(first_pair<1.or.nlocal<0.or.(nlocal>0.and.first_pair+nlocal-1>na*na)) local_error=1
    if(eta<0d0.or.head_tol<=0d0.or.head_tol>=1d0.or.eta_tol<=0d0.or.eta_tol>=1d0) local_error=1
    if(.not.screen_finite([cutoff,eta,head_tol,eta_tol]))local_error=1
    stage=20
    call screen_geometry(cell,rec,vol,err);local_error=max(local_error,err)
    matrix_free=present(apply_h)
    if(.not.matrix_free.and.nthreads>1)then
      if(.not.comm_thread_funneled())then
        local_error=1
        if(comm_is_root(rank))write(error_unit,'(a)')'HSEX_W_REJECT: outer q team requires primary MPI_THREAD_FUNNELED'
      endif
    endif
    only_head=.false.
    if(present(head_only))only_head=head_only
    if(only_head.and..not.matrix_free)local_error=1
    spin_weight=dble(3-ncomponent)
    solver_temperature=1d0
    if(matrix_free.and.ncomponent==2)then
      if(.not.present(temperature)) local_error=1
      if(present(temperature))then
        if(.not.screen_finite([temperature]).or.temperature<=0d0) local_error=1
        solver_temperature=temperature
      endif
      if(any(pop(nb,:)>1d-12)) local_error=1
    endif
    pump=eta>0d0
    if(ncomponent==1)then
      call screen_population(levels,pop,nv,insulating,err)
      local_error=max(local_error,err)
    else
      ! A spinor band holds one electron. Equilibrium metallic occupations
      ! conserve charge over the full uniform mesh, rather than at each k.
      insulating=.false.
      if(.not.pump) local_error=1
      if(.not.screen_finite(levels).or..not.screen_finite(pop))then
        local_error=1
      else
        if(any(levels(2:,:)-levels(:nb-1,:)<-1d-10)) local_error=1
        if(minval(pop)<-1d-10.or.maxval(pop)>1d0+1d-7) local_error=1
        if(abs(sum(pop)/nk-nv)>1d-4) local_error=1
        if(nb>nv.and.present(apply_velocity))then
          insulating=maxval(abs(pop(:nv,:)-1d0))<1d-12.and.maxval(abs(pop(nv+1:,:)))<1d-12
          insulating=insulating.and.minval(levels(nv+1,:))-maxval(levels(nv,:))>1d-8
        endif
      endif
    endif
    if(.not.screen_finite(states))local_error=1
    if(.not.pump)then
      if(any(pop(:nv,:)/=2d0).or.any(pop(nv+1:,:)/=0d0))local_error=1
      insulating=.true.
    endif
    if(matrix_free.and.ncomponent==1)then
      ! Scalar CG currently represents the static, integer-filled insulator.
      ! An empty target band demonstrates a common positive gap; fractional
      ! scalar pump screening retains its existing dense/eta treatment.
      if(pump.or.nb<=nv)then
        local_error=1
      else
        if(any(pop(:nv,:)/=2d0).or.any(pop(nv+1:,:)/=0d0))local_error=1
        if(minval(levels(nv+1,:))-maxval(levels(nv,:))<=1d-8)local_error=1
      endif
    endif
    project_bands=nb
    if(matrix_free.and.ncomponent==1)project_bands=nv
    if(compressed.and.(ncomponent/=1.or.pump.or..not.insulating.or.only_head))local_error=1
    optical=insulating.and.present(apply_velocity)
    if(only_head.and.optical)local_error=1
    if(comm_is_root(rank).and.pump.and.ncomponent==1)then
      write(*,'(a,l1,a,es12.4,a,es12.4)')'HSEX_W_POP insulating=',insulating, &
        ' max per-k missing=',maxval(abs(sum(pop,dim=1)-2*nv)), &
        ' max valence departure=',maxval(abs(pop(:nv,:)-2d0))
      if(nb>nv)write(*,'(a,es12.4)')'HSEX_W_POP indirect gap Ha=',minval(levels(nv+1:,:))-maxval(levels(:nv,:))
    endif
    if(.not.insulating.and.any(qmesh<6))then
      local_error=1
      if(comm_is_root(rank))write(*,'(a)')'HSEX_W_REJECT: fractional/non-common manifold needs q mesh >=6 in every direction'
    endif
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    nc=merge(2,1,pump);nu=merge(0,nv,pump);dv=vol/nreal
    allocate(sites(3,nk),qsites(3,nq),order(3,nk),frac(3,nk),lookup(nk),blocks(nq),heads(nq),qheads(nq))
    call screen_sites(mesh,sites);call screen_sites(qmesh,qsites)
    stage=30
    frac=matmul(transpose(cell),k)/(2d0*acos(-1d0));lookup=0
    do i=1,nk
      point3=(frac(:,i)-frac(:,1))*mesh;order(:,i)=nint(point3)
      if(maxval(abs(point3-order(:,i)))>1d-8)local_error=1
      idx=modulo(order(:,i),mesh);slot=1+idx(1)+mesh(1)*(idx(2)+mesh(2)*idx(3))
      if(lookup(slot)/=0)local_error=1
      lookup(slot)=i
    enddo
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    allocate(active_q(nq))
    active_q=.true.
    if(.not.insulating)then
      allocate(radii(nq),samples(nq));nsel=0
      do iq=1,nq
        radii(iq)=huge(1d0)
        if(any(qsites(:,iq)/=0).and.all(2*abs(qsites(:,iq))<qmesh))then
          radii(iq)=sum(matmul(rec,dble(qsites(:,iq))/qmesh)**2);nsel=nsel+1
        endif
      enddo
      if(nsel<52)then
        status=71
        if(comm_is_root(rank))write(*,'(a)')'HSEX_W_REJECT: fewer than 52 metallic head samples'
        return
      endif
      block
        real(8)::sorted(nq),tmp
        sorted=radii
        do a=2,nq
          tmp=sorted(a);b=a-1
          do while(b>=1)
            if(sorted(b)<=tmp)exit
            sorted(b+1)=sorted(b);b=b-1
          enddo
          sorted(b+1)=tmp
        enddo
        radius=sorted(52)*(1d0+1d-10)
      end block
      nsel=0
      do iq=1,nq
        if(radii(iq)>radius)cycle
        nsel=nsel+1;samples(nsel)=iq
      enddo
      if(only_head)then
        active_q=.false.
        active_q(samples(:nsel))=.true.
        if(comm_is_root(rank))write(output_unit,'(a,i0,a,i0)')'HSEX_HEAD_ONLY q samples=',nsel,' of ',nq
      endif
    endif
    stage=40
    do a=1,3
      igmax(a)=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell(:,a)**2))/(2d0*acos(-1d0))+1d0)
    enddo
    allg=product(2*igmax+1);allocate(gall(3,allg));i=0
    ! Match Python product ordering: Gz varies fastest.
    do a=-igmax(1),igmax(1)
      do b=-igmax(2),igmax(2)
        do c=-igmax(3),igmax(3)
          i=i+1;gall(:,i)=[a,b,c]
        enddo
      enddo
    enddo
    do iq=1,nq
      if(.not.active_q(iq))cycle
      q=dble(qsites(:,iq))/qmesh;mg=0
      do i=1,allg
        if(sum(matmul(rec,gall(:,i)+q)**2)/2d0<=cutoff)mg=mg+1
      enddo
      if(mg<1)local_error=1
      allocate(blocks(iq)%g(3,mg),blocks(iq)%chi(mg,mg,nc),blocks(iq)%reverse(mg))
      blocks(iq)%chi=0d0;j=0
      do i=1,allg
        if(sum(matmul(rec,gall(:,i)+q)**2)/2d0>cutoff)cycle
        j=j+1;blocks(iq)%g(:,j)=gall(:,i)
      enddo
      idx=modulo(-qsites(:,iq),qmesh);jq=1+idx(1)+qmesh(1)*(idx(2)+qmesh(2)*idx(3))
      blocks(iq)%minus=jq;blocks(iq)%owner=mod(min(iq,jq)-1,nproc)
    enddo
    do iq=1,nq
      if(.not.active_q(iq))cycle
      jq=blocks(iq)%minus
      do ig=1,size(blocks(iq)%g,2)
        idx=-blocks(iq)%g(:,ig)-(qsites(:,iq)+qsites(:,jq))/qmesh
        blocks(iq)%reverse(ig)=0
        do jg=1,size(blocks(jq)%g,2)
          if(all(idx==blocks(jq)%g(:,jg)))blocks(iq)%reverse(ig)=jg
        enddo
        if(blocks(iq)%reverse(ig)==0)local_error=1
      enddo
    enddo
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    if(.not.matrix_free)then
      allocate(h(ng,ng),adj(ng,ng),d(ng),e(ng-1),u(ng,nu),tau(ng-1),work(64*ng),etas(nc))
      allocate(original_h(merge(ng,1,pump),merge(ng,1,pump)))
      etas(1)=eta
      if(pump)etas(2)=eta/2d0
    endif
    if(optical)then
      mg=size(blocks(1)%g,2)-1
      allocate(optical_g(3,mg),optical_reverse(mg+3),optical_chi(mg+3,mg+3),optical_total(mg+3,mg+3))
      mg=0
      do i=1,size(blocks(1)%g,2)
        if(all(blocks(1)%g(:,i)==0))cycle
        mg=mg+1;optical_g(:,mg)=blocks(1)%g(:,i)
      enddo
      optical_reverse(:3)=[1,2,3]
      do i=1,mg
        do a=1,mg
          if(all(optical_g(:,a)==-optical_g(:,i)))optical_reverse(i+3)=a+3
        enddo
      enddo
      optical_chi=0d0
      if(.not.matrix_free.and..not.present(apply_optical_h))allocate(optical_h(ng,ng))
    endif
    stage=50
    rmax=0d0
    k_slots=local_k
    batch_threads=nthreads
    if(distributed)then
      ! Every rank must split the q loop at the same batch boundaries.
      batch_threads=-nthreads
      call comm_get_max(batch_threads,icomm)
      batch_threads=-batch_threads
      allocate(ranges(2,nproc),local_range(2,nproc),k_owner(nk))
      local_range=0
      local_range(:,rank+1)=[ik_s,ik_e]
      call comm_summation(local_range,ranges,size(ranges),icomm)
      k_owner=0
      do peer=1,nproc
        if(ranges(1,peer)<1.or.ranges(2,peer)>nk)then
          status=stage
          return
        endif
        if(any(k_owner(ranges(1,peer):ranges(2,peer))/=0))then
          status=stage
          return
        endif
        k_owner(ranges(1,peer):ranges(2,peer))=peer
      enddo
      if(any(k_owner==0))then
        status=stage
        return
      endif
      k_slots=maxval(max(0,ranges(2,:)-ranges(1,:)+1))
    endif
    if(compressed)then
      allocate(occupied_all(ng,nv,nk),stat=err)
      local_error=merge(1,0,err/=0);call comm_get_max(local_error,icomm)
      if(local_error/=0)then
        status=50
        return
      endif
      if(distributed)then
        allocate(occupied_part(ng,nv,nk),stat=err)
        local_error=merge(1,0,err/=0);call comm_get_max(local_error,icomm)
        if(local_error/=0)then
          status=50
          return
        endif
        occupied_part=0d0
        if(ik_e>=ik_s)occupied_part(:,:,ik_s:ik_e)=states(:,:nv,:)
        call comm_summation(occupied_part,occupied_all,size(occupied_all),icomm)
        deallocate(occupied_part)
      else
        occupied_all=states(:,:nv,:)
      endif
    endif
    if(matrix_free)then
      do j=ik_s,ik_e
        gram=matmul(conjg(transpose(states(:,:,j-orbital_offset))),states(:,:,j-orbital_offset))
        do a=1,nb
          gram(a,a)=gram(a,a)-1d0
        enddo
        if(.not.screen_finite(gram).or.maxval(abs(gram))>1d-8)then
          write(error_unit,*)'HSEX_STERNHEIMER_ORTHONORMALITY rank/k/error',rank,j,maxval(abs(gram))
          local_error=1
          exit
        endif
        deallocate(gram)
        allocate(vertex(ng,1))
        do a=1,project_bands
          call apply_h(j,states(:,a,j-orbital_offset),vertex(:,1),err)
          if(err/=0.or..not.screen_finite(vertex))then
            write(error_unit,*)'HSEX_STERNHEIMER_ACTION rank/k/band/status',rank,j,a,err
            local_error=1
            exit
          endif
          error=norm2(abs(vertex(:,1)-states(:,a,j-orbital_offset)*levels(a,j)))
          if(error>1d-7)then
            write(error_unit,*)'HSEX_STERNHEIMER_EIGENPAIR rank/k/band/residual',rank,j,a,error
            local_error=1
            exit
          endif
        enddo
        deallocate(vertex)
        if(local_error/=0)exit
      enddo
      err=local_error
      call comm_get_max(err,icomm)
      if(err/=0)then
        status=stage
        return
      endif
      if(comm_is_root(rank))then
        write(output_unit,'(a)')'HSEX_STERNHEIMER_EIGENPAIRS_PASS'
        flush(output_unit)
      endif
    endif
    do k_slot=1,k_slots
      j=ik_s+k_slot-1
      if(.not.matrix_free)then
        if(j<=ik_e)call prepare_dense_target()
        if(distributed)call comm_get_max(local_error,icomm)
        if(local_error/=0)exit
      endif
      if(optical.and.matrix_free.and.j<=ik_e)call accumulate_optical()
      if(distributed)call comm_get_max(local_error,icomm)
      if(local_error/=0)exit
      if(compressed)then
        if(j<=ik_e)call compressed_target()
        if(distributed)call comm_get_max(local_error,icomm)
        if(local_error/=0)exit
        cycle
      endif
      if(.not.matrix_free)then
        call dense_q_batches()
        if(local_error/=0)exit
        cycle
      endif
      beg=1
      do while(beg<=nq)
        if(.not.active_q(beg))then
          beg=beg+1
          cycle
        endif
        ! Bound padded source and phase arrays, allowing one oversized q.
        mg=0;last=beg-1
        do iq=beg,min(nq,beg+2*batch_threads-1)
          if(.not.active_q(iq))exit
          mq=max(mg,size(blocks(iq)%g,2))
          if(iq>beg.and.16d0*ng*(nb+mq)*(iq-beg+1)>64d0*1024**2)exit
          mg=mq;last=iq
        enddo
        batch=last-beg+1
        allocate(psi(ng,nb,batch),phase(ng,mg,batch),lev(nb,batch),occ(nb,batch),counts(batch), &
          response(mg,mg,nc,batch),residual(batch),istat(batch))
        psi=0d0;phase=0d0;occ=0d0
        if(distributed)then
          call system_clock(clock_phase)
          call exchange_source_batch(err)
          call record_phase(4,clock_phase)
          if(err/=0)then
            local_error=err
            exit
          endif
        endif
        if(j<=ik_e)then
          call system_clock(clock_phase)
          do iq=beg,last
            l=iq-beg+1
            source=source_k(j,iq)
            if(.not.distributed)psi(:,:,l)=states(:,:,source)
            lev(:,l)=levels(:,source);occ(:,l)=pop(:,source)
            if(.not.pump)occ(nv+1:,l)=0d0
            shift=nint(frac(:,source)+dble(qsites(:,iq))/qmesh-frac(:,j));counts(l)=size(blocks(iq)%g,2)
            do ig=1,counts(l)
              idx=blocks(iq)%g(:,ig)+shift
              if(any(2*abs(idx)>=n))then
                if(matrix_free.and.local_error==0)write(error_unit,*)'HSEX_STERNHEIMER_VERTEX rank/k/q/G',rank,j,iq,idx
                local_error=1
              endif
              ! The charge vertex is diagonal and identical in both spinor components.
              do i=0,nreal-1
                point3=dble([mod(i,n(1)),mod(i/n(1),n(2)),i/(n(1)*n(2))])/n
                phase(i+1,ig,l)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(point3*modulo(idx,n)),8))
              enddo
              do component=2,ncomponent
                phase((component-1)*nreal+1:component*nreal,ig,l)=phase(:nreal,ig,l)
              enddo
            enddo
          enddo
          call record_phase(5,clock_phase)
          if(local_error==0)then
            call system_clock(clock_phase)
            do l=1,batch
              i=counts(l)
              allocate(vertex(ng,i),band_response(i,i))
              response(:,:,:,l)=0d0
              do band=1,nb
                ! Integer scalar screening uses Q onto the unoccupied
                ! complement, like the dense static reference. Empty source
                ! bands and approximate empty target eigenpairs are unused.
                if(ncomponent==1.and.occ(band,l)==0d0)cycle
                ! Empty source bands still contribute to the explicit occupation difference.
                do a=1,i
                  vertex(:,a)=phase(:,a,l)*psi(:,band,l)
                enddo
                ! The shared solver uses one-spin occupations (0..1).
                ! Restore scalar spin degeneracy exactly once in chi.
                call sternheimer_band(apply_at_target,states(:,:project_bands,j-orbital_offset), &
                  levels(:project_bands,j),pop(:project_bands,j)/spin_weight, &
                  lev(band,l),occ(band,l)/spin_weight,solver_temperature,vertex,band_response,err)
                if(err/=0)then
                  write(error_unit,*)'HSEX_STERNHEIMER_SOLVE rank/k/q/band/status',rank,j,beg+l-1,band,err
                  local_error=1
                  exit
                endif
                response(:i,:i,1,l)=response(:i,:i,1,l)+spin_weight*band_response
              enddo
              ! The static solution has no eta dependence. The downstream
              ! metallic closure currently expects two response channels.
              response(:,:,:,l)=spread(response(:,:,1,l),3,nc)
              deallocate(vertex,band_response)
              if(local_error/=0)exit
            enddo
            call record_phase(6,clock_phase)
            do iq=beg,last
              l=iq-beg+1;i=counts(l);blocks(iq)%chi=blocks(iq)%chi+response(:i,:i,:,l)
            enddo
          endif
        endif
        if(distributed)call comm_get_max(local_error,icomm)
        deallocate(psi,phase,lev,occ,counts,response,residual,istat);beg=last+1
        if(local_error/=0)exit
      enddo
      if(local_error/=0)exit
      if(matrix_free.and.comm_is_root(rank).and.j<=ik_e)then
        write(output_unit,'(a,i0,a,i0)')'HSEX_STERNHEIMER_K completed=',j-ik_s+1,' local_total=',ik_e-ik_s+1
        flush(output_unit)
      endif
    enddo
    if(allocated(h))deallocate(h,adj,d,e,u,tau,work)
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    call system_clock(clock_response)
    timing(1)=rmax;call comm_get_max(timing,timing_max,1,icomm);max_residual=timing_max(1)
    ! Reduce one q block at a time to the common q/-q owner.
    call system_clock(clock_phase)
    do iq=1,nq
      if(.not.active_q(iq))cycle
      mg=size(blocks(iq)%g,2);allocate(reduced(mg,mg,nc))
      call comm_summation(blocks(iq)%chi,reduced,size(reduced),icomm,blocks(iq)%owner)
      if(rank==blocks(iq)%owner)then
        blocks(iq)%chi=reduced
      else
        deallocate(blocks(iq)%chi)
      endif
      deallocate(reduced)
    enddo
    call record_phase(7,clock_phase)
    call system_clock(clock_phase)
    stage=60
    heads=0d0;qheads=0d0
    do iq=1,nq
      if(.not.active_q(iq))cycle
      if(rank/=blocks(iq)%owner)cycle
      jq=blocks(iq)%minus;mg=size(blocks(iq)%g,2)
      allocate(reduced(mg,mg,nc),Qcoords(3,mg),blocks(iq)%box(na,na,nc))
      do c=1,nc
        do b=1,mg
          do a=1,mg
            reduced(a,b,c)=blocks(iq)%chi(a,b,c)+ &
              conjg(blocks(jq)%chi(blocks(iq)%reverse(a),blocks(iq)%reverse(b),c))
          enddo
        enddo
        do a=1,mg
          Qcoords(:,a)=blocks(iq)%g(:,a)+dble(qsites(:,iq))/qmesh
        enddo
        call screen_assemble(n,boxes,Qcoords,rec,vol,nk,reduced(:,:,c),pump, &
          blocks(iq)%box(:,:,c),head,qhead,err)
        local_error=max(local_error,err)
      enddo
      heads(iq)=head;qheads(iq)=qhead
      if(pump.and.(insulating.or.iq/=1))then
        call screen_channels(blocks(iq)%box(:,:,1:1),blocks(iq)%box(:,:,2:2),error,.false.)
        if(error>eta_tol)local_error=1
      endif
      deallocate(reduced,Qcoords)
    enddo
    do iq=1,nq
      if(allocated(blocks(iq)%chi))deallocate(blocks(iq)%chi)
    enddo
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    block
      real(8)::total(nq)
      call comm_summation(heads,total,nq,icomm);heads=total
      call comm_summation(qheads,total,nq,icomm);qheads=total
    end block
    if(comm_is_root(rank))then
      write(output_unit,'(a)')'HSEX_W_FINITE_Q_COLUMNS iq qfrac(1:3) min_abs_q_plus_G[bohr^-1] inverse_head'
      do iq=1,nq
        if(qheads(iq)<=0d0)cycle
        ! HEADS includes local fields; QHEADS is the minimum |q+G| squared.
        write(output_unit,'(a,1x,i0,5(1x,es24.16))')'HSEX_W_FINITE_Q',iq, &
          dble(qsites(:,iq))/qmesh,sqrt(qheads(iq)),heads(iq)
      enddo
    endif
    qmin=minval(qheads,mask=qheads>0d0);einv=0d0;i=0
    do iq=1,nq
      if(abs(qheads(iq)-qmin)<=1d-8+1d-5*abs(qmin))then
        einv=einv+heads(iq);i=i+1
      endif
    enddo
    if(i==0)then
      status=61
      if(comm_is_root(rank))write(*,'(a)')'HSEX_W_REJECT: no finite-q dielectric head'
      return
    endif
    einv=einv/i
    allocate(point(na,na),origin(na,na));point=0d0;origin=0d0
    stage=70
    if(optical)then
      call comm_summation(optical_chi,optical_total,size(optical_chi),icomm)
      do b=1,size(optical_chi,2)
        do a=1,size(optical_chi,1)
          optical_chi(a,b)=optical_total(a,b)+merge(-1d0,1d0,(a<=3).neqv.(b<=3))* &
            conjg(optical_total(optical_reverse(a),optical_reverse(b)))
        enddo
      enddo
      if(comm_is_root(rank))then
        call longwave_cell(n,boxes,optical_g,rec,vol,nk,optical_chi,macro,point,origin,local_error)
        if(local_error==0)then
          write(output_unit,'(a,3es18.9)')'HSEX_W_EPS_LONGWAVE diagonal ',(macro(a,a),a=1,3)
          write(output_unit,'(a,es18.9)')'HSEX_W_EPS_FINITE_Q ',1d0/einv
        else
          write(error_unit,*)'HSEX_LONGWAVE_REJECT status=',local_error
        endif
      endif
      call comm_get_max(local_error,icomm)
      if(local_error/=0)then
        status=stage
        return
      endif
      call comm_bcast(macro,icomm,0)
      call comm_bcast(point,icomm,0)
      call comm_bcast(origin,icomm,0)
      if(rank==blocks(1)%owner)blocks(1)%box(:,:,nc)=point
    endif
    if(.not.insulating)then
      ! Gather only the finite-q head stencil, never the full Nq*Na^2 kernel.
      if(comm_is_root(rank))allocate(head_samples(na,na,nsel,2))
      allocate(buf(na,na,nc))
      do i=1,nsel
        iq=samples(i);buf=0d0
        if(rank==blocks(iq)%owner)buf=blocks(iq)%box
        call comm_bcast(buf,icomm,blocks(iq)%owner)
        if(comm_is_root(rank))head_samples(:,:,i,:)=buf
      enddo
      deallocate(buf)
      if(comm_is_root(rank))then
        if(matrix_free)then
          call save_head_samples(err)
          if(err/=0)local_error=8
        endif
        allocate(coarse_point(na,na),coarse_origin(na,na))
        call screen_metal_head(head_samples(:,:,:,2),qmesh,rec,mesh,head_tol,point,origin,err,samples(:nsel))
        if(err/=0)write(error_unit,*)'HSEX_HEAD_REJECT fine/status',err
        local_error=max(local_error,err)
        if(local_error==0.and..not.matrix_free)then
          call screen_metal_head(head_samples(:,:,:,1),qmesh,rec,mesh,head_tol, &
            coarse_point,coarse_origin,err,samples(:nsel))
          if(err/=0)write(error_unit,*)'HSEX_HEAD_REJECT coarse/status',err
          local_error=max(local_error,err)
          if(err==0)then
            call screen_channels(reshape(coarse_origin,[na,na,1]),reshape(origin,[na,na,1]),error,.false.)
            if(error>eta_tol)then
              write(error_unit,*)'HSEX_HEAD_REJECT eta comparison error/tolerance',error,eta_tol
              local_error=1
            endif
          endif
        endif
        deallocate(head_samples,coarse_point,coarse_origin)
      endif
      err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
      if(only_head)then
        status=0
        return
      endif
      call comm_bcast(point,icomm,0);call comm_bcast(origin,icomm,0)
      if(rank==blocks(1)%owner)blocks(1)%box(:,:,nc)=point
    endif
    stage=80
    call system_clock(clock_assembly)
    call record_phase(8,clock_phase)
    call system_clock(clock_phase)
    allocate(coarse(nlocal,nq),values(nlocal,nk),matrix(na,na))
    do iq=1,nq
      matrix=0d0
      if(rank==blocks(iq)%owner)matrix=blocks(iq)%box(:,:,nc)
      call comm_bcast(matrix,icomm,blocks(iq)%owner)
      do i=1,nlocal
        ia=mod(first_pair+i-2,na)+1;ib=(first_pair+i-2)/na+1
        coarse(i,iq)=matrix(ia,ib)
      enddo
      if(allocated(blocks(iq)%box))deallocate(blocks(iq)%box)
    enddo
    if(all(qmesh==mesh))then
      values=coarse
    else
      allocate(tail(nlocal,1))
      do iq=1,nq
        call long_range(dble(qsites(:,iq))/qmesh,tail(:,1));coarse(:,iq)=coarse(:,iq)-tail(:,1)
      enddo
      do iq=1,nk
        point3=dble(sites(:,iq))*qmesh/mesh;lo=floor(point3);fraction=point3-lo
        call long_range(dble(sites(:,iq))/mesh,values(:,iq))
        do dz=0,1
          do dy=0,1
            do dx=0,1
              idx=[dx,dy,dz];weight=product(merge(fraction,1d0-fraction,idx==1))
              if(weight==0d0)cycle
              vtx=modulo(lo+idx,qmesh);jq=1+vtx(1)+qmesh(1)*(vtx(2)+qmesh(2)*vtx(3))
              values(:,iq)=values(:,iq)+weight*coarse(:,jq)
            enddo
          enddo
        enddo
      enddo
      deallocate(tail)
    endif
    if(.not.insulating.or.optical)then
      do i=1,nlocal
        ia=mod(first_pair+i-2,na)+1;ib=(first_pair+i-2)/na+1;values(i,1)=origin(ia,ib)
      enddo
    endif
    call record_phase(9,clock_phase)
    call system_clock(clock_phase)
    call screen_transform(values,mesh)
    call record_phase(10,clock_phase)
    call system_clock(clock_phase)
    if(insulating.and..not.optical)then
      qc=(6d0*acos(-1d0)**2/(vol*nk))**(1d0/3d0)
      values=values+12d0*acos(-1d0)*einv/(qc*qc*vol*nk)
    endif
    if(.not.screen_finite(values))local_error=1
    if(size(values)>0)then
      if(maxval(abs(aimag(values)))>1d-7)local_error=1
    endif
    ! R/-R transpose reciprocity is checked with one block per rank.
    allocate(current(na,na),partner(na,na))
    do iq=1,nk
      matrix=0d0
      do i=1,nlocal
        ia=mod(first_pair+i-2,na)+1;ib=(first_pair+i-2)/na+1;matrix(ia,ib)=values(i,iq)
      enddo
      call comm_summation(matrix,current,na*na,icomm)
      idx=modulo(-sites(:,iq),mesh);jq=1+idx(1)+mesh(1)*(idx(2)+mesh(2)*idx(3));matrix=0d0
      do i=1,nlocal
        ia=mod(first_pair+i-2,na)+1;ib=(first_pair+i-2)/na+1;matrix(ia,ib)=values(i,jq)
      enddo
      call comm_summation(matrix,partner,na*na,icomm)
      if(maxval(abs(current-transpose(partner)))>1d-8)local_error=1
    enddo
    err=local_error;call comm_get_max(err,icomm)
    if(err/=0)then
      status=stage
      if(comm_is_root(rank))write(*,'(a,i0,a,i0)')'HSEX_W_REJECT stage=',stage,' error=',err
      return
    endif
    call record_phase(11,clock_phase)
    w=real(values,8);status=0
    call system_clock(clock_end)
    timing=dble([clock_response-clock_start,clock_assembly-clock_response,clock_end-clock_assembly])/clock_rate
    call comm_get_max(timing,timing_max,3,icomm)
    call comm_get_max(phase_seconds,phase_max,11,icomm)
    if(comm_is_root(rank))then
      write(*,'(a,3f12.4)')'HSEX_W_SECONDS response/assembly/transform: ',timing_max
      write(*,'(a,11es14.6)')'HSEX_W_PHASE_SECONDS build/factor/optical/exchange/vertices/solve/'// &
        'reduce/assemble/redistribute/fft/reciprocity: ',phase_max
      write(*,'(a,i0,a,es12.4)')'HSEX_W_OMP team=',nthreads,' maximum residual=',max_residual
      if(.not.matrix_free)write(*,'(a,i0,a,i0)')'HSEX_W_Q_PARALLEL regions=',q_regions,' batches=',q_batches
    endif
  contains
    subroutine record_phase(index,started)
      implicit none
      integer,intent(in) :: index
      integer(8),intent(in) :: started
      integer(8) :: finished
      call system_clock(finished)
      phase_seconds(index)=phase_seconds(index)+dble(finished-started)/dble(clock_rate)
    end subroutine record_phase

    subroutine compressed_target()
      implicit none
      type(s_response_stars) :: stars
      type(s_response_space_stats) :: statistics
      integer,allocatable :: all_g(:,:,:),all_counts(:),sources(:),representatives(:),batch_q(:),offsets(:)
      complex(8),allocatable :: rhs(:,:),solutions(:,:),wave(:),reps_chi(:,:,:),local_chi(:,:)
      real(8),allocatable :: qfrac(:,:),energies(:),checks(:)
      integer :: maximum,code,target_q,representative,rep_count,entry,finish,nrhs,p,source_band,g_index,point_index
      integer :: q_index,column,capacity_batches,capped,actions,validation,grid_wrap(3),q_column
      real(8) :: position(3),worst
      integer(8) :: started
      call system_clock(started)
      maximum=0
      do target_q=1,nq
        maximum=max(maximum,size(blocks(target_q)%g,2))
      enddo
      allocate(all_g(3,maximum,nq),all_counts(nq),sources(nq),representatives(nq), &
        qfrac(3,nq),batch_q(compressed_batch),offsets(compressed_batch), &
        wave(ng),reps_chi(maximum,maximum,nq),local_chi(maximum,maximum),stat=code)
      if(code/=0)then
        local_error=1
        return
      endif
      all_g=0;reps_chi=0d0
      do target_q=1,nq
        all_counts(target_q)=size(blocks(target_q)%g,2)
        all_g(:,:all_counts(target_q),target_q)=blocks(target_q)%g
        sources(target_q)=source_k(j,target_q)
        qfrac(:,target_q)=dble(qsites(:,target_q))/qmesh
      enddo
      if(compressed_symmetry)then
        call response_q_stars(n,mesh,qmesh,cell,rec,atoms,species,frac,levels,occupied_all,j, &
          qfrac,all_g,all_counts,sources,apply_at_target,stars,code)
        if(code/=0)then
          local_error=code
          return
        endif
      else
        allocate(stars%representative(nq),stars%operation(nq),stars%column(maximum,nq),stat=code)
        if(code/=0)then
          local_error=1
          return
        endif
        stars%nsym=1;stars%nrep=nq;stars%column=0
        do target_q=1,nq
          stars%representative(target_q)=target_q;stars%operation(target_q)=1
          do g_index=1,all_counts(target_q)
            stars%column(g_index,target_q)=g_index
          enddo
        enddo
      endif
      rep_count=0
      do target_q=1,nq
        if(stars%representative(target_q)/=target_q)cycle
        rep_count=rep_count+1;representatives(rep_count)=target_q
      enddo
      entry=1;capacity_batches=0;capped=0;actions=stars%h_probes;validation=0;worst=0d0
      do while(entry<=rep_count)
        finish=min(rep_count,entry+compressed_batch-1);nrhs=0
        do column=1,finish-entry+1
          q_index=representatives(entry+column-1);batch_q(column)=q_index
          offsets(column)=nrhs+1;nrhs=nrhs+nv*all_counts(q_index)
        enddo
        allocate(rhs(ng,nrhs),solutions(ng,nrhs),energies(nrhs),checks(nrhs),stat=code)
        if(code/=0)then
          local_error=1
          return
        endif
        do column=1,finish-entry+1
          q_index=batch_q(column)
          grid_wrap=nint(frac(:,sources(q_index))+qfrac(:,q_index)-frac(:,j))
          do g_index=1,all_counts(q_index)
            if(any(2*abs(all_g(:,g_index,q_index)+grid_wrap)>=n))then
              local_error=1
              return
            endif
            do point_index=0,nreal-1
              position=dble([mod(point_index,n(1)),mod(point_index/n(1),n(2)),point_index/(n(1)*n(2))])/n
              wave(point_index+1)=exp(cmplx(0d0,2d0*acos(-1d0)* &
                sum(position*(all_g(:,g_index,q_index)+grid_wrap)),8))
            enddo
            do source_band=1,nv
              p=offsets(column)+(source_band-1)*all_counts(q_index)+g_index-1
              rhs(:,p)=wave*occupied_all(:,source_band,sources(q_index))
              energies(p)=levels(source_band,sources(q_index))
            enddo
          enddo
        enddo
        call response_space_solve(apply_at_target,occupied_all(:,:,j),energies,rhs,compressed_capacity, &
          compressed_tolerance,solutions,checks,statistics,code)
        if(code/=0.and.code/=2)then
          local_error=code
          return
        endif
        if(.not.screen_finite(solutions).or..not.screen_finite(checks))then
          local_error=1
          return
        endif
        capacity_batches=capacity_batches+1
        if(code==2)capped=capped+1
        actions=actions+statistics%h_actions;validation=validation+statistics%validation_actions
        worst=max(worst,maxval(checks))
        do column=1,finish-entry+1
          q_index=batch_q(column);q_column=all_counts(q_index)
          do source_band=1,nv
            p=offsets(column)+(source_band-1)*q_column
            call zgemm('C','N',q_column,q_column,ng,(1d0,0d0),rhs(:,p:),ng,solutions(:,p:),ng, &
              (1d0,0d0),reps_chi(:,:,q_index),maximum)
          enddo
        enddo
        deallocate(rhs,solutions,energies,checks);entry=finish+1
      enddo
      do target_q=1,nq
        q_column=all_counts(target_q);representative=stars%representative(target_q)
        local_chi(:q_column,:q_column)=reps_chi(stars%column(:q_column,target_q), &
          stars%column(:q_column,target_q),representative)
        if(.not.screen_finite(local_chi(:q_column,:q_column)))then
          local_error=1
          return
        endif
        ! sternheimer_band contributes 0.5 per one-spin occupied source;
        ! scalar spin_weight=2 restores exactly one copy here. The existing
        ! downstream k normalization and +/-q assembly remain unchanged.
        blocks(target_q)%chi(:q_column,:q_column,1)= &
          blocks(target_q)%chi(:q_column,:q_column,1)+local_chi(:q_column,:q_column)
      enddo
      rmax=max(rmax,worst)
      write(output_unit,'(a,8i9,es15.6)')'HSEX_REDUCED_W_K rank/k/ops/q/rep_q/batches/unfinished/Hpsi/residual ', &
        rank,j,stars%nsym,nq,rep_count,capacity_batches,capped,actions,worst
      flush(output_unit)
      call record_phase(6,started)
    end subroutine compressed_target

    subroutine dense_q_batches()
!$    use omp_lib, only: omp_get_num_threads
      implicit none
      type(hsex_q_shared) :: shared
      type(hsex_q_workspace) :: workspace
      integer :: max_g,entry,code,release_error
      integer(8) :: started,finished
      real(8) :: vertex_seconds,solve_seconds
      logical :: done
      max_g=0
      do entry=1,nq
        if(active_q(entry))max_g=max(max_g,size(blocks(entry)%g,2))
      enddo
      if(j<=ik_e)then
        call hsex_q_prepare(ng,nb,max_g,nu,merge(1,0,pump),merge(1,0,pump),etas,original_h,shared,code)
        if(code/=0)local_error=1
      endif
      if(distributed)call comm_get_max(local_error,icomm)
      if(local_error/=0)then
        if(j<=ik_e)call hsex_q_release(shared)
        return
      endif
      beg=1;done=.false.;release_error=0;vertex_seconds=0d0;solve_seconds=0d0;team=1
      q_regions=q_regions+1
      !$omp parallel num_threads(nthreads) default(shared) private(workspace,entry,code,started,finished) &
      !$omp& reduction(max:release_error,vertex_seconds,solve_seconds)
      vertex_seconds=0d0;solve_seconds=0d0;release_error=0
      if(j<=ik_e)call hsex_q_worker_prepare(ng,nb,max_g,nu,merge(1,0,pump),merge(1,0,pump),shared,workspace)
      !$omp master
!$    team=omp_get_num_threads()
      if(team/=nthreads)local_error=1
      if(distributed)call comm_get_max(local_error,icomm)
      !$omp end master
      do
        !$omp master
        do while(beg<=nq)
          if(active_q(beg))exit
          beg=beg+1
        enddo
        done=beg>nq.or.local_error/=0
        if(.not.done)then
          mg=0;last=beg-1
          do iq=beg,min(nq,beg+2*batch_threads-1)
            if(.not.active_q(iq))exit
            mq=max(mg,size(blocks(iq)%g,2))
            if(iq>beg.and.16d0*ng*(nb+mq)*(iq-beg+1)>64d0*1024**2)exit
            mg=mq;last=iq
          enddo
          batch=last-beg+1
          allocate(psi(ng,nb,batch),phase(ng,mg,batch),lev(nb,batch),occ(nb,batch),counts(batch), &
            response(mg,mg,nc,batch),residual(batch),istat(batch))
          psi=0d0;phase=0d0;occ=0d0;response=0d0;residual=0d0;istat=0
          if(distributed)then
            call system_clock(clock_phase)
            call exchange_source_batch(code)
            if(code/=0)local_error=1
            call record_phase(4,clock_phase)
            call comm_get_max(local_error,icomm)
          endif
          done=local_error/=0
        endif
        !$omp end master
        !$omp barrier
        if(done)exit
        !$omp do schedule(dynamic,1)
        do entry=1,batch
          if(j>ik_e)cycle
          call system_clock(started)
          call dense_vertex(entry,code)
          call system_clock(finished)
          vertex_seconds=vertex_seconds+dble(finished-started)/clock_rate
          if(code/=0)then
            istat(entry)=1
            cycle
          endif
          call system_clock(started)
          call hsex_q_one(ng,nb,mg,nu,nc,merge(1,0,pump),merge(1,0,pump),adj,d,e,u,psi(:,:,entry), &
            phase(:,:,entry),counts(entry),lev(:,entry),occ(:,entry),etas,shared,workspace, &
            response(:,:,:,entry),residual(entry),istat(entry))
          call system_clock(finished)
          solve_seconds=solve_seconds+dble(finished-started)/clock_rate
        enddo
        !$omp end do
        !$omp master
        if(any(istat/=0))local_error=1
        rmax=max(rmax,maxval(residual))
        if(j<=ik_e.and.local_error==0)then
          do iq=beg,last
            l=iq-beg+1;i=counts(l)
            blocks(iq)%chi=blocks(iq)%chi+response(:i,:i,:,l)
          enddo
        endif
        if(distributed)call comm_get_max(local_error,icomm)
        deallocate(psi,phase,lev,occ,counts,response,residual,istat)
        q_batches=q_batches+1;beg=last+1;done=local_error/=0
        !$omp end master
        ! Workers next wait at the top barrier; only the master touches
        ! shared batch state before it. No extra bottom barrier is needed.
      enddo
      if(j<=ik_e)then
        call hsex_q_worker_release(workspace,code)
        release_error=max(release_error,abs(code))
      endif
      !$omp end parallel
      if(allocated(psi))deallocate(psi,phase,lev,occ,counts,response,residual,istat)
      if(j<=ik_e)call hsex_q_release(shared)
      local_error=max(local_error,release_error)
      phase_seconds(5)=phase_seconds(5)+vertex_seconds
      phase_seconds(6)=phase_seconds(6)+solve_seconds
    end subroutine dense_q_batches

    recursive subroutine dense_vertex(entry,code)
      implicit none
      integer,intent(in) :: entry
      integer,intent(out) :: code
      integer :: q_index,source_index,shift_q(3),idx_q(3),g_index,grid_index,spin
      real(8) :: location(3)
      code=0;q_index=beg+entry-1
      source_index=source_k(j,q_index)
      if(.not.distributed)psi(:,:,entry)=states(:,:,source_index)
      lev(:,entry)=levels(:,source_index);occ(:,entry)=pop(:,source_index)
      if(.not.pump)occ(nv+1:,entry)=0d0
      shift_q=nint(frac(:,source_index)+dble(qsites(:,q_index))/qmesh-frac(:,j))
      counts(entry)=size(blocks(q_index)%g,2)
      do g_index=1,counts(entry)
        idx_q=blocks(q_index)%g(:,g_index)+shift_q
        if(any(2*abs(idx_q)>=n))then
          code=1
          return
        endif
        do grid_index=0,nreal-1
          location=dble([mod(grid_index,n(1)),mod(grid_index/n(1),n(2)),grid_index/(n(1)*n(2))])/n
          phase(grid_index+1,g_index,entry)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(location*modulo(idx_q,n)),8))
        enddo
        do spin=2,ncomponent
          phase((spin-1)*nreal+1:spin*nreal,g_index,entry)=phase(:nreal,g_index,entry)
        enddo
      enddo
    end subroutine dense_vertex

    subroutine prepare_dense_target()
      implicit none
      integer(8) :: started
      real(8) :: defect
      call system_clock(started)
      call load_h(j,h,err)
      call record_phase(1,started)
      if(err/=0)then
        write(error_unit,*)'HSEX_DENSE_TARGET_REJECT load rank/k/status',rank,j,err
        local_error=1;return
      endif
      if(.not.screen_finite(h))then
        write(error_unit,*)'HSEX_DENSE_TARGET_REJECT nonfinite rank/k',rank,j
        local_error=1;return
      endif
      defect=maxval(abs(h-transpose(conjg(h))))
      if(defect>1d-9)then
        write(error_unit,*)'HSEX_DENSE_TARGET_REJECT Hermitian rank/k/error',rank,j,defect
        local_error=1;return
      endif
      if(optical)then
        if(allocated(optical_h))optical_h=h
        call accumulate_optical()
        if(local_error/=0)return
      endif
      call system_clock(started)
      allocate(hu(ng,merge(nb,nv,pump)),gram(merge(nb,nv,pump),merge(nb,nv,pump)))
      i=size(hu,2)
      call zgemm('N','N',ng,i,ng,cmplx(1d0,0d0,8),h,ng,states(:,:,j-orbital_offset),ng,cmplx(0d0,0d0,8),hu,ng)
      call zgemm('C','N',i,i,ng,cmplx(1d0,0d0,8),states(:,:,j-orbital_offset),ng,states(:,:,j-orbital_offset),ng, &
        cmplx(0d0,0d0,8),gram,i)
      do a=1,i
        gram(a,a)=gram(a,a)-1d0
        defect=sqrt(sum(abs(hu(:,a)-states(:,a,j-orbital_offset)*levels(a,j))**2))
        if(defect>1d-4)then
          write(error_unit,*)'HSEX_DENSE_TARGET_REJECT eigenpair rank/k/band/residual',rank,j,a,defect
          local_error=1
        endif
      enddo
      defect=maxval(abs(gram))
      if(defect>1d-8)then
        write(error_unit,*)'HSEX_DENSE_TARGET_REJECT orthogonality rank/k/error',rank,j,defect
        local_error=1
      endif
      if(.not.screen_finite(gram))local_error=1
      if(.not.screen_finite(hu))local_error=1
      if(.not.pump)then
        u=states(:,:nv,j-orbital_offset);allocate(huq(ng,nv));gram=0d0
        call zgemm('C','N',nv,nv,ng,cmplx(1d0,0d0,8),u,ng,hu,ng,cmplx(0d0,0d0,8),gram,nv)
        do a=1,nv
          gram(a,a)=gram(a,a)+maxval(levels(:nv,j))+2d0
        enddo
        huq=hu-.5d0*matmul(u,gram)
        call zgemm('N','C',ng,ng,nv,cmplx(-1d0,0d0,8),huq,ng,u,ng,cmplx(1d0,0d0,8),h,ng)
        call zgemm('N','C',ng,ng,nv,cmplx(-1d0,0d0,8),u,ng,huq,ng,cmplx(1d0,0d0,8),h,ng)
        deallocate(huq)
      endif
      deallocate(hu,gram)
      if(local_error/=0)return
      if(pump)then
        ! Preserve the Hermitian operator selected by ZHETRD('L'), including
        ! the real diagonal; exported HSP2 snapshots use the same convention.
        do a=1,ng
          h(a,a)=cmplx(real(h(a,a),8),0d0,8)
          do b=a+1,ng
            h(a,b)=conjg(h(b,a))
          enddo
        enddo
        original_h=h
      endif
      call zhetrd('L',ng,h,ng,d,e,tau,work,size(work),err)
      if(err/=0)then
        local_error=1;return
      endif
      adj=0d0;adj(1,1)=1d0;adj(2:,2:)=h(2:,:ng-1)
      call zungqr(ng-1,ng-1,ng-1,adj(2,2),ng,tau,work,size(work),err)
      if(err/=0)then
        local_error=1;return
      endif
      adj=transpose(conjg(adj))
      call record_phase(2,started)
    end subroutine prepare_dense_target

    subroutine apply_longwave_target(vector,action,ierr)
      implicit none
      complex(8),intent(in) :: vector(:)
      complex(8),intent(out) :: action(:)
      integer,intent(out) :: ierr
      if(matrix_free)then
        call apply_h(j,vector,action,ierr)
      else if(present(apply_optical_h))then
        call apply_optical_h(j,vector,action,ierr)
      else
        action=matmul(optical_h,vector)
        ierr=0
      endif
    end subroutine

    subroutine accumulate_optical()
      implicit none
      complex(8) :: body(ng,size(optical_g,2)),velocity(ng,3),response(size(optical_chi,1),size(optical_chi,2))
      integer :: band,g,site,comp,code
      integer(8) :: started
      real(8) :: pos(3)
      call system_clock(started)
      do band=1,nv
        call apply_velocity(j,states(:,band,j-orbital_offset),velocity,code)
        if(code/=0)then
          local_error=1;return
        endif
        do g=1,size(optical_g,2)
          do comp=1,ncomponent
            do site=0,nreal-1
              pos=dble([mod(site,n(1)),mod(site/n(1),n(2)),site/(n(1)*n(2))])/n
              body(site+1+(comp-1)*nreal,g)=states(site+1+(comp-1)*nreal,band,j-orbital_offset)* &
                exp(cmplx(0d0,2d0*acos(-1d0)*sum(pos*optical_g(:,g)),8))
            enddo
          enddo
        enddo
        call longwave_band(apply_longwave_target,states(:,:nv,j-orbital_offset),levels(band,j),velocity,body,response,code)
        if(code/=0)then
          write(error_unit,*)'HSEX_OPTICAL_REJECT rank/k/band/status',rank,j,band,code
          local_error=1;return
        endif
        optical_chi=optical_chi+.5d0*pop(band,j)*response
      enddo
      call record_phase(3,started)
    end subroutine

    recursive integer function source_k(destination,q_index) result(source_index)
      implicit none
      integer,intent(in) :: destination,q_index
      integer :: coord(3)
      coord=modulo(order(:,destination)-qsites(:,q_index)*(mesh/qmesh),mesh)
      source_index=lookup(1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3)))
    end function source_k

    subroutine exchange_source_batch(status)
      implicit none
      integer,intent(out) :: status
      integer :: send_counts(nproc),recv_counts(nproc),send_offsets(nproc),recv_offsets(nproc)
      integer :: destinations(batch),positions(nproc)
      integer :: p,kdest,ksrc,q_index,entry,owner_rank,block_size,offset,astat
      integer(8) :: send_size,recv_size,orbital_size
      complex(8),allocatable :: send(:),recv(:)

      send_counts=0
      recv_counts=0
      destinations=0
      do p=1,nproc
        kdest=ranges(1,p)+k_slot-1
        if(kdest>ranges(2,p))cycle
        do q_index=beg,last
          ksrc=source_k(kdest,q_index)
          owner_rank=k_owner(ksrc)
          if(owner_rank==rank+1.and.p/=rank+1)send_counts(p)=send_counts(p)+1
          if(p==rank+1)then
            entry=q_index-beg+1
            if(owner_rank==rank+1)then
              psi(:,:,entry)=states(:,:,ksrc-orbital_offset)
            else
              destinations(entry)=owner_rank
              recv_counts(owner_rank)=recv_counts(owner_rank)+1
            endif
          endif
        enddo
      enddo
      orbital_size=int(ng,8)*nb
      send_size=sum(int(send_counts,8))*orbital_size
      recv_size=sum(int(recv_counts,8))*orbital_size
      status=0
      if(max(send_size,recv_size,orbital_size)>huge(block_size))status=1
      call comm_get_max(status,icomm)
      if(status/=0)return
      block_size=int(orbital_size)
      send_counts=send_counts*block_size
      recv_counts=recv_counts*block_size
      do p=1,nproc
        send_offsets(p)=sum(send_counts(:p-1))
        recv_offsets(p)=sum(recv_counts(:p-1))
      enddo
      allocate(send(max(1,int(send_size))),recv(max(1,int(recv_size))),stat=astat)
      status=merge(1,0,astat/=0)
      call comm_get_max(status,icomm)
      if(status/=0)return
      do p=1,nproc
        kdest=ranges(1,p)+k_slot-1
        if(kdest>ranges(2,p))cycle
        if(p==rank+1)cycle
        offset=send_offsets(p)
        do q_index=beg,last
          ksrc=source_k(kdest,q_index)
          if(k_owner(ksrc)/=rank+1)cycle
          send(offset+1:offset+block_size)=reshape(states(:,:,ksrc-orbital_offset),[block_size])
          offset=offset+block_size
        enddo
      enddo
      call comm_alltoallv_complex(send,send_counts,send_offsets,recv,recv_counts,recv_offsets,icomm)
      positions=recv_offsets
      do entry=1,batch
        owner_rank=destinations(entry)
        if(owner_rank==0)cycle
        offset=positions(owner_rank)
        psi(:,:,entry)=reshape(recv(offset+1:offset+block_size),[ng,nb])
        positions(owner_rank)=offset+block_size
      enddo
    end subroutine exchange_source_batch

    subroutine save_head_samples(ios)
      use, intrinsic :: iso_c_binding, only: c_int,c_char,c_null_char
      implicit none
      integer,intent(out) :: ios
      integer :: unit,close_status
      character(*),parameter :: path='hsex-spinor-head-samples.bin', partial=path//'.partial'
      interface
        function link_file(old,new) bind(C,name='link') result(code)
          import c_int,c_char
          implicit none
          character(c_char),intent(in) :: old(*),new(*)
          integer(c_int) :: code
        end function link_file
        function unlink_file(name) bind(C,name='unlink') result(code)
          import c_int,c_char
          implicit none
          character(c_char),intent(in) :: name(*)
          integer(c_int) :: code
        end function unlink_file
      end interface
      open(newunit=unit,file=partial,status='new',access='stream',form='unformatted',iostat=ios)
      if(ios==0)then
        write(unit,iostat=ios)'HXQ1',qmesh,mesh,na,nsel,rec,head_tol,samples(:nsel),head_samples(:,:,:,2)
        close(unit,iostat=close_status)
        if(ios==0)ios=close_status
        if(ios==0)ios=link_file(partial//c_null_char,path//c_null_char)
        if(ios==0)ios=unlink_file(partial//c_null_char)
      endif
      if(ios/=0)write(error_unit,*)'HSEX_HEAD_SAMPLES_WRITE status',ios
    end subroutine save_head_samples

    subroutine apply_at_target(x,y,ierr)
      implicit none
      complex(8),intent(in)::x(:)
      complex(8),intent(out)::y(:)
      integer,intent(out)::ierr
      call apply_h(j,x,y,ierr)
    end subroutine

    subroutine long_range(qfrac,result)
      implicit none
      real(8),intent(in)::qfrac(3)
      complex(8),intent(out)::result(:)
      complex(8),allocatable::forms(:,:)
      real(8),allocatable::coords(:,:),bare(:)
      real(8)::norm2,cartesian_q(3)
      integer::g,count_g,p,left,right
      result=0d0
      if(.not.insulating)return
      count_g=0
      do g=1,allg
        norm2=sum(matmul(rec,gall(:,g)+qfrac)**2)
        if(norm2>1d-14.and.norm2<=2d0*cutoff)count_g=count_g+1
      enddo
      allocate(coords(3,count_g),bare(count_g),forms(na,count_g));count_g=0
      do g=1,allg
        norm2=sum(matmul(rec,gall(:,g)+qfrac)**2)
        if(norm2<=1d-14.or.norm2>2d0*cutoff)cycle
        count_g=count_g+1;coords(:,count_g)=gall(:,g)+qfrac
        if(optical)then
          cartesian_q=matmul(rec,coords(:,count_g))
          bare(count_g)=4d0*acos(-1d0)/(vol*dot_product(cartesian_q,matmul(macro,cartesian_q)))
        else
          bare(count_g)=4d0*acos(-1d0)*einv/(vol*norm2)
        endif
      enddo
      call screen_forms(n,boxes,coords,forms)
      do p=1,nlocal
        left=mod(first_pair+p-2,na)+1;right=(first_pair+p-2)/na+1
        result(p)=sum(conjg(forms(left,:))*forms(right,:)*bare)
      enddo
    end subroutine
  end subroutine
end module
