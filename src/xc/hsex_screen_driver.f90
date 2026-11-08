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
  end interface
contains
  subroutine screen_compute(n,mesh,qmesh,boxes,cutoff,cell,k,levels,pop,states,nv,eta,head_tol,eta_tol, &
      nthreads,ik_s,ik_e,icomm,rank,nproc,load_h,first_pair,nlocal,w,status,apply_h,temperature,head_only,distributed_states)
    use, intrinsic :: iso_fortran_env, only: error_unit,output_unit
    use communication, only: comm_summation,comm_bcast,comm_get_max,comm_is_root,comm_alltoallv_complex
    use hsex_screen_math, only: screen_sites,screen_geometry,screen_assemble,screen_population,screen_channels, &
      screen_metal_head,screen_forms,screen_transform,screen_finite
    use hsex_q_response, only: hsex_q_solve
    use hsex_sternheimer, only: sternheimer_band
    implicit none
    integer,intent(in)::n(3),mesh(3),qmesh(3),boxes,nv,nthreads,ik_s,ik_e,icomm,rank,nproc,first_pair,nlocal
    real(8),intent(in)::cutoff,cell(3,3),k(:,:),levels(:,:),pop(:,:),eta,head_tol,eta_tol
    complex(8),intent(in)::states(:,:,:)
    procedure(h_loader)::load_h
    procedure(h_action),optional::apply_h
    real(8),intent(in),optional::temperature
    complex(8),allocatable::vertex(:,:),band_response(:,:)
    logical,intent(in),optional::head_only,distributed_states
    logical,allocatable :: active_q(:)
    logical::matrix_free,only_head,distributed
    integer :: orbital_offset,local_k,k_slots,k_slot,peer,batch_threads
    integer,allocatable :: ranges(:,:),local_range(:,:),k_owner(:)
    integer::band
    real(8),intent(out)::w(nlocal,product(mesh))
    integer,intent(out)::status
    type(q_block),allocatable::blocks(:)
    integer,allocatable::sites(:,:),qsites(:,:),order(:,:),lookup(:),gall(:,:),counts(:),istat(:),samples(:)
    integer::stage
    integer::nk,nq,ng,nreal,ncomponent,component,nb,na,nc,igmax(3),allg,i,j,l,a,b,c,iq,jq,ig,jg,owner,local_error,err,team
    integer::idx(3),shift(3),source,beg,last,batch,mq,mg,nu,slot,nsel,ia,ib,dx,dy,dz,lo(3),vtx(3)
    real(8)::rec(3,3),vol,dv,q(3),q2,head,qhead,einv,qmin,error,rmax,radius,qc,weight,fraction(3),point3(3)
    real(8),allocatable::frac(:,:),d(:),e(:),lev(:,:),occ(:,:),etas(:),residual(:),heads(:),qheads(:),radii(:)
    complex(8),allocatable::h(:,:),original_h(:,:),adj(:,:),u(:,:),tau(:),work(:),hu(:,:),huq(:,:),gram(:,:)
    complex(8),allocatable::psi(:,:,:),phase(:,:,:),response(:,:,:,:),reduced(:,:,:),qvec(:,:),matrix(:,:)
    complex(8),allocatable::coarse(:,:),values(:,:),tail(:,:),origin(:,:),point(:,:),head_samples(:,:,:,:)
    complex(8),allocatable::buf(:,:,:),coarse_origin(:,:),coarse_point(:,:),current(:,:),partner(:,:)
    real(8),allocatable::Qcoords(:,:)
    logical::pump,insulating
    integer(8)::clock_start,clock_response,clock_assembly,clock_end,clock_rate
    real(8)::timing(3),timing_max(3),max_residual
    call system_clock(clock_start,clock_rate)
    stage=10
    status=1;w=0d0;nk=product(mesh);nq=product(qmesh);nreal=product(n);nb=size(levels,1);na=boxes**3
    ng=size(states,1)
    local_error=0
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
    only_head=.false.
    if(present(head_only))only_head=head_only
    if(only_head.and..not.matrix_free)local_error=1
    if(matrix_free)then
      if(.not.present(temperature).or.ncomponent/=2) local_error=1
      if(present(temperature))then
        if(.not.screen_finite([temperature]).or.temperature<=0d0) local_error=1
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
      endif
    endif
    if(.not.screen_finite(states))local_error=1
    if(.not.pump)then
      if(any(pop(:nv,:)/=2d0).or.any(pop(nv+1:,:)/=0d0))local_error=1
      insulating=.true.
    endif
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
        do a=1,nb
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
          call exchange_source_batch(err)
          if(err/=0)then
            local_error=err
            exit
          endif
        endif
        if(j<=ik_e)then
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
          if(local_error==0)then
            if(matrix_free)then
              do l=1,batch
                i=counts(l)
                allocate(vertex(ng,i),band_response(i,i))
                response(:,:,:,l)=0d0
                do band=1,nb
                  ! Empty source bands still contribute to the explicit occupation difference.
                  do a=1,i
                    vertex(:,a)=phase(:,a,l)*psi(:,band,l)
                  enddo
                  call sternheimer_band(apply_at_target,states(:,:,j-orbital_offset),levels(:,j),pop(:,j),lev(band,l),occ(band,l), &
                    temperature,vertex,band_response,err)
                  if(err/=0)then
                    write(error_unit,*)'HSEX_STERNHEIMER_SOLVE rank/k/q/band/status',rank,j,beg+l-1,band,err
                    local_error=1
                    exit
                  endif
                  response(:i,:i,1,l)=response(:i,:i,1,l)+band_response
                enddo
                ! The static solution has no eta dependence. The downstream
                ! metallic closure currently expects two response channels.
                response(:,:,:,l)=spread(response(:,:,1,l),3,nc)
                deallocate(vertex,band_response)
                if(local_error/=0)exit
              enddo
            else
              call hsex_q_solve(ng,nb,batch,mg,nu,nc,nthreads,merge(1,0,pump),adj,d,e,u,psi,phase,counts,lev,occ, &
                etas,response,residual,istat,team,merge(1,0,pump),original_h)
              if(any(istat/=0).or.team/=nthreads)local_error=1
              rmax=max(rmax,maxval(residual))
            endif
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
    if(.not.insulating)then
      do i=1,nlocal
        ia=mod(first_pair+i-2,na)+1;ib=(first_pair+i-2)/na+1;values(i,1)=origin(ia,ib)
      enddo
    endif
    call screen_transform(values,mesh)
    if(insulating)then
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
    w=real(values,8);status=0
    call system_clock(clock_end)
    timing=dble([clock_response-clock_start,clock_assembly-clock_response,clock_end-clock_assembly])/clock_rate
    call comm_get_max(timing,timing_max,3,icomm)
    if(comm_is_root(rank))then
      write(*,'(a,3f12.4)')'HSEX_W_SECONDS response/assembly/transform: ',timing_max
      write(*,'(a,i0,a,es12.4)')'HSEX_W_OMP team=',nthreads,' maximum residual=',max_residual
    endif
  contains
    subroutine prepare_dense_target()
      implicit none
      call load_h(j,h,err)
      if(err/=0)then
        local_error=1;return
      endif
      if(.not.screen_finite(h))then
        local_error=1;return
      endif
      if(maxval(abs(h-transpose(conjg(h))))>1d-9)then
        local_error=1;return
      endif
      allocate(hu(ng,merge(nb,nv,pump)),gram(merge(nb,nv,pump),merge(nb,nv,pump)))
      i=size(hu,2)
      call zgemm('N','N',ng,i,ng,cmplx(1d0,0d0,8),h,ng,states(:,:,j-orbital_offset),ng,cmplx(0d0,0d0,8),hu,ng)
      call zgemm('C','N',i,i,ng,cmplx(1d0,0d0,8),states(:,:,j-orbital_offset),ng,states(:,:,j-orbital_offset),ng, &
        cmplx(0d0,0d0,8),gram,i)
      do a=1,i
        gram(a,a)=gram(a,a)-1d0
        if(sqrt(sum(abs(hu(:,a)-states(:,a,j-orbital_offset)*levels(a,j))**2))>1d-4)local_error=1
      enddo
      if(maxval(abs(gram))>1d-8)local_error=1
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
    end subroutine prepare_dense_target

    integer function source_k(destination,q_index) result(source_index)
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
      real(8)::norm2
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
        bare(count_g)=4d0*acos(-1d0)*einv/(vol*norm2)
      enddo
      call screen_forms(n,boxes,coords,forms)
      do p=1,nlocal
        left=mod(first_pair+p-2,na)+1;right=(first_pair+p-2)/na+1
        result(p)=sum(conjg(forms(left,:))*forms(right,:)*bare)
      enddo
    end subroutine
  end subroutine
end module
