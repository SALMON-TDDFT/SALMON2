! Fixed-target-k diagnostic: source/q/G compression, never publishes W.
module hsex_reduced_q
  implicit none
  private
  public :: reduced_q_test
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
  subroutine reduced_q_test(n,mesh,qmesh,cutoff,cell,k,levels,states,nv,ks,ke,comm,rank, &
      destination,qcount,maxrank,tolerance,apply,status,use_symmetry,atoms,species)
    use hsex_screen_math, only: screen_geometry,screen_sites,screen_finite
    use hsex_sternheimer, only: sternheimer_solve,s_sternheimer_stats
    use hsex_response_space, only: response_space_solve,s_response_space_stats
    use hsex_response_symmetry, only: s_response_stars,response_q_stars
    use communication, only: comm_summation,comm_get_max
    implicit none
    integer,intent(in) :: n(3),mesh(3),qmesh(3),nv,ks,ke,comm,rank,destination,qcount,maxrank
    real(8),intent(in) :: cutoff,cell(3,3),k(:,:),levels(:,:),tolerance
    complex(8),intent(in) :: states(:,:,:)
    procedure(grid_action) :: apply
    integer,intent(out) :: status
    logical,optional,intent(in) :: use_symmetry
    real(8),optional,intent(in) :: atoms(:,:)
    integer,optional,intent(in) :: species(:)
    complex(8),allocatable :: local(:,:,:),occupied(:,:,:),rhs(:,:),reference(:,:),solution(:,:),phase(:)
    complex(8),allocatable :: projected(:,:),gram(:,:),work(:),reference_block(:,:),actual_block(:,:)
    real(8),allocatable :: frac(:,:),radius(:),energy(:),residual(:),spectrum(:),rwork(:)
    real(8) :: rec(3,3),vol,q(3),r(3),error,reference_norm,solution_error
    integer,allocatable :: sites(:,:),selected(:),lookup(:),order(:,:),counts(:),starts(:),sources(:),g(:,:,:)
    integer :: ng,nk,nq,igmax(3),i,j,a,b,c,iq,ig,source,wrap(3),idx(3),slot,code,astat,mg,mmax,nrhs,at,band
    complex(8),allocatable :: compact_rhs(:,:),ref_blocks(:,:,:),approx_blocks(:,:,:),restored(:,:)
    real(8),allocatable :: qfractions(:,:)
    integer,allocatable :: retained(:),compact_start(:)
    logical,allocatable :: representative_rhs(:)
    logical :: symmetry
    type(s_response_stars) :: stars
    integer :: p,info,min_dim,svd_rank(4),nret,jq,reference_h_before,sym_cg_h
    real(8) :: symmetry_error,max_symmetry_error,reference_time_before,sym_cg_seconds
    real(8),parameter :: thresholds(4)=[1d-2,1d-4,1d-6,1d-8]
    integer(8) :: started,prepared,ref_done,solved,spectral,restored_end,rate
    type(s_sternheimer_stats) :: cg
    type(s_response_space_stats) :: reduced
    status=1;ng=product(n);nk=product(mesh);nq=product(qmesh);code=0
    symmetry=.false.
    if(present(use_symmetry))symmetry=use_symmetry
    if(symmetry.and.(.not.present(atoms).or..not.present(species)))code=1
    if(any(n<1).or.any(mesh<1).or.any(qmesh<1).or.nv<1)code=1
    call comm_get_max(code,comm)
    if(code/=0)return
    if(any(mod(mesh,qmesh)/=0))code=1
    if(destination<1.or.destination>nk.or.qcount<1.or.qcount>=nq.or.maxrank<1)code=1
    if(size(states,1)/=ng.or.size(states,2)<nv.or.size(states,3)/=max(0,ke-ks+1))code=1
    if(size(levels,1)<nv+1.or.size(levels,2)/=nk.or.any(shape(k)/=[3,nk]))code=1
    if(.not.screen_finite(states).or..not.screen_finite(levels).or..not.screen_finite(k))code=1
    if(.not.screen_finite([cutoff,tolerance]).or.cutoff<=0d0.or.tolerance<=0d0.or.tolerance>=1d0)code=1
    call comm_get_max(code,comm)
    if(code/=0)return
    if(minval(levels(nv+1,:))-maxval(levels(nv,:))<=1d-8)return
    call screen_geometry(cell,rec,vol,code)
    call comm_get_max(code,comm)
    if(code/=0)return
    call system_clock(started,rate)
    allocate(local(ng,nv,nk),occupied(ng,nv,nk),sites(3,nq),selected(qcount),lookup(nk),order(3,nk), &
      frac(3,nk),radius(nq),counts(qcount),starts(qcount),sources(qcount),phase(ng),stat=astat)
    code=abs(astat);call comm_get_max(code,comm)
    if(code/=0)return
    local=0d0
    if(ke>=ks)local(:,:,ks:ke)=states(:,:nv,:)
    call comm_summation(local,occupied,size(occupied),comm)
    deallocate(local)
    call screen_sites(qmesh,sites)
    do iq=1,nq
      radius(iq)=sum(matmul(rec,dble(sites(:,iq))/qmesh)**2)
      if(all(sites(:,iq)==0))radius(iq)=huge(1d0)
    enddo
    do i=1,qcount
      selected(i)=minloc(radius,dim=1);radius(selected(i))=huge(1d0)
    enddo
    frac=matmul(transpose(cell),k)/(2d0*acos(-1d0));lookup=0
    do i=1,nk
      r=(frac(:,i)-frac(:,1))*mesh;order(:,i)=nint(r)
      if(maxval(abs(r-order(:,i)))>1d-8)code=1
      idx=modulo(order(:,i),mesh);slot=1+idx(1)+mesh(1)*(idx(2)+mesh(2)*idx(3))
      if(lookup(slot)/=0)code=1
      lookup(slot)=i
    enddo
    call comm_get_max(code,comm)
    if(code/=0)return
    ! Only the owner of this fixed target may use its real-space Hpsi buffers.
    ! Other ranks participate in collection and final failure consensus.
    if(destination>=ks.and.destination<=ke)then
      do a=1,3
        igmax(a)=ceiling(sqrt(2d0*cutoff)*sqrt(sum(cell(:,a)**2))/(2d0*acos(-1d0))+1d0)
      enddo
      counts=0
      do i=1,qcount
        q=dble(sites(:,selected(i)))/qmesh
        do a=-igmax(1),igmax(1)
          do b=-igmax(2),igmax(2)
            do c=-igmax(3),igmax(3)
              if(sum(matmul(rec,[a,b,c]+q)**2)/2d0<=cutoff)counts(i)=counts(i)+1
            enddo
          enddo
        enddo
      enddo
      if(any(counts<1))code=1
    endif
    call comm_get_max(code,comm)
    if(code/=0)return
    if(destination>=ks.and.destination<=ke)then
      mmax=maxval(counts);nrhs=nv*sum(counts);min_dim=min(ng,nrhs)
      allocate(g(3,mmax,qcount),rhs(ng,nrhs),reference(ng,nrhs),solution(ng,nrhs), &
        energy(nrhs),residual(nrhs),projected(ng,nrhs),spectrum(min_dim),rwork(5*min_dim), &
        reference_block(mmax,mmax),actual_block(mmax,mmax),gram(nv,nv),stat=astat)
      if(astat/=0)then
        code=1
      else
        call run_target()
      endif
    endif
    call comm_get_max(code,comm)
    status=code
  contains
    subroutine run_target()
      implicit none
      complex(8) :: dummy(1,1),query(1)
      real(8) :: setup_seconds,reference_seconds,svd_seconds
      integer :: offset,lwork,result_code
      gram=matmul(conjg(transpose(occupied(:,:,destination))),occupied(:,:,destination))
      do j=1,nv
        gram(j,j)=gram(j,j)-1d0
      enddo
      if(.not.screen_finite(gram).or.maxval(abs(gram))>1d-8)then
        code=1
        return
      endif
      do i=1,qcount
        q=dble(sites(:,selected(i)))/qmesh
        idx=modulo(order(:,destination)-nint(q*mesh),mesh)
        slot=1+idx(1)+mesh(1)*(idx(2)+mesh(2)*idx(3));source=lookup(slot);sources(i)=source
        wrap=nint(frac(:,source)+q-frac(:,destination));mg=0
        do a=-igmax(1),igmax(1)
          do b=-igmax(2),igmax(2)
            do c=-igmax(3),igmax(3)
              if(sum(matmul(rec,[a,b,c]+q)**2)/2d0>cutoff)cycle
              mg=mg+1;g(:,mg,i)=[a,b,c]
              if(any(2*abs(g(:,mg,i)+wrap)>=n))then
                code=1
                return
              endif
            enddo
          enddo
        enddo
        starts(i)=1
        if(i>1)starts(i)=starts(i-1)+nv*counts(i-1)
        do ig=1,counts(i)
          do at=0,ng-1
            r=dble([mod(at,n(1)),mod(at/n(1),n(2)),at/(n(1)*n(2))])/n
            phase(at+1)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(r*(g(:,ig,i)+wrap)),8))
          enddo
          do band=1,nv
            p=starts(i)+(band-1)*counts(i)+ig-1
            rhs(:,p)=phase*occupied(:,band,source);energy(p)=levels(band,source)
          enddo
        enddo
        write(*,'(a,4i6,3f12.6)')'HSEX_REDUCED_Q target/q/source/NG/qfrac ', &
          destination,selected(i),source,counts(i),q
      enddo
      allocate(qfractions(3,qcount),retained(nrhs),compact_start(qcount),representative_rhs(nrhs), &
        ref_blocks(mmax,mmax,qcount),approx_blocks(mmax,mmax,qcount),restored(mmax,mmax),stat=astat)
      if(astat/=0)then
        code=1
        return
      endif
      qfractions=dble(sites(:,selected))/spread(dble(qmesh),2,qcount)
      if(symmetry)then
        call response_q_stars(n,mesh,qmesh,cell,rec,atoms,species,frac,levels,occupied,destination, &
          qfractions,g,counts,sources,at_target,stars,code)
        if(code/=0)return
      else
        allocate(stars%representative(qcount),stars%operation(qcount),stars%column(mmax,qcount),stat=astat)
        if(astat/=0)then
          code=1
          return
        endif
        stars%nrep=qcount;stars%nsym=1;stars%column=0
        do i=1,qcount
          stars%representative(i)=i;stars%operation(i)=1
          do ig=1,counts(i)
            stars%column(ig,i)=ig
          enddo
        enddo
      endif
      nret=0;representative_rhs=.false.;compact_start=0
      do i=1,qcount
        if(stars%representative(i)/=i)cycle
        compact_start(i)=nret+1
        do j=starts(i),starts(i)+nv*counts(i)-1
          nret=nret+1;retained(nret)=j;representative_rhs(j)=.true.
        enddo
      enddo
      allocate(compact_rhs(ng,nret),stat=astat)
      if(astat/=0)then
        code=1
        return
      endif
      compact_rhs=rhs(:,retained(:nret))
      if(symmetry)write(*,'(a,5i8)')'HSEX_REDUCED_SYMMETRY_COUNTS operations/q/full_RHS/rep_q/rep_RHS ', &
        stars%nsym,qcount,nrhs,stars%nrep,nret
      projected=rhs-matmul(occupied(:,:,destination), &
        matmul(conjg(transpose(occupied(:,:,destination))),rhs))
      ! Singular values diagnose input redundancy; this work is included in
      ! setup cost, not hidden inside a claimed solver speedup.
      call zgesvd('N','N',ng,nrhs,projected,ng,spectrum,dummy,1,dummy,1,query,-1,rwork,info)
      if(info/=0)then
        code=1
        return
      endif
      lwork=max(1,nint(real(query(1),8)));allocate(work(lwork),stat=astat)
      if(astat/=0)then
        code=1
        return
      endif
      call system_clock(spectral)
      call zgesvd('N','N',ng,nrhs,projected,ng,spectrum,dummy,1,dummy,1,work,lwork,rwork,info)
      call system_clock(prepared)
      if(info/=0.or..not.screen_finite(spectrum))then
        code=1
        return
      endif
      do i=1,4
        svd_rank(i)=count(spectrum>thresholds(i)*max(spectrum(1),tiny(1d0)))
      enddo
      write(*,'(a,3i8)')'HSEX_REDUCED_SIZE grid/RHS/maxrank ',ng,nrhs,maxrank
      write(*,'(a,4i8)')'HSEX_REDUCED_RHS_RANK relative1e-2/1e-4/1e-6/1e-8 ',svd_rank
      write(*,'(a,8es15.6)')'HSEX_REDUCED_SINGULAR_LEADING ',spectrum(:min(8,min_dim))
      deallocate(projected,work,rwork)
      cg=s_sternheimer_stats();sym_cg_h=0;sym_cg_seconds=0d0
      do j=1,nrhs
        reference_h_before=cg%h_actions;reference_time_before=cg%solve_seconds
        call sternheimer_solve(at_target,occupied(:,:,destination),energy(j),rhs(:,j),reference(:,j),code,cg)
        if(code/=0)return
        if(representative_rhs(j))then
          sym_cg_h=sym_cg_h+cg%h_actions-reference_h_before
          sym_cg_seconds=sym_cg_seconds+cg%solve_seconds-reference_time_before
        endif
      enddo
      call system_clock(ref_done)
      call response_space_solve(at_target,occupied(:,:,destination),energy(retained(:nret)),compact_rhs,maxrank,tolerance, &
        solution(:,:nret),residual(:nret),reduced,result_code)
      call system_clock(solved)
      solution_error=norm2(reshape(abs(solution(:,:nret)-reference(:,retained(:nret))),[ng*nret]))/ &
        max(norm2(reshape(abs(reference(:,retained(:nret))),[ng*nret])),tiny(1d0))
      ref_blocks=0d0;approx_blocks=0d0
      do i=1,qcount
        mg=counts(i)
        do band=1,nv
          offset=starts(i)+(band-1)*mg
          call zgemm('C','N',mg,mg,ng,(1d0,0d0),rhs(:,offset:),ng,reference(:,offset:),ng, &
            (1d0,0d0),ref_blocks(:,:,i),mmax)
          if(stars%representative(i)/=i)cycle
          offset=compact_start(i)+(band-1)*mg
          call zgemm('C','N',mg,mg,ng,(1d0,0d0),compact_rhs(:,offset:),ng,solution(:,offset:),ng, &
            (1d0,0d0),approx_blocks(:,:,i),mmax)
        enddo
      enddo
      max_symmetry_error=0d0
      do i=1,qcount
        mg=counts(i);jq=stars%representative(i)
        reference_block=ref_blocks(:,:,i);actual_block=0d0;restored=0d0
        actual_block(:mg,:mg)=approx_blocks(stars%column(:mg,i),stars%column(:mg,i),jq)
        restored(:mg,:mg)=ref_blocks(stars%column(:mg,i),stars%column(:mg,i),jq)
        reference_norm=norm2(reshape(abs(reference_block(:mg,:mg)),[mg*mg]))
        error=norm2(reshape(abs(actual_block(:mg,:mg)-reference_block(:mg,:mg)),[mg*mg]))/ &
          max(reference_norm,tiny(1d0))
        symmetry_error=norm2(reshape(abs(restored(:mg,:mg)-reference_block(:mg,:mg)),[mg*mg]))/ &
          max(reference_norm,tiny(1d0))
        if(.not.screen_finite([reference_norm,error,symmetry_error,solution_error]))then
          code=8
          return
        endif
        max_symmetry_error=max(max_symmetry_error,symmetry_error)
        if(symmetry)write(*,'(a,i6,es20.10)')'HSEX_REDUCED_SYMMETRY_CHI_ERROR q/relative ',selected(i),symmetry_error
        write(*,'(a,i6,2es20.10)')'HSEX_REDUCED_TARGET_CHI_ERROR q/relative/norm ',selected(i),error,reference_norm
      enddo
      call system_clock(restored_end)
      reduced%build_seconds=reduced%build_seconds+dble(restored_end-solved)/rate
      setup_seconds=dble(prepared-started)/rate;svd_seconds=dble(prepared-spectral)/rate
      reference_seconds=dble(ref_done-prepared)/rate
      write(*,'(a,5i10)')'HSEX_REDUCED_COUNTS CG_RHS/CG_H/shared_rank/shared_H/validation_H ', &
        cg%rhs,cg%h_actions,reduced%dimension,reduced%h_actions,reduced%validation_actions
      write(*,'(a,2es20.10,i6)')'HSEX_REDUCED_ERROR solution/true_residual/status ', &
        solution_error,maxval(residual(:nret)),result_code
      write(*,'(a,6es16.7)')'HSEX_REDUCED_SECONDS setup/SVD/CG/shared_build/validation/total ', &
        setup_seconds,svd_seconds,reference_seconds,reduced%build_seconds,reduced%validation_seconds, &
        dble(restored_end-started)/rate
      if(symmetry)then
        write(*,'(a,i10,es20.10)')'HSEX_REDUCED_SYMMETRY_CG retained_H/seconds ',sym_cg_h,sym_cg_seconds
        write(*,'(a,i10,4es20.10)')'HSEX_REDUCED_SYMMETRY_VALIDATION H_probes/H_cov/projector/energy/chi ', &
          stars%h_probes,stars%max_h_covariance,stars%max_projector_error,stars%max_energy_error,max_symmetry_error
        if(max_symmetry_error>1d-6)then
          code=9
          return
        endif
      endif
      ! A rank cap is a measured negative result, never a publishable W.
      if(result_code/=0.and.result_code/=2)then
        code=result_code
      else
        code=0
        write(*,'(a,i0,a,i0)')'HSEX_REDUCED_DIAGNOSTIC_PASS owner=',rank,' converged=',merge(1,0,result_code==0)
      endif
    end subroutine
    subroutine at_target(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      call apply(destination,x,y,ierr)
    end subroutine
  end subroutine
end module
