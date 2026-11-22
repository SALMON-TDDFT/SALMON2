! Single-q validation coordinator. This intentionally is not a production W driver.
module hsex_frozen_q
  implicit none
  private
  public :: frozen_q_probe
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
  subroutine frozen_q_probe(n,mesh,qmesh,qp,cell,k,states,nv,ks,ke,icomm,rank,np,capacity,tol,eta,reference,apply,status, &
    summed_error,trace_error,probe_shifts)
    use communication, only: comm_bcast,comm_summation,comm_get_max,comm_is_root
    use hsex_screen_math, only: screen_geometry,screen_finite
    use hsex_frozen_response, only: frozen_sources,frozen_solve,frozen_trace,s_frozen_solution,s_frozen_stats
    use hsex_frozen_response, only: frozen_background_rate
    use hsex_frozen_response, only: frozen_density
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    integer,intent(in) :: n(3),mesh(3),qmesh(3),qp(3),nv,ks,ke,icomm,rank,np,capacity
    real(8),intent(in) :: cell(3,3),k(:,:),tol,eta
    complex(8),intent(in) :: states(:,:,:)
    logical,intent(in) :: reference
    procedure(grid_action) :: apply
    integer,intent(out) :: status
    real(8),optional,intent(out) :: summed_error(2,3),trace_error
    complex(8),optional,intent(in) :: probe_shifts(:)
    interface
      function hsex_peak_rss_mib() bind(C) result(rss)
        import c_double
        real(c_double) :: rss
      end function
    end interface
    type(s_frozen_solution) :: solution
    type(s_frozen_stats) :: stats
    integer :: nz,ng,nk,code,j,target,source,trial,shift(3),coord(3),idx,a,b,s,unit,ios,unfinished
    integer :: command(2),owner_rank,request_code,info,points(3,3),rsize(2,np),local_size(2,np)
    integer,allocatable :: owners(:),lookup(:),order(:,:)
    complex(8),allocatable :: buffer(:),answer(:),orbitals(:,:),ps(:,:),pt(:,:),vertices(:,:),u(:,:),v(:,:),c(:,:,:)
    complex(8),allocatable :: hpt(:,:)
    complex(8),allocatable :: induced(:,:),reference_density(:),reference_product(:,:)
    complex(8),allocatable :: sum_density(:,:),sum_reference(:,:),fold(:)
    complex(8),allocatable :: left_h(:,:),right_h(:,:),rot(:,:),density(:,:),weighted(:,:),matrix(:,:)
    complex(8),allocatable :: eigwork(:)
    real(8),allocatable :: le(:),re(:),rw(:),frac(:,:)
    real(8) :: reciprocal(3,3),volume,q(3),position(3),elapsed,seconds(3),worst,bg,bg_square
    real(8),allocatable :: relative(:),density_error(:)
    real(8) :: rss(2),rss_max(2)
    real(8),allocatable :: sum_error(:,:)
    real(8) :: closure_error
    character(64) :: output_file
    complex(8),allocatable :: shifts(:),values(:,:),ref_values(:,:),total(:,:),flat(:),density_trace(:,:)
    integer(8) :: started,finished,rate,phase_start,phase_end
    status=1;ng=product(n);nk=product(mesh);code=0;rss(1)=hsex_peak_rss_mib()
    nz=2
    if(present(probe_shifts))nz=size(probe_shifts)
    if(nz<2)return
    allocate(shifts(nz),values(nz,3),ref_values(nz,3),total(nz,3),flat(3*nz),density_trace(nz,3), &
      relative(nz),density_error(nz),sum_error(nz,3))
    shifts(:2)=[cmplx(0d0,eta,8),cmplx(0d0,eta/2d0,8)]
    if(present(probe_shifts))shifts=probe_shifts
    if(.not.screen_finite(shifts).or.any(aimag(shifts)<=0d0))code=1
    output_file='hsex-frozen-q.dat'
    if(all(modulo(qp,max(qmesh,1))==0))output_file='hsex-frozen-q0.dat'
    sum_error=-1d0;closure_error=-1d0
    if(present(summed_error))summed_error=-1d0
    if(present(trace_error))trace_error=-1d0
    if(any(qmesh<1).or.any(mesh<1).or.nv<1.or.capacity<1)code=1
    if(size(k,1)/=3.or.size(k,2)/=nk.or.size(states,1)/=ng.or.size(states,2)<nv)code=1
    if(size(states,3)/=max(0,ke-ks+1).or..not.screen_finite(states))code=1
    if(.not.screen_finite([tol,eta]).or..not.screen_finite(k))code=1
    if(reference.and.ng>512)code=1
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(any(mod(mesh,qmesh)/=0).or.tol<=0d0.or.eta<=0d0)return
    if(all(modulo(qp,qmesh)==0).and..not.present(probe_shifts))return
    call screen_geometry(cell,reciprocal,volume,code)
    if(code/=0)return
    allocate(owners(nk),lookup(nk),order(3,nk),frac(3,nk),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    local_size=0;local_size(:,rank+1)=[ks,ke]
    call comm_summation(local_size,rsize,size(rsize),icomm)
    owners=-1
    do j=1,np
      do a=rsize(1,j),rsize(2,j)
        if(a<1.or.a>nk)then
          code=1;cycle
        endif
        if(owners(a)/=-1)code=1
        owners(a)=j-1
      enddo
    enddo
    if(any(owners<0))code=1
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
    q=dble(modulo(qp,qmesh))/qmesh
    points(:,1)=[0,0,0];points(:,2)=[1,0,0];points(:,3)=[0,1,0]
    if(all(modulo(qp,qmesh)==0))then
      points(:,1)=[1,0,0];points(:,2)=[0,1,0];points(:,3)=[0,0,1]
    endif
    allocate(buffer(ng),answer(ng),orbitals(ng,nv),stat=code)
    call comm_get_max(code,icomm)
    if(code/=0)return
    if(comm_is_root(rank))then
      call system_clock(started,rate)
      total=0d0;unfinished=0;worst=0d0;seconds=0d0;bg_square=0d0
      write(*,'(a,3f12.6,a)')'HSEX_FROZEN_Q qfrac=',q,' coordinator diagnostic; no W publication'
      open(newunit=unit,file=trim(output_file),status='new',iostat=ios)
      status=merge(0,1,ios==0)
      if(status==0)then
        write(unit,'(a)',iostat=ios)'# source target trial shift reference_valid Re_chi Im_chi Re_ref Im_ref relerr residual' &
          //' eta left_basis right_basis shared_H_actions shared_pair_seconds density_relerr omega_Ha'
        if(ios/=0)status=1
        if(status==0)then
        allocate(ps(ng,nv),pt(ng,nv),vertices(ng,3),hpt(ng,nv),induced(ng,3*nz),reference_density(ng), &
          sum_density(ng,3*nz),sum_reference(ng,3*nz),fold(ng),stat=ios)
        if(ios/=0)status=1
        if(status==0)then
          sum_density=0d0;sum_reference=0d0
        endif
        pairs: do source=1,nk
          if(status/=0)exit pairs
          coord=modulo(order(:,source)+modulo(qp,qmesh)*(mesh/qmesh),mesh)
          target=lookup(1+coord(1)+mesh(1)*(coord(2)+mesh(2)*coord(3)))
          shift=nint(frac(:,source)+q-frac(:,target))
          do a=0,ng-1
            position=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
            fold(a+1)=exp(cmplx(0d0,-2d0*acos(-1d0)*sum(position*shift),8))
          enddo
          command=[2,source];call request();ps=orbitals
          command=[2,target];call request();pt=orbitals
          do b=1,nv
            call left_action(pt(:,b),hpt(:,b),code)
            if(code/=0)then
              status=1;exit pairs
            endif
          enddo
          call frozen_background_rate(pt,hpt,bg,code)
          if(code/=0)then
            status=1;exit pairs
          endif
          bg_square=bg_square+bg*bg/nk
          ref_values=0d0
          do trial=1,3
            if(any(2*abs(points(:,trial)+shift)>=n))then
              status=1;exit pairs
            endif
            do a=0,ng-1
              position=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
              vertices(a+1,trial)=exp(cmplx(0d0,2d0*acos(-1d0)*sum(position*(points(:,trial)+shift)),8))
            enddo
          enddo
          call frozen_sources(ps,pt,2d0*[(1d0,a=1,nv)],2d0*[(1d0,a=1,nv)],vertices,u,v,c,code)
          if(code/=0)then
            status=1;exit pairs
          endif
          call system_clock(phase_start)
          call frozen_solve(left_action,right_action,u,c(:,:,1),v,shifts,capacity,tol,solution,stats,code, &
            extra_sources=c(:,:,2:))
          call system_clock(phase_end)
          elapsed=dble(phase_end-phase_start)/rate;seconds(1)=seconds(1)+elapsed
          if(code/=0.and.code/=2)then
            status=1;exit pairs
          endif
          call frozen_density(solution,induced,code)
          if(code/=0)then
            status=1;exit pairs
          endif
          do j=1,3*nz
            sum_density(:,j)=sum_density(:,j)+fold*induced(:,j)*dble(ng)/(nk*volume)
          enddo
          do trial=1,3
            if(maxval(solution%residual(nz*(trial-1)+1:nz*trial))>tol)unfinished=unfinished+1
            call frozen_trace(solution,vertices(:,trial),flat,code)
            if(code/=0)then
              status=1;exit pairs
            endif
            values(:,trial)=flat(nz*(trial-1)+1:nz*trial)
            relative=-1d0;density_error=-1d0
            if(reference)then
              call system_clock(phase_start)
              if(trial==1)then
                call spectral_h(target,left_h,le,code)
                if(code/=0)then
                  status=1;exit pairs
                endif
                call spectral_h(source,right_h,re,code)
                if(code/=0)then
                  status=1;exit pairs
                endif
              endif
              rot=matmul(matmul(conjg(transpose(left_h)),u),matmul(c(:,:,trial),matmul(conjg(transpose(v)),right_h)))
              weighted=left_h
              do a=1,ng
                weighted(:,a)=conjg(vertices(:,trial))*weighted(:,a)
              enddo
              matrix=matmul(conjg(transpose(right_h)),weighted)
              do s=1,nz
                density=rot
                do b=1,ng
                  do a=1,ng
                    density(a,b)=density(a,b)/(shifts(s)-le(a)+re(b))
                  enddo
                enddo
                ref_values(s,trial)=sum(transpose(matrix)*density)
                relative(s)=abs(values(s,trial)-ref_values(s,trial))/max(1d-12,abs(ref_values(s,trial)))
                reference_product=matmul(left_h,density)
                reference_density=sum(reference_product*conjg(right_h),dim=2)
                j=nz*(trial-1)+s
                sum_reference(:,j)=sum_reference(:,j)+fold*reference_density*dble(ng)/(nk*volume)
                density_error(s)=sqrt(sum(abs(induced(:,nz*(trial-1)+s)-reference_density)**2))/ &
                  max(1d-14,sqrt(sum(abs(reference_density)**2)))
              enddo
              call system_clock(phase_end);seconds(2)=seconds(2)+dble(phase_end-phase_start)/rate
            endif
            worst=max(worst,maxval(solution%residual))
            do s=1,nz
              write(unit,'(5i8,7es21.11,3i9,3es21.11)',iostat=ios)source,target,trial,s,merge(1,0,reference), &
                real(values(s,trial),8),aimag(values(s,trial)),real(ref_values(s,trial),8), &
                aimag(ref_values(s,trial)),relative(s),solution%residual(nz*(trial-1)+s),aimag(shifts(s)), &
                stats%left_size,stats%right_size,stats%h_actions,elapsed,density_error(s),real(shifts(s),8)
              if(ios/=0)then
                status=1;exit pairs
              endif
            enddo
            flush(unit,iostat=ios)
            if(ios/=0)then
              status=1;exit pairs
            endif
          enddo
          total=total+values/(nk*volume)
          write(*,'(a,3i7,es16.7)')'HSEX_FROZEN_Q_PAIR source/target/unfinished/worst_residual ', &
            source,target,unfinished,worst
        enddo pairs
        if(status==0)then
          do trial=1,3
            do a=0,ng-1
              position=dble([mod(a,n(1)),mod(a/n(1),n(2)),a/(n(1)*n(2))])/n
              fold(a+1)=exp(cmplx(0d0,-2d0*acos(-1d0)*sum(position*points(:,trial)),8))
            enddo
            do s=1,nz
              j=nz*(trial-1)+s
              density_trace(s,trial)=sum(fold*sum_density(:,j))/ng
              if(reference)sum_error(s,trial)=sqrt(sum(abs(sum_density(:,j)-sum_reference(:,j))**2))/ &
                max(1d-14,sqrt(sum(abs(sum_reference(:,j))**2)))
            enddo
          enddo
          closure_error=maxval(abs(density_trace-total))/max(1d-12,maxval(abs(total)))
          if(.not.screen_finite(sum_error).or..not.screen_finite([closure_error]))status=1
        endif
        endif
        close(unit,iostat=ios)
        if(ios/=0)status=1
      endif
      command=0;call comm_bcast(command,icomm,0)
      call system_clock(finished);seconds(3)=dble(finished-started)/rate
      write(*,'(a,3es18.8)')'HSEX_FROZEN_Q_SECONDS solver/reference/total ',seconds
      write(*,'(a,2es18.8)')'HSEX_FROZEN_BACKGROUND rate_Ha/rate_over_eta_half ',sqrt(bg_square),sqrt(bg_square)/(eta/2d0)
      write(*,'(a,i8,es18.8)')'HSEX_FROZEN_Q_DIAGNOSTIC_DONE unfinished/worst_residual ',unfinished,worst
      write(*,'(a,i4)')'HSEX_FROZEN_Q_STATUS ',status
      if(status==0)then
        do trial=1,3
          write(*,'(a,i4,*(es22.12))')'HSEX_FROZEN_Q_CHI physical_retarded ',trial,total(:,trial)
          write(*,'(a,i4,*(es22.12))')'HSEX_FROZEN_Q_SUM_DENSITY_ERROR trial/eta/eta_half ',trial,sum_error(:,trial)
        enddo
        write(*,'(a,es22.12)')'HSEX_FROZEN_Q_TRACE_CLOSURE ',closure_error
        open(newunit=unit,file='summary-'//trim(output_file),status='new',iostat=ios)
        if(ios/=0)then
          status=1
        else
          write(unit,'(a)')'# trial omega_Ha eta_Ha Re_chi0 Im_chi0 worst_residual unfinished'
          do trial=1,3
            do s=1,nz
              write(unit,'(i6,5es24.15,i8)')trial,real(shifts(s),8),aimag(shifts(s)), &
                real(total(s,trial),8),aimag(total(s,trial)),worst,unfinished
            enddo
          enddo
          close(unit,iostat=ios)
          if(ios/=0)status=1
        endif
      endif
    else
      do
        call comm_bcast(command,icomm,0)
        if(command(1)==0)exit
        call service()
      enddo
    endif
    rss(2)=hsex_peak_rss_mib()
    call comm_get_max(rss,rss_max,2,icomm)
    if(comm_is_root(rank))write(*,'(a,2es18.8)')'HSEX_FROZEN_RSS_MIB before/after_rankmax ',rss_max
    call comm_bcast(status,icomm,0)
    call comm_bcast(sum_error,icomm,0)
    call comm_bcast(closure_error,icomm,0)
    if(present(summed_error))summed_error=sum_error(:2,:)
    if(present(trace_error))trace_error=closure_error
  contains
    subroutine request()
      implicit none
      call comm_bcast(command,icomm,0)
      call service()
    end subroutine
    subroutine service()
      implicit none
      owner_rank=owners(command(2));request_code=0
      if(command(1)==1)then
        call comm_bcast(buffer,icomm,0)
        if(rank==owner_rank)then
          call apply(command(2),buffer,answer,request_code)
          if(request_code==0)then
            if(.not.screen_finite(answer))request_code=1
          endif
          if(request_code/=0)answer=0d0
        endif
        call comm_bcast(request_code,icomm,owner_rank)
        call comm_bcast(answer,icomm,owner_rank)
      else
        if(rank==owner_rank)orbitals=states(:,:nv,command(2)-ks+1)
        call comm_bcast(orbitals,icomm,owner_rank)
      endif
    end subroutine
    subroutine left_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      command=[1,target];buffer=x;call request();y=answer;ierr=request_code
    end subroutine
    subroutine right_action(x,y,ierr)
      implicit none
      complex(8),intent(in) :: x(:)
      complex(8),intent(out) :: y(:)
      integer,intent(out) :: ierr
      command=[1,source];buffer=x;call request();y=answer;ierr=request_code
    end subroutine
    subroutine spectral_h(index,h,levels,ierr)
      implicit none
      integer,intent(in) :: index
      complex(8),allocatable,intent(out) :: h(:,:)
      real(8),allocatable,intent(out) :: levels(:)
      integer,intent(out) :: ierr
      integer :: column
      allocate(h(ng,ng),levels(ng),stat=ierr)
      if(ierr/=0)return
      do column=1,ng
        command=[1,index];buffer=0d0;buffer(column)=1d0;call request()
        if(request_code/=0)then
          ierr=1;return
        endif
        h(:,column)=answer
      enddo
      ierr=1
      if(maxval(abs(h-conjg(transpose(h))))>1d-9)return
      if(allocated(eigwork))deallocate(eigwork,rw)
      allocate(eigwork(2*ng),rw(3*ng),stat=ierr)
      if(ierr/=0)return
      call zheev('V','L',ng,h,ng,levels,eigwork,size(eigwork),rw,ierr)
      if(.not.screen_finite(h).or..not.screen_finite(levels))ierr=1
    end subroutine
  end subroutine
end module
