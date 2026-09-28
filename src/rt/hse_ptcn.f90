! Native SALMON adapter for the verified PT-CN/ACE iteration.
module hse_ptcn
  use iso_c_binding
  use structures
  use hse_ptcn_core
  use hse_native
  use communication, only: comm_summation,comm_is_root
  use density_matrix, only: calc_density
  use hartree_sub, only: hartree
  use salmon_xc, only: exchange_correlation
  use hamiltonian, only: hpsi,update_vlocal
  implicit none
  private
  include 'fftw3.f03'
  public :: native_hse_step
contains
  subroutine native_hse_step(dt,lg,mg,system,info,stencil,xc_func,srg,srg_scalar,pp,ppg,ppn, &
      input,output,rho,rho_s,vlocal,vh,vxc,vpsl,fg,poisson,energy)
    real(8),intent(in) :: dt
    type(s_rgrid),intent(in) :: lg,mg
    type(s_dft_system),intent(inout) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_xc_functional),intent(in) :: xc_func
    type(s_sendrecv_grid),intent(inout) :: srg,srg_scalar
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(in) :: ppg
    type(s_pp_nlcc),intent(in) :: ppn
    type(s_orbital),intent(in) :: input
    type(s_orbital),intent(inout) :: output
    type(s_scalar),intent(inout) :: rho,rho_s(:),vlocal(:),vh,vxc(:)
    type(s_scalar),intent(in) :: vpsl
    type(s_reciprocal_grid),intent(inout) :: fg
    type(s_poisson),intent(inout) :: poisson
    type(s_dft_energy),intent(inout) :: energy
    type(s_orbital),save :: trial,htrial
    complex(8),allocatable,save :: u(:,:,:),x(:,:,:),fftwork(:,:),denominator(:,:),gram(:,:)
    type(c_ptr),save :: forward=c_null_ptr,backward=c_null_ptr
    real(8) :: q(3),kv(3),symbol,error,ge,ne,ne_local,t0,t1,initial_times(4),local_seconds,precondition_seconds
    integer :: ng,no,nk,n,ik,ki,ix,iy,iz,j,a,d,index,builds,apps,ierr,dims(3),clock0,clock1,rate
    complex(8),parameter :: one=(1d0,0d0),zero=(0d0,0d0)
    external :: zgemm
    call system_clock(clock0,rate)
    initial_times=hse_timings;local_seconds=0d0;precondition_seconds=0d0
    ng=product(mg%num);no=info%numo;nk=info%numk;n=mg%num(1)
    if(allocated(u))then
      if(any(shape(u)/=[ng,no,nk]))then
        call fftw_destroy_plan(forward);call fftw_destroy_plan(backward)
        call deallocate_orbital(trial);call deallocate_orbital(htrial)
        deallocate(u,x,fftwork,denominator,gram)
      endif
    endif
    if(.not.allocated(u))then
      call allocate_orbital_complex(1,mg,info,trial)
      call allocate_orbital_complex(1,mg,info,htrial)
      allocate(u(ng,no,nk),x(ng,no,nk),fftwork(ng,no),denominator(ng,nk),gram(no,no))
      dims=n
      forward=fftw_plan_many_dft(3,dims,no,fftwork,dims,1,ng,fftwork,dims,1,ng, &
        FFTW_FORWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
      backward=fftw_plan_many_dft(3,dims,no,fftwork,dims,1,ng,fftwork,dims,1,ng, &
        FFTW_BACKWARD,ior(FFTW_ESTIMATE,FFTW_UNALIGNED))
    endif
    call hse_pack(input,mg,info,u)
    if(.not.c_associated(forward).or..not.c_associated(backward))error stop 'HSE PT-CN FFT plan failed'
    do ki=1,nk
      ik=ki+info%ik_s-1;kv=system%vec_k(:,ik)+system%vec_Ac
      index=0
      do iz=0,n-1;do iy=0,n-1;do ix=0,n-1
        index=index+1;q=2*acos(-1d0)*real([ix,iy,iz],8)/n
        symbol=stencil%coef_lap0+.5d0*sum(kv*kv)
        do a=1,3;do d=1,4
          symbol=symbol-stencil%coef_lap(d,a)*cos(d*q(a))+2*kv(a)*stencil%coef_nab(d,a)*sin(d*q(a))
        enddo;enddo
        denominator(index,ki)=one+(0d0,.5d0)*dt*symbol
      enddo;enddo;enddo
    enddo
    call hse_ptcn_solve(u,dt,system%hvol,action,precondition,global_norm,x,error,builds,apps,ierr)
    hse_freeze=.false.
    if(ierr/=0)then
      if(comm_is_root(info%id_rko))write(*,*)'HSE PT-CN rejected full residual:',error
      error stop 'HSE PT-CN nonlinear convergence failed'
    endif
    ! PT-CN is not exactly orthogonality preserving; never hide drift by rescaling.
    ge=0d0
    do ki=1,nk
      call zgemm('C','N',no,no,ng,one*system%hvol,x(1,1,ki),ng,x(1,1,ki),ng,zero,gram(1,1),no)
      do j=1,no;gram(j,j)=gram(j,j)-one;enddo
      ge=max(ge,maxval(abs(gram)))
    enddo
    ne_local=2*sum(abs(x)**2)*system%hvol/system%nk
    call comm_summation(ne_local,ne,info%icomm_k)
    ierr=0
    if(ge>1d-8.or.abs(ne-2*system%no)>1d-7)ierr=1
    call comm_summation(ierr,j,info%icomm_k)
    if(j/=0)then
      write(*,*)'HSE PT-CN rejected: rank, electron number, local Gram error=',info%id_k,ne,ge
      error stop 'HSE PT-CN norm/orthogonality gate exceeded'
    endif
    call hse_unpack(x,output,mg,info)
    call system_clock(clock1)
    if(comm_is_root(info%id_rko))write(*,'(a,es14.6,a,i0,a,i0,a,f12.6)') &
      'HSE_PT_CN residual=',error,' builds=',builds,' applications=',apps,' seconds=',real(clock1-clock0,8)/rate
    if(comm_is_root(info%id_rko))write(*,'(a,6f12.6)') &
      'HSE_TIMING EXX ACE_build ACE_apply EXX_collectives local_H precondition=', &
      hse_timings-initial_times,local_seconds,precondition_seconds
  contains
    subroutine action(a,b,refresh)
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
      logical,intent(in) :: refresh
      real(8) :: start,counts(4)
      start=hse_walltime();counts=hse_timings
      call hse_unpack(a,trial,mg,info)
      call calc_density(system,rho_s,trial,info,mg)
      rho%f=rho_s(1)%f
      call hartree(lg,mg,info,system,fg,poisson,srg_scalar,stencil,rho,vh)
      hse_freeze=.not.refresh
      call exchange_correlation(system,xc_func,mg,srg_scalar,srg,rho_s,pp,ppn,info,trial,stencil,vxc,energy%E_xc)
      call update_vlocal(mg,1,vh,vpsl,vxc,vlocal)
      call hpsi(trial,htrial,info,mg,vlocal,system,stencil,srg,ppg)
      call hse_pack(htrial,mg,info,b)
      local_seconds=local_seconds+hse_walltime()-start-sum(hse_timings-counts)
    end subroutine
    subroutine precondition(a,b)
      complex(8),intent(in) :: a(:,:,:)
      complex(8),intent(out) :: b(:,:,:)
      integer :: k,j
      real(8) :: start
      start=hse_walltime()
      do k=1,nk
        fftwork=a(:,:,k)
        call fftw_execute_dft(forward,fftwork,fftwork)
        do j=1,no;fftwork(:,j)=fftwork(:,j)/denominator(:,k);enddo
        call fftw_execute_dft(backward,fftwork,fftwork)
        b(:,:,k)=fftwork/ng
      enddo
      precondition_seconds=precondition_seconds+hse_walltime()-start
    end subroutine
    real(8) function global_norm(a)
      complex(8),intent(in) :: a(:,:,:)
      real(8) :: local,total
      local=sum(abs(a)**2)
      call comm_summation(local,total,info%icomm_k)
      global_norm=sqrt(total)
    end function
  end subroutine
end module
