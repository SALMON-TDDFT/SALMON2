! Mesh HSE RT: occupied rotations compress exchange sources, never the state space.
! Mesh rows remain distributed over spatial MPI; orbital groups share local factors.
module hse_grid_rt
 use structures,only:s_dft_system,s_rgrid,s_parallel_info,s_orbital
 use salmon_global,only:yn_hse_realspace_rt,yn_hse_wannier,hse_omega,ae_shape1, &
   hse_rt_wf_radius,hse_rt_ace_interval,hse_rt_u_interval,hse_rt_fft_batch, &
   yn_hse_rt_fft_measure,yn_hse_rt_seed_distributed,hse_mlwf_maxiter,hse_mlwf_tolerance
 use communication,only:comm_summation,comm_bcast
 use hse_grid_exchange
 use hse_grid_wannier
 use hse_ace,only:hse_ace_state,hse_ace_average
 use lcfo_dist_dense,only:lcfo_distributed_ace_build
 implicit none
 private
 public :: grid_rt_enabled,grid_rt_refresh,grid_rt_add_action,grid_rt_stage
 type(s_hse_grid_exchange),save :: exchange
 type(s_hse_grid_wannier),save :: gauge
 type(hse_ace_state),save :: ace,initial_ace,midpoint_ace
 real(8),allocatable,save :: position(:,:)
 real(8),save :: lengths(3),dv
 integer,save :: spatial_comm,orbital_comm,orbital_rank,spatial_rank,step=0,builds=0
 logical,save :: ready=.false.,in_step=.false.,midpoint=.false.,scheduled=.true.
contains
 logical function grid_rt_enabled()
   grid_rt_enabled=yn_hse_realspace_rt=='y'
 end function
 subroutine initialize(system,mg,info)
   type(s_dft_system),intent(in) :: system
   type(s_rgrid),intent(in) :: mg
   type(s_parallel_info),intent(in) :: info
   integer :: n(3),x,y,z,g,status,j
   integer,allocatable :: indices(:)
   real(8) :: offdiag(3,3)
   if(ready)return
   if(system%nk/=1.or.info%numk/=1.or.info%numm/=1.or.system%nspin/=1) &
     error stop 'Grid HSE RT requires Gamma, unpolarized, single image'
   if(maxval(abs(system%vec_k))>1d-12.or.maxval(abs(system%rocc-2d0))>1d-12) &
     error stop 'Grid HSE RT requires Gamma and occupied spin pairs'
   offdiag=system%primitive_a
   do j=1,3
     lengths(j)=offdiag(j,j);offdiag(j,j)=0d0
   enddo
   if(maxval(abs(offdiag))>1d-12.or.any(lengths<=0d0))error stop 'Grid HSE RT requires orthogonal cell'
   n=nint(lengths/system%hgs);dv=system%hvol
   if(product(n)/=system%ngrid)error stop 'Grid HSE RT inconsistent global grid'
   spatial_comm=info%icomm_r;orbital_comm=info%icomm_o
   orbital_rank=info%id_o;spatial_rank=info%id_r
   allocate(indices(product(mg%num)),position(3,product(mg%num)))
   g=0
   do z=mg%is(3),mg%ie(3);do y=mg%is(2),mg%ie(2);do x=mg%is(1),mg%ie(1)
     g=g+1
     indices(g)=x+n(1)*((y-1)+n(2)*(z-1))
     position(:,g)=real([x-1,y-1,z-1],8)*system%hgs
   enddo;enddo;enddo
   if(orbital_rank==0)then
     call grid_exchange_init(exchange,n,system%hgs,hse_omega,indices,spatial_comm, &
       hse_rt_fft_batch,yn_hse_rt_fft_measure=='y',status)
     if(status/=0)error stop 'Grid HSE RT exchange initialization failed'
   endif
   ready=.true.
   if(spatial_rank==0.and.orbital_rank==0)then
     write(*,'(a)')'Real-space HSE RT active: mesh Taylor4; no fixed LCFO projection'
     write(*,'(a,3i9)')'Grid HSE global mesh: ',n
     write(*,'(a)')'Grid exchange reference: full periodic FFT, distributed source columns'
   endif
 end subroutine
 subroutine pack_all(psi,system,mg,info,c)
   type(s_orbital),intent(in) :: psi
   type(s_dft_system),intent(in) :: system
   type(s_rgrid),intent(in) :: mg
   type(s_parallel_info),intent(in) :: info
   complex(8),allocatable,intent(out) :: c(:,:)
   complex(8),allocatable :: part(:,:)
   integer :: io,ng
   ng=product(mg%num);allocate(c(ng,system%no),part(ng,system%no));part=0d0
   do io=info%io_s,info%io_e
     part(:,io)=reshape(psi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
       mg%is(3):mg%ie(3),1,io,1,1),[ng])
   enddo
   call comm_summation(part,c,size(c),info%icomm_o,0)
 end subroutine
 subroutine factor_action(factors,c,w,weight,comm)
   complex(8),intent(in) :: factors(:,:),c(:,:)
   complex(8),intent(out) :: w(:,:)
   real(8),intent(in) :: weight
   integer,intent(in) :: comm
   complex(8),allocatable :: overlap(:,:),total(:,:)
   overlap=matmul(conjg(transpose(factors)),c)*weight
   allocate(total(size(overlap,1),size(overlap,2)))
   call comm_summation(overlap,total,size(total),comm)
   w=-matmul(factors,total)
 end subroutine
 subroutine grid_rt_refresh(system,mg,info,psi,energy,force)
   type(s_dft_system),intent(in) :: system
   type(s_rgrid),intent(in) :: mg
   type(s_parallel_info),intent(in) :: info
   type(s_orbital),intent(in) :: psi
   real(8),intent(out) :: energy
   logical,optional,intent(in) :: force
   complex(8),allocatable :: c(:,:),source(:,:),w(:,:)
   integer :: status,dims(3),j
   real(8) :: local_energy
   logical :: rebuild
   call initialize(system,mg,info)
   call pack_all(psi,system,mg,info,c)
   rebuild=.not.allocated(ace%factors).or.scheduled
   if(present(force))rebuild=rebuild.or.force
   if(orbital_rank==0)then
     call grid_wannier_source(gauge,c,position,lengths,dv,spatial_comm,yn_hse_wannier=='y', &
       hse_rt_wf_radius,hse_mlwf_maxiter,hse_mlwf_tolerance,hse_rt_u_interval, &
       yn_hse_rt_seed_distributed=='y',source)
     allocate(w(size(c,1),size(c,2)))
     if(rebuild)then
       call grid_exchange_set_source(exchange,source,status)
       if(status/=0)error stop 'Grid HSE exchange source failed'
       call grid_exchange_apply(exchange,c,w,status)
       if(status/=0)error stop 'Grid HSE exchange build action failed'
       w=.25d0*w
       call lcfo_distributed_ace_build(ace,c,w,dv,spatial_comm,status)
       if(status/=0)error stop 'Grid HSE ACE metric failed (no clipping or fixed-basis fallback)'
       builds=builds+1
       if(spatial_rank==0)write(*,'(a,3i9)')'Grid HSE ACE build count/step/interval: ',builds,step,hse_rt_ace_interval
     else
       call factor_action(ace%factors(:,:,1),c,w,dv,spatial_comm)
     endif
     local_energy=0d0
     do j=1,system%no
       local_energy=local_energy+.5d0*system%rocc(j,1,1)*dv*real(sum(conjg(c(:,j))*w(:,j)),8)
     enddo
     call comm_summation(local_energy,energy,spatial_comm)
   endif
   call comm_bcast(energy,orbital_comm,0)
   if(rebuild)then
     if(orbital_rank==0)dims=shape(ace%factors)
     call comm_bcast(dims,orbital_comm,0)
     if(orbital_rank/=0)then
       if(allocated(ace%factors))deallocate(ace%factors)
       allocate(ace%factors(dims(1),dims(2),dims(3)))
     endif
     call comm_bcast(ace%factors,orbital_comm,0)
     call comm_bcast(ace%dv,orbital_comm,0)
   endif
 end subroutine
 subroutine grid_rt_add_action(psi,hpsi,system,mg,info)
   type(s_orbital),intent(in) :: psi
   type(s_orbital),intent(inout) :: hpsi
   type(s_dft_system),intent(in) :: system
   type(s_rgrid),intent(in) :: mg
   type(s_parallel_info),intent(in) :: info
   complex(8),allocatable :: c(:,:),w(:,:)
   integer :: io,j,ng
   if(.not.ready.or..not.allocated(ace%factors))error stop 'Grid HSE ACE uninitialized'
   ng=product(mg%num);allocate(c(ng,info%numo),w(ng,info%numo))
   do io=info%io_s,info%io_e
     c(:,io-info%io_s+1)=reshape(psi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
       mg%is(3):mg%ie(3),1,io,1,1),[ng])
   enddo
   if(in_step.and.midpoint)then
     call factor_action(midpoint_ace%factors(:,:,1),c,w,dv,info%icomm_r)
   else if(in_step)then
     call factor_action(initial_ace%factors(:,:,1),c,w,dv,info%icomm_r)
   else
     call factor_action(ace%factors(:,:,1),c,w,dv,info%icomm_r)
   endif
   do io=info%io_s,info%io_e
     j=io-info%io_s+1
     ! Add exchange only. Kinetic/local/nonlocal mesh action is never projected.
     hpsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,io,1,1)= &
       hpsi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),1,io,1,1)+reshape(w(:,j),mg%num)
   enddo
 end subroutine
 subroutine grid_rt_stage(stage,system,mg,info,psi)
   integer,intent(in) :: stage
   type(s_dft_system),optional,intent(in) :: system
   type(s_rgrid),optional,intent(in) :: mg
   type(s_parallel_info),optional,intent(in) :: info
   type(s_orbital),optional,intent(in) :: psi
   integer :: status,origin
   real(8) :: ignored_energy
   select case(stage)
   case(0)
     if(.not.present(system).or..not.present(mg).or..not.present(info).or..not.present(psi)) &
       error stop 'Grid HSE Taylor start requires state'
     step=step+1;origin=0
     if(ae_shape1=='impulse')origin=1
     scheduled=mod(step-origin,hse_rt_ace_interval)==0
     if(orbital_rank==0)call grid_wannier_stage(gauge,0)
     if(ae_shape1=='impulse'.and.step==1)call grid_rt_refresh(system,mg,info,psi,ignored_energy,force=.true.)
     initial_ace=ace;in_step=.true.;midpoint=.false.
   case(1)
     call hse_ace_average(initial_ace,ace,midpoint_ace,status)
     if(status/=0)error stop 'Grid HSE midpoint ACE failed'
     midpoint=.true.
     if(orbital_rank==0)call grid_wannier_stage(gauge,1)
   case(2)
     in_step=.false.;midpoint=.false.
     if(allocated(initial_ace%factors))deallocate(initial_ace%factors)
     if(allocated(midpoint_ace%factors))deallocate(midpoint_ace%factors)
     if(orbital_rank==0)call grid_wannier_stage(gauge,2)
   case default
     error stop 'Invalid grid HSE Taylor stage'
   end select
 end subroutine
end module
