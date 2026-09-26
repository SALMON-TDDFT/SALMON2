! Fixed LCFO subspace adapter for native real-space RT routines.
module lcfo_rt_basis
  use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_orbital
  use communication, only: comm_summation
  implicit none
  logical,save :: lcfo_rt_active=.false.
  complex(8),allocatable :: lcfo_basis(:,:)
  integer,allocatable :: lcfo_counts(:),lcfo_offsets(:),lcfo_origins(:,:)
  integer :: lcfo_grid(3),lcfo_core(3),lcfo_buffer(3),lcfo_rank,lcfo_comm,lcfo_orb_rank,lcfo_orb_comm
  real(8) :: lcfo_h(3),lcfo_dv
contains
  logical function lcfo_rt_requested()
    character(16) :: value
    integer :: status
    call get_environment_variable('SALMON_LCFO_RT',value,status=status)
    lcfo_rt_requested=status==0.and.trim(value)=='1'
  end function

  subroutine lcfo_rt_configure(basis,jxyz,meta,counts,system,mg,info)
    use salmon_global, only: yn_restart,write_rt_wfn_k,checkpoint_interval,time_shutdown
    complex(8),intent(in) :: basis(:,:)
    integer,intent(in) :: jxyz(:,:),meta(20),counts(:)
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    integer :: f,a,j,n,origins(3,size(counts)),bad,total_bad
    complex(8),allocatable :: overlap(:,:)
    real(8) :: offdiag(3,3)
    if(yn_restart=='y'.or.write_rt_wfn_k=='y'.or.checkpoint_interval>0.or.time_shutdown>0d0) &
      error stop 'LCFO RT: restart/checkpoint is not supported yet; disable restart, checkpoint and shutdown output'
    if(system%nspin/=1.or.system%nk/=1.or.maxval(abs(system%vec_k))>1d-12) &
      error stop 'LCFO RT: currently requires Gamma and one spin'
    if(info%isize_r/=size(counts).or.info%isize_ro/=size(counts)*info%isize_o.or.info%numm/=1) &
      error stop 'LCFO RT: one real-space rank per fragment required'
    if(any(mg%num/=meta(7:9)).or.any(mod(meta(4:6)-meta(7:9),2)/=0)) &
      error stop 'LCFO RT: incompatible core or buffer mesh'
    lcfo_rank=info%id_r;lcfo_comm=info%icomm_r
    lcfo_orb_rank=info%id_o;lcfo_orb_comm=info%icomm_o
    if(info%id_ro/=lcfo_rank+info%isize_r*lcfo_orb_rank)error stop 'LCFO RT: rank ordering mismatch'
    lcfo_grid=meta(1:3);lcfo_core=meta(7:9);lcfo_buffer=(meta(4:6)-lcfo_core)/2
    lcfo_h=system%hgs;lcfo_dv=system%hvol
    offdiag=system%primitive_a
    do a=1,3
      offdiag(a,a)=0d0
    enddo
    if(maxval(abs(offdiag))>1d-12)error stop 'LCFO RT: orthogonal cell required'

    if(any(lcfo_buffer<0))error stop 'LCFO RT: negative buffer'
    bad=0
    do a=1,3
      do j=1,lcfo_core(a)
        if(jxyz(j,a)/=mg%is(a)+j-1)bad=1
      enddo
    enddo
    call comm_summation(bad,total_bad,lcfo_comm)
    if(total_bad/=0)error stop 'LCFO RT: rank grid must coincide with its fragment core'
    lcfo_basis=basis;lcfo_counts=counts
    allocate(lcfo_offsets(size(counts)+1),lcfo_origins(3,size(counts)))
    lcfo_offsets(1)=0
    do f=1,size(counts)
      lcfo_offsets(f+1)=lcfo_offsets(f)+counts(f)
    enddo
    origins=0;origins(:,lcfo_rank+1)=mg%is-1
    call comm_summation(origins,lcfo_origins,size(origins),lcfo_comm)
    n=size(basis,2)
    overlap=matmul(conjg(transpose(basis)),basis)*lcfo_dv
    do j=1,n
      overlap(j,j)=overlap(j,j)-1d0
    enddo
    bad=0
    if(maxval(abs(overlap))>1d-10)bad=1
    call comm_summation(bad,total_bad,lcfo_comm)
    if(total_bad/=0)error stop 'LCFO RT: core basis is not orthonormal'
    lcfo_rt_active=.true.
    if(lcfo_rank==0.and.lcfo_orb_rank==0)write(*,*) &
      'Native LCFO RT active: fragments, basis, orbital groups =',size(counts),sum(counts),info%isize_o
  end subroutine

  subroutine lcfo_collect_coefficients(local_grid,global)
    complex(8),intent(in) :: local_grid(:,:)
    complex(8),intent(out) :: global(:,:)
    complex(8),allocatable :: work(:,:)
    allocate(work(size(global,1),size(global,2)));work=0d0
    work(lcfo_offsets(lcfo_rank+1)+1:lcfo_offsets(lcfo_rank+2),:)= &
      matmul(conjg(transpose(lcfo_basis)),local_grid)*lcfo_dv
    call comm_summation(work,global,size(work),lcfo_comm)
  end subroutine

  subroutine lcfo_project_array(local_grid)
    complex(8),intent(inout) :: local_grid(:,:)
    local_grid=matmul(lcfo_basis,matmul(conjg(transpose(lcfo_basis)),local_grid))*lcfo_dv
  end subroutine

  subroutine lcfo_project_orbital(psi,mg,info)
    type(s_orbital),intent(inout) :: psi
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    complex(8),allocatable :: work(:,:)
    integer :: io,is(3),ie(3)
    if(.not.lcfo_rt_active)return
    is=mg%is;ie=mg%ie
    allocate(work(product(mg%num),info%numo))
    do io=info%io_s,info%io_e
      work(:,io-info%io_s+1)=reshape(psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[product(mg%num)])
    enddo
    call lcfo_project_array(work)
    do io=info%io_s,info%io_e
      psi%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)=reshape(work(:,io-info%io_s+1),mg%num)
    enddo
    psi%update_zwf_overlap=.false.
  end subroutine
end module
