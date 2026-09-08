module structures
  implicit none
  type s_parallel_info
    logical::flag_eigenexa_init=.false.
    integer::nrow_local=0,ncol_local=0,nprow=0,npcol=0,myrow=0,mycol=0
  end type
end module

module eigen_eigenexa
  use structures,only:s_parallel_info
  use eigen_libs_mod
  use mpi
  implicit none
contains
  subroutine eigen_pdsyevd_ex_distributed_blocks(info,n,local_matrix,eigenvalues,&
      local_eigenvectors,ok,message)
    type(s_parallel_info),intent(in)::info
    integer,intent(in)::n
    real(8),intent(inout)::local_matrix(:,:)
    real(8),intent(out)::eigenvalues(:),local_eigenvectors(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(8),allocatable::full_local(:,:),full_global(:,:)
    integer::lr,lc,gr,gc,error
    local_eigenvectors=0d0;eigenvalues=0d0
    call eigen_sx(n,n,local_matrix,info%nrow_local,eigenvalues,local_eigenvectors,info%nrow_local)
    allocate(full_local(n,n),full_global(n,n));full_local=0d0
    do lc=eigen_loop_start(1,info%npcol,info%mycol),eigen_loop_end(n,info%npcol,info%mycol)
      gc=eigen_translate_l2g(lc,info%npcol,info%mycol)
      do lr=eigen_loop_start(1,info%nprow,info%myrow),eigen_loop_end(n,info%nprow,info%myrow)
        gr=eigen_translate_l2g(lr,info%nprow,info%myrow)
        full_local(gr,gc)=local_eigenvectors(lr,lc)
      enddo
    enddo
    call MPI_Allreduce(full_local,full_global,n*n,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,error)
    if(maxval(abs(full_global))==0d0)then;ok=.false.;message='test wrapper failed to reconstruct EigenExa vectors';return;endif
    ok=.true.;message=''
  end subroutine
end module
