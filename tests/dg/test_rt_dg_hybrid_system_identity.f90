program test_rt_dg_hybrid_system_identity
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  use structures,only:s_dft_system
  use rt_dg_hybrid_system_identity,only:fingerprint_rt_dg_hybrid_system
  use dg_portable_sha256,only:dg_sha256_mix_int64
  implicit none
  type(s_dft_system)::system,changed
  integer(int64)::reference,collision_left,collision_right,collision_word
  integer::ierr,rank,nproc
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  allocate(system%Rion(3,2),system%kion(2),system%vec_k(3,2),system%wtk(2),system%rocc(3,2,1))
  system%nion=2;system%nspin=1;system%ngrid=120;system%no=3;system%nk=2
  system%hgs=[0.2d0,0.25d0,0.3d0];system%primitive_a=0d0
  system%primitive_a(1,1)=4d0;system%primitive_a(2,2)=5d0;system%primitive_a(3,3)=6d0
  system%kion=[1,2];system%Rion=reshape([0d0,0d0,0d0,1d0,2d0,3d0],[3,2])
  system%vec_k=reshape([0d0,0d0,0d0,0.25d0,0d0,0d0],[3,2])
  system%wtk=[0.5d0,0.5d0];system%rocc=reshape([2d0,2d0,0d0,2d0,1d0,0d0],[3,2,1])
  reference=value(system,9919_int64,[1])
  call require(reference/=0_int64,'valid physical system has zero identity')
  changed=system;changed%Rion(1,2)=changed%Rion(1,2)+1d-12
  call require(value(changed,9919_int64,[1])/=reference,'atom position omitted')
  changed=system;changed%kion(2)=1
  call require(value(changed,9919_int64,[1])/=reference,'species omitted')
  call require(value(system,9920_int64,[1])/=reference,'pseudopotential omitted')
  changed=system;changed%primitive_a(1,1)=4.1d0
  call require(value(changed,9919_int64,[1])/=reference,'lattice omitted')
  call require(fingerprint_rt_dg_hybrid_system(system,[4,5,7],.true.,9,9919_int64,[1],&
    .false.,.false.,.false.,.false.,.false.)/=reference,'grid omitted')
  call require(value(system,9919_int64,[2])/=reference,'XC identity omitted')
  changed=system;changed%no=2;changed%rocc=changed%rocc(:2,:,:)
  call require(value(changed,9919_int64,[1])/=reference,'orbital count omitted')
  changed=system;changed%vec_k(1,2)=changed%vec_k(1,2)+1d-12
  call require(value(changed,9919_int64,[1])/=reference,'k-vector omitted')
  changed=system;changed%wtk(1)=0.4d0;changed%wtk(2)=0.6d0
  call require(value(changed,9919_int64,[1])/=reference,'k weights omitted')
  changed=system;changed%rocc(2,2,1)=changed%rocc(2,2,1)+1d-12
  call require(value(changed,9919_int64,[1])/=reference,'occupations/electron count omitted')
  collision_word=int(z'0000000000008001',int64)
  collision_left=int(z'243F6A8885A308D3',int64);collision_right=collision_left
  call dg_sha256_mix_int64(collision_left,0_int64)
  call dg_sha256_mix_int64(collision_right,collision_word)
  call require(collision_left/=collision_right,'portable SHA-256 receipt preserves rotate/XOR two-bit collision')
  if(rank==0)print '(a,i0)','PASS authoritative Hybrid system identity mutations ranks=',nproc
  call MPI_Finalize(ierr)
contains
  integer(int64) function value(candidate,pp,xctype)
    type(s_dft_system),intent(in)::candidate
    integer(int64),intent(in)::pp
    integer,intent(in)::xctype(:)
    value=fingerprint_rt_dg_hybrid_system(candidate,[4,5,6],.true.,9,pp,xctype,&
      .false.,.false.,.false.,.false.,.false.)
  end function value
  subroutine require(condition,message)
    logical,intent(in)::condition;character(*),intent(in)::message
    integer::bad,global_bad
    bad=merge(0,1,condition)
    call MPI_Allreduce(bad,global_bad,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(.not.condition)write(0,'(a)')trim(message)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_rt_dg_hybrid_system_identity
