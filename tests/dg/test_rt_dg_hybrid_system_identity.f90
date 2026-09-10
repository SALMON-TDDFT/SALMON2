program test_rt_dg_hybrid_system_identity
  use mpi
  use,intrinsic::iso_fortran_env,only:int64
  use structures,only:s_dft_system
  use rt_dg_hybrid_system_identity,only:fingerprint_rt_dg_hybrid_system
  use dg_portable_sha256,only:s_dg_sha256_context,dg_sha256_init,dg_sha256_update_bytes,dg_sha256_final
  implicit none
  type(s_dft_system)::system,changed,malformed
  integer(int64)::reference(4),digest(4),abc(3),pp_digest(4)
  type(s_dg_sha256_context)::sha
  integer::ierr,rank,nproc
  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr);call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  allocate(system%Rion(3,2),system%kion(2),system%vec_k(3,2),system%wtk(2),system%rocc(3,2,1))
  system%nion=2;system%nspin=1;system%ngrid=120;system%no=3;system%nk=2
  system%hgs=[0.2d0,0.25d0,0.3d0];system%primitive_a=0d0
  system%primitive_a(1,1)=4d0;system%primitive_a(2,2)=5d0;system%primitive_a(3,3)=6d0
  system%kion=[1,2];system%Rion=reshape([0d0,0d0,0d0,1d0,2d0,3d0],[3,2])
  system%vec_k=reshape([0d0,0d0,0d0,0.25d0,0d0,0d0],[3,2])
  system%wtk=[0.5d0,0.5d0];system%rocc=reshape([2d0,2d0,0d0,2d0,1d0,0d0],[3,2,1])
  pp_digest=[9919_int64,9920_int64,9921_int64,9922_int64]
  reference=value(system,pp_digest,[1],5,[0,0])
  call require(.not.all(reference==0_int64),'valid physical system has zero identity')
  changed=system;changed%Rion(1,2)=changed%Rion(1,2)+1d-12
  call require(any(value(changed,pp_digest,[1],5,[0,0])/=reference),'atom position omitted')
  changed=system;changed%kion(2)=1
  call require(any(value(changed,pp_digest,[1],5,[0,0])/=reference),'species omitted')
  call require(any(value(system,pp_digest+[1_int64,0_int64,0_int64,0_int64],[1],5,[0,0])/=reference),'pseudopotential omitted')
  changed=system;changed%primitive_a(1,1)=4.1d0
  call require(any(value(changed,pp_digest,[1],5,[0,0])/=reference),'lattice omitted')
  call require(any(fingerprint_rt_dg_hybrid_system(system,[4,5,7],.true.,9,pp_digest,[1],&
    .false.,.false.,.false.,.false.,.false.,5,[0,0])/=reference),'grid omitted')
  call require(any(value(system,pp_digest,[2],5,[0,0])/=reference),'XC identity omitted')
  changed=system;changed%no=4;deallocate(changed%rocc);allocate(changed%rocc(4,2,1))
  changed%rocc=0d0;changed%rocc(1:3,:,:)=system%rocc
  call require(all(value(changed,pp_digest,[1],5,[0,0])==reference),'trailing unoccupied solver state changed physical identity')
  changed=system;changed%vec_k(1,2)=changed%vec_k(1,2)+1d-12
  call require(any(value(changed,pp_digest,[1],5,[0,0])/=reference),'k-vector omitted')
  changed=system;changed%wtk(1)=0.4d0;changed%wtk(2)=0.6d0
  call require(any(value(changed,pp_digest,[1],5,[0,0])/=reference),'k weights omitted')
  changed=system;changed%rocc=0d0;changed%rocc(1,1,1)=1.45d0;changed%rocc(2,1,1)=1.05d0
  changed%rocc(1,2,1)=1.55d0;changed%rocc(2,2,1)=0.95d0
  call require(all(value(changed,pp_digest,[1],5,[0,0])==reference),&
    '300K fractional occupations changed the physical electron identity')
  call require(any(value(system,pp_digest,[1],4,[0,0])/=reference),'total electron count omitted')
  call require(any(value(system,pp_digest,[1],5,[3,2])/=reference),'spin-resolved electron specification omitted')
  malformed%nion=0;malformed%nspin=0;malformed%no=0;malformed%nk=0
  call require(all(value(malformed,pp_digest,[1],5,[0,0])==0_int64),'unallocated malformed system did not return invalid digest')
  abc=[97_int64,98_int64,99_int64];call dg_sha256_init(sha)
  call dg_sha256_update_bytes(sha,abc);call dg_sha256_final(sha,digest)
  call require(all(digest==[int(z'BA7816BF8F01CFEA',int64),int(z'414140DE5DAE2223',int64),&
    int(z'B00361A396177A9C',int64),int(z'B410FF61F20015AD',int64)]),'SHA-256 NIST abc vector mismatch')
  call dg_sha256_init(sha);call dg_sha256_final(sha,digest)
  call require(all(digest==[int(z'E3B0C44298FC1C14',int64),int(z'9AFBF4C8996FB924',int64),&
    int(z'27AE41E4649B934C',int64),int(z'A495991B7852B855',int64)]),'SHA-256 NIST empty vector mismatch')
  if(rank==0)print '(a,i0)','PASS authoritative Hybrid system identity mutations ranks=',nproc
  call MPI_Finalize(ierr)
contains
  function value(candidate,pp,xctype,electrons,spin_electrons) result(result_digest)
    type(s_dft_system),intent(in)::candidate
    integer(int64),intent(in)::pp(4)
    integer,intent(in)::xctype(:)
    integer,intent(in)::electrons,spin_electrons(2)
    integer(int64)::result_digest(4)
    result_digest=fingerprint_rt_dg_hybrid_system(candidate,[4,5,6],.true.,9,pp,xctype,&
      .false.,.false.,.false.,.false.,.false.,electrons,spin_electrons)
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
