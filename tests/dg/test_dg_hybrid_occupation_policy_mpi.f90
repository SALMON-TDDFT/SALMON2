#include "config.h"
program test_dg_hybrid_occupation_policy_mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use occupation_kernel,only:solve_spectrum_occupations
  use dg_hybrid_occupation_policy,only:s_dg_hybrid_occupation_result,derive_dg_hybrid_occupation_policy
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer::comm,rank,nproc,ierr
  real(real64),allocatable::spectrum(:,:,:),weights(:),kernel_occupations(:,:,:),eigenvalues(:)
  real(real64)::chemical_potential,electron_count
  type(s_dg_hybrid_occupation_result)::result
  integer(int64)::reference_fingerprint
  logical::ok
  character(256)::message
#ifdef USE_MPI
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
#else
  comm=0;rank=0;nproc=1
#endif

  allocate(spectrum(4,1,1),weights(1));weights=1d0
  spectrum(:,1,1)=[-1d0,-0.2d0,0.4d0,1d0]
  call solve_spectrum_occupations(spectrum,weights,4d0,0d0,.false.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok,'zero-temperature shared occupation kernel failed: '//trim(message))
  call require(all(abs(kernel_occupations(:,1,1)-[2d0,2d0,0d0,0d0])<1d-14).and.&
    abs(electron_count-4d0)<1d-10,'zero-temperature shared occupation convention changed')
  call require(abs(chemical_potential+0.1875d0)<1d-14,&
    'zero-temperature chemical-potential convergence changed')

  spectrum(:,1,1)=[-1d0,0d0,0d0,1d0]
  call solve_spectrum_occupations(spectrum,weights,6d0,0d0,.false.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok.and.all(abs(kernel_occupations(:,1,1)-[2d0,2d0,2d0,0d0])<1d-14),&
    'filled Fermi-degenerate shell changed the shared occupation convention')

  call solve_spectrum_occupations(spectrum,weights,4d0,0.2d0,.false.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok.and.abs(kernel_occupations(2,1,1)-kernel_occupations(3,1,1))<1d-14.and.&
    abs(electron_count-4d0)<1d-9,'finite-temperature Fermi degeneracy is inconsistent')

  deallocate(spectrum);allocate(spectrum(2,1,1));spectrum(:,1,1)=[-1d0,-0.5d0]
  call solve_spectrum_occupations(spectrum,weights,1d0,0.02d0,.false.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok.and.abs(electron_count-1d0)<1d-9,&
    'finite-temperature bracket near the lowest eigenvalue did not converge')

  deallocate(spectrum,weights);allocate(spectrum(1,2,1),weights(2))
  spectrum(1,:,1)=[-1d0,1d0];weights=[1d0-1d-10,1d-10]
  call solve_spectrum_occupations(spectrum,weights,2d0,0.02d0,.false.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok.and.abs(chemical_potential-0.9999999990686774d0)<1d-12.and.&
    abs(kernel_occupations(1,2,1)-0.9999999767169357d0)<1d-8,&
    'finite-temperature small-weight k-point occupation changed')

  deallocate(spectrum,weights);allocate(spectrum(1,1,2),weights(1));spectrum=-0.5d0;weights=1d0
  call solve_spectrum_occupations(spectrum,weights,1d0,0d0,.true.,kernel_occupations,&
    chemical_potential,electron_count,ok,message)
  call require(ok.and.all(abs(kernel_occupations-1d0)<1d-14).and.abs(electron_count-1d0)<1d-10,&
    'spin-orbit duplicate-component convention changed')

  allocate(eigenvalues(4));eigenvalues=[-1d0,-0.2d0,0.4d0,1d0]
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,4d0,0d0,1d-10,result,ok,message)
  call require(ok,'zero-temperature Hybrid occupation policy failed: '//trim(message))
  call require(result%valid.and.result%noccupied==2.and.abs(result%e_homo+0.2d0)<1d-14.and.&
    all(abs(result%occupations-[2d0,2d0,0d0,0d0])<1d-14),&
    'Hybrid occupied count or HOMO does not follow the final spectrum')
  reference_fingerprint=result%fingerprint

  eigenvalues=[-1d0,0d0,0d0,1d0]
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,6d0,0d0,1d-10,result,ok,message)
  call require(ok.and.result%noccupied==3.and.result%e_homo==0d0,&
    'Hybrid Fermi-degenerate HOMO selection is inconsistent')

  eigenvalues=[0d0,-1d0,1d0,2d0]
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,2d0,0d0,1d-10,result,ok,message)
  call require(.not.ok.and.index(message,'ascending')>0,&
    'Hybrid policy accepted reordered final eigenstates')

  deallocate(eigenvalues);allocate(eigenvalues(1));eigenvalues=-1d0
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,3d0,0d0,1d-10,result,ok,message)
  call require(.not.ok.and.index(message,'capacity')>0,&
    'Hybrid policy accepted insufficient spectrum capacity')

  deallocate(eigenvalues);allocate(eigenvalues(2));eigenvalues=[-1d0,1d0]
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,2d0,1d0/33d0,1d-12,result,ok,message)
  call require(ok.and.result%noccupied==1.and.result%e_homo==-1d0.and.&
    result%omitted_occupation_tail>0d0.and.result%omitted_occupation_tail<1d-12,&
    'Hybrid occupation threshold or omitted-tail receipt is incorrect')
  call derive_dg_hybrid_occupation_policy(comm,eigenvalues,2d0,1d0/33d0,1d-15,result,ok,message)
  call require(.not.ok.and.index(message,'tail')>0,&
    'Hybrid policy accepted an omitted occupation tail above tolerance')

  if(nproc>1)then
    eigenvalues=[-1d0,merge(1.1d0,1d0,rank==0)]
    call derive_dg_hybrid_occupation_policy(comm,eigenvalues,2d0,0d0,1d-10,result,ok,message)
    call require(.not.ok.and.index(message,'rank')>0,&
      'Hybrid policy accepted rank-dependent final eigenvalues')
  endif
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'HYBRID_OCCUPATION_POLICY ranks=',nproc,' fingerprint=',reference_fingerprint
    write(*,'(a,i0,a)')'PASS hybrid occupation policy on ',nproc,' ranks'
  endif
#ifdef USE_MPI
  call MPI_Finalize(ierr)
#endif
contains
  subroutine require(condition,label)
    logical,intent(in)::condition
    character(*),intent(in)::label
    integer::local_bad,global_bad
    local_bad=merge(0,1,condition)
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
#else
    global_bad=local_bad
#endif
    if(global_bad/=0)then
      if(rank==0)write(0,'(a)')trim(label)
#ifdef USE_MPI
      call MPI_Abort(comm,1,ierr)
#else
      error stop label
#endif
    endif
  end subroutine require
end program test_dg_hybrid_occupation_policy_mpi
