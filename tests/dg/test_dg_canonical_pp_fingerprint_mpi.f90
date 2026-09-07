#include "config.h"
program test_dg_canonical_pp_fingerprint_mpi
  use iso_fortran_env,only:int64
  use mpi
  use structures,only:s_pp_info
  use dg_canonical_pp_fingerprint,only:canonical_pp_fingerprint,canonical_pp_valence_sum
  implicit none
  type(s_pp_info)::reference,candidate
  integer::ierr,rank,nproc
  integer(int64)::reference_fingerprint,candidate_fingerprint,minimum_fingerprint,maximum_fingerprint

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)

  call populate_pp(reference)
  reference_fingerprint=canonical_pp_fingerprint(reference)
  call require(reference_fingerprint/=0_int64,'valid pseudopotential has a zero fingerprint')

  ! zion is root-only scratch for several readers; the per-species zps array is authoritative.
  candidate=reference
  candidate%zion=1000d0+real(rank,8)
  call require(canonical_pp_fingerprint(candidate)==reference_fingerprint,&
    'root-only zion scratch changed the canonical fingerprint')

  candidate=reference
  candidate%vpp=11d0+real(rank,8)
  candidate%upp=12d0+real(rank,8)
  candidate%dvpp=13d0+real(rank,8)
  candidate%dupp=14d0+real(rank,8)
  candidate%vpp_f=15d0+real(rank,8)
  candidate%upp_f=16d0+real(rank,8)
  candidate%dupptbl_ao=21d0+real(rank,8)
  candidate_fingerprint=canonical_pp_fingerprint(candidate)
  call require(candidate_fingerprint==reference_fingerprint,&
    'root-only pseudopotential scratch changed the canonical fingerprint')

  candidate=reference
  candidate%vloctbl(candidate%nrloc(1)+1:,1)=19d0+real(rank,8)
  candidate%dvloctbl(candidate%nrloc(1)+1:,1)=20d0+real(rank,8)
  call require(canonical_pp_fingerprint(candidate)==reference_fingerprint,&
    'inactive local-table tail changed the canonical fingerprint')

  candidate=reference
  candidate%udvtbl(candidate%nrps(1)+1,0,1)=17d0+real(rank,8)
  candidate%dudvtbl(candidate%nrps(1)+1,0,1)=18d0+real(rank,8)
  call require(canonical_pp_fingerprint(candidate)==reference_fingerprint,&
    'inactive nonlocal radial tail changed the canonical fingerprint')

  candidate=reference
  candidate%udvtbl(2,sum(candidate%nproj(:,1)),1)=23d0+real(rank,8)
  candidate%dudvtbl(2,sum(candidate%nproj(:,1)),1)=24d0+real(rank,8)
  call require(canonical_pp_fingerprint(candidate)==reference_fingerprint,&
    'first undeclared nonlocal channel changed the canonical fingerprint')

  candidate=reference
  candidate%upptbl_ao(candidate%nrps_ao(1)+1,0,1)=25d0+real(rank,8)
  candidate%upptbl_ao(2,sum(candidate%nproj(:,1)),1)=26d0+real(rank,8)
  call require(canonical_pp_fingerprint(candidate)==reference_fingerprint,&
    'inactive atomic-orbital projector tail changed the canonical fingerprint')

  candidate=reference
  candidate%rho_nlcc_tbl(6,1)=27d0+real(rank,8)
  candidate%tau_nlcc_tbl(6,1)=28d0+real(rank,8)
  candidate_fingerprint=canonical_pp_fingerprint(candidate)
  call require(candidate_fingerprint==reference_fingerprint,&
    'NLCC tail after the interpolation endpoint changed the canonical fingerprint')
  call MPI_Allreduce(candidate_fingerprint,minimum_fingerprint,1,MPI_INTEGER8,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_Allreduce(candidate_fingerprint,maximum_fingerprint,1,MPI_INTEGER8,MPI_MAX,MPI_COMM_WORLD,ierr)
  call require(minimum_fingerprint==maximum_fingerprint,&
    'rank-local scratch changed the canonical fingerprint across ranks')
  call require(abs(canonical_pp_valence_sum(candidate)-10d0)<1d-14,&
    'the valence receipt did not use authoritative per-species charges')

  candidate=reference
  candidate%vloctbl(2,1)=candidate%vloctbl(2,1)+1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'an active local production-operator value was not authenticated')
  candidate=reference
  candidate%udvtbl(2,0,1)=candidate%udvtbl(2,0,1)-1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'an active nonlocal production-operator value was not authenticated')
  candidate=reference
  candidate%upptbl_ao(2,0,1)=candidate%upptbl_ao(2,0,1)+1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'an active Hybrid atomic-orbital projector value was not authenticated')
  candidate=reference
  candidate%udvtbl_so(2,0,1)=candidate%udvtbl_so(2,0,1)-1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'an active spin-orbit projector value was not authenticated')
  candidate=reference
  candidate%rho_nlcc_tbl(2,1)=candidate%rho_nlcc_tbl(2,1)+1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'an active NLCC value was not authenticated')
  candidate=reference
  candidate%rho_nlcc_tbl(4,1)=1d-7
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'the first NLCC cutoff sentinel value was not authenticated')
  candidate=reference
  candidate%rho_nlcc_tbl(5,1)=candidate%rho_nlcc_tbl(5,1)+1d-12
  call require(canonical_pp_fingerprint(candidate)/=reference_fingerprint,&
    'the NLCC interpolation endpoint was not authenticated')

  if(rank==0)write(*,'(a,i0,a)')'PASS canonical PP fingerprint on ',nproc,' ranks'
  call MPI_Finalize(ierr)
contains
  subroutine populate_pp(pp)
    type(s_pp_info),intent(out)::pp
    integer::i,j,k
    pp%zion=-999d0;pp%lmax=2;pp%lmax0=8;pp%nrmax=6;pp%nrmax0=12;pp%flag_nlcc=.true.
    allocate(pp%atom_symbol(2),pp%rmass(2),pp%mr(2),pp%lref(2),pp%nrps(2),pp%mlps(2),&
      pp%nproj(0:2,2),pp%num_orb(2),pp%zps(2),pp%nrloc(2),pp%rloc(2),pp%rps(2))
    pp%atom_symbol=['Si','O '];pp%rmass=[28d0,16d0];pp%mr=[5,4];pp%lref=[1,2]
    pp%nrps=[4,3];pp%mlps=[1,2];pp%nproj=0;pp%nproj(:,1)=[1,1,0];pp%nproj(:,2)=[2,1,1]
    pp%num_orb=0;pp%zps=[4,6];pp%nrloc=pp%nrps;pp%rloc=[0.4d0,0.3d0];pp%rps=pp%rloc
    allocate(pp%anorm(0:5,2),pp%inorm(0:5,2),pp%anorm_so(0:5,2),pp%inorm_so(0:5,2))
    allocate(pp%rad(6,2),pp%radnl(6,2),pp%vloctbl(6,2),pp%dvloctbl(6,2))
    allocate(pp%udvtbl(6,0:5,2),pp%dudvtbl(6,0:5,2))
    allocate(pp%rho_pp_tbl(6,2),pp%rho_nlcc_tbl(6,2),pp%tau_nlcc_tbl(6,2))
    pp%anorm=0d0;pp%inorm=0;pp%anorm_so=0d0;pp%inorm_so=0
    pp%anorm_so(0,1)=0.75d0;pp%inorm_so(0,1)=1
    do k=1,2
      do j=0,sum(pp%nproj(:,k))-1
        pp%anorm(j,k)=0.5d0+0.01d0*real(j+10*k,8);pp%inorm(j,k)=1
      enddo
      do i=1,6
        pp%rad(i,k)=0.1d0*real(i,8);pp%radnl(i,k)=pp%rad(i,k)
        pp%vloctbl(i,k)=-real(k,8)/real(i+1,8);pp%dvloctbl(i,k)=0.1d0*real(i+k,8)
        pp%rho_pp_tbl(i,k)=0.01d0*real(i*k,8)
        pp%rho_nlcc_tbl(i,k)=0d0;pp%tau_nlcc_tbl(i,k)=0d0
      enddo
    enddo
    pp%rho_nlcc_tbl(1:3,1)=[0.08d0,0.04d0,0.01d0]
    pp%tau_nlcc_tbl(1:3,1)=[0.02d0,0.01d0,0.002d0]
    pp%rho_nlcc_tbl(1:2,2)=[0.07d0,0.02d0]
    pp%tau_nlcc_tbl(1:2,2)=[0.01d0,0.003d0]
    pp%udvtbl=0d0;pp%dudvtbl=0d0
    allocate(pp%udvtbl_so(6,0:5,2),pp%dudvtbl_so(6,0:5,2))
    pp%udvtbl_so=0d0;pp%dudvtbl_so=0d0
    do k=1,2;do j=0,sum(pp%nproj(:,k))-1;do i=1,pp%nrps(k)
      pp%udvtbl(i,j,k)=0.001d0*real(i+10*j+100*k,8)
      pp%dudvtbl(i,j,k)=-0.002d0*real(i+10*j+100*k,8)
    enddo;enddo;enddo
    pp%udvtbl_so(1:pp%nrps(1),0,1)=[0.11d0,0.12d0,0.13d0,0.14d0]
    pp%dudvtbl_so(1:pp%nrps(1),0,1)=[-0.21d0,-0.22d0,-0.23d0,-0.24d0]
    allocate(pp%nrps_ao(2),pp%rps_ao(2),pp%upptbl_ao(6,0:5,2),pp%dupptbl_ao(6,0:5,2))
    pp%nrps_ao=[4,3];pp%rps_ao=[0.4d0,0.3d0];pp%upptbl_ao=0d0;pp%dupptbl_ao=0d0
    do k=1,2;do j=0,sum(pp%nproj(:,k))-1;do i=1,pp%nrps_ao(k)
      pp%upptbl_ao(i,j,k)=0.004d0*real(i+10*j+100*k,8)
    enddo;enddo;enddo
    allocate(pp%vpp(0:12,0:6),pp%upp(0:12,0:5),pp%dvpp(0:12,0:6),pp%dupp(0:12,0:5))
    allocate(pp%vpp_f(0:12,0:6,2),pp%upp_f(0:12,0:5,2))
    pp%vpp=0d0;pp%upp=0d0;pp%dvpp=0d0;pp%dupp=0d0;pp%vpp_f=0d0;pp%upp_f=0d0
  end subroutine populate_pp

  subroutine require(condition,message)
    logical,intent(in)::condition
    character(*),intent(in)::message
    if(.not.condition)then
      write(*,'(a)')trim(message)
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
    endif
  end subroutine require
end program test_dg_canonical_pp_fingerprint_mpi
