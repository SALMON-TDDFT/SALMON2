#include "config.h"
program test_dg_nonlocal_projector_range_mpi
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_nonlocal_projector_range,only:s_dg_nonlocal_range_receipt,&
    analyze_dg_nonlocal_projector_range
  implicit none
  integer,parameter::ngrid=8,nwann=2,natom=4,nprojector=4
  integer::ierr,rank,nproc,p,nlocal,ilma,nnz
  integer(int64),allocatable::grid_ids(:),projector_grid_ids(:)
  integer(int64),allocatable::row_ids(:)
  integer,allocatable::projector_offsets(:)
  complex(real64),allocatable::wannier(:,:),projector(:,:)
  complex(real64),allocatable::projector_values(:)
  complex(real64),allocatable::operator_rows(:,:)
  real(real64)::lattice(3,3),centers(3,nwann),atoms(3,natom),strength(nprojector)
  integer::species(natom),projector_atom(nprojector),rotation(3,3),fragment_shape(3)
  real(real64)::translation(3)
  complex(real64)::wannier_representation(nwann,nwann),projector_representation(nprojector,nprojector)
  type(s_dg_nonlocal_range_receipt)::receipt,receipt_wide
  logical::ok
  character(256)::message

  call MPI_Init(ierr);call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  nlocal=count([(mod(p-1,nproc)==rank,p=1,ngrid)])
  allocate(grid_ids(nlocal),wannier(nwann,nlocal),projector(nprojector,nlocal))
  nlocal=0
  do p=1,ngrid
    if(mod(p-1,nproc)/=rank)cycle
    nlocal=nlocal+1;grid_ids(nlocal)=p
    wannier(:,nlocal)=[cmplx(merge(1d0,0d0,p==1)+merge(0.5d0,0d0,p==3)+&
      merge(0.25d0,0d0,p==5),0d0,real64),cmplx(merge(1d0,0d0,p==5),0d0,real64)]
    projector(:,nlocal)=cmplx([merge(1d0,0d0,p==1),merge(2d0,0d0,p==3),&
      merge(1d0,0d0,p==5),merge(2d0,0d0,p==7)],0d0,real64)
  enddo
  lattice=0d0;lattice(1,1)=8d0;lattice(2,2)=2d0;lattice(3,3)=2d0
  centers=reshape([0d0,0d0,0d0,0.5d0,0d0,0d0],[3,nwann])
  atoms=reshape([0d0,0d0,0d0,2d0,0d0,0d0,4d0,0d0,0d0,6d0,0d0,0d0],[3,natom])
  species=1;projector_atom=[1,2,3,4];strength=1d0
  rotation=0;rotation(1,1)=1;rotation(2,2)=-1;rotation(3,3)=-1
  translation=0d0;fragment_shape=[4,1,1]
  wannier_representation=0d0;projector_representation=0d0
  do p=1,nwann;wannier_representation(p,p)=1d0;enddo
  do p=1,nprojector;projector_representation(p,p)=1d0;enddo
  call pack_sparse_projectors
  allocate(row_ids(count([(mod(p-1,nproc)==rank,p=1,nwann)])))
  row_ids=pack([(int(p,int64),p=1,nwann)],[(mod(p-1,nproc)==rank,p=1,nwann)])
  call analyze_dg_nonlocal_projector_range(MPI_COMM_WORLD,grid_ids,wannier,centers,lattice,&
    atoms,species,projector_atom,strength,projector_offsets,projector_grid_ids,projector_values,&
    rotation,translation,fragment_shape,1,&
    wannier_representation,projector_representation,receipt,ok,message,row_ids,operator_rows)
  call require(ok,trim(message))
  call require(all(shape(operator_rows)==[size(row_ids),nwann]),'row-owned nonlocal operator shape mismatch')
  call require(all(abs(aimag(operator_rows))<1d-14),'real synthetic nonlocal operator became complex')
  call require(receipt%local_contribution>0d0,'local projector contribution missing')
  call require(receipt%adjacent_contribution>0d0,'adjacent projector contribution missing')
  call require(receipt%remote_contribution>0d0,'broad Wannier remote contribution missing')
  call require(receipt%unmatched_channel_count==0,'exact projector partners not matched')
  call require(receipt%symmetry_pair_defect<1d-14,'exact twofold projector pairing failed')
  call analyze_dg_nonlocal_projector_range(MPI_COMM_WORLD,grid_ids,wannier,centers,lattice,&
    atoms,species,projector_atom,strength,projector_offsets,projector_grid_ids,projector_values,&
    rotation,translation,fragment_shape,2,&
    wannier_representation,projector_representation,receipt_wide,ok,message)
  call require(ok,trim(message))
  call require(maxval(abs([receipt%local_contribution-receipt_wide%local_contribution,&
    receipt%adjacent_contribution-receipt_wide%adjacent_contribution,&
    receipt%remote_contribution-receipt_wide%remote_contribution,&
    receipt%symmetry_pair_defect-receipt_wide%symmetry_pair_defect]))<1d-14,&
    'range receipts depend on tile width')
  do p=1,size(grid_ids)
    wannier(:,p)=[cmplx(merge(1d0,0d0,grid_ids(p)==1),0d0,real64),&
      cmplx(merge(1d0,0d0,grid_ids(p)==5),0d0,real64)]
    projector(:,p)=0d0
    if(grid_ids(p)==1)projector(1,p)=1d0
    if(grid_ids(p)==5)projector(3,p)=1d0
  enddo
  wannier_representation=0d0;wannier_representation(1,2)=1d0;wannier_representation(2,1)=1d0
  projector_representation=0d0
  projector_representation(1,3)=1d0;projector_representation(3,1)=1d0
  projector_representation(2,2)=1d0;projector_representation(4,4)=1d0
  call pack_sparse_projectors
  call analyze_dg_nonlocal_projector_range(MPI_COMM_WORLD,grid_ids,wannier,centers,lattice,&
    atoms,species,projector_atom,strength,projector_offsets,projector_grid_ids,projector_values,&
    rotation,translation,fragment_shape,2,&
    wannier_representation,projector_representation,receipt_wide,ok,message)
  call require(ok.and.receipt_wide%symmetry_pair_defect<1d-14,&
    'nontrivial Wannier/projector permutation covariance failed')
  projector_representation(1,1)=2d0
  call analyze_dg_nonlocal_projector_range(MPI_COMM_WORLD,grid_ids,wannier,centers,lattice,&
    atoms,species,projector_atom,strength,projector_offsets,projector_grid_ids,projector_values,&
    rotation,translation,fragment_shape,1,&
    wannier_representation,projector_representation,receipt_wide,ok,message)
  call require(.not.ok,'nonunitary projector representation was accepted')
  if(rank==0)write(*,'(a)')'PASS nonlocal projector range diagnostic'
  call MPI_Finalize(ierr)
contains
  subroutine pack_sparse_projectors
    if(allocated(projector_offsets))deallocate(projector_offsets,projector_grid_ids,projector_values)
    nnz=count(abs(projector)>0d0)
    allocate(projector_offsets(nprojector+1),projector_grid_ids(nnz),projector_values(nnz))
    nnz=0;projector_offsets(1)=1
    do ilma=1,nprojector
      do p=1,size(grid_ids)
        if(abs(projector(ilma,p))==0d0)cycle
        nnz=nnz+1;projector_grid_ids(nnz)=grid_ids(p);projector_values(nnz)=projector(ilma,p)
      enddo
      projector_offsets(ilma+1)=nnz+1
    enddo
  end subroutine

  subroutine require(condition,label)
    logical,intent(in)::condition;character(*),intent(in)::label
    if(.not.condition)then;write(0,'(a,i0,2a)')'rank ',rank,': ',trim(label);call MPI_Abort(MPI_COMM_WORLD,1,ierr);endif
  end subroutine
end program
