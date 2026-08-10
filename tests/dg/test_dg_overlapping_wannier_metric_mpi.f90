#include "config.h"
program test_dg_overlapping_wannier_metric_mpi
  use mpi
  use,intrinsic::ieee_arithmetic,only:ieee_get_halting_mode,ieee_set_halting_mode,ieee_overflow
  use dg_overlapping_wannier_metric,only:assemble_dg_overlapping_wannier_metric,&
    assemble_dg_overlapping_wannier_metric_rows,assemble_dg_eigenexa_cyclic_metric_block,&
    assemble_dg_stitched_overlap_density_rows
  implicit none
  integer::comm,rank,nproc,ierr,i,j,nlocal,owned,rejected,reference_owned
  integer::nprow,npcol,myrow,mycol,nrowlocal,ncollocal,ilocal,jlocal
  integer(8),allocatable::ids(:),row_ids(:),stitched_row_ids(:)
  real(8),allocatable::weights(:)
  complex(8),allocatable::values(:,:),metric(:,:),vectors(:,:)
  complex(8),allocatable::base_values(:,:),reference_metric(:,:),metric_rows(:,:)
  complex(8),allocatable::gamma_values(:,:)
  complex(8),allocatable::stitched_values(:,:),stitched_srows(:,:),stitched_rhorows(:,:)
  complex(8)::stitched_s_reference(2,2),stitched_rho_reference(2,2)
  integer(8)::stitched_ids(4),stitched_peak_elements
  real(8)::stitched_weights(4),stitched_density(4),stitched_electron_count,&
    stitched_s_hermiticity,stitched_rho_hermiticity,stitched_minimum_pivot,stitched_pivot_condition
  real(8),allocatable::cyclic_metric(:,:),gamma_local_metric(:,:),gamma_metric(:,:)
  real(8),allocatable::spectrum(:),reference_spectrum(:)
  complex(8)::rotation(3,3)
  logical,allocatable::pairs(:,:)
  real(8)::minimum,condition,reference_minimum,reference_condition
  integer::reference_rejected
  integer(8)::cyclic_peak_elements
  logical::ok,overflow_halting
  character(256)::message
  call MPI_Init(ierr);comm=MPI_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr);call MPI_Comm_size(comm,nproc,ierr)
  nlocal=count([(mod(i-1,nproc)==rank,i=1,4)])
  allocate(ids(nlocal),weights(nlocal),values(3,nlocal),pairs(3,nlocal))
  nlocal=0
  do i=1,4
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;ids(nlocal)=i;weights(nlocal)=0.5d0+0.1d0*i
    values(:,nlocal)=[cmplx(1d0+0.1d0*i,0.05d0*i,8),&
      cmplx((-1d0)**i*0.3d0,0.02d0*i,8),cmplx(0.2d0*i,-0.04d0*i,8)]
  enddo
  pairs=.true.
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok,trim(message));call require(owned==4,'unique ownership')
  call require(maxval(abs(metric-conjg(transpose(metric))))<1d-13,'Hermitian metric')
  call require(abs(metric(1,2))>1d-8,'off-fragment block')
  call require(minimum>0d0.and.condition>=1d0.and.rejected==0,'positive rank revelation')
  call require_same_matrix(metric)
  reference_spectrum=spectrum;base_values=values;reference_metric=metric;reference_owned=owned
  nprow=int(sqrt(real(nproc,8)))
  do while(nprow>1.and.mod(nproc,nprow)/=0);nprow=nprow-1;enddo
  npcol=nproc/nprow;myrow=mod(rank,nprow)+1;mycol=rank/nprow+1
  nrowlocal=count([(mod(i-1,nprow)==myrow-1,i=1,3)])
  ncollocal=count([(mod(i-1,npcol)==mycol-1,i=1,3)])
  gamma_values=cmplx(real(base_values),0d0,8)
  allocate(gamma_local_metric(3,3),gamma_metric(3,3));gamma_local_metric=0d0
  do j=1,3;do i=1,3
    gamma_local_metric(i,j)=sum(weights*real(conjg(gamma_values(i,:))*gamma_values(j,:)))
  enddo;enddo
  call MPI_Allreduce(gamma_local_metric,gamma_metric,9,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call assemble_dg_eigenexa_cyclic_metric_block(comm,nprow,npcol,myrow,mycol,&
    nrowlocal+2,ncollocal+1,gamma_values,weights,cyclic_metric,cyclic_peak_elements,ok,message)
  call require(ok.and.all(shape(cyclic_metric)==[nrowlocal+2,ncollocal+1]),trim(message))
  call require(all(cyclic_metric(nrowlocal+1:,:)==0d0).and.&
    all(cyclic_metric(:,ncollocal+1:)==0d0),'EigenExa padding remains zero')
  minimum=0d0
  do jlocal=1,ncollocal;do ilocal=1,nrowlocal
    i=myrow+(ilocal-1)*nprow;j=mycol+(jlocal-1)*npcol
    minimum=max(minimum,abs(cyclic_metric(ilocal,jlocal)-gamma_metric(i,j)))
  enddo;enddo
  call require(minimum<1d-13,'direct cyclic metric block matches dense Gamma reference')
  call require(cyclic_peak_elements>=int((nrowlocal+2)*(ncollocal+1),8).and.&
    cyclic_peak_elements<=int((nrowlocal+2)*(ncollocal+1)+18,8),&
    'cyclic metric peak storage includes padding and one bounded row tile')
  allocate(row_ids(count([(mod(i-1,nproc)==rank,i=1,3)])))
  nlocal=0
  do i=1,3
    if(mod(i-1,nproc)/=rank)cycle
    nlocal=nlocal+1;row_ids(nlocal)=i
  enddo
  call assemble_dg_overlapping_wannier_metric_rows(comm,3,row_ids,ids,weights,values,pairs,4_8,&
    1d-10,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok,trim(message))
  minimum=0d0
  if(size(row_ids)>0)minimum=maxval(abs(metric_rows-reference_metric(int(row_ids),:)))
  call require(minimum<1d-13,'row-owned metric reference')
  call require(maxval(abs(spectrum-reference_spectrum))<1d-12,'row-owned metric spectrum')

  stitched_ids=[1_8,2_8,3_8,4_8];stitched_weights=1d0/real(nproc,8)
  stitched_density=[1d0,2d0,3d0,4d0]
  allocate(stitched_values(2,4));stitched_values(1,:)=1d0
  stitched_values(2,:)=[1d0,-1d0,2d0,-2d0]
  stitched_s_reference=(0d0,0d0);stitched_rho_reference=(0d0,0d0)
  do i=1,4
    do j=1,2
      stitched_s_reference(:,j)=stitched_s_reference(:,j)+&
        conjg(stitched_values(:,i))*stitched_values(j,i)
      stitched_rho_reference(:,j)=stitched_rho_reference(:,j)+stitched_density(i)*&
        conjg(stitched_values(:,i))*stitched_values(j,i)
    enddo
  enddo
  stitched_row_ids=pack(row_ids,row_ids<=2_8)
  call assemble_dg_stitched_overlap_density_rows(comm,2,stitched_row_ids,stitched_ids,&
    stitched_weights,stitched_values,stitched_density,1d0,4_8,10d0,1d-12,stitched_srows,&
    stitched_rhorows,stitched_electron_count,stitched_s_hermiticity,stitched_rho_hermiticity,&
    stitched_minimum_pivot,stitched_pivot_condition,stitched_peak_elements,ok,message)
  call require(ok,trim(message))
  minimum=0d0
  if(size(stitched_srows,1)>0)minimum=maxval(abs(stitched_srows-&
    stitched_s_reference(int(stitched_row_ids),:)))
  call require(minimum<1d-12,'stitched overlap matches duplicate-buffer dense reference')
  minimum=0d0
  if(size(stitched_rhorows,1)>0)minimum=maxval(abs(stitched_rhorows-&
    stitched_rho_reference(int(stitched_row_ids),:)))
  call require(minimum<1d-12,'stitched density matches duplicate-buffer dense reference')
  call require(abs(stitched_electron_count-10d0)<1d-12.and.stitched_s_hermiticity<1d-12.and.&
    stitched_rho_hermiticity<1d-12.and.stitched_minimum_pivot>0d0.and.&
    stitched_pivot_condition>=1d0,'stitched charge, Hermiticity, and positive-rank receipts')
  call require(stitched_peak_elements>0_8.and.stitched_peak_elements<=int(128+4*size(stitched_row_ids),8),&
    'stitched assembly storage is row-tiled and affine-order independent')
  stitched_values(2,:)=stitched_values(1,:)
  call assemble_dg_stitched_overlap_density_rows(comm,2,stitched_row_ids,stitched_ids,&
    stitched_weights,stitched_values,stitched_density,1d0,4_8,10d0,1d-12,stitched_srows,&
    stitched_rhorows,stitched_electron_count,stitched_s_hermiticity,stitched_rho_hermiticity,&
    stitched_minimum_pivot,stitched_pivot_condition,stitched_peak_elements,ok,message)
  call require(.not.ok,'stitched overlap rejects rank loss')
  stitched_values(2,:)=[1d0,-1d0,2d0,-2d0]
  if(rank==0)stitched_weights(1)=stitched_weights(1)+0.25d0
  call assemble_dg_stitched_overlap_density_rows(comm,2,stitched_row_ids,stitched_ids,&
    stitched_weights,stitched_values,stitched_density,1d0,4_8,10d0,1d-12,stitched_srows,&
    stitched_rhorows,stitched_electron_count,stitched_s_hermiticity,stitched_rho_hermiticity,&
    stitched_minimum_pivot,stitched_pivot_condition,stitched_peak_elements,ok,message)
  call require(.not.ok.and.index(message,'coverage')>0,'stitched assembly rejects nonunit point coverage')
  stitched_weights=1d0/real(nproc,8)

  values=base_values;values(2,:)=-values(2,:)
  call check_invariant('sign invariance')
  values=base_values([2,1,3],:)
  call check_invariant('permutation invariance')
  rotation=(0d0,0d0);rotation(1,1)=sqrt(0.5d0);rotation(1,2)=sqrt(0.5d0)
  rotation(2,1)=-sqrt(0.5d0);rotation(2,2)=sqrt(0.5d0);rotation(3,3)=1d0
  values=matmul(rotation,base_values)
  call check_invariant('unitary candidate-window invariance')

  values=base_values;values(3,:)=values(1,:)
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok.and.rejected==1,'positive-metric null-space rank revelation')
  reference_rejected=rejected;reference_minimum=minimum;reference_condition=condition
  reference_spectrum=spectrum
  call assemble_dg_overlapping_wannier_metric_rows(comm,3,row_ids,ids,weights,values,pairs,4_8,&
    1d-10,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok.and.rejected==reference_rejected,'row-owned null-space rank revelation')
  call require(abs(minimum-reference_minimum)<1d-12,'row-owned retained minimum')
  call require(abs(condition-reference_condition)<1d-10*reference_condition,'row-owned condition')
  call require(maxval(abs(spectrum-reference_spectrum))<1d-12,'row-owned retained spectrum')

  values=base_values
  values(3,:)=values(1,:)+1d-4*values(2,:)+1d-5*cmplx(real(ids,8)**2,0d0,8)
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-12,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok,'ill-conditioned legacy metric')
  reference_rejected=rejected;reference_minimum=minimum;reference_condition=condition
  reference_spectrum=spectrum
  call assemble_dg_overlapping_wannier_metric_rows(comm,3,row_ids,ids,weights,values,pairs,4_8,&
    1d-12,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(ok.and.rejected==reference_rejected,'row-owned ill-conditioned rank')
  call require(abs(minimum-reference_minimum)<1d-12*max(1d0,reference_minimum),&
    'row-owned ill-conditioned minimum')
  call require(abs(condition-reference_condition)<1d-4*reference_condition,&
    'row-owned ill-conditioned condition')
  call require(maxval(abs(spectrum-reference_spectrum))<1d-10,'row-owned ill-conditioned spectrum')

  if(nproc>1)then
    call ieee_get_halting_mode(ieee_overflow,overflow_halting)
    call ieee_set_halting_mode(ieee_overflow,.false.)
    values=(1d0,0d0);weights=huge(1d0)/3d0
    call assemble_dg_overlapping_wannier_metric_rows(comm,3,row_ids,ids,weights,values,pairs,4_8,&
      1d-10,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(.not.ok,'nonfinite assembled row-owned metric rejection')
    call ieee_set_halting_mode(ieee_overflow,overflow_halting)
  endif

  values=base_values
  if(rank==0.and.size(ids)>0)pairs(1,1)=.false.
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(.not.ok,'missing owner-pair rejection')
  pairs=.true.

  if(nproc>1)then
    call assemble_dg_overlapping_wannier_metric_rows(comm,merge(4,3,rank==0),row_ids,ids,weights,&
      values,pairs,4_8,1d-10,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(.not.ok,'rank-inconsistent row-owned metric nwann rejection')
    call assemble_dg_overlapping_wannier_metric_rows(comm,3,row_ids,ids,weights,values,pairs,&
      merge(5_8,4_8,rank==0),1d-10,metric_rows,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(.not.ok,'rank-inconsistent row-owned metric core-count rejection')
    call assemble_dg_overlapping_wannier_metric(comm,merge(4,3,rank==0),ids,weights,values,pairs,&
      4_8,1d-10,metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(.not.ok,'rank-inconsistent metric contract rejection')
  endif

  if(rank==0.and.size(ids)>0)then
    ids=[ids,ids(1)];weights=[weights,weights(1)]
    values=reshape([values,values(:,1)],[3,size(ids)])
    pairs=reshape([pairs,pairs(:,1)],[3,size(ids)])
  endif
  call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
    metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
  call require(.not.ok,'duplicate core quadrature rejection')
  if(rank==0)then
    write(*,'(a,i0,a,i0)')'METRIC ranks=',nproc,' ownership=',reference_owned
    write(*,'(*(es24.16,1x))')[(real(reference_metric(i,1)),aimag(reference_metric(i,1)),&
      i=1,size(reference_metric,1))]
    write(*,'(a,i0,a)')'PASS overlapping-Wannier metric on ',nproc,' ranks'
  endif
  call MPI_Finalize(ierr)
contains
  subroutine check_invariant(label)
    character(*),intent(in)::label
    call assemble_dg_overlapping_wannier_metric(comm,3,ids,weights,values,pairs,4_8,1d-10,&
      metric,vectors,spectrum,minimum,condition,rejected,owned,ok,message)
    call require(ok,label)
    call require(size(spectrum)==size(reference_spectrum),label)
    call require(maxval(abs(spectrum-reference_spectrum))<1d-12,label)
  end subroutine
  subroutine require(c,label)
    logical,intent(in)::c;character(*),intent(in)::label
    integer::l,g
    l=merge(0,1,c);call MPI_Allreduce(l,g,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(g/=0)error stop label
  end subroutine
  subroutine require_same_matrix(a)
    complex(8),intent(in)::a(:,:)
    complex(8)::reference(size(a,1),size(a,2))
    reference=a
    call MPI_Bcast(reference,size(a),MPI_DOUBLE_COMPLEX,0,comm,ierr)
    call require(maxval(abs(reference-a))<1d-14,'deterministic matrix agreement')
  end subroutine
end program
