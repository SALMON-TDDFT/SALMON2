#include "config.h"
module dg_hybrid_broken_volume
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::assemble_dg_hybrid_broken_volume_rows,assemble_dg_hybrid_exact_nonlocal_rows,&
    collect_dg_hybrid_exact_projector_overlaps
contains
  subroutine collect_dg_hybrid_exact_projector_overlaps(comm,global_basis_count,atom_ids,ordinals,strength,&
      partial_overlap,projector_ids,projector_owner,complete_strength,complete_overlap,ok,message)
    integer,intent(in)::comm,global_basis_count,atom_ids(:),ordinals(:)
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::partial_overlap(:,:)
    integer(int64),allocatable,intent(out)::projector_ids(:)
    integer,allocatable,intent(out)::projector_owner(:)
    real(real64),allocatable,intent(out)::complete_strength(:)
    complex(real64),allocatable,intent(out)::complete_overlap(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nprojector,nproc,ierr,local_bad,global_bad,p,q
    integer,allocatable::atom_min(:),atom_max(:),ordinal_min(:),ordinal_max(:)
    real(real64),allocatable::strength_min(:),strength_max(:)
    ok=.false.;message='';local_bad=0;nprojector=size(atom_ids)
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    if(global_basis_count<1.or.nprojector<1.or.size(ordinals)/=nprojector.or.size(strength)/=nprojector.or.&
        any(shape(partial_overlap)/=[global_basis_count,nprojector]).or.any(atom_ids<1).or.any(ordinals<1).or.&
        .not.all(ieee_is_finite(strength)).or..not.finite_matrix(partial_overlap))local_bad=1
    do p=1,nprojector;do q=p+1,nprojector
      if(atom_ids(p)==atom_ids(q).and.ordinals(p)==ordinals(q))local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid exact projector overlap catalog';return;endif
    allocate(atom_min(nprojector),atom_max(nprojector),ordinal_min(nprojector),ordinal_max(nprojector),&
      strength_min(nprojector),strength_max(nprojector))
    call MPI_Allreduce(atom_ids,atom_min,nprojector,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(atom_ids,atom_max,nprojector,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(ordinals,ordinal_min,nprojector,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(ordinals,ordinal_max,nprojector,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(strength,strength_min,nprojector,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(strength,strength_max,nprojector,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(atom_min/=atom_max).or.any(ordinal_min/=ordinal_max).or.&
        any(strength_min/=strength_max))then;message='rank-disagreeing exact projector catalog';return;endif
    allocate(projector_ids(nprojector),projector_owner(nprojector),complete_strength(nprojector),&
      complete_overlap(global_basis_count,nprojector))
    projector_ids=[(int(p,int64),p=1,nprojector)]
    projector_owner=[(mod(p-1,nproc),p=1,nprojector)]
    complete_strength=strength
    call MPI_Allreduce(partial_overlap,complete_overlap,size(partial_overlap),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.finite_matrix(complete_overlap)
    if(ok)then;message='';else;message='exact projector overlap reduction failed';endif
#else
    ok=.false.;message='MPI is required for exact projector overlap collection'
#endif
  end subroutine collect_dg_hybrid_exact_projector_overlaps

  subroutine assemble_dg_hybrid_exact_nonlocal_rows(comm,global_basis_count,row_ids,projector_ids,&
      projector_owner,strength,overlap,nonlocal_rows,diagnostics,ok,message)
    integer,intent(in)::comm,global_basis_count,projector_owner(:)
    integer(int64),intent(in)::row_ids(:),projector_ids(:)
    real(real64),intent(in)::strength(:)
    complex(real64),intent(in)::overlap(:,:)
    complex(real64),allocatable,intent(out)::nonlocal_rows(:,:)
    real(real64),intent(out)::diagnostics(2)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,nprojector,nowned,local_bad,global_bad,r,i,j,p,row,nrows,offset
    integer,allocatable::counts(:),displs(:),all_rows(:),row_owner(:),owner_min(:),owner_max(:)
    integer(int64),allocatable::all_row_ids(:),id_min(:),id_max(:)
    real(real64),allocatable::strength_min(:),strength_max(:),real_min(:,:),real_max(:,:),&
      imag_min(:,:),imag_max(:,:)
    complex(real64),allocatable::partial(:,:),reduced(:,:),remote(:)
    real(real64)::local_defect,global_defect,local_scale,global_scale

    ok=.false.;message='';diagnostics=huge(1d0);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    nprojector=size(projector_ids);nowned=size(row_ids)
    if(global_basis_count<1.or.nprojector<1.or.size(projector_owner)/=nprojector.or.&
        size(strength)/=nprojector.or.any(shape(overlap)/=[global_basis_count,nprojector]).or.&
        any(projector_ids<1_int64).or.any(projector_owner<0).or.any(projector_owner>=nproc).or.&
        any(row_ids<1_int64).or.any(row_ids>int(global_basis_count,int64)).or.&
        .not.all(ieee_is_finite(strength)).or..not.finite_matrix(overlap))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid exact nonlocal contract';return;endif

    allocate(owner_min(nprojector),owner_max(nprojector),id_min(nprojector),id_max(nprojector),&
      strength_min(nprojector),strength_max(nprojector),real_min(global_basis_count,nprojector),&
      real_max(global_basis_count,nprojector),imag_min(global_basis_count,nprojector),&
      imag_max(global_basis_count,nprojector))
    call MPI_Allreduce(projector_owner,owner_min,nprojector,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(projector_owner,owner_max,nprojector,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(projector_ids,id_min,nprojector,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(projector_ids,id_max,nprojector,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(strength,strength_min,nprojector,MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(strength,strength_max,nprojector,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(real(overlap),real_min,size(overlap),MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(real(overlap),real_max,size(overlap),MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(aimag(overlap),imag_min,size(overlap),MPI_DOUBLE_PRECISION,MPI_MIN,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(aimag(overlap),imag_max,size(overlap),MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(owner_min/=owner_max).or.any(id_min/=id_max).or.&
        any(strength_min/=strength_max).or.any(real_min/=real_max).or.any(imag_min/=imag_max))then
      message='rank-disagreeing exact nonlocal projector payload';return
    endif
    do p=1,nprojector
      if(count(projector_ids==projector_ids(p))/=1)then
        message='duplicate exact nonlocal projector ID';return
      endif
    enddo

    allocate(counts(nproc),displs(nproc));call MPI_Allgather(nowned,1,MPI_INTEGER,counts,1,MPI_INTEGER,comm,ierr)
    displs(1)=0;do r=2,nproc;displs(r)=displs(r-1)+counts(r-1);enddo
    allocate(all_row_ids(sum(counts)),all_rows(sum(counts)),row_owner(global_basis_count))
    call MPI_Allgatherv(row_ids,nowned,MPI_INTEGER8,all_row_ids,counts,displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.size(all_row_ids)/=global_basis_count)then
      message='exact nonlocal basis rows are incomplete';return
    endif
    all_rows=int(all_row_ids);row_owner=-1
    do r=1,nproc;do i=1,counts(r)
      row=all_rows(displs(r)+i)
      if(row<1.or.row>global_basis_count.or.row_owner(row)/=-1)local_bad=1
      if(row>=1.and.row<=global_basis_count)row_owner(row)=r-1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0.or.any(row_owner<0))then
      message='exact nonlocal basis rows are not owned exactly once';return
    endif

    allocate(nonlocal_rows(nowned,global_basis_count));nonlocal_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=counts(r+1);offset=displs(r+1)
      allocate(partial(nrows,global_basis_count),reduced(nrows,global_basis_count));partial=(0d0,0d0)
      do p=1,nprojector
        if(projector_owner(p)/=rank)cycle
        do i=1,nrows
          row=all_rows(offset+i)
          do j=1,global_basis_count
            partial(i,j)=partial(i,j)+strength(p)*conjg(overlap(row,p))*overlap(j,p)
          enddo
        enddo
      enddo
      call MPI_Reduce(partial,reduced,nrows*global_basis_count,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='exact nonlocal row reduction failed';return;endif
      if(rank==r)nonlocal_rows=reduced
      deallocate(partial,reduced)
    enddo

    allocate(remote(global_basis_count));local_defect=0d0;local_scale=1d0
    if(nowned>0)local_scale=max(local_scale,maxval(abs(nonlocal_rows)))
    do row=1,global_basis_count
      remote=(0d0,0d0)
      if(rank==row_owner(row))then;i=findloc(row_ids,int(row,int64),dim=1);remote=nonlocal_rows(i,:);endif
      call MPI_Bcast(remote,global_basis_count,MPI_DOUBLE_COMPLEX,row_owner(row),comm,ierr)
      do i=1,nowned
        local_defect=max(local_defect,abs(nonlocal_rows(i,row)-conjg(remote(int(row_ids(i))))))
      enddo
    enddo
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr==MPI_SUCCESS)call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    diagnostics=[global_defect,global_scale]
    if(ierr/=MPI_SUCCESS.or.global_defect>1d-12*global_scale)then
      message='exact nonlocal rows are not Hermitian';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='MPI is required for exact nonlocal assembly';diagnostics=huge(1d0)
#endif
  end subroutine assemble_dg_hybrid_exact_nonlocal_rows

  subroutine assemble_dg_hybrid_broken_volume_rows(comm,global_basis_count,row_ids,basis_fragment,&
      interior_ids,interior_fragment,weights,basis_values,basis_gradients,local_potential,&
      kinetic_rows,local_rows,diagnostics,ok,message)
    integer,intent(in)::comm,global_basis_count,basis_fragment(:),interior_fragment(:)
    integer(int64),intent(in)::row_ids(:),interior_ids(:)
    real(real64),intent(in)::weights(:),local_potential(:)
    complex(real64),intent(in)::basis_values(:,:),basis_gradients(:,:,:)
    complex(real64),allocatable,intent(out)::kinetic_rows(:,:),local_rows(:,:)
    real(real64),intent(out)::diagnostics(4)
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,local_bad,global_bad,nlocal,nowned,total_points,max_point
    integer::i,j,p,r,row,nrows,offset
    integer,allocatable::row_counts(:),row_displs(:),all_rows(:),row_owner(:),point_ownership(:)
    integer(int64),allocatable::all_row_ids(:)
    complex(real64),allocatable::partial_t(:,:),partial_v(:,:),reduced_t(:,:),reduced_v(:,:),&
      remote_t(:),remote_v(:)
    real(real64)::local_defects(4),global_defects(4),t_scale,v_scale

    ok=.false.;message='';diagnostics=huge(1d0);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    nlocal=size(interior_ids);nowned=size(row_ids)
    if(global_basis_count<1.or.size(basis_fragment)/=global_basis_count.or.&
        size(interior_fragment)/=nlocal.or.size(weights)/=nlocal.or.size(local_potential)/=nlocal.or.&
        any(shape(basis_values)/=[global_basis_count,nlocal]).or.&
        any(shape(basis_gradients)/=[3,global_basis_count,nlocal]).or.&
        any(row_ids<1_int64).or.any(row_ids>int(global_basis_count,int64)).or.any(interior_ids<1_int64).or.&
        any(basis_fragment<1).or.any(interior_fragment<1).or.any(weights<=0d0))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(local_potential)).or.&
        .not.finite_matrix(basis_values).or..not.finite_tensor(basis_gradients))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid broken-volume contract';return;endif

    allocate(row_counts(nproc),row_displs(nproc))
    call MPI_Allgather(nowned,1,MPI_INTEGER,row_counts,1,MPI_INTEGER,comm,ierr)
    row_displs(1)=0
    do r=2,nproc;row_displs(r)=row_displs(r-1)+row_counts(r-1);enddo
    allocate(all_row_ids(sum(row_counts)),all_rows(sum(row_counts)),row_owner(global_basis_count))
    call MPI_Allgatherv(row_ids,nowned,MPI_INTEGER8,all_row_ids,row_counts,row_displs,MPI_INTEGER8,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.size(all_row_ids)/=global_basis_count)then
      message='broken-volume basis rows are incomplete';return
    endif
    all_rows=int(all_row_ids);row_owner=-1
    do r=1,nproc;do i=1,row_counts(r)
      row=all_rows(row_displs(r)+i)
      if(row<1.or.row>global_basis_count.or.row_owner(row)/=-1)local_bad=1
      if(row>=1.and.row<=global_basis_count)row_owner(row)=r-1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0.or.any(row_owner<0))then
      message='broken-volume basis rows are not owned exactly once';return
    endif

    call MPI_Allreduce(nlocal,total_points,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    max_point=0;if(nlocal>0)max_point=int(maxval(interior_ids))
    call MPI_Allreduce(MPI_IN_PLACE,max_point,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.max_point/=total_points)then
      message='broken-volume interior points are incomplete';return
    endif
    allocate(point_ownership(max_point));point_ownership=0
    do p=1,nlocal;point_ownership(int(interior_ids(p)))=point_ownership(int(interior_ids(p)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,point_ownership,max_point,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(point_ownership/=1))then
      message='broken-volume interior points are not owned exactly once';return
    endif

    allocate(kinetic_rows(nowned,global_basis_count),local_rows(nowned,global_basis_count))
    kinetic_rows=(0d0,0d0);local_rows=(0d0,0d0)
    do r=0,nproc-1
      nrows=row_counts(r+1);offset=row_displs(r+1)
      allocate(partial_t(nrows,global_basis_count),partial_v(nrows,global_basis_count),&
        reduced_t(nrows,global_basis_count),reduced_v(nrows,global_basis_count))
      partial_t=(0d0,0d0);partial_v=(0d0,0d0)
      do i=1,nrows
        row=all_rows(offset+i)
        do p=1,nlocal
          if(basis_fragment(row)/=interior_fragment(p))cycle
          do j=1,global_basis_count
            if(basis_fragment(j)/=interior_fragment(p))cycle
            partial_t(i,j)=partial_t(i,j)+0.5d0*weights(p)*&
              sum(conjg(basis_gradients(:,row,p))*basis_gradients(:,j,p))
            partial_v(i,j)=partial_v(i,j)+weights(p)*local_potential(p)*&
              conjg(basis_values(row,p))*basis_values(j,p)
          enddo
        enddo
      enddo
      call MPI_Reduce(partial_t,reduced_t,nrows*global_basis_count,MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Reduce(partial_v,reduced_v,nrows*global_basis_count,&
        MPI_DOUBLE_COMPLEX,MPI_SUM,r,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='broken-volume row reduction failed';return;endif
      if(rank==r)then;kinetic_rows=reduced_t;local_rows=reduced_v;endif
      deallocate(partial_t,partial_v,reduced_t,reduced_v)
    enddo

    allocate(remote_t(global_basis_count),remote_v(global_basis_count))
    local_defects=0d0;t_scale=1d0;v_scale=1d0
    if(nowned>0)then
      t_scale=max(t_scale,maxval(abs(kinetic_rows)));v_scale=max(v_scale,maxval(abs(local_rows)))
    endif
    do row=1,global_basis_count
      remote_t=(0d0,0d0);remote_v=(0d0,0d0)
      if(rank==row_owner(row))then
        i=findloc(row_ids,int(row,int64),dim=1);remote_t=kinetic_rows(i,:);remote_v=local_rows(i,:)
      endif
      call MPI_Bcast(remote_t,global_basis_count,MPI_DOUBLE_COMPLEX,row_owner(row),comm,ierr)
      if(ierr==MPI_SUCCESS)call MPI_Bcast(remote_v,global_basis_count,MPI_DOUBLE_COMPLEX,row_owner(row),comm,ierr)
      do i=1,nowned
        local_defects(1)=max(local_defects(1),abs(kinetic_rows(i,row)-conjg(remote_t(int(row_ids(i))))))
        local_defects(2)=max(local_defects(2),abs(local_rows(i,row)-conjg(remote_v(int(row_ids(i))))))
      enddo
    enddo
    local_defects(3)=t_scale;local_defects(4)=v_scale
    call MPI_Allreduce(local_defects,global_defects,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    diagnostics=global_defects
    if(ierr/=MPI_SUCCESS.or.diagnostics(1)>1d-12*diagnostics(3).or.&
        diagnostics(2)>1d-12*diagnostics(4))then
      message='broken-volume rows are not Hermitian';return
    endif
    ok=.true.;message=''
#else
    ok=.false.;message='MPI is required for broken-volume assembly';diagnostics=huge(1d0)
#endif
  end subroutine assemble_dg_hybrid_broken_volume_rows

  logical function finite_matrix(values)
    complex(real64),intent(in)::values(:,:)
    finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_matrix

  logical function finite_tensor(values)
    complex(real64),intent(in)::values(:,:,:)
    finite_tensor=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_tensor
end module dg_hybrid_broken_volume
