module dg_hybrid_lcfo_support_redistribution
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  private
  public::redistribute_dg_hybrid_lcfo_support_tile
contains
  subroutine redistribute_dg_hybrid_lcfo_support_tile(comm,bases,global_point_count,global_basis_count,&
      spatial_ids,first_column,column_count,tile_values,peak_elements,fingerprint,ok,message)
    integer,intent(in)::comm,global_point_count,global_basis_count,first_column,column_count
    type(s_dg_hybrid_fragment_basis),intent(in)::bases(:)
    integer(int64),intent(in)::spatial_ids(:)
    complex(real64),allocatable,intent(out)::tile_values(:,:)
    integer(int64),intent(out)::peak_elements,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    complex(real64),allocatable::local_values(:),global_values(:)
    integer,allocatable::basis_ownership(:),point_ownership(:)
    integer::b,i,j,column,slot,target,ierr,local_bad,global_bad
    integer(int64)::bits

    ok=.false.;message='';peak_elements=0_int64;fingerprint=0_int64;local_bad=0
    if(global_point_count<1.or.global_basis_count<1.or.first_column<1.or.column_count<1.or.&
        first_column+column_count-1>global_basis_count)local_bad=1
    if(any(spatial_ids<1_int64).or.any(spatial_ids>int(max(0,global_point_count),int64)))local_bad=1
    allocate(basis_ownership(max(1,global_basis_count)),point_ownership(max(1,global_point_count)))
    basis_ownership=0;point_ownership=0
    do i=1,size(spatial_ids);point_ownership(int(spatial_ids(i)))=point_ownership(int(spatial_ids(i)))+1;enddo
    do b=1,size(bases)
      if(.not.allocated(bases(b)%global_ids).or..not.allocated(bases(b)%buffer_point_ids).or.&
          .not.allocated(bases(b)%buffer_values))then;local_bad=1;cycle;endif
      if(size(bases(b)%buffer_values,1)/=size(bases(b)%buffer_point_ids).or.&
          size(bases(b)%buffer_values,2)/=size(bases(b)%global_ids))local_bad=1
      if(any(bases(b)%global_ids<1_int64).or.any(bases(b)%global_ids>int(global_basis_count,int64)))local_bad=1
      if(any(bases(b)%buffer_point_ids<1_int64).or.&
          any(bases(b)%buffer_point_ids>int(global_point_count,int64)))local_bad=1
      if(.not.all(ieee_is_finite(real(bases(b)%buffer_values))).or.&
          .not.all(ieee_is_finite(aimag(bases(b)%buffer_values))))local_bad=1
      do i=1,size(bases(b)%global_ids)
        basis_ownership(int(bases(b)%global_ids(i)))=basis_ownership(int(bases(b)%global_ids(i)))+1
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,basis_ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Allreduce(MPI_IN_PLACE,point_ownership,global_point_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(basis_ownership/=1).or.any(point_ownership/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid LCFO support redistribution contract';return;endif
    allocate(local_values(column_count),global_values(column_count),&
      tile_values(column_count,size(spatial_ids)),source=(0d0,0d0))
    fingerprint=int(z'6A09E667F3BCC909',int64)
    do target=1,global_point_count
      local_values=(0d0,0d0)
      do b=1,size(bases);do i=1,size(bases(b)%global_ids)
        column=int(bases(b)%global_ids(i));if(column<first_column.or.column>=first_column+column_count)cycle
        slot=column-first_column+1
        do j=1,size(bases(b)%buffer_point_ids)
          if(bases(b)%buffer_point_ids(j)==int(target,int64))local_values(slot)=bases(b)%buffer_values(j,i)
        enddo
      enddo;enddo
      call MPI_Allreduce(local_values,global_values,column_count,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;message='LCFO support point reduction failed';return;endif
      do i=1,size(spatial_ids)
        if(spatial_ids(i)==int(target,int64))tile_values(:,i)=global_values
      enddo
      do j=1,column_count
        bits=transfer(real(global_values(j)),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
        bits=transfer(aimag(global_values(j)),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
      enddo
    enddo
    peak_elements=int(size(tile_values)+size(local_values)+size(global_values),int64)
    fingerprint=ieor(fingerprint,int(first_column,int64));if(fingerprint==0_int64)fingerprint=1427_int64
    ok=.true.;message=''
  end subroutine redistribute_dg_hybrid_lcfo_support_tile
end module dg_hybrid_lcfo_support_redistribution
