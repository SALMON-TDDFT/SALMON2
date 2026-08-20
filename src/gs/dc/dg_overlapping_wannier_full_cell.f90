#include "config.h"
module dg_overlapping_wannier_full_cell
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::project_dg_full_cell_hamiltonian_tiles

  abstract interface
    subroutine dg_full_cell_tile_operator(tile_in,tile_out,ok)
      import real64
      complex(real64),intent(in)::tile_in(:,:)
      complex(real64),intent(out)::tile_out(:,:)
      logical,intent(out)::ok
    end subroutine
  end interface
contains
  subroutine project_dg_full_cell_hamiltonian_tiles(comm,global_spatial_count,spatial_ids,weights,&
      basis_values,row_ids,tile_width,apply_tile,matrix_rows,workspace_peak_bytes,ok,message)
    integer,intent(in)::comm,global_spatial_count,tile_width
    integer(int64),intent(in)::spatial_ids(:),row_ids(:)
    real(real64),intent(in)::weights(:)
    complex(real64),intent(in)::basis_values(:,:)
    procedure(dg_full_cell_tile_operator)::apply_tile
    complex(real64),allocatable,intent(out)::matrix_rows(:,:)
    integer(int64),intent(out)::workspace_peak_bytes
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::i,j,j0,j1,width,nstate,nlocal,rank,nproc,ierr,local_bad,global_bad,allocation_status
    integer::minimum_value,maximum_value,root,local_position
    integer,allocatable::spatial_count(:),row_count(:),row_owner(:),row_position(:)
    integer(int64)::complex_elements,integer_elements
    complex(real64),allocatable::tile_in(:,:),tile_out(:,:),local_row(:),reduced_row(:)
    logical::callback_ok
    ok=.false.;message='';workspace_peak_bytes=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    nstate=size(basis_values,1);nlocal=size(spatial_ids)
    call agree_integer(global_spatial_count,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell spatial extent';return
    endif
    call agree_integer(nstate,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell orbital extent';return
    endif
    call agree_integer(tile_width,minimum_value,maximum_value,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_value/=maximum_value)then
      message='inconsistent full-cell tile width';return
    endif
    if(global_spatial_count<=0.or.nstate<=0.or.tile_width<=0)local_bad=1
    if(size(weights)/=nlocal.or.size(basis_values,2)/=nlocal)local_bad=1
    if(any(spatial_ids<1_int64).or.any(spatial_ids>int(global_spatial_count,int64)))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(nstate,int64)))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.finite_complex(basis_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid full-cell tiled projection contract';return
    endif
    allocate(spatial_count(global_spatial_count),row_count(nstate),row_owner(nstate),row_position(nstate),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      if(allocated(spatial_count))deallocate(spatial_count)
      if(allocated(row_count))deallocate(row_count)
      if(allocated(row_owner))deallocate(row_owner)
      if(allocated(row_position))deallocate(row_position)
      message='cannot allocate full-cell ownership workspace';return
    endif
    spatial_count=0;row_count=0;row_owner=-1;row_position=0
    do i=1,nlocal;spatial_count(int(spatial_ids(i)))=spatial_count(int(spatial_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,spatial_count,global_spatial_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    do i=1,size(row_ids)
      row_count(int(row_ids(i)))=row_count(int(row_ids(i)))+1
      row_owner(int(row_ids(i)))=rank;row_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,row_count,nstate,MPI_INTEGER,MPI_SUM,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,row_owner,nstate,MPI_INTEGER,MPI_MAX,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,row_position,nstate,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=merge(0,1,all(spatial_count==1).and.all(row_count==1))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='duplicate or missing full-cell spatial/orbital owner';return
    endif
    ! Peak owned workspace: output rows, two bounded complex tiles, two row buffers,
    ! and four integer ownership arrays.  Accumulate in int64 before allocating.
    complex_elements=int(size(row_ids),int64)*int(nstate,int64)
    if(2_int64*int(tile_width,int64)>huge(complex_elements)/max(1_int64,int(nlocal,int64)))local_bad=1
    if(local_bad==0)complex_elements=complex_elements+2_int64*int(tile_width,int64)*int(nlocal,int64)+&
      2_int64*int(tile_width,int64)
    integer_elements=2_int64*int(global_spatial_count,int64)+2_int64*int(nstate,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64)local_bad=1
    if(integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='full-cell tiled projection workspace overflow';return
    endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(matrix_rows(size(row_ids),nstate),tile_in(tile_width,nlocal),tile_out(tile_width,nlocal),&
      local_row(tile_width),reduced_row(tile_width),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup_outputs();message='cannot allocate full-cell tiled projection arrays';return
    endif
    matrix_rows=(0d0,0d0)
    do j0=1,nstate,tile_width
      j1=min(nstate,j0+tile_width-1);width=j1-j0+1
      tile_in(1:width,:)=basis_values(j0:j1,:)
      call apply_tile(tile_in(1:width,:),tile_out(1:width,:),callback_ok)
      local_bad=merge(0,1,callback_ok.and.finite_complex(tile_out(1:width,:)))
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
        call cleanup_outputs();message='full-cell tile Hamiltonian callback failed';return
      endif
      do i=1,nstate
        do j=1,width
          local_row(j)=sum(weights*conjg(basis_values(i,:))*tile_out(j,:))
        enddo
        root=row_owner(i);reduced_row(1:width)=(0d0,0d0)
        call MPI_Reduce(local_row,reduced_row,width,MPI_DOUBLE_COMPLEX,MPI_SUM,root,comm,ierr)
        if(ierr/=MPI_SUCCESS)then
          call cleanup_outputs();message='full-cell projected-row reduction failed';return
        endif
        if(rank==root)then
          local_position=row_position(i)
          matrix_rows(local_position,j0:j1)=reduced_row(1:width)
        endif
      enddo
    enddo
    ok=.true.
#else
    ok=.false.;message='full-cell tiled projection requires MPI';workspace_peak_bytes=0_int64
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup_outputs()
      if(allocated(matrix_rows))deallocate(matrix_rows)
      if(allocated(tile_in))deallocate(tile_in)
      if(allocated(tile_out))deallocate(tile_out)
      if(allocated(local_row))deallocate(local_row)
      if(allocated(reduced_row))deallocate(reduced_row)
    end subroutine
#endif
  end subroutine

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function

#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine
#endif
end module dg_overlapping_wannier_full_cell
