#include "config.h"
module dg_hybrid_density
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  abstract interface
    subroutine dg_hybrid_density_basis_provider(first_column,column_count,tile_values,ok)
      import real64
      integer,intent(in)::first_column,column_count
      complex(real64),intent(out)::tile_values(:,:)
      logical,intent(out)::ok
    end subroutine dg_hybrid_density_basis_provider
  end interface
  public::reconstruct_dg_hybrid_density,reconstruct_dg_hybrid_occupied_state
contains
  subroutine reconstruct_dg_hybrid_occupied_state(comm,global_basis_count,row_ids,metric_rows,basis_values,&
      weights,coefficients,occupations,density,gamma_rows,projector_rows,s_coefficients,electron_count,ok,message)
    integer,intent(in)::comm,global_basis_count
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::metric_rows(:,:),basis_values(:,:),coefficients(:,:)
    real(real64),intent(in)::weights(:),occupations(:)
    real(real64),allocatable,intent(out)::density(:)
    complex(real64),allocatable,intent(out)::gamma_rows(:,:),projector_rows(:,:),s_coefficients(:,:)
    real(real64),intent(out)::electron_count
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    complex(real64),allocatable::global_coefficients(:,:),global_s_coefficients(:,:),local_coefficients(:,:),&
      local_s_coefficients(:,:),weighted_coefficients(:,:),spatial_states(:,:)
    integer,allocatable::ownership(:)
    integer::i,p,nowned,nocc,npoint,ierr,local_bad,global_bad
    real(real64)::local_electron_count
    ok=.false.;message='';electron_count=0d0;nowned=size(row_ids);nocc=size(occupations);npoint=size(weights)
    local_bad=0
    if(global_basis_count<1.or.nocc<1.or.any(shape(metric_rows)/=[nowned,global_basis_count]).or.&
        any(shape(basis_values)/=[global_basis_count,npoint]).or.any(shape(coefficients)/=[nowned,nocc]).or.&
        any(row_ids<1_int64).or.any(row_ids>int(global_basis_count,int64)).or.any(weights<=0d0).or.&
        any(occupations<0d0).or..not.all(ieee_is_finite(weights)).or..not.all(ieee_is_finite(occupations)).or.&
        .not.finite_complex_state(metric_rows).or..not.finite_complex_state(basis_values).or.&
        .not.finite_complex_state(coefficients))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid concrete occupied-state reconstruction';return;endif
    allocate(ownership(global_basis_count));ownership=0
    do i=1,nowned;ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))then;message='occupied-state rows are not owned exactly once';return;endif
    allocate(local_coefficients(global_basis_count,nocc),global_coefficients(global_basis_count,nocc),&
      local_s_coefficients(global_basis_count,nocc),global_s_coefficients(global_basis_count,nocc),&
      s_coefficients(nowned,nocc),weighted_coefficients(nowned,nocc),spatial_states(nocc,npoint),&
      density(npoint),gamma_rows(nowned,global_basis_count),projector_rows(nowned,global_basis_count))
    local_coefficients=(0d0,0d0)
    do i=1,nowned;local_coefficients(int(row_ids(i)),:)=coefficients(i,:);enddo
    call MPI_Allreduce(local_coefficients,global_coefficients,size(local_coefficients),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='occupied coefficient row assembly failed';return;endif
    s_coefficients=matmul(metric_rows,global_coefficients);local_s_coefficients=(0d0,0d0)
    do i=1,nowned;local_s_coefficients(int(row_ids(i)),:)=s_coefficients(i,:);enddo
    call MPI_Allreduce(local_s_coefficients,global_s_coefficients,size(local_s_coefficients),&
      MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='metric occupied coefficient row assembly failed';return;endif
    weighted_coefficients=coefficients
    do i=1,nocc;weighted_coefficients(:,i)=occupations(i)*weighted_coefficients(:,i);enddo
    gamma_rows=matmul(weighted_coefficients,conjg(transpose(global_coefficients)))
    projector_rows=matmul(coefficients,conjg(transpose(global_s_coefficients)))
    spatial_states=matmul(transpose(global_coefficients),basis_values);density=0d0
    do i=1,nocc;density=density+occupations(i)*abs(spatial_states(i,:))**2;enddo
    local_electron_count=sum(weights*density)
    call MPI_Allreduce(local_electron_count,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    ok=ierr==MPI_SUCCESS.and.all(ieee_is_finite(density)).and.ieee_is_finite(electron_count).and.&
      finite_complex_state(gamma_rows).and.finite_complex_state(projector_rows).and.&
      finite_complex_state(s_coefficients)
    if(ok)then;message='';else;message='nonfinite concrete occupied-state reconstruction';endif
#else
    ok=.false.;message='MPI is required for concrete occupied-state reconstruction';electron_count=0d0
#endif
  end subroutine reconstruct_dg_hybrid_occupied_state

  subroutine reconstruct_dg_hybrid_density(comm,global_spatial_count,spatial_ids,weights,global_basis_count,row_ids,&
      coefficients,occupations,materialize_basis,basis_tile_width,occupied_tile_width,basis_fingerprint,tolerance,&
      density,electron_count,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_spatial_count,global_basis_count,basis_tile_width,occupied_tile_width
    integer(int64),intent(in)::spatial_ids(:),row_ids(:),basis_fingerprint
    real(real64),intent(in)::weights(:),occupations(:),tolerance
    complex(real64),intent(in)::coefficients(:,:)
    procedure(dg_hybrid_density_basis_provider)::materialize_basis
    real(real64),allocatable,intent(out)::density(:)
    real(real64),intent(out)::electron_count
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,ierr,nlocal,nowned,noccupied,local_bad,global_bad,minimum_integer,maximum_integer
    integer::i,j,row,first_basis,basis_count,first_occupied,occupied_count,allocation_status
    integer,allocatable::spatial_ownership(:),row_ownership(:),owner(:),position(:)
    integer(int64)::bits,minimum_bits,maximum_bits,complex_elements,real_elements,integer_elements,local_hash,global_hash,&
      entry_hash,quantized
    complex(real64),allocatable::basis_tile(:,:),coefficient_tile(:,:),spatial_states(:,:),coefficient_stream(:)
    real(real64)::local_maximum,global_maximum,safe_magnitude,local_electrons,quantization_scale,quantization_limit
    logical::callback_ok
    ok=.false.;message='';electron_count=0d0;workspace_peak_bytes=0_int64;fingerprint=0_int64
    nlocal=size(spatial_ids);nowned=size(row_ids);noccupied=size(occupations);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='hybrid density communicator failed';return;endif
    call agree_integer(global_spatial_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing density grid extent';return;endif
    call agree_integer(global_basis_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing density basis extent';return;endif
    call agree_integer(noccupied,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing density occupation count';return;endif
    call agree_integer(basis_tile_width,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing density basis tile';return;endif
    call agree_integer(occupied_tile_width,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing density occupied tile';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing density tolerance';return;endif
    call agree_int64(basis_fingerprint,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing density basis receipt';return;endif
    if(global_spatial_count<1.or.global_basis_count<1.or.noccupied<1.or.&
      basis_tile_width<1.or.occupied_tile_width<1)local_bad=1
    if(size(weights)/=nlocal.or.size(coefficients,1)/=nowned.or.size(coefficients,2)/=noccupied)local_bad=1
    if(any(spatial_ids<1_int64).or.any(spatial_ids>int(max(0,global_spatial_count),int64)))local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_basis_count),int64)))local_bad=1
    if(basis_fingerprint==0_int64.or..not.ieee_is_finite(tolerance).or.tolerance<1d-15.or.tolerance>1d-2)local_bad=1
    if(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or..not.all(ieee_is_finite(occupations)).or.&
      any(occupations<0d0).or..not.finite_complex(coefficients))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid hybrid density contract';return;endif
    do i=1,noccupied
      bits=transfer(occupations(i),bits);call agree_int64(bits,minimum_bits,maximum_bits,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing density occupations';return;endif
    enddo
    complex_elements=0_int64
    call add_product(complex_elements,int(basis_tile_width,int64),int(nlocal,int64),local_bad)
    call add_product(complex_elements,int(basis_tile_width,int64),int(occupied_tile_width,int64),local_bad)
    call add_product(complex_elements,int(occupied_tile_width,int64),int(nlocal,int64),local_bad)
    call add_count(complex_elements,int(occupied_tile_width,int64),local_bad)
    real_elements=int(nlocal,int64);integer_elements=0_int64
    call add_product(integer_elements,4_int64,int(global_spatial_count,int64),local_bad)
    call add_product(integer_elements,4_int64,int(global_basis_count,int64),local_bad)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.real_elements>huge(workspace_peak_bytes)/8_int64.or.&
      integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0)then
      workspace_peak_bytes=16_int64*complex_elements+8_int64*real_elements
      if(workspace_peak_bytes>huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
      if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+4_int64*integer_elements
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='hybrid density workspace receipt overflow';return;endif
    allocate(spatial_ownership(global_spatial_count),row_ownership(global_basis_count),owner(global_basis_count),&
      position(global_basis_count),density(nlocal),basis_tile(basis_tile_width,nlocal),&
      coefficient_tile(basis_tile_width,occupied_tile_width),spatial_states(occupied_tile_width,nlocal),&
      coefficient_stream(occupied_tile_width),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate hybrid density workspace';return;endif
    spatial_ownership=0;row_ownership=0;owner=-1;position=0
    do i=1,nlocal;spatial_ownership(int(spatial_ids(i)))=spatial_ownership(int(spatial_ids(i)))+1;enddo
    do i=1,nowned;row=int(row_ids(i));row_ownership(row)=row_ownership(row)+1;owner(row)=rank;position(row)=i;enddo
    call MPI_Allreduce(MPI_IN_PLACE,spatial_ownership,global_spatial_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='density grid ownership reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,row_ownership,global_basis_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='density basis ownership reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_basis_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='density basis owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,global_basis_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(spatial_ownership/=1).or.any(row_ownership/=1))then
      call cleanup();message='hybrid density ownership is not exactly once';return
    endif
    local_maximum=0d0;if(nowned>0)local_maximum=maxval(abs(coefficients))
    call MPI_Allreduce(local_maximum,global_maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    safe_magnitude=sqrt(sqrt(huge(1d0)/(64d0*real(global_basis_count,real64))))
    if(ierr/=MPI_SUCCESS.or.global_maximum>safe_magnitude)then;call cleanup();message='unsafe hybrid density coefficient magnitude';return;endif
    density=0d0
    do first_occupied=1,noccupied,occupied_tile_width
      occupied_count=min(occupied_tile_width,noccupied-first_occupied+1);spatial_states(1:occupied_count,:)=(0d0,0d0)
      do first_basis=1,global_basis_count,basis_tile_width
        basis_count=min(basis_tile_width,global_basis_count-first_basis+1)
        basis_tile(1:basis_count,:)=(0d0,0d0);call materialize_basis(first_basis,basis_count,basis_tile(1:basis_count,:),callback_ok)
        local_bad=merge(0,1,callback_ok.and.finite_complex(basis_tile(1:basis_count,:)))
        if(local_bad==0.and.size(basis_tile(1:basis_count,:))>0)then
          if(maxval(abs(basis_tile(1:basis_count,:)))>safe_magnitude)local_bad=1
        endif
        call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid density basis callback failed';return;endif
        do i=1,basis_count
          row=first_basis+i-1;coefficient_stream(1:occupied_count)=(0d0,0d0)
          if(rank==owner(row))coefficient_stream(1:occupied_count)=&
            coefficients(position(row),first_occupied:first_occupied+occupied_count-1)
          call MPI_Bcast(coefficient_stream,occupied_count,MPI_DOUBLE_COMPLEX,owner(row),comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='density coefficient broadcast failed';return;endif
          coefficient_tile(i,1:occupied_count)=coefficient_stream(1:occupied_count)
        enddo
        spatial_states(1:occupied_count,:)=spatial_states(1:occupied_count,:)+&
          matmul(transpose(coefficient_tile(1:basis_count,1:occupied_count)),basis_tile(1:basis_count,:))
      enddo
      do i=1,occupied_count
        density=density+occupations(first_occupied+i-1)*abs(spatial_states(i,:))**2
      enddo
    enddo
    local_bad=merge(0,1,all(ieee_is_finite(density)).and.all(density>=0d0))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='nonfinite reconstructed hybrid density';return;endif
    local_electrons=sum(weights*density);call MPI_Allreduce(local_electrons,electron_count,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or..not.ieee_is_finite(electron_count))then;call cleanup();message='hybrid electron-count reduction failed';return;endif
    quantization_scale=1000d0*tolerance;quantization_limit=0.25d0*real(huge(0_int64),real64)*quantization_scale
    local_hash=0_int64
    do i=1,nlocal
      if(density(i)>quantization_limit)then;local_bad=1;cycle;endif
      quantized=nint(density(i)/quantization_scale,int64)
      entry_hash=ieor(spatial_ids(i),ishftc(quantized,17));local_hash=ieor(local_hash,ishftc(entry_hash,mod(int(spatial_ids(i)),63)))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='hybrid density fingerprint range is unsafe';return;endif
    call MPI_Allreduce(local_hash,global_hash,1,MPI_INTEGER8,MPI_BXOR,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid density fingerprint reduction failed';return;endif
    fingerprint=ieor(global_hash,basis_fingerprint);if(fingerprint==0_int64)fingerprint=1237_int64
    ok=.true.;message='';call cleanup(.true.)
  contains
    subroutine agree_integer(value,minimum_value,maximum_value,status)
      integer,intent(in)::value;integer,intent(out)::minimum_value,maximum_value,status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine agree_integer
    subroutine agree_int64(value,minimum_value,maximum_value,status)
      integer(int64),intent(in)::value;integer(int64),intent(out)::minimum_value,maximum_value;integer,intent(out)::status
      call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status);if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
    end subroutine agree_int64
    logical function finite_complex(values)
      complex(real64),intent(in)::values(:,:)
      finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
    end function finite_complex
    subroutine add_product(total,left,right,bad)
      integer(int64),intent(inout)::total
      integer(int64),intent(in)::left,right
      integer,intent(inout)::bad
      integer(int64)::term
      if(bad/=0)return
      if(left<0_int64.or.right<0_int64)then;bad=1;return;endif
      if(left/=0_int64)then;if(right>huge(term)/left)then;bad=1;return;endif;endif
      term=left*right;call add_count(total,term,bad)
    end subroutine add_product
    subroutine add_count(total,term,bad)
      integer(int64),intent(inout)::total
      integer(int64),intent(in)::term
      integer,intent(inout)::bad
      if(bad/=0)return
      if(term<0_int64.or.total>huge(total)-term)then;bad=1;return;endif
      total=total+term
    end subroutine add_count
    subroutine cleanup(keep_density)
      logical,intent(in),optional::keep_density;logical::keep
      keep=.false.;if(present(keep_density))keep=keep_density
      if(allocated(spatial_ownership))deallocate(spatial_ownership)
      if(allocated(row_ownership))deallocate(row_ownership)
      if(allocated(owner))deallocate(owner)
      if(allocated(position))deallocate(position)
      if(.not.keep.and.allocated(density))deallocate(density)
      if(allocated(basis_tile))deallocate(basis_tile)
      if(allocated(coefficient_tile))deallocate(coefficient_tile)
      if(allocated(spatial_states))deallocate(spatial_states)
      if(allocated(coefficient_stream))deallocate(coefficient_stream)
    end subroutine cleanup
#else
    ok=.false.;message='MPI is required for hybrid density reconstruction';electron_count=0d0
    workspace_peak_bytes=0_int64;fingerprint=0_int64
#endif
  end subroutine reconstruct_dg_hybrid_density

  logical function finite_complex_state(values) result(finite)
    complex(real64),intent(in)::values(:,:)
    finite=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex_state
end module dg_hybrid_density
