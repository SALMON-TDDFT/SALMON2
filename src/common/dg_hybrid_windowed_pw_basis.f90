#include "config.h"
module dg_hybrid_windowed_pw_basis
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use dg_hybrid_windowed_pw_types,only:s_dg_hybrid_basis_catalog
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_dg_hybrid_windowed_pw_basis,materialize_dg_hybrid_windowed_pw_columns
contains
  subroutine build_dg_hybrid_windowed_pw_basis(comm,global_row_count,row_ids,coordinates,raw_windows,&
      fragment_action,row_action,g_vectors,reciprocal_rotation,g_action,g_star,g_conjugate,tile_width,tolerance,&
      windows,catalog,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count,tile_width
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::coordinates(:,:),raw_windows(:,:),g_vectors(:,:),reciprocal_rotation(:,:,:),tolerance
    integer,intent(in)::fragment_action(:,:),row_action(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable,intent(out)::windows(:,:)
    type(s_dg_hybrid_basis_catalog),intent(out)::catalog
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::rank,nproc,ierr,nlocal,nfragment,nwindow_operation,nreciprocal_operation,ng,nstar,i,j,f,op,target,root,position
    integer::minimum_integer,maximum_integer,local_bad,global_bad,allocation_status,packet_index
    integer,allocatable::ownership_count(:),owner(:),owner_position(:),star_count(:)
    integer(int64)::bits,minimum_bits,maximum_bits,real_elements,complex_elements,integer_elements,quantized,&
      packet_count,packet_descriptor_bytes
    real(real64)::scale,norm_scaled,covariance_defect
    real(real64),allocatable::remote_window(:),scaled(:)
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)return
    nlocal=size(row_ids);nfragment=size(raw_windows,1);nwindow_operation=size(fragment_action,2)
    nreciprocal_operation=size(g_action,2);ng=size(g_star)
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid spatial row count';return
    endif
    call agree_integer(nfragment,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid fragment count';return
    endif
    call agree_integer(nwindow_operation,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid window operation count';return
    endif
    call agree_integer(nreciprocal_operation,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid reciprocal operation count';return
    endif
    call agree_integer(ng,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid reciprocal mode count';return
    endif
    call agree_integer(tile_width,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid PW tile width';return
    endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent hybrid basis tolerance';return
    endif
    if(global_row_count<1.or.nfragment<1.or.nwindow_operation<1.or.nreciprocal_operation<1.or.ng<1.or.tile_width<1)&
      local_bad=1
    if(size(raw_windows,2)/=nlocal.or.any(shape(coordinates)/=[3,nlocal]))local_bad=1
    if(any(shape(fragment_action)/=[nfragment,nwindow_operation]))local_bad=1
    if(any(shape(row_action)/=[global_row_count,nwindow_operation]))local_bad=1
    if(any(shape(g_vectors)/=[3,ng]).or.any(shape(reciprocal_rotation)/=[3,3,nreciprocal_operation]).or.&
      any(shape(g_action)/=[ng,nreciprocal_operation]))local_bad=1
    if(size(g_conjugate)/=ng)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance))local_bad=1
    if(.not.all(ieee_is_finite(coordinates)).or..not.all(ieee_is_finite(raw_windows)))local_bad=1
    if(.not.all(ieee_is_finite(g_vectors)).or..not.all(ieee_is_finite(reciprocal_rotation)))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid windowed PW basis shape or finite contract';return
    endif
    if(tolerance<1d-15.or.tolerance>1d-2)local_bad=1
    do op=1,nwindow_operation
      do f=1,nfragment
        call agree_integer(fragment_action(f,op),minimum_integer,maximum_integer,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
          message='inconsistent hybrid fragment action';return
        endif
      enddo
      do i=1,global_row_count
        call agree_integer(row_action(i,op),minimum_integer,maximum_integer,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
          message='inconsistent hybrid spatial action';return
        endif
      enddo
    enddo
    do op=1,nreciprocal_operation
      do i=1,ng
        call agree_integer(g_action(i,op),minimum_integer,maximum_integer,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
          message='inconsistent hybrid reciprocal action';return
        endif
      enddo
      do j=1,3;do i=1,3
        bits=transfer(reciprocal_rotation(i,j,op),bits)
        call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
          message='inconsistent hybrid reciprocal rotation';return
        endif
      enddo;enddo
    enddo
    do i=1,ng
      call agree_integer(g_star(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid G-star catalog';return
      endif
      call agree_integer(g_conjugate(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid reciprocal conjugate catalog';return
      endif
      do j=1,3
        bits=transfer(g_vectors(j,i),bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
        if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
          message='inconsistent hybrid reciprocal vector';return
        endif
      enddo
    enddo
    if(any(fragment_action<1).or.any(fragment_action>nfragment))local_bad=1
    if(any(row_action<1).or.any(row_action>global_row_count))local_bad=1
    if(any(g_action<1).or.any(g_action>ng))local_bad=1
    if(any(g_star<1).or.any(g_conjugate<1).or.any(g_conjugate>ng))local_bad=1
    if(local_bad==0)then
      do op=1,nwindow_operation
        if(.not.is_permutation(fragment_action(:,op),nfragment))local_bad=1
        if(.not.is_permutation(row_action(:,op),global_row_count))local_bad=1
      enddo
      do op=1,nreciprocal_operation
        if(.not.is_permutation(g_action(:,op),ng))local_bad=1
      enddo
      do i=1,ng
        if(g_conjugate(g_conjugate(i))/=i)local_bad=1
        if(g_star(g_conjugate(i))/=g_star(i))local_bad=1
        if(maxval(abs(g_vectors(:,g_conjugate(i))+g_vectors(:,i)))>tolerance)local_bad=1
        do op=1,nreciprocal_operation
          if(g_star(g_action(i,op))/=g_star(i))local_bad=1
          if(maxval(abs(g_vectors(:,g_action(i,op))-&
            matmul(reciprocal_rotation(:,:,op),g_vectors(:,i))))>100d0*tolerance)local_bad=1
        enddo
      enddo
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid symmetry action or incomplete G star';return
    endif
    nstar=maxval(g_star)
    packet_count=int(nfragment,int64)*int(nstar,int64)
    if(packet_count>int(huge(0),int64))local_bad=1
    if(int(nfragment,int64)*int(ng,int64)>int(huge(0),int64))local_bad=1
    if(packet_count>huge(packet_descriptor_bytes)/64_int64)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='hybrid packet extent is not representable';return
    endif
    packet_descriptor_bytes=64_int64*packet_count
    allocate(ownership_count(global_row_count),owner(global_row_count),owner_position(global_row_count),&
      star_count(nstar),remote_window(nfragment),scaled(nfragment),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid basis ownership workspace';return
    endif
    star_count=0
    do i=1,ng;star_count(g_star(i))=star_count(g_star(i))+1;enddo
    if(any(star_count==0))local_bad=1
    ownership_count=0;owner=-1;owner_position=0
    do i=1,nlocal
      ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1
      owner(int(row_ids(i)))=rank;owner_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid row ownership count failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid row owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid row position reduction failed';return;endif
    if(any(ownership_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='duplicate or missing hybrid spatial row';return
    endif
    real_elements=int(nfragment,int64)*int(nlocal,int64)+int(nfragment,int64)
    complex_elements=int(tile_width,int64)*int(nlocal,int64)
    integer_elements=3_int64*int(global_row_count,int64)+int(nstar,int64)+&
      int(nfragment,int64)*int(ng,int64)+3_int64*packet_count
    if(real_elements>huge(workspace_peak_bytes)/8_int64)local_bad=1
    if(complex_elements>huge(workspace_peak_bytes)/16_int64)local_bad=1
    if(integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0)then
      workspace_peak_bytes=8_int64*real_elements
      if(16_int64*complex_elements>huge(workspace_peak_bytes)-workspace_peak_bytes)local_bad=1
      if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+16_int64*complex_elements
      if(4_int64*integer_elements>huge(workspace_peak_bytes)-workspace_peak_bytes)local_bad=1
      if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+4_int64*integer_elements
      if(packet_descriptor_bytes>huge(workspace_peak_bytes)-workspace_peak_bytes)local_bad=1
      if(local_bad==0)workspace_peak_bytes=workspace_peak_bytes+packet_descriptor_bytes
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid basis workspace receipt overflow';return
    endif
    allocate(windows(nfragment,nlocal),catalog%packets(int(packet_count)),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate normalized hybrid basis outputs';return
    endif
    do i=1,nlocal
      scale=maxval(abs(raw_windows(:,i)))
      if(scale<=0d0)then;local_bad=1;cycle;endif
      scaled=raw_windows(:,i)/scale;norm_scaled=sqrt(sum(scaled**2))
      if(.not.ieee_is_finite(norm_scaled).or.norm_scaled<=0d0)then;local_bad=1;cycle;endif
      windows(:,i)=scaled/norm_scaled
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot normalize hybrid partition windows';return
    endif
    covariance_defect=0d0
    do target=1,global_row_count
      root=owner(target);remote_window=0d0
      if(rank==root)remote_window=windows(:,owner_position(target))
      call MPI_Bcast(remote_window,nfragment,MPI_DOUBLE_PRECISION,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid window row broadcast failed';return;endif
      do i=1,nlocal;do op=1,nwindow_operation
        if(row_action(int(row_ids(i)),op)/=target)cycle
        do f=1,nfragment
          covariance_defect=max(covariance_defect,&
            abs(remote_window(fragment_action(f,op))-windows(f,i)))
        enddo
      enddo;enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,covariance_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.covariance_defect>100d0*tolerance)then
      call cleanup();message='hybrid partition windows are not symmetry covariant';return
    endif
    packet_index=0
    do f=1,nfragment;do j=1,nstar
      packet_index=packet_index+1
      catalog%packets(packet_index)%fragment_id=f
      catalog%packets(packet_index)%star_id=j
      catalog%packets(packet_index)%owner_rank=mod(f-1,nproc)
      allocate(catalog%packets(packet_index)%g_indices(star_count(j)),stat=allocation_status)
      if(allocation_status/=0)then;local_bad=1;cycle;endif
      position=0
      do i=1,ng
        if(g_star(i)==j)then;position=position+1;catalog%packets(packet_index)%g_indices(position)=i;endif
      enddo
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid packet membership';return
    endif
    fingerprint=ieor(int(z'BB67AE8584CAA73B',int64),transfer(tolerance,bits))
    fingerprint=ieor(ishftc(fingerprint,7),int(nwindow_operation,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(nreciprocal_operation,int64))
    do op=1,nwindow_operation
      do f=1,nfragment
        fingerprint=ieor(ishftc(fingerprint,7),int(fragment_action(f,op),int64))
      enddo
      do i=1,global_row_count
        fingerprint=ieor(ishftc(fingerprint,7),int(row_action(i,op),int64))
      enddo
    enddo
    do op=1,nreciprocal_operation
      do i=1,ng
        fingerprint=ieor(ishftc(fingerprint,7),int(g_action(i,op),int64))
      enddo
      do j=1,3;do i=1,3
        fingerprint=ieor(ishftc(fingerprint,7),transfer(reciprocal_rotation(i,j,op),bits))
      enddo;enddo
    enddo
    do i=1,ng
      fingerprint=ieor(ishftc(fingerprint,7),int(g_star(i),int64))
      fingerprint=ieor(ishftc(fingerprint,7),int(g_conjugate(i),int64))
      do j=1,3
        fingerprint=ieor(ishftc(fingerprint,7),transfer(g_vectors(j,i),bits))
      enddo
    enddo
    do target=1,global_row_count
      root=owner(target);remote_window=0d0
      if(rank==root)remote_window=windows(:,owner_position(target))
      call MPI_Bcast(remote_window,nfragment,MPI_DOUBLE_PRECISION,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid window fingerprint broadcast failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,7),int(target,int64))
      do f=1,nfragment
        quantized=nint(remote_window(f)/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,7),quantized)
      enddo
    enddo
    do f=1,nfragment*nstar
      fingerprint=ieor(ishftc(fingerprint,7),int(catalog%packets(f)%fragment_id,int64))
      fingerprint=ieor(ishftc(fingerprint,7),int(catalog%packets(f)%star_id,int64))
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    catalog%window_fingerprint=fingerprint
    catalog%packet_fingerprint=int(z'3C6EF372FE94F82B',int64)
    do f=1,size(catalog%packets)
      catalog%packet_fingerprint=ieor(ishftc(catalog%packet_fingerprint,11),&
        int(catalog%packets(f)%fragment_id,int64))
      catalog%packet_fingerprint=ieor(ishftc(catalog%packet_fingerprint,11),&
        int(catalog%packets(f)%star_id,int64))
      do i=1,size(catalog%packets(f)%g_indices)
        catalog%packet_fingerprint=ieor(ishftc(catalog%packet_fingerprint,11),&
          int(catalog%packets(f)%g_indices(i),int64))
      enddo
    enddo
    catalog%catalog_fingerprint=ieor(catalog%window_fingerprint,ishftc(catalog%packet_fingerprint,17))
    if(catalog%catalog_fingerprint==0_int64)catalog%catalog_fingerprint=1_int64
    catalog%window_operation_count=nwindow_operation
    catalog%operation_count=nreciprocal_operation
    catalog%pw_mode_count=ng
    catalog%valid=.true.;ok=.true.
#else
    ok=.false.;message='hybrid windowed PW basis requires MPI';workspace_peak_bytes=0_int64;fingerprint=0_int64
    allocate(windows(0,0))
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup()
      if(allocated(windows))deallocate(windows)
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position)
      if(allocated(star_count))deallocate(star_count)
      if(allocated(remote_window))deallocate(remote_window)
      if(allocated(scaled))deallocate(scaled)
      if(allocated(catalog%packets))deallocate(catalog%packets)
      catalog%valid=.false.
    end subroutine cleanup
#endif
  end subroutine build_dg_hybrid_windowed_pw_basis

  subroutine materialize_dg_hybrid_windowed_pw_columns(catalog,g_vectors,coordinates,windows,&
      first_column,column_count,values,ok,message)
    type(s_dg_hybrid_basis_catalog),intent(in)::catalog
    real(real64),intent(in)::g_vectors(:,:),coordinates(:,:),windows(:,:)
    integer,intent(in)::first_column,column_count
    complex(real64),allocatable,intent(out)::values(:,:)
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::packet,g,output,fragment,allocation_status
    integer(int64)::flat,total_columns,last_column
    real(real64)::phase(size(coordinates,2)),coordinate_scale,g_scale
    ok=.false.;message='';flat=0_int64;total_columns=0_int64
    if(.not.catalog%valid.or.column_count<1.or.first_column<1.or.size(coordinates,1)/=3.or.&
      size(g_vectors,1)/=3.or.size(windows,2)/=size(coordinates,2))then
      message='invalid hybrid PW materialization contract';return
    endif
    do packet=1,size(catalog%packets)
      if(int(size(catalog%packets(packet)%g_indices),int64)>huge(total_columns)-total_columns)then
        message='hybrid PW column extent overflow';return
      endif
      total_columns=total_columns+int(size(catalog%packets(packet)%g_indices),int64)
    enddo
    if(int(column_count,int64)>huge(last_column)-int(first_column,int64))then
      message='hybrid PW requested tile extent overflow';return
    endif
    last_column=int(first_column,int64)+int(column_count,int64)-1_int64
    if(last_column>total_columns)then
      message='hybrid PW requested tile lies outside the catalog';return
    endif
    if(.not.all(ieee_is_finite(coordinates)).or..not.all(ieee_is_finite(g_vectors)).or.&
      .not.all(ieee_is_finite(windows)))then
      message='nonfinite hybrid PW materialization payload';return
    endif
    coordinate_scale=maxval(abs(coordinates));g_scale=maxval(abs(g_vectors))
    if(g_scale>0d0)then
      if(coordinate_scale>(huge(coordinate_scale)/4d0)/g_scale)then
        message='hybrid PW phase magnitude is unsafe';return
      endif
    endif
    allocate(values(column_count,size(coordinates,2)),stat=allocation_status)
    if(allocation_status/=0)then;message='cannot allocate hybrid PW tile';return;endif
    values=(0d0,0d0);output=0
    do packet=1,size(catalog%packets)
      fragment=catalog%packets(packet)%fragment_id
      if(fragment<1.or.fragment>size(windows,1))then
        deallocate(values);message='invalid hybrid packet fragment';return
      endif
      do g=1,size(catalog%packets(packet)%g_indices)
        flat=flat+1_int64
        if(flat<int(first_column,int64).or.flat>last_column)cycle
        if(catalog%packets(packet)%g_indices(g)<1.or.&
          catalog%packets(packet)%g_indices(g)>size(g_vectors,2))then
          deallocate(values);message='invalid hybrid packet reciprocal index';return
        endif
        output=output+1
        phase=matmul(g_vectors(:,catalog%packets(packet)%g_indices(g)),coordinates)
        values(output,:)=windows(fragment,:)*exp(cmplx(0d0,phase,real64))
      enddo
    enddo
    if(output/=column_count.or..not.all(ieee_is_finite(real(values))).or.&
      .not.all(ieee_is_finite(aimag(values))))then
      deallocate(values);message='hybrid PW tile range or finite check failed';return
    endif
    ok=.true.
  end subroutine materialize_dg_hybrid_windowed_pw_columns

  logical function is_permutation(values,n)
    integer,intent(in)::values(:),n
    integer::i
    integer::counts(n)
    counts=0
    do i=1,size(values)
      if(values(i)<1.or.values(i)>n)then;is_permutation=.false.;return;endif
      counts(values(i))=counts(values(i))+1
    enddo
    is_permutation=all(counts==1)
  end function is_permutation

#ifdef USE_MPI
  subroutine agree_integer(value,minimum_value,maximum_value,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum_value,maximum_value,ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer

  subroutine agree_int64(value,minimum_value,maximum_value,comm,ierr)
    integer(int64),intent(in)::value
    integer,intent(in)::comm
    integer(int64),intent(out)::minimum_value,maximum_value
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64
#endif
end module dg_hybrid_windowed_pw_basis
