#include "config.h"
module dg_hybrid_wannier_complement
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::project_dg_hybrid_wannier_complement,compute_dg_hybrid_wannier_projection_tile
contains
  subroutine compute_dg_hybrid_wannier_projection_tile(comm,global_row_count,row_ids,weights,wannier_values,&
      pw_tile,wannier_fingerprint,packet_fingerprint,first_column,tolerance,coefficients,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count,first_column
    integer(int64),intent(in)::row_ids(:),wannier_fingerprint,packet_fingerprint
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(in)::wannier_values(:,:),pw_tile(:,:)
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::nw,width,nlocal,i,j,k,ierr,local_bad,global_bad,minimum_integer,maximum_integer
    integer,allocatable::ownership(:)
    integer(int64)::minimum_bits,maximum_bits,bits
    complex(real64),allocatable::local_coefficients(:,:),gram_local(:,:),gram_global(:,:)
    real(real64)::gram_defect
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;local_bad=0
    nlocal=size(row_ids);nw=size(wannier_values,1);width=size(pw_tile,1)
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent projection row extent';return;endif
    call agree_integer(first_column,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='inconsistent projection tile origin';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='inconsistent projection tolerance';return;endif
    call agree_int64(wannier_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.wannier_fingerprint==0_int64)then
      message='invalid projection Wannier provenance';return
    endif
    call agree_int64(packet_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.packet_fingerprint==0_int64)then
      message='invalid projection packet provenance';return
    endif
    if(global_row_count<1.or.first_column<1.or.nw<1.or.width<1.or.&
        size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal.or.size(pw_tile,2)/=nlocal)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(max(0,global_row_count),int64)))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or.any(weights<=0d0).or..not.ieee_is_finite(tolerance).or.&
        tolerance<1d-15.or.tolerance>1d-2.or..not.finite_complex(wannier_values).or.&
        .not.finite_complex(pw_tile))local_bad=1
    allocate(ownership(max(1,global_row_count)));ownership=0
    do i=1,nlocal
      if(row_ids(i)>=1_int64.and.row_ids(i)<=int(max(0,global_row_count),int64))&
        ownership(int(row_ids(i)))=ownership(int(row_ids(i)))+1
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,max(0,global_row_count),MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      local_bad=1
    elseif(global_row_count>0)then
      if(any(ownership(:global_row_count)/=1))local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid projection tile contract';return;endif
    allocate(local_coefficients(nw,width),coefficients(nw,width),gram_local(nw,nw),gram_global(nw,nw))
    do j=1,width;do i=1,nw
      local_coefficients(i,j)=sum(weights*conjg(wannier_values(i,:))*pw_tile(j,:))
    enddo;enddo
    do j=1,nw;do i=1,nw
      gram_local(i,j)=sum(weights*conjg(wannier_values(i,:))*wannier_values(j,:))
    enddo;enddo
    call MPI_Allreduce(local_coefficients,coefficients,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projection coefficient reduction failed';return;endif
    call MPI_Allreduce(gram_local,gram_global,nw*nw,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='projection Gram reduction failed';return;endif
    gram_defect=0d0
    do j=1,nw;do i=1,nw
      if(i==j)then;gram_defect=max(gram_defect,abs(gram_global(i,j)-1d0))
      else;gram_defect=max(gram_defect,abs(gram_global(i,j)));endif
    enddo;enddo
    if(gram_defect>100d0*tolerance)then;message='projection Wannier frame is not orthonormal';return;endif
    workspace_peak_bytes=16_int64*int(2*nw*width+2*nw*nw,int64)+4_int64*int(global_row_count,int64)
    fingerprint=ieor(wannier_fingerprint,ishftc(packet_fingerprint,13))
    fingerprint=ieor(fingerprint,int(first_column,int64))
    do j=1,width;do i=1,nw
      bits=transfer(real(coefficients(i,j)),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=transfer(aimag(coefficients(i,j)),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
    enddo;enddo
    if(fingerprint==0_int64)fingerprint=1877_int64
    ok=.true.;message=''
#else
    ok=.false.;message='hybrid projection tile requires MPI';workspace_peak_bytes=0_int64;fingerprint=0_int64
    allocate(coefficients(0,0))
#endif
  end subroutine compute_dg_hybrid_wannier_projection_tile

  subroutine project_dg_hybrid_wannier_complement(comm,global_row_count,row_ids,weights,wannier_values,&
      pw_values,wannier_fingerprint,packet_fingerprint,packet_ids,near_offsets,near_wannier_ids,&
      diagnose_full_tail,tolerance,projected_values,&
      omitted_tail,workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,global_row_count
    integer(int64),intent(in)::row_ids(:)
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(in)::wannier_values(:,:),pw_values(:,:)
    integer(int64),intent(in)::wannier_fingerprint,packet_fingerprint
    integer,intent(in)::packet_ids(:),near_offsets(:),near_wannier_ids(:)
    logical,intent(in)::diagnose_full_tail
    complex(real64),allocatable,intent(out)::projected_values(:,:)
    real(real64),intent(out)::omitted_tail
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer,parameter::block_width=8
    integer::rank,ierr,nlocal,nw,np,i,j,j0,j1,width,p,k,local_bad,global_bad,allocation_status
    integer::minimum_integer,maximum_integer,root,position
    integer,allocatable::ownership_count(:),owner(:),owner_position(:)
    integer(int64)::bits,minimum_bits,maximum_bits,complex_elements,integer_elements,quantized
    complex(real64),allocatable::local_block(:,:),global_block(:,:),remote_row(:)
    complex(real64)::local_scalar,global_scalar
    real(real64)::local_max,global_value_scale,global_weight_scale,safe_weight,gram_defect,tail_square,&
      packet_tail_square,cross_defect,quantization_limit,value_bound_scale,global_quantization_scale
    logical::diagnose_min,diagnose_max
    ok=.false.;message='';omitted_tail=huge(1d0);workspace_peak_bytes=0_int64;fingerprint=0_int64
    local_bad=0;nlocal=size(row_ids);nw=size(wannier_values,1);np=size(pw_values,1)
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)return
    call agree_integer(global_row_count,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement spatial extent';return
    endif
    call agree_integer(nw,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement Wannier extent';return
    endif
    call agree_integer(np,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent hybrid complement PW extent';return
    endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent hybrid complement tolerance';return
    endif
    call agree_int64(wannier_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.wannier_fingerprint==0_int64)then
      message='invalid or inconsistent retained-Wannier provenance';return
    endif
    call agree_int64(packet_fingerprint,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits.or.packet_fingerprint==0_int64)then
      message='invalid or inconsistent windowed-PW packet provenance';return
    endif
    call agree_logical(diagnose_full_tail,diagnose_min,diagnose_max,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.diagnose_min.neqv.diagnose_max)then
      message='inconsistent hybrid omitted-tail diagnostic mode';return
    endif
    if(global_row_count<1.or.nw<1.or.np<1)local_bad=1
    if(size(weights)/=nlocal.or.size(wannier_values,2)/=nlocal.or.size(pw_values,2)/=nlocal)local_bad=1
    if(size(packet_ids)/=np.or.size(near_offsets)/=np+1)local_bad=1
    if(any(row_ids<1_int64).or.any(row_ids>int(global_row_count,int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance))local_bad=1
    if(.not.all(ieee_is_finite(weights)).or..not.finite_complex(wannier_values).or.&
      .not.finite_complex(pw_values))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid complement shape or finite contract';return
    endif
    if(tolerance<1d-15.or.tolerance>1d-2.or.any(weights<=0d0).or.any(packet_ids<1))local_bad=1
    do i=1,np
      call agree_integer(packet_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid PW packet IDs';return
      endif
    enddo
    do i=1,np+1
      call agree_integer(near_offsets(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid complement neighbor offsets';return
      endif
    enddo
    do i=1,size(near_wannier_ids)
      call agree_integer(near_wannier_ids(i),minimum_integer,maximum_integer,comm,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
        message='inconsistent hybrid complement neighbor IDs';return
      endif
    enddo
    if(near_offsets(1)/=1.or.near_offsets(np+1)/=size(near_wannier_ids)+1)local_bad=1
    if(any(near_offsets(2:np+1)<near_offsets(1:np)))local_bad=1
    if(any(near_wannier_ids<1).or.any(near_wannier_ids>nw))local_bad=1
    do p=1,np
      do k=near_offsets(p)+1,near_offsets(p+1)-1
        if(near_wannier_ids(k)<=near_wannier_ids(k-1))local_bad=1
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid hybrid complement sparse-neighbor catalog';return
    endif
    allocate(ownership_count(global_row_count),owner(global_row_count),owner_position(global_row_count),&
      local_block(nw,min(block_width,nw)),global_block(nw,min(block_width,nw)),&
      remote_row(np),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate hybrid complement workspace';return
    endif
    ownership_count=0;owner=-1;owner_position=0
    do i=1,nlocal
      ownership_count(int(row_ids(i)))=ownership_count(int(row_ids(i)))+1
      owner(int(row_ids(i)))=rank;owner_position(int(row_ids(i)))=i
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_row_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement ownership count failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner_position,global_row_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement position reduction failed';return;endif
    if(any(ownership_count/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='duplicate or missing hybrid complement spatial row';return
    endif
    local_max=0d0
    if(nlocal>0)local_max=max(maxval(abs(wannier_values)),maxval(abs(pw_values)))
    call MPI_Allreduce(local_max,global_value_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement value scale reduction failed';return;endif
    local_max=0d0;if(nlocal>0)local_max=maxval(weights)
    call MPI_Allreduce(local_max,global_weight_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement weight scale reduction failed';return;endif
    if(global_value_scale>0d0)then
      value_bound_scale=max(1d0,global_value_scale)
      safe_weight=huge(1d0)/64d0
      safe_weight=safe_weight/real(global_row_count,real64)/real(max(1,nw),real64)
      safe_weight=safe_weight/value_bound_scale/value_bound_scale/value_bound_scale
      if(global_weight_scale>safe_weight)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid complement multiplication magnitude is unsafe';return
    endif
    gram_defect=0d0
    do j0=1,nw,block_width
      j1=min(nw,j0+block_width-1);width=j1-j0+1;local_block(:,1:width)=(0d0,0d0)
      do j=1,width;do i=1,nw
        local_block(i,j)=sum(weights*conjg(wannier_values(i,:))*wannier_values(j0+j-1,:))
      enddo;enddo
      call MPI_Allreduce(local_block,global_block,nw*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid Wannier Gram reduction failed';return;endif
      do j=1,width;do i=1,nw
        if(i==j0+j-1)then
          gram_defect=max(gram_defect,abs(global_block(i,j)-1d0))
        else
          gram_defect=max(gram_defect,abs(global_block(i,j)))
        endif
      enddo;enddo
    enddo
    if(gram_defect>100d0*tolerance)then
      call cleanup();message='retained Wannier frame is not orthonormal';return
    endif
    complex_elements=int(np,int64)*int(nlocal,int64)+&
      2_int64*int(nw,int64)*int(min(block_width,nw),int64)+&
      int(np,int64)
    integer_elements=3_int64*int(global_row_count,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.&
      integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>&
      huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='hybrid complement workspace receipt overflow';return
    endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(projected_values(np,nlocal),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate projected hybrid PW values';return
    endif
    projected_values=pw_values;tail_square=0d0
    do p=1,np
      packet_tail_square=0d0
      if(diagnose_full_tail)then
        do j0=1,nw,block_width
          j1=min(nw,j0+block_width-1);width=j1-j0+1
          do j=1,width
            local_block(j,1)=sum(weights*conjg(wannier_values(j0+j-1,:))*pw_values(p,:))
          enddo
          call MPI_Allreduce(local_block(1:width,1),global_block(1:width,1),width,&
            MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement overlap reduction failed';return;endif
          do j=1,width
            i=j0+j-1
            if(is_near(i,p,near_offsets,near_wannier_ids))then
              projected_values(p,:)=projected_values(p,:)-wannier_values(i,:)*global_block(j,1)
            else
              packet_tail_square=packet_tail_square+abs(global_block(j,1))**2
            endif
          enddo
        enddo
        tail_square=max(tail_square,packet_tail_square)
      else
        width=near_offsets(p+1)-near_offsets(p)
        do j=1,width
          i=near_wannier_ids(near_offsets(p)+j-1)
          local_scalar=sum(weights*conjg(wannier_values(i,:))*pw_values(p,:))
          call MPI_Allreduce(local_scalar,global_scalar,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
          if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid local overlap reduction failed';return;endif
          projected_values(p,:)=projected_values(p,:)-wannier_values(i,:)*global_scalar
        enddo
      endif
    enddo
    omitted_tail=sqrt(tail_square)
    if(diagnose_full_tail.and.omitted_tail>tolerance)then
      call cleanup();message='omitted Wannier projection tail exceeds tolerance';return
    endif
    cross_defect=0d0
    do p=1,np;do i=1,nw
      local_scalar=sum(weights*conjg(wannier_values(i,:))*projected_values(p,:))
      call MPI_Allreduce(local_scalar,global_scalar,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement residual reduction failed';return;endif
      cross_defect=max(cross_defect,abs(global_scalar))
    enddo;enddo
    if(cross_defect>2d0*tolerance)then
      call cleanup();message='projected PW is not orthogonal to retained Wanniers';return
    endif
    quantization_limit=0.25d0*real(huge(0_int64),real64)*100d0*tolerance
    local_max=0d0
    if(nlocal>0)local_max=max(maxval(abs(real(projected_values))),maxval(abs(aimag(projected_values))))
    call MPI_Allreduce(local_max,global_quantization_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then
      call cleanup();message='hybrid complement quantization scale reduction failed';return
    endif
    if(global_quantization_scale>quantization_limit)then
      call cleanup();message='hybrid complement fingerprint range is unsafe';return
    endif
    fingerprint=ieor(int(z'A54FF53A5F1D36F1',int64),wannier_fingerprint)
    fingerprint=ieor(ishftc(fingerprint,9),packet_fingerprint)
    do i=1,np
      fingerprint=ieor(ishftc(fingerprint,9),int(packet_ids(i),int64))
      fingerprint=ieor(ishftc(fingerprint,9),int(near_offsets(i+1)-near_offsets(i),int64))
      do k=near_offsets(i),near_offsets(i+1)-1
        fingerprint=ieor(ishftc(fingerprint,9),int(near_wannier_ids(k),int64))
      enddo
    enddo
    do i=1,global_row_count
      root=owner(i);remote_row=(0d0,0d0)
      if(rank==root)remote_row=projected_values(:,owner_position(i))
      call MPI_Bcast(remote_row,np,MPI_DOUBLE_COMPLEX,root,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup();message='hybrid complement fingerprint broadcast failed';return;endif
      fingerprint=ieor(ishftc(fingerprint,9),int(i,int64))
      do p=1,np
        quantized=nint(real(remote_row(p))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
        quantized=nint(aimag(remote_row(p))/(100d0*tolerance),int64)
        fingerprint=ieor(ishftc(fingerprint,9),quantized)
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    ok=.true.
#else
    ok=.false.;message='hybrid Wannier complement requires MPI';omitted_tail=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;allocate(projected_values(0,0))
#endif
  contains
#ifdef USE_MPI
    subroutine cleanup()
      if(allocated(projected_values))deallocate(projected_values)
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(owner))deallocate(owner)
      if(allocated(owner_position))deallocate(owner_position)
      if(allocated(local_block))deallocate(local_block)
      if(allocated(global_block))deallocate(global_block)
      if(allocated(remote_row))deallocate(remote_row)
    end subroutine cleanup
#endif
  end subroutine project_dg_hybrid_wannier_complement

  logical function is_near(wannier_id,pw_id,offsets,ids)
    integer,intent(in)::wannier_id,pw_id,offsets(:),ids(:)
    integer::k
    is_near=.false.
    do k=offsets(pw_id),offsets(pw_id+1)-1
      if(ids(k)==wannier_id)then;is_near=.true.;return;endif
    enddo
  end function is_near

  logical function finite_complex(values)
    complex(real64),intent(in)::values(:,:)
    finite_complex=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
  end function finite_complex

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
  subroutine agree_logical(value,minimum_value,maximum_value,comm,ierr)
    logical,intent(in)::value
    logical,intent(out)::minimum_value,maximum_value
    integer,intent(in)::comm
    integer,intent(out)::ierr
    integer::input,minimum_integer,maximum_integer
    input=merge(1,0,value)
    call agree_integer(input,minimum_integer,maximum_integer,comm,ierr)
    minimum_value=minimum_integer==1;maximum_value=maximum_integer==1
  end subroutine agree_logical
#endif
end module dg_hybrid_wannier_complement
