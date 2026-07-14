module lcfo_wannier_sawf_dmn
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: int64, file_storage_size
  use lcfo_wannier_sawf, only: t_sawf_symop
  use lcfo_wannier_sawf_band, only: validate_sawf_dmn_band
  implicit none
  private

  type, public :: t_sawf_dmn_writer
    character(:), allocatable :: final_path, text_path, wann_scratch_path, band_scratch_path
    integer :: text_unit=0, wann_unit=0, band_unit=0
    integer :: num_bands=0, num_wann=0, num_symmetry=0, appended=0
    real(8) :: tolerance=0.0d0
    real(8) :: hamiltonian_tolerance=0.0d0
    real(8) :: max_unitarity=0.0d0, max_hamiltonian=0.0d0, max_amn=0.0d0
    real(8) :: max_group_wann=0.0d0, max_group_band=0.0d0
    logical :: active=.false.
    logical :: owns_text=.false., owns_wann_scratch=.false., owns_band_scratch=.false.
  end type t_sawf_dmn_writer

  type, public :: t_sawf_operation_index
    integer, allocatable :: slots(:)
    integer(int64), allocatable :: fingerprints(:)
    real(8) :: tolerance=0.0d0
  end type t_sawf_operation_index

  integer, parameter :: max_sawf_index_operations=4096
  integer, parameter :: max_temp_attempts=128
  integer(int64), parameter :: hash_modulus=2147483629_int64
  integer(int64), parameter :: hash_factor=65599_int64

  public :: begin_sawf_dmn, append_sawf_dmn_operation, finish_sawf_dmn, abort_sawf_dmn
  public :: validate_sawf_dmn_covariances
  public :: build_sawf_operation_index, lookup_sawf_operation_product

  interface
    integer(c_int) function c_rename(old_path,new_path) bind(C,name='rename')
      import :: c_char,c_int
      character(c_char), intent(in) :: old_path(*),new_path(*)
    end function c_rename
    integer(c_int) function c_remove(path) bind(C,name='remove')
      import :: c_char,c_int
      character(c_char), intent(in) :: path(*)
    end function c_remove
  end interface

contains

  subroutine begin_sawf_dmn(writer,final_path,num_bands,num_wann,num_symmetry,tolerance,ok,message, &
      temp_nonce,hamiltonian_tolerance)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    character(*), intent(in) :: final_path
    integer, intent(in) :: num_bands,num_wann,num_symmetry
    real(8), intent(in) :: tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, intent(in), optional :: temp_nonce
    real(8), intent(in), optional :: hamiltonian_tolerance
    integer :: clock_count,io_status,attempt
    integer(int64) :: nonce
    character(256) :: suffix,io_message
    integer(int64) :: matrix_elements

    call abort_sawf_dmn(writer)
    ok=.false.; message=''
    if(len_trim(final_path)==0 .or. num_bands<=0 .or. num_wann<=0 .or. num_symmetry<=0 .or. &
        .not.ieee_is_finite(tolerance) .or. tolerance<=0.0d0) then
      message='SAWF DMN writer dimensions, path, or tolerance are invalid'
      return
    end if
    matrix_elements=int(num_bands,int64)*int(num_bands,int64)
    if(matrix_elements>int(huge(0),int64) .or. &
        int(num_wann,int64)*int(num_wann,int64)>int(huge(0),int64)) then
      message='SAWF DMN matrix dimensions exceed addressable default-integer indexing'
      return
    end if
    call system_clock(clock_count)
    nonce=int(abs(clock_count),int64)
    if(present(temp_nonce)) nonce=int(abs(temp_nonce),int64)
    writer%final_path=trim(final_path)
    do attempt=0,max_temp_attempts-1
      write(suffix,'(a,i0,a,i0)') '.tmp.',nonce,'.',attempt
      writer%text_path=trim(final_path)//trim(suffix)
      writer%wann_scratch_path=trim(final_path)//trim(suffix)//'.wann.scratch'
      writer%band_scratch_path=trim(final_path)//trim(suffix)//'.band.scratch'
      writer%text_unit=0
      open(newunit=writer%text_unit,file=writer%text_path,status='new',action='write', &
        form='formatted',iostat=io_status,iomsg=io_message)
      if(io_status/=0) then
        writer%text_unit=0
        cycle
      end if
      writer%owns_text=.true.
      writer%wann_unit=0
      open(newunit=writer%wann_unit,file=writer%wann_scratch_path,status='new',action='readwrite', &
        access='stream',form='unformatted',iostat=io_status,iomsg=io_message)
      if(io_status==0) writer%owns_wann_scratch=.true.
      if(io_status/=0) writer%wann_unit=0
      if(io_status==0) then
        writer%band_unit=0
        open(newunit=writer%band_unit,file=writer%band_scratch_path,status='new',action='readwrite', &
          access='stream',form='unformatted',iostat=io_status,iomsg=io_message)
        if(io_status==0) writer%owns_band_scratch=.true.
        if(io_status/=0) writer%band_unit=0
      end if
      if(io_status==0) exit
      call release_owned_temps(writer)
    end do
    if(io_status/=0) then
      message='SAWF DMN could not exclusively create temporary files after bounded retries: '//trim(io_message)
      call abort_sawf_dmn(writer)
      return
    end if
    writer%num_bands=num_bands; writer%num_wann=num_wann
    writer%num_symmetry=num_symmetry; writer%tolerance=tolerance; writer%appended=0
    writer%hamiltonian_tolerance=tolerance
    if(present(hamiltonian_tolerance)) then
      if(.not.ieee_is_finite(hamiltonian_tolerance) .or. hamiltonian_tolerance<tolerance) then
        message='SAWF DMN Hamiltonian tolerance must be finite and no smaller than closure tolerance'
        call abort_sawf_dmn(writer)
        return
      end if
      writer%hamiltonian_tolerance=hamiltonian_tolerance
    end if
    writer%max_unitarity=0.0d0; writer%max_hamiltonian=0.0d0; writer%max_amn=0.0d0
    writer%max_group_wann=0.0d0; writer%max_group_band=0.0d0; writer%active=.true.
    write(writer%text_unit,'(a)',iostat=io_status,iomsg=io_message) &
      'SALMON SAWF Gamma-only symmetry data'
    if(io_status==0) write(writer%text_unit,'(4i9)',iostat=io_status,iomsg=io_message) &
      num_bands,num_symmetry,1,1
    if(io_status==0) write(writer%text_unit,*,iostat=io_status,iomsg=io_message)
    if(io_status==0) write(writer%text_unit,'(10i9)',iostat=io_status,iomsg=io_message) 1
    if(io_status==0) write(writer%text_unit,*,iostat=io_status,iomsg=io_message)
    if(io_status==0) write(writer%text_unit,'(10i9)',iostat=io_status,iomsg=io_message) 1
    if(io_status==0) write(writer%text_unit,*,iostat=io_status,iomsg=io_message)
    if(io_status==0) write(writer%text_unit,'(10i9)',iostat=io_status,iomsg=io_message) &
      [(1,clock_count=1,num_symmetry)]
    if(io_status/=0) then
      message='SAWF DMN header write failed: '//trim(io_message)
      call abort_sawf_dmn(writer)
      return
    end if
    ok=.true.
  end subroutine begin_sawf_dmn


  subroutine append_sawf_dmn_operation(writer,operation_index,d_wann,d_band,eigenvalues,amn, &
      is_identity,ok,message,singular_min,singular_max,closure_residual)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    integer, intent(in) :: operation_index
    complex(8), intent(in) :: d_wann(:,:),d_band(:,:),amn(:,:)
    real(8), intent(in) :: eigenvalues(:)
    logical, intent(in) :: is_identity
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), intent(out), optional :: singular_min,singular_max,closure_residual
    real(8) :: smin,smax,closure,unitarity,ham_residual,amn_residual,identity_residual
    integer :: io_status
    character(256) :: io_message
    logical :: covariance_ok

    ok=.false.; message=''
    if(.not.writer%active .or. operation_index/=writer%appended+1 .or. &
        operation_index<1 .or. operation_index>writer%num_symmetry) then
      message='SAWF DMN operations must be appended once in normalized order'
      return
    end if
    if(any(shape(d_wann)/=[writer%num_wann,writer%num_wann]) .or. &
        any(shape(d_band)/=[writer%num_bands,writer%num_bands]) .or. &
        size(eigenvalues)/=writer%num_bands .or. &
        any(shape(amn)/=[writer%num_bands,writer%num_wann])) then
      message='SAWF DMN operation dimensions are inconsistent with the header'
      return
    end if
    call validate_sawf_dmn_band(d_band,writer%tolerance,smin,smax,closure,ok,message)
    if(present(singular_min)) singular_min=smin
    if(present(singular_max)) singular_max=smax
    if(present(closure_residual)) closure_residual=closure
    if(.not.ok) return
    call validate_sawf_dmn_covariances(d_wann,d_band,eigenvalues,amn,writer%tolerance, &
      covariance_ok,message,unitarity,ham_residual,amn_residual, &
      writer%hamiltonian_tolerance)
    if(.not.covariance_ok) then
      ok=.false.
      return
    end if
    identity_residual=0.0d0
    if(is_identity) then
      identity_residual=max(identity_matrix_residual(d_wann),identity_matrix_residual(d_band))
      if(.not.ieee_is_finite(identity_residual) .or. identity_residual>writer%tolerance) then
        write(message,'(a,es13.5)') 'SAWF DMN identity operation residual=',identity_residual
        return
      end if
      if(operation_index/=1) then
        message='SAWF DMN identity operation must be first'
        return
      end if
    else if(operation_index==1) then
      message='SAWF DMN first operation must be identity'
      return
    end if
    write(writer%text_unit,*,iostat=io_status,iomsg=io_message)
    if(io_status==0) write(writer%text_unit,"(1p,(' (',e18.10,',',e18.10,')'))", &
      iostat=io_status,iomsg=io_message) d_wann
    if(io_status==0) write(writer%wann_unit,iostat=io_status,iomsg=io_message) d_wann
    if(io_status==0) write(writer%band_unit,iostat=io_status,iomsg=io_message) d_band
    if(io_status/=0) then
      message='SAWF DMN streamed operation write failed: '//trim(io_message)
      return
    end if
    writer%appended=operation_index
    writer%max_unitarity=max(writer%max_unitarity,unitarity,closure)
    writer%max_hamiltonian=max(writer%max_hamiltonian,ham_residual)
    writer%max_amn=max(writer%max_amn,amn_residual)
    ok=.true.
  end subroutine append_sawf_dmn_operation


  subroutine validate_sawf_dmn_covariances(d_wann,d_band,eigenvalues,amn,tolerance,ok,message, &
      unitarity_residual,hamiltonian_residual,amn_residual,hamiltonian_tolerance)
    complex(8), intent(in) :: d_wann(:,:),d_band(:,:),amn(:,:)
    real(8), intent(in) :: eigenvalues(:),tolerance
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8), intent(out), optional :: unitarity_residual,hamiltonian_residual,amn_residual
    real(8), intent(in), optional :: hamiltonian_tolerance
    complex(8), allocatable :: metric(:,:),left(:,:),right(:,:),weighted(:,:)
    real(8) :: ures,hres,ares,energy_scale,amn_scale,energy_center,energy_min,energy_max
    real(8) :: hamiltonian_tolerance_effective
    real(8), allocatable :: centered_energy(:)
    integer :: allocation_status,i,nb,nw

    ok=.false.; message=''; ures=huge(0.0d0); hres=huge(0.0d0); ares=huge(0.0d0)
    hamiltonian_tolerance_effective=tolerance
    if(present(hamiltonian_tolerance))hamiltonian_tolerance_effective=hamiltonian_tolerance
    nb=size(d_band,1); nw=size(d_wann,1)
    if(nb<=0 .or. nw<=0 .or. size(d_band,2)/=nb .or. size(d_wann,2)/=nw .or. &
        size(eigenvalues)/=nb .or. any(shape(amn)/=[nb,nw]) .or. &
        .not.ieee_is_finite(tolerance) .or. tolerance<=0.0d0 .or. &
        .not.ieee_is_finite(hamiltonian_tolerance_effective) .or. &
        hamiltonian_tolerance_effective<tolerance .or. &
        .not.all(ieee_is_finite(eigenvalues)) .or. .not.complex_finite(d_wann) .or. &
        .not.complex_finite(d_band) .or. .not.complex_finite(amn)) then
      message='SAWF DMN covariance inputs are invalid or non-finite'
      call set_optional_residuals()
      return
    end if
    allocate(metric(max(nb,nw),max(nb,nw)),left(nb,nw),right(nb,nw), &
      weighted(nb,nb),centered_energy(nb),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF DMN covariance work allocation failed'
      call set_optional_residuals()
      return
    end if
    metric=(0.0d0,0.0d0)
    metric(1:nw,1:nw)=matmul(conjg(transpose(d_wann)),d_wann)
    do i=1,nw; metric(i,i)=metric(i,i)-(1.0d0,0.0d0); enddo
    ures=maxval(abs(metric(1:nw,1:nw)))
    metric(1:nb,1:nb)=matmul(conjg(transpose(d_band)),d_band)
    do i=1,nb; metric(i,i)=metric(i,i)-(1.0d0,0.0d0); enddo
    ures=max(ures,maxval(abs(metric(1:nb,1:nb))))
    if(ures>2.0d0*tolerance) then
      write(message,'(a,es13.5)') 'SAWF DMN unitarity residual=',ures
      deallocate(metric,left,right,weighted,centered_energy); call set_optional_residuals(); return
    end if
    energy_min=minval(eigenvalues); energy_max=maxval(eigenvalues)
    energy_center=energy_min+0.5d0*(energy_max-energy_min)
    centered_energy=eigenvalues-energy_center
    energy_scale=maxval(abs(centered_energy))
    weighted=(0.0d0,0.0d0)
    do i=1,nb; weighted(i,:)=centered_energy(i)*d_band(i,:); enddo
    weighted=matmul(conjg(transpose(d_band)),weighted)
    do i=1,nb; weighted(i,i)=weighted(i,i)-centered_energy(i); enddo
    if(energy_scale==0.0d0) then
      hres=0.0d0
    else
      hres=maxval(abs(weighted))/energy_scale
    end if
    if(hres>hamiltonian_tolerance_effective) then
      write(message,'(a,es13.5,a,es13.5)') 'SAWF DMN Hamiltonian covariance residual=', &
        hres,' tolerance=',hamiltonian_tolerance_effective
      deallocate(metric,left,right,weighted,centered_energy); call set_optional_residuals(); return
    end if
    ! Wannier90 sitesym equation at Gamma: U=d_band U D_wann^dag,
    ! equivalently d_band*A=A*D_wann for the AMN trial matrix A.
    left=matmul(d_band,amn); right=matmul(amn,d_wann)
    amn_scale=max(1.0d0,maxval(abs(amn)))
    ares=maxval(abs(left-right))/amn_scale
    if(ares>2.0d0*tolerance) then
      write(message,'(a,es13.5,a,es13.5)') 'SAWF DMN AMN covariance residual=', &
        ares,' tolerance=',2.0d0*tolerance
      deallocate(metric,left,right,weighted,centered_energy); call set_optional_residuals(); return
    end if
    deallocate(metric,left,right,weighted,centered_energy)
    ok=.true.; call set_optional_residuals()
  contains
    subroutine set_optional_residuals()
      if(present(unitarity_residual)) unitarity_residual=ures
      if(present(hamiltonian_residual)) hamiltonian_residual=hres
      if(present(amn_residual)) amn_residual=ares
    end subroutine
  end subroutine validate_sawf_dmn_covariances


  subroutine finish_sawf_dmn(writer,operations,ok,message)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    type(t_sawf_symop), intent(in) :: operations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: matrix(:,:)
    integer :: iop,io_status,close_status,rename_status,allocation_status
    character(256) :: io_message

    ok=.false.; message=''
    if(.not.writer%active .or. writer%appended/=writer%num_symmetry .or. &
        size(operations)/=writer%num_symmetry) then
      message='SAWF DMN cannot publish before every normalized operation is appended'
      return
    end if
    call validate_streamed_group(writer,operations,ok,message)
    if(.not.ok) return
    rewind(writer%band_unit)
    allocate(matrix(writer%num_bands,writer%num_bands),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF DMN final stream buffer allocation failed'
      return
    end if
    do iop=1,writer%num_symmetry
      read(writer%band_unit,iostat=io_status,iomsg=io_message) matrix
      if(io_status==0) write(writer%text_unit,*,iostat=io_status,iomsg=io_message)
      if(io_status==0) write(writer%text_unit,"(1p,(' (',e18.10,',',e18.10,')'))", &
        iostat=io_status,iomsg=io_message) matrix
      if(io_status/=0) then
        message='SAWF DMN final band stream failed: '//trim(io_message)
        deallocate(matrix)
        return
      end if
    end do
    deallocate(matrix)
    flush(writer%text_unit,iostat=io_status,iomsg=io_message)
    close_status=0
    close(writer%text_unit,iostat=close_status)
    writer%text_unit=0
    if(io_status/=0 .or. close_status/=0) then
      message='SAWF DMN temporary file flush/close failed: '//trim(io_message)
      return
    end if
    close(writer%wann_unit); writer%wann_unit=0
    close(writer%band_unit); writer%band_unit=0
    rename_status=rename_path(writer%text_path,writer%final_path)
    if(rename_status/=0) then
      write(message,'(a,i0)') 'SAWF DMN atomic same-directory rename failed; status=',rename_status
      call remove_path(writer%text_path); call remove_path(writer%wann_scratch_path)
      call remove_path(writer%band_scratch_path); writer%active=.false.
      return
    end if
    if(writer%owns_wann_scratch) call remove_path(writer%wann_scratch_path)
    if(writer%owns_band_scratch) call remove_path(writer%band_scratch_path)
    writer%owns_wann_scratch=.false.; writer%owns_band_scratch=.false.; writer%owns_text=.false.
    writer%active=.false.; ok=.true.
  end subroutine finish_sawf_dmn


  subroutine validate_streamed_group(writer,operations,ok,message)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    type(t_sawf_symop), intent(in) :: operations(:)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer, allocatable :: generators(:)
    logical, allocatable :: reached(:)
    integer :: ngen,iop,igen,product,allocation_status,io_status
    logical :: changed
    real(8) :: residual
    type(t_sawf_operation_index) :: operation_index

    ok=.false.; message=''
    flush(writer%wann_unit,iostat=io_status)
    if(io_status==0) flush(writer%band_unit,iostat=io_status)
    if(io_status/=0) then
      message='SAWF DMN group scratch flush failed'
      return
    end if
    allocate(generators(size(operations)),reached(size(operations)),stat=allocation_status)
    if(allocation_status/=0) then; message='SAWF DMN group traversal allocation failed'; return; endif
    call build_sawf_operation_index(operations,writer%tolerance,operation_index,ok,message)
    if(.not.ok) return
    reached=.false.; reached(1)=.true.; ngen=0
    do while(.not.all(reached))
      do iop=2,size(operations)
        if(.not.reached(iop)) exit
      end do
      ngen=ngen+1; generators(ngen)=iop; changed=.true.
      do while(changed)
        changed=.false.
        do iop=1,size(operations)
          if(.not.reached(iop)) cycle
          do igen=1,ngen
            call lookup_sawf_operation_product(operation_index,operations,generators(igen),iop, &
              product,ok,message)
            if(.not.ok) return
            if(.not.reached(product)) then; reached(product)=.true.; changed=.true.; endif
            call lookup_sawf_operation_product(operation_index,operations,iop,generators(igen), &
              product,ok,message)
            if(.not.ok) return
            if(.not.reached(product)) then; reached(product)=.true.; changed=.true.; endif
          end do
        end do
      end do
    end do
    ! Exact generator edges include spanning-tree edges and every cycle relation.
    ! Since the selected generators reach every operation, these equations establish
    ! the representation on the whole finite group without an O(nsym^2) matrix table.
    do igen=1,ngen
      do iop=1,size(operations)
        call lookup_sawf_operation_product(operation_index,operations,generators(igen),iop, &
          product,ok,message)
        if(.not.ok) return
        call validate_scratch_relation(writer,.true.,generators(igen),iop,product,residual,ok,message)
        writer%max_group_wann=max(writer%max_group_wann,residual)
        if(.not.ok) return
        call validate_scratch_relation(writer,.false.,generators(igen),iop,product,residual,ok,message)
        writer%max_group_band=max(writer%max_group_band,residual)
        if(.not.ok) return
        call lookup_sawf_operation_product(operation_index,operations,iop,generators(igen), &
          product,ok,message)
        if(.not.ok) return
        call validate_scratch_relation(writer,.true.,iop,generators(igen),product,residual,ok,message)
        writer%max_group_wann=max(writer%max_group_wann,residual)
        if(.not.ok) return
        call validate_scratch_relation(writer,.false.,iop,generators(igen),product,residual,ok,message)
        writer%max_group_band=max(writer%max_group_band,residual)
        if(.not.ok) return
      end do
    end do
    deallocate(generators,reached); ok=.true.
  end subroutine validate_streamed_group


  subroutine validate_scratch_relation(writer,wannier_space,left,right,product,residual,ok,message)
    type(t_sawf_dmn_writer), intent(in) :: writer
    logical, intent(in) :: wannier_space
    integer, intent(in) :: left,right,product
    real(8), intent(out) :: residual
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    complex(8), allocatable :: a(:,:),b(:,:),c(:,:),work(:,:)
    integer :: n,unit,allocation_status,io_status
    integer(int64) :: position
    character(256) :: io_message

    ok=.false.; message=''; residual=huge(0.0d0)
    if(product<1) then; message='SAWF DMN group product has no normalized operation'; return; endif
    if(wannier_space) then; n=writer%num_wann; unit=writer%wann_unit
    else; n=writer%num_bands; unit=writer%band_unit; endif
    allocate(a(n,n),b(n,n),c(n,n),work(n,n),stat=allocation_status)
    if(allocation_status/=0) then; message='SAWF DMN group work allocation failed'; return; endif
    call checked_stream_position(n,left,position,ok,message); if(.not.ok) return
    read(unit,pos=position,iostat=io_status,iomsg=io_message) a
    if(io_status==0) then
      call checked_stream_position(n,right,position,ok,message); if(.not.ok) return
      read(unit,pos=position,iostat=io_status,iomsg=io_message) b
    endif
    if(io_status==0) then
      call checked_stream_position(n,product,position,ok,message); if(.not.ok) return
      read(unit,pos=position,iostat=io_status,iomsg=io_message) c
    endif
    if(io_status/=0) then; message='SAWF DMN group scratch read failed: '//trim(io_message); deallocate(a,b,c,work); return; endif
    work=matmul(a,b)-c
    residual=maxval(abs(work))
    ok=.false.
    if(.not.ieee_is_finite(residual) .or. residual>4.0d0*writer%tolerance) then
      write(message,'(a,3(i0,a),es13.5)') 'SAWF DMN group relation left=',left, &
        ' right=',right,' product=',product,' residual=',residual
      deallocate(a,b,c,work); return
    end if
    deallocate(a,b,c,work); ok=.true.
  end subroutine validate_scratch_relation


  subroutine build_sawf_operation_index(operations,tolerance,index,ok,message)
    type(t_sawf_symop), intent(in) :: operations(:)
    real(8), intent(in) :: tolerance
    type(t_sawf_operation_index), intent(inout) :: index
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: table_size,allocation_status,iop,slot,step,match_count,matched
    integer(int64) :: fingerprint,qtau(3)

    ok=.false.; message=''
    if(allocated(index%slots)) deallocate(index%slots)
    if(allocated(index%fingerprints)) deallocate(index%fingerprints)
    if(size(operations)<1 .or. size(operations)>max_sawf_index_operations .or. &
        .not.ieee_is_finite(tolerance) .or. tolerance<=0.0d0) then
      message='SAWF operation index size or tolerance is invalid or exceeds its cap'
      return
    end if
    table_size=1
    do while(table_size<4*size(operations))
      if(table_size>huge(table_size)/2) then
        message='SAWF operation hash table size overflows default integer'
        return
      end if
      table_size=2*table_size
    end do
    allocate(index%slots(table_size),index%fingerprints(table_size),stat=allocation_status)
    if(allocation_status/=0) then
      message='SAWF operation hash table allocation failed'
      return
    end if
    index%slots=0; index%fingerprints=0_int64; index%tolerance=tolerance
    do iop=1,size(operations)
      if(.not.allocated(operations(iop)%atom_map) .or. &
          .not.all(ieee_is_finite(operations(iop)%tau))) then
        write(message,'(a,i0,a)') 'SAWF operation ',iop,' has no finite full atom-map key'
        return
      end if
      call quantize_periodic_tau(operations(iop)%tau,tolerance,qtau,ok,message)
      if(.not.ok) return
      call find_index_matches(index,operations,operations(iop)%W,operations(iop)%tau, &
        operations(iop)%atom_map,matched,match_count,ok,message)
      if(.not.ok) return
      if(match_count/=0) then
        write(message,'(a,i0,a,i0)') 'SAWF operation index duplicate full key: operation=', &
          iop,' existing=',matched
        ok=.false.
        return
      end if
      fingerprint=operation_fingerprint(operations(iop)%W,operations(iop)%atom_map,qtau)
      slot=1+int(modulo(fingerprint,int(table_size,int64)))
      do step=0,table_size-1
        if(index%slots(slot)==0) then
          index%slots(slot)=iop; index%fingerprints(slot)=fingerprint
          exit
        end if
        slot=1+modulo(slot,table_size)
      end do
      if(step==table_size) then
        message='SAWF operation hash table unexpectedly filled'
        return
      end if
    end do
    ok=.true.
  end subroutine build_sawf_operation_index


  subroutine lookup_sawf_operation_product(index,operations,left,right,product,ok,message)
    type(t_sawf_operation_index), intent(in) :: index
    type(t_sawf_symop), intent(in) :: operations(:)
    integer, intent(in) :: left,right
    integer, intent(out) :: product
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer :: i,j,k,match_count,natom
    integer(int64) :: term,sum64,new_sum
    integer :: product_w(3,3)
    integer, allocatable :: product_map(:)
    real(8) :: product_tau(3)

    ok=.false.; message=''; product=0
    if(.not.allocated(index%slots) .or. left<1 .or. right<1 .or. &
        left>size(operations) .or. right>size(operations)) then
      message='SAWF operation product lookup index or operands are invalid'
      return
    end if
    if(.not.allocated(operations(left)%atom_map) .or. .not.allocated(operations(right)%atom_map) .or. &
        size(operations(left)%atom_map)/=size(operations(right)%atom_map)) then
      message='SAWF operation product atom-map dimensions are inconsistent'
      return
    end if
    do i=1,3
      do j=1,3
        sum64=0_int64
        do k=1,3
          call checked_int64_multiply(int(operations(left)%W(i,k),int64), &
            int(operations(right)%W(k,j),int64),term,ok)
          if(.not.ok) then; message='SAWF operation product rotation multiplication overflows int64'; return; endif
          call checked_int64_add(sum64,term,new_sum,ok)
          if(.not.ok) then; message='SAWF operation product rotation sum overflows int64'; return; endif
          sum64=new_sum
        end do
        if(sum64<int(-huge(0),int64) .or. sum64>int(huge(0),int64)) then
          message='SAWF operation product rotation exceeds default integer'
          return
        end if
        product_w(i,j)=int(sum64)
      end do
    end do
    product_tau=matmul(real(operations(left)%W,8),operations(right)%tau)+operations(left)%tau
    if(.not.all(ieee_is_finite(product_tau))) then
      message='SAWF operation product translation is non-finite'
      return
    end if
    natom=size(operations(left)%atom_map)
    allocate(product_map(natom))
    do i=1,natom
      j=operations(right)%atom_map(i)
      if(j<1 .or. j>natom) then
        message='SAWF right operation atom map is out of range'; deallocate(product_map); return
      end if
      product_map(i)=operations(left)%atom_map(j)
      if(product_map(i)<1 .or. product_map(i)>natom) then
        message='SAWF composed operation atom map is out of range'; deallocate(product_map); return
      end if
    end do
    call find_index_matches(index,operations,product_w,product_tau,product_map, &
      product,match_count,ok,message)
    deallocate(product_map)
    if(.not.ok) return
    if(match_count==0) then
      message='SAWF normalized operation set is not closed for the full W/tau/atom-map key'
      ok=.false.; return
    else if(match_count>1) then
      message='SAWF operation product lookup is ambiguous for the full key'
      ok=.false.; return
    end if
    ok=.true.
  end subroutine lookup_sawf_operation_product


  subroutine find_index_matches(index,operations,w,tau,atom_map,matched,match_count,ok,message)
    type(t_sawf_operation_index), intent(in) :: index
    type(t_sawf_symop), intent(in) :: operations(:)
    integer, intent(in) :: w(3,3),atom_map(:)
    real(8), intent(in) :: tau(3)
    integer, intent(out) :: matched,match_count
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: qtau(3),trial_q(3),fingerprint
    integer :: dx,dy,dz,slot,step,candidate,table_size

    ok=.false.; message=''; matched=0; match_count=0
    if(.not.allocated(index%slots) .or. .not.allocated(index%fingerprints)) then
      message='SAWF operation hash index is not allocated'
      return
    end if
    call quantize_periodic_tau(tau,index%tolerance,qtau,ok,message)
    if(.not.ok) return
    table_size=size(index%slots)
    do dz=-1,1; do dy=-1,1; do dx=-1,1
      trial_q=qtau+[int(dx,int64),int(dy,int64),int(dz,int64)]
      fingerprint=operation_fingerprint(w,atom_map,trial_q)
      slot=1+int(modulo(fingerprint,int(table_size,int64)))
      do step=0,table_size-1
        candidate=index%slots(slot)
        if(candidate==0) exit
        if(index%fingerprints(slot)==fingerprint) then
          if(operation_key_equal(operations(candidate),w,tau,atom_map,index%tolerance)) then
            if(match_count==0) then
              matched=candidate; match_count=1
            else if(candidate/=matched) then
              match_count=2
              ok=.true.
              return
            end if
          end if
        end if
        slot=1+modulo(slot,table_size)
      end do
    end do; end do; end do
    ok=.true.
  end subroutine find_index_matches


  logical function operation_key_equal(operation,w,tau,atom_map,tolerance) result(equal)
    type(t_sawf_symop), intent(in) :: operation
    integer, intent(in) :: w(3,3),atom_map(:)
    real(8), intent(in) :: tau(3),tolerance
    real(8) :: delta(3)
    equal=.false.
    if(any(operation%W/=w) .or. .not.allocated(operation%atom_map) .or. &
        size(operation%atom_map)/=size(atom_map)) return
    if(any(operation%atom_map/=atom_map)) return
    delta=operation%tau-tau; delta=delta-anint(delta)
    equal=maxval(abs(delta))<=tolerance
  end function operation_key_equal


  subroutine quantize_periodic_tau(tau,tolerance,qtau,ok,message)
    real(8), intent(in) :: tau(3),tolerance
    integer(int64), intent(out) :: qtau(3)
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    real(8) :: wrapped,scaled
    integer :: i
    ok=.false.; message=''; qtau=0_int64
    if(.not.all(ieee_is_finite(tau)) .or. .not.ieee_is_finite(tolerance) .or. tolerance<=0.0d0) then
      message='SAWF operation translation quantization input is invalid'
      return
    end if
    do i=1,3
      wrapped=modulo(tau(i),1.0d0)
      if(1.0d0-wrapped<=tolerance) wrapped=0.0d0
      scaled=wrapped/tolerance
      if(.not.ieee_is_finite(scaled) .or. scaled>dble(huge(0_int64)-2_int64)) then
        message='SAWF operation translation quantization exceeds int64'
        return
      end if
      qtau(i)=nint(scaled,kind=int64)
    end do
    ok=.true.
  end subroutine quantize_periodic_tau


  integer(int64) function operation_fingerprint(w,atom_map,qtau) result(hash)
    integer, intent(in) :: w(3,3),atom_map(:)
    integer(int64), intent(in) :: qtau(3)
    integer :: i,j
    hash=17_int64
    do j=1,3; do i=1,3; call hash_mix(hash,int(w(i,j),int64)); end do; end do
    do i=1,3; call hash_mix(hash,qtau(i)); end do
    call hash_mix(hash,int(size(atom_map),int64))
    do i=1,size(atom_map); call hash_mix(hash,int(atom_map(i),int64)); end do
  end function operation_fingerprint

  subroutine hash_mix(hash,value)
    integer(int64), intent(inout) :: hash
    integer(int64), intent(in) :: value
    hash=modulo(hash*hash_factor+modulo(value,hash_modulus),hash_modulus)
  end subroutine hash_mix

  subroutine checked_int64_multiply(a,b,result,ok)
    integer(int64), intent(in) :: a,b
    integer(int64), intent(out) :: result
    logical, intent(out) :: ok
    result=0_int64; ok=.false.
    if(a==0_int64 .or. b==0_int64) then; ok=.true.; return; endif
    if(a==-1_int64 .and. b==-huge(0_int64)) return
    if(b==-1_int64 .and. a==-huge(0_int64)) return
    if(abs(a)>huge(0_int64)/abs(b)) return
    result=a*b; ok=.true.
  end subroutine checked_int64_multiply

  subroutine checked_int64_add(a,b,result,ok)
    integer(int64), intent(in) :: a,b
    integer(int64), intent(out) :: result
    logical, intent(out) :: ok
    result=0_int64; ok=.false.
    if((b>0_int64 .and. a>huge(0_int64)-b) .or. &
        (b<0_int64 .and. a< -huge(0_int64)-b)) return
    result=a+b; ok=.true.
  end subroutine checked_int64_add

  subroutine checked_stream_position(matrix_order,operation_index,position,ok,message)
    integer, intent(in) :: matrix_order,operation_index
    integer(int64), intent(out) :: position
    logical, intent(out) :: ok
    character(*), intent(out) :: message
    integer(int64) :: element_units,matrix_units,offset
    position=0_int64; ok=.false.; message=''
    if(matrix_order<=0 .or. operation_index<=0 .or. file_storage_size<=0 .or. &
        modulo(storage_size(cmplx(0.0d0,0.0d0,kind=8)),file_storage_size)/=0) then
      message='SAWF scratch stream storage-unit assumptions are invalid'
      return
    end if
    element_units=int(storage_size(cmplx(0.0d0,0.0d0,kind=8))/file_storage_size,int64)
    if(int(matrix_order,int64)>huge(0_int64)/int(matrix_order,int64)) then
      message='SAWF scratch matrix element count overflows int64'; return
    end if
    matrix_units=int(matrix_order,int64)*int(matrix_order,int64)
    if(matrix_units>huge(0_int64)/element_units) then
      message='SAWF scratch matrix storage size overflows int64'; return
    end if
    matrix_units=matrix_units*element_units
    if(int(operation_index-1,int64)>huge(0_int64)/matrix_units) then
      message='SAWF scratch operation offset overflows int64'; return
    end if
    offset=int(operation_index-1,int64)*matrix_units
    if(offset==huge(0_int64)) then
      message='SAWF scratch stream position overflows int64'; return
    end if
    position=1_int64+offset; ok=.true.
  end subroutine checked_stream_position


  subroutine abort_sawf_dmn(writer)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    integer :: io_status
    if(writer%text_unit/=0) close(writer%text_unit,iostat=io_status)
    if(writer%wann_unit/=0) close(writer%wann_unit,iostat=io_status)
    if(writer%band_unit/=0) close(writer%band_unit,iostat=io_status)
    writer%text_unit=0; writer%wann_unit=0; writer%band_unit=0
    if(allocated(writer%text_path) .and. writer%owns_text) call remove_path(writer%text_path)
    if(allocated(writer%wann_scratch_path) .and. writer%owns_wann_scratch) &
      call remove_path(writer%wann_scratch_path)
    if(allocated(writer%band_scratch_path) .and. writer%owns_band_scratch) &
      call remove_path(writer%band_scratch_path)
    writer%owns_text=.false.; writer%owns_wann_scratch=.false.; writer%owns_band_scratch=.false.
    writer%active=.false.; writer%appended=0
  end subroutine abort_sawf_dmn

  subroutine release_owned_temps(writer)
    type(t_sawf_dmn_writer), intent(inout) :: writer
    integer :: io_status
    if(writer%text_unit/=0) close(writer%text_unit,iostat=io_status)
    if(writer%wann_unit/=0) close(writer%wann_unit,iostat=io_status)
    if(writer%band_unit/=0) close(writer%band_unit,iostat=io_status)
    writer%text_unit=0; writer%wann_unit=0; writer%band_unit=0
    if(writer%owns_text) call remove_path(writer%text_path)
    if(writer%owns_wann_scratch) call remove_path(writer%wann_scratch_path)
    if(writer%owns_band_scratch) call remove_path(writer%band_scratch_path)
    writer%owns_text=.false.; writer%owns_wann_scratch=.false.; writer%owns_band_scratch=.false.
  end subroutine release_owned_temps


  real(8) function identity_matrix_residual(matrix) result(residual)
    complex(8), intent(in) :: matrix(:,:)
    integer :: i
    complex(8), allocatable :: work(:,:)
    allocate(work(size(matrix,1),size(matrix,2))); work=matrix
    do i=1,min(size(work,1),size(work,2)); work(i,i)=work(i,i)-(1.0d0,0.0d0); enddo
    residual=maxval(abs(work)); deallocate(work)
  end function identity_matrix_residual

  logical function complex_finite(matrix) result(finite)
    complex(8), intent(in) :: matrix(:,:)
    finite=all(ieee_is_finite(real(matrix,8))) .and. all(ieee_is_finite(aimag(matrix)))
  end function complex_finite

  integer function rename_path(old_path,new_path) result(status)
    character(*), intent(in) :: old_path,new_path
    character(c_char), allocatable :: cold(:),cnew(:)
    call make_c_string(old_path,cold); call make_c_string(new_path,cnew)
    status=int(c_rename(cold,cnew)); deallocate(cold,cnew)
  end function rename_path

  subroutine remove_path(path)
    character(*), intent(in) :: path
    character(c_char), allocatable :: cpath(:)
    integer(c_int) :: ignored
    call make_c_string(path,cpath); ignored=c_remove(cpath); deallocate(cpath)
  end subroutine remove_path

  subroutine make_c_string(text,c_text)
    character(*), intent(in) :: text
    character(c_char), allocatable, intent(out) :: c_text(:)
    integer :: i,n
    n=len_trim(text); allocate(c_text(n+1))
    do i=1,n; c_text(i)=text(i:i); enddo
    c_text(n+1)=c_null_char
  end subroutine make_c_string

end module lcfo_wannier_sawf_dmn
