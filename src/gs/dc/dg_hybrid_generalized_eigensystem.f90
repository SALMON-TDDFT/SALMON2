#include "config.h"
module dg_hybrid_generalized_eigensystem
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero,ieee_overflow
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::solve_dg_hybrid_generalized_scalapack
contains
  subroutine solve_dg_hybrid_generalized_scalapack(comm,global_count,nstate,row_ids,hrows,srows,tolerance,&
      coefficients,eigenvalues,maximum_residual,orthogonality_defect,projector_defect,workspace_peak_bytes,&
      fingerprint,ok,message)
    integer,intent(in)::comm,global_count,nstate
    integer(int64),intent(in)::row_ids(:)
    complex(real64),intent(in)::hrows(:,:),srows(:,:)
    real(real64),intent(in)::tolerance
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::eigenvalues(:),maximum_residual,orthogonality_defect,projector_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#if defined(USE_MPI) && defined(USE_SCALAPACK)
    integer::rank,nproc,ierr,nowned,local_bad,global_bad,minimum_integer,maximum_integer
    integer::nprow,npcol,myrow,mycol,ictxt,mb,nb,local_rows,local_columns,allocation_status
    integer::i,j,k,row,owner_rank,owner_position,local_i,local_j,process_row,process_column,info
    integer::group_comm,group_world,minimum_workspace,maximum_workspace,lwork,lrwork,liwork
    integer,allocatable::ownership_count(:),owner(:),position(:),gridmap(:,:),comm_ranks(:),world_ranks(:),iwork(:)
    integer::desca(9),descb(9),descz(9)
    integer(int64)::bits,minimum_bits,maximum_bits,entry_hash,complex_elements,integer_elements,quantized,extra_bytes
    complex(real64),allocatable::adiv(:,:),bdiv(:,:),zdiv(:,:),work(:),remote_row(:),stream(:),hc(:,:),sc(:,:)
    real(real64),allocatable::rwork(:),all_eigenvalues(:)
    complex(real64)::local_dot,global_dot
    real(real64)::scale,local_value,global_value,quantization_scale,quantization_limit,hermitian_defect,metric_defect,&
      matrix_scale
    logical::halt_invalid,halt_zero,halt_overflow,halting_disabled
    integer,external::NUMROC
    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;halting_disabled=.false.
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);projector_defect=huge(1d0)
    eigenvalues=0d0;nowned=size(row_ids);local_bad=0
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='ScaLAPACK communicator rank failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='ScaLAPACK communicator size failed';return;endif
    call agree_integer(global_count,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing generalized extent';return;endif
    call agree_integer(nstate,minimum_integer,maximum_integer,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then;message='rank-disagreeing generalized state count';return;endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then;message='rank-disagreeing generalized tolerance';return;endif
    if(global_count<1.or.nstate<1.or.nstate>global_count)local_bad=1
    if(size(hrows,1)/=nowned.or.size(hrows,2)/=global_count.or.any(shape(srows)/=shape(hrows)))local_bad=1
    if(size(eigenvalues)/=nstate.or.any(row_ids<1_int64).or.any(row_ids>int(max(0,global_count),int64)))local_bad=1
    if(.not.ieee_is_finite(tolerance).or.tolerance<1d-15.or.tolerance>1d-2)local_bad=1
    if(.not.finite_matrix(hrows).or..not.finite_matrix(srows))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid complex generalized eigensystem contract';return;endif
    nprow=1
    do i=1,int(sqrt(real(nproc,real64)));if(mod(nproc,i)==0)nprow=i;enddo
    npcol=nproc/nprow;mb=1;nb=1
    local_rows=NUMROC(global_count,mb,mod(rank,nprow),0,nprow)
    local_columns=NUMROC(global_count,nb,rank/nprow,0,npcol)
    complex_elements=3_int64*int(max(1,local_rows),int64)*int(max(1,local_columns),int64)+&
      int(global_count,int64)+int(nstate,int64)+3_int64*int(nowned,int64)*int(nstate,int64)
    integer_elements=3_int64*int(global_count,int64)+3_int64*int(nproc,int64)
    if(complex_elements>huge(workspace_peak_bytes)/16_int64.or.integer_elements>huge(workspace_peak_bytes)/4_int64)local_bad=1
    if(local_bad==0.and.16_int64*complex_elements>huge(workspace_peak_bytes)-4_int64*integer_elements)local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='generalized ScaLAPACK base receipt overflow';return;endif
    workspace_peak_bytes=16_int64*complex_elements+4_int64*integer_elements
    allocate(ownership_count(global_count),owner(global_count),position(global_count),gridmap(nprow,npcol),&
      comm_ranks(nproc),world_ranks(nproc),adiv(max(1,local_rows),max(1,local_columns)),&
      bdiv(max(1,local_rows),max(1,local_columns)),zdiv(max(1,local_rows),max(1,local_columns)),&
      remote_row(global_count),stream(nstate),hc(nowned,nstate),sc(nowned,nstate),all_eigenvalues(global_count),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup();message='cannot allocate generalized ScaLAPACK workspace';return;endif
    ownership_count=0;owner=-1;position=0
    do i=1,nowned;row=int(row_ids(i));ownership_count(row)=ownership_count(row)+1;owner(row)=rank;position(row)=i;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership_count,global_count,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='generalized ownership count failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,owner,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup();message='generalized owner reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,position,global_count,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership_count/=1))then;call cleanup();message='generalized rows are not owned exactly once';return;endif
    do i=1,nproc;comm_ranks(i)=i-1;enddo
    call MPI_Comm_group(comm,group_comm,ierr);if(ierr/=MPI_SUCCESS)then;call cleanup();message='communicator group failed';return;endif
    call MPI_Comm_group(MPI_COMM_WORLD,group_world,ierr)
    call MPI_Group_translate_ranks(group_comm,nproc,comm_ranks,group_world,world_ranks,ierr)
    call MPI_Group_free(group_world,ierr);call MPI_Group_free(group_comm,ierr)
    if(any(world_ranks==MPI_UNDEFINED))then;call cleanup();message='BLACS rank translation failed';return;endif
    k=0
    do j=1,npcol;do i=1,nprow;k=k+1;gridmap(i,j)=world_ranks(k);enddo;enddo
    call BLACS_GET(0,0,ictxt);call BLACS_GRIDMAP(ictxt,gridmap,nprow,nprow,npcol)
    call BLACS_GRIDINFO(ictxt,nprow,npcol,myrow,mycol)
    call DESCINIT(desca,global_count,global_count,mb,nb,0,0,ictxt,max(1,local_rows),info)
    descb=desca;descz=desca
    if(info/=0)then;call cleanup_grid();message='ScaLAPACK descriptor initialization failed';return;endif
    adiv=(0d0,0d0);bdiv=(0d0,0d0);hermitian_defect=0d0;metric_defect=0d0;matrix_scale=1d0
    do row=1,global_count
      remote_row=(0d0,0d0);if(rank==owner(row))remote_row=hrows(position(row),:)
      call MPI_Bcast(remote_row,global_count,MPI_DOUBLE_COMPLEX,owner(row),comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup_grid();message='Hamiltonian row broadcast failed';return;endif
      matrix_scale=max(matrix_scale,maxval(abs(remote_row)))
      do i=1,nowned
        hermitian_defect=max(hermitian_defect,abs(hrows(i,row)-conjg(remote_row(int(row_ids(i))))))
      enddo
      do j=1,global_count
        call INFOG2L(row,j,desca,nprow,npcol,myrow,mycol,local_i,local_j,process_row,process_column)
        if(myrow==process_row.and.mycol==process_column)adiv(local_i,local_j)=remote_row(j)
      enddo
      remote_row=(0d0,0d0);if(rank==owner(row))remote_row=srows(position(row),:)
      call MPI_Bcast(remote_row,global_count,MPI_DOUBLE_COMPLEX,owner(row),comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup_grid();message='metric row broadcast failed';return;endif
      matrix_scale=max(matrix_scale,maxval(abs(remote_row)))
      do i=1,nowned
        metric_defect=max(metric_defect,abs(srows(i,row)-conjg(remote_row(int(row_ids(i))))))
      enddo
      do j=1,global_count
        call INFOG2L(row,j,descb,nprow,npcol,myrow,mycol,local_i,local_j,process_row,process_column)
        if(myrow==process_row.and.mycol==process_column)bdiv(local_i,local_j)=remote_row(j)
      enddo
    enddo
    call MPI_Allreduce(MPI_IN_PLACE,hermitian_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;call cleanup_grid();message='Hamiltonian Hermiticity reduction failed';return;endif
    call MPI_Allreduce(MPI_IN_PLACE,metric_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.max(hermitian_defect,metric_defect)>100d0*tolerance*matrix_scale)then
      call cleanup_grid();message='complex generalized pencil is not Hermitian';return
    endif
    call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
    call ieee_get_halting_mode(ieee_overflow,halt_overflow)
    call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.);halting_disabled=.true.
    call PZPOTRF('L',global_count,bdiv,1,1,descb,info)
    call MPI_Allreduce(info,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='complex hybrid metric is not positive definite';return;endif
    call PZHEGST(1,'L',global_count,adiv,1,1,desca,bdiv,1,1,descb,scale,info)
    call MPI_Allreduce(info,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0.or..not.ieee_is_finite(scale))then
      call cleanup_grid();message='complex generalized reduction failed';return
    endif
    lwork=-1;lrwork=-1;liwork=-1
    allocate(work(1),rwork(1),iwork(1));call PZHEEVD('V','L',global_count,adiv,1,1,desca,&
      all_eigenvalues,zdiv,1,1,descz,work,lwork,rwork,lrwork,iwork,liwork,info)
    local_bad=merge(0,1,info==0.and.ieee_is_finite(real(work(1))).and.ieee_is_finite(rwork(1)))
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='PZHEEVD workspace query failed';return;endif
    lwork=max(1,int(real(work(1))));lrwork=max(1,int(rwork(1)));liwork=max(1,iwork(1))
    local_bad=0;extra_bytes=0_int64
    if(int(lwork,int64)>huge(extra_bytes)/16_int64.or.int(lrwork,int64)>huge(extra_bytes)/8_int64.or.&
      int(liwork,int64)>huge(extra_bytes)/4_int64.or.int(global_count,int64)>huge(extra_bytes)/8_int64)local_bad=1
    if(local_bad==0)then
      extra_bytes=16_int64*int(lwork,int64)+8_int64*int(lrwork,int64)+4_int64*int(liwork,int64)+&
        8_int64*int(global_count,int64)
      if(workspace_peak_bytes>huge(workspace_peak_bytes)-extra_bytes)local_bad=1
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='PZHEEVD workspace receipt overflow';return;endif
    workspace_peak_bytes=workspace_peak_bytes+extra_bytes
    deallocate(work,rwork,iwork);allocate(work(lwork),rwork(lrwork),iwork(liwork),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0);call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='cannot allocate PZHEEVD workspace';return;endif
    call PZHEEVD('V','L',global_count,adiv,1,1,desca,all_eigenvalues,zdiv,1,1,descz,&
      work,lwork,rwork,lrwork,iwork,liwork,info)
    call MPI_Allreduce(info,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='complex ScaLAPACK eigensystem failed';return;endif
    call PZTRSM('L','L','C','N',global_count,global_count,(1d0,0d0),bdiv,1,1,descb,zdiv,1,1,descz)
    call restore_halting()
    allocate(coefficients(nowned,nstate),stat=allocation_status);local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='cannot allocate occupied ScaLAPACK rows';return;endif
    coefficients=(0d0,0d0)
    do row=1,global_count
      stream=(0d0,0d0)
      do j=1,nstate
        call INFOG2L(row,j,descz,nprow,npcol,myrow,mycol,local_i,local_j,process_row,process_column)
        if(myrow==process_row.and.mycol==process_column)stream(j)=zdiv(local_i,local_j)
      enddo
      call MPI_Allreduce(MPI_IN_PLACE,stream,nstate,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(ierr/=MPI_SUCCESS)then;call cleanup_grid();message='occupied row collection failed';return;endif
      if(rank==owner(row))coefficients(position(row),:)=stream
    enddo
    eigenvalues=all_eigenvalues(1:nstate)
    hc=(0d0,0d0);sc=(0d0,0d0)
    do row=1,global_count
      stream=(0d0,0d0);if(rank==owner(row))stream=coefficients(position(row),:)
      call MPI_Bcast(stream,nstate,MPI_DOUBLE_COMPLEX,owner(row),comm,ierr)
      do i=1,nowned;hc(i,:)=hc(i,:)+hrows(i,row)*stream;sc(i,:)=sc(i,:)+srows(i,row)*stream;enddo
    enddo
    local_value=0d0
    do j=1,nstate;local_value=max(local_value,maxval(abs(hc(:,j)-eigenvalues(j)*sc(:,j))));enddo
    call MPI_Allreduce(local_value,maximum_residual,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    orthogonality_defect=0d0
    do j=1,nstate;do k=1,nstate
      local_dot=sum(conjg(coefficients(:,j))*sc(:,k));call MPI_Allreduce(local_dot,global_dot,1,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      if(j==k)global_dot=global_dot-(1d0,0d0);orthogonality_defect=max(orthogonality_defect,abs(global_dot))
    enddo;enddo
    projector_defect=orthogonality_defect
    global_value=max(1d0,maxval(abs(eigenvalues)))
    local_bad=merge(0,1,maximum_residual<=100d0*tolerance*global_value.and.&
      orthogonality_defect<=100d0*tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='complex generalized eigensystem residual gate failed';return;endif
    quantization_scale=1000d0*tolerance;quantization_limit=0.25d0*real(huge(0_int64),real64)*quantization_scale
    fingerprint=1709_int64
    do i=1,global_count
      stream=(0d0,0d0);if(rank==owner(i))stream=coefficients(position(i),:)
      call MPI_Bcast(stream,nstate,MPI_DOUBLE_COMPLEX,owner(i),comm,ierr)
      do j=1,global_count
        remote_row(1:nstate)=(0d0,0d0)
        if(rank==owner(j))remote_row(1:nstate)=coefficients(position(j),:)
        call MPI_Bcast(remote_row,nstate,MPI_DOUBLE_COMPLEX,owner(j),comm,ierr)
        local_dot=sum(stream*conjg(remote_row(1:nstate)))
        if(abs(real(local_dot))>quantization_limit.or.abs(aimag(local_dot))>quantization_limit)local_bad=1
        if(local_bad==0)then
          quantized=nint(real(local_dot)/quantization_scale,int64)
          entry_hash=ieor(int(i,int64),ishftc(int(j,int64),7));entry_hash=ieor(entry_hash,ishftc(quantized,17))
          quantized=nint(aimag(local_dot)/quantization_scale,int64);entry_hash=ieor(entry_hash,ishftc(quantized,31))
          fingerprint=ieor(ishftc(fingerprint,1),entry_hash)
        endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;call cleanup_grid();message='projector fingerprint range is unsafe';return;endif
    do j=1,nstate
      if(abs(eigenvalues(j))>quantization_limit)then;call cleanup_grid();message='eigenvalue fingerprint range is unsafe';return;endif
      quantized=nint(eigenvalues(j)/quantization_scale,int64)
      fingerprint=ieor(fingerprint,ishftc(quantized,mod(11*j,63)))
    enddo
    if(fingerprint==0_int64)fingerprint=991_int64
    ok=.true.;message='';call cleanup_grid(.false.)
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
    logical function finite_matrix(values)
      complex(real64),intent(in)::values(:,:)
      finite_matrix=all(ieee_is_finite(real(values))).and.all(ieee_is_finite(aimag(values)))
    end function finite_matrix
    subroutine cleanup_grid(clear_output)
      logical,intent(in),optional::clear_output
      logical::remove_output
      remove_output=.true.;if(present(clear_output))remove_output=clear_output
      call restore_halting()
      call BLACS_GRIDEXIT(ictxt)
      if(remove_output.and.allocated(coefficients))deallocate(coefficients)
      call cleanup()
    end subroutine cleanup_grid
    subroutine restore_halting()
      if(.not.halting_disabled)return
      call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
      call ieee_set_flag(ieee_overflow,.false.);call ieee_set_halting_mode(ieee_invalid,halt_invalid)
      call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero);call ieee_set_halting_mode(ieee_overflow,halt_overflow)
      halting_disabled=.false.
    end subroutine restore_halting
    subroutine cleanup()
      if(allocated(ownership_count))deallocate(ownership_count)
      if(allocated(owner))deallocate(owner)
      if(allocated(position))deallocate(position)
      if(allocated(gridmap))deallocate(gridmap)
      if(allocated(comm_ranks))deallocate(comm_ranks)
      if(allocated(world_ranks))deallocate(world_ranks)
      if(allocated(adiv))deallocate(adiv)
      if(allocated(bdiv))deallocate(bdiv)
      if(allocated(zdiv))deallocate(zdiv)
      if(allocated(work))deallocate(work)
      if(allocated(rwork))deallocate(rwork)
      if(allocated(all_eigenvalues))deallocate(all_eigenvalues)
      if(allocated(iwork))deallocate(iwork)
      if(allocated(remote_row))deallocate(remote_row)
      if(allocated(stream))deallocate(stream)
      if(allocated(hc))deallocate(hc)
      if(allocated(sc))deallocate(sc)
    end subroutine cleanup
#else
    ok=.false.;message='complex ScaLAPACK is required for the hybrid reference eigensystem'
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);projector_defect=huge(1d0)
    workspace_peak_bytes=0_int64;fingerprint=0_int64;eigenvalues=0d0
#endif
  end subroutine solve_dg_hybrid_generalized_scalapack
end module dg_hybrid_generalized_eigensystem
