module dg_hybrid_production_fragment_basis
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis,build_dg_hybrid_fragment_basis
  use dg_hybrid_wannier_complement,only:project_dg_hybrid_wannier_complement
  implicit none
  private
  public::build_dg_hybrid_production_fragment_basis
contains
  subroutine build_dg_hybrid_production_fragment_basis(comm,global_point_count,point_ids,weights,&
      wannier_values,pw_values,packet_ids,near_offsets,near_wannier_ids,wannier_fingerprint,&
      packet_fingerprint,global_basis_offset,fragment_id,tolerance,basis,workspace_peak_bytes,&
      fingerprint,ok,message,physical_point_ids)
    integer,intent(in)::comm,global_point_count,fragment_id
    integer(int64),intent(in)::point_ids(:),wannier_fingerprint,packet_fingerprint,global_basis_offset
    real(real64),intent(in)::weights(:),tolerance
    complex(real64),intent(in)::wannier_values(:,:),pw_values(:,:)
    integer,intent(in)::packet_ids(:),near_offsets(:),near_wannier_ids(:)
    type(s_dg_hybrid_fragment_basis),intent(out)::basis
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(in),optional::physical_point_ids(:)
    complex(real64),allocatable::projected_pw(:,:),local_full(:,:),global_full(:,:),wf_owned(:,:),pw_owned(:,:)
    integer(int64),allocatable::wf_ids(:),pw_ids(:)
    integer(int64),allocatable::local_physical_ids(:),global_physical_ids(:)
    integer::rank,nproc,ierr,nw,np,nbasis,nwf_owned,npw_owned,i,p,iw,ip
    integer(int64)::complement_workspace,complement_fingerprint,transpose_bytes
    real(real64)::omitted_tail
    logical::project_ok
    character(256)::project_message

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)then;message='fragment basis rank failed';return;endif
    call MPI_Comm_size(comm,nproc,ierr);if(ierr/=MPI_SUCCESS)then;message='fragment basis size failed';return;endif
    nw=size(wannier_values,1);np=size(pw_values,1);nbasis=nw+np
    if(present(physical_point_ids))then
      if(size(physical_point_ids)/=size(point_ids))then;message='fragment physical point ID shape mismatch';return;endif
    endif
    call project_dg_hybrid_wannier_complement(comm,global_point_count,point_ids,weights,wannier_values,pw_values,&
      wannier_fingerprint,packet_fingerprint,packet_ids,near_offsets,near_wannier_ids,.false.,tolerance,&
      projected_pw,omitted_tail,complement_workspace,complement_fingerprint,project_ok,project_message)
    if(.not.project_ok)then;message=trim(project_message);return;endif
    allocate(local_full(global_point_count,nbasis),global_full(global_point_count,nbasis),source=(0d0,0d0))
    local_full=(0d0,0d0)
    do i=1,size(point_ids)
      local_full(int(point_ids(i)),1:nw)=wannier_values(:,i)
      local_full(int(point_ids(i)),nw+1:nbasis)=projected_pw(:,i)
    enddo
    call MPI_Allreduce(local_full,global_full,size(global_full),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment spatial-to-column redistribution failed';return;endif
    allocate(local_physical_ids(global_point_count),global_physical_ids(global_point_count));local_physical_ids=0_int64
    do i=1,size(point_ids)
      if(present(physical_point_ids))then
        local_physical_ids(int(point_ids(i)))=physical_point_ids(i)
      else
        local_physical_ids(int(point_ids(i)))=point_ids(i)
      endif
    enddo
    call MPI_Allreduce(local_physical_ids,global_physical_ids,global_point_count,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(global_physical_ids<=0_int64))then
      message='fragment physical point ID redistribution failed';return
    endif
    nwf_owned=count([(mod(p-1,nproc)==rank,p=1,nw)])
    npw_owned=count([(mod(nw+p-1,nproc)==rank,p=1,np)])
    allocate(wf_ids(nwf_owned),pw_ids(npw_owned),wf_owned(global_point_count,nwf_owned),&
      pw_owned(global_point_count,npw_owned));iw=0;ip=0
    do p=1,nw
      if(mod(p-1,nproc)/=rank)cycle
      iw=iw+1;wf_ids(iw)=global_basis_offset+int(p,int64);wf_owned(:,iw)=global_full(:,p)
    enddo
    do p=1,np
      if(mod(nw+p-1,nproc)/=rank)cycle
      ip=ip+1;pw_ids(ip)=global_basis_offset+int(nw+p,int64);pw_owned(:,ip)=global_full(:,nw+p)
    enddo
    call build_dg_hybrid_fragment_basis(comm,fragment_id,wf_ids,wf_owned,pw_ids,pw_owned,0,0,basis,ok,message)
    if(.not.ok)return
    basis%buffer_point_ids=global_physical_ids
    transpose_bytes=32_int64*int(global_point_count,int64)*int(nbasis,int64)
    workspace_peak_bytes=max(complement_workspace,transpose_bytes)
    fingerprint=ieor(basis%provenance_fingerprint,ishftc(complement_fingerprint,17))
    if(fingerprint==0_int64)fingerprint=811_int64
  end subroutine build_dg_hybrid_production_fragment_basis
end module dg_hybrid_production_fragment_basis
