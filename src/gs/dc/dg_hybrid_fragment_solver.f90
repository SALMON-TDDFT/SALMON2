module dg_hybrid_fragment_solver
  use mpi
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite,ieee_get_halting_mode,ieee_set_halting_mode,ieee_set_flag,&
    ieee_invalid,ieee_divide_by_zero,ieee_overflow
  use dg_hybrid_fragment_basis,only:s_dg_hybrid_fragment_basis
  implicit none
  private
  abstract interface
    subroutine fragment_apply(input,output,ok)
      import real64
      complex(real64),intent(in)::input(:,:);complex(real64),intent(out)::output(:,:);logical,intent(out)::ok
    end subroutine fragment_apply
  end interface
  public::solve_dg_hybrid_fragment_basis
contains
  subroutine solve_dg_hybrid_fragment_basis(comm,basis,nstate,occupations,core_mask,apply_h,apply_s,tolerance,&
      coefficients,eigenvalues,core_density,core_electron_count,maximum_residual,orthogonality_defect,&
      workspace_peak_bytes,fingerprint,ok,message)
    integer,intent(in)::comm,nstate
    type(s_dg_hybrid_fragment_basis),intent(in)::basis
    real(real64),intent(in)::occupations(:),tolerance
    logical,intent(in)::core_mask(:)
    procedure(fragment_apply)::apply_h,apply_s
    complex(real64),allocatable,intent(out)::coefficients(:,:)
    real(real64),intent(out)::eigenvalues(:),core_density(:),core_electron_count,maximum_residual,orthogonality_defect
    integer(int64),intent(out)::workspace_peak_bytes,fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer::i,j,k,nowned,nbasis,npoint,ierr,local_bad,global_bad,info,lwork,first_id,last_id,position
    integer,allocatable::ownership(:)
    complex(real64),allocatable::vectors(:,:),hvectors(:,:),svectors(:,:),hmat(:,:),smat(:,:),hcopy(:,:),scopy(:,:),&
      wavefunctions(:,:),work(:)
    real(real64),allocatable::all_eigenvalues(:),rwork(:)
    complex(real64)::query(1),value
    real(real64)::scale
    integer(int64)::bits
    logical::callback_ok,halt_invalid,halt_zero,halt_overflow,halting_disabled
    external::zhegv

    ok=.false.;message='';workspace_peak_bytes=0_int64;fingerprint=0_int64;halting_disabled=.false.
    maximum_residual=huge(1d0);orthogonality_defect=huge(1d0);core_electron_count=0d0
    if(.not.allocated(basis%global_ids).or..not.allocated(basis%buffer_values))then
      message='fragment solver basis is not allocated';return
    endif
    nowned=size(basis%global_ids);npoint=size(basis%buffer_values,1)
    call MPI_Allreduce(nowned,nbasis,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.nbasis<1.or.nstate<1.or.nstate>nbasis.or.size(occupations)/=nstate.or.&
        size(eigenvalues)/=nstate.or.size(core_mask)/=npoint.or.size(core_density)/=npoint)then
      message='invalid fragment eigensystem shape';return
    endif
    first_id=huge(0);last_id=-huge(0)
    if(nowned>0)then;first_id=int(minval(basis%global_ids));last_id=int(maxval(basis%global_ids));endif
    call MPI_Allreduce(MPI_IN_PLACE,first_id,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    call MPI_Allreduce(MPI_IN_PLACE,last_id,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    local_bad=0
    if(last_id-first_id+1/=nbasis.or.any(occupations<0d0).or..not.ieee_is_finite(tolerance).or.&
        tolerance<1d-15.or.tolerance>1d-2)local_bad=1
    allocate(ownership(nbasis));ownership=0
    do i=1,nowned;position=int(basis%global_ids(i))-first_id+1;if(position<1.or.position>nbasis)then
      local_bad=1;else;ownership(position)=ownership(position)+1;endif;enddo
    call MPI_Allreduce(MPI_IN_PLACE,ownership,nbasis,MPI_INTEGER,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.any(ownership/=1))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='invalid fragment basis ownership';return;endif
    allocate(vectors(npoint,nbasis),source=(0d0,0d0))
    do i=1,nowned;position=int(basis%global_ids(i))-first_id+1;vectors(:,position)=basis%buffer_values(:,i);enddo
    call MPI_Allreduce(MPI_IN_PLACE,vectors,size(vectors),MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='fragment basis column collection failed';return;endif
    allocate(hvectors(npoint,nbasis),svectors(npoint,nbasis))
    call apply_h(vectors,hvectors,callback_ok);call callback_agreement(callback_ok,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment Hamiltonian application failed';return;endif
    call apply_s(vectors,svectors,callback_ok);call callback_agreement(callback_ok,global_bad,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;message='fragment metric application failed';return;endif
    allocate(hmat(nbasis,nbasis),smat(nbasis,nbasis),hcopy(nbasis,nbasis),scopy(nbasis,nbasis))
    hmat=matmul(conjg(transpose(vectors)),hvectors);smat=matmul(conjg(transpose(vectors)),svectors)
    hcopy=hmat;scopy=smat;allocate(all_eigenvalues(nbasis),rwork(max(1,3*nbasis-2)))
    call ieee_get_halting_mode(ieee_invalid,halt_invalid);call ieee_get_halting_mode(ieee_divide_by_zero,halt_zero)
    call ieee_get_halting_mode(ieee_overflow,halt_overflow)
    call ieee_set_halting_mode(ieee_invalid,.false.);call ieee_set_halting_mode(ieee_divide_by_zero,.false.)
    call ieee_set_halting_mode(ieee_overflow,.false.);halting_disabled=.true.
    lwork=-1;call zhegv(1,'V','U',nbasis,hmat,nbasis,smat,nbasis,all_eigenvalues,query,lwork,rwork,info)
    if(info/=0)then;call restore_halting();message='fragment ZHEGV workspace query failed';return;endif
    lwork=max(1,int(real(query(1))));allocate(work(lwork));hmat=hcopy;smat=scopy
    call zhegv(1,'V','U',nbasis,hmat,nbasis,smat,nbasis,all_eigenvalues,work,lwork,rwork,info)
    call restore_halting()
    if(info/=0)then;message='fragment generalized eigensolve failed';return;endif
    allocate(coefficients(nowned,nstate),wavefunctions(npoint,nstate))
    do i=1,nowned;position=int(basis%global_ids(i))-first_id+1;coefficients(i,:)=hmat(position,1:nstate);enddo
    eigenvalues=all_eigenvalues(1:nstate);wavefunctions=matmul(vectors,hmat(:,1:nstate));core_density=0d0
    do j=1,nstate;core_density=core_density+occupations(j)*abs(wavefunctions(:,j))**2;enddo
    where(.not.core_mask)core_density=0d0
    core_electron_count=sum(core_density)
    maximum_residual=0d0;orthogonality_defect=0d0
    do j=1,nstate
      maximum_residual=max(maximum_residual,maxval(abs(matmul(hcopy,hmat(:,j))-&
        eigenvalues(j)*matmul(scopy,hmat(:,j)))))
      do k=1,nstate
        value=dot_product(hmat(:,j),matmul(scopy,hmat(:,k)))
        if(j==k)value=value-(1d0,0d0)
        orthogonality_defect=max(orthogonality_defect,abs(value))
      enddo
    enddo
    scale=max(1d0,maxval(abs(eigenvalues)))
    if(maximum_residual>100d0*tolerance*scale.or.orthogonality_defect>100d0*tolerance)then
      message='fragment generalized eigensystem receipt failed';return
    endif
    workspace_peak_bytes=16_int64*int(3*npoint*nbasis+4*nbasis*nbasis+npoint*nstate+lwork,int64)+&
      8_int64*int(nbasis+size(rwork)+npoint,int64)+4_int64*int(nbasis,int64)
    fingerprint=ieor(basis%provenance_fingerprint,int(z'9E3779B97F4A7C15',int64))
    do j=1,nstate
      bits=transfer(eigenvalues(j),bits);fingerprint=ieor(ishftc(fingerprint,11),bits)
    enddo
    if(fingerprint==0_int64)fingerprint=1223_int64
    ok=.true.;message=''
  contains
    subroutine callback_agreement(local_ok,bad,status)
      logical,intent(in)::local_ok;integer,intent(out)::bad,status;integer::local_value
      local_value=merge(0,1,local_ok);call MPI_Allreduce(local_value,bad,1,MPI_INTEGER,MPI_MAX,comm,status)
    end subroutine callback_agreement
    subroutine restore_halting()
      if(.not.halting_disabled)return
      call ieee_set_flag(ieee_invalid,.false.);call ieee_set_flag(ieee_divide_by_zero,.false.)
      call ieee_set_flag(ieee_overflow,.false.);call ieee_set_halting_mode(ieee_invalid,halt_invalid)
      call ieee_set_halting_mode(ieee_divide_by_zero,halt_zero);call ieee_set_halting_mode(ieee_overflow,halt_overflow)
      halting_disabled=.false.
    end subroutine restore_halting
  end subroutine solve_dg_hybrid_fragment_basis
end module dg_hybrid_fragment_solver
