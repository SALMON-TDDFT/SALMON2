#include "config.h"
module dg_nonlocal_projector_range
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::s_dg_nonlocal_range_receipt,analyze_dg_nonlocal_projector_range
  type s_dg_nonlocal_range_receipt
    real(real64)::local_contribution=0d0,adjacent_contribution=0d0,remote_contribution=0d0
    real(real64)::maximum_remote_fraction=0d0,symmetry_pair_defect=huge(1d0)
    integer::unmatched_channel_count=0
    integer(int64)::workspace_peak_bytes=0_int64
  end type
contains
  subroutine analyze_dg_nonlocal_projector_range(comm,grid_ids,wannier,centers_fractional,lattice,&
      atoms_cartesian,species,projector_atom,strength,projector_offsets,projector_grid_ids,projector_values,&
      rotation,translation,fragment_shape,&
      tile_width,wannier_representation,projector_representation,receipt,ok,message,row_ids,operator_rows)
    integer,intent(in)::comm,species(:),projector_atom(:),projector_offsets(:),rotation(3,3),&
      fragment_shape(3),tile_width
    integer(int64),intent(in)::grid_ids(:),projector_grid_ids(:)
    complex(real64),intent(in)::wannier(:,:),projector_values(:),wannier_representation(:,:),&
      projector_representation(:,:)
    real(real64),intent(in)::centers_fractional(:,:),lattice(3,3),atoms_cartesian(:,:),strength(:),translation(3)
    type(s_dg_nonlocal_range_receipt),intent(out)::receipt
    logical,intent(out)::ok
    character(*),intent(out)::message
    integer(int64),intent(in),optional::row_ids(:)
    complex(real64),allocatable,intent(out),optional::operator_rows(:,:)
#ifdef USE_MPI
    integer::ierr,nw,np,na,i,p,q,g,j0,j1,width,local_bad,global_bad,fc(3),fa(3),df(3),class
    real(real64)::inverse(3,3),determinant,atom_fractional(3),value,total(3),remote_by_w(size(wannier,1)),&
      total_by_w(size(wannier,1)),local_scale,global_scale,local_defect,global_defect
    complex(real64),allocatable::local_overlap(:,:),overlap(:,:),transformed(:,:),gram(:,:)
    receipt=s_dg_nonlocal_range_receipt();ok=.false.;message='';nw=size(wannier,1);np=size(projector_atom);na=size(species)
    local_bad=0
    if(nw<1.or.np<1.or.na<1.or.tile_width<1.or.size(grid_ids)/=size(wannier,2).or.&
      any(shape(centers_fractional)/=[3,nw]).or.&
      any(shape(atoms_cartesian)/=[3,na]).or.size(projector_atom)/=np.or.size(strength)/=np.or.&
      size(projector_offsets)/=np+1.or.size(projector_grid_ids)/=size(projector_values).or.&
      any(shape(wannier_representation)/=[nw,nw]).or.any(shape(projector_representation)/=[np,np]).or.&
      any(projector_atom<1).or.any(projector_atom>na).or.any(fragment_shape<1))local_bad=1
    if(present(row_ids).neqv.present(operator_rows))local_bad=1
    if(present(row_ids))then
      if(any(row_ids<1_int64).or.any(row_ids>int(nw,int64)))local_bad=1
    endif
    if(local_bad==0)then
      if(projector_offsets(1)/=1.or.projector_offsets(np+1)/=size(projector_values)+1.or.&
        any(projector_offsets(2:np+1)<projector_offsets(1:np)))local_bad=1
    endif
    if(.not.all(ieee_is_finite(real(wannier))).or..not.all(ieee_is_finite(aimag(wannier))).or.&
      .not.all(ieee_is_finite(real(projector_values))).or..not.all(ieee_is_finite(aimag(projector_values))).or.&
      .not.all(ieee_is_finite(strength)).or.&
      .not.all(ieee_is_finite(real(wannier_representation))).or.&
      .not.all(ieee_is_finite(aimag(wannier_representation))).or.&
      .not.all(ieee_is_finite(real(projector_representation))).or.&
      .not.all(ieee_is_finite(aimag(projector_representation))))local_bad=1
    if(local_bad==0)then
      allocate(gram(max(nw,np),max(nw,np)));gram=0d0
      gram(1:nw,1:nw)=matmul(conjg(transpose(wannier_representation)),wannier_representation)
      do i=1,nw;gram(i,i)=gram(i,i)-1d0;enddo
      if(maxval(abs(gram(1:nw,1:nw)))>1d-10)local_bad=1
      gram=0d0
      gram(1:np,1:np)=matmul(conjg(transpose(projector_representation)),projector_representation)
      do i=1,np;gram(i,i)=gram(i,i)-1d0;enddo
      if(maxval(abs(gram(1:np,1:np)))>1d-10)local_bad=1
      deallocate(gram)
    endif
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='invalid nonlocal projector range contract';return;endif
    call invert3(lattice,inverse,determinant)
    if(abs(determinant)<1d-14)then;message='singular nonlocal diagnostic lattice';return;endif
    allocate(overlap(np,nw));overlap=0d0
    do j0=1,nw,tile_width
      j1=min(nw,j0+tile_width-1);width=j1-j0+1
      allocate(local_overlap(np,width));local_overlap=0d0
      do p=1,np
        do q=projector_offsets(p),projector_offsets(p+1)-1
          g=find_grid_position(grid_ids,projector_grid_ids(q))
          if(g<1)then
            local_bad=1;cycle
          endif
          do i=1,width
            local_overlap(p,i)=local_overlap(p,i)+conjg(projector_values(q))*wannier(j0+i-1,g)
          enddo
        enddo
      enddo
      call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
      if(global_bad/=0.or.ierr/=MPI_SUCCESS)then;message='projector support is not present on its owning grid rank';return;endif
      call MPI_Allreduce(local_overlap,overlap(:,j0:j1),np*width,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,ierr)
      deallocate(local_overlap)
      if(ierr/=MPI_SUCCESS)then;message='nonlocal overlap reduction failed';return;endif
    enddo
    total=0d0;remote_by_w=0d0;total_by_w=0d0
    do i=1,nw
      fc=floor(modulo(centers_fractional(:,i),1d0)*real(fragment_shape,real64))
      do p=1,np
        atom_fractional=modulo(matmul(inverse,atoms_cartesian(:,projector_atom(p))),1d0)
        fa=floor(atom_fractional*real(fragment_shape,real64));df=abs(fc-fa);df=min(df,fragment_shape-df)
        if(maxval(df)==0)then;class=1;else if(maxval(df)<=1)then;class=2;else;class=3;endif
        value=abs(strength(p))*abs(overlap(p,i))**2;total(class)=total(class)+value
        total_by_w(i)=total_by_w(i)+value;if(class==3)remote_by_w(i)=remote_by_w(i)+value
      enddo
    enddo
    receipt%local_contribution=total(1);receipt%adjacent_contribution=total(2);receipt%remote_contribution=total(3)
    do i=1,nw
      if(total_by_w(i)>tiny(1d0))receipt%maximum_remote_fraction=max(receipt%maximum_remote_fraction,&
        remote_by_w(i)/total_by_w(i))
    enddo
    allocate(transformed(np,nw));transformed=matmul(conjg(transpose(projector_representation)),&
      matmul(overlap,wannier_representation))
    local_defect=maxval(abs(transformed-overlap));local_scale=max(1d0,maxval(abs(overlap)))
    call MPI_Allreduce(local_defect,global_defect,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call MPI_Allreduce(local_scale,global_scale,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    receipt%symmetry_pair_defect=global_defect/global_scale
    if(present(operator_rows))then
      allocate(operator_rows(size(row_ids),nw));operator_rows=0d0
      do i=1,size(row_ids)
        do p=1,np
          operator_rows(i,:)=operator_rows(i,:)+strength(p)*conjg(overlap(p,int(row_ids(i))))*overlap(p,:)
        enddo
      enddo
    endif
    receipt%workspace_peak_bytes=16_int64*int(2*np*min(tile_width,nw)+np*nw,int64)
    if(present(operator_rows))receipt%workspace_peak_bytes=receipt%workspace_peak_bytes+&
      16_int64*int(size(operator_rows),int64)
    ok=ierr==MPI_SUCCESS.and.ieee_is_finite(receipt%symmetry_pair_defect)
    if(.not.ok)message='nonlocal projector range result is not finite'
#else
    receipt=s_dg_nonlocal_range_receipt();ok=.false.;message='nonlocal projector range diagnostic requires MPI'
#endif
  contains
    integer function find_grid_position(ids,target) result(position)
      integer(int64),intent(in)::ids(:),target
      integer::lo,hi,mid
      lo=1;hi=size(ids);position=0
      do while(lo<=hi)
        mid=lo+(hi-lo)/2
        if(ids(mid)==target)then;position=mid;return
        else if(ids(mid)<target)then;lo=mid+1
        else;hi=mid-1
        endif
      enddo
    end function

    subroutine invert3(a,b,d)
      real(real64),intent(in)::a(3,3);real(real64),intent(out)::b(3,3),d
      d=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+&
        a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1));b=0d0;if(abs(d)<1d-14)return
      b(1,:)=[a(2,2)*a(3,3)-a(2,3)*a(3,2),a(1,3)*a(3,2)-a(1,2)*a(3,3),a(1,2)*a(2,3)-a(1,3)*a(2,2)]/d
      b(2,:)=[a(2,3)*a(3,1)-a(2,1)*a(3,3),a(1,1)*a(3,3)-a(1,3)*a(3,1),a(1,3)*a(2,1)-a(1,1)*a(2,3)]/d
      b(3,:)=[a(2,1)*a(3,2)-a(2,2)*a(3,1),a(1,2)*a(3,1)-a(1,1)*a(3,2),a(1,1)*a(2,2)-a(1,2)*a(2,1)]/d
    end subroutine
  end subroutine
end module
