#include "config.h"
module dg_hybrid_reciprocal_catalog
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  public::build_dg_hybrid_reciprocal_catalog
contains
  subroutine build_dg_hybrid_reciprocal_catalog(comm,reciprocal_lattice,reciprocal_rotation,&
      cutoff,tolerance,g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::reciprocal_lattice(3,3),reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,allocatable,intent(out)::g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable,intent(out)::g_vectors(:,:)
    integer(int64),intent(out)::fingerprint
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,noperation,local_bad,global_bad,i,j,k,op,ng,nmax(3),candidate_count,target,best
    integer::minimum_integer,maximum_integer,allocation_status,label,new_label
    integer,allocatable::candidate_integer(:,:)
    integer(int64)::bits,minimum_bits,maximum_bits
    real(real64)::inverse_lattice(3,3),determinant,gmax,g(3),mapped(3),distance,best_distance
    ok=.false.;message='';fingerprint=0_int64;noperation=size(reciprocal_rotation,3)
    call agree_integer(noperation,minimum_integer,maximum_integer,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='inconsistent reciprocal operation count';return
    endif
    bits=transfer(cutoff,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent reciprocal cutoff';return
    endif
    bits=transfer(tolerance,bits);call agree_int64(bits,minimum_bits,maximum_bits,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='inconsistent reciprocal tolerance';return
    endif
    local_bad=merge(0,1,noperation>0.and.cutoff>=0d0.and.tolerance>=1d-15.and.tolerance<=1d-2.and.&
      all(ieee_is_finite(reciprocal_lattice)).and.all(ieee_is_finite(reciprocal_rotation)).and.&
      ieee_is_finite(cutoff).and.ieee_is_finite(tolerance))
    call agree_real_array(reciprocal_lattice,comm,local_bad,ierr)
    call agree_real_array(reciprocal_rotation,comm,local_bad,ierr)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid reciprocal catalog input';return
    endif
    call invert_3x3(reciprocal_lattice,inverse_lattice,determinant)
    local_bad=merge(0,1,ieee_is_finite(determinant).and.abs(determinant)>tolerance)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='singular reciprocal lattice';return
    endif
    gmax=sqrt(2d0*cutoff)
    do i=1,3
      nmax(i)=ceiling(gmax*sqrt(sum(inverse_lattice(i,:)**2))+10d0*tolerance)
    enddo
    candidate_count=0
    do i=-nmax(1),nmax(1);do j=-nmax(2),nmax(2);do k=-nmax(3),nmax(3)
      g=matmul(reciprocal_lattice,real([i,j,k],real64))
      if(0.5d0*sum(g*g)<=cutoff+tolerance)candidate_count=candidate_count+1
    enddo;enddo;enddo
    if(candidate_count<1)then;message='empty reciprocal catalog';return;endif
    allocate(candidate_integer(3,candidate_count),g_integer(3,candidate_count),&
      g_vectors(3,candidate_count),g_action(candidate_count,noperation),g_star(candidate_count),&
      g_conjugate(candidate_count),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate reciprocal catalog';return
    endif
    ng=0
    do i=-nmax(1),nmax(1);do j=-nmax(2),nmax(2);do k=-nmax(3),nmax(3)
      g=matmul(reciprocal_lattice,real([i,j,k],real64))
      if(0.5d0*sum(g*g)>cutoff+tolerance)cycle
      ng=ng+1;candidate_integer(:,ng)=[i,j,k];g_vectors(:,ng)=g
    enddo;enddo;enddo
    g_integer=candidate_integer
    local_bad=0
    do i=1,ng
      call find_mode(-g_vectors(:,i),g_vectors,tolerance,best,best_distance)
      if(best==0)then;local_bad=1;else;g_conjugate(i)=best;endif
      do op=1,noperation
        mapped=matmul(reciprocal_rotation(:,:,op),g_vectors(:,i))
        call find_mode(mapped,g_vectors,100d0*tolerance,best,best_distance)
        if(best==0)then;local_bad=1;else;g_action(i,op)=best;endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='reciprocal cutoff is not closed under symmetry and conjugation';return
    endif
    do i=1,ng;g_star(i)=i;enddo
    do
      local_bad=0
      do i=1,ng
        label=min(g_star(i),g_star(g_conjugate(i)))
        do op=1,noperation;label=min(label,g_star(g_action(i,op)));enddo
        if(label/=g_star(i))then;g_star(i)=label;local_bad=1;endif
      enddo
      do i=ng,1,-1
        new_label=g_star(i)
        do j=1,ng
          if(g_star(j)==i)new_label=min(new_label,g_star(j))
        enddo
        if(new_label/=g_star(i))then;g_star(i)=new_label;local_bad=1;endif
      enddo
      if(local_bad==0)exit
    enddo
    label=0
    do i=1,ng
      if(g_star(i)==i)then
        label=label+1
        do j=1,ng;if(g_star(j)==i)g_star(j)=label;enddo
      endif
    enddo
    fingerprint=int(z'6A09E667F3BCC909',int64)
    do j=1,3;do i=1,3
      fingerprint=ieor(ishftc(fingerprint,7),transfer(reciprocal_lattice(i,j),bits))
    enddo;enddo
    fingerprint=ieor(ishftc(fingerprint,7),transfer(cutoff,bits))
    do i=1,ng
      do j=1,3;fingerprint=ieor(ishftc(fingerprint,7),int(g_integer(j,i),int64));enddo
      fingerprint=ieor(ishftc(fingerprint,7),int(g_star(i),int64))
      fingerprint=ieor(ishftc(fingerprint,7),int(g_conjugate(i),int64))
      do op=1,noperation
        fingerprint=ieor(ishftc(fingerprint,7),int(g_action(i,op),int64))
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(candidate_integer);ok=.true.
  contains
    subroutine cleanup()
      if(allocated(candidate_integer))deallocate(candidate_integer)
      if(allocated(g_integer))deallocate(g_integer)
      if(allocated(g_vectors))deallocate(g_vectors)
      if(allocated(g_action))deallocate(g_action)
      if(allocated(g_star))deallocate(g_star)
      if(allocated(g_conjugate))deallocate(g_conjugate)
    end subroutine cleanup
#else
    ok=.false.;message='reciprocal catalog requires MPI';fingerprint=0_int64
#endif
  end subroutine build_dg_hybrid_reciprocal_catalog

#ifdef USE_MPI
  subroutine agree_integer(value,minimum,maximum,comm,ierr)
    integer,intent(in)::value,comm
    integer,intent(out)::minimum,maximum,ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER,MPI_MAX,comm,ierr)
  end subroutine agree_integer

  subroutine agree_int64(value,minimum,maximum,comm,ierr)
    integer(int64),intent(in)::value
    integer(int64),intent(out)::minimum,maximum
    integer,intent(in)::comm
    integer,intent(out)::ierr
    call MPI_Allreduce(value,minimum,1,MPI_INTEGER8,MPI_MIN,comm,ierr);if(ierr/=MPI_SUCCESS)return
    call MPI_Allreduce(value,maximum,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
  end subroutine agree_int64

  subroutine agree_real_array(values,comm,bad,ierr)
    real(real64),intent(in)::values(..)
    integer,intent(in)::comm
    integer,intent(inout)::bad
    integer,intent(out)::ierr
    integer::i,j,k
    integer(int64)::bits,minimum,maximum
    select rank(values)
    rank(2)
      do j=1,size(values,2);do i=1,size(values,1)
        bits=transfer(values(i,j),bits);call agree_int64(bits,minimum,maximum,comm,ierr)
        if(ierr/=MPI_SUCCESS)return
        if(minimum/=maximum)bad=1
      enddo;enddo
    rank(3)
      do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
        bits=transfer(values(i,j,k),bits);call agree_int64(bits,minimum,maximum,comm,ierr)
        if(ierr/=MPI_SUCCESS)return
        if(minimum/=maximum)bad=1
      enddo;enddo;enddo
    end select
  end subroutine agree_real_array

  subroutine invert_3x3(matrix,inverse,determinant)
    real(real64),intent(in)::matrix(3,3)
    real(real64),intent(out)::inverse(3,3),determinant
    determinant=matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))-&
      matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+&
      matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
    inverse(1,:)=[matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2),&
      matrix(1,3)*matrix(3,2)-matrix(1,2)*matrix(3,3),matrix(1,2)*matrix(2,3)-matrix(1,3)*matrix(2,2)]
    inverse(2,:)=[matrix(2,3)*matrix(3,1)-matrix(2,1)*matrix(3,3),&
      matrix(1,1)*matrix(3,3)-matrix(1,3)*matrix(3,1),matrix(1,3)*matrix(2,1)-matrix(1,1)*matrix(2,3)]
    inverse(3,:)=[matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1),&
      matrix(1,2)*matrix(3,1)-matrix(1,1)*matrix(3,2),matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)]
    if(determinant/=0d0)inverse=inverse/determinant
  end subroutine invert_3x3

  subroutine find_mode(vector,vectors,tolerance,index,distance)
    real(real64),intent(in)::vector(3),vectors(:,:),tolerance
    integer,intent(out)::index
    real(real64),intent(out)::distance
    integer::i
    index=0;distance=huge(1d0)
    do i=1,size(vectors,2)
      if(maxval(abs(vector-vectors(:,i)))<distance)then
        distance=maxval(abs(vector-vectors(:,i)));index=i
      endif
    enddo
    if(distance>tolerance)index=0
  end subroutine find_mode
#endif
end module dg_hybrid_reciprocal_catalog
