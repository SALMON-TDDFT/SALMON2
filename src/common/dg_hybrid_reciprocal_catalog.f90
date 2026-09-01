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
      cutoff,tolerance,g_integer,g_vectors,g_action,g_star,g_conjugate,fingerprint,&
      effective_cutoff,shell_added,orbit_added,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::reciprocal_lattice(3,3),reciprocal_rotation(:,:,:),cutoff,tolerance
    integer,allocatable,intent(out)::g_integer(:,:),g_action(:,:),g_star(:),g_conjugate(:)
    real(real64),allocatable,intent(out)::g_vectors(:,:)
    integer(int64),intent(out)::fingerprint
    real(real64),intent(out)::effective_cutoff
    integer,intent(out)::shell_added,orbit_added
    logical,intent(out)::ok
    character(*),intent(out)::message
#ifdef USE_MPI
    integer::ierr,noperation,local_bad,global_bad,i,j,k,op,ng,nmax(3),candidate_count,best
    integer::minimum_integer,maximum_integer,allocation_status,label,new_label,a,b,c,identity_count,&
      candidate_index,target_candidate
    integer,allocatable::candidate_integer(:,:),integer_rotation(:,:,:),selected_position(:)
    integer(int64)::bits,minimum_bits,maximum_bits
    real(real64)::inverse_lattice(3,3),determinant,gmax,g(3),distance,best_distance,&
      base_tau,search_cutoff,boundary_energy,operation_scale
    real(real64)::integer_transform(3,3),identity_real(3,3)
    real(real64),allocatable::candidate_vectors(:,:),candidate_energy(:)
    integer::identity_integer(3,3),product_integer(3,3),mapped_integer(3)
    logical,allocatable::selected(:)
    logical::changed,found
    ok=.false.;message='';fingerprint=0_int64;effective_cutoff=0d0;shell_added=0;orbit_added=0
    noperation=size(reciprocal_rotation,3)
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
    call agree_real_matrix(reciprocal_lattice,comm,local_bad,ierr)
    call agree_real_cube(reciprocal_rotation,comm,local_bad,ierr)
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

    allocate(integer_rotation(3,3,noperation),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate reciprocal operation workspace';return
    endif
    identity_integer=0;identity_real=0d0
    do i=1,3;identity_integer(i,i)=1;identity_real(i,i)=1d0;enddo
    local_bad=0;identity_count=0
    do op=1,noperation
      integer_transform=matmul(inverse_lattice,matmul(reciprocal_rotation(:,:,op),reciprocal_lattice))
      integer_rotation(:,:,op)=nint(integer_transform)
      operation_scale=max(1d0,maxval(abs(integer_transform)))
      if(maxval(abs(integer_transform-real(integer_rotation(:,:,op),real64)))>&
          100d0*tolerance*operation_scale)local_bad=1
      operation_scale=max(1d0,maxval(abs(reciprocal_rotation(:,:,op))))
      if(maxval(abs(matmul(transpose(reciprocal_rotation(:,:,op)),reciprocal_rotation(:,:,op))-&
          identity_real))>100d0*tolerance*operation_scale)local_bad=1
      if(all(integer_rotation(:,:,op)==identity_integer))identity_count=identity_count+1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='reciprocal operation does not preserve the lattice metric';return
    endif
    if(identity_count/=1)then
      call cleanup();message='reciprocal operation catalog requires exactly one identity';return
    endif
    local_bad=0
    do a=1,noperation-1;do b=a+1,noperation
      if(all(integer_rotation(:,:,a)==integer_rotation(:,:,b)))local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='reciprocal operation catalog contains duplicate actions';return
    endif
    local_bad=0
    do a=1,noperation;do b=1,noperation
      product_integer=matmul(integer_rotation(:,:,a),integer_rotation(:,:,b));found=.false.
      do c=1,noperation
        if(all(product_integer==integer_rotation(:,:,c)))then;found=.true.;exit;endif
      enddo
      if(.not.found)local_bad=1
    enddo;enddo
    do a=1,noperation
      found=.false.
      do b=1,noperation
        if(all(matmul(integer_rotation(:,:,a),integer_rotation(:,:,b))==identity_integer).and.&
          all(matmul(integer_rotation(:,:,b),integer_rotation(:,:,a))==identity_integer))then
          found=.true.;exit
        endif
      enddo
      if(.not.found)local_bad=1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='reciprocal operation catalog is not closed';return
    endif

    base_tau=64d0*epsilon(1d0)*max(1d0,abs(cutoff))
    search_cutoff=cutoff+max(4d0*base_tau,256d0*tolerance*max(1d0,abs(cutoff)))
    gmax=sqrt(2d0*search_cutoff)
    do i=1,3
      nmax(i)=ceiling(gmax*sqrt(sum(inverse_lattice(i,:)**2))+10d0*tolerance)
    enddo
    candidate_count=0
    do i=-nmax(1),nmax(1);do j=-nmax(2),nmax(2);do k=-nmax(3),nmax(3)
      g=matmul(reciprocal_lattice,real([i,j,k],real64))
      if(0.5d0*sum(g*g)<=search_cutoff)candidate_count=candidate_count+1
    enddo;enddo;enddo
    if(candidate_count<1)then;message='empty reciprocal catalog';return;endif
    allocate(candidate_integer(3,candidate_count),candidate_vectors(3,candidate_count),&
      candidate_energy(candidate_count),selected(candidate_count),selected_position(candidate_count),&
      stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate reciprocal catalog';return
    endif
    candidate_index=0;selected=.false.;selected_position=0
    do i=-nmax(1),nmax(1);do j=-nmax(2),nmax(2);do k=-nmax(3),nmax(3)
      g=matmul(reciprocal_lattice,real([i,j,k],real64))
      if(0.5d0*sum(g*g)>search_cutoff)cycle
      candidate_index=candidate_index+1
      candidate_integer(:,candidate_index)=[i,j,k]
      candidate_vectors(:,candidate_index)=g
      candidate_energy(candidate_index)=0.5d0*sum(g*g)
      if(candidate_energy(candidate_index)<=cutoff+pw_tolerance(cutoff,candidate_energy(candidate_index)))&
        selected(candidate_index)=.true.
    enddo;enddo;enddo
    if(candidate_index/=candidate_count.or..not.any(selected))then
      call cleanup();message='reciprocal cutoff selection is empty or inconsistent';return
    endif
    boundary_energy=maxval(candidate_energy,mask=selected)
    do i=1,candidate_count
      if(selected(i))cycle
      if(abs(candidate_energy(i)-boundary_energy)<=&
          pw_tolerance(cutoff,candidate_energy(i))+pw_tolerance(cutoff,boundary_energy))then
        selected(i)=.true.;shell_added=shell_added+1
      endif
    enddo
    changed=.true.
    do while(changed)
      changed=.false.
      do i=1,candidate_count
        if(.not.selected(i))cycle
        mapped_integer=-candidate_integer(:,i)
        target_candidate=find_integer_mode(mapped_integer,candidate_integer)
        if(target_candidate==0)then
          call cleanup();message='conjugate reciprocal partner lies outside the safe search catalog';return
        endif
        if(.not.selected(target_candidate))then
          selected(target_candidate)=.true.;orbit_added=orbit_added+1;changed=.true.
        endif
        do op=1,noperation
          mapped_integer=matmul(integer_rotation(:,:,op),candidate_integer(:,i))
          target_candidate=find_integer_mode(mapped_integer,candidate_integer)
          if(target_candidate==0)then
            call cleanup();message='authoritative reciprocal orbit lies outside the safe search catalog';return
          endif
          if(.not.selected(target_candidate))then
            selected(target_candidate)=.true.;orbit_added=orbit_added+1;changed=.true.
          endif
        enddo
      enddo
    enddo
    ng=count(selected)
    allocate(g_integer(3,ng),g_vectors(3,ng),g_action(ng,noperation),g_star(ng),&
      g_conjugate(ng),stat=allocation_status)
    local_bad=merge(0,1,allocation_status==0)
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='cannot allocate completed reciprocal catalog';return
    endif
    j=0
    do i=1,candidate_count
      if(.not.selected(i))cycle
      j=j+1;selected_position(i)=j
      g_integer(:,j)=candidate_integer(:,i);g_vectors(:,j)=candidate_vectors(:,i)
    enddo
    effective_cutoff=maxval(candidate_energy,mask=selected)
    local_bad=0
    do i=1,ng
      call find_mode(-g_vectors(:,i),g_vectors,100d0*tolerance,best,best_distance)
      if(best==0)then;local_bad=1;else;g_conjugate(i)=best;endif
      do op=1,noperation
        mapped_integer=matmul(integer_rotation(:,:,op),g_integer(:,i))
        candidate_index=find_integer_mode(mapped_integer,candidate_integer)
        if(candidate_index==0)then
          local_bad=1
        else
          best=selected_position(candidate_index)
          if(best==0)then;local_bad=1;else;g_action(i,op)=best;endif
        endif
      enddo
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      call cleanup();message='completed reciprocal cutoff is not closed under symmetry and conjugation';return
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
    fingerprint=ieor(ishftc(fingerprint,7),transfer(effective_cutoff,bits))
    fingerprint=ieor(ishftc(fingerprint,7),int(shell_added,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(orbit_added,int64))
    do op=1,noperation;do j=1,3;do i=1,3
      fingerprint=ieor(ishftc(fingerprint,7),int(integer_rotation(i,j,op),int64))
    enddo;enddo;enddo
    do i=1,ng
      do j=1,3;fingerprint=ieor(ishftc(fingerprint,7),int(g_integer(j,i),int64));enddo
      fingerprint=ieor(ishftc(fingerprint,7),int(g_star(i),int64))
      fingerprint=ieor(ishftc(fingerprint,7),int(g_conjugate(i),int64))
      do op=1,noperation
        fingerprint=ieor(ishftc(fingerprint,7),int(g_action(i,op),int64))
      enddo
    enddo
    if(fingerprint==0_int64)fingerprint=1_int64
    deallocate(candidate_integer,candidate_vectors,candidate_energy,selected,selected_position,&
      integer_rotation);ok=.true.
  contains
    subroutine cleanup()
      if(allocated(candidate_integer))deallocate(candidate_integer)
      if(allocated(candidate_vectors))deallocate(candidate_vectors)
      if(allocated(candidate_energy))deallocate(candidate_energy)
      if(allocated(selected))deallocate(selected)
      if(allocated(selected_position))deallocate(selected_position)
      if(allocated(integer_rotation))deallocate(integer_rotation)
      if(allocated(g_integer))deallocate(g_integer)
      if(allocated(g_vectors))deallocate(g_vectors)
      if(allocated(g_action))deallocate(g_action)
      if(allocated(g_star))deallocate(g_star)
      if(allocated(g_conjugate))deallocate(g_conjugate)
    end subroutine cleanup

    integer function find_integer_mode(vector,vectors) result(index)
      integer,intent(in)::vector(3),vectors(:,:)
      integer::mode
      index=0
      do mode=1,size(vectors,2)
        if(all(vectors(:,mode)==vector))then;index=mode;return;endif
      enddo
    end function find_integer_mode

    real(real64) function pw_tolerance(requested,energy) result(value)
      real(real64),intent(in)::requested,energy
      value=64d0*epsilon(1d0)*max(1d0,abs(requested),abs(energy))
    end function pw_tolerance
#else
    ok=.false.;message='reciprocal catalog requires MPI';fingerprint=0_int64
    effective_cutoff=0d0;shell_added=0;orbit_added=0
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

  subroutine agree_real_matrix(values,comm,bad,ierr)
    real(real64),intent(in)::values(:,:)
    integer,intent(in)::comm
    integer,intent(inout)::bad
    integer,intent(out)::ierr
    integer::i,j
    integer(int64)::bits,minimum,maximum
    do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(values(i,j),bits);call agree_int64(bits,minimum,maximum,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(minimum/=maximum)bad=1
    enddo;enddo
  end subroutine agree_real_matrix

  subroutine agree_real_cube(values,comm,bad,ierr)
    real(real64),intent(in)::values(:,:,:)
    integer,intent(in)::comm
    integer,intent(inout)::bad
    integer,intent(out)::ierr
    integer::i,j,k
    integer(int64)::bits,minimum,maximum
    do k=1,size(values,3);do j=1,size(values,2);do i=1,size(values,1)
      bits=transfer(values(i,j,k),bits);call agree_int64(bits,minimum,maximum,comm,ierr)
      if(ierr/=MPI_SUCCESS)return
      if(minimum/=maximum)bad=1
    enddo;enddo;enddo
  end subroutine agree_real_cube

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
