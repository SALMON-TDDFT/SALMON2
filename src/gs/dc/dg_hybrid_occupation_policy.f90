#include "config.h"
module dg_hybrid_occupation_policy
  use,intrinsic::iso_fortran_env,only:int64,real64
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use occupation_kernel,only:solve_spectrum_occupations
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  private
  type,public::s_dg_hybrid_occupation_result
    logical::valid=.false.
    integer::state_count=0,noccupied=0
    real(real64)::chemical_potential=0d0,electron_count=0d0,e_homo=0d0
    real(real64)::occupation_threshold=64d0*epsilon(1d0),omitted_occupation_tail=0d0
    integer(int64)::fingerprint=0_int64
    real(real64),allocatable::occupations(:)
  end type s_dg_hybrid_occupation_result
  public::derive_dg_hybrid_occupation_policy
contains
  subroutine derive_dg_hybrid_occupation_policy(comm,eigenvalues,expected_electron_count,&
      electronic_temperature,electron_count_tolerance,result,ok,message)
    integer,intent(in)::comm
    real(real64),intent(in)::eigenvalues(:),expected_electron_count,electronic_temperature,&
      electron_count_tolerance
    type(s_dg_hybrid_occupation_result),intent(out)::result
    logical,intent(out)::ok
    character(*),intent(out)::message
    real(real64),allocatable::spectrum(:,:,:),weights(:),kernel_occupations(:,:,:)
    real(real64)::chemical_potential,electron_count
    integer::i,ierr,minimum_integer,maximum_integer,local_bad,global_bad
    integer(int64)::bits,minimum_bits,maximum_bits,fingerprint
    logical::kernel_ok
    character(256)::kernel_message

    result=s_dg_hybrid_occupation_result();ok=.false.;message=''
#ifdef USE_MPI
    call MPI_Allreduce(size(eigenvalues),minimum_integer,1,MPI_INTEGER,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Hybrid occupation spectrum rank agreement failed';return;endif
    call MPI_Allreduce(size(eigenvalues),maximum_integer,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_integer/=maximum_integer)then
      message='Hybrid occupation spectrum rank agreement failed';return
    endif
    call agree_real(expected_electron_count,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='Hybrid electron target rank agreement failed';return
    endif
    call agree_real(electronic_temperature,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='Hybrid temperature rank agreement failed';return
    endif
    call agree_real(electron_count_tolerance,minimum_bits,maximum_bits,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='Hybrid electron tolerance rank agreement failed';return
    endif
    do i=1,size(eigenvalues)
      call agree_real(eigenvalues(i),minimum_bits,maximum_bits,ierr)
      if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
        message='Hybrid final eigenvalues are rank-dependent';return
      endif
    enddo
#endif
    local_bad=merge(0,1,size(eigenvalues)>0.and.all(ieee_is_finite(eigenvalues)).and.&
      ieee_is_finite(expected_electron_count).and.expected_electron_count>=0d0.and.&
      ieee_is_finite(electronic_temperature).and.electronic_temperature>=0d0.and.&
      ieee_is_finite(electron_count_tolerance).and.electron_count_tolerance>0d0)
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='invalid Hybrid occupation policy input';return
    endif
#else
    global_bad=local_bad
    if(global_bad/=0)then;message='invalid Hybrid occupation policy input';return;endif
#endif
    if(size(eigenvalues)>1)then
      local_bad=merge(0,1,all(eigenvalues(2:)>=eigenvalues(:size(eigenvalues)-1)))
    else
      local_bad=0
    endif
#ifdef USE_MPI
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then
      message='Hybrid final LCFO spectrum must be ascending';return
    endif
#else
    if(local_bad/=0)then;message='Hybrid final LCFO spectrum must be ascending';return;endif
#endif
    allocate(spectrum(size(eigenvalues),1,1),weights(1));spectrum(:,1,1)=eigenvalues;weights=1d0
    call solve_spectrum_occupations(spectrum,weights,expected_electron_count,electronic_temperature,.false.,&
      kernel_occupations,chemical_potential,electron_count,kernel_ok,kernel_message)
    if(.not.kernel_ok)then;message=trim(kernel_message);return;endif
    result%state_count=size(eigenvalues)
    result%noccupied=count(kernel_occupations(:,1,1)>result%occupation_threshold)
    if(result%noccupied<1)then;message='Hybrid final spectrum has no occupied state';return;endif
    result%omitted_occupation_tail=sum(kernel_occupations(result%noccupied+1:,1,1))
    if(result%omitted_occupation_tail>electron_count_tolerance)then
      write(message,'(a,2(es24.16,a))')'Hybrid omitted occupation tail exceeds tolerance: tail=',&
        result%omitted_occupation_tail,' tolerance=',electron_count_tolerance,''
      return
    endif
    if(abs(electron_count-expected_electron_count)>electron_count_tolerance)then
      message='Hybrid occupation electron count exceeds tolerance';return
    endif
    allocate(result%occupations(size(eigenvalues)))
    result%occupations=kernel_occupations(:,1,1)
    result%chemical_potential=chemical_potential;result%electron_count=electron_count
    result%e_homo=eigenvalues(result%noccupied)
    fingerprint=int(z'5BE0CD19137E2179',int64)
    fingerprint=ieor(ishftc(fingerprint,7),int(result%state_count,int64))
    fingerprint=ieor(ishftc(fingerprint,7),int(result%noccupied,int64))
    do i=1,size(eigenvalues)
      bits=transfer(eigenvalues(i),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
      bits=transfer(result%occupations(i),bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
    enddo
    bits=transfer(result%chemical_potential,bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
    bits=transfer(result%e_homo,bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
    bits=transfer(result%omitted_occupation_tail,bits);fingerprint=ieor(ishftc(fingerprint,7),bits)
    if(fingerprint==0_int64)fingerprint=1_int64
#ifdef USE_MPI
    call MPI_Allreduce(fingerprint,minimum_bits,1,MPI_INTEGER8,MPI_MIN,comm,ierr)
    if(ierr/=MPI_SUCCESS)then;message='Hybrid occupation fingerprint reduction failed';return;endif
    call MPI_Allreduce(fingerprint,maximum_bits,1,MPI_INTEGER8,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.minimum_bits/=maximum_bits)then
      message='Hybrid occupation fingerprint is rank-dependent';return
    endif
#endif
    result%fingerprint=fingerprint;result%valid=.true.;ok=.true.
  contains
#ifdef USE_MPI
    subroutine agree_real(value,minimum_value,maximum_value,status)
      real(real64),intent(in)::value
      integer(int64),intent(out)::minimum_value,maximum_value
      integer,intent(out)::status
      integer(int64)::value_bits
      value_bits=transfer(value,value_bits)
      call MPI_Allreduce(value_bits,minimum_value,1,MPI_INTEGER8,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS)return
      call MPI_Allreduce(value_bits,maximum_value,1,MPI_INTEGER8,MPI_MAX,comm,status)
    end subroutine agree_real
#endif
  end subroutine derive_dg_hybrid_occupation_policy
end module dg_hybrid_occupation_policy
