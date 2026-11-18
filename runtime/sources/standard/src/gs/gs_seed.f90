! Value-only initial orbitals/occupations and inputs, without MPI/Libxc owners.
module salmonpy_gs_seed
use structures
use salmonpy_model_inputs, only: s_model_inputs
implicit none
private
public :: s_gs_seed,load_gs_seed
type :: s_gs_seed
  type(s_model_inputs) :: inputs
  type(s_orbital) :: psi
  real(8),allocatable :: occupations(:,:,:)
  integer :: grid(3)=0,core_start(3)=0,core_end(3)=0,decomposition(5)=0
  integer :: orbital_start=0,orbital_end=-1,k_start=0,k_end=-1
end type
contains
subroutine load_gs_seed(seed,system,mg,info,psi)
  use communication, only: comm_summation
  type(s_gs_seed),intent(in) :: seed
  type(s_dft_system),intent(inout) :: system
  type(s_rgrid),intent(in) :: mg
  type(s_parallel_info),intent(in) :: info
  type(s_orbital),intent(inout) :: psi
  complex(8),allocatable :: source(:,:,:),destination(:,:,:)
  integer :: o,k,s,m
  ! Inputs and the subgroup/decomposition are inherited from the source model.
  system%rocc=seed%occupations(:system%no,:system%nk,:system%nspin)
  psi%zwf=(0d0,0d0)
  allocate(source(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  allocate(destination(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
  do m=info%im_s,info%im_e
  do k=info%ik_s,info%ik_e
  do s=1,system%nspin
  do o=1,system%no
    source=(0d0,0d0)
    if(o>=seed%orbital_start.and.o<=seed%orbital_end)then
      if(allocated(seed%psi%rwf))then
        source=cmplx(seed%psi%rwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2), &
          mg%is(3):mg%ie(3),s,o,k,m),0d0,kind=8)
      else
        source=seed%psi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),s,o,k,m)
      endif
    endif
    ! GS may carry empty bands and thus a different orbital ownership split.
    call comm_summation(source,destination,size(source),info%icomm_o)
    if(o>=info%io_s.and.o<=info%io_e) &
      psi%zwf(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3),s,o,k,m)=destination
  enddo
  enddo
  enddo
  enddo
end subroutine
end module
