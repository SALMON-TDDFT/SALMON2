! Value-only initial orbitals/occupations and inputs, without MPI/Libxc owners.
module salmonpy_gs_seed
use iso_c_binding
use structures
use salmonpy_model_inputs, only: s_model_inputs
implicit none
private
public :: s_gs_seed,load_gs_seed,s_array_guess,read_array_guess
type :: s_gs_seed
  type(s_model_inputs) :: inputs
  type(s_orbital) :: psi
  real(8),allocatable :: occupations(:,:,:)
  integer :: grid(3)=0,core_start(3)=0,core_end(3)=0,decomposition(5)=0
  integer :: orbital_start=0,orbital_end=-1,k_start=0,k_end=-1
end type
type :: s_array_guess
  type(c_funptr) :: reader=c_null_funptr
  type(c_ptr) :: cookie=c_null_ptr
end type
abstract interface
  integer(c_int) function array_reader(cookie,bounds,shape,waves,occupations,points,weights) bind(C)
    import c_int,c_ptr
    type(c_ptr),value :: cookie,bounds,shape,waves,occupations,points,weights
  end function
end interface
contains
subroutine read_array_guess(data,system,lg,mg,info,psi,seed,status)
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use math_constants,only:pi
 type(s_array_guess),intent(in)::data
 type(s_dft_system),intent(in)::system
 type(s_rgrid),intent(in)::lg,mg
 type(s_parallel_info),intent(in)::info
 type(s_orbital),intent(in)::psi
 type(s_gs_seed),intent(out)::seed
 integer,intent(out)::status
 procedure(array_reader),pointer::reader
 integer(c_int64_t),target::bounds(10),shape5(5)
 complex(c_double_complex),allocatable,target::waves(:)
 real(c_double),allocatable,target::occ(:),points(:,:),weights(:)
 real(8)::expected
 integer::nx,ny,nz,nb,nk,nlocal,i,j,k,o,x,y,z,idx,allocation
 status=11
 if(.not.c_associated(data%reader))return
 if(system%nspin/=1.or.info%im_s/=1.or.info%im_e/=1)return
 nx=mg%ie(1)-mg%is(1)+1;ny=mg%ie(2)-mg%is(2)+1;nz=mg%ie(3)-mg%is(3)+1
 nb=system%no;nk=system%nk
 bounds=[int(mg%is,c_int64_t),int(mg%ie,c_int64_t), &
         int([info%io_s,info%io_e,info%ik_s,info%ik_e],c_int64_t)]
 shape5=[int([nx,ny,nz,nb,nk],c_int64_t)]
 if(any(shape5<1).or.info%numo<1.or.info%numk<1)return
 if(product(real(shape5(:3),8))*info%numo*info%numk>real(huge(nlocal),8))return
 nlocal=nx*ny*nz*info%numo*info%numk
 allocate(waves(nlocal),occ(nb*nk),points(3,nk),weights(nk),stat=allocation)
 if(allocation/=0)then
  status=12;return
 endif
 ! Poison outputs so an incomplete provider cannot silently initialize a state.
 block
  use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan
  waves=cmplx(ieee_value(0d0,ieee_quiet_nan),0d0,kind=8)
  occ=ieee_value(0d0,ieee_quiet_nan);points=ieee_value(0d0,ieee_quiet_nan)
  weights=ieee_value(0d0,ieee_quiet_nan)
 end block
 call c_f_procpointer(data%reader,reader)
 status=reader(data%cookie,c_loc(bounds),c_loc(shape5),c_loc(waves),c_loc(occ),c_loc(points),c_loc(weights))
 if(status/=0)then
  status=11;return
 endif
 status=11
 do i=1,nlocal
  if(.not.ieee_is_finite(real(waves(i),8)).or..not.ieee_is_finite(aimag(waves(i))))return
  if(allocated(psi%rwf))then
   if(aimag(waves(i))/=0d0)return
  endif
 enddo
 do k=1,nk
  if(.not.ieee_is_finite(weights(k)))return
  if(abs(weights(k)-system%wtk(k))>1d-12)return
  do j=1,3
   if(.not.ieee_is_finite(points(j,k)))return
   expected=dot_product(system%primitive_a(:,j),system%vec_k(:,k))/(2d0*pi)
   if(abs(points(j,k)-expected)>1d-12)return
  enddo
  do o=1,nb
   if(.not.ieee_is_finite(occ((k-1)*nb+o)))return
   if(occ((k-1)*nb+o)<0d0)return
  enddo
 enddo
 call seed%inputs%capture()
 seed%grid=lg%num;seed%core_start=mg%is;seed%core_end=mg%ie
 seed%decomposition=[info%npk,info%nporbital,info%nprgrid]
 seed%orbital_start=info%io_s;seed%orbital_end=info%io_e
 seed%k_start=info%ik_s;seed%k_end=info%ik_e
 allocate(seed%occupations(nb,nk,1),stat=allocation)
 if(allocation/=0)then
  status=12;return
 endif
 seed%occupations(:,:,1)=reshape(occ,[nb,nk])
 if(allocated(psi%rwf))then
  allocate(seed%psi%rwf(lbound(psi%rwf, &
    1):ubound(psi%rwf, &
    1), &
    lbound(psi%rwf, &
    2):ubound(psi%rwf, &
    2), &
    lbound(psi%rwf, &
    3):ubound(psi%rwf, &
    3), &
    lbound(psi%rwf, &
    4):ubound(psi%rwf, &
    4), &
    lbound(psi%rwf, &
    5):ubound(psi%rwf, &
    5), &
    lbound(psi%rwf, &
    6):ubound(psi%rwf, &
    6), &
    lbound(psi%rwf, &
    7):ubound(psi%rwf, &
    7)), &
    stat=allocation)
  if(allocation/=0)then
   status=12;return
  endif
  seed%psi%rwf=0d0
 endif
 if(allocated(psi%zwf))then
  allocate(seed%psi%zwf(lbound(psi%zwf, &
    1):ubound(psi%zwf, &
    1), &
    lbound(psi%zwf, &
    2):ubound(psi%zwf, &
    2), &
    lbound(psi%zwf, &
    3):ubound(psi%zwf, &
    3), &
    lbound(psi%zwf, &
    4):ubound(psi%zwf, &
    4), &
    lbound(psi%zwf, &
    5):ubound(psi%zwf, &
    5), &
    lbound(psi%zwf, &
    6):ubound(psi%zwf, &
    6), &
    lbound(psi%zwf, &
    7):ubound(psi%zwf, &
    7)), &
    stat=allocation)
  if(allocation/=0)then
   status=12;return
  endif
  seed%psi%zwf=0d0
 endif
 idx=0
 do k=info%ik_s,info%ik_e
 do o=info%io_s,info%io_e
 do z=mg%is(3),mg%ie(3)
 do y=mg%is(2),mg%ie(2)
 do x=mg%is(1),mg%ie(1)
  idx=idx+1
  if(allocated(seed%psi%rwf))seed%psi%rwf(x,y,z,1,o,k,1)=real(waves(idx),8)
  if(allocated(seed%psi%zwf))seed%psi%zwf(x,y,z,1,o,k,1)=waves(idx)
 enddo
 enddo
 enddo
 enddo
 enddo
 status=0
end subroutine
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
