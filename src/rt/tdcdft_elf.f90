! Copyright 2026 SALMON developers. Apache License, Version 2.0.
module tdcdft_elf
  implicit none
  private
  public :: measure_elf
contains
  subroutine measure_elf(system,mg,info,stencil,srg,psi,measure)
    use structures, only: s_dft_system,s_rgrid,s_parallel_info,s_stencil,s_sendrecv_grid,s_orbital
    use stencil_sub, only: calc_gradient_psi
    use sendrecv_grid, only: update_overlap_complex8
    use communication, only: comm_summation
    use tdcdft_lrc, only: elf_value
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(s_dft_system),intent(in) :: system
    type(s_rgrid),intent(in) :: mg
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_orbital),intent(inout) :: psi
    real(8),intent(out) :: measure
    real(8),allocatable :: local(:,:,:,:),moments(:,:,:,:)
    complex(8),allocatable :: orbital(:,:,:),gradient(:,:,:,:)
    complex(8) :: val,g(3),overlap(3)
    real(8) :: weight,n,elf,local_sum(2),total(2)
    integer :: ik,io,ix,iy,iz
    if(system%nspin/=1.or.info%im_s/=1.or.info%im_e/=1) error stop 'TDCDFT ELF: unsupported spin or image'
    if(.not.allocated(psi%zwf)) error stop 'TDCDFT ELF: complex orbitals required'
    allocate(local(8,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))

    allocate(moments(8,mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3)))
    allocate(orbital(mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
                     mg%is_array(3):mg%ie_array(3)))
    allocate(gradient(3,mg%is_array(1):mg%ie_array(1),mg%is_array(2):mg%ie_array(2), &
                       mg%is_array(3):mg%ie_array(3)))
    if(info%if_divide_rspace) call update_overlap_complex8(srg,mg,psi%zwf)
    local=0d0
    do ik=info%ik_s,info%ik_e
    do io=info%io_s,info%io_e
      weight=0.5d0*system%rocc(io,ik,1)*system%wtk(ik)
      if(weight==0d0) cycle
      orbital=psi%zwf(:,:,:,1,io,ik,1)
      call calc_gradient_psi(orbital,gradient,mg%is_array,mg%ie_array,mg%is,mg%ie, &
                            mg%idx,mg%idy,mg%idz,stencil%coef_nab,system%rmatrix_B)
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        val=orbital(ix,iy,iz)
        g=gradient(:,ix,iy,iz)+cmplx(0d0,1d0,8)*system%vec_k(:,ik)*val
        overlap=conjg(val)*g
        local(1,ix,iy,iz)=local(1,ix,iy,iz)+weight*abs(val)**2
        local(2,ix,iy,iz)=local(2,ix,iy,iz)+weight*sum(abs(g)**2)
        local(3:5,ix,iy,iz)=local(3:5,ix,iy,iz)+weight*real(overlap,8)
        local(6:8,ix,iy,iz)=local(6:8,ix,iy,iz)+weight*aimag(overlap)
      end do
      end do
      end do
    end do
    end do
    call comm_summation(local,moments,size(local),info%icomm_ko)
    if(.not.all(ieee_is_finite(moments))) error stop 'TDCDFT ELF: nonfinite orbital moments'
    local_sum=0d0
    do iz=mg%is(3),mg%ie(3)
    do iy=mg%is(2),mg%ie(2)
    do ix=mg%is(1),mg%ie(1)
      n=moments(1,ix,iy,iz)
      elf=elf_value(n,moments(2,ix,iy,iz),moments(3:5,ix,iy,iz),moments(6:8,ix,iy,iz))
      local_sum=local_sum+[n*(elf-0.5d0)**2,n]
    end do
    end do
    end do
    call comm_summation(local_sum,total,2,info%icomm_r)
    if(total(2)<=0d0) error stop 'TDCDFT ELF: empty density'
    measure=total(1)/total(2) ! Cell volume and spin multiplicity cancel.
  end subroutine measure_elf
end module tdcdft_elf
