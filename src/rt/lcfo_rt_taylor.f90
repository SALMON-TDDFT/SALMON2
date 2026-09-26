! Untruncated coefficient-space Taylor reference; native Hamiltonian is unchanged.
module lcfo_rt_taylor
  implicit none
  private
  public :: lcfo_coefficient_taylor
contains
  subroutine lcfo_coefficient_taylor(mg,system,info,stencil,srg,input,output,work,ppg,vlocal,rt)
    use structures
    use lcfo_rt_basis, only: lcfo_basis,lcfo_dv,lcfo_rank,lcfo_orb_rank
    use hamiltonian, only: hpsi
    use salmon_global, only: n_hamil
    implicit none
    type(s_rgrid),intent(in) :: mg
    type(s_dft_system),intent(in) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_sendrecv_grid),intent(inout) :: srg
    type(s_orbital),intent(inout) :: input,output,work
    type(s_pp_grid),intent(in) :: ppg
    type(s_scalar),intent(in) :: vlocal(system%nspin)
    type(s_rt),intent(in) :: rt
    complex(8),allocatable :: power(:,:),result(:,:),grid(:,:),next_power(:,:)
    integer :: n,j,io,ng,is(3),ie(3)
    if(n_hamil/=4)error stop 'LCFO direct WF requires Taylor4'
    ng=product(mg%num);is=mg%is;ie=mg%ie
    allocate(grid(ng,info%numo))
    do j=1,info%numo
      io=info%io_s+j-1
      grid(:,j)=reshape(input%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1),[ng])
    enddo
    power=matmul(conjg(transpose(lcfo_basis)),grid)*lcfo_dv
    result=power
    allocate(next_power(size(power,1),size(power,2)))
    do n=1,4
      grid=matmul(lcfo_basis,power)
      do j=1,info%numo
        io=info%io_s+j-1
        input%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)=reshape(grid(:,j),mg%num)
      enddo
      input%update_zwf_overlap=.false.
      call hpsi(input,work,info,mg,vlocal,system,stencil,srg,ppg,lcfo_coeff=power,lcfo_action=next_power)
      power=next_power
      result=result+rt%zc(n)*power
    enddo
    grid=matmul(lcfo_basis,result)
    do j=1,info%numo
      io=info%io_s+j-1
      output%zwf(is(1):ie(1),is(2):ie(2),is(3):ie(3),1,io,1,1)=reshape(grid(:,j),mg%num)
    enddo
    output%update_zwf_overlap=.false.
    if(lcfo_rank==0.and.lcfo_orb_rank==0)write(*,'(a)')'LCFO direct WF coefficient Taylor4'
  end subroutine
end module
