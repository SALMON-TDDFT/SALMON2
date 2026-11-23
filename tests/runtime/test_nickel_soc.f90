! Small Ni SOC execution/ownership test. No GS orbital/density export or RT.
program test_nickel_soc
 use salmon_global
 use parallelization
 use inputoutput
 use math_constants
 use timer
 use plusU_global,only:read_Hubbard_parameters
 use sym_sub,only:read_sw_symmetry
 use salmonpy_gs_context,only:s_gs_context,gs_run,gs_close
 use noncollinear_module,only:s_noncollinear_state,snapshot_noncollinear_state,swap_noncollinear_state
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 implicit none
 type(s_gs_context),target::context
 type(s_noncollinear_state)::saved,changed,restored
 integer::status,i,j,k
 real(8)::volume,electrons,magnetization(3),hermiticity
 logical::exists
 call set_math_constants
 call setup_parallel
 call read_input
 calc_mode='GS'
 call read_Hubbard_parameters
 call read_sw_symmetry(yn_symmetry)
 call timer_initialize
 if(yn_spinorbit/='y'.or.spin/='noncollinear')error stop 'Ni test requires SOC spinors'
 if(any(num_kgrid/=[2,2,2]))error stop 'Ni test requires 2x2x2'
 if(yn_gamma_centered/='n'.or.any(dk_shift/=0d0))error stop 'SALMON standard k division required'
 if(xc/='libxc_pbe'.or.natom/=1.or.nelem/=1.or.izatom(1)/=28)error stop 'Ni PBE input required'
 if(nelec/=18.or.nstate/=2*nelec.or.temperature_k/=300d0)error stop 'Ni metal occupations required'
 if(method_init_wf/='random'.or.method_mixing/='simple_dm')error stop 'explicit spinor initialization required'
 if(write_gs_restart_data/='no'.or.yn_self_checkpoint/='n')error stop 'CI must not export GS data'
 if(nscf/=3)error stop 'CI has exactly three SCF iterations'
 call gs_run(context,status)
 if(status/=0)error stop 'Ni SCF execution failed'
 call snapshot_noncollinear_state(saved)
 if(.not.allocated(saved%den_mat).or..not.allocated(saved%vxc_mat))error stop 'SOC state missing'
 if(.not.allocated(saved%rot_ang).or..not.allocated(saved%dmat_old))error stop 'SOC mixing state missing'
 if(.not.finite5(saved%den_mat).or..not.finite5(saved%vxc_mat))error stop 'nonfinite SOC state'
 if(.not.finite5(saved%dmat_old))error stop 'nonfinite SOC mixing history'
 hermiticity=maxval(abs(saved%den_mat(:,:,:,1,2)-conjg(saved%den_mat(:,:,:,2,1))))
 if(hermiticity>1d-12)error stop 'density matrix is not Hermitian'
 volume=abs(al_vec1(1)*(al_vec2(2)*al_vec3(3)-al_vec2(3)*al_vec3(2)) &
           -al_vec1(2)*(al_vec2(1)*al_vec3(3)-al_vec2(3)*al_vec3(1)) &
           +al_vec1(3)*(al_vec2(1)*al_vec3(2)-al_vec2(2)*al_vec3(1)))
 electrons=0d0;magnetization=0d0
 do k=1,size(saved%den_mat,3)
 do j=1,size(saved%den_mat,2)
 do i=1,size(saved%den_mat,1)
  electrons=electrons+dble(saved%den_mat(i,j,k,1,1)+saved%den_mat(i,j,k,2,2))
  magnetization(1)=magnetization(1)+dble(saved%den_mat(i,j,k,1,2))
  magnetization(2)=magnetization(2)-aimag(saved%den_mat(i,j,k,1,2))
  magnetization(3)=magnetization(3)+0.5d0*dble(saved%den_mat(i,j,k,1,1)-saved%den_mat(i,j,k,2,2))
 enddo
 enddo
 enddo
 electrons=electrons*volume/dble(product(num_rgrid))
 magnetization=magnetization*volume/dble(product(num_rgrid))
 if(.not.ieee_is_finite(electrons).or.abs(electrons-18d0)>1d-7)error stop 'electron-number mismatch'
 do i=1,3
  if(.not.ieee_is_finite(magnetization(i)))error stop 'nonfinite magnetization'
 enddo
 changed=saved
 changed%den_mat(1,1,1,1,1)=changed%den_mat(1,1,1,1,1)+0.125d0
 call swap_noncollinear_state(changed)
 call snapshot_noncollinear_state(restored)
 if(all(restored%den_mat==saved%den_mat))error stop 'changed owner was not selected'
 call swap_noncollinear_state(saved)
 call snapshot_noncollinear_state(restored)
 if(any(restored%den_mat/=changed%den_mat))error stop 'density restore differs'
 if(any(restored%vxc_mat/=changed%vxc_mat))error stop 'XC restore differs'
 if(any(restored%rot_ang/=changed%rot_ang))error stop 'rotation restore differs'
 if(any(restored%dmat_old/=changed%dmat_old))error stop 'mixing restore differs'
 if(allocated(restored%old_mat).neqv.allocated(changed%old_mat))error stop 'old-matrix allocation differs'
 if(allocated(restored%old_mat))then
  if(any(restored%old_mat/=changed%old_mat))error stop 'old-matrix restore differs'
 endif
 inquire(file='data_for_restart/wfn.bin',exist=exists)
 if(exists)error stop 'GS wavefunction was exported'
 inquire(file='data_for_restart/rho.bin',exist=exists)
 if(exists)error stop 'GS density was exported'
 print *, 'NI_SOC_ELECTRONS',electrons
 print *, 'NI_SOC_SPIN_EXPECTATION',magnetization
 print *, 'NI_SOC_HERMITICITY',hermiticity
 print *, 'NI_SOC_RESTORE_EXACT_PASS'
 print *, 'NI_SOC_SHORT_SCF_PASS'
 call gs_close(context)
 call end_parallel
contains
 logical function finite5(value)result(ok)
  complex(8),intent(in)::value(:,:,:,:,:)
  integer::i1,i2,i3,i4,i5
  ok=.false.
  do i5=1,size(value,5)
  do i4=1,size(value,4)
  do i3=1,size(value,3)
  do i2=1,size(value,2)
  do i1=1,size(value,1)
   if(.not.ieee_is_finite(dble(value(i1,i2,i3,i4,i5))))return
   if(.not.ieee_is_finite(aimag(value(i1,i2,i3,i4,i5))))return
  enddo
  enddo
  enddo
  enddo
  enddo
  ok=.true.
 end function
end program
