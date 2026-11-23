program test_soc_resources
 use noncollinear_module,only:s_noncollinear_state,swap_noncollinear_state,snapshot_noncollinear_state,noncollinear_finite
 use,intrinsic::ieee_arithmetic,only:ieee_value,ieee_quiet_nan,ieee_positive_inf,ieee_negative_inf
 use structures,only:s_kpoint_config,s_dft_system,s_unfold
 use lattice,only:init_kvector_config
 use parallelization,only:setup_parallel,end_parallel
 use sym_sub,only:use_symmetry,DISPLAY
 implicit none
 type(s_noncollinear_state)::left,right,snapshot,observed,finite_seed
 type(s_kpoint_config)::config
 type(s_dft_system)::system
 type(s_unfold)::unfold
 integer::i,j,k
 real(8)::bad
 call setup_parallel
 config%num_kgrid=[2,2,1];config%dk_shift=0d0
 config%yn_gamma_centered='n';config%file_kw='none';config%dm_unfold_option='no'
 system%primitive_b=0d0
 do i=1,3
  system%primitive_b(i,i)=1d0
 enddo
 use_symmetry=.false.;DISPLAY=.false.
 call init_kvector_config(system,unfold,config)
 if(system%nk/=4)error stop 'expected four k points'
 do i=1,system%nk
  if(maxval(abs(system%vec_k(:,i)))<1d-12)error stop 'Gamma point is present'
  print *, 'K_FRACTIONAL',i,system%vec_k(:,i)
 enddo
 print *, 'GAMMA_FREE_2X2X1_PASS'
 allocate(left%den_mat(-1:0,1:1,1:1,1:2,1:2))
 allocate(left%vxc_mat(-1:0,1:1,1:1,1:2,1:2))
 allocate(left%old_mat(-1:0,1:1,1:1))
 allocate(left%rot_ang(-1:0,1:1,1:1,1:2))
 allocate(left%dmat_old(-1:0,1:1,1:1,1:2,1:2))
 left%den_mat=(1d0,2d0);left%vxc_mat=(3d0,4d0);left%old_mat=(5d0,6d0)
 left%rot_ang=7d0;left%dmat_old=(8d0,9d0)
 right=left;right%den_mat=(11d0,12d0)
 call swap_noncollinear_state(left)
 if(allocated(left%den_mat))error stop 'owner transfer did not empty source'
 call snapshot_noncollinear_state(snapshot)
 if(any(lbound(snapshot%den_mat)/=[-1,1,1,1,1]))error stop 'bounds lost'
 call swap_noncollinear_state(right)
 if(right%den_mat(-1,1,1,1,1)/=(1d0,2d0))error stop 'first owner lost'
 call snapshot_noncollinear_state(observed)
 if(observed%den_mat(-1,1,1,1,1)/=(11d0,12d0))error stop 'second owner not selected'
 call swap_noncollinear_state(left)
 if(left%den_mat(-1,1,1,1,1)/=(11d0,12d0))error stop 'second owner lost'
 call swap_noncollinear_state(snapshot)
 call snapshot_noncollinear_state(observed)
 if(observed%den_mat(-1,1,1,1,1)/=(1d0,2d0))error stop 'snapshot restore failed'
 if(observed%vxc_mat(-1,1,1,1,1)/=(3d0,4d0))error stop 'xc state lost'
 if(observed%old_mat(-1,1,1)/=(5d0,6d0))error stop 'old matrix lost'
 if(observed%rot_ang(-1,1,1,1)/=7d0)error stop 'rotation state lost'
 if(observed%dmat_old(-1,1,1,1,1)/=(8d0,9d0))error stop 'mixing history lost'
 observed%den_mat=(99d0,0d0)
 call snapshot_noncollinear_state(observed)
 if(observed%den_mat(-1,1,1,1,1)/=(1d0,2d0))error stop 'snapshot aliases live state'
 if(.not.noncollinear_finite())error stop 'finite SOC workspace rejected'
 finite_seed=observed
 do j=1,3
  select case(j)
  case(1);bad=ieee_value(0d0,ieee_quiet_nan)
  case(2);bad=ieee_value(0d0,ieee_positive_inf)
  case(3);bad=ieee_value(0d0,ieee_negative_inf)
  end select
  do k=1,5
   observed=finite_seed
   select case(k)
   case(1);observed%den_mat(0,1,1,2,2)=cmplx(bad,0d0,8)
   case(2);observed%vxc_mat(0,1,1,2,2)=cmplx(0d0,bad,8)
   case(3);observed%old_mat(0,1,1)=cmplx(bad,0d0,8)
   case(4);observed%rot_ang(0,1,1,2)=bad
   case(5);observed%dmat_old(0,1,1,2,2)=cmplx(0d0,bad,8)
   end select
   call swap_noncollinear_state(observed)
   if(noncollinear_finite())error stop 'nonfinite SOC workspace accepted'
   call swap_noncollinear_state(observed)
   if(.not.noncollinear_finite())error stop 'finite SOC workspace was not preserved'
  enddo
 enddo
 print *, 'SOC finite/NaN/Inf scalar inquiry PASS'
 print *, 'SOC resource ownership PASS'
 call end_parallel
end program
