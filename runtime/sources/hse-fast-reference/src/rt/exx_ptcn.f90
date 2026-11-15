module exx_ptcn
 implicit none
 private
 public::native_exx_step
contains
  subroutine native_exx_step(dt,lg,mg,system,info,stencil,xc_func,srg,srg_scalar,pp,ppg,ppn, &
      input,output,rho,rho_s,vlocal,vh,vxc,vpsl,fg,poisson,energy,ac_start,ac_end)
    use structures, only: s_rgrid,s_dft_system,s_parallel_info,s_stencil, &
      s_xc_functional,s_sendrecv_grid,s_pp_info,s_pp_grid, &
      s_pp_nlcc,s_orbital,s_scalar,s_reciprocal_grid, &
      s_poisson,s_dft_energy,allocate_orbital_complex,deallocate_orbital
    use ptcn_native,only:native_ptcn_step
    implicit none
    real(8),intent(in) :: dt
    type(s_rgrid),intent(in) :: lg,mg
    type(s_dft_system),intent(inout) :: system
    type(s_parallel_info),intent(in) :: info
    type(s_stencil),intent(in) :: stencil
    type(s_xc_functional),intent(in) :: xc_func
    type(s_sendrecv_grid),intent(inout) :: srg,srg_scalar
    type(s_pp_info),intent(in) :: pp
    type(s_pp_grid),intent(inout) :: ppg
    type(s_pp_nlcc),intent(in) :: ppn
    type(s_orbital),intent(in) :: input
    type(s_orbital),intent(inout) :: output
    type(s_scalar),intent(inout) :: rho,rho_s(:),vlocal(:),vh,vxc(:)
    type(s_scalar),intent(in) :: vpsl
    type(s_reciprocal_grid),intent(inout) :: fg
    type(s_poisson),intent(inout) :: poisson
    type(s_dft_energy),intent(inout) :: energy
    real(8),intent(in),optional::ac_start(3),ac_end(3)
    real(8)::a0(3),a1(3)
    a0=system%vec_Ac;a1=a0
    if(present(ac_start))a0=ac_start
    if(present(ac_end))a1=ac_end
    call native_ptcn_step(dt,a0,a1,lg,mg,system,info,stencil,xc_func,srg,srg_scalar,pp,ppg,ppn, &
      input,output,rho,rho_s,vlocal,vh,vxc,vpsl,fg,poisson,energy)
  end subroutine
end module
