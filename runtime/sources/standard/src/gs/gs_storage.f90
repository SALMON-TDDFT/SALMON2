! Stable owned GS numerical storage; never copy a live owner.
module salmonpy_gs_storage
use structures
use salmon_xc, only: finalize_xc
use communication, only: comm_release_children
implicit none
private
public :: s_gs_storage,gs_storage_clear
type :: s_gs_storage
  logical :: active=.false.,xc_live=.false.,converged=.false.,imported=.false.
  integer :: mpi_child_owner=0,iterations=0,report_status=0
  real(8) :: residual=0d0,threshold=0d0
type(s_rgrid) :: lg
type(s_rgrid) :: mg
type(s_parallel_info) :: info
type(s_sendrecv_grid) :: srg, srg_scalar
type(s_orbital) :: spsi,shpsi,sttpsi
type(s_dft_system) :: system
type(s_poisson) :: poisson
type(s_stencil) :: stencil
type(s_xc_functional) :: xc_func
type(s_scalar) :: rho,rho_jm,Vh,Vpsl
type(s_scalar),allocatable :: V_local(:),rho_s(:),Vxc(:)
type(s_reciprocal_grid) :: fg
type(s_pp_info) :: pp
type(s_pp_grid) :: ppg
type(s_pp_nlcc) :: ppn
type(s_dft_energy) :: energy
type(s_ewald_ion_ion) :: ewald
type(s_cg)     :: cg
type(s_mixing) :: mixing
type(s_ofile)  :: ofl
type(s_band_dft) ::band
type(s_opt) :: opt
type(s_dcdft) :: dc
type(s_unfold) :: unfold
end type
contains
subroutine gs_storage_clear(ctx)
  type(s_gs_storage),intent(inout) :: ctx
  type(s_gs_storage) :: empty
  if(ctx%xc_live)call finalize_xc(ctx%xc_func)
  ! Request caches have already retired before gs_solve returns.
  call comm_release_children(ctx%mpi_child_owner)
  ctx=empty
end subroutine
end module
