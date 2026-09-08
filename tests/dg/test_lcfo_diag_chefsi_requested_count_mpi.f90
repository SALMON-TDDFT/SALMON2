program test_lcfo_diag_chefsi_requested_count_mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi
  use structures, only: s_dcdft
  use lcfo_diag_chefsi, only: diag_chefsi
  implicit none
  integer, parameter :: basis_count=3,nspin=1
  type(s_dcdft) :: dc
  integer :: ierr,rank,nproc,i,j
  integer :: n_basis(1,nspin),n_mat(nspin)
  integer :: halo_src(0),halo_dst(0),halo_root_src(0),halo_dvec(3,0)
  real(8) :: h_diag(basis_count,basis_count,nspin)
  real(8) :: h_halo(basis_count,basis_count,nspin,0)
  real(8) :: esp_expanded(basis_count,nspin),esp_conventional(basis_count,nspin)
  real(8), allocatable :: coef_expanded(:,:,:),coef_conventional(:,:,:)
  real(8) :: overlap

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=1) error stop 'runtime fixture currently supports exactly one MPI rank'

  dc%i_frag=1
  dc%id_frag=0
  dc%id_tot=rank
  dc%n_frag=1
  dc%isize_frag=1
  dc%isize_tot=nproc
  dc%nstate_tot=2
  dc%icomm_frag=MPI_COMM_WORLD
  dc%icomm_tot=MPI_COMM_WORLD
  n_basis(1,1)=basis_count
  n_mat(1)=basis_count
  h_diag=0d0
  h_diag(1,1,1)=1d0
  h_diag(2,2,1)=2d0
  h_diag(3,3,1)=3d0
  esp_expanded=-huge(1d0)
  esp_conventional=-huge(1d0)

  allocate(coef_expanded(basis_count,basis_count,nspin))
  coef_expanded=0d0
#ifdef OLD_CHEFSI_INTERFACE
  call diag_chefsi(dc,nspin,8,0,4,1d-11,n_basis,n_mat,0,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_expanded,coef_expanded)
#else
  call diag_chefsi(dc,nspin,8,0,4,1d-11,basis_count,n_basis,n_mat,0,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_expanded,coef_expanded)
#endif

  if(any(.not.ieee_is_finite(esp_expanded(:,1)))) &
    error stop 'expanded eigenvalues are not all finite'
  if(maxval(abs(esp_expanded(:,1)-[1d0,2d0,3d0]))>1d-10) &
    error stop 'expanded eigenvalues do not match the diagonal fixture'
  if(any(.not.ieee_is_finite(coef_expanded(:,:,1)))) &
    error stop 'expanded coefficients are not all finite'
  do j=1,basis_count
    if(abs(sum(coef_expanded(:,j,1)**2)-1d0)>1d-10) &
      error stop 'an expanded coefficient column is missing or not normalized'
    do i=1,basis_count
      overlap=sum(coef_expanded(:,i,1)*coef_expanded(:,j,1))
      if(i==j) overlap=overlap-1d0
      if(abs(overlap)>1d-10) error stop 'expanded coefficients are not orthonormal'
    end do
  end do
  if(maxval(abs(coef_expanded(:,3,1)))<0.5d0) &
    error stop 'the state beyond nstate_tot was not produced'

  allocate(coef_conventional(basis_count,dc%nstate_tot,nspin))
  coef_conventional=0d0
#ifdef OLD_CHEFSI_INTERFACE
  call diag_chefsi(dc,nspin,8,0,4,1d-11,n_basis,n_mat,0,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_conventional,coef_conventional)
#else
  call diag_chefsi(dc,nspin,8,0,4,1d-11,dc%nstate_tot,n_basis,n_mat,0,halo_src, &
  & halo_dst,halo_root_src,halo_dvec,h_diag,h_halo,esp_conventional,coef_conventional)
#endif
  if(any(.not.ieee_is_finite(esp_conventional(1:dc%nstate_tot,1)))) &
    error stop 'conventional eigenvalues are not finite'
  if(maxval(abs(esp_conventional(1:dc%nstate_tot,1)- &
    & esp_expanded(1:dc%nstate_tot,1)))>1d-10) &
    error stop 'conventional leading eigenvalues changed'
  do j=1,dc%nstate_tot
    overlap=abs(sum(coef_conventional(:,j,1)*coef_expanded(:,j,1)))
    if(abs(overlap-1d0)>1d-10) error stop 'conventional leading coefficient changed'
  end do

  if(rank==0) write(*,'(a)') 'PASS actual CheFSI expanded and conventional runtime fixture'
  deallocate(coef_conventional,coef_expanded)
  call MPI_Finalize(ierr)
end program test_lcfo_diag_chefsi_requested_count_mpi
