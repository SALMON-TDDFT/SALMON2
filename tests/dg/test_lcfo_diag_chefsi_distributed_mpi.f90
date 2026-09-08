program test_lcfo_diag_chefsi_distributed_mpi
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use mpi
  use structures, only: s_dcdft
  use lcfo_diag_chefsi, only: diag_chefsi
  implicit none
  integer, parameter :: max_basis=40,nfrag=2,nspin=1,nrequested=66,nmatrix=80
  type(s_dcdft) :: dc
  integer :: ierr,rank,nproc,i,j
  integer :: n_basis(nfrag,nspin),n_mat(nspin)
  integer :: halo_src(1),halo_dst(1),halo_root_src(1),halo_dvec(3,1)
  real(8) :: h_diag(max_basis,max_basis,nspin)
  real(8) :: h_halo(max_basis,max_basis,nspin,1)
  real(8) :: esp_expanded(nmatrix,nspin),esp_conventional(nmatrix,nspin)
  real(8) :: expected(nmatrix),coupling,delta
  real(8), allocatable :: gram_local(:,:),gram_global(:,:),cross_local(:),cross_global(:)
  real(8), allocatable :: coef_expanded(:,:,:),coef_conventional(:,:,:)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc/=nfrag) error stop 'distributed fixture requires two MPI ranks'
  dc%i_frag=rank+1
  dc%id_frag=0
  dc%id_tot=rank
  dc%n_frag=nfrag
  dc%isize_frag=1
  dc%isize_tot=nproc
  dc%nstate_tot=64
  dc%icomm_frag=MPI_COMM_SELF
  dc%icomm_tot=MPI_COMM_WORLD
  n_basis(:,1)=max_basis
  n_mat(1)=nmatrix
  h_diag=0d0
  coupling=0.2d0
  if(rank==0) then
    do i=1,max_basis
      h_diag(i,i,1)=real(i,8)
    end do
    halo_src(1)=2
    halo_dst(1)=2
    halo_root_src(1)=1
    halo_dvec(:,1)=[1,0,0]
  else
    do i=1,max_basis
      h_diag(i,i,1)=real(i+max_basis,8)
    end do
    halo_src(1)=1
    halo_dst(1)=1
    halo_root_src(1)=0
    halo_dvec(:,1)=[-1,0,0]
  end if
  h_halo=0d0
  do i=1,max_basis
    h_halo(i,i,1,1)=coupling
  end do
  delta=sqrt((0.5d0*real(max_basis,8))**2+coupling**2)
  do i=1,max_basis
    expected(i)=real(i,8)+0.5d0*real(max_basis,8)-delta
    expected(max_basis+i)=real(i,8)+0.5d0*real(max_basis,8)+delta
  end do
  esp_expanded=-huge(1d0)
  esp_conventional=-huge(1d0)
  allocate(coef_expanded(max_basis,nrequested,nspin))
  coef_expanded=0d0
  call diag_chefsi(dc,nspin,20,8,80,1d-10,nrequested,n_basis,n_mat,1, &
    & halo_src,halo_dst,halo_root_src,halo_dvec,h_diag,h_halo, &
    & esp_expanded,coef_expanded)
  if(any(.not.ieee_is_finite(esp_expanded(1:nrequested,1)))) &
    error stop 'distributed expanded eigenvalues are not finite'
  if(maxval(abs(esp_expanded(1:nrequested,1)-expected(1:nrequested)))>1d-8) &
    error stop 'distributed expanded eigenvalues differ from dense reference'
  if(any(.not.ieee_is_finite(coef_expanded(:,:,1)))) &
    error stop 'distributed expanded coefficients are not finite'
  allocate(gram_local(nrequested,nrequested),gram_global(nrequested,nrequested))
  gram_local=matmul(transpose(coef_expanded(:,:,1)),coef_expanded(:,:,1))
  call MPI_Allreduce(gram_local,gram_global,size(gram_local),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  do i=1,nrequested
    gram_global(i,i)=gram_global(i,i)-1d0
  end do
  if(maxval(abs(gram_global))>1d-8) error stop 'distributed expanded coefficients are not orthonormal'

  allocate(coef_conventional(max_basis,dc%nstate_tot,nspin))
  coef_conventional=0d0
  call diag_chefsi(dc,nspin,20,8,80,1d-10,dc%nstate_tot,n_basis,n_mat,1, &
    & halo_src,halo_dst,halo_root_src,halo_dvec,h_diag,h_halo, &
    & esp_conventional,coef_conventional)
  if(maxval(abs(esp_conventional(1:dc%nstate_tot,1)-expected(1:dc%nstate_tot)))>1d-8) &
    error stop 'distributed conventional eigenvalues differ from dense reference'
  allocate(cross_local(dc%nstate_tot),cross_global(dc%nstate_tot))
  do j=1,dc%nstate_tot
    cross_local(j)=sum(coef_conventional(:,j,1)*coef_expanded(:,j,1))
  end do
  call MPI_Allreduce(cross_local,cross_global,size(cross_local),MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(abs(cross_global)-1d0))>1d-7) error stop 'distributed conventional leading vector changed'
  if(rank==0) write(*,'(a)') 'PASS actual distributed CheFSI expanded and conventional fixture'
  deallocate(cross_global,cross_local,gram_global,gram_local)
  deallocate(coef_conventional,coef_expanded)
  call MPI_Finalize(ierr)
end program test_lcfo_diag_chefsi_distributed_mpi
