! A symmetric pair of occupied combinations has zero spread gradient but is
! delocalized. Gamma Jacobi and the optional projected seed must both escape it.
program probe
 use mpi
 use exx_spatial
 implicit none
 type(spatial_exx_state) :: seeded,plain,orbital,canonical
 integer :: ierr,np,rank,n(3)=[16,8,8],m(3),g,x,y,z,status
 complex(8),allocatable :: psi(:,:,:)
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 m=[16,8/np,8];allocate(psi(product(m),2,1));psi=0;g=0
 do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
  g=g+1
  if(z/=3.or.y+rank*m(2)/=3)cycle
  if(x==2)psi(g,:,1)=1d0/sqrt(2d0)
  if(x==8)psi(g,:,1)=[1d0,-1d0]/sqrt(2d0)
 enddo;enddo;enddo
 call spatial_exx_refresh(plain,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,20,1d-10,status)
 if(rank==0)print *,'PLAIN status/spread/gradient: ',status,plain%spread,plain%gradient
 if(status/=0.or.plain%spread>1d-8.or.plain%gradient>1d-10)error stop 'Gamma Jacobi failed to escape saddle'
 canonical=plain
 call spatial_exx_canonical_source(canonical,psi,reshape([2d0,.5d0],[2,1]),MPI_COMM_WORLD,status)
 if(status/=0.or.allocated(canonical%gauge).or.allocated(canonical%previous))error stop 'canonical gauge allocation'
 if(maxval(abs(canonical%source(:,1)-psi(:,1,1)))>1d-14.or. &
    maxval(abs(canonical%source(:,2)-.5d0*psi(:,2,1)))>1d-14)error stop 'canonical occupation weighting'
 call spatial_exx_canonical_source(canonical,psi,reshape([0d0,2d0],[2,1]),MPI_COMM_WORLD,status,comm_o=MPI_COMM_SELF)
 if(status/=0.or.any(canonical%source(:,1)/=0d0).or. &
    maxval(abs(canonical%source(:,2)-psi(:,2,1)))>1d-14)error stop 'canonical refresh'
 seeded%seed_localized=.true.
 call spatial_exx_refresh(seeded,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,20,1d-10,status)
 if(status/=0.or.seeded%last_localization_status/=0.or.seeded%spread>1d-8)error stop 'seed failed to localize'
 orbital%seed_localized=.true.
 call spatial_exx_refresh(orbital,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,20,1d-10,status,comm_o=MPI_COMM_SELF)
 if(status/=0.or.orbital%last_localization_status/=0.or.orbital%spread>1d-8)error stop 'orbital path seed failed'
 ! Force a rank-independent transport singularity and defer localization once.
 seeded%previous=0
 call spatial_exx_refresh(seeded,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,0,1d-10,status)
 if(status/=0.or..not.seeded%seed_needed)error stop 'lost reset seed requirement'
 call spatial_exx_refresh(seeded,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,20,1d-10,status)
 if(status/=0.or.seeded%seed_needed.or.seeded%spread>1d-8)error stop 'deferred reseed failed'
 seeded%retain_accepted_gauge=.true.;orbital%retain_accepted_gauge=.true.
 ! Perturb the occupied subspace, then deliberately exhaust minimization.
 do g=1,size(psi,1)
  psi(g,1,1)=psi(g,1,1)+1d-3*cmplx(sin(real(g+rank,8)),cos(real(2*g+rank,8)),8)
  psi(g,2,1)=psi(g,2,1)+1d-3*cmplx(cos(real(3*g+rank,8)),sin(real(g+rank,8)),8)
 enddo
 call spatial_exx_refresh(seeded,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,1,1d-30,status)
 if(status/=0.or.seeded%localization_status==0)error stop 'failed minimization not exercised'
 if(.not.seeded%retained_gauge.or.seeded%last_localization_status/=0)error stop 'accepted gauge lost'
 call spatial_exx_refresh(orbital,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,1,1d-30,status,comm_o=MPI_COMM_SELF)
 if(status/=0.or.orbital%localization_status==0)error stop 'orbital failure not exercised'
 if(.not.orbital%retained_gauge.or.orbital%last_localization_status/=0)error stop 'orbital accepted gauge lost'
 ! SCF/default mode must retain the old retry/fallback contract.
 seeded%retain_accepted_gauge=.false.
 call spatial_exx_refresh(seeded,n,[1d0,1d0,1d0],[np,1],[rank,0],[MPI_COMM_WORLD,MPI_COMM_SELF], &
  MPI_COMM_WORLD,psi,1,1d-30,status)
 if(status/=0.or.seeded%last_localization_status==0.or.seeded%retained_gauge)error stop 'SCF retention enabled'
 if(rank==0)print *,'SPREAD plain/seeded: ',plain%spread,seeded%spread
 call MPI_Finalize(ierr)
end program
