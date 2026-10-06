program rotation_probe
 use mpi, only: MPI_COMM_SELF,MPI_COMM_WORLD
 use communication, only: comm_init,comm_finalize,comm_get_groupinfo,comm_is_root
 use omp_lib, only: omp_set_num_threads
 use exx_orbitals, only: orbital_rotate
 use exx_distributed_gauge, only: s_exx_gauge,gauge_tiles_rotate
 implicit none
 integer,parameter :: n=3,ng=17
 type(s_exx_gauge) :: gauge
 complex(8) :: matrix(n,n),wave(ng,n),weighted(ng,n),reference(ng,n)
 complex(8),allocatable :: local(:,:),output(:,:),baseline(:,:)
 real(8) :: weights(n),angle
 integer,allocatable :: counts(:)
 integer :: rank,np,first,last,i,j,g,threads,status,mode
 call comm_init
 call comm_get_groupinfo(MPI_COMM_WORLD,rank,np)
 allocate(counts(0:np-1))
 do i=0,np-1
  counts(i)=(i+1)*n/np-i*n/np
 enddo
 first=rank*n/np+1;last=(rank+1)*n/np
 allocate(local(ng,counts(rank)),output(ng,counts(rank)),baseline(ng,counts(rank)))
 gauge%n=n;gauge%comm=MPI_COMM_WORLD
 allocate(gauge%rows(n),gauge%cols(counts(rank)),gauge%matrix(n,counts(rank)))
 do j=1,n
  weights(j)=0.2d0*j
  do i=1,n
   angle=2d0*acos(-1d0)*(i-1)*(j-1)/n
   matrix(i,j)=exp(cmplx(0d0,angle,8))/sqrt(dble(n))
  enddo
  do g=1,ng
   wave(g,j)=cmplx(sin(dble(g*j)),cos(dble(g+j)),8)
  enddo
  weighted(:,j)=wave(:,j)*weights(j)
 enddo
 gauge%rows=[(i,i=1,n)];gauge%cols=[(i,i=first,last)]
 gauge%matrix=matrix(:,first:last);local=wave(:,first:last)
 do mode=1,3
  if(mode==1)reference=matmul(wave,matrix)
  if(mode==2)reference=matmul(weighted,matrix)
  if(mode==3)reference=matmul(wave,conjg(transpose(matrix)))
  do threads=1,4
   call omp_set_num_threads(threads)
   if(mode==2)then
    call gauge_tiles_rotate(gauge,local,MPI_COMM_SELF,MPI_COMM_WORLD,output,status,weights=weights(first:last))
   else
    call gauge_tiles_rotate(gauge,local,MPI_COMM_SELF,MPI_COMM_WORLD,output,status,adjoint=mode==3)
   endif
   call check(status==0,'gauge status')
   call check(all(abs(output-reference(:,first:last))<2d-14),'gauge reference')
   if(threads==1)baseline=output
   call check(all(output==baseline),'gauge OMP exact equality')
   if(mode==2)then
    call orbital_rotate(local,matrix,MPI_COMM_WORLD,counts,first,output,weights(first:last))
   elseif(mode==3)then
    call orbital_rotate(local,conjg(transpose(matrix)),MPI_COMM_WORLD,counts,first,output)
   else
    call orbital_rotate(local,matrix,MPI_COMM_WORLD,counts,first,output)
   endif
   call check(all(abs(output-reference(:,first:last))<2d-14),'orbital reference')
   call check(all(output==baseline),'orbital/gauge equality')
  enddo
 enddo
 if(comm_is_root(rank))print *, 'PASS weighted/forward/adjoint rotations, OMP 1/2/3/4 and empty orbital ranks'
 call comm_finalize
contains
 subroutine check(ok,message)
  use mpi, only: MPI_Abort
  implicit none
  logical,intent(in) :: ok
  character(*),intent(in) :: message
  integer :: err
  if(ok)return
  print *,rank,message
  call MPI_Abort(MPI_COMM_WORLD,1,err)
 end subroutine
end program
