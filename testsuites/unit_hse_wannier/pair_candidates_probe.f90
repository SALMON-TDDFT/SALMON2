program pair_candidates_probe
 use mpi
 use iso_fortran_env, only:int64
 use exx_pair_candidates
 implicit none
 type(exx_pair_catalog) :: catalog
 complex(8),allocatable :: target(:,:),source(:)
 integer,allocatable :: selected(:)
 integer :: ierr,rank,peers,n(3),m(3),lo(3),norb,g,x,y,z,j,k,st,ns,low(3),high(3),trial,world_rank,group,color,split
 integer(int64) :: total,visited
 real(8) :: threshold,tail
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,world_rank,ierr)
 do split=0,1
 color=mod(world_rank,1+split)
 call MPI_Comm_split(MPI_COMM_WORLD,color,world_rank,group,ierr)
 call MPI_Comm_rank(group,rank,ierr)
 call MPI_Comm_size(group,peers,ierr)
 do trial=1,3
  norb=4*2**(trial-1)+color;n=[8*norb,8,8];m=[n(1),8/peers,8];lo=[0,rank*m(2),0]
  allocate(target(product(m),norb),source(product(m)),selected(norb))
  target=0;g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1
   if(y+lo(2)/=2.or.z/=2)cycle
   do j=1,norb
    if(x==8*(j-1)+1)target(g,j)=cmplx(1d0,.25d0*j,8)
   enddo
  enddo;enddo;enddo
  call pair_catalog_build(catalog,n,lo,m,group,target,0d0,st)
  if(st/=0)error stop 'catalog build'
  total=0;visited=0
  do j=1,norb
   source=target(:,j)
   call pair_source_box(n,lo,m,group,source,low,high,st)
   if(st/=0)error stop 'source box'
   call pair_catalog_query(catalog,low,high,0d0,selected,ns,st)
   if(st/=0.or.ns/=1)error stop 'chain candidate count'
   if(selected(1)/=j)error stop 'chain candidate identity'
   total=total+ns;visited=visited+catalog%last_visited
  enddo
  if(total/=norb.or.visited/=norb)error stop 'quadratic candidate generation'
  if(catalog%entries/=norb)error stop 'dense catalogue storage'
  if(rank==0)print *,'CHAIN orbitals/candidates/entries/visited',norb,total,catalog%entries,visited
  ! Nonzero tails may be excluded only with a positive certified threshold.
  tail=1d-12;target=target+cmplx(tail,-tail,8);threshold=1d-10
  call pair_catalog_build(catalog,n,lo,m,group,target,threshold,st)
  if(st/=0)error stop 'tail catalogue'
  do j=1,norb
   low=[8*(j-1),0,0];high=low+[7,7,7]
   call pair_catalog_query(catalog,low,high,threshold,selected,ns,st)
   if(st/=0.or.ns/=1.or.selected(1)/=j)error stop 'tail candidate'
   ! Independent omitted-target bound over the entire queried mesh block.
   do k=1,norb
    if(k==j)cycle
    g=0
    do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
     g=g+1
     if(x>=low(1).and.x<=high(1).and.abs(target(g,k))>threshold)error stop 'unsafe omission'
    enddo;enddo;enddo
   enddo
  enddo
  call pair_catalog_build(catalog,n,lo,m,group,target,0d0,st)
  low=[0,0,0];high=[7,7,7]
  call pair_catalog_query(catalog,low,high,0d0,selected,ns,st)
  if(st/=0.or.ns/=norb)error stop 'zero threshold dropped tiny tail'
  ! Periodic boundary lobes: a loose global box must retain both ends.
  target=0;g=0
  do z=0,m(3)-1;do y=0,m(2)-1;do x=0,m(1)-1
   g=g+1
   if(y+lo(2)==2.and.z==2.and.(x==0.or.x==n(1)-1))target(g,1)=(1d0,1d0)
  enddo;enddo;enddo
  call pair_catalog_build(catalog,n,lo,m,group,target,0d0,st)
  call pair_source_box(n,lo,m,group,target(:,1),low,high,st)
  if(st/=0.or.low(1)/=0.or.high(1)/=n(1)-1)error stop 'periodic bounding box'
  call pair_catalog_query(catalog,low,high,0d0,selected,ns,st)
  if(st/=0.or.ns/=1.or.selected(1)/=1)error stop 'duplicates or lost periodic lobe'
  source=0
  call pair_source_box(n,lo,m,group,source,low,high,st)
  call pair_catalog_query(catalog,low,high,0d0,selected,ns,st)
  if(st/=0.or.ns/=0)error stop 'empty source'
  call pair_catalog_build(catalog,n,lo,m,group,target(:,:0),0d0,st)
  low=0;high=n-1
  call pair_catalog_query(catalog,low,high,0d0,selected(:0),ns,st)
  if(st/=0.or.ns/=0.or.catalog%entries/=0)error stop 'empty target partition'
  deallocate(target,source,selected)
 enddo
 call MPI_Comm_free(group,ierr)
 enddo
 call MPI_Finalize(ierr)
end program
