program probe
 use mpi
 use lcfo_dist_rows
 implicit none
 type(s_lcfo_halo) :: rows
 type(s_lcfo_column_halo) :: cols
 integer :: rank,np,ierr,trial,i,j,n,nc
 integer,allocatable :: counts(:),chosen(:),ids(:)
 complex(8),allocatable :: local(:,:),full(:,:),compact(:,:)
 call MPI_Init(ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
 allocate(counts(np));counts=2;counts(1)=0
 allocate(local(counts(rank+1),7))
 do trial=1,4
  if(trial==1)then
   chosen=[sum(counts),1,sum(counts)]
   if(rank==0)then
    ids=[integer::]
   else
    ids=[7,rank,1]
   endif
  elseif(trial==2)then
   chosen=[integer::];ids=[3,1]
  elseif(trial==3)then
   chosen=[1,sum(counts)];ids=[(j,j=1,7)]
  else
   ! Rank 1 owns rows requested by others but receives no columns itself.
   chosen=[1,sum(counts)]
   ids=[6,2]
   if(rank==1)ids=[integer::]
  endif
  call lcfo_halo_init(rows,counts,chosen,MPI_COMM_WORLD)
  call lcfo_column_halo_init(cols,rows,ids,7)
  do n=1,2
   do j=1,7;do i=1,size(local,1)
    local(i,j)=cmplx(100*rank+10*j+i,n*j,8)
   enddo;enddo
   call lcfo_halo_get(rows,local,full)
   call lcfo_column_halo_get(cols,rows,local,compact)
   if(any(shape(compact)/=[size(chosen),size(ids)]))error stop 'shape'
   if(any(compact/=full(:,ids)))error stop 'column transport mismatch'
  enddo
 enddo
 if(rank==0)print *, 'selected column halo PASS',np
 call MPI_Finalize(ierr)
end program
