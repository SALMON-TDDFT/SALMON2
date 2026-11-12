#include "config.h"
! Exact scalar/SOC orbit exchange: k ownership is separate from r-block computation.
module exx_spinor_wave_workers
 use iso_c_binding
 use iso_fortran_env,only:int64
 !$ use omp_lib,only:omp_get_max_threads,omp_get_thread_num
 use mpi,only:MPI_Comm_size,MPI_Allreduce,MPI_Bcast,MPI_Allgather,MPI_Comm_split,MPI_Allgatherv, &
 MPI_Gatherv,MPI_Comm_free,MPI_Send,MPI_Recv,MPI_COMM_NULL,MPI_INTEGER,MPI_MAX, &
 MPI_DOUBLE_COMPLEX,MPI_STATUS_SIZE,MPI_STATUS_IGNORE
 use exx_spinor_wave_fft,only:wave_fft_plans,wave_fft_create,wave_fft_execute,wave_fft_destroy
 use exx_k_exchange,only:exx_k_kernel
 use exx_symmetry,only:s_exx_symmetry_map,symmetry_spin_rotation
 use exx_worker_layout,only:choose_worker_layout,choose_balanced_strips
 use exx_blas_threads,only:exx_blas_thread_control
 implicit none
 private
 include 'fftw3.f03'
 public::apply_wave_workers,release_worker_fft
 type(wave_fft_plans),save::fft_plans
contains
 subroutine release_worker_fft
 call wave_fft_destroy(fft_plans)
 end subroutine
 subroutine apply_wave_workers(op,map,es,target,action,starts,counts,rep_start,rank,comm,ierr, &
  payload,pairwork,densitywork,contractwork,scratchbytes,tile_divisor,rprime_group_count,balanced_pairs)
 implicit none
 type(exx_k_kernel),intent(in)::op
 type(s_exx_symmetry_map),intent(in)::map
 complex(8),intent(in)::es(:,:,:),target(:,:,:)
 complex(8),intent(out)::action(:,:,:)
 integer,intent(in)::starts(:),counts(:),rep_start,rank,comm
 integer,intent(out)::ierr
 integer(int64),intent(out)::payload,pairwork,densitywork,contractwork,scratchbytes
 integer,optional,intent(in)::tile_divisor,rprime_group_count
 logical,optional,intent(in)::balanced_pairs
 logical::balanced
 character(32)::balanced_env
 integer::env_status,env_length,balanced_flag,row_stride,column_stride,row_capacity
 integer,allocatable::balanced_row(:),balanced_col(:),balanced_count(:),balanced_group(:),leaders(:),group_sizes(:)
 integer::nc,np,ng,nk,nr,nt,nsym,b,chunk,divisor,width,widthmax,base,rows,nfirst,lo,nbatch
 integer::s,g,i,j,k,p,r,c,sl,sr,so,st,gi,gj,kg,q,owner_rank,jlocal,opid,d,nd,ndmax,totaldest
 integer::e,bad,globalbad,col,ncols,column_comm,fft_threads
 integer::contract_threads,contract_team,thread,requested_column_groups
 integer::density_threads
#ifdef HAVE_EXX_OPENBLAS_THREADS
 integer::density_blas_state(3)
#endif
 integer::row_comm,row_group,row_leader,row_members,jfirst,jcount,message_count,k_batch
 integer,allocatable::row_leaders(:)
 integer::settings(15),expected_settings(15)
 integer::column_groups,row_groups,row_rank,column_rank,group_col,column_base,col_g,ncols_g, &
 ksource,delta(3),idx(3),cc(3),source_g,global_k,ii
 integer(int64)::peak_pairs,active_ranks
 integer,allocatable::forward(:,:),wrap(:,:,:),owner(:),operation(:),first(:),kperm(:,:),lookup(:),inverse(:)
 integer,allocatable::repstarts(:),repcounts(:),widths(:),dest(:),dr(:),destcounts(:),displs(:),alldest(:)
 integer,allocatable::datacounts(:),datadispls(:),expected_layout(:)
 integer,allocatable::rep_owner(:)
 complex(8),allocatable::gauge(:,:),spin(:,:,:),rowbuf(:,:,:),rowwf(:,:,:,:),colbuf(:,:,:),tbuf(:,:,:)
 complex(8),allocatable::target_owned(:,:,:,:)
 complex(8),allocatable::row_message(:,:,:)
 complex(8),allocatable::partial(:,:,:,:),packed(:),gathered(:),row(:,:),contracted(:,:)
 complex(8),allocatable,target::rho(:,:,:,:)
 complex(8)::coefficient,output_phase
 complex(8),parameter::one=(1d0,0d0),zero=(0d0,0d0)
 real(8)::length(3),pi,k0(3),v(3),scaled(3),angle
 interface
 subroutine zgemm(ta,tb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
 implicit none
 character,intent(in)::ta,tb
 integer,intent(in)::m,n,k,lda,ldb,ldc
 complex(8),intent(in)::alpha,beta,a(lda,*),b(ldb,*)
 complex(8),intent(inout)::c(ldc,*)
 end subroutine
 end interface
 column_comm=MPI_COMM_NULL
 row_comm=MPI_COMM_NULL
 ierr=1;payload=0;pairwork=0;densitywork=0;contractwork=0;scratchbytes=0;action=zero
 call MPI_Comm_size(comm,np,e)
 bad=0
 if(size(starts)/=np.or.size(counts)/=np.or.rank<0.or.rank>=np)bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 ng=op%ng;nk=op%nk;nr=map%nrep;nt=size(target,2);nsym=map%nsym;b=op%block
 bad=0
 if(min(ng,nk,nr,nt,nsym,b)<1)bad=1
 nc=0
 if(ng>0)nc=size(target,1)/ng
 if(nc/=1.and.nc/=2)bad=1
 if(size(starts)/=np.or.size(counts)/=np.or.size(target,1)/=nc*ng.or.nt<1.or. &
  size(es,1)/=nc*ng.or.size(es,3)/=counts(rank+1).or.any(shape(action)/=shape(target)))bad=1
 if(op%nq/=nk.or.any(op%n/=map%n).or.sum(counts)/=nk.or.starts(1)/=1)bad=1
 do p=2,np
 if(starts(p)/=starts(p-1)+counts(p-1))bad=1
 enddo
 if(any(counts<0).or.any(starts<1))bad=1
 if(allocated(op%phase))then
 if(op%phase_start/=starts(rank+1).or.size(op%phase,2)/=counts(rank+1))bad=1
 else
 bad=1
 endif
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 allocate(expected_layout(2*np));expected_layout=[starts,counts]
 call MPI_Bcast(expected_layout,2*np,MPI_INTEGER,0,comm,e)
 if(any(expected_layout/=[starts,counts]))bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 allocate(repstarts(np),repcounts(np),widths(np))
 call MPI_Allgather(rep_start,1,MPI_INTEGER,repstarts,1,MPI_INTEGER,comm,e)
 i=size(target,3);call MPI_Allgather(i,1,MPI_INTEGER,repcounts,1,MPI_INTEGER,comm,e)
 width=size(es,2);call MPI_Allgather(width,1,MPI_INTEGER,widths,1,MPI_INTEGER,comm,e)
 widthmax=maxval(widths)
 allocate(rep_owner(nr))
 rep_owner=-1
 do p=1,np
 if(repcounts(p)==0)cycle
 i=repstarts(p);j=i+repcounts(p)-1
 if(i<1.or.j>nr.or.j<i)then
 bad=1
 cycle
 endif
 if(any(rep_owner(i:j)/=-1))bad=1
 rep_owner(i:j)=p-1
 enddo
 if(any(rep_owner<0))bad=1
 if(rep_start<1.or.rep_start+size(target,3)>nr+1)then
 bad=1
 else
 if(starts(rank+1)/=map%first(rep_start))bad=1
 if(counts(rank+1)/=map%first(rep_start+size(target,3))-map%first(rep_start))bad=1
 endif
 divisor=0
 if(present(tile_divisor))divisor=tile_divisor
 if(divisor<0)bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 column_groups=0
 if(present(rprime_group_count))column_groups=rprime_group_count
 if(column_groups<0)bad=1
 if(column_groups>0)then
 if(modulo(np,column_groups)/=0)bad=1
 endif
 balanced_env=''
 call get_environment_variable('SALMON_EXX_WORKER_BALANCED_PAIRS',balanced_env,length=env_length,status=env_status)
 balanced=.false.
 if(env_status/=0.and.env_status/=1)bad=1
 if(env_length>len(balanced_env))bad=1
 select case(trim(balanced_env))
 case('','0')
 case('1')
  balanced=.true.
 case default
  bad=1
 end select
 if(present(balanced_pairs))balanced=balanced_pairs
 balanced_flag=merge(1,0,balanced)
 settings=[op%n,op%mesh,nk,nr,nt,b,divisor,nsym,column_groups,nc,balanced_flag];expected_settings=settings
 call MPI_Bcast(expected_settings,size(settings),MPI_INTEGER,0,comm,e)
 if(any(settings/=expected_settings))bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 allocate(forward(ng,nsym),wrap(3,ng,nsym),owner(ng),operation(ng),first(ng),gauge(ng,nsym),spin(nc,nc,nsym))
 length=map%n*map%h;pi=acos(-1d0);k0=map%full_k(:,1)
 do s=1,nsym
 if(nc==2)then
 call symmetry_spin_rotation(map%rotation(:,:,s),spin(:,:,s),bad)
 if(bad/=0)exit
 else
 spin(1,1,s)=one
 endif
 do g=1,ng
 v=matmul(map%rotation(:,:,s),real(op%point(:,g),8)/op%n)+map%translation(:,s)
 cc=nint(v*op%n);wrap(:,g,s)=floor(real(cc,8)/op%n);cc=modulo(cc,op%n)
 forward(g,s)=1+cc(1)+op%n(1)*(cc(2)+op%n(2)*cc(3))
 source_g=map%grid(g,s)
 angle=sum(k0*map%h*real(op%point(:,source_g)-op%point(:,g),8))
 gauge(g,s)=exp(cmplx(0d0,angle,8))
 enddo
 enddo
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 owner=0;operation=0;nfirst=0
 do g=1,ng
 if(owner(g)/=0)cycle
 nfirst=nfirst+1;first(nfirst)=g
 do s=1,nsym
 d=forward(g,s)
 if(owner(d)/=0)cycle
 owner(d)=nfirst;operation(d)=s
 enddo
 enddo
 allocate(kperm(nk,nsym),lookup(nk),inverse(nk));lookup=0
 do i=1,nk
 scaled=(map%full_k(:,i)-k0)*length*op%mesh/(2*pi);idx=modulo(nint(scaled),op%mesh)
 j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));lookup(j)=i
 inverse(op%order(i))=i
 enddo
 do s=1,nsym;do i=1,nk
 v=matmul(transpose(map%rotation(:,:,s)),map%full_k(:,i))
 scaled=(v-k0)*length*op%mesh/(2*pi);idx=modulo(nint(scaled),op%mesh)
 j=1+idx(1)+op%mesh(1)*(idx(2)+op%mesh(2)*idx(3));kperm(i,s)=lookup(j)
 enddo;enddo
 ! Select only after counting actual symmetry-reduced r rows. K ownership is independent.
 allocate(leaders(np),group_sizes(np));leaders=0;group_sizes=0
 if(balanced)then
 requested_column_groups=column_groups
 allocate(balanced_row(np),balanced_col(np),balanced_count(np),balanced_group(np))
 call choose_balanced_strips(nfirst,ng,nk,nc,np,balanced_row,balanced_col,balanced_count,balanced_group, &
  leaders,column_groups,chunk,peak_pairs,active_ranks,bad)
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)then
 balanced=.false.
 column_groups=requested_column_groups
 endif
 endif
 if(balanced)then
 b=1;row_groups=nfirst;column_rank=balanced_group(rank+1);row_rank=balanced_row(rank+1)-1
 do p=1,np
 group_sizes(balanced_group(p)+1)=group_sizes(balanced_group(p)+1)+1
 enddo
 if(rank==0)print *,'WAVE_BALANCED_PAIRS_MIN_MAX_ACTIVE',minval(balanced_count),peak_pairs,active_ranks
 else
 call choose_worker_layout(nfirst,ng,nk,nc,np,column_groups,divisor, &
  row_groups,ii,b,chunk,peak_pairs,active_ranks,bad)
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 column_groups=ii;row_rank=rank/column_groups;column_rank=modulo(rank,column_groups)
 do p=1,column_groups
 leaders(p)=p-1;group_sizes(p)=row_groups
 enddo
 endif
 row_stride=row_groups*b;column_stride=column_groups*chunk;row_capacity=row_groups*b
 if(balanced)then
 row_stride=nfirst;column_stride=ng;row_capacity=nfirst
 endif
 if(rank==0)print *,'WAVE_LAYOUT_R_RPRIME_PEAK_ACTIVE',nfirst,ng,peak_pairs,active_ranks
 if(rank==0)print *,'WAVE_BATCH_ROWS_COLUMNS_GROUPS',b,chunk,row_groups,column_groups
 ndmax=b*nsym
 contract_threads=1
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 !$ contract_threads=min(ndmax,omp_get_max_threads())
#endif
 k_batch=min(16,maxval(counts))
 allocate(rowbuf(row_capacity,widthmax,nc),rowwf(b,widthmax,nk,nc),colbuf(chunk,widthmax,nc*k_batch),tbuf(chunk,nt,nc))
 allocate(target_owned(chunk,nt,nc,nsym))
 allocate(rho(b,chunk,nk,nc*nc),partial(ndmax,nt,nc,nr),dest(ndmax),dr(ndmax))
 allocate(destcounts(np),displs(np),datacounts(np),datadispls(np),alldest(np*ndmax))
 allocate(packed(ndmax*nt*nc),gathered(np*ndmax*nt*nc),row(chunk,contract_threads),contracted(nt,contract_threads))
 scratchbytes=16_int64*(size(rowbuf,kind=int64)+size(rowwf,kind=int64)+size(colbuf,kind=int64)+ &
 size(tbuf,kind=int64)+size(target_owned,kind=int64)+size(rho,kind=int64)+ &
 size(partial,kind=int64)+size(packed,kind=int64)+ &
 size(gathered,kind=int64)+size(row,kind=int64)+size(contracted,kind=int64))
 call wave_fft_create(fft_plans,op%mesh,rho,bad)
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)goto 900
 call MPI_Comm_split(comm,column_rank,row_rank,column_comm,e)
 if(balanced)then
 allocate(row_leaders(row_groups))
 do p=np,1,-1
 row_leaders(balanced_row(p))=p-1
 enddo
 call MPI_Comm_split(comm,row_rank,rank,row_comm,e)
 call MPI_Comm_size(row_comm,row_members,e)
 allocate(row_message(widthmax,nc,k_batch))
 scratchbytes=scratchbytes+16_int64*size(row_message,kind=int64)
 endif
 density_threads=1
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_OPENBLAS_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 !$ density_threads=min(k_batch,omp_get_max_threads())
#endif
#ifdef HAVE_EXX_OPENBLAS_THREADS
 ! OpenBLAS uses process-wide control; set and restore it outside the teams.
 density_blas_state=0
 if(max(density_threads,contract_threads)>1)then
 call exx_blas_thread_control(1,density_blas_state)
 if(density_blas_state(1)==0)then
 density_threads=1
 contract_threads=1
 endif
 endif
#endif
 do base=1,nfirst,row_stride
 lo=base+row_rank*b;rows=max(0,min(b,nfirst-lo+1));nbatch=min(row_groups*b,nfirst-base+1)
 rowwf=zero
 if(balanced)then
 do p=0,np-1
 do jfirst=1,counts(p+1),size(row_message,3)
 jcount=min(size(row_message,3),counts(p+1)-jfirst+1)
 message_count=widthmax*nc*jcount
 do row_group=1,row_groups
 row_leader=row_leaders(row_group)
 if(rank==p)then
 gi=first(base+row_group-1)
 row_message(:,:,:jcount)=zero
 do j=1,jcount;do c=1,nc
 row_message(:width,c,j)=es(gi+(c-1)*ng,:,jfirst+j-1)*op%phase(gi,jfirst+j-1)
 enddo;enddo
 if(p/=row_leader)then
 call MPI_Send(row_message,message_count,MPI_DOUBLE_COMPLEX,row_leader,700,comm,e)
 payload=payload+16_int64*message_count
 endif
 else if(rank==row_leader)then
 call MPI_Recv(row_message,message_count,MPI_DOUBLE_COMPLEX,p,700,comm,MPI_STATUS_IGNORE,e)
 endif
 if(row_rank+1==row_group)then
 call MPI_Bcast(row_message,message_count,MPI_DOUBLE_COMPLEX,0,row_comm,e)
 if(rank==row_leader)payload=payload+16_int64*message_count*(row_members-1)
 do j=1,jcount
 kg=starts(p+1)+jfirst+j-2
 rowwf(1,:,kg,:)=row_message(:,:,j)
 enddo
 endif
 enddo
 enddo;enddo
 else
 do p=0,np-1;do j=1,counts(p+1)
 kg=starts(p+1)+j-1;rowbuf=zero
 if(rank==p)then
 do c=1,nc;do g=1,nbatch
 gi=first(base+g-1)
 rowbuf(g,:width,c)=es(gi+(c-1)*ng,:,j)*op%phase(gi,j)
 enddo;enddo
 endif
 call MPI_Bcast(rowbuf,size(rowbuf),MPI_DOUBLE_COMPLEX,p,comm,e)
 if(rank==p)payload=payload+int(size(rowbuf),int64)*16*(np-1)
 if(rows>0)rowwf(:rows,:,kg,:)=rowbuf(row_rank*b+1:row_rank*b+rows,:,:)
 enddo;enddo
 endif
 nd=0
 do d=1,ng
 if(owner(d)<lo.or.owner(d)>=lo+rows)cycle
 nd=nd+1;dest(nd)=d;dr(nd)=owner(d)-lo+1
 enddo
 call MPI_Allgather(nd,1,MPI_INTEGER,destcounts,1,MPI_INTEGER,comm,e)
 displs(1)=0
 do p=2,np;displs(p)=displs(p-1)+destcounts(p-1);enddo
 totaldest=sum(destcounts)
 call MPI_Allgatherv(dest,nd,MPI_INTEGER,alldest,destcounts,displs,MPI_INTEGER,comm,e)
 datacounts=destcounts*nt*nc;datadispls=displs*nt*nc
 partial=zero
 do column_base=1,ng,column_stride
 col=column_base+column_rank*chunk;ncols=max(0,min(chunk,ng-col+1))
 if(balanced)then
 col=balanced_col(rank+1);ncols=balanced_count(rank+1)
 endif
 rho=zero
 ! Stream each full-k source column tile, never all source columns/full WFs.
 do p=0,np-1;do jfirst=1,counts(p+1),k_batch
 jcount=min(k_batch,counts(p+1)-jfirst+1)
 do group_col=0,column_groups-1
 col_g=column_base+group_col*chunk;ncols_g=max(0,min(chunk,ng-col_g+1))
 if(balanced)then
 col_g=balanced_col(leaders(group_col+1)+1);ncols_g=balanced_count(leaders(group_col+1)+1)
 endif
 if(rank==p)then
 colbuf(:,:,:nc*jcount)=zero
 do j=1,jcount;do c=1,nc;do g=1,ncols_g
 gi=col_g+g-1
 colbuf(g,:width,(j-1)*nc+c)=es(gi+(c-1)*ng,:,jfirst+j-1)*op%phase(gi,jfirst+j-1)
 enddo;enddo;enddo
 endif
 call route_column(colbuf(:,:,:nc*jcount),p,group_col,701)
 if(group_col/=column_rank)cycle
 if(rows==0.or.ncols==0)cycle
 ! Each k point writes to a distinct rho plane. MPI runs on the caller.
 !$omp parallel if(min(jcount,density_threads)>1) num_threads(min(jcount,density_threads)) private(kg,c,sl,sr)
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 block
 integer::density_blas_state(3)
 call exx_blas_thread_control(1,density_blas_state)
#endif
 !$omp do schedule(static)
 do j=1,jcount
 kg=starts(p+1)+jfirst+j-2
 ki_block:block
 integer::ki
 ki=inverse(kg)
 do c=1,nc*nc
 sl=1+modulo(c-1,nc);sr=1+(c-1)/nc
 call zgemm('N','C',rows,ncols,widthmax,one,rowwf(1,1,kg,sl),b, &
 colbuf(1,1,(j-1)*nc+sr),chunk,zero,rho(1,1,ki,c),b)
 enddo
 end block ki_block
 enddo
 !$omp end do nowait
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 call exx_blas_thread_control(0,density_blas_state)
 end block
#endif
 !$omp end parallel
 densitywork=densitywork+int(rows,int64)*ncols*jcount*nc*nc
 enddo
 enddo;enddo
 call wave_fft_execute(fft_plans,rho,op%n,op%point,op%shift,op%kernel,first,lo,col,rows,ncols,fft_threads)
 if(base==1.and.column_base==1.and.rank==0)print *,'WAVE_FFT_OMP_THREADS',fft_threads
 pairwork=pairwork+4_int64*rows*ncols
 ! Route all operations for one representative before threading over destinations.
 do gj=1,nr
 owner_rank=rep_owner(gj)
 global_k=map%first(gj)
 do opid=1,nsym
 do group_col=0,column_groups-1
 col_g=column_base+group_col*chunk;ncols_g=max(0,min(chunk,ng-col_g+1))
 if(balanced)then
 col_g=balanced_col(leaders(group_col+1)+1);ncols_g=balanced_count(leaders(group_col+1)+1)
 endif
 tbuf=zero
 if(rank==owner_rank)then
 jlocal=gj-rep_start+1;ii=global_k-starts(rank+1)+1
 do st=1,nc;do g=1,ncols_g
 gi=forward(col_g+g-1,opid)
 tbuf(g,:,st)=target(gi+(st-1)*ng,:,jlocal)*op%phase(gi,ii)
 enddo;enddo
 endif
 call route_column(tbuf,owner_rank,group_col,702)
 if(group_col==column_rank)target_owned(:,:,:,opid)=tbuf
 enddo
 enddo
 if(nd==0.or.ncols==0)cycle
 contract_team=min(nd,contract_threads)
 !$omp parallel if(contract_team>1) num_threads(contract_team) &
 !$omp private(d,opid,r,gi,ksource,q,c,sl,sr,g,source_g,ii,delta,angle,so,st,coefficient,output_phase,thread) &
 !$omp reduction(+:contractwork)
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 block
 integer::contract_blas_state(3)
 contract_blas_state=0
 if(contract_team>1)call exx_blas_thread_control(1,contract_blas_state)
#endif
 !$omp do schedule(static)
 do i=1,nd
 thread=1
 !$ thread=omp_get_thread_num()+1
 d=dest(i);opid=operation(d)
 r=dr(i);gi=first(lo+r-1);ksource=kperm(global_k,opid);q=inverse(ksource)
 angle=sum((map%full_k(:,global_k)-k0)*map%h*op%point(:,d))
 output_phase=exp(cmplx(0d0,-angle,8))
 do c=1,nc*nc
 sl=1+modulo(c-1,nc);sr=1+(c-1)/nc
 do g=1,ncols
 source_g=col+g-1;ii=forward(source_g,opid);delta=wrap(:,source_g,opid)-wrap(:,gi,opid)
 angle=sum(map%full_k(:,global_k)*length*delta)
 row(g,thread)=rho(r,g,q,c)*exp(cmplx(0d0,angle,8))*gauge(d,opid)*conjg(gauge(ii,opid))
 enddo
 do so=1,nc;do st=1,nc
 coefficient=spin(so,sl,opid)*conjg(spin(st,sr,opid))
 ! Spin rotations can contain exact structural zeros; no numerical cutoff is applied.
 if(coefficient==zero)cycle
 call zgemm('N','N',1,nt,ncols,-coefficient/real(nk,8),row(1,thread),1, &
 target_owned(1,1,st,opid),chunk,zero,contracted(1,thread),1)
 partial(i,:,so,gj)=partial(i,:,so,gj)+contracted(:,thread)*output_phase
 contractwork=contractwork+int(ncols,int64)*nt
 enddo;enddo
 enddo
 enddo
 !$omp end do nowait
#if defined(HAVE_EXX_NVPL_THREADS) || defined(HAVE_EXX_FUJITSU_THREADS)
 call exx_blas_thread_control(0,contract_blas_state)
 end block
#endif
 !$omp end parallel
 enddo
 enddo ! columns
 ! Return bounded partial actions only to their k owner, no action replica on helpers.
 do gj=1,nr
 owner_rank=rep_owner(gj)
 packed(:nd*nt*nc)=reshape(partial(:nd,:,:,gj),[nd*nt*nc])
 call MPI_Gatherv(packed,nd*nt*nc,MPI_DOUBLE_COMPLEX,gathered,datacounts,datadispls, &
 MPI_DOUBLE_COMPLEX,owner_rank,comm,e)
 if(rank==owner_rank)payload=payload+int(sum(datacounts)-datacounts(owner_rank+1),int64)*16
 if(rank==owner_rank)then
 jlocal=gj-rep_start+1
 do p=1,np
 do st=1,nc;do j=1,nt;do i=1,destcounts(p)
 ii=datadispls(p)+(st-1)*destcounts(p)*nt+(j-1)*destcounts(p)+i
 d=alldest(displs(p)+i)
 action(d+(st-1)*ng,j,jlocal)=action(d+(st-1)*ng,j,jlocal)+gathered(ii)
 enddo;enddo;enddo
 enddo
 endif
 enddo
 enddo ! first row batches
#ifdef HAVE_EXX_OPENBLAS_THREADS
 call exx_blas_thread_control(0,density_blas_state)
#endif
 ierr=0
900 continue
 if(row_comm/=MPI_COMM_NULL)call MPI_Comm_free(row_comm,e)
 if(column_comm/=MPI_COMM_NULL)call MPI_Comm_free(column_comm,e)
contains
 subroutine route_column(tile,source,group,tag)
 implicit none
 complex(8),intent(inout)::tile(:,:,:)
 integer,intent(in)::source,group,tag
 integer::leader,stat(MPI_STATUS_SIZE),rc
 if(column_groups==1)then
 call MPI_Bcast(tile,size(tile),MPI_DOUBLE_COMPLEX,source,comm,rc)
 if(rank==source)payload=payload+int(size(tile),int64)*16*(np-1)
 return
 endif
 leader=leaders(group+1)
 ! One copy crosses to the selected column group; other column groups receive nothing.
 if(source/=leader)then
 if(rank==source)then
 call MPI_Send(tile,size(tile),MPI_DOUBLE_COMPLEX,leader,tag,comm,rc)
 payload=payload+int(size(tile),int64)*16
 else if(rank==leader)then
 call MPI_Recv(tile,size(tile),MPI_DOUBLE_COMPLEX,source,tag,comm,stat,rc)
 endif
 endif
 if(column_rank==group)then
 call MPI_Bcast(tile,size(tile),MPI_DOUBLE_COMPLEX,0,column_comm,rc)
 if(rank==leader)payload=payload+int(size(tile),int64)*16*(group_sizes(group+1)-1)
 endif
 end subroutine
 end subroutine
end module
