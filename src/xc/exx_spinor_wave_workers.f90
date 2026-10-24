! Exact scalar/SOC orbit exchange: k ownership is separate from r-block computation.
module exx_spinor_wave_workers
 use iso_c_binding
 use iso_fortran_env,only:int64
 use mpi,only:MPI_Comm_size,MPI_Allreduce,MPI_Bcast,MPI_Allgather,MPI_Comm_split,MPI_Allgatherv, &
 MPI_Abort,MPI_Gatherv,MPI_Comm_free,MPI_Send,MPI_Recv,MPI_COMM_NULL,MPI_INTEGER,MPI_MAX, &
 MPI_DOUBLE_COMPLEX,MPI_STATUS_SIZE
 use,intrinsic::ieee_arithmetic,only:ieee_is_finite
 use exx_spinor_wave_fft,only:wave_fft_plans,wave_fft_create,wave_fft_execute,wave_fft_destroy
 use exx_k_exchange,only:exx_k_kernel
 use exx_symmetry,only:s_exx_symmetry_map,symmetry_spin_rotation
 implicit none
 private
 include 'fftw3.f03'
 public::apply_wave_workers
contains
 subroutine apply_wave_workers(op,map,es,target,action,starts,counts,rep_start,rank,comm,ierr, &
  payload,pairwork,densitywork,contractwork,scratchbytes,tile_divisor,rprime_group_count)
 implicit none
 type(exx_k_kernel),intent(in)::op
 type(s_exx_symmetry_map),intent(in)::map
 complex(8),intent(in)::es(:,:,:),target(:,:,:)
 complex(8),intent(out)::action(:,:,:)
 integer,intent(in)::starts(:),counts(:),rep_start,rank,comm
 integer,intent(out)::ierr
 integer(int64),intent(out)::payload,pairwork,densitywork,contractwork,scratchbytes
 integer,optional,intent(in)::tile_divisor,rprime_group_count
 integer::nc,np,ng,nk,nr,nt,nsym,b,chunk,divisor,width,widthmax,base,rows,nfirst,lo,nbatch
 integer::s,g,i,j,k,p,r,c,sl,sr,so,st,gi,gj,kg,q,owner_rank,jlocal,opid,d,nd,ndmax,totaldest
 integer::e,bad,globalbad,col,ncols,column_comm,fft_threads
 integer::settings(14),expected_settings(14)
 integer::column_groups,row_groups,row_rank,column_rank,group_col,column_base,col_g,ncols_g, &
 ksource,delta(3),idx(3),cc(3),source_g,global_k,ii
 integer,allocatable::forward(:,:),wrap(:,:,:),owner(:),operation(:),first(:),kperm(:,:),lookup(:),inverse(:)
 integer,allocatable::repstarts(:),repcounts(:),widths(:),dest(:),dr(:),destcounts(:),displs(:),alldest(:)
 integer,allocatable::datacounts(:),datadispls(:),expected_layout(:)
 complex(8),allocatable::gauge(:,:),spin(:,:,:),rowbuf(:,:,:),rowwf(:,:,:,:),colbuf(:,:,:),tbuf(:,:,:)
 complex(8),allocatable::col_owned(:,:,:),target_owned(:,:,:)
 complex(8),allocatable::partial(:,:,:,:),packed(:),gathered(:),row(:,:),contracted(:,:),tm(:,:)
 complex(8),allocatable,target::rho(:,:,:,:)
 complex(8)::coefficient
 complex(8),parameter::one=(1d0,0d0),zero=(0d0,0d0)
 real(8)::length(3),pi,k0(3),v(3),scaled(3),angle
 type(wave_fft_plans)::fft_plans
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
 if(.not.spinor_values_are_finite(target).or..not.spinor_values_are_finite(es))bad=1
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
 do j=1,nr
 i=0
 do p=1,np
 if(j>=repstarts(p).and.j<repstarts(p)+repcounts(p))i=i+1
 enddo
 if(i/=1)bad=1
 enddo
 if(rep_start<1.or.rep_start+size(target,3)>nr+1)then
 bad=1
 else
 if(starts(rank+1)/=map%first(rep_start))bad=1
 if(counts(rank+1)/=map%first(rep_start+size(target,3))-map%first(rep_start))bad=1
 endif
 divisor=16
 if(present(tile_divisor))divisor=tile_divisor
 if(divisor<1)bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 column_groups=1
 if(present(rprime_group_count))column_groups=rprime_group_count
 if(column_groups<1)bad=1
 if(column_groups>0)then
 if(modulo(np,column_groups)/=0)bad=1
 endif
 settings=[op%n,op%mesh,nk,nr,nt,b,divisor,nsym,column_groups,nc];expected_settings=settings
 call MPI_Bcast(expected_settings,size(settings),MPI_INTEGER,0,comm,e)
 if(any(settings/=expected_settings))bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 row_groups=np/column_groups;row_rank=rank/column_groups;column_rank=modulo(rank,column_groups)
 chunk=1+(ng-1)/divisor
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
 ! Batch more orbit rows to reuse column/target transfers, while bounding four-spin rho.
 b=max(1,min(16,1+(nfirst-1)/row_groups))
 ! Refuse if even the minimum single-pair four-spin rho cannot fit the declared cap.
 if(int(nk,int64)*(16*nc*nc)>67108864_int64)bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)return
 b=min(b,int(67108864_int64/(int(nk,int64)*(16*nc*nc))))
 chunk=max(1,min(chunk,int(67108864_int64/(int(b,int64)*nk*nc*nc*16))))
 if(rank==0)print *,'WAVE_BATCH_ROWS_COLUMNS_GROUPS',b,chunk,row_groups,column_groups
 ndmax=b*nsym
 allocate(rowbuf(row_groups*b,widthmax,nc),rowwf(b,widthmax,nk,nc),colbuf(chunk,widthmax,nc),tbuf(chunk,nt,nc))
 allocate(col_owned(chunk,widthmax,nc),target_owned(chunk,nt,nc))
 allocate(rho(b,chunk,nk,nc*nc),partial(ndmax,nt,nc,nr),dest(ndmax),dr(ndmax))
 allocate(destcounts(np),displs(np),datacounts(np),datadispls(np),alldest(np*ndmax))
 allocate(packed(ndmax*nt*nc),gathered(np*ndmax*nt*nc),row(1,chunk),tm(chunk,nt),contracted(1,nt))
 scratchbytes=16_int64*(size(rowbuf,kind=int64)+size(rowwf,kind=int64)+size(colbuf,kind=int64)+ &
 size(tbuf,kind=int64)+size(col_owned,kind=int64)+size(target_owned,kind=int64)+size(rho,kind=int64)+ &
 size(partial,kind=int64)+size(packed,kind=int64)+ &
 size(gathered,kind=int64)+size(row,kind=int64)+size(tm,kind=int64)+size(contracted,kind=int64))
 call wave_fft_create(fft_plans,op%mesh,rho,bad)
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad/=0)goto 900
 call MPI_Comm_split(comm,column_rank,row_rank,column_comm,e)
 do base=1,nfirst,row_groups*b
 lo=base+row_rank*b;rows=max(0,min(b,nfirst-lo+1));nbatch=min(row_groups*b,nfirst-base+1)
 rowwf=zero
 ! Only a small row tile is broadcast; retain only this worker's rows.
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
 do column_base=1,ng,column_groups*chunk
 col=column_base+column_rank*chunk;ncols=max(0,min(chunk,ng-col+1));rho=zero
 ! Stream each full-k source column tile, never all source columns/full WFs.
 do p=0,np-1;do j=1,counts(p+1)
 kg=starts(p+1)+j-1
 do group_col=0,column_groups-1
 col_g=column_base+group_col*chunk;ncols_g=max(0,min(chunk,ng-col_g+1));colbuf=zero
 if(rank==p)then
 do c=1,nc;do g=1,ncols_g
 gi=col_g+g-1;colbuf(g,:width,c)=es(gi+(c-1)*ng,:,j)*op%phase(gi,j)
 enddo;enddo
 endif
 call route_column(colbuf,p,group_col,701)
 if(group_col==column_rank)col_owned=colbuf
 enddo
 ki_block:block
 integer::ki
 ki=inverse(kg)
 do c=1,nc*nc
 sl=1+modulo(c-1,nc);sr=1+(c-1)/nc
 if(rows>0)call zgemm('N','C',rows,ncols,widthmax,one,rowwf(1,1,kg,sl),b, &
 col_owned(1,1,sr),chunk,zero,rho(1,1,ki,c),b)
 densitywork=densitywork+int(rows,int64)*ncols
 enddo
 end block ki_block
 enddo;enddo
 call wave_fft_execute(fft_plans,rho,op%n,op%point,op%shift,op%kernel,first,lo,col,rows,ncols,fft_threads)
 if(base==1.and.column_base==1.and.rank==0)print *,'WAVE_FFT_OMP_THREADS',fft_threads
 pairwork=pairwork+4_int64*rows*ncols
 ! Stream one representative target/operation tile. All workers contract their own destinations.
 do gj=1,nr
 owner_rank=-1
 do p=0,np-1
 if(gj>=repstarts(p+1).and.gj<repstarts(p+1)+repcounts(p+1))owner_rank=p
 enddo
 if(owner_rank<0)call MPI_Abort(comm,91,e)
 global_k=map%first(gj)
 do opid=1,nsym
 do group_col=0,column_groups-1
 col_g=column_base+group_col*chunk;ncols_g=max(0,min(chunk,ng-col_g+1))
 tbuf=zero
 if(rank==owner_rank)then
 jlocal=gj-rep_start+1;ii=global_k-starts(rank+1)+1
 do st=1,nc;do g=1,ncols_g
 gi=forward(col_g+g-1,opid)
 tbuf(g,:,st)=target(gi+(st-1)*ng,:,jlocal)*op%phase(gi,ii)
 enddo;enddo
 endif
 call route_column(tbuf,owner_rank,group_col,702)
 if(group_col==column_rank)target_owned=tbuf
 enddo
 do i=1,nd
 d=dest(i);if(operation(d)/=opid)cycle
 r=dr(i);gi=first(lo+r-1);ksource=kperm(global_k,opid);q=inverse(ksource)
 do c=1,nc*nc
 sl=1+modulo(c-1,nc);sr=1+(c-1)/nc
 do g=1,ncols
 source_g=col+g-1;ii=forward(source_g,opid);delta=wrap(:,source_g,opid)-wrap(:,gi,opid)
 angle=sum(map%full_k(:,global_k)*length*delta)
 row(1,g)=rho(r,g,q,c)*exp(cmplx(0d0,angle,8))*gauge(d,opid)*conjg(gauge(ii,opid))
 enddo
 do so=1,nc;do st=1,nc
 coefficient=spin(so,sl,opid)*conjg(spin(st,sr,opid))
 ! Spin rotations can contain exact structural zeros; no numerical cutoff is applied.
 if(coefficient==zero)cycle
 tm=target_owned(:,:,st)
 call zgemm('N','N',1,nt,ncols,-coefficient/real(nk,8),row,1,tm,chunk,zero,contracted,1)
 angle=sum((map%full_k(:,global_k)-k0)*map%h*op%point(:,d))
 partial(i,:,so,gj)=partial(i,:,so,gj)+contracted(1,:)*exp(cmplx(0d0,-angle,8))
 contractwork=contractwork+int(ncols,int64)*nt
 enddo;enddo
 enddo
 enddo
 enddo
 enddo
 enddo ! columns
 ! Return bounded partial actions only to their k owner, no action replica on helpers.
 do gj=1,nr
 owner_rank=-1
 do p=0,np-1
 if(gj>=repstarts(p+1).and.gj<repstarts(p+1)+repcounts(p+1))owner_rank=p
 enddo
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
 bad=0
 if(.not.spinor_values_are_finite(action))bad=1
 call MPI_Allreduce(bad,globalbad,1,MPI_INTEGER,MPI_MAX,comm,e)
 if(globalbad==0)ierr=0
900 continue
 call wave_fft_destroy(fft_plans)
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
 leader=group
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
 if(rank==leader)payload=payload+int(size(tile),int64)*16*(row_groups-1)
 endif
 end subroutine
 end subroutine
 logical function spinor_values_are_finite(a) result(ok)
 implicit none
 complex(8),intent(in)::a(:,:,:)
 integer::i,j,k
 ok=.true.
 do k=1,size(a,3)
 do j=1,size(a,2)
 do i=1,size(a,1)
  if(.not.ieee_is_finite(real(a(i,j,k))).or..not.ieee_is_finite(aimag(a(i,j,k))))ok=.false.
 enddo
 enddo
 enddo
 end function
end module
