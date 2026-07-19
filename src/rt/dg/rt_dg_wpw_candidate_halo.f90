module rt_dg_wpw_candidate_halo
  use,intrinsic::ieee_arithmetic,only:ieee_is_finite
  use mpi,only:MPI_Comm_rank,MPI_Comm_size,MPI_Allreduce,MPI_Irecv,MPI_Isend,MPI_Waitall,&
    MPI_INTEGER,MPI_MAX,MPI_DOUBLE_COMPLEX,MPI_STATUSES_IGNORE,MPI_SUCCESS,MPI_Wtime
  implicit none
  private
  integer,parameter,public::wpw_candidate_kind_wp=1,wpw_candidate_kind_pp=2
  type,public::s_dg_wpw_staged_candidate
    integer::kind=0,image_id=0,source_fragment=0
    integer::wp_w=0,wp_p=0,pp_r=0,pp_c=0
    complex(8)::wp_h=(0d0,0d0),wp_s=(0d0,0d0),pp_h=(0d0,0d0),pp_s=(0d0,0d0)
  end type
  type,public::s_dg_wpw_owned_candidates
    integer::epoch=-1
    integer,allocatable::wp_w(:),wp_p(:),pp_r(:),pp_c(:)
    complex(8),allocatable::wp_h(:),wp_s(:),pp_h(:),pp_s(:)
  end type
  type::s_peer_payload
    integer::peer=-1,count=0
    integer,allocatable::ibuf(:)
    complex(8),allocatable::zbuf(:)
  end type
  public::route_dg_wpw_candidate_halo
contains
  subroutine route_dg_wpw_candidate_halo(comm,epoch,n_frag,n_g,support_fragments,staged,owned,info,&
      max_candidates_per_peer)
    integer,intent(in)::comm,epoch,n_frag,n_g,support_fragments(:),max_candidates_per_peer
    type(s_dg_wpw_staged_candidate),intent(in)::staged(:)
    type(s_dg_wpw_owned_candidates),intent(out)::owned
    integer,intent(out)::info
    type(s_peer_payload),allocatable::send(:),recv(:)
    type(s_dg_wpw_staged_candidate),allocatable::work(:)
    integer,allocatable::requests(:),headers_send(:,:),headers_recv(:,:),wp_order(:),pp_order(:)
    integer::rank,nrank,ierr,local_bad,global_bad,npeer,i,j,k,owner,total,nrequest
    integer::global_candidate_bound
    real(8)::route_started

    info=1;local_bad=0;route_started=MPI_Wtime()
    call MPI_Comm_rank(comm,rank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    call MPI_Comm_size(comm,nrank,ierr);if(ierr/=MPI_SUCCESS)local_bad=1
    if(epoch<0.or.n_frag<=0.or.n_g<=0.or.max_candidates_per_peer<=0.or.nrank/=n_frag.or.&
       .not.strictly_increasing(support_fragments).or.&
       .not.any(support_fragments==rank+1))local_bad=1
    do i=1,size(staged)
      if(staged(i)%image_id<=0.or.candidate_owner(staged(i),n_g)<0.or.&
         .not.any(support_fragments==candidate_owner(staged(i),n_g)+1))local_bad=1
      select case(staged(i)%kind)
      case(wpw_candidate_kind_wp)
        if(staged(i)%wp_w<=0.or.staged(i)%wp_p<=0.or.staged(i)%wp_p>n_frag*n_g.or.&
           .not.finite_complex(staged(i)%wp_h).or..not.finite_complex(staged(i)%wp_s))local_bad=1
      case(wpw_candidate_kind_pp)
        if(staged(i)%pp_r<=0.or.staged(i)%pp_c<=0.or.staged(i)%pp_r>n_frag*n_g.or.&
           staged(i)%pp_c>n_frag*n_g.or..not.finite_complex(staged(i)%pp_h).or.&
           .not.finite_complex(staged(i)%pp_s))local_bad=1
      case default;local_bad=1
      end select
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call MPI_Allreduce(max_candidates_per_peer,global_candidate_bound,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS)return

    npeer=size(support_fragments)-1
    allocate(send(npeer),recv(npeer),headers_send(2,npeer),headers_recv(2,npeer))
    k=0
    do i=1,size(support_fragments)
      if(support_fragments(i)==rank+1)cycle
      k=k+1;send(k)%peer=support_fragments(i)-1;recv(k)%peer=send(k)%peer
      send(k)%count=0
      do j=1,size(staged)
        if(candidate_owner(staged(j),n_g)==send(k)%peer)send(k)%count=send(k)%count+1
      enddo
      if(send(k)%count>global_candidate_bound)local_bad=1
      headers_send(:,k)=[epoch,send(k)%count];headers_recv(:,k)=-1
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    allocate(requests(max(1,2*npeer)));nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(headers_recv(:,i),2,MPI_INTEGER,recv(i)%peer,9101,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(headers_send(:,i),2,MPI_INTEGER,send(i)%peer,9101,comm,requests(nrequest),ierr)
      if(ierr/=MPI_SUCCESS)local_bad=1
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS)local_bad=1
    do i=1,npeer
      if(headers_recv(1,i)/=epoch.or.headers_recv(2,i)<0.or.&
         headers_recv(2,i)>global_candidate_bound)local_bad=1
      recv(i)%count=max(0,headers_recv(2,i))
    enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return

    do i=1,npeer
      allocate(send(i)%ibuf(6*send(i)%count),send(i)%zbuf(4*send(i)%count))
      k=0
      do j=1,size(staged)
        if(candidate_owner(staged(j),n_g)/=send(i)%peer)cycle
        k=k+1
        if(staged(j)%kind==wpw_candidate_kind_wp)then
          send(i)%ibuf(6*k-5:6*k)=[staged(j)%kind,staged(j)%image_id,rank+1,&
            staged(j)%wp_w,staged(j)%wp_p,0]
          send(i)%zbuf(4*k-3:4*k)=[staged(j)%wp_h,staged(j)%wp_s,(0d0,0d0),(0d0,0d0)]
        else
          send(i)%ibuf(6*k-5:6*k)=[staged(j)%kind,staged(j)%image_id,rank+1,&
            staged(j)%pp_r,staged(j)%pp_c,0]
          send(i)%zbuf(4*k-3:4*k)=[staged(j)%pp_h,staged(j)%pp_s,(0d0,0d0),(0d0,0d0)]
        endif
      enddo
      allocate(recv(i)%ibuf(6*recv(i)%count),recv(i)%zbuf(4*recv(i)%count))
    enddo
    deallocate(requests);allocate(requests(max(1,4*npeer)));nrequest=0
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Irecv(recv(i)%ibuf,size(recv(i)%ibuf),MPI_INTEGER,recv(i)%peer,9102,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Irecv(recv(i)%zbuf,size(recv(i)%zbuf),MPI_DOUBLE_COMPLEX,recv(i)%peer,9103,comm,requests(nrequest),ierr)
    enddo
    do i=1,npeer
      nrequest=nrequest+1;call MPI_Isend(send(i)%ibuf,size(send(i)%ibuf),MPI_INTEGER,send(i)%peer,9102,comm,requests(nrequest),ierr)
      nrequest=nrequest+1;call MPI_Isend(send(i)%zbuf,size(send(i)%zbuf),MPI_DOUBLE_COMPLEX,send(i)%peer,9103,comm,requests(nrequest),ierr)
    enddo
    if(nrequest>0)call MPI_Waitall(nrequest,requests,MPI_STATUSES_IGNORE,ierr)
    local_bad=merge(0,1,ierr==MPI_SUCCESS)

    total=count([(candidate_owner(staged(i),n_g)==rank,i=1,size(staged))])
    do i=1,npeer;total=total+recv(i)%count;enddo
    allocate(work(total));k=0
    do i=1,size(staged)
      if(candidate_owner(staged(i),n_g)/=rank)cycle
      k=k+1;work(k)=staged(i);work(k)%source_fragment=rank+1
    enddo
    do i=1,npeer;do j=1,recv(i)%count
      k=k+1
      work(k)%kind=recv(i)%ibuf(6*j-5);work(k)%image_id=recv(i)%ibuf(6*j-4)
      work(k)%source_fragment=recv(i)%ibuf(6*j-3)
      if(work(k)%kind==wpw_candidate_kind_wp)then
        work(k)%wp_w=recv(i)%ibuf(6*j-2);work(k)%wp_p=recv(i)%ibuf(6*j-1)
        work(k)%wp_h=recv(i)%zbuf(4*j-3);work(k)%wp_s=recv(i)%zbuf(4*j-2)
      else if(work(k)%kind==wpw_candidate_kind_pp)then
        work(k)%pp_r=recv(i)%ibuf(6*j-2);work(k)%pp_c=recv(i)%ibuf(6*j-1)
        work(k)%pp_h=recv(i)%zbuf(4*j-3);work(k)%pp_s=recv(i)%zbuf(4*j-2)
      else;local_bad=1;endif
      if(work(k)%source_fragment/=recv(i)%peer+1.or.candidate_owner(work(k),n_g)/=rank)local_bad=1
    enddo;enddo
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)return
    call build_candidate_order(work,.true.,wp_order)
    if(has_duplicate_image_record(work,wp_order,.true.))local_bad=1
    call build_candidate_order(work,.false.,pp_order)
    if(has_duplicate_image_record(work,pp_order,.false.))local_bad=1
    call MPI_Allreduce(local_bad,global_bad,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS.or.global_bad/=0)then;info=2;return;endif
    call publish_aggregated(work,wp_order,pp_order,epoch,owned);info=0
    if(rank==0)write(*,'(1x,a,i0,3(a,i0),a,f12.6)')'[DG-WPW-ROUTE] staged_local=',size(staged),&
      ' routed_local=',size(work),' wp_owned=',size(owned%wp_w),' pp_owned=',size(owned%pp_r),&
      ' seconds=',MPI_Wtime()-route_started
  end subroutine

  integer pure function owner_of(column_id,n_g)result(owner)
    integer,intent(in)::column_id,n_g
    owner=(column_id-1)/n_g
  end function
  integer pure function candidate_owner(candidate,n_g)result(owner)
    type(s_dg_wpw_staged_candidate),intent(in)::candidate
    integer,intent(in)::n_g
    select case(candidate%kind)
    case(wpw_candidate_kind_wp);owner=owner_of(candidate%wp_p,n_g)
    case(wpw_candidate_kind_pp);owner=owner_of(candidate%pp_r,n_g)
    case default;owner=-1
    end select
  end function
  elemental logical function finite_complex(value)result(ok)
    complex(8),intent(in)::value
    ok=ieee_is_finite(real(value,8)).and.ieee_is_finite(aimag(value))
  end function
  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function
  subroutine build_candidate_order(a,wp_kind,order)
    type(s_dg_wpw_staged_candidate),intent(in)::a(:)
    logical,intent(in)::wp_kind
    integer,allocatable,intent(out)::order(:)
    integer,allocatable::scratch(:)
    integer::i,left,middle,right,ia,ib,k,width,n

    n=count([(wpw_kind_matches(a(i),wp_kind),i=1,size(a))])
    allocate(order(n),scratch(n));k=0
    do i=1,size(a)
      if(.not.wpw_kind_matches(a(i),wp_kind))cycle
      k=k+1;order(k)=i
    enddo
    width=1
    do while(width<n)
      left=1
      do while(left<=n)
        middle=min(left+width-1,n);right=min(left+2*width-1,n)
        ia=left;ib=middle+1;k=left
        do while(ia<=middle.and.ib<=right)
          if(candidate_le(a(order(ia)),a(order(ib)),wp_kind))then
            scratch(k)=order(ia);ia=ia+1
          else
            scratch(k)=order(ib);ib=ib+1
          endif
          k=k+1
        enddo
        if(ia<=middle)scratch(k:k+middle-ia)=order(ia:middle)
        if(ib<=right)scratch(k:k+right-ib)=order(ib:right)
        left=left+2*width
      enddo
      order=scratch
      if(width>n/2)exit
      width=2*width
    enddo
  end subroutine build_candidate_order
  logical pure function candidate_le(a,b,wp_kind)result(le)
    type(s_dg_wpw_staged_candidate),intent(in)::a,b
    logical,intent(in)::wp_kind
    integer::ar,ac,br,bc
    if(wpw_kind_matches(a,wp_kind).neqv.wpw_kind_matches(b,wp_kind))then
      le=wpw_kind_matches(a,wp_kind);return
    endif
    if(wp_kind)then;ar=a%wp_p;ac=a%wp_w;br=b%wp_p;bc=b%wp_w
    else;ar=a%pp_r;ac=a%pp_c;br=b%pp_r;bc=b%pp_c;endif
    le=ar<br.or.(ar==br.and.(ac<bc.or.(ac==bc.and.(a%source_fragment<b%source_fragment.or.&
      (a%source_fragment==b%source_fragment.and.a%image_id<=b%image_id)))))
  end function
  logical function has_duplicate_image_record(a,order,wp_kind)result(duplicate)
    type(s_dg_wpw_staged_candidate),intent(in)::a(:)
    integer,intent(in)::order(:)
    logical,intent(in)::wp_kind
    integer::i,icurr,iprev
    duplicate=.false.
    do i=2,size(order)
      icurr=order(i);iprev=order(i-1)
      if(a(icurr)%source_fragment/=a(iprev)%source_fragment.or.&
         a(icurr)%image_id/=a(iprev)%image_id)cycle
      if(wp_kind)then
        duplicate=a(icurr)%wp_p==a(iprev)%wp_p.and.a(icurr)%wp_w==a(iprev)%wp_w
      else
        duplicate=a(icurr)%pp_r==a(iprev)%pp_r.and.a(icurr)%pp_c==a(iprev)%pp_c
      endif
      if(duplicate)return
    enddo
  end function
  logical pure function wpw_kind_matches(candidate,wp_kind)result(matches)
    type(s_dg_wpw_staged_candidate),intent(in)::candidate
    logical,intent(in)::wp_kind
    matches=(wp_kind.and.candidate%kind==wpw_candidate_kind_wp).or.&
      (.not.wp_kind.and.candidate%kind==wpw_candidate_kind_pp)
  end function
  subroutine publish_aggregated(work,wp_order,pp_order,epoch,owned)
    type(s_dg_wpw_staged_candidate),intent(in)::work(:)
    integer,intent(in)::wp_order(:),pp_order(:)
    integer,intent(in)::epoch
    type(s_dg_wpw_owned_candidates),intent(out)::owned
    integer::i,j,n,owned_key_row,owned_key_col
    owned%epoch=epoch
    n=0;owned_key_row=-1;owned_key_col=-1
    do i=1,size(wp_order)
      j=wp_order(i)
      if(n==0.or.work(j)%wp_p/=owned_key_row.or.work(j)%wp_w/=owned_key_col)n=n+1
      owned_key_row=work(j)%wp_p;owned_key_col=work(j)%wp_w
    enddo
    allocate(owned%wp_w(n),owned%wp_p(n),owned%wp_h(n),owned%wp_s(n));owned%wp_h=0;owned%wp_s=0;n=0
    do i=1,size(wp_order)
      j=wp_order(i)
      if(n==0)then
        n=n+1;owned%wp_w(n)=work(j)%wp_w;owned%wp_p(n)=work(j)%wp_p
      else if(work(j)%wp_p/=owned%wp_p(n).or.work(j)%wp_w/=owned%wp_w(n))then
        n=n+1;owned%wp_w(n)=work(j)%wp_w;owned%wp_p(n)=work(j)%wp_p
      endif
      owned%wp_h(n)=owned%wp_h(n)+work(j)%wp_h;owned%wp_s(n)=owned%wp_s(n)+work(j)%wp_s
    enddo
    n=0;owned_key_row=-1;owned_key_col=-1
    do i=1,size(pp_order)
      j=pp_order(i)
      if(n==0.or.work(j)%pp_r/=owned_key_row.or.work(j)%pp_c/=owned_key_col)n=n+1
      owned_key_row=work(j)%pp_r;owned_key_col=work(j)%pp_c
    enddo
    allocate(owned%pp_r(n),owned%pp_c(n),owned%pp_h(n),owned%pp_s(n));owned%pp_h=0;owned%pp_s=0;n=0
    do i=1,size(pp_order)
      j=pp_order(i)
      if(n==0)then
        n=n+1;owned%pp_r(n)=work(j)%pp_r;owned%pp_c(n)=work(j)%pp_c
      else if(work(j)%pp_r/=owned%pp_r(n).or.work(j)%pp_c/=owned%pp_c(n))then
        n=n+1;owned%pp_r(n)=work(j)%pp_r;owned%pp_c(n)=work(j)%pp_c
      endif
      owned%pp_h(n)=owned%pp_h(n)+work(j)%pp_h;owned%pp_s(n)=owned%pp_s(n)+work(j)%pp_s
    enddo
  end subroutine publish_aggregated
end module
