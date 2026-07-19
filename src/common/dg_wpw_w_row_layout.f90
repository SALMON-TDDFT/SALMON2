module dg_wpw_w_row_layout
  implicit none
  private
  public::build_dg_wpw_w_row_layout
contains
  subroutine build_dg_wpw_w_row_layout(owned_fragment,support_fragments,n_basis,index_basis,&
      owned_ids,support_ids,info)
    integer,intent(in)::owned_fragment,support_fragments(:),n_basis(:),index_basis(:,:)
    integer,allocatable,intent(out)::owned_ids(:),support_ids(:)
    integer,intent(out)::info
    integer::i,j,k,n,key
    info=1
    if(owned_fragment<1.or.owned_fragment>size(n_basis).or.size(index_basis,2)/=size(n_basis).or.&
       .not.strictly_increasing(support_fragments).or.any(support_fragments<1).or.&
       any(support_fragments>size(n_basis)).or..not.any(support_fragments==owned_fragment).or.&
       any(n_basis<0).or.any(n_basis>size(index_basis,1)).or.n_basis(owned_fragment)<=0)return
    allocate(owned_ids(n_basis(owned_fragment)))
    owned_ids=index_basis(1:n_basis(owned_fragment),owned_fragment)
    n=sum(n_basis(support_fragments));if(n<=0.or..not.strictly_increasing(owned_ids))then
      deallocate(owned_ids);return
    endif
    allocate(support_ids(n));k=0
    do i=1,size(support_fragments)
      do j=1,n_basis(support_fragments(i))
        k=k+1;support_ids(k)=index_basis(j,support_fragments(i))
      enddo
    enddo
    do i=2,n
      key=support_ids(i);j=i-1
      do while(j>=1)
        if(support_ids(j)<=key)exit
        support_ids(j+1)=support_ids(j);j=j-1
      enddo
      support_ids(j+1)=key
    enddo
    if(.not.strictly_increasing(support_ids).or.any(support_ids<=0))then
      deallocate(owned_ids,support_ids);return
    endif
    info=0
  end subroutine

  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function
end module dg_wpw_w_row_layout
