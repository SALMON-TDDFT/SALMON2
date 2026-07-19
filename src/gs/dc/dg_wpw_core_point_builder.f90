module dg_wpw_core_point_builder
  use dg_wpw_windows,only:evaluate_dg_wpw_normalized_windows
  use rt_dg_wpw_point_evaluator,only:evaluate_windowed_kg_point
  implicit none
  private
  public::evaluate_dg_wpw_core_point
contains
  subroutine evaluate_dg_wpw_core_point(core_lo,core_hi,fragment_ids,total_grid,hgs,buffer,width,&
      gvec,grid,omega_cell,column_ids,chi,grad_chi,p,grad_p,info)
    integer,intent(in)::core_lo(:,:),core_hi(:,:),fragment_ids(:),total_grid(3),buffer(3),width(3),grid(3)
    real(8),intent(in)::hgs(3),gvec(:,:),omega_cell
    integer,intent(out)::column_ids(:),info
    real(8),intent(out)::chi(:),grad_chi(:,:)
    complex(8),intent(out)::p(:),grad_p(:,:)
    integer::k,g,index,point_info,nk,ng
    real(8)::position(3)

    column_ids=0;chi=0d0;grad_chi=0d0;p=(0d0,0d0);grad_p=(0d0,0d0);info=1
    nk=size(fragment_ids);ng=size(gvec,2)
    if(nk<=0.or.ng<=0.or.size(gvec,1)/=3.or..not.strictly_increasing(fragment_ids).or.&
       any(shape(core_lo)/=[3,nk]).or.any(shape(core_hi)/=[3,nk]).or.size(chi)/=nk.or.&
       any(shape(grad_chi)/=[3,nk]).or.size(column_ids)/=nk*ng.or.size(p)/=nk*ng.or.&
       any(shape(grad_p)/=[3,nk*ng]))return
    call evaluate_dg_wpw_normalized_windows(core_lo,core_hi,total_grid,hgs,buffer,width,grid,chi,grad_chi,info)
    if(info/=0)return
    position=dble(grid)*hgs
    do k=1,nk
      do g=1,ng
        index=(k-1)*ng+g
        column_ids(index)=(fragment_ids(k)-1)*ng+g
        call evaluate_windowed_kg_point(chi(k),grad_chi(:,k),gvec(:,g),position,omega_cell,&
          p(index),grad_p(:,index),point_info)
        if(point_info/=0)then;info=2;column_ids=0;p=(0d0,0d0);grad_p=(0d0,0d0);return;endif
      enddo
    enddo
    info=0
  end subroutine

  logical function strictly_increasing(ids)result(ok)
    integer,intent(in)::ids(:);integer::i
    ok=size(ids)>0
    do i=2,size(ids);if(ids(i)<=ids(i-1))then;ok=.false.;return;endif;enddo
  end function
end module dg_wpw_core_point_builder
