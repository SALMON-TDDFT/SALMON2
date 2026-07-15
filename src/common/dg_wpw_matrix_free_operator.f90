module dg_wpw_matrix_free_operator
 implicit none
 private
 public::apply_h_batch,apply_s_batch,global_gram_batch
 abstract interface
  subroutine apply_h_batch(context,xw_owned,xp_owned,yw_owned,yp_owned,info)
   class(*), intent(inout) :: context
   complex(8),intent(in)::xw_owned(:,:),xp_owned(:,:)
   complex(8),intent(out)::yw_owned(:,:),yp_owned(:,:)
   integer,intent(out)::info
  end subroutine
  subroutine apply_s_batch(context,xw_owned,xp_owned,yw_owned,yp_owned,info)
   class(*), intent(inout) :: context
   complex(8),intent(in)::xw_owned(:,:),xp_owned(:,:)
   complex(8),intent(out)::yw_owned(:,:),yp_owned(:,:)
   integer,intent(out)::info
  end subroutine
  subroutine global_gram_batch(x,y,n_local_rows,nx,ny,g,info)
   integer,intent(in)::n_local_rows,nx,ny;complex(8),intent(in)::x(n_local_rows,nx),y(n_local_rows,ny)
   complex(8),intent(out)::g(nx,ny);integer,intent(out)::info
  end subroutine
 end interface
end module
