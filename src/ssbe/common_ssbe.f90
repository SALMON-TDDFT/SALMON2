module common_ssbe
    implicit none
contains

# define DX(dt) kdx(kx+(dt)),ky,kz
# define DY(dt) kx,kdy(ky+(dt)),kz
# define DZ(dt) kx,ky,kdz(kz+(dt))    

subroutine grad_k_array_nb1d_double(nb,nk,b_matrix,  &
                                &   array_nb1d,grad_k_array_nb1d)
    use salmon_global, only: num_kgrid
    use initialization_sub, only: set_bN
    implicit none
    integer, intent(in) :: nb
    integer, intent(in) :: nk
    real(8), intent(in) :: b_matrix(3,3)
    real(8), intent(in) :: array_nb1d(1:nb,1:nk)
    real(8), intent(out) :: grad_k_array_nb1d(1:nb,3,1:nk)
    integer :: ik, ib
    integer :: j
    integer :: kx, ky, kz
    integer :: iik
    integer :: kdx(-3:num_kgrid(1)+4)
    integer :: kdy(-3:num_kgrid(2)+4)
    integer :: kdz(-3:num_kgrid(3)+4)
    real(8) :: work(num_kgrid(1),num_kgrid(2),num_kgrid(3))
    real(8) :: grad_k_work(3,num_kgrid(1),num_kgrid(2),num_kgrid(3))
    real(8) :: bnmat(4,4)
    real(8) :: nabt(4,3)
    real(8) :: w

    if (nk /= num_kgrid(1)*num_kgrid(2)*num_kgrid(3)) then
        write(*,*) "ERROR: nk mismatch", nk, num_kgrid, num_kgrid(1)*num_kgrid(2)*num_kgrid(3)
        stop
    end if
    
    do j=-3,num_kgrid(1)+4
        kdx(j) = mod(j+num_kgrid(1)-1,num_kgrid(1))+1
    end do
    do j=-3,num_kgrid(2)+4
        kdy(j) = mod(j+num_kgrid(2)-1,num_kgrid(2))+1
    end do
    do j=-3,num_kgrid(3)+4
        kdz(j) = mod(j+num_kgrid(3)-1,num_kgrid(3))+1
    end do
    
    call set_bN(bnmat)
    do j=1,4
        nabt(j,1) = bnmat(j,4)/(b_matrix(1,1)/num_kgrid(1))
        nabt(j,2) = bnmat(j,4)/(b_matrix(2,2)/num_kgrid(2))
        nabt(j,3) = bnmat(j,4)/(b_matrix(3,3)/num_kgrid(3))
    end do 

    !$omp parallel do default(shared) &
    !$omp private(ib,kx,ky,kz,iik,work,w,grad_k_work)
    do ib=1,nb
        do kz=1,num_kgrid(3)
        do ky=1,num_kgrid(2)
        do kx=1,num_kgrid(1)
            iik=(kz-1)*num_kgrid(2)*num_kgrid(1)+(ky-1)*num_kgrid(1)+kx
            work(kx,ky,kz)=array_nb1d(ib,iik)
        end do
        end do
        end do

        do kz=1,num_kgrid(3)
        do ky=1,num_kgrid(2)
        do kx=1,num_kgrid(1)
            w =  nabt(1,1)*(work(DX(1)) - work(DX(-1))) &
              & +nabt(2,1)*(work(DX(2)) - work(DX(-2))) &
              & +nabt(3,1)*(work(DX(3)) - work(DX(-3))) &
              & +nabt(4,1)*(work(DX(4)) - work(DX(-4)))
            grad_k_work(1,kx,ky,kz) = w

            w =  nabt(1,2)*(work(DY(1)) - work(DY(-1))) &
              & +nabt(2,2)*(work(DY(2)) - work(DY(-2))) &
              & +nabt(3,2)*(work(DY(3)) - work(DY(-3))) &
              & +nabt(4,2)*(work(DY(4)) - work(DY(-4)))
            grad_k_work(2,kx,ky,kz) = w

            w =  nabt(1,3)*(work(DZ(1)) - work(DZ(-1))) &
              & +nabt(2,3)*(work(DZ(2)) - work(DZ(-2))) &
              & +nabt(3,3)*(work(DZ(3)) - work(DZ(-3))) &
              & +nabt(4,3)*(work(DZ(4)) - work(DZ(-4)))
            grad_k_work(3,kx,ky,kz) = w
        end do
        end do
        end do
            
        do kz=1,num_kgrid(3)
        do ky=1,num_kgrid(2)
        do kx=1,num_kgrid(1)
            iik=(kz-1)*num_kgrid(2)*num_kgrid(1)+(ky-1)*num_kgrid(1)+kx
            grad_k_array_nb1d(ib,1,iik) = grad_k_work(1,kx,ky,kz)
            grad_k_array_nb1d(ib,2,iik) = grad_k_work(2,kx,ky,kz)
            grad_k_array_nb1d(ib,3,iik) = grad_k_work(3,kx,ky,kz)
        end do
        end do
        end do
    end do

end subroutine grad_k_array_nb1d_double

subroutine grad_k_array_nb2d_dcomplex(nb,nk,b_matrix,  &
                                 &    array_nb2d,grad_k_array_nb2d)
    use salmon_global, only: num_kgrid
    use initialization_sub, only: set_bN
    implicit none
    integer, intent(in) :: nb
    integer, intent(in) :: nk
    real(8), intent(in) :: b_matrix(3,3)
    complex(8), intent(in) :: array_nb2d(1:nb,1:nb,1:nk)
    complex(8), intent(out) :: grad_k_array_nb2d(1:nb,1:nb,3,1:nk)
    integer :: ik, ib, jb
    integer :: j
    integer :: kx, ky, kz
    integer :: iik
    integer :: kdx(-3:num_kgrid(1)+4)
    integer :: kdy(-3:num_kgrid(2)+4)
    integer :: kdz(-3:num_kgrid(3)+4)
    complex(8) :: work(num_kgrid(1),num_kgrid(2),num_kgrid(3))
    complex(8) :: grad_k_work(3,num_kgrid(1),num_kgrid(2),num_kgrid(3))
    real(8) :: bnmat(4,4)
    real(8) :: nabt(4,3)
    real(8) :: w

    do j=-3,num_kgrid(1)+4
        kdx(j) = mod(j+num_kgrid(1)-1,num_kgrid(1))+1
    end do
    do j=-3,num_kgrid(2)+4
        kdy(j) = mod(j+num_kgrid(2)-1,num_kgrid(2))+1
    end do
    do j=-3,num_kgrid(3)+4
        kdz(j) = mod(j+num_kgrid(3)-1,num_kgrid(3))+1
    end do
    
    call set_bN(bnmat)
    do j=1,4
        nabt(j,1) = bnmat(j,4)/(b_matrix(1,1)/num_kgrid(1))
        nabt(j,2) = bnmat(j,4)/(b_matrix(2,2)/num_kgrid(2))
        nabt(j,3) = bnmat(j,4)/(b_matrix(3,3)/num_kgrid(3))
    end do 

    !$omp parallel do default(shared) &
    !$omp private(ib,jb,kx,ky,kz,iik,work,w,grad_k_work) collapse(2)
    do ib=1,nb
        do jb=1,nb
            do kz=1,num_kgrid(3)
            do ky=1,num_kgrid(2)
            do kx=1,num_kgrid(1)
                iik=(kz-1)*num_kgrid(2)*num_kgrid(1)+(ky-1)*num_kgrid(1)+kx
                work(kx,ky,kz)=array_nb2d(ib,jb,iik)
            end do
            end do
            end do

            do kz=1,num_kgrid(3)
            do ky=1,num_kgrid(2)
            do kx=1,num_kgrid(1)
                w =  nabt(1,1)*(work(DX(1)) - work(DX(-1))) &
                  & +nabt(2,1)*(work(DX(2)) - work(DX(-2))) &
                  & +nabt(3,1)*(work(DX(3)) - work(DX(-3))) &
                  & +nabt(4,1)*(work(DX(4)) - work(DX(-4)))
                grad_k_work(1,kx,ky,kz) = w
    
                w =  nabt(1,2)*(work(DY(1)) - work(DY(-1))) &
                  & +nabt(2,2)*(work(DY(2)) - work(DY(-2))) &
                  & +nabt(3,2)*(work(DY(3)) - work(DY(-3))) &
                  & +nabt(4,2)*(work(DY(4)) - work(DY(-4)))
                grad_k_work(2,kx,ky,kz) = w
    
                w =  nabt(1,3)*(work(DZ(1)) - work(DZ(-1))) &
                  & +nabt(2,3)*(work(DZ(2)) - work(DZ(-2))) &
                  & +nabt(3,3)*(work(DZ(3)) - work(DZ(-3))) &
                  & +nabt(4,3)*(work(DZ(4)) - work(DZ(-4)))
                grad_k_work(3,kx,ky,kz) = w
            end do
            end do
            end do
                
            do kz=1,num_kgrid(3)
            do ky=1,num_kgrid(2)
            do kx=1,num_kgrid(1)
                iik=(kz-1)*num_kgrid(2)*num_kgrid(1)+(ky-1)*num_kgrid(1)+kx
                grad_k_array_nb2d(ib,jb,1,iik) = grad_k_work(1,kx,ky,kz)
                grad_k_array_nb2d(ib,jb,2,iik) = grad_k_work(2,kx,ky,kz)
                grad_k_array_nb2d(ib,jb,3,iik) = grad_k_work(3,kx,ky,kz)
            end do
            end do
            end do
        end do
    end do

end subroutine grad_k_array_nb2d_dcomplex

end module common_ssbe

