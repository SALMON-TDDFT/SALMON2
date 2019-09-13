!
!  Copyright 2019 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!
module structure_opt_sub
  implicit none
  ! Global variables like below is not allowed by our policy
  ! These should be moved (do it later...)
  real(8),allocatable :: a_dRion(:), dFion(:)
  real(8),allocatable :: Hess_mat(:,:), Hess_mat_last(:,:)
contains

  !==============================================================initilize
  subroutine structure_opt_ini(natom)
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none
    integer,intent(in) :: natom

    allocate(a_dRion(3*natom),dFion(3*natom))
    allocate(Hess_mat(3*natom,3*natom),Hess_mat_last(3*natom,3*natom))

    a_dRion(:)=0d0
    dFion(:)  =0d0
    Hess_mat(:,:)     =0d0
    Hess_mat_last(:,:)=0d0
    if(comm_is_root(nproc_id_global))then
      write(*,*) "===== Grand State Optimization Start ====="
      write(*,*) "       (Quasi-Newton method using Force only)       "
    end if

  end subroutine structure_opt_ini

  !======================================================convergence check
  subroutine structure_opt_check(natom,iopt,flag_opt_conv,Force)
    use structures
    use salmon_global, only: convrg_opt_fmax,unit_system,flag_opt_atom
    use salmon_parallel, only: nproc_id_global,nproc_group_global
    use salmon_communication, only: comm_is_root,comm_bcast
    implicit none
    integer,intent(in) :: natom,iopt
    real(8),intent(in) :: Force(3,natom)
    logical,intent(inout) :: flag_opt_conv
    real(8),parameter :: a_B=0.529177d0,Ry=13.6058d0
    integer :: iatom,iatom_count
    real(8) :: fabs,fmax,fave

    fmax =0d0
    fave =0d0
    iatom_count=0
    do iatom=1,natom
      if(flag_opt_atom(iatom)=='y') then
        iatom_count = iatom_count+1
        fabs = sum(Force(:,iatom)**2d0)
        fave = fave + fabs
        if(fabs>=fmax) fmax = fabs
      end if
    enddo

    select case(unit_system)
    case('au','a.u.')
      fmax = sqrt(fmax)
      fave = sqrt(fave/iatom_count)
    case('A_eV_fs')
      fmax = sqrt(fmax)*2d0*Ry/a_B
      fave = sqrt(fave/iatom_count)*2d0*Ry/a_B
    end select

    if(comm_is_root(nproc_id_global))then
      write(*,*) " Max-force=",fmax, "  Mean-force=",fave
      write(*,*) "==================================================="
      write(*,*) "Quasi-Newton Optimization Step = ", iopt
      if(fmax<=convrg_opt_fmax) flag_opt_conv=.true.
    end if
    call comm_bcast(flag_opt_conv,nproc_group_global)

  end subroutine structure_opt_check

  !===========================================================optimization
  subroutine structure_opt(natom,iopt,system)
    use structures, only: s_dft_system
    use salmon_global, only: flag_opt_atom
    use salmon_communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(inout) :: system
    integer,intent(in) :: natom,iopt
    !theta_opt=0.0d0:DFP,theta_opt=1.0d0:BFGS in Quasi_Newton method
    real(8), parameter :: alpha=1.0d0,theta_opt=1.0d0
    integer :: ii,ij,icount,iatom, NA3
    real(8) :: const1,const2
    real(8) :: dRion(3,natom)
    real(8) :: force_1d(3*natom),dRion_1d(3*natom),optmat_1d(3*natom)
    real(8) :: optmat1_2d(3*natom,3*natom),optmat2_2d(3*natom,3*natom),optmat3_2d(3*natom,3*natom)

    NA3 = 3*natom

    icount=1
    do iatom=1,natom
    do ii=1,3
       force_1d(icount) = system%Force(ii,iatom)
       icount = icount+1
    end do
    end do

    if(iopt==1)then
      !update Hess_mat
      do ii=1,NA3
      do ij=1,NA3
         if(ii==ij)then
            Hess_mat(ii,ij) = 1d0
         else
            Hess_mat(ii,ij) = 0d0
         end if
         Hess_mat_last(ii,ij) = Hess_mat(ii,ij)
      end do
      end do
    else
      !update dFion
      dFion=-(force_1d-dFion)
      !prepare const and matrix
      call dgemm('n','n',1,1,NA3,1.0d0,a_dRion,1,dFion,NA3,0d0,const1,1)
      call dgemm('n','n',NA3,1,NA3,1d0,Hess_mat,NA3,dFion,NA3,0d0,optmat_1d,NA3)
      call dgemm('n','n',1,1,NA3,1d0,dFion,1,optmat_1d,NA3,0d0,const2,1)
      call dgemm('n','n',NA3,NA3,1,1d0,a_dRion,NA3,a_dRion,1,0d0,optmat1_2d,NA3)
      !update Hess_mat
      Hess_mat = Hess_mat_last + ((const1+theta_opt*const2)/(const1**2d0))*optmat1_2d
      if(theta_opt==0.0d0)then
        !theta_opt=0.0d0:DFP
        call dgemm('n','n',NA3,NA3,1,1d0,optmat_1d,NA3,optmat_1d,1,0d0,optmat2_2d,NA3)
        Hess_mat = Hess_mat-(1d0/const2)*optmat2_2d
      elseif(theta_opt==1.0d0)then
        !theta_opt=1.0d0:BFGS
        call dgemm('n','n',NA3,NA3,1,1d0,optmat_1d,NA3,a_dRion,1,0d0,optmat2_2d,NA3)
        call dgemm('n','n',NA3,NA3,1,1d0,a_dRion,NA3,optmat_1d,1,0d0,optmat3_2d,NA3)
        Hess_mat = Hess_mat-(theta_opt/const1)*(optmat2_2d+optmat3_2d)
      endif
      !update Hess_mat_last
      Hess_mat_last = Hess_mat
    end if

    !update dRion_1d and dRion
    dRion_1d(:) = 0d0
    call dgemm('n','n',NA3,1,NA3,1d0,Hess_mat,NA3,force_1d,NA3,0d0,dRion_1d,NA3)
    do iatom=1,natom
      dRion(1:3,iatom) = dRion_1d((1+3*(iatom-1)):(3+3*(iatom-1)))
    end do

    !update a_dRion,dFion
    a_dRion = alpha * dRion_1d
    dFion   = force_1d

    !update Rion
    do iatom=1,natom
      if(flag_opt_atom(iatom)=='y') then
        system%Rion(1:3,iatom) = system%Rion(1:3,iatom) +alpha*dRion(1:3,iatom)
      end if
    end do
   !call comm_bcast(system%Rion,nproc_group_global) !<-- need?

  end subroutine structure_opt

  !===============================================================finilize
  subroutine structure_opt_fin
    use salmon_parallel, only: nproc_id_global
    use salmon_communication, only: comm_is_root
    implicit none
    deallocate(a_dRion,dFion)
    deallocate(Hess_mat,Hess_mat_last)
    if(comm_is_root(nproc_id_global)) write(*,*) "Optimization Converged"
  end subroutine structure_opt_fin

end module structure_opt_sub
