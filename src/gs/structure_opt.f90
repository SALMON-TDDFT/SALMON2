!
!  Copyright 2019-2020 SALMON developers
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
contains

  !==============================================================initilize
  subroutine  initialization_opt(Miopt,opt,system,flag_opt_conv,nopt_max,ofl)
    use structures, only: s_opt,s_dft_system,s_ofile
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    use salmon_global, only: natom,nopt,yn_restart
    use write_sub, only: write_xyz
    use checkpoint_restart_sub, only: restart_opt
    implicit none
    type(s_opt) :: opt
    type(s_dft_system) :: system
    integer :: nopt_max,Miopt,NA3
    logical :: flag_opt_conv
    type(s_ofile) :: ofl
    character(100) :: comment_line

    flag_opt_conv = .false.
    nopt_max = nopt

    NA3 = 3*natom
    allocate( opt%a_dRion(NA3), opt%dFion(NA3) )
    allocate( opt%Hess_mat(NA3,NA3), opt%Hess_mat_last(NA3,NA3) )

    if(yn_restart == 'y') then
       call restart_opt(Miopt,opt)
    else
       opt%a_dRion(:) = 0d0
       opt%dFion(:)   = 0d0
       opt%Hess_mat(:,:)      = 0d0
       opt%Hess_mat_last(:,:) = 0d0
    endif

    if(comm_is_root(nproc_id_global))then
      write(*,*) "===== Grand State Optimization Start ====="
      write(*,*) "       (Quasi-Newton method using Force only)       "
    end if

    write(comment_line,10) 0
    call write_xyz(comment_line,"new","r  ",system,ofl)
10  format("#opt iteration step=",i5)

  end subroutine  initialization_opt

  !======================================================convergence check
  subroutine structure_opt_check(iopt,flag_opt_conv,Force)
    use structures
    use salmon_global, only: natom,convrg_opt_fmax,unit_system,flag_opt_atom
    use inputoutput, only: au_length_aa, au_energy_ev
    use parallelization, only: nproc_id_global,nproc_group_global
    use communication, only: comm_is_root,comm_bcast
    implicit none
    integer,intent(in) :: iopt
    real(8),intent(in) :: Force(3,natom)
    logical,intent(inout) :: flag_opt_conv
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
      fmax = sqrt(fmax)*au_energy_ev/au_length_aa
      fave = sqrt(fave/iatom_count)*au_energy_ev/au_length_aa
    end select

    if(comm_is_root(nproc_id_global))then
      write(*,*) " Max-force=",fmax, "  Mean-force=",fave
      write(*,*) "==================================================="
      write(*,*) "Quasi-Newton Optimization Step = ", iopt
      if( fmax <= convrg_opt_fmax ) flag_opt_conv=.true.
    end if
    call comm_bcast(flag_opt_conv,nproc_group_global)

  end subroutine structure_opt_check

  !===========================================================optimization
  subroutine structure_opt(opt,iopt,system)
    use structures, only: s_dft_system, s_opt
    use salmon_global, only: natom, flag_opt_atom, max_step_len_adjust, &
                             method_opt, step_steep, step_fire, Kion
    use const, only: umass
    use communication, only: comm_bcast
    implicit none
    type(s_dft_system),intent(inout) :: system
    type(s_opt) :: opt
    integer,intent(in) :: iopt
    !theta_opt=0.0d0:DFP,theta_opt=1.0d0:BFGS in Quasi_Newton method
    real(8), parameter :: theta_opt=1.0d0  !alpha=1.0d0 -- this is now from input
    real(8) :: alpha
    integer :: ii,ij,jj,icount,iatom, NA3,ixyz
    real(8) :: const1,const2,rtmp
    real(8) :: dRion(3,natom), dRabs(natom), dRabs_max
    real(8) :: force_1d(3*natom),dRion_1d(3*natom),optmat_1d(3*natom)
    real(8) :: optmat1_2d(3*natom,3*natom),optmat2_2d(3*natom,3*natom),optmat3_2d(3*natom,3*natom)
    real(8) :: f_inc,f_dec,alpha_start,f_alpha,dt_max
    integer :: n_min
    real(8) :: p_value,f_norm,v_norm
    real(8) :: mass_au

    NA3 = 3*natom

    icount=1
    do iatom=1,natom
    do ixyz=1,3
       force_1d(icount) = system%Force(ixyz,iatom)
       icount = icount+1
    end do
    end do
    
    select case(method_opt)
    case('bfgs')

      alpha = 1.0d0

      if(iopt==1)then
        !update opt%Hess_mat
        do ii=1,NA3
        do ij=1,NA3
           if(ii==ij)then
              opt%Hess_mat(ii,ij) = 1d0
           else
              opt%Hess_mat(ii,ij) = 0d0
           end if
           opt%Hess_mat_last(ii,ij) = opt%Hess_mat(ii,ij)
        end do
        end do
      else
        !update opt%dFion
        opt%dFion=-(force_1d-opt%dFion)
        !prepare const and matrix
        call dgemm('n','n',1,1,NA3,1.0d0,opt%a_dRion,1,opt%dFion,NA3,0d0,const1,1)
        call dgemm('n','n',NA3,1,NA3,1d0,opt%Hess_mat,NA3,opt%dFion,NA3,0d0,optmat_1d,NA3)
        call dgemm('n','n',1,1,NA3,1d0,opt%dFion,1,optmat_1d,NA3,0d0,const2,1)
        call dgemm('n','n',NA3,NA3,1,1d0,opt%a_dRion,NA3,opt%a_dRion,1,0d0,optmat1_2d,NA3)
        !update opt%Hess_mat
        rtmp = (const1+theta_opt*const2)/(const1**2d0)
        !$omp parallel do collapse(2) private(ii,jj)
        do ii=1,NA3
        do jj=1,NA3
           opt%Hess_mat(ii,jj) = opt%Hess_mat_last(ii,jj) + rtmp * optmat1_2d(ii,jj)
        enddo
        enddo
        !$omp end parallel do
        if(theta_opt==0.0d0)then
          !theta_opt=0.0d0:DFP
          call dgemm('n','n',NA3,NA3,1,1d0,optmat_1d,NA3,optmat_1d,1,0d0,optmat2_2d,NA3)
          !$omp parallel do collapse(2) private(ii,jj)
          do ii=1,NA3
          do jj=1,NA3
             opt%Hess_mat(ii,jj) = opt%Hess_mat(ii,jj)-(1d0/const2)*optmat2_2d(ii,jj)
          enddo
          enddo
          !$omp end parallel do
        elseif(theta_opt==1.0d0)then
          !theta_opt=1.0d0:BFGS
          call dgemm('n','n',NA3,NA3,1,1d0,optmat_1d,NA3,opt%a_dRion,1,0d0,optmat2_2d,NA3)
          call dgemm('n','n',NA3,NA3,1,1d0,opt%a_dRion,NA3,optmat_1d,1,0d0,optmat3_2d,NA3)
          rtmp = theta_opt/const1
          !$omp parallel do collapse(2) private(ii,jj)
          do ii=1,NA3
          do jj=1,NA3
             opt%Hess_mat(ii,jj) = opt%Hess_mat(ii,jj)- rtmp *(optmat2_2d(ii,jj)+optmat3_2d(ii,jj))
          enddo
          enddo
          !$omp end parallel do
        endif
        !update opt%Hess_mat_last
        !$omp parallel do collapse(2) private(ii,jj)
        do ii=1,NA3
        do jj=1,NA3
           opt%Hess_mat_last(ii,jj) = opt%Hess_mat(ii,jj)
        enddo
        enddo
        !$omp end parallel do
      end if
  
      !update dRion_1d and dRion
      dRion_1d(:) = 0d0
      call dgemm('n','n',NA3,1,NA3,1d0,opt%Hess_mat,NA3,force_1d,NA3,0d0,dRion_1d,NA3)
      !$omp parallel do collapse(2) private(iatom,ixyz)
      do iatom=1,natom
      do ixyz=1,3
         dRion(ixyz,iatom) = dRion_1d(ixyz+3*(iatom-1))       
  !      dRion(1:3,iatom) = dRion_1d((1+3*(iatom-1)):(3+3*(iatom-1)))
      end do
      end do
      !$omp end parallel do
  
      !adjust alpha if required from input
      if(max_step_len_adjust .gt. 0d0) then
         !$omp parallel do private(iatom)
         do iatom=1,natom
            if(flag_opt_atom(iatom)=='y') then
                 dRabs(iatom) = sqrt( sum(dRion(:,iatom)**2d0) )
            else
                 dRabs(iatom) = 0d0
            endif
         enddo
         !$omp end parallel do
         dRabs_max = maxval(dRabs(:))
         alpha = max_step_len_adjust / dRabs_max
      endif
  
      !update opt%a_dRion,opt%dFion
      !$omp parallel do private(ii)
      do ii=1,NA3
         opt%a_dRion(ii) = alpha * dRion_1d(ii)
         opt%dFion(ii)   = force_1d(ii)
      enddo
      !$omp end parallel do

      !update Rion
      !$omp parallel do private(iatom)
      do iatom=1,natom
         if(flag_opt_atom(iatom)=='y') &
            system%Rion(1:3,iatom) = system%Rion(1:3,iatom) +alpha*dRion(1:3,iatom)
      end do
      !$omp end parallel do

    case('steep')
      
      alpha = step_steep

      !update dRion_1d and dRion
      !$omp parallel do private(ii)
      do ii=1,NA3
        dRion_1d(ii) = force_1d(ii)
      enddo
      !$omp end parallel do
      !$omp parallel do collapse(2) private(iatom,ixyz)
      do iatom=1,natom
      do ixyz=1,3
         dRion(ixyz,iatom) = dRion_1d(ixyz+3*(iatom-1))
      end do
      end do
      !$omp end parallel do

      !adjust alpha if required from input
      if(max_step_len_adjust .gt. 0d0) then
         !$omp parallel do private(iatom)
         do iatom=1,natom
            if(flag_opt_atom(iatom)=='y') then
                 dRabs(iatom) = sqrt( sum(dRion(:,iatom)**2d0) )
            else
                 dRabs(iatom) = 0d0
            endif
         enddo
         !$omp end parallel do
         dRabs_max = maxval(dRabs(:))
         alpha = max_step_len_adjust / dRabs_max * 0.1d0
      endif

      !update opt%a_dRion,opt%dFion
      !$omp parallel do private(ii)
      do ii=1,NA3
         opt%a_dRion(ii) = alpha * 10.d0 * dRion_1d(ii)
         opt%dFion(ii)   = force_1d(ii)
      enddo
      !$omp end parallel do
  
      !update Rion
      !$omp parallel do private(iatom)
      do iatom=1,natom
         if(flag_opt_atom(iatom)=='y') &
            system%Rion(1:3,iatom) = system%Rion(1:3,iatom) +alpha*10.d0*dRion(1:3,iatom)
      end do
      !$omp end parallel do

    case('fire')
      f_inc = 1.1d0
      f_dec = 0.5d0
      alpha_start = 0.1d0
      f_alpha = 0.99d0
      n_min = 5
      dt_max = 1.d0/0.024188843d0
  
      if(iopt==1)then
        opt%v_fire(:) = 0.d0
        opt%alpha_fire = alpha_start
        opt%step_fire = step_fire
        opt%p_times = 0
      end if

      p_value = 0.d0
      !$omp parallel do reduction(+:p_value) private(ii)
      do ii=1,NA3
        p_value = p_value + force_1d(ii)*opt%v_fire(ii)
      end do
      !$omp end parallel do
  
      f_norm = 0.d0
      !$omp parallel do reduction(+:f_norm) private(ii)
      do ii=1,NA3
        f_norm = f_norm + force_1d(ii)**2
      end do
      !$omp end parallel do
      f_norm = sqrt(f_norm)
  
      v_norm = 0.d0
      !$omp parallel do reduction(+:v_norm) private(ii)
      do ii=1,NA3
        v_norm = v_norm + opt%v_fire(ii)**2
      end do
      !$omp end parallel do
      v_norm = sqrt(v_norm)

      !$omp parallel do private(ii)
      do ii=1,NA3
        opt%v_fire(ii) = (1.d0-opt%alpha_fire)*opt%v_fire(ii) &
                         + opt%alpha_fire*force_1d(ii)/f_norm*v_norm
      end do
      !$omp end parallel do
  
      if(p_value > 0.d0)then
        opt%p_times = opt%p_times + 1
        if(opt%p_times > n_min)then
          opt%step_fire = min(opt%step_fire*f_inc,dt_max)
          opt%alpha_fire = f_alpha*opt%alpha_fire
        end if
      else
        opt%p_times = 0
        opt%step_fire = f_dec * opt%step_fire
        opt%v_fire(:) = 0.d0
        opt%alpha_fire = alpha_start
      end if

      ! Euler method
      do iatom=1,natom
        if(flag_opt_atom(iatom)=='y') then
          mass_au = umass * system%Mass(Kion(iatom))
          do ixyz=1,3
            opt%v_fire(ixyz+3*(iatom-1)) = opt%v_fire(ixyz+3*(iatom-1)) &
                                         + force_1d(ixyz+3*(iatom-1))*opt%step_fire/mass_au
              system%Rion(ixyz,iatom) = system%Rion(ixyz,iatom) + opt%v_fire(ixyz+3*(iatom-1))*opt%step_fire
          end do
        end if
      end do
    end select
  
  end subroutine structure_opt

  !===============================================================finilize
  subroutine structure_opt_fin(opt)
    use structures, only: s_opt
    use parallelization, only: nproc_id_global
    use communication, only: comm_is_root
    implicit none
    type(s_opt) :: opt
    deallocate(opt%a_dRion,opt%dFion)
    deallocate(opt%Hess_mat,opt%Hess_mat_last)
    if(comm_is_root(nproc_id_global)) write(*,*) "Optimization Converged"
  end subroutine structure_opt_fin

end module structure_opt_sub
