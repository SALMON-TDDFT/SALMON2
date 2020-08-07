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
module jellium
  implicit none

contains
  
  !===========================================================================================
  != check condition =========================================================================
  subroutine check_condition_jm
    use salmon_global,   only: yn_md, yn_opt, yn_out_pdos, yn_out_tm, yn_out_rvf_rt, nelem, natom, nelec, spin, xc, &
                               yn_periodic, layout_multipole, shape_file_jm, num_jm, nelec_jm, method_singlescale
    use parallelization, only: nproc_id_global
    use communication,   only: comm_is_root
    implicit none
    
    call condition_yn_jm(yn_md,        'yn_md',        'n')
    call condition_yn_jm(yn_opt,       'yn_opt',       'n')
    call condition_yn_jm(yn_out_pdos,  'yn_out_pdos',  'n')
    call condition_yn_jm(yn_out_tm,    'yn_out_tm',    'n')
    call condition_yn_jm(yn_out_rvf_rt,'yn_out_rvf_rt','n')

    call condition_int_jm(nelem,'nelem',1)
    call condition_int_jm(natom,'natom',1)
    
    if(yn_periodic=='n'.and.layout_multipole/=1) then
      if(comm_is_root(nproc_id_global)) &
        write(*,'("For yn_jm = y and yn_periodic = n, layout_multipole must be 1.")')
      stop
    end if
    
    if(mod(nelec,2)/=0) then
      if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, nelec must be even number.")')
      stop
    end if
    
    if(trim(spin)/='unpolarized') then
      if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, spin must be even unpolarized.")')
      stop
    end if
    
    if(trim(xc)/='pz') then
      if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, xc must be pz.")')
      stop
    end if
    
    if(trim(method_singlescale)/='3d') then
      if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, method_singlescale must be 3d.")')
      stop
    end if
    
    if (trim(shape_file_jm)=='none' .and. nelec/=sum(nelec_jm(:)))then
      if(comm_is_root(nproc_id_global)) &
        write(*,'("For yn_jm = y and shape_file_jm = none, nelec must be sum(nelec_jm).")')
      stop
    end if
    
    if(num_jm<1) then
      if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, num_jm must be larger than 0.")')
      stop
    end if
    
    return
  contains
    
    !+ CONTAINED IN check_condition_jm +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check condition for y/n +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine condition_yn_jm(tar,name,ans)
      implicit none
      character(1),intent(in) :: tar
      character(*),intent(in) :: name
      character(1),intent(in) :: ans
      
      if (tar/=ans) then
        if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, ",A," must be ",A,".")') name,ans
        stop
      end if
      
      return
    end subroutine condition_yn_jm
    
    !+ CONTAINED IN check_condition_jm +++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check condition for integer +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine condition_int_jm(tar,name,ans)
      implicit none
      integer,     intent(in) :: tar
      character(*),intent(in) :: name
      integer,     intent(in) :: ans
      
      if (tar/=ans) then
        if(comm_is_root(nproc_id_global)) write(*,'("For yn_jm = y, ",A," must be ",I4,".")') name,ans
        stop
      end if
      
      return
    end subroutine condition_int_jm
    
  end subroutine check_condition_jm
  
  !===========================================================================================
  != meke positive back ground charge density ================================================
  subroutine make_rho_jm(lg,mg,info,system,rho_jm)
    use salmon_global,   only: shape_file_jm, num_jm, nelec_jm, rs_bohr_jm, sphere_loc_jm, &
                               yn_charge_neutral_jm, yn_output_dns_jm, yn_periodic, nelec, unit_system
    use inputoutput,     only: ulength_from_au
    use structures,      only: s_rgrid, s_dft_system, s_parallel_info, s_scalar, allocate_scalar
    use parallelization, only: nproc_id_global, nproc_group_global
    use communication,   only: comm_is_root, comm_summation
    use common_maxwell,  only: input_shape_em
    use write_file3d,    only: write_cube
    use math_constants,  only: pi
    implicit none
    type(s_rgrid),         intent(in)    :: lg, mg
    type(s_parallel_info), intent(in)    :: info
    type(s_dft_system),    intent(in)    :: system
    type(s_scalar),        intent(inout) :: rho_jm
    type(s_scalar)                       :: work_l1,work_l2
    integer,allocatable :: imedia(:,:,:)
    integer             :: ii, ix, iy, iz, mod_nelec
    real(8),allocatable :: dens(:), radi(:), mod_rs_bohr_jm(:)
    real(8)             :: rab, charge_sum, charge_error 
    character(60)       :: suffix
    character(30)       :: phys_quantity
    
    !set density
    allocate(dens(num_jm)); dens(:)=0.0d0;
    do ii=1,num_jm
      dens(ii) = 1.0d0/(4.0d0*pi/3.0*(rs_bohr_jm(ii)**3.0d0))
    end do
    
    !make rho_jm
    if (trim(shape_file_jm)=='none')then
      !**************************************************************************************!
      !*** rho_jm is generated by spherecal shapes ******************************************!
      !**************************************************************************************!
      !allocate radius
      allocate(radi(num_jm)); radi(:)=0.0d0;
      
      !make spheres
      do ii=1,num_jm
        !set radius
        radi(ii) = ( dble(nelec_jm(ii))/dens(ii)/(4.0d0*pi/3.0) )**(1.0d0/3.0d0)
        
        !make ii-th sphere
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          rab = sqrt( (lg%coordinate(ix,1)-sphere_loc_jm(ii,1))**2.0d0  &
                     +(lg%coordinate(iy,2)-sphere_loc_jm(ii,2))**2.0d0  &
                     +(lg%coordinate(iz,3)-sphere_loc_jm(ii,3))**2.0d0 )
          if(rab<=radi(ii)) rho_jm%f(ix,iy,iz)=dens(ii)
        end do
        end do
        end do
      end do
      
      !check charge neutrality
      call check_neutral_jm(charge_sum,charge_error,sum(nelec_jm(:)))
    else
      !**************************************************************************************!
      !*** rho_jm is generated by cube file *************************************************!
      !**************************************************************************************!
      !input shape
      allocate(imedia(mg%is(1):mg%ie(1),mg%is(2):mg%ie(2),mg%is(3):mg%ie(3))); imedia(:,:,:)=0;
      if(comm_is_root(nproc_id_global)) write(*,*)
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
      if(index(shape_file_jm,".cube", back=.true.)/=0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "shape file is inputed by .cube format."
        end if
        call input_shape_em(shape_file_jm,600,mg%is,mg%ie,lg%is,lg%ie,0,imedia,'cu')
      elseif(index(shape_file_jm,".mp", back=.true.)/=0) then
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "shape file is inputed by .mp format."
          write(*,*) "This version works for only .cube format.."
        end if
        stop
      else
        if(comm_is_root(nproc_id_global)) then
          write(*,*) "shape file must be .cube or .mp formats."
        end if
        stop
      end if
      if(comm_is_root(nproc_id_global)) write(*,*) "**************************"
      
      !make rho_jm from shape file
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        if(imedia(ix,iy,iz)>0) rho_jm%f(ix,iy,iz)=dens(imedia(ix,iy,iz))
      end do
      end do
      end do
      
      !check charge neutrality
      call check_neutral_jm(charge_sum,charge_error,nelec)
    end if
    
    
    !propose modified parameter & stop
    !or modify parameters
    allocate(mod_rs_bohr_jm(num_jm)); mod_rs_bohr_jm(:)=0.0d0;
    if(charge_error>=2.0d0/dble(nelec))then
      !stop & propose modified parameter
      mod_nelec = int(charge_sum)
      if(mod(mod_nelec,2)/=0) mod_nelec=mod_nelec+1
      if(comm_is_root(nproc_id_global))then
        write(*,*)
        write(*,'("Charge nertrality error is",E23.15E3,".")') charge_error
        write(*,'("To improve charge nertrality, change nelec to",I9,".")') mod_nelec
      end if
      stop
    else
      !modify parameters and recheck neutrality
      if(yn_charge_neutral_jm=='y')then
        rho_jm%f(:,:,:)   = rho_jm%f(:,:,:) * ( dble(nelec)/charge_sum )
        dens(:)           = dens(:)         * ( dble(nelec)/charge_sum )
        mod_rs_bohr_jm(:) = ( 1.0d0/(4.0d0*pi/3.0*dens(:)) )**(1.0d0/3.0d0)
        if (trim(shape_file_jm)=='none')then
          call check_neutral_jm(charge_sum,charge_error,sum(nelec_jm(:)))
        else
          call check_neutral_jm(charge_sum,charge_error,nelec)
        end if
      end if
    end if
    
    !output cube file
    if(yn_output_dns_jm=='y') then
      call allocate_scalar(lg,work_l1); call allocate_scalar(lg,work_l2);
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        work_l1%f(ix,iy,iz) = rho_jm%f(ix,iy,iz)
      end do
      end do
      end do
      call comm_summation(work_l1%f,work_l2%f,lg%num(1)*lg%num(2)*lg%num(3),nproc_group_global)
      suffix = "dns_jellium"; phys_quantity = "pbcd";
      call write_cube(lg,103,suffix,phys_quantity,work_l2%f,system)
    end if
    
    !write information
    if(comm_is_root(nproc_id_global))then
      write(*,*)
      write(*,*) '****************** Jellium information ******************'
      if(trim(shape_file_jm)=='none')then
        write(*,'("  Positive background charge density is generated by spherecal shapes:")')
        do ii=1,num_jm
          write(*, '(A,I3,A,E23.15E3)') '    Radius of sphere(',ii,') =', radi(ii)*ulength_from_au
        end do
        write(*,*) " in the unit system, ",trim(unit_system),"."
      else
        write(*,'("  Positive background charge density is generated by shape file.")')
      end if
      if(sum(mod_rs_bohr_jm(:))==0.0d0) then
        write(*,'("  Wigner-Seitz radius is set as follows:")')
        do ii=1,num_jm
          write(*, '(A,I3,A,E23.15E3)') '    rs_bohr_jm(',ii,') =', rs_bohr_jm(ii)
        end do
      else
        write(*,'("  To keep charge neutrality, Wigner-Seitz radius is modified as follows:")')
        do ii=1,num_jm
          write(*, '(A,I3,A,E23.15E3)') '    mod_rs_bohr_jm(',ii,') =', mod_rs_bohr_jm(ii)
        end do
      end if
      write(*,*) " in the atomic unit(Bohr)."
      write(*,'(A,E23.15E3," %")') '  Chrge neutrality error =', charge_error
      if(yn_periodic=='y') then
        write(*,*)
        write(*,'("  For yn_jm = y and yn_periodic=y, this version still cannot output Total Energy.")')
        write(*,*)
      end if
      write(*,*) '*********************************************************'
      write(*,*)
    end if
    
    !change sign in view of electron density
    rho_jm%f(:,:,:) = -rho_jm%f(:,:,:)
    
    return
   contains
    
    !+ CONTAINED IN make_rho_jm ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !+ check charge neutrality +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine check_neutral_jm(sum_c,err,num_e)
      implicit none
      real(8), intent(inout) :: sum_c
      real(8), intent(out)   :: err
      integer, intent(in)    :: num_e
      real(8)                :: sum_tmp
      
      sum_tmp = 0.0d0; sum_c = 0.0d0;
!$omp parallel
!$omp do private(ix,iy,iz) reduction( + : sum_tmp )
      do iz=mg%is(3),mg%ie(3)
      do iy=mg%is(2),mg%ie(2)
      do ix=mg%is(1),mg%ie(1)
        sum_tmp = sum_tmp + rho_jm%f(ix,iy,iz)
      end do
      end do
      end do
!$omp end do
!$omp end parallel
      call comm_summation(sum_tmp,sum_c,info%icomm_r)
      sum_c = sum_c * system%hvol
      err   = abs( (sum_c - dble(num_e)) / dble(num_e) )
      
      return
    end subroutine check_neutral_jm
    
  end subroutine make_rho_jm
  
end module jellium
