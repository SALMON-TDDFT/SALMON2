!
!  Copyright 2020 SALMON developers
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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130


module fdtd_weyl
    use structures, only: s_fdtd_system
    implicit none

    type ls_fdtd_weyl
        type(s_scalar) :: phi, rho_em
        type(s_vector) :: vec_e, vec_h, vec_a, vec_j_em
        type(s_vector) :: vec_Ac, vec_Ac_old
    end type ls_fdtd_weyl


contains


    subroutine weyl_init(fs, fw)
        implicit none
        type(s_fdtd_system), intent(inout) :: fs
        type(ls_fdtd_weyl), intent(inout) :: fw

        call allocate_vector(fs%mg, fw%vec_e)
        call allocate_vector(fs%mg, fw%vec_h)
        call allocate_vector(fs%mg, fw%vec_j_em)
        call allocate_vector_with_ovlp(fs%mg, fw%vec_Ac)
        call allocate_vector_with_ovlp(fs%mg, fw%vec_Ac_old)
        call allocate_vector_with_ovlp(fs%mg, fw%vec_Ac_tmp)
        call allocate_scalar(fs%mg, fw%edensity_emfield)
        call allocate_scalar(fs%mg, fw%edensity_absorb)


        contains

        subroutine weyl_allocate
        end subroutine weyl_allocate
                
        subroutine weyl_check_inc
        end subroutine weyl_check_inc
        
        
    end subroutine



    subroutine weyl_calc(fs, fw)
        implicit none
        type(s_fdtd_system), intent(inout) :: fs
        type(s_fdtd_field), intent(inout) :: fw
    
    end subroutine weyl_calc


    subroutine weyl_finalize(fs, fw)
        implicit none
        type(s_fdtd_system), intent(inout) :: fs
        type(s_fdtd_field), intent(inout) :: fw

    end subroutine weyl_finalize

    !   subroutine weyl_mpi_grid_sr(fs, fw)
    !   subroutine weyl_input_shape(ifn, ng_is, ng_ie, lg_is, lg_ie, Nd, imat, format)
    !   subroutine weyl_sendrecv(fs, fw, var)
    !   subroutine weyl_fd(ista, iend, ng_is, ng_ie, Nd, c1, c2, f1, f2, f3, var, dir)
    !   subroutine weyl_save_plane(id, ipl, conv, ng_is, ng_ie, lg_is, lg_ie, Nd, ifn, iobs, iter, f, var)
    !   subroutine weyl_fourier(nt, ne, dt, de, ti, ft, fr, fi)


end module fdtd_weyl
