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

!=======================================================================


module input_checker_ms
    use salmon_global
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global, nproc_size_global
    implicit none
contains

function check_input_variables_ms() result(flag)
    implicit none
    logical :: flag

    flag = .true.
    
    if (nx_m < 1) &
        call raise("ERROR! 'nx_m' must be larger than 1!")
        
    if (ny_m /= 1) &
        call raise("ERROR! 'ny_m' must be 1!")
        
    if (nz_m /= 1) &
        call raise("ERROR! 'nz_m' must be 1!")
        
    if (nproc_size_global < nx_m * ny_m * nz_m) &
        call raise("ERROR! MPI procs is too small!")
        
    if (mod(nproc_size_global, nx_m * ny_m * nz_m) > 0) &
        call raise("ERROR! MPI procs must be an integer multiple of macropoints!")
        
    if (hx_m < 1d-6 .and. dl_em(1) < 1d-6) &
        call raise("ERROR! 'hx_m' or 'dl_em(1)' must be specified!")
        
    if (hy_m < 1d-6 .and. dl_em(2) < 1d-6) &
        call raise("ERROR! 'hy_m' or 'dl_em(2)' must be specified!")
        
    if (hz_m < 1d-6 .and. dl_em(3) < 1d-6) &
        call raise("ERROR! 'hz_m' or 'dl_em(3)' must be specified!")
        
    if (dt < 1d-6) &
        call raise("ERROR! 'dt' must be specified!")
        
    if (dt_em > 1d-6) &
        call raise("ERROR! 'dt_em' must not be specified!")
        
    if (abs(nxvacl_m) < 2) &
        call raise("ERROR! 'nxvacl_m' must not larger than 2!")
        
    if (abs(nxvacr_m) < 2) &
        call raise("ERROR! 'nxvacr_m' must not larger than 2!")
        
    if (trim(boundary_em(1,1)) .ne. 'periodic' &
        & .and. trim(boundary_em(1,1)) .ne. 'pec' &
        & .and. trim(boundary_em(1,1)) .ne. 'abc') &
            call raise("ERROR! 'boundary_em(1,1)' unknown boundary condition!")

            
    if (trim(boundary_em(1,2)) .ne. 'periodic' &
        & .and. trim(boundary_em(1,2)) .ne. 'pec' &
        & .and. trim(boundary_em(1,2)) .ne. 'abc') &
            call raise("ERROR! 'boundary_em(1,2)' unknown boundary condition!")

    if (nmacro_chunk < 1) &
        call raise("ERROR! 'nmacro_chunk' must not larger than 1!")

    return

    contains

    subroutine raise(msg)
        implicit none
        character(*), intent(in) :: msg
        if (comm_is_root(nproc_id_global)) then
            write(*, '(a)') msg
        end if
        flag = .false.
    end subroutine raise

end function check_input_variables_ms
end module input_checker_ms


