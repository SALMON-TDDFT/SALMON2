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


module input_checker_sbe
    use salmon_global
    use communication, only: comm_is_root
    use parallelization, only: nproc_id_global, nproc_size_global
    implicit none
contains

function check_input_variables_sbe() result(flag)
    implicit none
    logical :: flag

    flag = .true.

    if ((norm2(al_vec1) <= 1d-9) .and. (norm2(al_vec2) <= 1d-9) .and. (norm2(al_vec3) <= 1d-9)) then
        if (norm2(al) > 1d-9) then
            al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
            al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
            al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
        else
            call raise("ERROR! 'al' or 'al_vec1..3' must be specified!")
        end if
    end if

    if (nstate_sbe <= 0) nstate_sbe = nstate
    if (nstate_sbe > nstate) then
        call raise("ERROR! 'nstate_sbe' must be smaller than 'nstate'!")
    end if

    if (trim(theory) .ne. "maxwell_sbe") return

    if (nx_m < 1) &
        call raise("ERROR! 'nx_m' must be larger than 1!")
    if (ny_m < 1) &
        call raise("ERROR! 'ny_m' must be larger than 1!")
    if (nz_m < 1) &
        call raise("ERROR! 'nz_m' must be larger than 1!")

    if ((nxvac_m(1) > 0) .and. (nxvacl_m > 0)) &
        call raise("ERROR! 'nxvac_m' and 'nxvacl_m' must not be specified simultaneously!")
    if ((nxvac_m(2) > 0) .and. (nxvacr_m > 0)) &
        call raise("ERROR! 'nxvac_m' and 'nxvacr_m' must not be specified simultaneously!")

    if ((hx_m > 0) .and. (dl_em(1) > 0)) &
        call raise("ERROR! 'hx_m' and 'dl_em' must not be specified simultaneously!")
    if ((hy_m > 0) .and. (dl_em(2) > 0)) &
        call raise("ERROR! 'hy_m' and 'dl_em' must not be specified simultaneously!")
    if ((hz_m > 0) .and. (dl_em(3) > 0)) &
        call raise("ERROR! 'hz_m' and 'dl_em' must not be specified simultaneously!")
    if ((hx_m <= 0) .and. (dl_em(1) <= 0)) &
        call raise("ERROR! 'hx_m' or 'dl_em(1)' must be specified!")
    if ((hy_m <= 0) .and. (dl_em(2) <= 0)) &
        call raise("ERROR! 'hy_m' or 'dl_em(2)' must be specified!")
    if ((hz_m <= 0) .and. (dl_em(3) <= 0)) &
        call raise("ERROR! 'hz_m' or 'dl_em(1)' must be specified!")
        
    select case(trim(fdtddim))
    case("1d", "1D")
    case("3d", "3D")
    case default
        call raise("ERROR! 'fdtddim' must be specified ('1d' or '3d')")
    end select

    if (abs(nxvacl_m) >= 1) nxvac_m(1) = abs(nxvacl_m)
    if (abs(nxvacr_m) >= 1) nxvac_m(2) = abs(nxvacr_m)

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

end function check_input_variables_sbe
end module input_checker_sbe


