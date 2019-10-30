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
!-----------------------------------------------------------------------------------------
subroutine weyl_init(fs,ff,fw)
  use structures,     only: s_fdtd_system, s_fdtd_field, allocate_vector
  use salmon_maxwell, only: ls_fdtd_work
  implicit none
  type(s_fdtd_system) :: fs
  type(s_fdtd_field)  :: ff
  type(ls_fdtd_work)  :: fw
  
  ! Allocation of the Vector field components:
  call allocate_vector(fs%mg, ff%vec_e)
  call allocate_vector(fs%mg, ff%vec_h)
  call allocate_vector(fs%mg, ff%vec_j_em)
  call allocate_vector_array(fs%mg, ff%vec_Ac)
  call allocate_vector_array(fs%mg, ff%vec_Ac_old)
  call allocate_vector_array(fs%mg, fw%vec_Ac_tmp)
  call allocate_scalar(fs%mg, fw%edensity_emfield)
  call allocate_scalar(fs%mg, fw%edensity_absorb)
  return
end subroutine weyl_init
