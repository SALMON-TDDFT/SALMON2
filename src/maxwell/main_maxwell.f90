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
!-----------------------------------------------------------------------------------------
subroutine main_maxwell
  use parallelization, only: nproc_id_global
  use communication,   only: comm_is_root
  use structures,      only: s_fdtd_system
  use fdtd_eh,         only: ls_fdtd_eh,eh_init,eh_calc,eh_finalize
  use misc_routines,   only: get_wtime
  implicit none
  type(s_fdtd_system) :: fs
  type(ls_fdtd_eh)    :: fe
  real(8)             :: elapsed_time
  
  elapsed_time=get_wtime()
  call eh_init(fs,fe)
  call eh_calc(fs,fe)
  call eh_finalize(fs,fe)
  elapsed_time=get_wtime()-elapsed_time
  if(comm_is_root(nproc_id_global)) then
    write(*,'(A,f16.8)') " elapsed time [s] = ", elapsed_time
    write(*,*)
  end if
  
  return
end subroutine main_maxwell
