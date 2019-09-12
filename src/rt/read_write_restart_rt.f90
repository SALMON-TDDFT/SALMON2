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
module read_write_restart_rt_sub
  implicit none

contains

  subroutine init_dir_out_restart(ofile)
    use structures,  only: s_ofile
    use inputoutput, only: theory
    use filesystem,  only: create_directory
    use salmon_communication, only: comm_is_root,comm_sync_all
    use salmon_parallel,      only: nproc_id_global
    implicit none
    type(s_ofile), intent(inout) :: ofile

    select case(theory)
    case('DFT','DFT_MD', &
         'TDDFT_response','TDDFT_pulse','Single_scale_Maxwell_TDDFT')
       ofile%dir_out_restart = 'data_for_restart/'
    end select

    !!! currently, set "./" : change later
    ofile%dir_out_restart = './'

    if(ofile%dir_out_restart(1:3).ne."./ ") then
      if(comm_is_root(nproc_id_global)) then
        if(.not. create_directory(ofile%dir_out_restart)) then
          stop 'fail: init_dir_out_restart::create_directory'
        end if
      end if
    endif

    call comm_sync_all ! sync until directory created

  end subroutine init_dir_out_restart
  
  subroutine write_checkpoint_rt(itt,ng,info)
    use structures, only: s_rgrid,s_orbital_parallel
    implicit none
    type(s_rgrid), intent(in) :: ng
    type(s_orbital_parallel), intent(in) :: info

    integer :: itt
    character(256) :: dir_checkpoint

    call create_checkpoint_dir(itt,dir_checkpoint)
    call write_rt_bin(dir_checkpoint,ng,info)

  end subroutine write_checkpoint_rt

  subroutine create_checkpoint_dir(itt,dir_checkpoint)
    use filesystem, only: create_directory
    use salmon_communication, only: comm_is_root,comm_sync_all
    use salmon_parallel,      only: nproc_id_global
    implicit none
    integer :: itt
    character(256) :: dir_checkpoint

    write(dir_checkpoint,'(A,I6.6,A)') "checkpoint_rt_",itt,"/"
    if(comm_is_root(nproc_id_global)) then
      if(.not. create_directory(dir_checkpoint)) then
        stop 'fail: create_checkpoint_dir::create_directory'
      end if
    end if

    call comm_sync_all ! sync until directory created

  end subroutine create_checkpoint_dir

end module read_write_restart_rt_sub
