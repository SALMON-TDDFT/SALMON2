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
MODULE scf_data
use salmon_global
implicit none

!-------------------- Global variables

integer :: itotNtime

real(8),allocatable :: vloc_t(:,:,:,:)
real(8),allocatable :: vloc_new(:,:,:,:)
real(8),allocatable :: vloc_old(:,:,:,:,:)

integer :: itt

character(100):: file_Projection


!real(8),allocatable :: curr(:,:)


integer :: CONTEXT, IAM, MYCOL, MYROW, NPCOL, NPROCS2, NPROW
integer :: DESCA( 50 ), DESCZ( 50 )

real(8),allocatable :: vonf_sd(:,:,:),eonf_sd(:,:,:,:)

!filename
character(100) :: file_eigen
character(100) :: file_RT_q
character(100) :: file_alpha_q
character(100) :: file_RT_e
character(100) :: file_RT_dip2
character(100) :: file_alpha_dip2
character(100) :: file_RT_dip2_q
character(100) :: file_alpha_dip2_q
character(100) :: file_RT_dip2_e
character(100) :: file_external
character(100) :: file_out_rt_bin
character(100) :: file_in_rt_bin

END MODULE scf_data

