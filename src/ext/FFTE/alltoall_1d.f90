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
!======================================================================
subroutine alltoall_1d(vec_a,vec_b,n,icomm,npu)
  use communication, only: comm_alltoall
  implicit none
  integer,intent(in)      :: n
  complex(8),intent(in)   :: vec_a(n)
  complex(8),intent(out)  :: vec_b(n)
  integer,intent(in)      :: icomm
  integer,intent(in)      :: npu

  call comm_alltoall(vec_a,vec_b,icomm,n/npu)

end subroutine alltoall_1d
!======================================================================
