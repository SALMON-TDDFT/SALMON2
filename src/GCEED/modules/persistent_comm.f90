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
module persistent_comm
  use pack_unpack, only: array_shape

  integer,public,allocatable :: nreqs_rorbital(:),nreqs_corbital(:)

  type(array_shape),public,allocatable :: nshape_orbital(:)
  type(array_shape),public,allocatable :: nshape_groupob(:),nrange_groupob(:,:)

  public :: init_persistent_requests

private
contains
  subroutine init_persistent_requests(info)
    use structures, only: s_orbital_parallel
    implicit none
    type(s_orbital_parallel),intent(in) :: info

    call init_comm_korbital(info)
  end subroutine

  subroutine init_comm_korbital(info)
    use structures, only: s_orbital_parallel
    use init_sendrecv_sub,    only: iup_array,idw_array,jup_array,jdw_array,kup_array,kdw_array
    use salmon_communication, only: comm_send_init, comm_recv_init, comm_proc_null
    use pack_unpack,          only: create_array_shape
    use scf_data
    implicit none
    type(s_orbital_parallel),intent(in) :: info
    integer :: icomm
    integer :: iup,idw,jup,jdw,kup,kdw

    icomm=info%icomm_r

#ifdef SALMON_USE_MPI
    iup=iup_array(1)
    idw=idw_array(1)
    jup=jup_array(1)
    jdw=jdw_array(1)
    kup=kup_array(1)
    kdw=kdw_array(1)
#else
    iup=comm_proc_null
    idw=comm_proc_null
    jup=comm_proc_null
    jdw=comm_proc_null
    kup=comm_proc_null
    kdw=comm_proc_null
#endif

    allocate(nreqs_rorbital(12))
    nreqs_rorbital( 1) = comm_send_init(srmatbox1_x_3d,iup,3,icomm)
    nreqs_rorbital( 2) = comm_recv_init(srmatbox2_x_3d,idw,3,icomm)
    nreqs_rorbital( 3) = comm_send_init(srmatbox3_x_3d,idw,4,icomm)
    nreqs_rorbital( 4) = comm_recv_init(srmatbox4_x_3d,iup,4,icomm)
    nreqs_rorbital( 5) = comm_send_init(srmatbox1_y_3d,jup,5,icomm)
    nreqs_rorbital( 6) = comm_recv_init(srmatbox2_y_3d,jdw,5,icomm)
    nreqs_rorbital( 7) = comm_send_init(srmatbox3_y_3d,jdw,6,icomm)
    nreqs_rorbital( 8) = comm_recv_init(srmatbox4_y_3d,jup,6,icomm)
    nreqs_rorbital( 9) = comm_send_init(srmatbox1_z_3d,kup,7,icomm)
    nreqs_rorbital(10) = comm_recv_init(srmatbox2_z_3d,kdw,7,icomm)
    nreqs_rorbital(11) = comm_send_init(srmatbox3_z_3d,kdw,8,icomm)
    nreqs_rorbital(12) = comm_recv_init(srmatbox4_z_3d,kup,8,icomm)

    allocate(nreqs_corbital(12))
    nreqs_corbital( 1) = comm_send_init(scmatbox1_x_3d,iup,3,icomm)
    nreqs_corbital( 2) = comm_recv_init(scmatbox2_x_3d,idw,3,icomm)
    nreqs_corbital( 3) = comm_send_init(scmatbox3_x_3d,idw,4,icomm)
    nreqs_corbital( 4) = comm_recv_init(scmatbox4_x_3d,iup,4,icomm)
    nreqs_corbital( 5) = comm_send_init(scmatbox1_y_3d,jup,5,icomm)
    nreqs_corbital( 6) = comm_recv_init(scmatbox2_y_3d,jdw,5,icomm)
    nreqs_corbital( 7) = comm_send_init(scmatbox3_y_3d,jdw,6,icomm)
    nreqs_corbital( 8) = comm_recv_init(scmatbox4_y_3d,jup,6,icomm)
    nreqs_corbital( 9) = comm_send_init(scmatbox1_z_3d,kup,7,icomm)
    nreqs_corbital(10) = comm_recv_init(scmatbox2_z_3d,kdw,7,icomm)
    nreqs_corbital(11) = comm_send_init(scmatbox3_z_3d,kdw,8,icomm)
    nreqs_corbital(12) = comm_recv_init(scmatbox4_z_3d,kup,8,icomm)

    allocate(nshape_orbital(3))
    nshape_orbital(1) = create_array_shape(mg_sta(1)-Nd, mg_end(1)+Nd)
    nshape_orbital(2) = create_array_shape(mg_sta(2)-Nd, mg_end(2)+Nd)
    nshape_orbital(3) = create_array_shape(mg_sta(3)-Nd, mg_end(3)+Nd)
  end subroutine
end module
