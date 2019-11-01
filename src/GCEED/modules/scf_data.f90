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
!-------------------- Parameters
integer, parameter :: maxntmg=10

! Physical constants
real(8),parameter :: E2=14.4d0, H2M=3.81d0, a_B=0.529177d0
real(8),parameter :: Ry=13.6058d0
real(8),parameter :: fs2eVinv = 1.51925d0

! Coefficients of finite-difference
integer,parameter :: Nd=4       ! (2*Nd+1)-points finite-difference used in Taylor expansion
integer,parameter :: Ndh=4       ! (2*Nd+1)-points finite-difference used in Hartree routine
real(8),parameter :: cN0=-205.d0/72.d0                 ! for Laplacian
real(8),parameter :: cN1=8.d0/5.d0  , cN2=-1.d0/5.d0
real(8),parameter :: cN3=8.d0/315.d0, cN4=-1.d0/560.d0
real(8),parameter :: bN1=4.d0/5.d0  , bN2=-1.d0/5.d0   ! for Nabla
real(8),parameter :: bN3=4.d0/105.d0, bN4=-1.d0/280.d0

real(8),parameter :: cN0_Nd1=-2.d0
real(8),parameter :: cN1_Nd1=1.d0
real(8),parameter :: bN1_Nd1=1.d0/2.d0

real(8),parameter :: cN0_Nd2=-5.d0/2.d0
real(8),parameter :: cN1_Nd2=4.d0/3.d0, cN2_Nd2=-1.d0/12.d0
real(8),parameter :: bN1_Nd2=2.d0/3.d0, bN2_Nd2=-1.d0/12.d0

real(8),parameter :: cN0_Nd3=-49.d0/18.d0 
real(8),parameter :: cN1_Nd3=3.d0/2.d0, cN2_Nd3=-3.d0/20.d0
real(8),parameter :: cN3_Nd3=1.d0/90.d0
real(8),parameter :: bN1_Nd3=3.d0/4.d0, bN2_Nd3=-3.d0/20.d0 
real(8),parameter :: bN3_Nd3=1.d0/60.d0

real(8),parameter :: cN0_Nd5=-5269.d0/1800.d0 
real(8),parameter :: cN1_Nd5=5.d0/3.d0, cN2_Nd5=-5.d0/21.d0
real(8),parameter :: cN3_Nd5=5.d0/126.d0, cN4_Nd5=-5.d0/1008.d0
real(8),parameter :: cN5_Nd5=1.d0/3150.d0

real(8),parameter :: cN0_Nd6=-5369.d0/1800.d0 
real(8),parameter :: cN1_Nd6=12.d0/7.d0, cN2_Nd6=-15.d0/56.d0
real(8),parameter :: cN3_Nd6=10.d0/189.d0, cN4_Nd6=-1.d0/112.d0
real(8),parameter :: cN5_Nd6=2.d0/1925.d0, cN6_Nd6=-1.d0/16632.d0

real(8),parameter :: cN0_Nd7=-266681.d0/88200.d0
real(8),parameter :: cN1_Nd7=7.d0/4.d0, cN2_Nd7=-7.d0/24.d0
real(8),parameter :: cN3_Nd7=7.d0/108.d0, cN4_Nd7=-7.d0/528.d0
real(8),parameter :: cN5_Nd7=7.d0/3300.d0, cN6_Nd7=-7.d0/30888.d0
real(8),parameter :: cN7_Nd7=1.d0/84084.d0

real(8),parameter :: cN0_Nd8=-1077749.d0/352800.d0 
real(8),parameter :: cN1_Nd8=16.d0/9.d0, cN2_Nd8=-14.d0/45.d0
real(8),parameter :: cN3_Nd8=112.d0/1485.d0, cN4_Nd8=-7.d0/396.d0
real(8),parameter :: cN5_Nd8=112.d0/32175.d0, cN6_Nd8=-2.d0/3861.d0
real(8),parameter :: cN7_Nd8=16.d0/315315.d0, cN8_Nd8=-1.d0/411840.d0

real(8),parameter :: cN0_Nd9=-9778141.d0/3175200.d0
real(8),parameter :: cN1_Nd9=9.d0/5.d0, cN2_Nd9=-18.d0/55.d0
real(8),parameter :: cN3_Nd9=14.d0/165.d0, cN4_Nd9=-63.d0/2860.d0
real(8),parameter :: cN5_Nd9=18.d0/3575.d0, cN6_Nd9=-2.d0/2145.d0
real(8),parameter :: cN7_Nd9=9.d0/70070.d0, cN8_Nd9=-9.d0/777920.d0
real(8),parameter :: cN9_Nd9=1.d0/1969110.d0

real(8),parameter :: cN0_Nd10=-1968329.d0/635040.d0
real(8),parameter :: cN1_Nd10=20.d0/11.d0, cN2_Nd10=-15.d0/44.d0
real(8),parameter :: cN3_Nd10=40.d0/429.d0, cN4_Nd10=-15.d0/572.d0
real(8),parameter :: cN5_Nd10=24.d0/3575.d0, cN6_Nd10=-5.d0/3432.d0
real(8),parameter :: cN7_Nd10=30.d0/119119.d0, cN8_Nd10=-5.d0/155584.d0
real(8),parameter :: cN9_Nd10=10.d0/3741309.d0, cN10_Nd10=-1.d0/9237800.d0

real(8),parameter :: cN0_Nd11=-239437889.d0/76839840.d0
real(8),parameter :: cN1_Nd11=11.d0/6.d0, cN2_Nd11=-55.d0/156.d0
real(8),parameter :: cN3_Nd11=55.d0/546.d0, cN4_Nd11=-11.d0/364.d0
real(8),parameter :: cN5_Nd11=11.d0/1300.d0, cN6_Nd11=-11.d0/5304.d0
real(8),parameter :: cN7_Nd11=55.d0/129948.d0, cN8_Nd11=-55.d0/806208.d0
real(8),parameter :: cN9_Nd11=11.d0/1360476.d0, cN10_Nd11=-11.d0/17635800.d0
real(8),parameter :: cN11_Nd11=1.d0/42678636.d0

real(8),parameter :: cN0_Nd12=-240505109.d0/76839840.d0
real(8),parameter :: cN1_Nd12=24.d0/13.d0, cN2_Nd12=-33.d0/91.d0
real(8),parameter :: cN3_Nd12=88.d0/819.d0, cN4_Nd12=-99.d0/2912.d0
real(8),parameter :: cN5_Nd12=396.d0/38675.d0, cN6_Nd12=-11.d0/3978.d0
real(8),parameter :: cN7_Nd12=132.d0/205751.d0, cN8_Nd12=-33.d0/268736.d0
real(8),parameter :: cN9_Nd12=44.d0/2380833.d0, cN10_Nd12=-3.d0/1469650.d0
real(8),parameter :: cN11_Nd12=12.d0/81800719.d0, cN12_Nd12=-1.d0/194699232.d0

real(8) :: cnmat(0:12,12),bnmat(4,4)

!-------------------- Global variables
integer,allocatable :: ista_Mxin(:,:),iend_Mxin(:,:),inum_Mxin(:,:)
integer,allocatable :: ista_Mxin_s(:,:),iend_Mxin_s(:,:),inum_Mxin_s(:,:)

integer :: Miter       ! Total number of Iteration for SCF calculation
integer :: Miter_rt    ! Total number of Iteration for RT calculation

integer :: iflag_diisjump

integer :: ilsda

integer :: MST(2),ifMST(2),itotMST

real(8) :: rnetot

real(8) :: Hgs(3)        ! Grid spacing
real(8) :: Hvol


!Nonlinear core correction
logical :: flag_nlcc = .false.


real(8),allocatable :: rho(:,:,:)       ! Single particle density
real(8),allocatable :: rho0(:,:,:)      ! Single particle density
integer :: ihpsieff
integer :: iSCFRT


integer, allocatable :: idiis_sd(:)


complex(8), allocatable :: zc(:)
real(8), allocatable :: Dp(:,:), Qp(:,:,:), rIe(:)
real(8), allocatable :: tene(:)
real(8) :: vecDs(3)

integer :: itotNtime

integer :: lg_sta(3),lg_end(3),lg_num(3)
integer :: mg_sta(3),mg_end(3),mg_num(3)
integer :: ng_sta(3),ng_end(3),ng_num(3)

integer :: iobnum

integer :: k_sta,k_end,k_num

integer :: num_kpoints_3d(3)
integer :: num_kpoints_rd

real(8),allocatable :: wtk(:)

real(8),allocatable :: norm_diff_psi_stock(:,:)

real(8) :: Etot
real(8) :: Exc

integer :: iflag_subspace_diag

real(8),allocatable :: vloc_t(:,:,:,:)
real(8),allocatable :: vloc_new(:,:,:,:)
real(8),allocatable :: vloc_old(:,:,:,:,:)

integer :: iDiter(maxntmg)

integer :: itt
integer :: ikind_eext   !0:No external field, 1: dipoleApprox


!some of following variables can be removed by small changes
real(8) :: Fst
real(8) :: romega, romega2(2)
real(8) :: pulse_T, pulse_T2(2) 
real(8) :: rlaser_I, rlaser_I2(2) 
real(8) :: tau, tau2(2), delay, rcycle

character(2)  :: denplane  ! plane for writing density (xy, yz, xz)
integer       :: idensum   ! whether density is summed up along direction
                           ! perpendicular to the plane
                           ! (0: not summed, 1: summed)
real(8)       :: posplane  ! position of the plane
                           ! (only for idensum = 0)

character(100):: file_Projection


real(8),allocatable :: curr(:,:)


integer :: ilasbound_sta(3),ilasbound_end(3)
real(8) :: rlaser_center(3)

integer :: iblacsinit
integer :: CONTEXT, IAM, MYCOL, MYROW, NPCOL, NPROCS2, NPROW
integer :: DESCA( 50 ), DESCZ( 50 )

real(8),allocatable :: A_ext(:,:)
real(8),allocatable :: A_ind(:,:)
real(8),allocatable :: A_tot(:,:)
real(8),allocatable :: E_ext(:,:)
real(8),allocatable :: E_ind(:,:)
real(8),allocatable :: E_tot(:,:)

real(8),allocatable :: vonf_sd(:,:,:),eonf_sd(:,:,:,:)

!filename
character(100) :: file_eigen
character(100) :: file_alpha_lr
character(100) :: file_alpha_pulse
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


CONTAINS

!======================================================================
subroutine old_mesh(lg,mg,ng,system,info)
use salmon_parallel, only: nproc_size_global
use structures
implicit none
type(s_rgrid),intent(in) :: lg,mg,ng
type(s_dft_system),intent(in) :: system
type(s_orbital_parallel),intent(in) :: info

lg_sta(1:3) = lg%is(1:3)
lg_end(1:3) = lg%ie(1:3)
lg_num(1:3) = lg%num(1:3)

mg_sta(1:3) = mg%is(1:3)
mg_end(1:3) = mg%ie(1:3)
mg_num(1:3) = mg%num(1:3)

ng_sta(1:3) = ng%is(1:3)
ng_end(1:3) = ng%ie(1:3)
ng_num(1:3) = ng%num(1:3)

allocate(ista_Mxin(3,0:nproc_size_global-1),iend_Mxin(3,0:nproc_size_global-1))
allocate(inum_Mxin(3,0:nproc_size_global-1))

allocate(ista_Mxin_s(3,0:nproc_size_global-1),iend_Mxin_s(3,0:nproc_size_global-1))
allocate(inum_Mxin_s(3,0:nproc_size_global-1))

ista_mxin = mg%is_all
iend_mxin = mg%ie_all
inum_mxin = mg%ie_all - mg%is_all + 1

ista_mxin_s = ng%is_all
iend_mxin_s = ng%ie_all
inum_mxin_s = ng%ie_all - ng%is_all + 1

Hvol = system%Hvol ! future work: remove this line
Hgs = system%Hgs ! future work: remove this line

k_sta = info%ik_s ! future work: remove this line
k_end = info%ik_e ! future work: remove this line
k_num = info%numk ! future work: remove this line
iobnum = info%numo ! future work: remove this line

end subroutine old_mesh

function check_rion_update() result(rion_update)
  implicit none
  logical :: rion_update
  if (iscfrt == 1) then
    rion_update = (yn_opt == 'y')
  else if (iscfrt == 2) then
    rion_update = (yn_md == 'y')
  else
    rion_update = .true.
  end if
end function

END MODULE scf_data

