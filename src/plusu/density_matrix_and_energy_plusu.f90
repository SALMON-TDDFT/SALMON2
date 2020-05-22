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
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130

module density_matrix_and_energy_plusU_sub

  use plusU_global, only: dm_mms_nla, prep_Hubbard_parameters &
                         ,PLUS_U_ON, U_eff, V_eff

  implicit none

  private
  public :: calc_density_matrix_and_energy_plusU
  public :: PLUS_U_ON

contains

  subroutine calc_density_matrix_and_energy_plusU( psi, ppg, info, system, E_U )
    use structures, only: s_orbital, s_pp_grid, s_parallel_info, s_dft_system
    implicit none
    type(s_orbital),intent(in) :: psi
    type(s_pp_grid),intent(in) :: ppg
    type(s_parallel_info),intent(in) :: info
    type(s_dft_system),intent(in) :: system
    real(8),intent(out) :: E_U
    integer :: Nlma,Nspin
    integer :: im,ik,io,ispin,ilma,jlma,iprj,Nproj_pairs
    integer :: ix,iy,iz,j,ia,io_s,io_e,ik_s,ik_e,im_s,im_e
    integer :: a,l,n,m1,m2
    complex(8) :: phipsi,ztmp,ztmp2
    complex(8),allocatable :: phipsi_lma(:)
    complex(8),parameter :: zero=(0.0d0,0.0d0)

    Nspin=system%nspin
    Nproj_pairs=size(ppg%proj_pairs_ao,2)
    Nlma=size(ppg%ia_tbl_ao)

    im_s=info%im_s
    im_e=info%im_e
    ik_s=info%ik_s
    ik_e=info%ik_e
    io_s=info%io_s
    io_e=info%io_e

    if( .not.allocated(V_eff) )then
      allocate( V_eff(Nproj_pairs,Nspin) )
      V_eff=0.0d0
    end if

    if( .not.allocated(U_eff) )then
      a=maxval( ppg%proj_pairs_info_ao(1,:) )
      l=maxval( ppg%proj_pairs_info_ao(2,:) )
      n=maxval( ppg%proj_pairs_info_ao(3,:) )
      allocate( U_eff(n,0:l,a) ); U_eff=0.0d0
      call prep_Hubbard_parameters( U_eff )
    end if

    if( .not.allocated(dm_mms_nla) )then
      a=maxval( ppg%proj_pairs_info_ao(1,:) )
      l=maxval( ppg%proj_pairs_info_ao(2,:) )
      n=maxval( ppg%proj_pairs_info_ao(3,:) )
      allocate( dm_mms_nla(-l:l,-l:l,Nspin,n,0:l,a) )
    end if
    dm_mms_nla=zero

    allocate( phipsi_lma(Nlma) ); phipsi_lma=zero

    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      if ( abs(system%rocc(io,ik,ispin) * system%wtk(ik)) < 1.d-10 ) cycle

      do ilma=1,Nlma
        ia = ppg%ia_tbl_ao(ilma)
        phipsi = zero
        do j=1,ppg%mps_ao(ia)
          ix = ppg%jxyz_ao(1,j,ia)
          iy = ppg%jxyz_ao(2,j,ia)
          iz = ppg%jxyz_ao(3,j,ia)
          phipsi = phipsi + conjg( ppg%zekr_phi_ao(j,ilma,ik) ) &
                          * psi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        phipsi_lma(ilma) = phipsi * system%Hvol
      end do !ilma

      do iprj=1,Nproj_pairs
        ilma = ppg%proj_pairs_ao(1,iprj)
        jlma = ppg%proj_pairs_ao(2,iprj)
        a = ppg%proj_pairs_info_ao(1,iprj)
        l = ppg%proj_pairs_info_ao(2,iprj)
        n = ppg%proj_pairs_info_ao(3,iprj)
        m1= ppg%proj_pairs_info_ao(4,iprj)
        m2= ppg%proj_pairs_info_ao(5,iprj)
        dm_mms_nla(m1,m2,ispin,n,l,a) = dm_mms_nla(m1,m2,ispin,n,l,a) &
           + system%rocc(io,ik,ispin) * system%wtk(ik) * phipsi_lma(ilma)*phipsi_lma(jlma)
      end do

    end do !ispin
    end do !io
    end do !ik
    end do !im

    E_U=0.0d0
    do a=1,size(dm_mms_nla,6)
    do l=0,size(dm_mms_nla,5)-1
    do n=1,size(dm_mms_nla,4)
      ztmp2=zero
      do ispin=1,Nspin
      do m2=-l,l
        ztmp = dm_mms_nla(m2,m2,ispin,n,l,a)
        do m1=-l,l
          ztmp = ztmp - dm_mms_nla(m2,m1,ispin,n,l,a) &
                      * dm_mms_nla(m1,m2,ispin,n,l,a)
        end do !m1
        ztmp2 = ztmp2 + ztmp
      end do !m2
      end do !ispin
      ztmp2 = 0.5d0*U_eff(n,l,a)*ztmp2
      !write(*,'(1x,3i4,2g25.15)') a,l,n,real(ztmp2),aimag(ztmp2)
      E_U = E_U + ztmp2
    end do !n
    end do !l
    end do !a

    V_eff=0.0d0
    do ispin=1,Nspin
      do iprj=1,Nproj_pairs
        ilma = ppg%proj_pairs_ao(1,iprj)
        jlma = ppg%proj_pairs_ao(2,iprj)
        a = ppg%proj_pairs_info_ao(1,iprj)
        l = ppg%proj_pairs_info_ao(2,iprj)
        n = ppg%proj_pairs_info_ao(3,iprj)
        m1= ppg%proj_pairs_info_ao(4,iprj)
        m2= ppg%proj_pairs_info_ao(5,iprj)
        if( m1 == m2 )then
           V_eff(iprj,ispin) = U_eff(n,l,a)*( 0.5d0 - dm_mms_nla(m1,m2,ispin,n,l,a) )
        else
           V_eff(iprj,ispin) = U_eff(n,l,a)*( -dm_mms_nla(m1,m2,ispin,n,l,a) )
        end if
      end do !iprj
    end do !ispin

    return
  end subroutine calc_density_matrix_and_energy_plusU

end module density_matrix_and_energy_plusU_sub
