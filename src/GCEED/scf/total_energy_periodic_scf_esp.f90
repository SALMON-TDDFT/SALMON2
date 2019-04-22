!=======================================================================

SUBROUTINE total_energy_periodic_scf_esp(psi_in)
use salmon_parallel, only: nproc_group_global, nproc_group_korbital, nproc_group_h
use salmon_communication, only: comm_summation
use misc_routines, only: get_wtime
use calc_allob_sub
use scf_data
use hpsi2_sub
implicit none

complex(8) :: psi_in(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3),  &
                1:iobnum,k_sta:k_end)
complex(8) :: tpsi(mg_sta(1)-Nd:mg_end(1)+Nd,mg_sta(2)-Nd:mg_end(2)+Nd,mg_sta(3)-Nd:mg_end(3)+Nd)
complex(8) :: htpsi(mg_sta(1):mg_end(1),mg_sta(2):mg_end(2),mg_sta(3):mg_end(3))
real(8) :: esp3(itotMST,num_kpoints_rd)
integer :: iob,iik
integer :: ix,iy,iz
complex(8) :: cbox
integer :: iob_allob



iwk_size=2
call make_iwksta_iwkend

ihpsieff=0

esp3=0.d0

!$OMP parallel do private(iz,iy,ix) 
do iz=mg_sta(3)-Nd,mg_end(3)+Nd
do iy=mg_sta(2)-Nd,mg_end(2)+Nd
do ix=mg_sta(1)-Nd,mg_end(1)+Nd
  tpsi(ix,iy,iz)=0.d0
end do
end do
end do



do iik=k_sta,k_end
do iob=1,iobnum
  call calc_allob(iob,iob_allob,iparaway_ob,itotmst,mst,iobnum)


!$OMP parallel do private(iz,iy,ix) 
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    tpsi(ix,iy,iz)=psi_in(ix,iy,iz,iob,iik)
  end do
  end do
  end do

  call hpsi2(tpsi,htpsi,iob_allob,iik,0,0)




  cbox=0.d0
!$OMP parallel do reduction ( + : cbox ) private(iz,iy,ix) 
  do iz=mg_sta(3),mg_end(3)
  do iy=mg_sta(2),mg_end(2)
  do ix=mg_sta(1),mg_end(1)
    cbox=cbox+conjg(psi_in(ix,iy,iz,iob,iik))*htpsi(ix,iy,iz)
  end do
  end do
  end do
  
  esp3(iob_allob,iik)=dble(cbox)*Hvol



end do
end do


call comm_summation(esp3,esp,itotMST*num_kpoints_rd,nproc_group_global)




end subroutine total_energy_periodic_scf_esp
