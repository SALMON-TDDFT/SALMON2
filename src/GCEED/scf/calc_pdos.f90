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
subroutine calc_pdos(lg,mg,system,info,pp,energy,tpsi)
use structures
use math_constants, only: pi,zi
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use salmon_global, only: out_dos_start, out_dos_end, out_dos_function, &
out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, nelec, kion, natom, nstate, unit_energy
use inputoutput, only: uenergy_from_au
use prep_pp_sub, only: bisection
implicit none
type(s_rgrid)           ,intent(in) :: lg,mg
type(s_dft_system)      ,intent(in) :: system
type(s_orbital_parallel),intent(in) :: info
type(s_pp_info)         ,intent(in) :: pp
type(s_dft_energy)      ,intent(in) :: energy
type(s_orbital)         ,intent(in) :: tpsi
!
integer :: iob,iatom,L,ix,iy,iz,iik,ispin
integer :: ikoa
integer :: intr
real(8) :: phi_r
real(8) :: rr
real(8) :: ratio1,ratio2
real(8) :: xx,yy,zz
real(8) :: Ylm
integer :: lm
real(8) :: rbox_pdos(25,natom)
real(8) :: rbox_pdos2(25,natom)
real(8) :: pdos_l_tmp(out_dos_nenergy,0:4,natom)
real(8) :: pdos_l(out_dos_nenergy,0:4,natom)
character(100) :: Outfile
real(8) :: fk,ww,dw
integer :: iw
real(8) :: ene_homo,ene_lumo,ene_min,ene_max,efermi,eshift
character(20) :: fileNumber

if( all(pp%upp_f==0.0d0) )then
  write(*,*) "@calc_pdos: Pseudoatom wave function is not available"
  return
end if

ene_min = minval(energy%esp(:,:,:))
ene_max = maxval(energy%esp(:,:,:))
if(yn_out_dos_set_fe_origin=='y'.and.nstate>nelec/2) then 
  ene_homo = energy%esp(nelec/2,1,1)
  ene_lumo = energy%esp(nelec/2+1,1,1)
  efermi = (ene_homo+ene_lumo)*0.5d0 
  eshift = efermi 
else 
  eshift = 0d0 
endif 
out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min))
out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min))
dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1) 

pdos_l_tmp=0.d0

do ispin=1,system%nspin
do iik=info%ik_s,info%ik_e
do iob=info%io_s,info%io_e
  rbox_pdos=0.d0
  do iatom=1,natom
    ikoa=Kion(iatom)
    do L=0,pp%mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        do iz=mg%is(3),mg%ie(3)
        do iy=mg%is(2),mg%ie(2)
        do ix=mg%is(1),mg%ie(1)
          xx=lg%coordinate(ix,1)-system%Rion(1,iatom)
          yy=lg%coordinate(iy,2)-system%Rion(2,iatom)
          zz=lg%coordinate(iz,3)-system%Rion(3,iatom)
          rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
          call bisection(rr,intr,ikoa,pp%nrmax,pp%rad)
          if(intr==1) intr=2
          ratio1=(rr-pp%rad(intr,ikoa))/(pp%rad(intr+1,ikoa)-pp%rad(intr,ikoa)) ; ratio2=1.d0-ratio1
          phi_r= ratio1*pp%upp_f(intr,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi)) +  &
                 ratio2*pp%upp_f(intr-1,pp%lref(ikoa),ikoa)/rr**(pp%lref(ikoa)+1)*sqrt((2*pp%lref(ikoa)+1)/(4*Pi))
                                        !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
          call Ylm_sub(xx,yy,zz,lm,Ylm)
          rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+tpsi%zwf(ix,iy,iz,ispin,iob,iik,1)*phi_r*Ylm*system%Hvol
        end do
        end do
        end do
      end do
    end do
  end do
  call comm_summation(rbox_pdos,rbox_pdos2,25*natom,info%icomm_r)
  do iatom=1,natom
    ikoa=Kion(iatom)
    do L=0,pp%mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        select case (out_dos_function)
        case('lorentzian') 
          fk=2.d0*out_dos_width/pi
          do iw=1,out_dos_nenergy 
            ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
            pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
              +abs(rbox_pdos2(lm,iatom))**2*fk/(ww**2+out_dos_width**2) 
          end do 
        case('gaussian')
          fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
          do iw=1,out_dos_nenergy 
            ww=out_dos_start+dble(iw-1)*dw+eshift-energy%esp(iob,iik,ispin)
            pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
              +abs(rbox_pdos2(lm,iatom))**2*fk*exp(-(0.5d0/out_dos_width**2)*ww**2) 
          end do
        end select
      end do
    end do
  end do
end do
end do
end do
call comm_summation(pdos_l_tmp,pdos_l,out_dos_nenergy*5*natom,info%icomm_ko)

if(comm_is_root(nproc_id_global))then
  do iatom=1,natom
    ikoa=Kion(iatom)
    write(fileNumber, '(i8)') iatom
    OutFile = "pdos"//trim(adjustl(fileNumber))//".data"
    open(101,file=OutFile)
    write(101,'("# Projected Density of States")') 
    select case(unit_energy)
    case('au','a.u.')
      if(pp%mlps(ikoa)==0)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.]")') 
      else if(pp%mlps(ikoa)==1)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.]")') 
      else if(pp%mlps(ikoa)==2)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.]")') 
      else if(pp%mlps(ikoa)==3)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.] PDOS(l=3)[a.u.]")')
      end if 
    case('ev','eV')
      if(pp%mlps(ikoa)==0)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV]")')
      else if(pp%mlps(ikoa)==1)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV]")')
      else if(pp%mlps(ikoa)==2)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV]")') 
      else if(pp%mlps(ikoa)==3)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV] PDOS(l=3)[1/eV]")')
      end if 
    end select
    write(101,'("#-----------------------")') 
    if(pp%mlps(ikoa)==0)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
      end do
    else if(pp%mlps(ikoa)==1)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,2f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
      end do
    else if(pp%mlps(ikoa)==2)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,3f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
      end do
    else if(pp%mlps(ikoa)==3)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,4f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,pp%mlps(ikoa))
      end do
    end if
    close(101)
  end do
end if

end subroutine calc_pdos
