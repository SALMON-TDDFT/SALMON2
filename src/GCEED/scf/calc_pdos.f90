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
subroutine calc_pdos(lg,info,pp)
use structures
use salmon_parallel, only: nproc_id_global
use salmon_communication, only: comm_is_root, comm_summation
use inputoutput, only: out_dos_start, out_dos_end, out_dos_function, &
                       out_dos_width, out_dos_nenergy, yn_out_dos_set_fe_origin, uenergy_from_au
use prep_pp_sub, only: bisection
use calc_allob_sub
use scf_data
implicit none
type(s_rgrid)           ,intent(in) :: lg
type(s_orbital_parallel),intent(in) :: info
type(s_pp_info)         ,intent(in) :: pp
integer :: iob,iobmax,iob_allob,iatom,L,ix,iy,iz,iik
integer :: ikoa
integer :: intr
real(8) :: phi_r
real(8) :: rr
real(8) :: ratio1,ratio2
real(8) :: xx,yy,zz
real(8) :: Ylm
integer :: lm
real(8) :: rbox_pdos(25,MI)
real(8) :: rbox_pdos2(25,MI)
real(8) :: pdos_l_tmp(out_dos_nenergy,0:4,MI)
real(8) :: pdos_l(out_dos_nenergy,0:4,MI)
character(100) :: Outfile
real(8) :: fk,ww,dw
integer :: iw
real(8) :: ene_homo,ene_lumo,ene_min,ene_max,efermi,eshift

if( all(pp%upp_f==0.0d0) )then
  write(*,*) "@calc_pdos: Pseudoatom wave function is not available"
  return
end if

call calc_pmax(iobmax)

ene_min = minval(esp(:,:))
ene_max = maxval(esp(:,:))
if(yn_out_dos_set_fe_origin=='y'.and.nstate>nelec/2) then 
  ene_homo = esp(nelec/2,1)
  ene_lumo = esp(nelec/2+1,1)
  efermi = (ene_homo+ene_lumo)*0.5d0 
  eshift = efermi 
else 
  eshift = 0d0 
endif 
out_dos_start = max(out_dos_start,ene_min-0.25d0*(ene_max-ene_min))
out_dos_end = min(out_dos_end,ene_max+0.25d0*(ene_max-ene_min))
dw=(out_dos_end-out_dos_start)/dble(out_dos_nenergy-1) 

pdos_l_tmp=0.d0

do iik=k_sta,k_end
do iob=1,iobmax
  call calc_allob(iob,info,iob_allob,itotmst,mst,iobnum)
  rbox_pdos=0.d0
  do iatom=1,MI
    ikoa=Kion(iatom)
    do L=0,Mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        do iz=mg_sta(3),mg_end(3)
        do iy=mg_sta(2),mg_end(2)
        do ix=mg_sta(1),mg_end(1)
          xx=lg%coordinate(ix,1)-Rion(1,iatom)
          yy=lg%coordinate(iy,2)-Rion(2,iatom)
          zz=lg%coordinate(iz,3)-Rion(3,iatom)
          rr=sqrt(xx**2+yy**2+zz**2)+1.d-50
          call bisection(rr,intr,ikoa,pp%nrmax,pp%rad)
          if(intr==1) intr=2
          ratio1=(rr-pp%rad(intr,ikoa))/(pp%rad(intr+1,ikoa)-pp%rad(intr,ikoa)) ; ratio2=1.d0-ratio1
          phi_r= ratio1*pp%upp_f(intr,Lref(ikoa),ikoa)/rr**(Lref(ikoa)+1)*sqrt((2*Lref(ikoa)+1)/(4*Pi)) +  &
                 ratio2*pp%upp_f(intr-1,Lref(ikoa),ikoa)/rr**(Lref(ikoa)+1)*sqrt((2*Lref(ikoa)+1)/(4*Pi))
                                        !Be carefull for upp(i,l)/vpp(i,l) reffering rad(i+1) as coordinate
          call Ylm_sub(xx,yy,zz,lm,Ylm)
          rbox_pdos(lm,iatom)=rbox_pdos(lm,iatom)+psi(ix,iy,iz,iob,iik)*phi_r*Ylm*Hvol
        end do
        end do
        end do
      end do
    end do
  end do
  call comm_summation(rbox_pdos,rbox_pdos2,25*MI,info%icomm_r) 
  do iatom=1,MI
    ikoa=Kion(iatom)
    do L=0,Mlps(ikoa)
      do lm=L**2+1,(L+1)**2
        select case (out_dos_function)
        case('lorentzian') 
          fk=2.d0*out_dos_width/pi
          do iw=1,out_dos_nenergy 
            ww=out_dos_start+dble(iw-1)*dw+eshift-esp(iob_allob,iik)  
            pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
              +abs(rbox_pdos2(lm,iatom))**2*fk/(ww**2+out_dos_width**2) 
          end do 
        case('gaussian')
          fk=2.d0/(sqrt(2.d0*pi)*out_dos_width)
          do iw=1,out_dos_nenergy 
            ww=out_dos_start+dble(iw-1)*dw+eshift-esp(iob_allob,iik)  
            pdos_l_tmp(iw,L,iatom)=pdos_l_tmp(iw,L,iatom)  &
              +abs(rbox_pdos2(lm,iatom))**2*fk*exp(-(0.5d0/out_dos_width**2)*ww**2) 
          end do
        end select
      end do
    end do
  end do
end do
end do
call comm_summation(pdos_l_tmp,pdos_l,out_dos_nenergy*5*MI,info%icomm_ko) 

if(comm_is_root(nproc_id_global))then
  do iatom=1,MI
    ikoa=Kion(iatom)
    write(fileNumber, '(i8)') iatom
    OutFile = "pdos"//trim(adjustl(fileNumber))//".data"
    open(101,file=OutFile)
    write(101,'("# Projected Density of States")') 
    select case(unit_energy)
    case('au','a.u.')
      if(Mlps(ikoa)==0)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.]")') 
      else if(Mlps(ikoa)==1)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.]")') 
      else if(Mlps(ikoa)==2)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.]")') 
      else if(Mlps(ikoa)==3)then
        write(101,'("# Energy[a.u.] PDOS(l=0)[a.u.] PDOS(l=1)[a.u.] PDOS(l=2)[a.u.] PDOS(l=3)[a.u.]")')
      end if 
    case('ev','eV')
      if(Mlps(ikoa)==0)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV]")')
      else if(Mlps(ikoa)==1)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV]")')
      else if(Mlps(ikoa)==2)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV]")') 
      else if(Mlps(ikoa)==3)then
        write(101,'("# Energy[eV]  PDOS(l=0)[1/eV] PDOS(l=1)[1/eV] PDOS(l=2)[1/eV] PDOS(l=3)[1/eV]")')
      end if 
    end select
    write(101,'("#-----------------------")') 
    if(Mlps(ikoa)==0)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==1)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,2f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==2)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,3f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    else if(Mlps(ikoa)==3)then
      do iw=1,out_dos_nenergy 
        ww=out_dos_start+dble(iw-1)*dw+eshift
        write(101,'(f10.5,4f14.8)') ww*uenergy_from_au,(pdos_l(iw,L,iatom)/uenergy_from_au,L=0,Mlps(ikoa))
      end do
    end if
    close(101)
  end do
end if

end subroutine calc_pdos
