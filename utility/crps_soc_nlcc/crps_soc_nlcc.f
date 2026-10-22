c  Program for Soft Pseudo-potential
c  with generalization by Hamann
c
c    References: N. Troullier and J.L. Martins,
c                  Phys. Rev. B43(1991)1993.
c                D.R. Hamann, Phys. Rev. B40(1989)2980.
c                S.G. Louie et.al, Phys. Rev. B26(1982)1738.
c                  (partial core correction)
c
c    1. Solve all-electron Dirac LDA problem (sub. All_Electron)
c    2. Construct pseudopotential (sub. Pseudo_Potential)
c       Pseudo-potential is obtained.
c    3. Solve pseudo-atom LDA problem (sub. Psuedo_Atom)
c       to construct pseudo-atom wave function necessary for
c       Kleinman-Bylander treatment of pseudo-potential in LDA.
c
c    To specify element, the following parameters in the main program
c    need to be given.
c
c       ZA       Atomic number of the element
c       NA       Number of orbitals in Dirac LDA calculation
c       occA(n)  occupation number in Dirac LDA calculation
c       Lmax     Maximum angular momentum of pseudo-potential
c       PtoA(m)=n
c                m-th orbital in pseudo-atom is n-th orbital
c                is n-th orbital in all-electron calculation
c       irRC(m)  pseudo-potential's connection radius
c                (int(radius*100.) for hr=0.01\AA)
c       irPC     radius of partial core (int(radius*100.))
c       ZPA      charge number of pseudo-atom
c       occPA(l) occupation number in pseudo-atom
c
      use crps_output, only: crps_write_v1
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      implicit none

      common /asy/ZA,NA,occA
      common /pob/lP,PtoA,AtoP
      common /psy/ZPA,Lmax,NP,occPA
      common /con/irRC,irPC
      common /lcr/x0,xN,hx,Rx
      common /par/rhoPA,VPA,spePA,uP
      common /vpp/VCP,VLSP,rhoPC
c parameter
      integer NAdim,NPdim,Ldim,Nx,Nr
      parameter(NAdim=26,NPdim=7,Ldim=3)
      parameter(Nx=640,Nr=2000)
      real*8 hr
      parameter(hr=0.01d0)
c common
      integer NA
      real*8 ZA,occA(NAdim)
      integer lP(NPdim),PtoA(NPdim),AtoP(NAdim)
      integer Lmax,NP
      real*8 ZPA,occPA(0:Ldim)
      integer irRC(NPdim),irPC
      real*8 x0,xN,hx,Rx(0:2*Nx)
      real*8 rhoPA(0:Nr),VPA(0:Nr),spePA(0:Ldim),uP(0:Nr,0:Ldim)
      real*8 VCP(0:Nr,0:Ldim),VLSP(0:Nr,0:Ldim),rhoPC(0:Nr)
c local
      integer ix,ir,mP,l,ios,jtwo(NPdim),lloc,nset,first,last
      character*256 config,output_file
      real*8 core_radius,c2lo(NPdim),c2hi(NPdim)
      real*8 cutoff_radius_bohr(NPdim),default_cutoff(NPdim)
      real*8,parameter::unset=-huge(1.d0)
      integer nodeA(NAdim),lA(NAdim),jA(NAdim),kA(NAdim)
      common /aob/nodeA,lA,jA,kA
      real*8 jpot(0:Nr,NPdim),jwave(0:Nr,NPdim)
      real*8 joc(NPdim),jenergy(NPdim),vlocal(0:Nr)
      logical write_legacy
      namelist /crps/ core_radius,c2lo,c2hi,output_file,lloc,
     & write_legacy,ZA,NA,occA,PtoA,cutoff_radius_bohr,Lmax

      interface
         subroutine Pseudo_Potential(c2lo,c2hi,jpot,jwave,jenergy)
            import :: NPdim,Nr
            implicit none
            real(8),intent(in)::c2lo(NPdim),c2hi(NPdim)
            real(8),intent(out)::jpot(0:Nr,NPdim),jwave(0:Nr,NPdim)
            real(8),intent(out)::jenergy(NPdim)
         end subroutine Pseudo_Potential
      end interface

      AtoP=0
      call init

c log coordinate
      x0=-10.d0
      xN=4.d0
      hx=(xN-x0)/Nx
      do 10 ix=0,2*Nx
         Rx(ix)=exp(x0+ix*hx/2.d0)
 10   continue

c all electron atom
      ZA=79.d0
c number of orbitals
      NA=26

c 1s1/2 (2)
      occA(1)=2.d0
c 2s1/2
      occA(2)=2.d0
c 2p1/2
      occA(3)=2.d0
c 2p3/2 (10)
      occA(4)=4.d0
c 3s1/2
      occA(5)=2.d0
c 3p1/2
      occA(6)=2.d0
c 3p3/2 (18)
      occA(7)=4.d0
c 3d3/2
      occA(8)=4.d0
c 3d5/2 (28)
      occA(9)=6.d0
c 4s1/2
      occA(10)=2.d0
c 4p1/2
      occA(11)=2.d0
c 4p3/2
      occA(12)=4.d0
c 4d3/2
      occA(13)=4.d0
c 4d5/2
      occA(14)=6.d0
c 4f5/2
      occA(15)=6.d0
c 4f7/2 (60)
      occA(16)=8.d0
c 5s1/2
      occA(17)=2.d0
c 5p1/2
      occA(18)=2.d0
c 5p3/2
      occA(19)=4.d0
c 5d3/2
      occA(20)=4.d0
c 5d5/2 (78)
      occA(21)=6.d0
c 6s1/2
      occA(22)=1.d0
c 6p1/2
      occA(23)=0.d0
c 6p3/2
      occA(24)=0.d0
c 5f5/2
      occA(25)=0.d0
c 5f7/2
      occA(26)=0.d0

c pseudo potential and atom
c the same occupation is assumed as that of all electron
      Lmax=3
      NP=2*Lmax+1

c  s1/2
      lP(1)=0
      PtoA(1)=22
      AtoP(22)=1
      irRC(1)=int(2.35d0*0.529d0*100.d0)
c  p1/2
      lP(2)=1
      PtoA(2)=23
      AtoP(23)=2
      irRC(2)=int(2.35d0*0.529d0*100.d0)
c  p3/2
      lP(3)=1
      PtoA(3)=24
      AtoP(24)=3
      irRC(3)=int(2.35d0*0.529d0*100.d0)
c  d3/2
      lP(4)=2
      PtoA(4)=20
      AtoP(20)=4
      irRC(4)=int(2.35d0*0.529d0*100.d0)
c  d5/2
      lP(5)=2
      PtoA(5)=21
      AtoP(21)=5
      irRC(5)=int(2.35d0*0.529d0*100.d0)
c  f5/2
      lP(6)=3
      PtoA(6)=25
      AtoP(25)=6
      irRC(6)=int(2.35d0*0.529d0*100.d0)
c  f7/2
      lP(7)=3
      PtoA(7)=26
      AtoP(26)=7
      irRC(7)=int(2.35d0*0.529d0*100.d0)

c partial core radius
c no partial core if irPC=0
      irPC=0
      core_radius=0.d0
      c2lo=-10.5d0
      c2hi=10.5d0
      output_file='Au_soc_nlcc_rps.dat'
      write_legacy=.false.
      lloc=Lmax
      cutoff_radius_bohr=irRC*hr/0.529177210903d0
      call get_command_argument(1,config)
      if(len_trim(config).gt.0) then
         NA=-1
         Lmax=-1
         lloc=-1
         output_file=''
         default_cutoff=cutoff_radius_bohr
         cutoff_radius_bohr=unset
         open(99,file=trim(config),status='old',iostat=ios)
         if(ios.ne.0) stop 1
         read(99,nml=crps,iostat=ios)
         close(99)
         if(ios.ne.0) stop 1
         nset=count(cutoff_radius_bohr.ne.unset)
         if(Lmax.eq.-1) then
            Lmax=3
            if(nset.ne.0) then
               if(nset.ne.3.and.nset.ne.5.and.nset.ne.7) then
                  write(*,*)'Set Lmax for a partial radius list'
                  stop 1
               endif
               if(any(cutoff_radius_bohr(1:nset).eq.unset))then
                  write(*,*)'Set Lmax for a sparse radius list'
                  stop 1
               endif
               Lmax=(nset-1)/2
            endif
         endif
         where(cutoff_radius_bohr.eq.unset)
            cutoff_radius_bohr=default_cutoff
         endwhere
         if(NA.eq.-1.and.Lmax.ge.1.and.Lmax.le.Ldim)
     &      NA=maxval(PtoA(1:2*Lmax+1))
         if(lloc.eq.-1)lloc=Lmax
         if(len_trim(output_file).eq.0) then
            first=1
            do ix=1,len_trim(config)
               if(config(ix:ix).eq.'/')first=ix+1
            enddo
            last=len_trim(config)
            do ix=last,first,-1
               if(config(ix:ix).eq.'.') then
                  last=ix-1
                  exit
               endif
            enddo
            if(last.lt.first.or.last-first.gt.220)stop 1
            output_file=config(first:last)//'_soc_nlcc_rps.dat'
         endif
      endif
      if(Lmax.lt.1.or.Lmax.gt.Ldim) stop 1
      NP=2*Lmax+1
      if(.not.ieee_is_finite(ZA)) stop 1
      if(ZA.lt.1.d0.or.ZA.gt.118.d0.or.NA.lt.1.or.NA.gt.NAdim)
     & stop 1
      if(.not.all(ieee_is_finite(occA))) stop 1
      if(any(occA.lt.0.d0)) stop 1
      if(abs(sum(occA(1:NA))-ZA).gt.1.d-8) stop 1
      if(.not.all(ieee_is_finite(cutoff_radius_bohr))) stop 1
      irRC=nint(cutoff_radius_bohr*0.529177210903d0/hr)
      if(any(irRC.lt.3).or.any(irRC.gt.Nr-3)) stop 1
      AtoP=0
      do mP=1,NP
         if(PtoA(mP).lt.1.or.PtoA(mP).gt.NA) stop 1
         if(AtoP(PtoA(mP)).ne.0) stop 1
         if(lA(PtoA(mP)).ne.lP(mP)) stop 1
         if(occA(PtoA(mP)).gt.jA(PtoA(mP))) stop 1
         AtoP(PtoA(mP))=mP
      enddo
      write(*,*) 'SPECIES_INPUT_VALID ',ZA,NA
      if(.not.ieee_is_finite(core_radius)) stop 1
      if(.not.all(ieee_is_finite(c2lo))) stop 1
      if(.not.all(ieee_is_finite(c2hi))) stop 1
      if(core_radius.lt.0.d0.or.core_radius.ge.(Nr-2)*hr)
     &   stop 1
      irPC=nint(core_radius/hr)
      if(core_radius.gt.0.d0.and.irPC.lt.2) stop 1
      if(lloc.lt.0.or.lloc.gt.Lmax) stop 1
      do mP=1,NP
         if(c2lo(mP).ge.c2hi(mP)) stop 1
         joc(mP)=occA(PtoA(mP))
         jtwo(mP)=2*lP(mP)-1
         if(mod(mP,2).eq.1) jtwo(mP)=2*lP(mP)+1
      enddo

c all electron calculation
      call All_Electron

      call Pseudo_Potential(c2lo,c2hi,jpot,jwave,jenergy)

c Full j-resolved interchange before scalar pseudo-atom iteration.
      ZPA=sum(joc(1:NP))
c Local reference must retain Coulomb tail, even for empty channels.
      if(lloc.eq.0) then
         vlocal=jpot(:,1)
      else
         vlocal=(lloc*jpot(:,2*lloc)
     &          +(lloc+1)*jpot(:,2*lloc+1))/(2*lloc+1)
      endif
      write(*,*) 'CORE_RANGE ',minval(rhoPC),maxval(rhoPC),
     & minloc(rhoPC)
      call crps_write_v1(trim(output_file),Nr,hr,Lmax,ZPA,
     & NP,irPC,lloc,lP,jtwo,irRC,joc,jenergy,rhoPC,vlocal,
     & jpot,jwave)
      write(*,*) 'CRPS_SOC_NLCC_WRITTEN ',trim(output_file)
      if(.not.write_legacy) stop
c non-relativisitic pseudo-atom calculation
c
      ZPA=sum(joc(1:NP))
      occPA=0.d0
      do mP=1,NP
         occPA(lP(mP))=occPA(lP(mP))+joc(mP)
      enddo

      call Pseudo_Atom

c write pseudo-potential file

      open(7,file='Au_crps.dat',status='new')
      write(7,1000) Nr,hr,Lmax,ZPA
      write(7,1010) irPC,(irRC(mP)*hr,mP=1,NP,2)
      do 100 ir=0,Nr
         write(7,1020) ir*hr,rhoPC(ir),(VCP(ir,l),l=0,Lmax)
 100  continue
      do 110 ir=0,Nr
         write(7,1030) ir*hr,(uP(ir,l),l=0,Lmax)
 110  continue
      close(7)

 1000 format(1x,i12,f12.5,i12,f12.5)
 1010 format(1x,i12,4f12.5)
 1020 format(1x,f10.4,f14.5,4f13.5)
 1030 format(1x,f10.4,14x,4f13.5)

      stop
      end

C---------------------------------------------------------------------
      SUBROUTINE All_Electron
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /aob/nodeA,lA,jA,kA
      common /asy/ZA,NA,occA
      common /lcr/x0,xN,hx,Rx
      common /aer/rhoA,VA,speA,FA,GA
c parameter
      integer NAdim,Nx,Niter
      parameter(NAdim=26)
      parameter(Nx=640)
      parameter(Niter=100)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer nodeA(NAdim),lA(NAdim),jA(NAdim),kA(NAdim)
      integer NA
      real*8 ZA,occA(NAdim)
      real*8 x0,xN,hx,Rx(0:2*Nx)
      real*8 rhoA(0:2*Nx),VA(0:2*Nx),speA(NAdim)
      real*8 FA(0:Nx,NAdim),GA(0:Nx,NAdim)
c local
      integer mA,ix,jx,iter
      real*8 V_xc,E_xc,Etot,sum,energy,rhoA_old(0:2*Nx)

c initial density
      do 10 ix=0,Nx
         rhoA(2*ix)=ZA*(0.5d0/Pi)**1.5d0*exp(-0.5d0*Rx(2*ix)**2)
 10   continue

c start all electron iteration

      do 20 iter=1,Niter
         write(*,*) 'iter=',iter

c interpolate density
      rhoA(1)=0.375d0*rhoA(0)+0.75d0*rhoA(2)-0.125d0*rhoA(4)
      rhoA(2*Nx-1)=-0.125d0*rhoA(2*Nx-4)+0.75d0*rhoA(2*Nx-2)
     &             +0.375d0*rhoA(2*Nx)
      do 30 ix=1,Nx-2
         rhoA(2*ix+1)=(-rhoA(2*ix-2)+9*rhoA(2*ix)
     &                 +9*rhoA(2*ix+2)-rhoA(2*ix+4))/16d0
 30   continue

c calculate Hartree potential
      do 40 ix=0,2*Nx
         sum=0.d0
         do 50 jx=0,ix
            sum=sum+Rx(jx)**3*rhoA(jx)
 50      continue
         sum=sum/Rx(ix)
         do 60 jx=ix+1,2*Nx
            sum=sum+Rx(jx)**2*rhoA(jx)
 60      continue
         VA(ix)=sum*4*Pi*E2*(hx/2)
         VA(ix)=VA(ix)+V_xc(rhoA(ix))
 40   continue

c Total energy calculation
      Etot=0.d0
      do 70 mA=1,NA
         Etot=Etot+occA(mA)*speA(mA)
 70   continue
      sum=0.d0
      do 80 ix=0,2*Nx
         sum=sum+Rx(ix)**3*rhoA(ix)
     &           *(E_xc(rhoA(ix))-VA(ix)/2-V_xc(rhoA(ix))/2)
 80   continue
      sum=sum*(hx/2)*4*Pi
      Etot=Etot+sum
      write(*,*) 'total energy = ',Etot

      if(iter.eq.Niter) go to 130

      do 90 ix=0,2*Nx
         rhoA_old(ix)=rhoA(ix)
         rhoA(ix)=0.d0
 90   continue

c solve dirac equation
      do 100 mA=1,NA
         if(mA.le.2) energy=-ZA**2*13.6d0/(nodeA(mA)+lA(mA)+1)**2
         call Solve_Dirac(mA,energy)

         write(*,*) mA,nodeA(mA),lA(mA),(jA(mA)-1),'/2   ',energy
         speA(mA)=energy
         if(abs(int(energy+1.d-15)-energy).lt.1d-14.or.
     &      abs(int(energy-1.d-15)-energy).lt.1d-14) speA(mA)=1.d0
         if(occA(mA).gt.1.d-5.and.speA(mA).gt.0.5d0) then
            write(*,*) 'occupied orbitals energy is positive'
c            stop
         endif
         do 110 ix=0,Nx
            rhoA(2*ix)=rhoA(2*ix)+(FA(ix,mA)**2+GA(ix,mA)**2)
     &                            /(4*Pi*Rx(2*ix)**2)*occA(mA)
 110     continue

 100  continue

      if(iter.ne.Niter-1) then
      do 120 ix=0,Nx
         rhoA(2*ix)=0.15d0*rhoA(2*ix)+0.85d0*rhoA_old(2*ix)
 120  continue
      endif

 20   continue

c construct positive energy wave function for unbound orbitals
 130  do 140 mA=1,NA
         if(speA(mA).gt.0.5d0) go to 150
 140  continue
      return

 150  energy=-1d6
      do 160 mA=1,NA
         if(speA(mA).lt.0.d0.and.speA(mA).gt.energy)
     &      energy=speA(mA)
 160  continue

      do 170 mA=1,NA
         if(speA(mA).gt.0.5d0) then
            call Cont_wave(mA,energy)
            speA(mA)=energy
         endif
 170  continue

      return
      end

C---------------------------------------------------------------------
c solve radial Dirac equation by Runge-Kutta method
      SUBROUTINE Solve_Dirac(mA,energy)
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /aob/nodeA,lA,jA,kA
      common /asy/ZA,NA,occA
      common /lcr/x0,xN,hx,Rx
      common /aer/rhoA,VA,speA,FA,GA
c parameter
      integer NAdim,Nx,Nadd
      parameter(NAdim=26)
      parameter(Nx=640)
      parameter(Nadd=80)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer nodeA(NAdim),lA(NAdim),jA(NAdim),kA(NAdim)
      integer NA
      real*8 ZA,occA(NAdim)
      real*8 x0,xN,hx,Rx(0:2*Nx)
      real*8 rhoA(0:2*Nx),VA(0:2*Nx),speA(NAdim)
      real*8 FA(0:Nx,NAdim),GA(0:Nx,NAdim)
c local
      integer mA,kappa,node,ix,node0,Ntp,ncount,iter,Mx
      real*8 eps,eps0,s,energy,rke0,rke1,epslo,epshi
      real*8 FM,GM,ratio,sum,deltaE
      real*8 Fk1,Gk1,Fk2,Gk2,Fk3,Gk3,Fk4,Gk4

      kappa=kA(mA)
      node=nodeA(mA)

      eps=energy
      s=sqrt(kappa**2-ZA**2*ALP**2)

      node0=-1
 10   continue
      if(eps.gt.-0.1d0) then
         write(*,*) 'No solution in Solve_Dirac'
         energy=int(energy)
         return
      endif
c classical turning point
         rke0=H2M*kappa*(kappa+1)/Rx(0)**2
     &        -ZA*E2/Rx(0)+VA(0)-eps
      do 20 ix=1,Nx
         rke1=H2M*kappa*(kappa+1)/Rx(2*ix)**2
     &        -ZA*E2/Rx(2*ix)+VA(2*ix)-eps
         if(rke1.gt.rke0.and.rke1.gt.0.d0)  go to 30
         rke0=rke1
 20   continue
c when classical turning point cannot be found
      eps=eps*0.8d0
      if(eps.gt.-0.1d0) then
         write(*,*) 'failed to find classical turning point'
         stop
      endif
      go to 10

 30   Ntp=ix+Nadd

      if(Ntp.ge.Nx-10) then
         write(*,*) 'too large Ntp'
         stop
      endif

c find approximate energy

      FA(0,mA)=Rx(0)**s
      GA(0,mA) =(s+kappa)/(ZA*ALP)*FA(0,mA)

      ncount=0
      do 40 ix=0,Ntp-1
      Fk1=hx*(-kappa*FA(ix,mA)
     &   +(ZA*ALP+Rx(2*ix)*(2*MC2+eps-VA(2*ix))/HC)*GA(ix,mA))
      Gk1=hx*(kappa*GA(ix,mA)
     &   -(ZA*ALP+Rx(2*ix)*(eps-VA(2*ix))/HC)*FA(ix,mA))
      Fk2=hx*(-kappa*(FA(ix,mA)+Fk1/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk1/2))
      Gk2=hx*(kappa*(GA(ix,mA)+Gk1/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk1/2))
      Fk3=hx*(-kappa*(FA(ix,mA)+Fk2/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk2/2))
      Gk3=hx* (kappa* (GA(ix,mA) +Gk2/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk2/2))
      Fk4=hx*(-kappa*(FA(ix,mA)+Fk3)
     &   +(ZA*ALP+Rx(2*ix+2 )*(2*MC2+eps-VA(2*ix+2))/HC)
     &   *(GA(ix,mA)+Gk3))
      Gk4=hx*(kappa*(GA(ix,mA)+Gk3)
     &   -(ZA*ALP+Rx(2*ix+2)*(eps-VA(2*ix+2))/HC)
     &   *(FA(ix,mA)+Fk3))
      FA(ix+1,mA)=FA(ix,mA)+(Fk1+2*Fk2+2*Fk3+Fk4)/6
      GA(ix+1,mA)=GA(ix,mA)+(Gk1+2*Gk2+2*Gk3+Gk4)/6
      if(FA(ix+1,mA)*FA(ix,mA).lt.0.d0) ncount=ncount+1
 40   continue

      if(node0.eq.-1) then
         if(ncount.gt.node) then
            eps0=eps
            eps=eps*1.1d0
         else
            eps0=eps
            eps=eps*0.9d0
         endif
         node0=ncount
         go to 10
      else
         if(abs(ncount-node0).ge.2) then
            write(*,*) 'ncount changed more than two unit'
            stop
         endif
         if(abs(ncount-node0).eq.1.and.ncount.eq.node) then
            epslo=eps
            epshi=eps0
            go to 50
         endif
         if(abs(ncount-node0).eq.1.and.node0.eq.node) then
            epslo=eps0
            epshi=eps
            go to 50
         endif
         if(ncount.gt.node) then
            eps0=eps
            eps=eps*1.1d0
         else
            eps0=eps
            eps=eps*0.9d0
         endif
         node0=ncount
         go to 10
      endif

 50   continue

c more accurate value by bi-section method (inner part only)

      do 60 iter=1,5
      eps=(epslo+epshi)/2

      FA(0,mA)=Rx(0)**s
      GA(0,mA) =(s+kappa)/(ZA*ALP)*FA(0,mA)

      ncount=0
      do 70 ix=0,Ntp-1
      Fk1=hx*(-kappa*FA(ix,mA)
     &   +(ZA*ALP+Rx(2*ix)*(2*MC2+eps-VA(2*ix))/HC)*GA(ix,mA))
      Gk1=hx*(kappa*GA(ix,mA)
     &   -(ZA*ALP+Rx(2*ix)*(eps-VA(2*ix))/HC)*FA(ix,mA))
      Fk2=hx*(-kappa*(FA(ix,mA)+Fk1/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk1/2))
      Gk2=hx*(kappa*(GA(ix,mA)+Gk1/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk1/2))
      Fk3=hx*(-kappa*(FA(ix,mA)+Fk2/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk2/2))
      Gk3=hx* (kappa* (GA(ix,mA) +Gk2/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk2/2))
      Fk4=hx*(-kappa*(FA(ix,mA)+Fk3)
     &   +(ZA*ALP+Rx(2*ix+2 )*(2*MC2+eps-VA(2*ix+2))/HC)
     &   *(GA(ix,mA)+Gk3))
      Gk4=hx*(kappa*(GA(ix,mA)+Gk3)
     &   -(ZA*ALP+Rx(2*ix+2)*(eps-VA(2*ix+2))/HC)*(FA(ix,mA)+Fk3))
      FA(ix+1,mA)=FA(ix,mA)+(Fk1+2*Fk2+2*Fk3+Fk4)/6
      GA(ix+1,mA)=GA(ix,mA)+(Gk1+2*Gk2+2*Gk3+Gk4)/6
      if(FA(ix+1,mA)*FA(ix,mA).lt.0.d0) ncount=ncount+1
 70   continue

      if(ncount.le.node) then
         epslo=eps
      else
         epshi=eps
      endif
 60   continue

      Ntp=Ntp-Nadd+5

c accurate eigenvalue by newton method
      do 80 iter=1,50

      FA(0,mA)=Rx(0)**s
      GA(0,mA)=(s+kappa)/(ZA*ALP)*FA(0,mA)

      do 90 ix=0,Ntp-1
      Fk1=hx*(-kappa*FA(ix,mA)
     &   +(ZA*ALP+Rx(2*ix)*(2*MC2+eps-VA(2*ix))/HC)*GA(ix,mA))
      Gk1=hx*(kappa*GA(ix,mA)
     &   -(ZA*ALP+Rx(2*ix)*(eps-VA(2*ix))/HC)*FA(ix,mA))
      Fk2=hx*(-kappa*(FA(ix,mA)+Fk1/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk1/2))
      Gk2=hx*(kappa*(GA(ix,mA)+Gk1/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk1/2))
      Fk3=hx*(-kappa*(FA(ix,mA)+Fk2/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk2/2))
      Gk3=hx* (kappa* (GA(ix,mA) +Gk2/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk2/2))
      Fk4=hx*(-kappa*(FA(ix,mA)+Fk3)
     &   +(ZA*ALP+Rx(2*ix+2 )*(2*MC2+eps-VA(2*ix+2))/HC)
     &   *(GA(ix,mA)+Gk3))
      Gk4=hx*(kappa*(GA(ix,mA)+Gk3)
     &   -(ZA*ALP+Rx(2*ix+2)*(eps-VA(2*ix+2))/HC)*(FA(ix,mA)+Fk3))
      FA(ix+1,mA)=FA(ix,mA)+(Fk1+2*Fk2+2*Fk3+Fk4)/6
      GA(ix+1,mA)=GA(ix,mA)+(Gk1+2*Gk2+2*Gk3+Gk4)/6
 90   continue
      FM=FA(Ntp,mA)
      GM=GA(Ntp,mA)

      Mx=Nx
 100  FA(Mx,mA)=exp(-sqrt(-eps/H2M)*Rx(2*Mx))
      GA(Mx,mA)=(-sqrt(-eps/H2M)+kappa/Rx(2*Mx))*FA(Mx,mA)
     &       /((2*MC2+eps)/HC)
      if(abs(FA(Mx,mA)).lt.1.d-50) then
         Mx=Mx-10
         go to 100
      else
         do 110 ix=Nx,Mx+1,-1
            FA(ix,mA)=0.d0
            GA(ix,mA)=0.d0
 110     continue
      endif

      do 120 ix=Mx,Ntp+1,-1
      Fk1=-hx*(-kappa*FA(ix,mA)
     &   +(ZA*ALP+Rx(2*ix)*(2*MC2+eps-VA(2*ix))/HC)*GA(ix,mA))
      Gk1=-hx*(kappa*GA(ix,mA)
     &   -(ZA*ALP+Rx(2*ix)*(eps-VA(2*ix))/HC)*FA(ix,mA))
      Fk2=-hx*(-kappa*(FA(ix,mA)+Fk1/2)
     &   +(ZA*ALP+Rx(2*ix-1)*(2*MC2+eps-VA(2*ix-1))/HC)
     &   *(GA(ix,mA)+Gk1/2))
      Gk2=-hx*(kappa*(GA(ix,mA)+Gk1/2)
     &   -(ZA*ALP+Rx(2*ix-1)*(eps-VA(2*ix-1))/HC)
     &   *(FA(ix,mA)+Fk1/2))
      Fk3=-hx*(-kappa*(FA(ix,mA)+Fk2/2)
     &   +(ZA*ALP+Rx(2*ix-1)*(2*MC2+eps-VA(2*ix-1))/HC)
     &   *(GA(ix,mA)+Gk2/2))
      Gk3=-hx*(kappa* (GA(ix,mA) +Gk2/2)
     &   -(ZA*ALP+Rx(2*ix-1)*(eps-VA(2*ix-1))/HC)
     &   *(FA(ix,mA)+Fk2/2))
      Fk4=-hx*(-kappa*(FA(ix,mA)+Fk3)
     &   +(ZA*ALP+Rx(2*ix-2)*(2*MC2+eps-VA(2*ix-2))/HC)
     &   *(GA(ix,mA)+Gk3))
      Gk4=-hx*(kappa*(GA(ix,mA)+Gk3)
     &   -(ZA*ALP+Rx(2*ix-2)*(eps-VA(2*ix-2))/HC)*(FA(ix,mA)+Fk3))
      FA(ix-1,mA)=FA(ix,mA)+(Fk1+2*Fk2+2*Fk3+Fk4)/6
      GA(ix-1,mA)=GA(ix,mA)+(Gk1+2*Gk2+2*Gk3+Gk4)/6
 120  continue

      ratio=FM/FA(Ntp,mA)
      do 130 ix=Ntp,Nx
         FA(ix,mA)=FA(ix,mA)*ratio
         GA(ix,mA)=GA(ix,mA)*ratio
 130  continue

      sum=0.d0
      do 140 ix=0,Nx
         sum=sum+Rx(2*ix)*(FA(ix,mA)**2+GA(ix,mA)**2)
 140  continue
      sum=sum*hx

      deltaE=FA(Ntp,mA)*(GA(Ntp,mA)-GM)/sum*HC
      if(abs(deltaE).lt.1d-5) then
         energy=eps
         do 150 ix=0,Nx
            FA(ix,mA)=FA(ix,mA)/sqrt(sum)
            GA(ix,mA)=GA(ix,mA)/sqrt(sum)
 150     continue
         return
      endif

      eps=eps-deltaE

 80   continue
      write(*,*) 'Solve iteration did not converge'
      stop
      end

C---------------------------------------------------------------------
c for waves without bound state
      SUBROUTINE Cont_wave(mA,energy)
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /aob/nodeA,lA,jA,kA
      common /asy/ZA,NA,occA
      common /pob/lP,PtoA,AtoP
      common /con/irRC,irPC
      common /lcr/x0,xN,hx,Rx
      common /aer/rhoA,VA,speA,FA,GA
c parameter
      integer NAdim,NPdim,Nx,Nadd
      parameter(NAdim=26,NPdim=7)
      parameter(Nx=640)
      parameter(Nadd=80)
      real*8 hr
      parameter(hr=0.01d0)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer nodeA(NAdim),lA(NAdim),jA(NAdim),kA(NAdim)
      integer NA
      real*8 ZA,occA(NAdim)
      integer lP(NPdim),PtoA(NPdim),AtoP(NAdim)
      integer irRC(NPdim),irPC
      real*8 x0,xN,hx,Rx(0:2*Nx)
      real*8 rhoA(0:2*Nx),VA(0:2*Nx),speA(NAdim)
      real*8 FA(0:Nx,NAdim),GA(0:Nx,NAdim)
c local
      integer mA,kappa,node,ix,Ntp
      real*8 eps,s,energy,sum
      real*8 Fk1,Gk1,Fk2,Gk2,Fk3,Gk3,Fk4,Gk4

      kappa=kA(mA)
      node=nodeA(mA)

      eps=energy
      s=sqrt(kappa**2-ZA**2*ALP**2)

      FA(0,mA)=Rx(0)**s
      GA(0,mA) =(s+kappa)/(ZA*ALP)*FA(0,mA)

c  integrate up to 2.5 times connection radius
      Ntp=int((log(irRC(AtoP(mA))*hr*2.5d0)-x0)/hx)

      do 10 ix=0,Ntp-1
      Fk1=hx*(-kappa*FA(ix,mA)
     &   +(ZA*ALP+Rx(2*ix)*(2*MC2+eps-VA(2*ix))/HC)*GA(ix,mA))
      Gk1=hx*(kappa*GA(ix,mA)
     &   -(ZA*ALP+Rx(2*ix)*(eps-VA(2*ix))/HC)*FA(ix,mA))
      Fk2=hx*(-kappa*(FA(ix,mA)+Fk1/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk1/2))
      Gk2=hx*(kappa*(GA(ix,mA)+Gk1/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk1/2))
      Fk3=hx*(-kappa*(FA(ix,mA)+Fk2/2)
     &   +(ZA*ALP+Rx(2*ix+1)*(2*MC2+eps-VA(2*ix+1))/HC)
     &   *(GA(ix,mA)+Gk2/2))
      Gk3=hx* (kappa* (GA(ix,mA) +Gk2/2)
     &   -(ZA*ALP+Rx(2*ix+1)*(eps-VA(2*ix+1))/HC)
     &   *(FA(ix,mA)+Fk2/2))
      Fk4=hx*(-kappa*(FA(ix,mA)+Fk3)
     &   +(ZA*ALP+Rx(2*ix+2 )*(2*MC2+eps-VA(2*ix+2))/HC)
     &   *(GA(ix,mA)+Gk3))
      Gk4=hx*(kappa*(GA(ix,mA)+Gk3)
     &   -(ZA*ALP+Rx(2*ix+2)*(eps-VA(2*ix+2))/HC)
     &   *(FA(ix,mA)+Fk3))
      FA(ix+1,mA)=FA(ix,mA)+(Fk1+2*Fk2+2*Fk3+Fk4)/6
      GA(ix+1,mA)=GA(ix,mA)+(Gk1+2*Gk2+2*Gk3+Gk4)/6
 10   continue

      do 20 ix=Ntp+1,Nx
         FA(ix,mA)=0.d0
         GA(ix,mA)=0.d0
 20   continue

      sum=0.d0
      do 30 ix=0,Nx
         sum=sum+Rx(2*ix)*(FA(ix,mA)**2+GA(ix,mA)**2)
 30   continue
      sum=sum*hx

      do 40 ix=0,Nx
         FA(ix,mA)=FA(ix,mA)/sqrt(sum)
         GA(ix,mA)=GA(ix,mA)/sqrt(sum)
 40   continue

      return
      end

C---------------------------------------------------------------------
      SUBROUTINE Pseudo_Potential(c2lo,c2hi,jpot,jwave,jenergy)
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /asy/ZA,NA,occA
      common /pob/lP,PtoA,AtoP
      common /psy/ZPA,Lmax,NP,occPA
      common /con/irRC,irPC
      common /lcr/x0,xN,hx,Rx
      common /aer/rhoA,VA,speA,FA,GA
      common /ppr/rhoP,VP,FP
      common /vpp/VCP,VLSP,rhoPC
      common /ctb/c0,c4,c6,c8,c10,c12,scr_int
c parameter
      integer NAdim,NPdim,Ldim,Nx,Nr
      parameter(NAdim=26,NPdim=7,Ldim=3)
      parameter(Nx=640,Nr=2000)
      real*8 hr
      parameter(hr=0.01d0)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer NA
      real*8 ZA,occA(NAdim)
      integer lP(NPdim),PtoA(NPdim),AtoP(NAdim)
      integer Lmax,NP
      real*8 ZPA,occPA(0:Ldim)
      integer irRC(NPdim),irPC
      real*8 x0,xN,hx,Rx(0:2*Nx)
      real*8 rhoA(0:2*Nx),VA(0:2*Nx),speA(NAdim)
      real*8 FA(0:Nx,NAdim),GA(0:Nx,NAdim)
      real*8 rhoP(0:Nr),VP(0:Nr,NPdim),FP(0:Nr,NPdim)
      real*8 VCP(0:Nr,0:Ldim),VLSP(0:Nr,0:Ldim),rhoPC(0:Nr)
      real*8 c0,c4,c6,c8,c10,c12,scr_int
c local
      integer ix,mP,mA,khi,klo,ir,k,ixc,ns,ic2,jr
      real*8 ya(0:Nx),xa(0:Nx),y2a(0:Nx),x,y
      real*8 sum,r,rho0,drho0,rhs,xhi,xlo,xm,Ca,Cb,rc,xc
      real*8 c2,fc2,c2min,fc2min,c2max,fc2max,c2try,fc2try
      real*8 p,pd,pdd,rr,vhP,vxcP,V_xc
      real*8 VPion(0:Nr,NPdim)

      real(8),intent(in)::c2lo(NPdim),c2hi(NPdim)
      real(8),intent(out)::jpot(0:Nr,NPdim),jwave(0:Nr,NPdim)
      real(8),intent(out)::jenergy(NPdim)
      real(8)::scanlo,scanhi,flo,fhi

      do 5 ix=0,Nx
         xa(ix)=Rx(2*ix)
 5    continue

c core density
      if(irPC.eq.0) then
         do 10 ir=0,Nr
            rhoPC(ir)=0.d0
 10      continue
         go to 20
      endif

c Sum core orbitals directly; rhoA-valence loses positivity in tails.
      do 30 ix=0,Nx
         sum=0.d0
         do 40 mA=1,NA
            if(AtoP(mA).eq.0) then
               sum=sum+(FA(ix,mA)**2+GA(ix,mA)**2)
     &              /(4*Pi*Rx(2*ix)**2)*occA(mA)
            endif
 40      continue
         ya(ix)=sum
 30   continue

      call spline(xa,ya,y2a)
      khi=Nx
      klo=0
      do 50 ir=1,Nr
         r=ir*hr
         x=r
         call splint(xa,ya,y2a,khi,klo,x,y)
c Cubic splines may overshoot below zero in exponentially small tails.
c Use the positive linear interpolant only on those intervals.
         if(y.lt.0.d0) y=(ya(klo)*(xa(khi)-x)
     &                    +ya(khi)*(x-xa(klo)))
     &                    /(xa(khi)-xa(klo))
         rhoPC(ir)=y
 50   continue
      rhoPC(0)=rhoPC(1)

c construct partial core
      rho0=rhoPC(irPC)
      drho0=(rhoPC(irPC+1)-rhoPC(irPC-1))/(2*hr)
      rhs=1.d0+(irPC*hr)*drho0/rho0
      if(rhs.ge.1.d0) then
         write(*,*) 'No solution in partial core'
         stop
      endif

      xhi=Pi-1.d-5
      xlo=1.d-5
      do 60 k=1,30
         xm=(xhi+xlo)/2
         if(xm/tan(xm)-rhs.lt.0.d0) then
            xhi=xm
         else
            xlo=xm
         endif
  60  continue
      write(*,*) xm

      Cb=xm/(irPC*hr)
      Ca=rho0*(irPC*hr)/sin(xm)

      rhoPC(0)=Ca*Cb
      do 70 ir=1,irPC
         rhoPC(ir)=Ca*sin(Cb*ir*hr)/(ir*hr)
 70   continue

 20   continue
c potential interpolation
      do 80 ix=0,Nx
         ya(ix)=VA(2*ix)
 80   continue
      call spline(xa,ya,y2a)
      khi=Nx
      klo=0
      do 90 ir=1,Nr
         r=ir*hr
         x=r
         call splint(xa,ya,y2a,khi,klo,x,y)
         do 100 mP=1,NP
            VP(ir,mP)=-ZA*E2/r+y
 100     continue
 90   continue
      do 110 mP=1,NP
         VP(0,mP)=VP(1,mP)
 110  continue

c pseudo atom wave function
      do 120 mP=1,NP
         mA=PtoA(mP)

         rc=irRC(mP)*hr
         xc=log(rc)
         ixc=int((xc-x0)/hx)

c inner integral
         sum=0.d0
         do 130 ix=0,ixc
            sum=sum+FA(ix,mA)**2*Rx(2*ix)*hx
 130     continue
         sum=sum-0.5d0*FA(ixc,mA)**2*Rx(2*ixc)*hx
     &           +FA(ixc,mA)**2*Rx(2*ixc)*((xc-x0)-ixc*hx)
         scr_int=sum

c wave function interpolation
         do 140 ix=0,Nx
            ya(ix)=FA(ix,mA)
 140     continue
         call spline(xa,ya,y2a)
         khi=Nx
         klo=0
         do 150 ir=1,Nr
            r=ir*hr
            x=r
            call splint(xa,ya,y2a,khi,klo,x,y)
            FP(ir,mP)=y
 150     continue
         FP(0,mP)=FP(1,mP)

c search for c0-c12 coefficients

         ns=0
         do 160 ic2=0,20
            scanlo=c2lo(mP)+(c2hi(mP)-c2lo(mP))*ic2/21.d0
            scanhi=c2lo(mP)+(c2hi(mP)-c2lo(mP))
     &             *(ic2+1)/21.d0
            flo=fc2(mP,scanlo)
            fhi=fc2(mP,scanhi)
            write(*,*) 'C2_SCAN ',mP,scanlo,scanhi,flo,fhi
            if(flo*fhi.lt.0.d0) then
               ns=ns+1
               c2min=scanlo
               c2max=scanhi
               fc2min=flo
               fc2max=fhi
            endif
 160     continue
         if(ns.ne.1) then
            write(*,*) 'C2_RANGE_REJECT channel,roots=',mP,ns
            stop 1
         endif
         do 170 k=1,50
            c2try=(c2min+c2max)/2.d0
            fc2try=fc2(mP,c2try)
            if(fc2min*fc2try.lt.0.d0)  then
               c2max=c2try
               fc2max=fc2try
            else
               c2min=c2try
               fc2min=fc2try
            endif
 170     continue

         c2=c2try
         write(*,*) 'C2_ROOT ',mP,c2,fc2try
         if(abs(fc2try).gt.1.d-9) stop 1

c pseudo-atom wave function and potential
         do 180 ir=1,irRC(mP)
            r=ir*hr
            p=c0+c2*r**2+c4*r**4+c6*r**6+c8*r**8
     &       +c10*r**10+c12*r**12
            pd=2*c2*r+4*c4*r**3+6*c6*r**5
     &        +8*c8*r**7+10*c10*r**9+12*c12*r**11
            pdd=2*c2+12*c4*r**2+30*c6*r**4
     &        +56*c8*r**6+90*c10*r**8+132*c12*r**10
            FP(ir,mP)=r**(lP(mP)+1)*exp(p)
            VP(ir,mP)=speA(mA)+H2M*(2*(lP(mP)+1)/r*pd+pd**2+pdd)
 180     continue
         FP(0,mP)=0.d0
         VP(0,mP)=VP(1,mP)

 120  continue

c pseudo atom density
      do 190 ir=0,Nr
         rhoP(ir)=0.d0
 190  continue
      do 200 mP=1,NP
         mA=PtoA(mP)
      do 200 ir=0,Nr
         rhoP(ir)=rhoP(ir)+FP(ir,mP)**2*occA(mA)
 200  continue
      do 210 ir=1,Nr
         r=ir*hr
         rhoP(ir)=rhoP(ir)/(4*Pi*r**2)
 210  continue
      rhoP(0)=rhoP(1)

c unscreening
      do 220 ir=1,Nr
         r=ir*hr
         sum=0.d0
         do 230 jr=0,ir
            rr=jr*hr
            sum=sum+rr**2*rhoP(jr)/r
 230     continue
         do 240 jr=ir+1,Nr
            rr=jr*hr
            sum=sum+rr*rhoP(jr)
 240     continue
         vhP=sum*hr*4*Pi*E2
         VxcP=V_xc(rhoP(ir)+rhoPC(ir))
         do 250 mP=1,NP
            VPion(ir,mP)=VP(ir,mP)-vhP-VxcP
 250     continue
 220  continue
      do 260 mP=1,NP
         VPion(0,mP)=VPion(1,mP)
 260  continue

c Preserve j-resolved potentials/waves before scalar averaging.
      jpot=VPion
      jwave=FP
      do mP=1,NP
         jenergy(mP)=speA(PtoA(mP))
      enddo

c central and LS pseudo-potentials
      do 270 ir=0,Nr
         VCP(ir,0)=VPion(ir,1)
         VCP(ir,1)=(VPion(ir,2)+2*VPion(ir,3))/3.d0
         VLSP(ir,1)=(VPion(ir,3)-VPion(ir,2))*2/3.d0
         if(NP.ge.5) then
            VCP(ir,2)=(2*VPion(ir,4)+3*VPion(ir,5))/5.d0
            VLSP(ir,2)=(VPion(ir,5)-VPion(ir,4))*2/5.d0
         endif
         if(NP.ge.7) then
            VCP(ir,3)=(3*VPion(ir,6)+4*VPion(ir,7))/7.d0
            VLSP(ir,3)=(VPion(ir,7)-VPion(ir,6))*2.d0/7.d0
         endif
 270  continue

      do 280 mP=1,NP
         mA=PtoA(mP)
         if(speA(mA).gt.-0.01d0) then
            do 290 ir=0,Nr
               VCP(ir,lP(mP))=0.d0
               VLSP(ir,lP(mP))=0.d0
 290        continue
         endif
 280  continue

      return
      end

C---------------------------------------------------------------------
      SUBROUTINE Pseudo_Atom
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /psy/ZPA,Lmax,NP,occPA
      common /con/irRC,irPC
      common /vpp/VCP,VLSP,rhoPC
      common /par/rhoPA,VPA,spePA,uP
c parameter
      integer Ldim,Nr,Niter,NPdim
      parameter(Ldim=3,NPdim=7)
      parameter(Nr=2000)
      parameter(Niter=100)
      real*8 hr
      parameter(hr=0.01d0)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer Lmax,NP,irRC(NPdim),irPC
      real*8 ZPA,occPA(0:Ldim)
      real*8 VCP(0:Nr,0:Ldim),VLSP(0:Nr,0:Ldim),rhoPC(0:Nr)
      real*8 rhoPA(0:Nr),VPA(0:Nr),spePA(0:Ldim),uP(0:Nr,0:Ldim)
c local
      integer ir,iter,jr,l,k,nd,Nmax,Nct
      real*8 r,sum,rr,V_xc,Etot,E_xc,rhoPA_old(0:Nr)
      real*8 emax,emin,etrial,sx0,sx1,sx2

c initial density
      do 10 ir=0,Nr
         r=ir*hr
         rhoPA(ir)=ZPA*(1d0/Pi)**1.5d0*exp(-1d0*r**2)
 10   continue

c start pseudo-atom iteration
      do 20 iter=1,Niter
         write(*,*) 'iter=',iter

c calculate Hartree potential
      do 30 ir=1,Nr
         r=ir*hr
         sum=0.d0
         do 40 jr=0,ir
            rr=jr*hr
            sum=sum+rr**2*rhoPA(jr)
 40      continue
         sum=sum/r
         do 50 jr=ir+1,Nr
            rr=jr*hr
            sum=sum+rr*rhoPA(jr)
 50      continue
         VPA(ir)=sum*4*Pi*E2*hr

         VPA(ir)=VPA(ir)+V_xc(rhoPA(ir)+rhoPC(ir))
 30   continue
         VPA(0)=VPA(1)

c Total energy calculation

      Etot=0.d0
      do 60 l=0,Lmax
         Etot=Etot+occPA(l)*spePA(l)
 60   continue
      sum=0.d0
      do 70 ir=0,Nr
         r=ir*hr
         sum=sum+r**2*rhoPA(ir)
     &                 *(-VPA(ir)/2-V_xc(rhoPA(ir)+rhoPC(ir))/2)
     &          +r**2*(rhoPA(ir)+rhoPC(ir))
     &                 *E_xc(rhoPA(ir)+rhoPC(ir))
 70   continue
      sum=sum*4*Pi*hr
      Etot=Etot+sum
      write(*,*) 'total energy = ',Etot
      write(*,*) (spePA(l),l=0,Lmax)

      do 80 ir=0,Nr
         rhoPA_old(ir)=rhoPA(ir)
         rhoPA(ir)=0.d0
 80   continue

c solve Schroedinger equation

      do 90 l=0,Lmax
c node is zero
         Nmax=Nr

 150     emin=-100.d0
         emax=0.d0
         do 100 k=1,50
            etrial=(emin+emax)/2.d0
            uP(1,l)=hr**(l+1)
            sx0=l*(l+1)/(1*hr)**2
     &         +(VCP(1,l)+VPA(1)-etrial)/H2M
            uP(2,l)=(2*hr)**(l+1)
            sx1=l*(l+1)/(2*hr)**2
     &         +(VCP(2,l)+VPA(2)-etrial)/H2M
            nd=0
            do 110 ir=2,Nmax-1
               sx2=l*(l+1)/((ir+1)*hr)**2
     &            +(VCP(ir+1,l)+VPA(ir+1)-etrial)/H2M
               uP(ir+1,l)=((2+sx1*5.d0/6.d0*hr**2)*uP(ir,l)
     &                       -(1-sx0/12.d0*hr**2)*uP(ir-1,l))
     &                       /(1-sx2/12.d0*hr**2)
               sx0=sx1
               sx1=sx2
              if(abs(uP(ir+1,l)).gt.1d10) then
                 do 120 jr=1,ir+1
                    uP(jr,l)=uP(jr,l)/uP(ir+1,l)
 120             continue
              endif
              if(uP(ir+1,l)*uP(ir,l).lt.0.d0) nd=nd+1
 110       continue
           if(nd.gt.0) then
               emax=etrial
           else
               emin=etrial
           endif
 100     continue

         sum=0.d0
         do 130 ir=1,Nmax
            sum=sum+uP(ir,l)**2*hr
 130     continue
         do 140 ir=1,Nmax
            uP(ir,l)=uP(ir,l)/sqrt(sum)
 140     continue
         if(abs(uP(Nmax-5,l)).gt.abs(uP(Nmax-10,l))) then
            Nmax=Nmax/2
            if(Nmax.lt.200) then
               write(*,*) 'too small Nmax   l=',l
               etrial=1.d0
               spePA(l)=1.d0
               do 145 ir=1,Nr
                  uP(ir,l)=0.d0
 145           continue
               go to 90
            endif
            go to 150
         endif

         uP(0,l)=0.d0
         spePA(l)=etrial

         if(Nmax.ne.Nr) then
            do 160 ir=Nmax+1,Nr
               uP(ir,l)=0.d0
 160        continue
         endif
 90   continue

      do 170 l=0,Lmax
      do 170 ir=1,Nr
         r=ir*hr
         rhoPA(ir)=rhoPA(ir)+uP(ir,l)**2/(4*Pi*r**2)*occPA(l)
 170  continue
      rhoPA(0)=rhoPA(1)

      do 180 ir=0,Nr
         rhoPA(ir)=0.15d0*rhoPA(ir)+0.85d0*rhoPA_old(ir)
 180  continue

 20   continue

      etrial=-1d6
      do 200 l=0,Lmax
      if(spePA(l).lt.0.d0.and.spePA(l).gt.etrial) etrial=spePA(l)
 200  continue
      if(etrial.lt.-1d5) return

      do 210 l=0,Lmax
         Nct=int(irRC(2*l+1)*2.5d0)
         if(spePA(l).lt.0.d0) go to 210
            uP(1,l)=hr**(l+1)
            sx0=l*(l+1)/(1*hr)**2
     &         +(VCP(1,l)+VPA(1)-etrial)/H2M
            uP(2,l)=(2*hr)**(l+1)
            sx1=l*(l+1)/(2*hr)**2
     &         +(VCP(2,l)+VPA(2)-etrial)/H2M
            do 220 ir=2,Nct-1
               sx2=l*(l+1)/((ir+1)*hr)**2
     &            +(VCP(ir+1,l)+VPA(ir+1)-etrial)/H2M
               uP(ir+1,l)=((2+sx1*5.d0/6.d0*hr**2)*uP(ir,l)
     &                       -(1-sx0/12.d0*hr**2)*uP(ir-1,l))
     &                       /(1-sx2/12.d0*hr**2)
               sx0=sx1
               sx1=sx2
              if(abs(uP(ir+1,l)).gt.1d10) then
                 do 230 jr=1,ir+1
                    uP(jr,l)=uP(jr,l)/uP(ir+1,l)
 230             continue
              endif
 220       continue

           sum=0.d0
           do 240 ir=1,Nct
              sum=sum+uP(ir,l)**2*hr
 240       continue
           do 250 ir=1,Nct
              uP(ir,l)=uP(ir,l)/sqrt(sum)
 250       continue
           do 260 ir=Nct+1,Nmax
              uP(ir,l)=0.d0
 260       continue

 210     continue

      return
      end

C---------------------------------------------------------------------
      FUNCTION fc2(mP,c2)
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /pob/lP,PtoA,AtoP
      common /con/irRC,irPC
      common /aer/rhoA,VA,speA,FA,GA
      common /ppr/rhoP,VP,FP
      common /ctb/c0,c4,c6,c8,c10,c12,scr_int
c parameter
      integer NAdim,NPdim,Nx,Nr
      parameter(NAdim=26,NPdim=7)
      parameter(Nx=640,Nr=2000)
      real*8 hr
      parameter(hr=0.01d0)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer lP(NPdim),PtoA(NPdim),AtoP(NAdim)
      integer irRC(NPdim),irPC
      real*8 rhoA(0:2*Nx),VA(0:2*Nx),speA(NAdim)
      real*8 FA(0:Nx,NAdim),GA(0:Nx,NAdim)
      real*8 rhoP(0:Nr),VP(0:Nr,NPdim),FP(0:Nr,NPdim)
      real*8 c0,c4,c6,c8,c10,c12,scr_int
c local
      integer mP,l,irc,ir
      real*8 c2,e,rc,Frc,dFrc,Vrc,dVrc,ddVrc,p0,p1,p2,p3,p4,fc2
      real*8 a(6,6),b(6,1),sum,fact,r

      l=lP(mP)
      e=speA(PtoA(mP))
      irc=irRC(mP)
      rc=irc*hr

      if(FP(irc,mP).lt.0.d0) then
          do 10 ir=0,Nr
             FP(ir,mP)=-FP(ir,mP)
 10       continue
      endif

      Frc=FP(irc,mP)
      dFrc=(-FP(irc+2,mP)+8*FP(irc+1,mP)
     &      -8*FP(irc-1,mP)+FP(irc-2,mP))/(12*hr)
      Vrc=VP(irc,mP)
      dVrc=(-VP(irc+2,mP)+8*VP(irc+1,mP)
     &      -8*VP(irc-1,mP)+VP(irc-2,mP))/(12*hr)
      ddVrc=(-VP(irc+2,mP)+16*VP(irc+1,mP)-30*VP(irc,mP)
     &        +16*VP(irc-1,mP)-VP(irc-2,mP))/(12*hr**2)

      p0=log(Frc/rc**(l+1))
      p1=dFrc/Frc-(l+1)/rc
      p2=(Vrc-e)/H2M-p1**2-2*p1*(l+1)/rc
      p3=dVrc/H2M-2*p1*p2-2*p2*(l+1)/rc
     &   +2*p1*(l+1)/rc**2
      p4=ddVrc/H2M-2*p2**2
     &   -2*p1*p3-2*p3*(l+1)/rc+4*p2*(l+1)/rc**2-4*p1*(l+1)/rc**3

      a(1,1)=1.d0
      a(1,2)=rc**4
      a(1,3)=rc**6
      a(1,4)=rc**8
      a(1,5)=rc**10
      a(1,6)=rc**12
      a(2,1)=0.d0
      a(2,2)=4*rc**3
      a(2,3)=6*rc**5
      a(2,4)=8*rc**7
      a(2,5)=10*rc**9
      a(2,6)=12*rc**11
      a(3,1)=0.d0
      a(3,2)=12*rc**2
      a(3,3)=30*rc**4
      a(3,4)=56*rc**6
      a(3,5)=90*rc**8
      a(3,6)=132*rc**10
      a(4,1)=0.d0
      a(4,2)=24*rc
      a(4,3)=120*rc**3
      a(4,4)=336*rc**5
      a(4,5)=720*rc**7
      a(4,6)=1320*rc**9
      a(5,1)=0.d0
      a(5,2)=24.d0
      a(5,3)=360*rc**2
      a(5,4)=1680*rc**4
      a(5,5)=5040*rc**6
      a(5,6)=11880*rc**8
      a(6,1)=0.d0
      a(6,2)=2*l+5
      a(6,3)=0.d0
      a(6,4)=0.d0
      a(6,5)=0.d0
      a(6,6)=0.d0
      b(1,1)=p0-c2*rc**2
      b(2,1)=p1-2*c2*rc
      b(3,1)=p2-2*c2
      b(4,1)=p3
      b(5,1)=p4
      b(6,1)=-c2**2
      call gaussj(a,6,6,b,1,1)
      c0=b(1,1)
      c4=b(2,1)
      c6=b(3,1)
      c8=b(4,1)
      c10=b(5,1)
      c12=b(6,1)

      fc2=2*c0
      sum=0.d0
      do 30 ir=1,irc
         fact=1.d0
         if(ir.eq.irc) fact=0.5d0
         r=ir*hr
         sum=sum+fact*r**(2*l+2)*exp(2*c2*r**2+2*c4*r**4
     &          +2*c6*r**6+2*c8*r**8+2*c10*r**10+2*c12*r**12)*hr
 30   continue
      fc2=fc2+log(sum)-log(scr_int)

      return
      end

c---------------------------------------------------------------------
      SUBROUTINE gaussj(a,n,np,b,m,mp)
      implicit none
      integer m,mp,n,np,NMAX
      real*8 a(np,np),b(np,mp)
      parameter(NMAX=50)
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),
     &        ipiv(NMAX)
      real*8 big,dum,pivinv
      do j=1,n
         ipiv(j)=0
      enddo
      do i=1,n
         big=0.d0
         do j=1,n
            if(ipiv(j).ne.1) then
               do k=1,n
                  if(ipiv(k).eq.0) then
                     if(abs(a(j,k)).ge.big) then
                        big=abs(a(j,k))
                        irow=j
                        icol=k
                     endif
                  else if(ipiv(k).gt.1) then
                     write(*,*) 'singular matrix in gaussj'
                     stop
                  endif
               enddo
            endif
         enddo
         ipiv(icol)=ipiv(icol)+1

         if(irow.ne.icol) then
            do l=1,n
               dum=a(irow,l)
               a(irow,l)=a(icol,l)
               a(icol,l)=dum
            enddo
            do l=1,m
               dum=b(irow,l)
               b(irow,l)=b(icol,l)
               b(icol,l)=dum
            enddo
         endif
         indxr(i)=irow
         indxc(i)=icol
         if(a(icol,icol).eq.0.d0) then
            write(*,*) 'singular matrix in gaussj'
            stop
         endif
         pivinv=1.d0/a(icol,icol)
         a(icol,icol)=1.d0
         do l=1,n
            a(icol,l)=a(icol,l)*pivinv
         enddo
         do l=1,m
            b(icol,l)=b(icol,l)*pivinv
         enddo
         do ll=1,n
            if(ll.ne.icol) then
               dum=a(ll,icol)
               a(ll,icol)=0.d0
               do l=1,n
                  a(ll,l)=a(ll,l)-a(icol,l)*dum
               enddo
               do l=1,m
                  b(ll,l)=b(ll,l)-b(icol,l)*dum
               enddo
            endif
         enddo
      enddo
      do l=n,1,-1
         if(indxr(l).ne.indxc(l)) then
            do k=1,n
               dum=a(k,indxr(l))
               a(k,indxr(l))=a(k,indxc(l))
               a(k,indxc(l))=dum
            enddo
         endif
      enddo
      return
      end

C---------------------------------------------------------------------
c exchange correlation potential with relativistic correction
      FUNCTION V_xc(rho)
      implicit none
c parameter
      real*8 Pi,E2,a_B,Ry
      real*8 gammaU,beta1U,beta2U,AU,BU,CU,DU
      parameter(Pi=3.141592653589793d0)
      parameter(E2=14.4d0,a_B=0.529177d0,Ry=13.6058d0)
      parameter(gammaU=-0.1423d0,beta1U=1.0529d0,beta2U=0.3334d0)
      parameter(AU=0.0311d0,BU=-0.048d0,CU=0.002d0,DU=-0.0116d0)
c local
      real*8 rs,rssq,rsln
      real*8 beta,eta,xi,psix
c in-out
      real*8 rho,V_xc

      if(rho.lt.1.d-80) rho=1.d-80
      rs=(3.d0/(4*Pi*rho))**(1.d0/3.d0)
      rs=rs/a_B

      beta=1.d0/(71.4036d0*rs)
      eta=sqrt(1+beta**2)
      xi=beta+eta
      psix=-0.5d0+1.5d0*log(xi)/beta/eta

      V_xc=-4.d0/3.d0*0.4582d0/rs*psix
      if(rs.gt.1.d0) then
         rssq=sqrt(rs)
         V_xc=V_xc+gammaU*(1.d0+7.d0/6.d0*beta1U*rssq
     &                   +4.d0/3.d0*beta2U*rs)
     &          /(1.d0+beta1U*rssq+beta2U*rs)**2
      else
         rsln=log(rs)
         V_xc=V_xc+AU*rsln+(BU-AU/3.d0)
     &          +2.d0/3.d0*CU*rs*rsln+(2.d0*DU-CU)/3.d0*rs
      endif
      V_xc=2.d0*Ry*V_xc
      return
      end

C---------------------------------------------------------------------
c exchange correlation energy density with relativistic correction
      FUNCTION E_xc(rho)
      implicit none
c parameter
      real*8 Pi,E2,a_B,Ry
      real*8 gammaU,beta1U,beta2U,AU,BU,CU,DU
      parameter(Pi=3.141592653589793d0)
      parameter(E2=14.4d0,a_B=0.529177d0,Ry=13.6058d0)
      parameter(gammaU=-0.1423d0,beta1U=1.0529d0,beta2U=0.3334d0)
      parameter(AU=0.0311d0,BU=-0.048d0,CU=0.002d0,DU=-0.0116d0)
c local
      real*8 rs,rssq,rsln
      real*8 beta,eta,xi,phix
c in-out
      real*8 rho,E_xc

      if(rho.lt.1.d-80) rho=1.d-80
      rs=(3.d0/(4*Pi*rho))**(1.d0/3.d0)
      rs=rs/a_B

      beta=1.d0/(71.4036d0*rs)
      eta=sqrt(1+beta**2)
      xi=beta+eta
      phix=1.d0-1.5d0*((beta*eta-log(xi))/beta**2)**2

      E_xc=-.4582d0/rs*phix
      if(rs.gt.1.d0) then
         rssq=sqrt(rs)
         E_xc=E_xc+gammaU/(1.d0+beta1U*rssq+beta2U*rs)
      else
         rsln=log(rs)
         E_xc=E_xc+AU*rsln+BU+CU*rs*rsln+DU*rs
      endif
      E_xc=2.d0*Ry*E_xc
      return
      end

C---------------------------------------------------------------------
c spline interpolation
c NOT for general purpose    (modified from Numerical Recipes)
      SUBROUTINE spline(x,y,y2)
      implicit none
c parameter
      integer Nx
      parameter(Nx=640)
c local
      integer i, k
      real*8 x(0:Nx),y(0:Nx),y2(0:Nx)
      real*8 p,qn,sig,un,u(0:Nx)

      y2(0)=0.d0
      u(0)=0.d0

      do 10 i=1,Nx-1
         sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
         p=sig*y2(i-1)+2
         y2(i)=(sig-1.d0)/p
         u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     &        /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
10    continue

      qn=0.d0
      un=0.d0
      y2(Nx)=(un-qn*u(Nx-1))/(qn*y2(Nx-1)+1.d0)
      do 20 k=Nx-1,0,-1
         y2(k)=y2(k)*y2(k+1)+u(k)
 20   continue
      return
      end

C---------------------------------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,klo,khi,x,y)
      implicit none
c parameter
      integer Nx
      parameter(Nx=640)
c local
      integer k,klo,khi
      real*8 x,y,xa(0:Nx),y2a(0:Nx),ya(0:Nx),a,b,h

      if(khi-klo.eq.1.and.xa(khi).gt.x.and.xa(klo).lt.x) go to 2

      klo=0
      khi=Nx
 1    if(khi-klo.gt.1) then
         k=(khi+klo)/2
         if(xa(k).gt.x) then
             khi=k
         else
             klo=k
         endif
      go to 1
      endif
      h=xa(khi)-xa(klo)
      if(h.eq.0.d0) then
          write(*,*)   'bad xa input in splint'
          stop
      endif

 2    a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     &   ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
      return
      end

C---------------------------------------------------------------------
      SUBROUTINE init
      implicit none

      common /prm/HC,MC2,H2M,E2,ALP,Pi
      common /aob/nodeA,lA,jA,kA
c parameter
      integer NAdim
      parameter(NAdim=26)
c common
      real*8 HC,MC2,H2M,E2,ALP,Pi
      integer nodeA(NAdim),lA(NAdim),jA(NAdim),kA(NAdim)
c local
      integer mA

      HC=1.97329d3
      MC2=0.511003d6
      H2M=HC*HC/(2*MC2)
      E2=14.39976d0
      ALP=E2/HC
      Pi=3.141592653589793d0

c 1s1/2 (2)
      nodeA(1)=0
      lA(1)=0
      jA(1)=2
c 2s1/2
      nodeA(2)=1
      lA(2)=0
      jA(2)=2
c 2p1/2
      nodeA(3)=0
      lA(3)=1
      jA(3)=2
c 2p3/2 (10)
      nodeA(4)=0
      lA(4)=1
      jA(4)=4
c 3s1/2
      nodeA(5)=2
      lA(5)=0
      jA(5)=2
c 3p1/2
      nodeA(6)=1
      lA(6)=1
      jA(6)=2
c 3p3/2 (18)
      nodeA(7)=1
      lA(7)=1
      jA(7)=4
c 3d3/2
      nodeA(8)=0
      lA(8)=2
      jA(8)=4
c 3d5/2 (28)
      nodeA(9)=0
      lA(9)=2
      jA(9)=6
c 4s1/2
      nodeA(10)=3
      lA(10)=0
      jA(10)=2
c 4p1/2
      nodeA(11)=2
      lA(11)=1
      jA(11)=2
c 4p3/2
      nodeA(12)=2
      lA(12)=1
      jA(12)=4
c 4d3/2
      nodeA(13)=1
      lA(13)=2
      jA(13)=4
C 4d5/2
      nodeA(14)=1
      lA(14)=2
      jA(14)=6
c 4f5/2
      nodeA(15)=0
      lA(15)=3
      jA(15)=6
c 4f7/2 (60)
      nodeA(16)=0
      lA(16)=3
      jA(16)=8
c 5s1/2
      nodeA(17)=4
      lA(17)=0
      jA(17)=2
c 5p1/2
      nodeA(18)=3
      lA(18)=1
      jA(18)=2
c 5p3/2
      nodeA(19)=3
      lA(19)=1
      jA(19)=4
c 5d3/2
      nodeA(20)=2
      lA(20)=2
      jA(20)=4
c 5d5/2 (78)
      nodeA(21)=2
      lA(21)=2
      jA(21)=6
c 6s1/2
      nodeA(22)=5
      lA(22)=0
      jA(22)=2
c 6p1/2
      nodeA(23)=4
      lA(23)=1
      jA(23)=2
c 6p3/2
      nodeA(24)=4
      lA(24)=1
      jA(24)=4
c 5f5/2
      nodeA(25)=1
      lA(25)=3
      jA(25)=6
c 5f7/2
      nodeA(26)=1
      lA(26)=3
      jA(26)=8

c kappa table
      do 10 mA=1,NAdim
         if(jA(mA).eq.2*lA(mA)+2) then
            kA(mA)=-jA(mA)/2
         else
            kA(mA)=jA(mA)/2
         endif
 10   continue

      return
      end
