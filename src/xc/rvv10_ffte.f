C
C  Copyright 2018-2020 SALMON developers
C
C  Licensed under the Apache License, Version 2.0 (the "License");
C  you may not use this file except in compliance with the License.
C  You may obtain a copy of the License at
C
C      http://www.apache.org/licenses/LICENSE-2.0
C
C  Unless required by applicable law or agreed to in writing, software
C  distributed under the License is distributed on an "AS IS" BASIS,
C  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C  See the License for the specific language governing permissions and
C  limitations under the License.
C
C======================================================================
C======================================================================
C
C     FFTE: A FAST FOURIER TRANSFORM PACKAGE
C
C     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
C                BY
C         DAISUKE TAKAHASHI
C         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
C         UNIVERSITY OF TSUKUBA
C         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
C         E-MAIL: daisuke@cs.tsukuba.ac.jp
C
C
C     Private-table wrapper for rVV10. Reuses native FFTE transpose/FFT
C     workers without overwriting the saved tables of Poisson's wrapper.
C     IOPT is -1 (forward) or +1 (normalized inverse).
      SUBROUTINE PZFFT3DV_RVV10(A,B,NX,NY,NZ,NPUY,NPUZ,IOPT,
     1                        ICOMMY,ICOMMZ)
      IMPLICIT REAL*8 (A-H,O-Z)
! Parameters written in param.h !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PARAMETER (MAXNPU=65536)
      PARAMETER (NDA2=65536)
      PARAMETER (NDA3=4096)
      PARAMETER (NBLK=16)
      PARAMETER (NP=8)
      PARAMETER (L2SIZE=2097152)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      COMPLEX*16 A(*),B(*)
      COMPLEX*16 C(NDA3)
      COMPLEX*16 WX(NDA3),WY(NDA3),WZ(NDA3)
      DIMENSION LNX(3),LNY(3),LNZ(3)
C
      NN=NX*(NY/NPUY)*(NZ/NPUZ)
C
      CALL FACTOR(NX,LNX)
      CALL FACTOR(NY,LNY)
      CALL FACTOR(NZ,LNZ)
C
        CALL SETTBL(WX,NX)
        CALL SETTBL(WY,NY)
        CALL SETTBL(WZ,NZ)
C
      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 10 I=1,NN
          A(I)=DCONJG(A(I))
   10   CONTINUE
      END IF
C
      IF (IOPT .EQ. -1 .OR. IOPT .EQ. -2) THEN
!$OMP PARALLEL PRIVATE(C)
        CALL PZFFT3DVF(A,A,A,A,A,A,A,B,B,B,B,B,B,B,C,WX,WY,WZ,NX,NY,NZ,
     1                 LNX,LNY,LNZ,NPUY,NPUZ,IOPT,ICOMMY,ICOMMZ)
!$OMP END PARALLEL
      ELSE
!$OMP PARALLEL PRIVATE(C)
        CALL PZFFT3DVB(A,A,A,A,A,A,A,B,B,B,B,B,B,B,C,WX,WY,WZ,NX,NY,NZ,
     1                 LNX,LNY,LNZ,NPUY,NPUZ,IOPT,ICOMMY,ICOMMZ)
!$OMP END PARALLEL
      END IF
C
      IF (IOPT .EQ. 1 .OR. IOPT .EQ. 2) THEN
        DN=1.0D0/(DBLE(NX)*DBLE(NY)*DBLE(NZ))
!$OMP PARALLEL DO
!DIR$ VECTOR ALIGNED
        DO 20 I=1,NN
          B(I)=DCONJG(B(I))*DN
   20   CONTINUE
      END IF
      RETURN
      END
