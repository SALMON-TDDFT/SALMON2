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
module eigen_subdiag_sub
  implicit none

contains

  subroutine eigen_real8_pdsyev(pinfo,info,h,e,v)
    ! scalapack is used: --enable-scalapack must be put in configure
    ! This is for test version for 4 nodes (assuming space grid is divided to 4 block)
    ! still incorrect calculation.....
    use structures, only: s_process_info, s_orbital_parallel
    use communication, only: comm_summation, comm_is_root
    use parallelization, only: nproc_id_global
    implicit none
    type(s_process_info) :: pinfo
    type(s_orbital_parallel),intent(in) :: info
    real(8), intent(in)  :: h(:,:)
    real(8), intent(out) :: e(:)
    real(8), intent(out) :: v(:,:)
    real(8), allocatable :: h_cp(:,:)
    real(8), allocatable :: h_div(:,:), v_div(:,:), v_tmp(:,:), work(:)
    integer :: n, nb, lwork, i,j,ii,jj, iia,jja
    integer :: nprow, npcol
    integer :: context, iam, ierr, mycol, myrow, nprocs
    integer :: lld_r,lld_c,numroc,indxg2p
    integer :: qrmem,mpc0,nqc0,iroffc,icoffc,icrow,iccol,sizemqrleft
    integer :: ic,jc,mb_c,nb_c,rsrc_c,csrc_c
    integer :: nn,np,nq,nrc,ldc
    integer,allocatable :: desca(:), descz(:)

    n  = ubound(h,1)
    nb = 1  !blocking factor -- probably parameter relating to efficiency

    !XXXX change here 
    if(pinfo%npdomain_orbital(1) > 1) then
       nprow = pinfo%npdomain_orbital(1)
       npcol = pinfo%npdomain_orbital(2)
    else
       nprow = pinfo%npdomain_orbital(2)
       npcol = pinfo%npdomain_orbital(3)
    endif

   !write(*,*) "nprow,npcol", nprow, npcol ; flush(6)

    allocate( h_cp(n,n) )
    h_cp = h
    allocate( v_tmp(n,n) )
    v_tmp= 0d0 !unnecessary?
    v    = 0d0 !unnecessary?


    if(.not.pinfo%flag_blacs_gridinit) then

       CALL BLACS_PINFO( pinfo%iam, pinfo%nprocs )
       write(*,*) "iam,nprocs=", pinfo%iam, pinfo%nprocs ; flush(6)

       IF( pinfo%nprocs .lt. 1 ) CALL BLACS_SETUP( pinfo%iam, NPROW*NPCOL )
       write(*,*) "iam=", pinfo%iam ; flush(6)

       CALL BLACS_GET( -1,0,pinfo%context )
       CALL BLACS_GRIDINIT( pinfo%context, 'R', NPROW, NPCOL )
       CALL BLACS_GRIDINFO( pinfo%context, NPROW,NPCOL, pinfo%myrow,pinfo%mycol )
       pinfo%flag_blacs_gridinit = .true.

       if(comm_is_root(nproc_id_global)) then
          write(*,*) "  scalapack:pdsyev is used"
          write(*,*) "    nprow, npcol = ", nprow, npcol
       endif
    endif

    context = pinfo%context
    iam     = pinfo%iam
    nprocs  = pinfo%nprocs
    myrow   = pinfo%myrow
    mycol   = pinfo%mycol

    !write(*,*) "    context = ", context
    !write(*,*) "myrow,mycol", myrow,mycol ; flush(6)

    lld_r = max(1,numroc(n, nb, myrow, 0, nprow))
    lld_c = max(1,numroc(n, nb, mycol, 0, npcol))
    allocate(desca(lld_r), descz(lld_r))
    allocate( h_div(lld_r,lld_c), v_div(lld_r,lld_c) )

    !write(*,*) "lld_r, lld_c=", lld_r, lld_c ; flush(6)

    !(set lwork) : by manual (correct?)
    qrmem = 2*n-2
    ic = 1
    jc = 1
    mb_c = nb
    nb_c = nb
    rsrc_c = 0
    csrc_c = 0
    iroffc = mod(ic-1, mb_c)
    icoffc = mod(jc-1, nb_c)
    icrow = indxg2p(ic, mb_c, MYROW, rsrc_c, NPROW)
    iccol = indxg2p(jc, nb_c, MYCOL, csrc_c, NPCOL)
    mpc0  = numroc(n+iroffc, mb_c, MYROW, icrow, NPROW)
    nqc0  = numroc(n+icoffc, nb_c, MYCOL, iccol, NPCOL)
    sizemqrleft = max( (nb*(nb-1))/2, (nqc0 + mpc0)*nb ) + nb*nb
    nn = max(n, nb, 2)
    np = numroc(nn, nb, 0, 0, NPROW)
    nq = numroc(max(n, nb, 2), nb, 0, 0, NPCOL)
    nrc = numroc(n, nb, myrow, 0, NPROCS)
    ldc = max(1, nrc)
    lwork = 5*n + n*ldc + max(sizemqrleft, qrmem) + 1
   !lwork = lwork * 5  !xxx larger in cases
    allocate(work(lwork))


    IF( MYROW.EQ.-1 ) GO TO 20

    CALL DESCINIT( DESCA, n, n, nb, nb, 0, 0, context, lld_r, ierr )
    if(ierr < 0) write(*,*) "error1 in pdsyev"
    CALL DESCINIT( DESCZ, n, n, nb, nb, 0, 0, context, lld_r, ierr )
    if(ierr < 0) write(*,*) "error2 in pdsyev"

    !cut out and put into small divided block for each process
    do i=1,lld_r
    do j=1,lld_c
       if(myrow+1==nprow) then
          ii = n-lld_r+i
       else
          ii = myrow*lld_r+i
       endif
       if(mycol+1==npcol) then
          jj = n-lld_c+j
       else
          jj = mycol*lld_c+j
       endif
       h_div(i,j) = h(ii,jj)
    enddo
    enddo

    if(myrow+1==nprow) then
       iia = n-lld_r+1
    else
       iia = myrow*lld_r+1
    endif
    if(mycol+1==npcol) then
       jja = n-lld_c+1
    else
       jja = mycol*lld_c+1
    endif


    CALL PDSYEV( 'V','U', n, h_div,1,1,DESCA,e,v_div,1,1,DESCZ,work,lwork,ierr )
!xxx    CALL PDSYEV( 'V','U', n, h_div,iia,jja,DESCA,e,v_div,iia,jja,DESCZ,work,lwork,ierr )
!xxx    CALL PDSYEV( 'V','U', n, h_cp,iia,jja,DESCA,e,v_tmp,iia,jja,DESCZ,work,lwork,ierr )
    if(ierr < 0) write(*,*) "error3 in pdsyev"

   !write(*,*) "hoge",iam, DESCZ

    !get together from each process to make n x n size
    do i=1,lld_r
    do j=1,lld_c
       if(myrow+1==nprow) then
          ii = n-lld_r+i
       else
          ii = myrow*lld_r+i
       endif
       if(mycol+1==npcol) then
          jj = n-lld_c+j
       else
          jj = mycol*lld_c+j
       endif
       v_tmp(ii,jj) = v_div(i,j) 
    enddo
    enddo
    call comm_summation(v_tmp, v, n*n, info%icomm_rko) !!!!xxxx


20  continue
!    CALL BLACS_GRIDEXIT( context )
!    CALL BLACS_EXIT( 1 )

    deallocate(work,h_div,v_div,v_tmp)

    return
  end subroutine eigen_real8_pdsyev

!----------------------------------------------------------------
subroutine eigen_subdiag(Rmat,evec,iter,ier2,pinfo)
  use parallelization, only: nproc_size_global
  use structures, only: s_process_info
  implicit none
  
  integer :: iter,ier2
  real(8) :: Rmat(iter,iter)
  real(8) :: evec(iter,iter)
  type(s_process_info),intent(in) :: pinfo
  
  character(1) :: JOBZ,UPLO
  integer :: N
  real(8) :: A(iter,iter)
  integer :: LDA
  real(8) :: W(iter)
  real(8) :: WORK(3*iter-1)
  integer :: LWORK
  
  ier2=0
  
  if(nproc_size_global==1)then
    JOBZ='V'
    UPLO='L'
    N=iter
    A=Rmat
    LDA=iter
    LWORK=3*iter-1
    call DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,ier2)
    evec=A
  else
    call SAMPLE_PDSYEV_CALL(Rmat,evec,iter,pinfo)
  end if

end subroutine eigen_subdiag

subroutine eigen_subdiag_periodic(Rmat,evec,iter,ier2)
  implicit none
  character :: JOBZ, UPLO
  integer :: LWORK
  integer :: iter,ier2
  real(8),allocatable :: RWORK(:)
  real(8) :: W(iter)
  complex(8) :: Rmat(iter,iter)
  complex(8),allocatable :: WORK(:)
  complex(8) :: evec(iter,iter)

  ier2=0

  JOBZ='V'
  UPLO='U'

  LWORK=2*iter-1
  allocate(WORK(LWORK))
  allocate(RWORK(3*iter-2))

  call ZHEEV(JOBZ,UPLO,iter,Rmat,iter,W,WORK,LWORK,RWORK,ier2)

  evec(:,:)=Rmat(:,:)

  deallocate(WORK,RWORK)

end subroutine eigen_subdiag_periodic
 
!
      subroutine SAMPLE_PDSYEV_CALL(Rmat,evec,iter,pinfo)
!
!
!  -- ScaLAPACK routine (version 1.2) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 10, 1996
!
!     This routine contains a sample call to PDSYEV.
!     When compiled and run, it produces output which can be
!     pasted directly into matlab.
!
!     .. Parameters ..
      use structures, only: s_process_info
      use parallelization, only: nproc_size_global
      integer :: iter
      real(8) :: Rmat(iter,iter)
      real(8) :: evec(iter,iter)
      type(s_process_info),intent(in) :: pinfo
      INTEGER            LWORK, MAXN
      INTEGER            LIWORK
      INTEGER,allocatable :: IWORK(:)
      DOUBLE PRECISION   ZERO
!      PARAMETER          ( LWORK = 264, MAXN = 100, ZERO = 0.0D+0 )
      PARAMETER          ( MAXN = 20000, ZERO = 0.0D+0 )
      INTEGER            LDA
      DOUBLE PRECISION   MONE
      INTEGER            MAXPROCS
!      PARAMETER          ( LDA = MAXN, MONE = -1.0D+0, MAXPROCS = 512 )
      PARAMETER          ( MONE = -1.0D+0, MAXPROCS = 2048 )
      integer :: ii,jj
      character(1) :: SCOPE,TOP
!     ..
!     .. Local Scalars ..
      INTEGER            INFO, N, NB 
!      INTEGER            CONTEXT, I, IAM, INFO, MYCOL, MYROW, N, NB, NPCOL, NPROCS, NPROW
      integer :: CONTEXT, IAM, MYCOL, MYROW, NPCOL, NPROCS2, NPROW
!     ..
!     .. Local Arrays ..
      INTEGER            DESCA( 50 ), DESCZ( 50 )
!      DOUBLE PRECISION   A( LDA, LDA ), W( MAXN ), WORK( LWORK ), Z( LDA, LDA )
      DOUBLE PRECISION   A( iter, iter ), W( iter ), Z( iter, iter )
      real(8), allocatable :: WORK( : )
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_EXIT, BLACS_GET, BLACS_GRIDEXIT,   &
                         BLACS_GRIDINFO, BLACS_GRIDINIT, BLACS_PINFO,&
                         BLACS_SETUP, DESCINIT, PDLAMODHILB, PDLAPRNT,&
                         PDSYEV
      INTEGER :: NP,NQ,TRILWMIN

!     ..
!     .. Executable Statements ..
!
!
!     Set up the problem
!
!      N = 4
      N = iter
      NB = 1
!      NPROW = 2
!      NPCOL = 2
      if(pinfo%npdomain_orbital(1)>1)then
        NPROW = pinfo%npdomain_orbital(1)
      else if(pinfo%npdomain_orbital(2)>1)then
        NPROW = pinfo%npdomain_orbital(2)
      else
        NPROW = pinfo%npdomain_orbital(3)
      end if
      NPCOL = nproc_size_global/NPROW
      LDA = iter
      
      NP = max(N,NPROW)
      NQ = max(N,NPCOL)
      TRILWMIN = 3*N + max( NP+1, 3 )
      LWORK=max(1+6*N+2*NP*NQ, TRILWMIN) + 2*N
      allocate(WORK(LWORK))
!
!
!      if(iblacsinit==0)then
!     Initialize the BLACS
!
        CALL BLACS_PINFO( IAM, NPROCS2 )
        IF( ( NPROCS2.LT.1 ) ) THEN
           CALL BLACS_SETUP( IAM, NPROW*NPCOL )
        END IF
!
!
!     Initialize a single BLACS context
!
        CALL BLACS_GET( -1, 0, CONTEXT )
        CALL BLACS_GRIDINIT( CONTEXT, 'R', NPROW, NPCOL )
        CALL BLACS_GRIDINFO( CONTEXT, NPROW, NPCOL, MYROW, MYCOL )
!
!      end if
!     Bail out if this process is not a part of this context.
!
      IF( MYROW.EQ.-1 ) GO TO 20
!
!
!     These are basic array descriptors
!
!      if(iblacsinit==0)then
        CALL DESCINIT( DESCA, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
        CALL DESCINIT( DESCZ, N, N, NB, NB, 0, 0, CONTEXT, LDA, INFO )
!        CALL BLACS_GRIDINFO( DESCA( 2 ), NPROW, NPCOL, MYROW, MYCOL )
!      end if
!
!     Build a matrix that you can create with
!     a one line matlab command:  hilb(n) + diag([1:-1/n:1/n])
!
!      CALL PDLAMODHILB( N, A, 1, 1, DESCA, INFO, Rmat, iter )
!
!     Ask PDSYEV to compute the entire eigendecomposition
!
      do jj= 1, N
         do ii= 1, N
           CALL PDELSET( A, ii, jj, DESCA, Rmat(ii,jj) ) 
         end do  
      end do

      LIWORK=7*N+8*NPCOL+2
      allocate(IWORK(LIWORK))

      CALL PDSYEVD( 'V', 'L', N, A, 1, 1, DESCA, W, Z, 1, 1, DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )

      deallocate(IWORK)
!
!     Print out the eigenvectors
!
!      CALL PDLAPRNT( N, N, Z, 1, 1, DESCZ, 0, 0, 'Z', 6, WORK )


!      evec2=0.d0

!      if(mod(N,NPROW)==0)then
!        maxii=N/NPROW
!      else
!        if(MYROW<mod(N,NPROW))then
!          maxii=N/NPROW+1
!        else
!          maxii=N/NPROW
!        end if
!      end if

!      if(mod(N,NPCOL)==0)then
!        maxjj=N/NPCOL
!      else
!        if(MYCOL<mod(N,NPCOL))then
!          maxjj=N/NPCOL+1
!        else
!          maxjj=N/NPCOL
!        end if
!      end if

!      do jj=0,maxjj-1
!      do ii=0,maxii-1
!        evec2(ii*NPROW+(MYROW+1),jj*NPCOL+(MYCOL+1))=Z(ii+1,jj+1)
!      end do
 !     end do

!      call comm_summation(evec2,evec,N*N,nproc_group_global)

      SCOPE='A'
      TOP=' '
      do jj=1,N
        do ii=1,N
          call PDELGET( SCOPE, TOP, evec(ii,jj), Z, ii, jj, DESCZ )
        end do
      end do

!      CALL BLACS_GRIDEXIT( CONTEXT )
!
   20 CONTINUE
!      iblacsinit=1
!
!      CALL BLACS_EXIT( 0 )
!      CALL BLACS_EXIT( 1 )
!
!
!     Uncomment this line on SUN systems to avoid the useless print out
!
!      CALL IEEE_FLAGS( 'clear', 'exception', 'underflow', '')
!
!
 9999 FORMAT( 'W=diag([', 4D16.12, ']);' )
!
!      STOP

      deallocate(WORK)

      END SUBROUTINE SAMPLE_PDSYEV_CALL
!
      SUBROUTINE PDLAMODHILB( N, A, IA, JA, DESCA, INFO, Rmat,iter )
!
!  -- ScaLAPACK routine (version 1.2) --
!     University of Tennessee, Knoxville, Oak Ridge National Laboratory,
!     and University of California, Berkeley.
!     May 10, 1996
!
!
!
!
!     .. Parameters ..
      integer :: iter
      real(8) :: Rmat(iter,iter)
      INTEGER            BLOCK_CYCLIC_2D, DLEN_, DT_, CTXT_, M_, N_,  &
                         MB_, NB_, RSRC_, CSRC_, LLD_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DT_ = 1,  &
                         CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6,  &
                         RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Scalar Arguments ..
      INTEGER            IA, INFO, JA, N
!     ..
!     .. Array Arguments ..
      INTEGER            DESCA( * )
      DOUBLE PRECISION   A( * )
!     ..
!     .. Local Scalars ..
      INTEGER            I, J, MYCOL, MYROW, NPCOL, NPROW
!     ..
!     .. External Subroutines ..
      EXTERNAL           BLACS_GRIDINFO, PDELSET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          DBLE
!     ..
!     .. Executable Statements ..
!
!
!     The matlab code for a real matrix is:
!         hilb(n) + diag( [ 1:-1/n:1/n ] )
!     The matlab code for a complex matrix is:
!         hilb(N) + toeplitz( [ 1 (1:(N-1))*i ] )
!
!       This is just to keep ftnchek happy
      IF( BLOCK_CYCLIC_2D*CSRC_*CTXT_*DLEN_*DT_*LLD_*MB_*M_*NB_*N_* RSRC_.LT.0 )RETURN
!
      INFO = 0
!
      CALL BLACS_GRIDINFO( DESCA( CTXT_ ), NPROW, NPCOL, MYROW, MYCOL )
!
!
      IF( IA.NE.1 ) THEN
         INFO = -3
      ELSE IF( JA.NE.1 ) THEN
         INFO = -4
      END IF
!
      DO 20 J = 1, N
         DO 10 I = 1, N
           CALL PDELSET( A, I, J, DESCA, Rmat(I,J) ) 
!            IF( I.EQ.J ) THEN
!               CALL PDELSET( A, I, J, DESCA,  &
!                             ( DBLE( N-I+1 ) ) / DBLE( N )+ONE /  &
!                             ( DBLE( I+J )-ONE ) )   &
!            ELSE
!               CALL PDELSET( A, I, J, DESCA, ONE / ( DBLE( I+J )-ONE ) )
!            END IF
   10    CONTINUE
   20 CONTINUE
!
!
      RETURN
!
!     End of PDLAMODHLIB
!
      END SUBROUTINE

end module eigen_subdiag_sub
