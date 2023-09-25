! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_discrete_ordinate_util

  use musica_constants, only : dk => musica_dk

  implicit none

  private
  public :: PSNDO, solver_constants_t

  integer, parameter :: NMUG = 10
  integer, parameter :: MAXSTR = 100
  integer, parameter :: MAXTRM = 100
  integer, parameter :: MAXSQT = 1000
  integer, parameter :: MXCLY = 151
  integer, parameter :: MXULV = 151
  integer, parameter :: MXCMU = 32
  integer, parameter :: MXUMU = 32
  integer, parameter :: MXPHI = 3
  real(dk), parameter    :: rZERO = 0.0_dk
  real(dk), parameter    :: ONEHALF = 0.5_dk
  real(dk), parameter    :: rONE  = 1.0_dk
  real(dk), parameter    :: rTWO  = 2.0_dk
  real(dk), parameter    :: rFOUR = 4.0_dk
  real(dk), parameter    :: rTEN  = 10.0_dk
  ! Discrete ordinate constants:
  ! For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;
  ! ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0
  integer, parameter :: NPHI  = 0
  integer, parameter :: IBCND = 0
  real(dk), parameter    :: ACCUR = 0.0001_dk
  real(dk), parameter    :: FBEAM = rONE
  real(dk), parameter    :: FISOT = rZERO
  real(dk), parameter    :: PHI0  = rZERO
  logical, parameter :: DELTAM = .true.
  logical, parameter :: LAMBER = .true.
  logical, parameter :: PLANK = .FALSE.
  logical, parameter :: USRANG = .FALSE.
  logical, parameter :: USRTAU = .FALSE.
  logical, parameter :: ONLYFL = .TRUE.
  logical, parameter :: PRNT(7) = .FALSE.

  ! This could probably be removed if something like this happens:
  ! https://github.com/j3-fortran/fortran_proposals/issues/214
  type :: solver_constants_t
    real(dk) :: DITHER
    real(dk) :: SQT(MAXSQT) ! Square roots of integers
    real(dk) :: GMU(NMUG)   ! Angle cosine quadrature points on (0,1)
    real(dk) :: GWT(NMUG)   ! Angle cosine quadrature weights on (0,1)
    real(dk) :: YLMG(0:MAXSTR,NMUG) ! Normalized associated Legendre polynomials
                                    !   and the NMUG quadrature angles
    real(dk) :: C(MAXTRM)   ! Integral from 0 to 1 of MU * P-sub-L(MU)
                            !   ( vanishes for L = 3, 5, 7, ... )
  end type solver_constants_t

  interface solver_constants_t
    module procedure :: constructor
  end interface solver_constants_t

contains

!=============================================================================*

  ! Constructs constant collection for Discrete Ordinate solver
  function constructor( ) result( this )

    type(solver_constants_t) :: this

    integer :: NS, K, JG, L
    real(dk) :: SGN, CL

    this%DITHER = rTEN*EPSILON( this%DITHER )
!                            ** Must dither more on Cray (14-digit prec)
    IF( this%DITHER < 1.E-10_dk ) this%DITHER = rTEN*this%DITHER

    DO NS = 1, MAXSQT
      this%SQT( NS ) = SQRT( REAL(NS,dk) )
    ENDDO

    CALL QGAUSN( NMUG, this%GMU, this%GWT )
    CALL LEPOLY( NMUG, 0, MAXSTR, this%GMU, this%YLMG, this )

!                       ** Convert Legendre polys. to negative GMU
    SGN  = -rONE

    DO K = 0, MAXSTR
      SGN  = - SGN
      DO JG = 1, NMUG
        this%YLMG(K,JG) = SGN*this%YLMG(K,JG)
      ENDDO
    ENDDO

    CL     = 0.125_dk
    this%C(2)   = rTEN*CL

    DO L = 4, MAXTRM, 2
      CL   = - CL*real((L - 3),dk) / real(L + 2,dk)
      this%C(L) = rTWO*real(2*L + 1,dk)*CL
    ENDDO

  end function constructor

!=============================================================================*

      SUBROUTINE PSNDO( dsdh, nid,                 &
                        NLYR, DTAUC, SSALB, PMOM,  &
                        ALBEDO, NSTR, UMU0,        &
                        RFLDIR, RFLDN, FLUP,       &
                        uavgso, uavgup, uavgdn,    &
                        solver_constants )

      use musica_constants, only : PI => kpi

      integer, intent(in) :: NLYR
      integer, intent(in) :: NSTR
      integer, intent(in) :: nid(0:)
      real(dk), intent(in)    :: ALBEDO
      real(dk), intent(in)    :: UMU0
      real(dk), intent(in)    :: dsdh(0:,:)
      real(dk), intent(in)    :: DTAUC(:)
      real(dk), intent(inout) :: SSALB(:)
      real(dk), intent(inout) :: PMOM(0:,:)
      real(dk), intent(out)   :: RFLDIR(:), RFLDN(:), FLUP(:)
      real(dk), intent(out)   :: uavgso(:)
      real(dk), intent(out)   :: uavgup(:)
      real(dk), intent(out)   :: uavgdn(:)
      type(solver_constants_t), intent(in) :: solver_constants

      real(dk), parameter     :: RPD  = PI / 180.0_dk
!     spherical geometry
      real(dk) tausla(0:NLYR), tauslau(0:NLYR), mu2(0:NLYR)
!     ..

!     local variables
      real(dk) ::   DFDT(NLYR+1),                                  &
                HL(0:NSTR), PHI(MXPHI),                        &
                TRNMED(NSTR), U0U(nstr,nlyr+1), UAVG(NLYR+1),  &
                UMU(NSTR), CWT(NSTR), UTAU(NLYR+1),            &
                UU(NSTR,NLYR,MXPHI)
!     ..
!     .. Local Scalars ..

      logical :: COMPAR, LYRCUT, PASS1
      integer :: IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL
      integer :: NCOS, NCUT, NN
      real(dk)    :: AZERR, AZTERM, BPLANK, COSPHI, DELM0, DUM, SGN, TPLANK
!     ..
!     .. Local Arrays ..

      integer :: NTAU, NUMU
      integer :: IPVT(NSTR*NLYR), LAYRU(NLYR+1)

      real(dk) ::   ANGCOS(1)
      real(dk) ::   AMB(NSTR/2,NSTR/2), APB(NSTR/2,NSTR/2),           &
                ARRAY(NSTR,NSTR),                                 &
                B(NSTR*NLYR), BDR(NSTR/2,0:NSTR/2), BEM(NSTR/2),  &
                CBAND(9*(NSTR/2)-2,NSTR*NLYR), CC(NSTR,NSTR),     &
                CMU(NSTR), DTAUCP(NLYR),                          &
                EMU(NSTR), EVAL(NSTR/2), EVECC(NSTR, NSTR),       &
                EXPBEA(0:NLYR), FLDIR(NLYR+1), FLDN(NLYR+1),      &
                FLYR(NLYR), GC(NSTR,NSTR,NLYR),                   &
                GL(0:NSTR,NLYR), GU(NSTR,NSTR,NLYR),              &
                HLPR(0:NSTR), KK(NSTR,NLYR), LL(NSTR,NLYR),       &
                OPRIM(NLYR), PHIRAD(MXPHI), PKAG(0:NLYR),         &
                PSI(NSTR), RMU(NSTR,0:NSTR/2), TAUC(0:NLYR),      &
                TAUCPR(0:NLYR), U0C(NSTR,NLYR+1), UTAUPR(NLYR+1), &
                UUM(NSTR,NLYR), WK(NSTR), XR0(NLYR),              &
                XR1(NLYR), YLM0(0:NSTR,1), YLMC(0:NSTR,NSTR),     &
                YLMU(0:NSTR,NSTR), Z(NSTR*NLYR), Z0(NSTR),        &
                Z0U(NSTR,NLYR), Z1(NSTR), Z1U(NSTR,NLYR),         &
                ZBEAM(NSTR,NLYR), ZJ(NSTR),                       &
                ZPLK0(NSTR,NLYR), ZPLK1(NSTR,NLYR), ZZ(NSTR,NLYR)

      real(dk) :: sindir(nlyr+1), sinup(nlyr+1), sindn(nlyr+1)

!gy added glsave and dgl to allow adjustable dimensioning in SOLVEC
      real(dk) GLSAVE(0:nstr), DGL(0:nstr)

      real(dk) :: AAD(NSTR/2,NSTR/2), EVECCD(NSTR/2,NSTR/2), &
                  EVALD(NSTR/2), WKD(NSTR)

      real(dk)    :: PLKAVG

!     ** Calculate cumulative optical depth
!     and dither single-scatter albedo
!     to improve numerical behavior of
!     eigen{value/vector} computation

      TAUC = rZERO

      DO LC = 1, NLYR
         IF( SSALB(LC) == rONE ) THEN
           SSALB(LC) = rONE - solver_constants%DITHER
         ENDIF
         TAUC(LC) = TAUC(LC - 1) + DTAUC(LC)
      ENDDO
!                                ** Check input dimensions and variables

      CALL CHEKIN( NLYR, DTAUC, SSALB, PMOM,          &
                   NTAU, UTAU, NSTR, NUMU, UMU, NPHI, &
                   PHI, UMU0, FISOT, ALBEDO,          &
                   HL, TAUC, solver_constants )

!                                 ** Zero internal and output arrays

      CALL  ZEROAL( EXPBEA(1:), FLYR, OPRIM, TAUCPR(1:), XR0, XR1, &
                    CMU, CWT, PSI, WK, Z0, Z1, ZJ,                 &
                    HLPR, YLM0(:,1), ARRAY, CC, EVECC,             &
                    GL, YLMC, YLMU,                                &
                    KK, LL, ZZ, ZPLK0, ZPLK1,                      &
                    GC, LAYRU, UTAUPR,                             &
                    GU, Z0U, Z1U, ZBEAM,                           &
                    EVAL, AMB, APB, IPVT, Z,                       &
                    RFLDIR, RFLDN, FLUP, UAVG, DFDT,               &
                    TRNMED, U0U, UU )

!                                 ** Perform various setup operations

      CALL SETDIS(                                              &
          dsdh, nid, tausla, tauslau, mu2,                      &
          CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR, &
          GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT,           &
          NCUT, NLYR, NTAU, NN, NSTR, PLANK,                    &
          NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU, &
          UTAUPR, UMU, UMU0, USRTAU, USRANG )

!                                 ** Print input information
      IF ( PRNT(1) ) THEN
           CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM, &
                        NTAU, UTAU, NSTR, NUMU, UMU,      &
                        PHI, UMU0, FISOT, ALBEDO, HL,     &
                        FLYR, LYRCUT,                     &
                        OPRIM, TAUC, TAUCPR, PRNT(7) )
      ENDIF

!                              ** Handle special case for getting albedo
!                                 and transmissivity of medium for many
!                                 beam angles at once
!                                   ** Calculate Planck functions

         BPLANK = rZERO
         TPLANK = rZERO
         PKAG   = rZERO

! ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
!           (EQ STWJ 5)

      KCONV = 0
!                                    ** Azimuth-independent case

      IF( FBEAM == rZERO .OR. (1. - UMU0) < 1.E-5 .OR. ONLYFL .OR. &
          (NUMU == 1 .AND. (1. - UMU(1)) < 1.E-5 ) ) THEN
        NAZ = 0
      ELSE
        NAZ  = NSTR - 1
      ENDIF

      AZIMUTH_LOOP: DO MAZIM = 0, NAZ

         IF( MAZIM == 0 ) THEN
           DELM0  = rONE
         ELSE
           DELM0  = rZERO
         ENDIF
!                             ** Get normalized associated Legendre
!                                polynomials for
!                                (a) incident beam angle cosine
!                                (b) computational and user polar angle
!                                    cosines
         IF( FBEAM > rZERO ) THEN
            NCOS   = 1
            ANGCOS = -UMU0
            CALL LEPOLY( NCOS, MAZIM, NSTR - 1, ANGCOS, YLM0, solver_constants )
         END IF


         IF( .NOT. ONLYFL .AND. USRANG ) THEN
            CALL LEPOLY( NUMU, MAZIM, NSTR-1, UMU, YLMU, solver_constants )
         ENDIF

         CALL LEPOLY( NN, MAZIM, NSTR-1, CMU, YLMC, solver_constants )

!                       ** Get normalized associated Legendre polys.
!                          with negative arguments from those with
!                          positive arguments; Dave/Armstrong Eq. (15)
         SGN  = -rONE

         DO L = MAZIM, NSTR - 1
            SGN  = - SGN
            DO IQ = NN + 1, NSTR
               YLMC(L,IQ) = SGN*YLMC(L,IQ - NN)
            ENDDO
         ENDDO
!     ** Specify users bottom reflectivity
!        and emissivity properties
         IF ( .NOT. LYRCUT ) THEN
           CALL  SURFAC(                           &
               ALBEDO, DELM0, FBEAM, HLPR, LAMBER, &
               MAZIM, NN, NUMU, NSTR, ONLYFL,      &
               UMU, USRANG, YLM0(:,1), YLMC, YLMU, &
               BDR, EMU, BEM, RMU, solver_constants )
         ENDIF


! ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO LC = 1, NCUT
            CALL SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL(0:,LC),    &
                 MAZIM, NN, NSTR, YLM0(:,1), YLMC, CC,            &
                 EVECC, EVAL, KK(:,LC ), GC(:,:,LC), AAD, EVECCD, &
                 EVALD, WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0,    &
                 ZJ, ZZ(:,LC), OPRIM(LC), LC, mu2(lc),            &
                 glsave, dgl, solver_constants)
         ENDDO


! ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


!                      ** Set coefficient matrix of equations combining
!                         boundary and layer interface conditions

         CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK, &
                      LAMBER, LYRCUT, NCOL, NCUT,                  &
                      NN, NSTR, TAUCPR, WK )

!                      ** Solve for constants of integration in homo-
!                         geneous solution (general boundary conditions)

         CALL SOLVE0(                                        &
             B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,   &
             FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,  &
             NCOL, NCUT, NN, NSTR, PI,                       &
             TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

!                                  ** Compute upward and downward fluxes

         IF ( MAZIM == 0 ) THEN
            CALL FLUXES( tausla, tauslau,                             &
                         CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,  &
                         NCUT, NN, NSTR, NTAU,                        &
                         PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, &
                         XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,      &
                         FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C,       &
                         uavgso, uavgup, uavgdn,                      &
                         sindir, sinup, sindn)
         ENDIF

         IF( ONLYFL ) THEN
            IF( MXUMU >= NSTR ) THEN
!     ** Save azimuthal-avg intensities
!        at quadrature angles
               DO LU = 1, NTAU
                  DO IQ = 1, NSTR
                     U0U(IQ,LU) = U0C(IQ,LU)
                  ENDDO
               ENDDO
            ENDIF
            EXIT AZIMUTH_LOOP
         ENDIF

         UUM = rZERO

         IF( MAZIM == 0 ) THEN
!     ** Save azimuthally averaged intensities
            DO LU = 1, NTAU
               DO IU = 1, NUMU
                  U0U(IU,LU) = UUM(IU,LU)
                  DO J = 1, NPHI
                     UU(IU,LU,J) = UUM(IU,LU)
                  ENDDO
               ENDDO
            ENDDO
!                              ** Print azimuthally averaged intensities
!                                 at user angles

            IF( PRNT( 4 ) ) THEN
               CALL PRAVIN( UMU, NUMU, UTAU, NTAU, U0U )
            ENDIF
            IF( NAZ > 0 ) THEN
               PHIRAD = rZERO
               DO J = 1, NPHI
                  PHIRAD(J) = RPD*(PHI(J) - PHI0)
               ENDDO
            END IF

         ELSE
!                                ** Increment intensity by current
!                                   azimuthal component (Fourier
!                                   cosine series);  Eq SD(2)
            AZERR  = rZERO

            DO J = 1, NPHI
               COSPHI = COS( MAZIM*PHIRAD(J) )
               DO LU = 1, NTAU
                  DO IU = 1, NUMU
                     AZTERM = UUM(IU,LU)*COSPHI
                     UU(IU,LU,J) = UU(IU,LU,J) + AZTERM
                     AZERR = MAX( AZERR, RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )
                  ENDDO
               ENDDO
            ENDDO

            IF( AZERR <= ACCUR ) KCONV  = KCONV + 1

            IF( KCONV >= 2 ) THEN
               EXIT AZIMUTH_LOOP
            ENDIF

         ENDIF

      ENDDO AZIMUTH_LOOP

! ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


!                                          ** Print intensities
      IF( PRNT( 5 ) .AND. .NOT. ONLYFL ) THEN
        CALL PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )
      ENDIF

      END SUBROUTINE PSNDO

      SUBROUTINE ASYMTX( AA, EVEC, EVAL, IER, WKD, AAD, EVECD, EVALD )

!    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

!       Solves eigenfunction problem for real asymmetric matrix
!       for which it is known a priori that the eigenvalues are real.

!       This is an adaptation of a subroutine EIGRF in the IMSL
!       library to use real instead of complex arithmetic, accounting
!       for the known fact that the eigenvalues and eigenvectors in
!       the discrete ordinate solution are real.  Other changes include
!       putting all the called subroutines in-line, deleting the
!       performance index calculation, updating many DO-loops
!       to Fortran77, and in calculating the machine precision
!       TOL instead of specifying it in a data statement.

!       EIGRF is based primarily on EISPACK routines.  The matrix is
!       first balanced using the Parlett-Reinsch algorithm.  Then
!       the Martin-Wilkinson algorithm is applied.

!       References:
!          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
!             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
!             Sources and Development of Mathematical Software,
!             Prentice-Hall, Englewood Cliffs, NJ
!         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
!             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
!         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
!             Clarendon Press, Oxford

!   I N P U T    V A R I A B L E S:

!       AA    :  input asymmetric matrix, destroyed after solved
!        M    :  order of  AA
!       IA    :  first dimension of  AA
!    IEVEC    :  first dimension of  EVEC

!   O U T P U T    V A R I A B L E S:

!       EVEC  :  (unnormalized) eigenvectors of  AA
!                   ( column J corresponds to EVAL(J) )

!       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

!       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
!                   in that case eigenvalues IER+1,IER+2,...,M  are
!                   correct but eigenvalues 1,...,IER are set to zero.

!   S C R A T C H   V A R I A B L E S:

!       WKD   :  work area ( dimension at least 2*M )
!       AAD   :  double precision stand-in for AA
!       EVECD :  double precision stand-in for EVEC
!       EVALD :  double precision stand-in for EVAL

!   Called by- SOLEIG
!   Calls- ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      integer, intent(out) :: IER
!     ..
!     .. Array Arguments ..

      real(dk), intent(in)  :: AA(:,:)
      real(dk), intent(out) :: EVAL(:), EVEC(:,:)
      real(dk), intent(inout) :: WKD(:)
      real(dk), intent(out)   :: AAD(:,:)
      real(dk), intent(out)   :: EVALD(:), EVECD(:,:)
!     ..
!     .. Local Scalars ..

      integer :: I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
      integer :: M, IA, IEVEC
      logical :: NOCONV, NOTLAS, CONV

      real(dk) :: COL, DISCRI, F, G, H,                        &
                 P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T, &
                 TOL, UU, VV, W, X, Y, Z

      real(dk), parameter :: C1 = 0.4375_dk
      real(dk), parameter :: C2 = ONEHALF
      real(dk), parameter :: C3 = 0.75_dk
      real(dk), parameter :: C4 = 0.95_dk
      real(dk), parameter :: C5 = 16._dk
      real(dk), parameter :: C6 = 256._dk

      IER  = 0
      TOL  = EPSILON( TOL )

      M     = size(EVAL)
      IA    = size(AA,dim=1)
      IEVEC = size(EVECD,dim=1)

      IF( M < 1 .OR. IA < M .OR. IEVEC < M ) &
          CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )

      ! ** Handle 1x1 and 2x2 special cases

      IF( M == 1 ) THEN
         EVAL( 1 )    = AA( 1, 1 )
         EVEC( 1, 1 ) = rONE
         RETURN
      ELSE IF( M == 2 ) THEN
         DISCRI = (AA(1,1) - AA(2,2))**2 + rFOUR*AA(1,2)*AA(2,1)

         IF( DISCRI.LT.rZERO )  &
             CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )

         SGN  = rONE

         IF( AA( 1,1 ) < AA( 2,2 ) ) SGN  = -rONE

         EVAL( 1 ) = ONEHALF*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
         EVAL( 2 ) = ONEHALF*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
         EVEC( 1, 1 ) = rONE
         EVEC( 2, 2 ) = rONE

         IF( AA( 1,1 ) == AA( 2,2 ) .AND.  &
             ( AA( 2,1 ) == rZERO .OR. AA( 1,2 ) == rZERO ) ) THEN

            RNORM  = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) + &
                     ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
            W  = TOL*RNORM
            EVEC( 2, 1 ) =   AA( 2, 1 ) / W
            EVEC( 1, 2 ) = - AA( 1, 2 ) / W
         ELSE
            EVEC( 2, 1 ) = AA( 2, 1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
            EVEC( 1, 2 ) = AA( 1, 2 ) / ( EVAL( 2 ) - AA( 1,1 ) )
         END IF

         RETURN
      END IF
!                               ** Put s.p. matrix into d.p. matrix
      DO J = 1, M
         DO K = 1, M
            AAD(J,K) = REAL(AA(J,K),kind=dk)
         ENDDO
      ENDDO

!                                ** Initialize output variables
      IER  = 0

      EVECD(1:M,1:M) = rZERO
      EVALD(1:M) = rZERO
      DO I = 1, M
         EVECD(I,I) = rONE
      ENDDO

      !  ** Balance the input matrix and reduce its norm by
      !     diagonal similarity transformation stored in WK;
      !     then search for rows isolating an eigenvalue
      !     and push them down
      RNORM  = rZERO
      L  = 1
      K  = M

OUTER_LOOP_A: &
      DO
         KKK  = K
         DO J = KKK, 1, -1
            ROW  = rZERO
            DO I = 1, K
               IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )
            ENDDO

            IF( ROW.EQ.rZERO ) THEN
               WKD( K ) = J
               IF( J.NE.K ) THEN
                  DO I = 1, K
                     REPL        = AAD( I, J )
                     AAD( I, J ) = AAD( I, K )
                     AAD( I, K ) = REPL
                  ENDDO

                  DO I = L, M
                     REPL        = AAD( J, I )
                     AAD( J, I ) = AAD( K, I )
                     AAD( K, I ) = REPL
                  ENDDO
               END IF
               K  = K - 1
               CYCLE OUTER_LOOP_A
            END IF
         ENDDO
         EXIT
      ENDDO OUTER_LOOP_A
      ! ** Search for columns isolating an
      !    eigenvalue and push them left
OUTER_LOOP_B: &
      DO
         LLL  = L
         DO J = LLL, K
            COL  = rZERO
            DO I = L, K
               IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )
            ENDDO

            IF( COL.EQ.rZERO ) THEN
               WKD( L ) = J
               IF( J.NE.L ) THEN
                  DO I = 1, K
                     REPL        = AAD( I, J )
                     AAD( I, J ) = AAD( I, L )
                     AAD( I, L ) = REPL
                  ENDDO

                  DO I = L, M
                     REPL        = AAD( J, I )
                     AAD( J, I ) = AAD( L, I )
                     AAD( L, I ) = REPL
                  ENDDO
               END IF

               L  = L + 1
               CYCLE OUTER_LOOP_B
            END IF
         END DO
         EXIT
      ENDDO OUTER_LOOP_B

      !  ** Balance the submatrix in rows L through K
      WKD( L:K ) = rONE

      DO
         CONV = .TRUE.
         DO I = L, K
            COL  = rZERO
            ROW  = rZERO

            DO J = L, K
               IF( J.NE.I ) THEN
                  COL  = COL + ABS( AAD( J,I ) )
                  ROW  = ROW + ABS( AAD( I,J ) )
               END IF
            ENDDO

            F  = rONE
            G  = ROW / C5
            H  = COL + ROW

            DO WHILE( COL.LT.G )
               F    = F*C5
               COL  = COL*C6
            ENDDO

            G  = ROW*C5

            DO WHILE( COL.GE.G )
               F    = F / C5
               COL  = COL / C6
            ENDDO
            !  ** Now balance
            IF( (COL + ROW )/F < C4*H ) THEN
               WKD( I ) = WKD(I)*F
               CONV = .FALSE.

               DO J = L, M
                  AAD( I, J ) = AAD( I, J ) / F
               ENDDO

               DO J = 1, K
                  AAD( J, I ) = AAD( J, I )*F
               ENDDO
            END IF
         ENDDO
         IF( CONV ) EXIT
      ENDDO
      !  ** Is A already in Hessenberg form?
      IF( K-1 >= L+1 ) THEN
         !  ** Transfer A to a Hessenberg form
         DO N = L + 1, K - 1
            H  = rZERO
            WKD( N + M ) = rZERO
            ! ** Scale column
            SCALE = SUM( ABS( AAD(N:K,N-1) ) )

            IF( SCALE /= rZERO ) THEN
               DO I = K, N, -1
                  WKD( I + M ) = AAD( I, N - 1 ) / SCALE
                  H  = H + WKD( I + M )**2
               ENDDO

               G    = - SIGN( SQRT( H ), WKD( N + M ) )
               H    = H - WKD( N + M )*G
               WKD( N + M ) = WKD( N + M ) - G
               ! ** Form (I-(U*UT)/H)*A
               DO J = N, M
                  F  = rZERO
                  DO I = K, N, -1
                     F  = F + WKD( I + M )*AAD( I, J )
                  ENDDO
                  DO I = N, K
                     AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
                  ENDDO
               ENDDO
               ! ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
               DO I = 1, K
                  F  = rZERO

                  DO J = K, N, -1
                     F  = F + WKD( J + M )*AAD( I, J )
                  ENDDO

                  DO J = N, K
                     AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
                  ENDDO
               ENDDO

               WKD(N + M) = SCALE*WKD(N + M)
               AAD(N, N - 1) = SCALE*G
            END IF
         ENDDO

         DO N = K - 2, L, -1
            N1   = N + 1
            N2   = N + 2
            F  = AAD( N + 1, N )

            IF( F /= rZERO ) THEN
               F  = F*WKD( N + 1 + M )
               DO I = N + 2, K
                  WKD( I + M ) = AAD( I, N )
               ENDDO

               IF( N + 1 <= K ) THEN
                  DO J = 1, M
                     G = DOT_PRODUCT( WKD(N+1+M:K+M),EVECD(N+1:K,J) )
                     G  = G / F
                     DO I = N + 1, K
                        EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
                     ENDDO
                  ENDDO
               END IF
            END IF
         ENDDO
      END IF

      N  = 1

      DO I = 1, M
         DO J = N, M
            RNORM  = RNORM + ABS( AAD( I,J ) )
         ENDDO

         N  = I

         IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )
      ENDDO

      N  = K
      T  = rZERO

      !  ** Search for next eigenvalues
OUTER_BLK: DO
      IF( N >= L ) THEN
         IN  = 0
         N1  = N - 1
         N2  = N - 2
      !  ** Look for single small sub-diagonal element
INNER_BLK: DO
         DO I = L, N
            LB  = N + L - I
            IF( LB == L ) EXIT
            S  = ABS(AAD(LB - 1,LB - 1)) + ABS(AAD(LB,LB))
            IF( S == rZERO ) S = RNORM
            IF(ABS( AAD(LB, LB-1)) <= TOL*S ) EXIT
         ENDDO

         X  = AAD( N, N )

         IF( LB == N ) THEN
         !  ** One eigenvalue found
            AAD( N, N ) = X + T
            EVALD( N ) = AAD( N, N )
            N  = N1
            CYCLE OUTER_BLK
         END IF

! next line has been included to avoid run time error caused by xlf

         IF( N1 <= 0 .OR. N <= 0 ) THEN
           WRITE(0,*) 'Subscript out of bounds in ASYMTX'
           STOP 9999
         ENDIF

         Y  = AAD( N1, N1 )
         W  = AAD( N, N1 )*AAD( N1, N )

         IF( LB == N1 ) THEN
!                                        ** Two eigenvalues found
            P  = ( Y - X )*C2
            Q  = P**2 + W
            Z  = SQRT( ABS( Q ) )
            AAD( N, N ) = X + T
            X  = AAD( N, N )
            AAD( N1, N1 ) = Y + T
!                                        ** Real pair
            Z  = P + SIGN( Z, P )
            EVALD( N1 ) = X + Z
            EVALD( N ) = EVALD( N1 )

            IF( Z.NE.rZERO ) EVALD( N ) = X - W / Z

            X  = AAD( N, N1 )
!                                  ** Employ scale factor in case
!                                     X and Z are very small
            R  = SQRT( X*X + Z*Z )
            P  = X / R
            Q  = Z / R
!                                             ** Row modification
            DO J = N1, M
               Z  = AAD( N1, J )
               AAD( N1, J ) = Q*Z + P*AAD( N, J )
               AAD( N, J ) = Q*AAD( N, J ) - P*Z
            ENDDO
!                                             ** Column modification
            DO I = 1, N
               Z  = AAD( I, N1 )
               AAD( I, N1 ) = Q*Z + P*AAD( I, N )
               AAD( I, N ) = Q*AAD( I, N ) - P*Z
            ENDDO
!                                          ** Accumulate transformations
            DO I = L, K
               Z  = EVECD( I, N1 )
               EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
               EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
            ENDDO

            N  = N2
            CYCLE OUTER_BLK
         END IF

         IF( IN == 30 ) THEN
      !   ** No convergence after 30 iterations; set error
      !      indicator to the index of the current eigenvalue
            IER  = N
            EVAL(1:M) = EVALD(1:M)
            DO K = 1, M
               EVEC(1:M,K) = EVECD(1:M,K)
            ENDDO
            RETURN
!                                                  ** Form shift
         ELSEIF( IN == 10 .OR. IN == 20 ) THEN
            T  = T + X
            DO I = L, N
               AAD( I, I ) = AAD( I, I ) - X
            ENDDO
            S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
            X  = C3*S
            Y  = X
            W  = -C1*S**2
         END IF

         IN  = IN + 1
      !  ** Look for two consecutive small sub-diagonal elements
         DO J = LB, N2
            I  = N2 + LB - J
            Z  = AAD( I, I )
            R  = X - Z
            S  = Y - Z
            P  = ( R*S - W ) / AAD( I + 1, I ) + AAD( I, I + 1 )
            Q  = AAD( I + 1, I + 1 ) - Z - R - S
            R  = AAD( I + 2, I + 1 )
            S  = ABS( P ) + ABS( Q ) + ABS( R )
            P  = P / S
            Q  = Q / S
            R  = R / S

            IF( I == LB ) EXIT
            UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
            VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) + &
                                ABS( AAD( I+1, I+1 ) ) )
            IF( UU <= TOL*VV ) EXIT
         ENDDO

         AAD( I+2, I ) = rZERO
         DO J = I + 3, N
            AAD( J, J - 2 ) = rZERO
            AAD( J, J - 3 ) = rZERO
         ENDDO

      ! ** Double QR step involving rows K to N and columns M to N
         DO KA = I, N1
            NOTLAS = KA /= N1
            IF( KA == I ) THEN
               S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
               IF( LB /= I ) AAD(KA,KA - 1) = -AAD(KA,KA - 1)
            ELSE
               P  = AAD( KA, KA - 1 )
               Q  = AAD( KA + 1, KA - 1 )
               R  = rZERO
               IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )
               X  = ABS( P ) + ABS( Q ) + ABS( R )
               IF( X == rZERO ) CYCLE
                  P  = P / X
                  Q  = Q / X
                  R  = R / X
                  S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )
                  AAD( KA, KA - 1 ) = -S*X
               END IF

               P  = P + S
               X  = P / S
               Y  = Q / S
               Z  = R / S
               Q  = Q / P
               R  = R / P
         ! ** Row modification
               DO J = KA, M
                  P  = AAD( KA, J ) + Q*AAD( KA + 1, J )
                  IF( NOTLAS ) THEN
                     P  = P + R*AAD( KA + 2, J )
                     AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z
                  END IF
                  AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
                  AAD( KA, J ) = AAD( KA, J ) - P*X
               ENDDO
         ! ** Column modification
               DO II = 1, MIN( N, KA + 3 )
                  P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )
                  IF( NOTLAS ) THEN
                     P  = P + Z*AAD( II, KA + 2 )
                     AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R
                  END IF
                  AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
                  AAD( II, KA ) = AAD( II, KA ) - P
               ENDDO
         ! ** Accumulate transformations
               DO II = L, K
                  P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )
                  IF( NOTLAS ) THEN
                     P  = P + Z*EVECD( II, KA + 2 )
                     EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R
                  END IF
                  EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
                  EVECD( II, KA ) = EVECD( II, KA ) - P
               ENDDO

            ENDDO

         ENDDO INNER_BLK
      ELSE
         EXIT
      ENDIF
      !  ** All evals found, now backsubstitute real vector
      ENDDO OUTER_BLK

      IF( RNORM /= rZERO ) THEN
         DO N = M, 1, -1
            N2   = N
            AAD( N, N ) = rONE
            DO I = N - 1, 1, -1
               W  = AAD( I, I ) - EVALD( N )
               IF( W.EQ.rZERO ) W  = TOL*RNORM
               R  = AAD( I, N )
               DO J = N2, N - 1
                  R  = R + AAD( I, J )*AAD( J, N )
               ENDDO
               AAD( I, N ) = -R / W
               N2   = I
            ENDDO
         ENDDO
!                      ** End backsubstitution vectors of isolated evals
         DO I = 1, M
            IF( I < L .OR. I > K ) THEN
               DO J = I, M
                  EVECD(I,J) = AAD(I,J)
               ENDDO
            END IF
         ENDDO
!                                   ** Multiply by transformation matrix
         IF( K /= 0 ) THEN
            DO J = M, L, -1
               DO I = L, K
                  Z  = rZERO
                  DO N = L, MIN( J, K )
                     Z  = Z + EVECD( I, N )*AAD( N, J )
                  ENDDO
                  EVECD( I, J ) = Z
               ENDDO
            ENDDO
         END IF
      END IF


      DO I = L, K
         DO J = 1, M
            EVECD(I,J) = EVECD(I,J)*WKD(I)
         ENDDO
      ENDDO

!                           ** Interchange rows if permutations occurred
      DO I = L-1, 1, -1
         J  = WKD( I )
         IF( I.NE.J ) THEN
            DO N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
            ENDDO
         END IF
      ENDDO

      DO I = K + 1, M
         J  = WKD( I )
         IF( I /= J ) THEN
            DO N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
            ENDDO
         END IF
      ENDDO
      !  ** Put results into output arrays
      EVAL(1:M) = EVALD(1:M)
      DO K = 1, M
         EVEC(1:M,K) = EVECD(1:M,K)
      ENDDO

      END SUBROUTINE ASYMTX

      SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, PMOM, &
                         NTAU, UTAU, NSTR, NUMU,   &
                         UMU, NPHI, PHI, UMU0,     &
                         FISOT, ALBEDO, HL, TAUC,  &
                         solver_constants )

!           Checks the input dimensions and variables

!   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
!   Called by- DISORT
! --------------------------------------------------------------------

!     .. Scalar Arguments ..

      integer, intent(in) ::   NLYR, NPHI, NSTR, NTAU, NUMU
      real(dk), intent(in)    ::   ALBEDO, FISOT, UMU0
!     ..
!     .. Array Arguments ..

      real(dk), intent(in) ::                 &
                DTAUC(:), HL(0:), PHI(:), &
                PMOM(0:,:), SSALB(:),     &
                TAUC(0:), UMU(:)
      real(dk), intent(inout) :: UTAU(:)
!     ..
!     .. Type Arguments ..

      type(solver_constants_t), intent(in) :: solver_constants

!     ..
!     .. Local Scalars ..

      logical  :: INPERR
      integer  :: IRMU, IU, J, K, LC, LU
      real(dk)     :: FLXALB, RMU

!     ..

      INPERR = .FALSE.

      IF( NLYR < 1 ) INPERR = WRTBAD( 'NLYR' )

      IF( NLYR > MXCLY ) INPERR = WRTBAD( 'MAXCLY' )

      IF( ANY( DTAUC(:) < rZERO ) ) THEN
         INPERR = WRTBAD( 'DTAUC' )
      ENDIF
      IF( ANY( SSALB(:) < rZERO ) .or. ANY( SSALB(:) > rONE ) ) THEN
         INPERR = WRTBAD( 'SSALB' )
      ENDIF
      DO LC = 1, NLYR
         IF( ANY( PMOM(:,LC) < -rONE ) .or.  &
             ANY( PMOM(:,LC) > rONE ) ) THEN
            INPERR = WRTBAD( 'PMOM' )
         ENDIF
      ENDDO

      IF( MXULV < NLYR + 1 ) INPERR = WRTBAD( 'MAXULV' )

      IF( NSTR < 2 .OR. MOD(NSTR,2) /= 0 ) INPERR = WRTBAD( 'NSTR' )


      IF( NSTR > MXCMU ) INPERR = WRTBAD( 'MAXCMU' )

      IF( USRANG ) THEN

         IF( NUMU < 0 ) INPERR = WRTBAD( 'NUMU' )

         IF( .NOT. ONLYFL .AND. NUMU == 0 ) INPERR = WRTBAD( 'NUMU' )

         IF( NUMU > MXUMU ) INPERR = WRTBAD( 'MXUMU' )

         IF( IBCND == 1 .AND. 2*NUMU > MXUMU ) INPERR = WRTBAD( 'MXUMU' )

         DO IU = 1, NUMU
            IF( UMU(IU) < -rONE .OR. UMU(IU) > rONE .OR. &
                UMU(IU) == rZERO ) INPERR = WRTBAD( 'UMU' )

            IF( IBCND == 1 .AND. UMU(IU) < rZERO ) INPERR = WRTBAD( 'UMU' )

            IF( IU > 1 ) THEN
               IF( UMU(IU) < UMU(IU-1) ) INPERR = WRTBAD( 'UMU' )
            END IF
         ENDDO

      ELSE

         IF( MXUMU < NSTR ) INPERR = WRTBAD( 'MAXUMU' )

      END IF


      IF( .NOT. ONLYFL .AND. IBCND /= 1 ) THEN

         IF( NPHI <= 0 ) INPERR = WRTBAD( 'NPHI' )

         IF( NPHI > MXPHI ) INPERR = WRTBAD( 'MAXPHI' )

         IF( ANY( PHI(:) < rZERO ) .OR. ANY( PHI(:) > 360.0 ) ) &
                INPERR = WRTBAD( 'PHI' )

      END IF


      IF( IBCND.LT.0 .OR. IBCND.GT.1 ) INPERR = WRTBAD( 'IBCND' )

      IF( IBCND.EQ.0 ) THEN

         IF( FBEAM.LT.rZERO ) INPERR = WRTBAD( 'FBEAM' )

         IF( FBEAM.GT.rZERO .AND. abs(UMU0).GT.rONE ) INPERR = WRTBAD( 'UMU0' )

         IF( FBEAM.GT.rZERO .AND. ( PHI0.LT.rZERO .OR.PHI0.GT.360.0 ) ) &
             INPERR = WRTBAD( 'PHI0' )

         IF( FISOT.LT.rZERO ) INPERR = WRTBAD( 'FISOT' )

!                    ** Make sure flux albedo at dense mesh of incident
!                       angles does not assume unphysical values
         IF( (.NOT. ONLYFL .AND. USRANG) .OR. .NOT. LAMBER ) THEN
           DO IRMU = 0, 100
             RMU  = REAL(IRMU,dk)*0.01
             FLXALB = DREF(RMU,HL,NSTR,solver_constants)
             IF( FLXALB < rZERO .OR. FLXALB > rONE ) INPERR = WRTBAD( 'HL' )
           ENDDO
         ENDIF

      ELSE IF( IBCND == 1 ) THEN

         IF( ALBEDO < rZERO .OR. ALBEDO > rONE ) INPERR = WRTBAD( 'ALBEDO' )

      END IF

      IF( ACCUR < rZERO .OR. ACCUR > 1.E-2 ) &
        INPERR = WRTBAD( 'ACCUR' )

      IF( MXCLY.LT.NLYR ) &
        INPERR = WRTDIM( 'MXCLY', NLYR )

      IF( IBCND /= 1 ) THEN
         IF( USRTAU .AND. MXULV.LT.NTAU ) &
           INPERR = WRTDIM( 'MXULV',NTAU )

         IF( .NOT.USRTAU .AND. MXULV .LT. NLYR + 1 ) &
           INPERR = WRTDIM( 'MXULV', NLYR + 1 )
      ELSE
         IF( MXULV.LT.2 ) INPERR = WRTDIM( 'MXULV', 2 )
      END IF

      IF( MXCMU.LT.NSTR ) &
        INPERR = WRTDIM( 'MXCMU', NSTR )

      IF( USRANG .AND. MXUMU.LT.NUMU ) &
        INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( USRANG .AND. IBCND.EQ.1 .AND.MXUMU.LT.2*NUMU ) &
        INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( .NOT.USRANG .AND. MXUMU.LT.NSTR ) &
        INPERR = WRTDIM( 'MXUMU', NSTR )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 .AND. MXPHI.LT.NPHI ) &
        INPERR = WRTDIM( 'MXPHI', NPHI )

      IF( INPERR ) &
        CALL ERRMSG( 'DISORT--input and/or dimension errors',.True.)

      END SUBROUTINE CHEKIN

      SUBROUTINE FLUXES( tausla, tauslau,                               &
                         CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,    &
                         NCUT, NN, NSTR, NTAU, PI,                      &
                         PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, XR0,  &
                         XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, FLDN, FLDIR,&
                         RFLDIR, RFLDN, UAVG, U0C,                      &
                         uavgso, uavgup, uavgdn,                        &
                         sindir, sinup, sindn)

!       Calculates the radiative fluxes, mean intensity, and flux
!       derivative with respect to optical depth from the m=0 intensity
!       components (the azimuthally-averaged intensity)

!    I N P U T     V A R I A B L E S:

!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
!       LAYRU    :  Layer number of user level UTAU
!       LL       :  Constants of integration in Eq. SC(1), obtained
!                     by solving scaled version of Eq. SC(5);
!                     exponential term of Eq. SC(12) not included
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
!       NCUT     :  Number of computational layer where absorption
!                     optical depth exceeds ABSCUT
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       UTAUPR   :  Optical depths of user output levels in delta-M
!                     coordinates;  equal to UTAU if no delta-M
!       XR0      :  Expansion of thermal source function in Eq. SS(14)
!       XR1      :  Expansion of thermal source function Eqs. SS(16)
!       ZZ       :  Beam source vectors in Eq. SS(19)
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
!       (remainder are DISORT input variables)


!                   O U T P U T     V A R I A B L E S:

!       U0C      :  Azimuthally averaged intensities
!                   ( at polar quadrature angles )
!       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)


!                   I N T E R N A L       V A R I A B L E S:

!       DIRINT   :  Direct intensity attenuated
!       FDNTOT   :  Total downward flux (direct + diffuse)
!       FLDIR    :  Direct-beam flux (delta-M scaled)
!       FLDN     :  Diffuse down-flux (delta-M scaled)
!       FNET     :  Net flux (total-down - diffuse-up)
!       FACT     :  EXP( - UTAUPR / UMU0 )
!       PLSORC   :  Planck source function (thermal)
!       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

!   Called by- DISORT
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      integer, intent(in) :: NCUT, NN, NSTR, NTAU
      real(dk), intent(in)    :: FBEAM, PI, UMU0
      logical, intent(in) :: LYRCUT
!     ..
!     .. Array Arguments ..

      logical, intent(in) :: PRNT(:)
      integer, intent(in) :: LAYRU(:)
      real(dk),    intent(in) :: CMU(:), CWT(:),        &
                GC(:,:,:), KK(:,:), LL(:,:),        &
                SSALB(:),  TAUCPR(0:),              &
                UTAU(:), UTAUPR(:), XR0(:), XR1(:), &
                ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)
      real(dk), intent(in)  :: tausla(0:), tauslau(0:)
      real(dk), intent(out) :: U0C(:,:)
      real(dk), intent(out) :: RFLDIR(:), RFLDN(:), FLUP(:)
      real(dk), intent(out) :: DFDT(:), UAVG(:)
      real(dk), intent(out) :: uavgso(:), uavgup(:), uavgdn(:)
      real(dk), intent(out) :: sindir(:), sinup(:), sindn(:)
      real(dk), intent(out) :: FLDIR(:), FLDN(:)
!     ..
!     .. Local Scalars ..

      integer :: IQ, JQ, LU, LYU
      real(dk)    :: ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT
!     ..

      IF( PRNT( 2 ) ) WRITE( *, 9000 )
!                                          ** Zero DISORT output arrays
      U0C   = 0.
      FLDIR = 0.
      FLDN  = 0.
      uavgso = 0.
      uavgup = 0.
      uavgdn = 0.
      sindir = 0.
      sinup  = 0.
      sindn  = 0.

!    ** Loop over user levels
      LEVEL_LOOP: DO LU = 1, NTAU

         LYU  = LAYRU( LU )

         IF( LYRCUT .AND. LYU > NCUT ) THEN
!                                                ** No radiation reaches
!                                                ** this level
            FDNTOT = rZERO
            FNET   = rZERO
            PLSORC = rZERO
            IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU,  &
               RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET, &
               UAVG( LU ), PLSORC, DFDT( LU )
            CYCLE
         END IF

         IF( FBEAM > rZERO ) THEN
 
            FACT  = EXP( - tausla(LU-1) )
            DIRINT       = FBEAM*FACT
            FLDIR( LU )  = UMU0*( FBEAM*FACT )
            RFLDIR( LU ) = UMU0*FBEAM * EXP( -tauslau(lu-1) )
            sindir( LU ) = SQRT(1.-UMU0*UMU0)*FBEAM * EXP( -tauslau(lu-1) )

         ELSE

            DIRINT       = rZERO
            FLDIR( LU )  = rZERO
            RFLDIR( LU ) = rZERO
            sindir( LU ) = rZERO

         END IF


         DO IQ = 1, NN

            ZINT   = rZERO

            DO JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU ) &
                        * EXP( -KK( JQ,LYU )*( UTAUPR( LU ) - TAUCPR( LYU ) ) )
            ENDDO

            DO JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )* &
                        EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -     &
                        TAUCPR( LYU - 1 ) ) )
            ENDDO

            U0C( IQ, LU ) = ZINT

            IF( FBEAM > rZERO ) THEN
              U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
            ENDIF

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) &
                            + ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( NN + 1 - IQ )*U0C( IQ, LU )
            uavgdn(lu) = uavgdn(lu) + cwt(nn+1-iq) * u0c( iq,lu )
            sindn(lu)  = sindn(lu)  + cwt(nn+1-iq)           &
                         *SQRT(1.-CMU(NN+1-IQ)*CMU(NN+1-IQ)) &
                         *U0C( IQ, LU )
            FLDN( LU ) = FLDN( LU ) + CWT( NN + 1 - IQ ) &
                         *CMU( NN + 1 - IQ )*U0C( IQ, LU )
         ENDDO


         DO IQ = NN + 1, NSTR

            ZINT   = rZERO

            DO JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU ) &
                        *EXP( -KK(JQ,LYU)*(UTAUPR(LU) - TAUCPR(LYU)) )
            ENDDO

            DO JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU ) &
                        *EXP( -KK( JQ,LYU )*(UTAUPR(LU) - TAUCPR(LYU - 1)) )
            ENDDO

            U0C( IQ, LU ) = ZINT

            IF( FBEAM > rZERO ) THEN
              U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
            ENDIF

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) &
                            + ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( IQ - NN )*U0C( IQ, LU )
            uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c( iq,lu )
            sinup (lu) = sinup(lu)  + cwt(iq-nn) & 
                         *SQRT(1. - CMU(IQ-NN)*CMU(IQ-NN))*U0C(IQ,LU)
            FLUP( LU ) = FLUP( LU ) + CWT( IQ - NN )*CMU( IQ - NN )* U0C( IQ, LU )
         ENDDO


         FLUP( LU )  = rTWO*PI*FLUP( LU )
         FLDN( LU )  = rTWO*PI*FLDN( LU )
         FDNTOT      = FLDN( LU ) + FLDIR( LU )
         FNET        = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU )  = ( rTWO*PI*UAVG( LU ) + DIRINT ) / ( rFOUR*PI )
         uavgso( lu ) = dirint / (rFOUR*pi)
         uavgup( lu ) = (rTWO * pi * uavgup(lu) )/ (rFOUR*pi)
         uavgdn( lu)  = (rTWO * pi * uavgdn(lu) )/ (rFOUR*pi)
         sindn ( lu ) = rTWO*PI*sindn ( LU )
         sinup ( lu ) = rTWO*PI*sinup ( LU )

         PLSORC      = XR0( LYU ) + XR1( LYU )*UTAUPR( LU )
         DFDT( LU )  = ( 1.- SSALB( LYU ) ) * rFOUR*PI *(UAVG(LU) - PLSORC)

         IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU,  &
             RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET, &
             UAVG( LU ), PLSORC, DFDT( LU )

      ENDDO LEVEL_LOOP


      IF( PRNT( 3 ) ) THEN

         WRITE( *, FMT = 9020 )

         DO LU = 1, NTAU

            WRITE( *, FMT = 9030 ) UTAU( LU )

            DO IQ = 1, NN
               ANG1   = 180./ PI* ACOS( CMU( 2*NN - IQ + 1 ) )
               ANG2   = 180./ PI* ACOS( CMU( IQ ) )
               WRITE( *, 9040 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),  &
                                ANG2, CMU(IQ),        U0C(IQ+NN,LU)
            ENDDO
         ENDDO

      END IF


 9000 FORMAT( //, 21X,                                                   &
       '<----------------------- FLUXES ----------------------->', /,    &
       '   Optical  Compu    Downward    Downward    Downward     ',     &
       ' Upward                    Mean      Planck   d(Net Flux)', /,   &
       '     Depth  Layer      Direct     Diffuse       Total     ',     &
       'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )
 9010 FORMAT( F10.4, I7, 1P, 7E12.3, E14.3 )
 9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES', &
            ' ( at polar quadrature angles ) *******' )
 9030 FORMAT( /, ' Optical depth =', F10.4, //,        &
        '     Angle (deg)   cos(Angle)     Intensity', &
        '     Angle (deg)   cos(Angle)     Intensity' )
 9040 FORMAT( 2( 0P,F16.4,F13.5,1P,E14.3 ) )

      END SUBROUTINE FLUXES

      SUBROUTINE LEPOLY( NMU, M, TWONM1, MU, YLM, solver_constants )

!       Computes the normalized associated Legendre polynomial,
!       defined in terms of the associated Legendre polynomial
!       Plm = P-sub-l-super-m as

!             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

!       for fixed order m and all degrees from l = m to TWONM1.
!       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
!       from a prior call to the routine.

!       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
!                  High-Order Associated Legendre Polynomials,
!                  J. Quant. Spectrosc. Radiat. Transfer 10,
!                  557-562, 1970.  (hereafter D/A)

!       METHOD: Varying degree recurrence relationship.

!       NOTE 1: The D/A formulas are transformed by
!               setting  M = n-1; L = k-1.
!       NOTE 2: Assumes that routine is called first with  M = 0,
!               then with  M = 1, etc. up to  M = TWONM1.
!       NOTE 3: Loops are written in such a way as to vectorize.

!  I N P U T     V A R I A B L E S:

!       NMU    :  Number of arguments of YLM
!       M      :  Order of YLM
!       MAXMU  :  First dimension of YLM
!       TWONM1 :  Max degree of YLM
!       MU(i)  :  Arguments of YLM (i = 1 to NMU)

!       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
!       from a prior call.

!  O U T P U T     V A R I A B L E:

!       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
!                   polynomials evaluated at argument MU(i)

!   Called by- DISORT, ALBTRN, SURFAC
!   Calls- ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      integer, intent(in) :: M, NMU, TWONM1
!     ..
!     .. Array Arguments ..

      real(dk), intent(in)  :: MU(:)
      real(dk), intent(out) :: YLM(0:,:)
!     ..
!     .. Type Arguments ..
      type(solver_constants_t), intent(in) :: solver_constants
!     ..
!     .. Local Scalars ..

      integer   :: I, L, NS
      real(dk)      :: TMP1, TMP2
!     ..


      IF( 2*TWONM1 > MAXSQT ) &
        CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)


      IF( M == 0 ) THEN
!                             ** Upward recurrence for ordinary
!                                Legendre polynomials
         DO I = 1, NMU
            YLM(0,I) = rONE
            YLM(1,I) = MU(I)
         ENDDO

         DO L = 2, TWONM1
            DO I = 1, NMU
               YLM(L,I) = ((2*L - 1)*MU(I)*YLM(L - 1,I) &
                                 - (L - 1)*YLM(L - 2,I)) / L
            ENDDO
         ENDDO
      ELSE
         DO I = 1, NMU
!                               ** Y-sub-m-super-m; derived from
!                               ** D/A Eqs. (11,12)
            YLM(M,I) = -solver_constants%SQT(2*M - 1) / &
              solver_constants%SQT(2*M )*SQRT(rONE - MU(I)**2)*YLM(M - 1,I)

!                              ** Y-sub-(m+1)-super-m; derived from
!                              ** D/A Eqs.(13,14) using Eqs.(11,12)
            YLM(M + 1,I) = solver_constants%SQT(2*M + 1)*MU(I)*YLM(M,I)
         ENDDO
!                                   ** Upward recurrence; D/A EQ.(10)
         DO L = M + 2, TWONM1
            TMP1 = solver_constants%SQT(L - M )*solver_constants%SQT(L + M)
            TMP2 = solver_constants%SQT(L - M - 1)*solver_constants%SQT(L + M - 1)
            DO I = 1, NMU
               YLM(L,I) = ((2*L - 1 )*MU(I)*YLM(L-1,I) - TMP2*YLM(L-2,I)) / TMP1
            ENDDO
         ENDDO
      END IF

      END SUBROUTINE LEPOLY

      SUBROUTINE PRAVIN( UMU, NUMU, UTAU, NTAU, U0U )

!        Print azimuthally averaged intensities at user angles

!   Called by- DISORT

!     LENFMT   Max number of polar angle cosines UMU that can be
!                printed on one line, as set in FORMAT statement
! --------------------------------------------------------------------

!     .. Scalar Arguments ..

      integer   NTAU, NUMU
!     ..
!     .. Array Arguments ..

      real(dk)      U0U( :, : ), UMU( : ), UTAU( : )
!     ..
!     .. Local Scalars ..

      integer   IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MIN
!     ..


      IF( NUMU.LT.1 )  RETURN

      WRITE( *, '(//,A)' )                                &
         ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' // &
         '(at user polar angles)  ********'

      LENFMT = 8
      NPASS  = 1 + (NUMU-1) / LENFMT

      WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines', &
                            '     Depth'

      DO NP = 1, NPASS
         IUMIN  = 1 + LENFMT * ( NP - 1 )
         IUMAX  = MIN( LENFMT*NP, NUMU )
         WRITE( *,'(/,10X,8F14.5)') ( UMU(IU), IU = IUMIN, IUMAX )

         DO LU = 1, NTAU
            WRITE( *, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ), &
                 ( U0U( IU,LU ), IU = IUMIN, IUMAX )
         ENDDO
      ENDDO

      END SUBROUTINE PRAVIN

      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM,  &
                         NTAU, UTAU, NSTR, NUMU, UMU,       &
                         PHI, UMU0, FISOT,                  &
                         ALBEDO, HL, FLYR, LYRCUT,          &
                         OPRIM, TAUC, TAUCPR, PRTMOM )

!        Print values of input variables

!   Called by- DISORT
! --------------------------------------------------------------------

!     .. Scalar Arguments ..

      logical, intent(in) :: LYRCUT, PRTMOM
      integer, intent(in) :: NLYR, NSTR, NTAU, NUMU
      real(dk), intent(in)    :: ALBEDO, FISOT, UMU0
!     ..
!     .. Array Arguments ..

      real(dk), intent(in) ::  DTAUC(:), DTAUCP(:), FLYR(:), HL(0:),  &
                OPRIM(:), PHI(:), PMOM(0:,:), SSALB(:),           &
                TAUC(0:), TAUCPR(0:), UMU(:), UTAU(:)
!     ..
!     .. Local Scalars ..

      integer :: IU, J, K, LC, LU
      real(dk)    :: YESSCT
!     ..


      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,  &
             '     No. computational layers =', NLYR

      IF( IBCND /= 1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' )  &
          NTAU,' User optical depths :', ( UTAU(LU), LU = 1, NTAU )

      IF( .NOT. ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' ) &
          NUMU,' User polar angle cosines :',( UMU(IU), IU = 1, NUMU )

      IF( .NOT. ONLYFL .AND. IBCND /= 1 ) &
          WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' ) &
                 NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )

      IF( .NOT. PLANK .OR. IBCND == 1 ) WRITE( *, '(A)' ) ' No thermal emission'


      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND

      IF( IBCND == 0 ) THEN

         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' ) &
                '    Incident beam with intensity =', FBEAM,      &
                ' and polar angle cosine = ', UMU0,               &
                '  and azimuth angle =', PHI0,                    &
                '    plus isotropic incident intensity =', FISOT

         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' ) '    Bottom albedo (Lambertian) =', ALBEDO

         IF( .NOT. LAMBER ) WRITE( *, '(A,/,(10X,10F9.5))' )              &
           '    Legendre coeffs of bottom bidirectional reflectivity :',  &
               ( HL( K ), K = 0, NSTR )

      ELSE IF( IBCND == 1 ) THEN

         WRITE(*,'(A)') '    Isotropic illumination from top and bottom'
         WRITE( *, '(A,0P,F8.4)' ) '    Bottom albedo (Lambertian) =', ALBEDO
      END IF


      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'
      IF( .NOT.DELTAM ) WRITE( *, '(A)' ) ' Does not use delta-M method'


      IF( IBCND == 1 ) THEN

         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of medium vs. incident beam angle'

      ELSE IF( ONLYFL ) THEN

         WRITE( *, '(A)' ) ' Calculate fluxes and azim-averaged intensities only'

      ELSE

         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

      END IF


      WRITE( *, '(A,1P,E11.2)' ) ' Relative convergence criterion for azimuth series =',ACCUR

      IF( LYRCUT ) WRITE( *, '(A)' ) ' Sets radiation = 0 below absorption optical depth 10'


!                                        ** Print layer variables
      IF( PLANK ) WRITE( *, FMT = 9180 )
      IF( .NOT. PLANK ) WRITE( *, FMT = 9190 )

      YESSCT = rZERO

      DO LC = 1, NLYR
         YESSCT = YESSCT + SSALB( LC )

         IF( PLANK ) &
             WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')   &
                   LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),  &
                   DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM(1,LC)

         IF( .NOT.PLANK )  &
             WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)')         &
                   LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),  &
                   DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
      ENDDO


      IF( PRTMOM .AND. YESSCT > rZERO ) THEN

         WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'

         DO LC = 1, NLYR
            IF( SSALB( LC ).GT.rZERO )  &
                WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' )  &
                       LC, ( PMOM( K, LC ), K = 0, NSTR )
         ENDDO

      END IF

!                ** (Read every other line in these formats)

 9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,     &
      '                   Total    Single                           ',  &
                     'Total    Single', /,                              &
      '       Optical   Optical   Scatter   Truncated   ',              &
         'Optical   Optical   Scatter    Asymm', /,                     &
      '         Depth     Depth    Albedo    Fraction     ',            &
           'Depth     Depth    Albedo   Factor   Temperature' )
 9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,     &
      '                   Total    Single                           ',  &
                     'Total    Single', /,                              &
      '       Optical   Optical   Scatter   Truncated   ',              &
         'Optical   Optical   Scatter    Asymm', /,                     &
      '         Depth     Depth    Albedo    Fraction     ',            &
           'Depth     Depth    Albedo   Factor' )

      END SUBROUTINE PRTINP

      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )

!         Prints the intensity at user polar and azimuthal angles

!     All arguments are DISORT input or output variables

!   Called by- DISORT

!     LENFMT   Max number of azimuth angles PHI that can be printed
!                on one line, as set in FORMAT statement
! +-------------------------------------------------------------------+


!     .. Scalar Arguments ..

      integer   NPHI, NTAU, NUMU
!     ..
!     .. Array Arguments ..

      real(dk)      PHI(:), UMU(:), UTAU(:), UU(:,:,:)
!     ..
!     .. Local Scalars ..

      integer   IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS
!     ..
!     .. Intrinsic Functions ..

      INTRINSIC MIN
!     ..


      IF( NPHI.LT.1 )  RETURN

      WRITE( *, '(//,A)' ) ' *********  I N T E N S I T I E S  *********'

      LENFMT = 10
      NPASS  = 1 + (NPHI-1) / LENFMT

      WRITE( *, '(/,A,/,A,/,A)' )  &
         '             Polar   Azimuth angles (degrees)',  &
         '   Optical   Angle', '    Depth   Cosine'

      DO LU = 1, NTAU
         DO NP = 1, NPASS

            JMIN   = 1 + LENFMT * ( NP - 1 )
            JMAX   = MIN( LENFMT*NP, NPHI )

            WRITE( *, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )

            IF( NP.EQ.1 ) &
              WRITE( *, '(F10.4,F8.4,1P,10E11.3)' )  &
                   UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            IF( NP.GT.1 ) &
              WRITE( *, '(10X,F8.4,1P,10E11.3)' )  &
                             UMU(1), (UU(1, LU, J), J = JMIN, JMAX)

            DO IU = 2, NUMU
               WRITE( *, '(10X,F8.4,1P,10E11.3)' )  &
                 UMU( IU ), ( UU( IU, LU, J ), J = JMIN, JMAX )
            END DO
         END DO
      END DO

      END SUBROUTINE PRTINT

      SUBROUTINE QGAUSN( M, GMU, GWT )

!       Compute weights and abscissae for ordinary Gaussian quadrature
!       on the interval (0,1);  that is, such that

!           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

!       is a good approximation to

!           integral(0 to 1) ( f(x) dx )

!   INPUT :    M       order of quadrature rule

!   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
!             GWT(I)   array of weights (I = 1 TO M)

!   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
!                   Integration, Academic Press, New York, pp. 87, 1975

!   METHOD:  Compute the abscissae as roots of the Legendre
!            polynomial P-sub-M using a cubically convergent
!            refinement of Newton's method.  Compute the
!            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
!            that Newton's method can very easily diverge; only a
!            very good initial guess can guarantee convergence.
!            The initial guess used here has never led to divergence
!            even for M up to 1000.

!   ACCURACY:  relative error no better than TOL or computer
!              precision (machine epsilon), whichever is larger

!   INTERNAL VARIABLES:

!    ITER      : number of Newton Method iterations
!    MAXIT     : maximum allowed iterations of Newton Method
!    PM2,PM1,P : 3 successive Legendre polynomials
!    PPR       : derivative of Legendre polynomial
!    P2PRI     : 2nd derivative of Legendre polynomial
!    TOL       : convergence criterion for Legendre poly root iteration
!    X,XI      : successive iterates in cubically-convergent version
!                of Newtons Method (seeking roots of Legendre poly.)

!   Called by- SETDIS, SURFAC
!   Calls- ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      integer, intent(in) :: M
!     ..
!     .. Array Arguments ..

      real(dk), intent(out)  :: GMU(:), GWT(:)
!     ..
!     .. Local Scalars ..

      integer, parameter ::  MAXIT = 1000

      integer ::  ITER, K, LIM, NN, NP1
      real(dk)    ::  CONA, PI, T
      real(dk) ::  EN, NNP1, P, P2PRI, PM1, PM2, PPR, PROD
      real(dk) ::  TMP, TOL, X, XI


      PI   = rTWO*ASIN( rONE )
      TOL  = rTEN*EPSILON( TOL )

      IF( M < 1 ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)

      IF( M == 1 ) THEN
         GMU( 1 ) = ONEHALF
         GWT( 1 ) = rONE
         RETURN
      END IF

      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = REAL( M - 1,dk ) / REAL( 8*M**3,dk )

      LIM  = M / 2

      DO K = 1, LIM
!                                        ** Initial guess for k-th root
!                                           of Legendre polynomial, from
!                                           Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
         X  = COS( T + CONA / TAN( T ) )
         ITER = 0
!                                        ** Upward recurrence for
!                                           Legendre polynomials
UNC_LOOP: &
         DO WHILE( .true. )
            ITER   = ITER + 1
            PM2    = rONE
            PM1    = X

            DO NN = 2, M
               P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / real(NN,dk)
               PM2  = PM1
               PM1  = P
            ENDDO
!                                              ** Newton Method
            TMP    = rONE / ( rONE - X**2 )
            PPR    = EN*( PM2 - X*P )*TMP
            P2PRI  = (rTWO*X*PPR - NNP1*P)*TMP
            XI     = X - ( P / PPR )*( rONE + (P / PPR)*P2PRI / (rTWO*PPR) )

!                                              ** Check for convergence
            IF( ABS( XI - X ) > TOL ) THEN
               IF( ITER.GT.MAXIT ) &
                  CALL ERRMSG( 'QGAUSN--max iteration count',.True.)
               X  = XI
               CYCLE UNC_LOOP
            END IF
            EXIT
         END DO UNC_LOOP
!                             ** Iteration finished--calculate weights,
!                                abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = rTWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
      ENDDO
!                                    ** Set middle abscissa and weight
!                                       for rules of odd order
      IF( MOD( M,2 ).NE.0 ) THEN

         GMU( LIM + 1 ) = rZERO
         PROD   = rONE

         DO K = 3, M, 2
            PROD   = PROD * K / ( K - 1 )
         ENDDO

         GWT( LIM + 1 ) = rTWO / PROD**2
      END IF

!                                        ** Convert from (-1,1) to (0,1)
      DO K = 1, M
         GMU( K ) = ONEHALF*GMU( K ) + ONEHALF
         GWT( K ) = ONEHALF*GWT( K )
      ENDDO

      END SUBROUTINE QGAUSN

      SUBROUTINE SETDIS( dsdh, nid, tausla, tauslau, mu2,                 &
                         CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM,  &
                         FLYR, GL, HL, HLPR, IBCND, LAMBER, LAYRU,        &
                         LYRCUT, NCUT, NLYR,                              &
                         NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, OPRIM,      &
                         PMOM, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU,    &
                         UMU0, USRTAU, USRANG )

!          Perform miscellaneous setting-up operations

!       INPUT :  all are DISORT input variables (see DOC file)

!       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
!                NUMU,UMU    if USRANG = FALSE
!                CMU,CWT     computational polar angles and
!                               corresponding quadrature weights
!                EXPBEA      transmission of direct beam
!                FLYR        truncated fraction in delta-M method
!                GL          phase function Legendre coefficients multi-
!                              plied by (2L+1) and single-scatter albedo
!                HLPR        Legendre moments of surface bidirectional
!                              reflectivity, times 2K+1
!                LAYRU       Computational layer in which UTAU falls
!                LYRCUT      flag as to whether radiation will be zeroed
!                              below layer NCUT
!                NCUT        computational layer where absorption
!                              optical depth first exceeds  ABSCUT
!                NN          NSTR / 2
!                OPRIM       delta-M-scaled single-scatter albedo
!                TAUCPR      delta-M-scaled optical depth
!                UTAUPR      delta-M-scaled version of  UTAU

!   Called by- DISORT
!   Calls- QGAUSN, ERRMSG
! ----------------------------------------------------------------------

!     use tuv_params, only : largest

!     .. Scalar Arguments ..

      integer, intent(in)  ::  IBCND, NLYR, NSTR
      integer, intent(out) ::  NCUT, NN, NTAU, NUMU
      logical, intent(in)  ::  DELTAM, LAMBER, ONLYFL
      logical, intent(in)  ::  PLANK, USRANG, USRTAU
      logical, intent(out) ::  LYRCUT
      real(dk), intent(in)     ::  FBEAM, UMU0

! geometry
      integer, intent(in) :: nid(0:)
      real(dk), intent(in)    :: dsdh(0:,:)
      real(dk), intent(out)   :: tausla(0:), tauslau(0:), mu2(0:)

      real(dk), parameter     :: largest = 1.e36_dk
      real(dk) :: sum, sumu
!     ..
!     .. Array Arguments ..

      integer, intent(out) :: LAYRU(:)
      real(dk), intent(in)   :: DTAUC(:), HL(0:), SSALB(:), TAUC(0:)
      real(dk), intent(out)  :: UTAU(:), UTAUPR(:)
      real(dk), intent(out)  :: CMU(:), CWT(:), DTAUCP(:), EXPBEA(0:)
      real(dk), intent(out)  :: FLYR(:), GL(0:,:), HLPR(0:), UMU(:)
      real(dk), intent(out)  :: OPRIM(:), PMOM(0:,:), TAUCPR(0:) 

!     ..
!     .. Local Scalars ..

      real(dk), parameter :: ABSCUT = 10000._dk

      integer   :: IQ, IU, K, LC, LU, I
      real(dk)      :: ABSTAU, F

      IF( .NOT. USRTAU ) THEN
!  ** Set output levels at computational layer boundaries
         NTAU  = NLYR + 1
         DO LC = 0, NTAU - 1
            UTAU(LC + 1) = TAUC(LC)
         ENDDO
      END IF
!                        ** Apply delta-M scaling and move description
!                           of computational layers to local variables
      EXPBEA(0) = rONE
      TAUCPR(0) = rZERO
      ABSTAU    = rZERO

      tausla  = rZERO
      tauslau = rZERO
      mu2 = rONE/largest

      DO LC = 1, NLYR

         PMOM( 0, LC ) = rONE

         IF( ABSTAU < ABSCUT ) NCUT  = LC

         ABSTAU = ABSTAU + (rONE - SSALB(LC))*DTAUC(LC)

         IF( .NOT. DELTAM ) THEN
            OPRIM( LC )  = SSALB( LC )
            DTAUCP( LC ) = DTAUC( LC )
            TAUCPR( LC ) = TAUC( LC )

            DO K = 0, NSTR - 1
               GL(K,LC) = REAL(2*K + 1,dk)*OPRIM(LC)*PMOM(K,LC)
            ENDDO

            F  = rZERO
         ELSE
!                                    ** Do delta-M transformation
            F  = PMOM(NSTR,LC)
            OPRIM(LC) = SSALB(LC) * (rONE - F) / ( rONE - F*SSALB(LC))
            DTAUCP(LC) = ( rONE - F*SSALB(LC))*DTAUC(LC)
            TAUCPR(LC) = TAUCPR(LC - 1) + DTAUCP(LC)

            DO K = 0, NSTR - 1
               GL(K,LC) = REAL(2*K + 1,dk) * OPRIM(LC)*(PMOM(K,LC) - F) / (rONE - F)
            ENDDO
         END IF

         FLYR(LC) = F
         EXPBEA(LC) = rZERO
      ENDDO
! 
! calculate slant optical depth
!              
         IF(umu0 <  rZERO) THEN
           IF(nid(0) < 0) THEN
             tausla(0) = largest
             tauslau(0) = largest
           ELSE
             sum = rZERO
             sumu = rZERO
             DO lc = 1, nid(0)
               sum = sum + rTWO*dtaucp(lc)*dsdh(0,lc)
               sumu = sumu + rTWO*dtauc(lc)*dsdh(0,lc)
             END DO
             tausla(0) = sum 
             tauslau(0) = sumu 
           END IF
         END IF

         expbea(0) = EXP( -tausla(0) )

!
         DO lc = 1, nlyr
          IF(nid(lc) < 0) THEN
            tausla(lc) = largest
            tauslau(lc) = largest
          ELSE
            sum = rZERO
            sumu = rZERO
            DO lu = 1, MIN(nid(lc),lc)
               sum = sum + dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + dtauc(lu)*dsdh(lc,lu)
            ENDDO
            DO lu = MIN(nid(lc),lc)+1,nid(lc)
               sum = sum + rTWO*dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + rTWO*dtauc(lu)*dsdh(lc,lu)
            ENDDO
            tausla(lc) = sum 
            tauslau(lc) = sumu 
            IF(tausla(lc) == tausla(lc-1)) THEN
              mu2(lc) = largest
            ELSE
              mu2(lc) = (taucpr(lc)-taucpr(lc-1)) / (tausla(lc)-tausla(lc-1))
              mu2(lc) = SIGN( MAX(ABS(mu2(lc)),1./largest),mu2(lc) )
            END IF
          END IF
          expbea(lc) = EXP( -tausla(lc) )
         ENDDO

!                      ** If no thermal emission, cut off medium below
!                         absorption optical depth = ABSCUT ( note that
!                         delta-M transformation leaves absorption
!                         optical depth invariant ).  Not worth the
!                         trouble for one-layer problems, though.

      LYRCUT = ABSTAU >= ABSCUT .AND. .NOT. PLANK .AND. IBCND /= 1 .AND. NLYR > 1

      IF( .NOT. LYRCUT ) NCUT = NLYR

!                             ** Set arrays defining location of user
!                             ** output levels within delta-M-scaled
!                             ** computational mesh
      DO LU = 1, NTAU
         DO LC = 1, NLYR
            IF( UTAU(LU) >= TAUC(LC - 1 ) .AND. UTAU(LU) <= TAUC(LC) ) EXIT
         ENDDO
         LC = MIN( NLYR,LC )

         UTAUPR(LU) = UTAU(LU)
         IF( DELTAM ) THEN
           UTAUPR( LU ) = TAUCPR(LC - 1)  &
                        + (rONE - SSALB(LC)*FLYR(LC))*(UTAU(LU) - TAUC(LC-1))
         ENDIF
         LAYRU(LU) = LC
      ENDDO
!                      ** Calculate computational polar angle cosines
!                         and associated quadrature weights for Gaussian
!                         quadrature on the interval (0,1) (upward)
      NN   = NSTR / 2

      CALL QGAUSN( NN, CMU, CWT )
!                                  ** Downward (neg) angles and weights
      DO IQ = 1, NN
         CMU(IQ + NN) = -CMU(IQ)
         CWT(IQ + NN) = CWT(IQ)
      ENDDO


         DO IQ = 1, NN
!                      ** Dither mu2 if it is close to one of the 
!                         quadrature angles.
           DO  lc = 1, nlyr
             IF (  ABS(mu2(lc)) < 1.E5 ) THEN
               IF( ABS(rONE - ABS(mu2(lc))/CMU(IQ)) < 0.05 ) mu2(lc) = mu2(lc)*0.999
             ENDIF
           END DO
         END DO

      IF( .NOT. USRANG .OR. ( ONLYFL .AND. MXUMU >= NSTR ) ) THEN
!                                   ** Set output polar angles to
!                                      computational polar angles
         NUMU = NSTR
         DO IU = 1, NN
            UMU(IU) = -CMU(NN + 1 - IU)
         ENDDO

         DO IU = NN + 1, NSTR
            UMU(IU) = CMU(IU - NN)
         ENDDO
      END IF


      IF( USRANG .AND. IBCND == 1 ) THEN
!                               ** Shift positive user angle cosines to
!                                  upper locations and put negatives
!                                  in lower locations
         DO IU = 1, NUMU
            UMU(IU+NUMU) = UMU(IU)
         ENDDO

         DO IU = 1, NUMU
            UMU(IU) = -UMU(2*NUMU + 1 - IU)
         ENDDO

         NUMU   = 2*NUMU
      END IF

      IF( .NOT. LYRCUT .AND. .NOT. LAMBER ) THEN
         DO K = 0, NSTR
            HLPR(K) = REAL((2*K + 1),dk)*HL(K)
         ENDDO
      END IF

      END SUBROUTINE SETDIS

      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,  &
                         LAMBER, LYRCUT, NCOL, NCUT,                   &
                         NN, NSTR, TAUCPR, WK )

!        Calculate coefficient matrix for the set of equations
!        obtained from the boundary conditions and the continuity-
!        of-intensity-at-layer-interface equations;  store in the
!        special banded-matrix format required by LINPACK routines

!     I N P U T      V A R I A B L E S:

!       BDR      :  Surface bidirectional reflectivity
!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       DELM0    :  Kronecker delta, delta-sub-m0
!       GC       :  Eigenvectors at polar quadrature angles, SC(1)
!       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       NN       :  Number of streams in a hemisphere (NSTR/2)
!       NCUT     :  Total number of computational layers considered
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       (remainder are DISORT input variables)

!   O U T P U T     V A R I A B L E S:

!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                      scaled by Eq. SC(12); in banded form required
!                      by LINPACK solution routines
!       NCOL     :  Counts of columns in CBAND

!   I N T E R N A L    V A R I A B L E S:

!       IROW     :  Points to row in CBAND
!       JCOL     :  Points to position in layer block
!       LDA      :  Row dimension of CBAND
!       NCD      :  Number of diagonals below or above main diagonal
!       NSHIFT   :  For positioning number of rows in band storage
!       WK       :  Temporary storage for EXP evaluations

!   Called by- DISORT, ALBTRN
! +--------------------------------------------------------------------+


!     .. Scalar Arguments ..

      logical, intent(in)  :: LAMBER, LYRCUT
      integer, intent(in)  :: NCUT, NN, NSTR
      integer, intent(out) :: NCOL
      real(dk), intent(in) :: DELM0
!     ..
!     .. Array Arguments ..

      real(dk), intent(in)  ::  &
                BDR(:,0:), CMU(:),  &
                CWT(:), DTAUCP(:), GC(:,:,:),   &
                KK(:,:), TAUCPR(0:)
      real(dk), intent(inout) :: WK(:)
      real(dk), intent(out)   :: CBAND(:,:)
!     ..
!     .. Local Scalars ..

      integer   :: IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      real(dk)  :: EXPA, SUM
!     ..

      CBAND = rZERO

      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
!  ** Use continuity conditions of Eq. STWJ(17)
!     to form coefficient matrix in STWJ(20);
!     employ scaling transformation STWJ(22)
LAYER_LOOP: &
      DO LC = 1, NCUT
         DO IQ = 1, NN
            WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
         ENDDO
         JCOL  = 0
         DO IQ = 1, NN
            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL
            DO JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
               IROW  = IROW + 1
            ENDDO
            JCOL  = JCOL + 1
         ENDDO

         DO IQ = NN + 1, NSTR
            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL
            DO JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )*WK(NSTR + 1 - IQ)
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
               IROW  = IROW + 1
            ENDDO
            JCOL  = JCOL + 1
         ENDDO
      ENDDO LAYER_LOOP

!  ** Use top boundary condition of STWJ(20a) for first layer
      JCOL  = 0

      DO IQ = 1, NN
         EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
         IROW  = NSHIFT - JCOL + NN
         DO JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
            IROW  = IROW + 1
         ENDDO
         JCOL  = JCOL + 1
      ENDDO


      DO IQ = NN + 1, NSTR
         IROW  = NSHIFT - JCOL + NN
         DO JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
            IROW  = IROW + 1
         ENDDO
         JCOL  = JCOL + 1
      ENDDO

!  ** Use bottom boundary condition of STWJ(20c) for last layer
      NNCOL = NCOL - NSTR
      JCOL  = 0

      DO IQ = 1, NN
         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         DO JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

!  ** No azimuthal-dependent intensity if Lambert surface;
!     no intensity component if truncated bottom layer
               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )
            ELSE
               SUM  = rZERO
               DO K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )  &
                             * GC( NN + 1 - K, IQ, NCUT )
               ENDDO
               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) - (1. + DELM0)*SUM
            END IF
            IROW  = IROW + 1
         ENDDO
         JCOL  = JCOL + 1
      ENDDO

      DO IQ = NN + 1, NSTR
         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         EXPA   = WK( NSTR + 1 - IQ )

         DO JQ = NN + 1, NSTR
            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN
               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA
            ELSE
               SUM  = rZERO
               DO K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )  &
                               * GC( NN + 1 - K, IQ, NCUT )
               ENDDO
               CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) - (1. + DELM0)*SUM )*EXPA
            END IF
            IROW  = IROW + 1
         ENDDO
         JCOL  = JCOL + 1
      ENDDO

      END SUBROUTINE SETMTX

      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MAZIM,     &
                         NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC,  &
                         AAD, EVECCD, EVALD, WKD )

!         Solves eigenvalue/vector problem necessary to construct
!         homogeneous part of discrete ordinate solution; STWJ(8b)
!         ** NOTE ** Eigenvalue problem is degenerate when single
!                    scattering albedo = 1;  present way of doing it
!                    seems numerically more stable than alternative
!                    methods that we tried

!   I N P U T     V A R I A B L E S:

!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                    (including factors 2l+1 and single-scatter albedo)
!       CMU    :  Computational polar angle cosines
!       CWT    :  Weights for quadrature over polar angle cosine
!       MAZIM  :  Order of azimuthal component
!       NN     :  Half the total number of streams
!       YLMC   :  Normalized associated Legendre polynomial
!                    at the quadrature angles CMU
!       (remainder are DISORT input variables)

!   O U T P U T    V A R I A B L E S:

!       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
!       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
!                    but then square roots taken
!       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
!                    from ASYMTX ( column j corresponds to EVAL(j) )
!                    but then  (G+) + (G-)  is calculated from SS(10),
!                    G+  and  G-  are separated, and  G+  is stacked on
!                    top of  G-  to form NSTR eigenvectors of SS(7)
!       GC     :  Permanent storage for all NSTR eigenvectors, but
!                    in an order corresponding to KK
!       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
!                    but re-ordered with negative values first ( square
!                    roots of EVAL taken and negatives added )

!   I N T E R N A L   V A R I A B L E S:

!       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
!                    eigenvalue problem
!       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
!                    problem: (alfa+beta)*(alfa-beta)
!       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
!       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
!       WKD     :  Scratch array required by ASYMTX

!   Called by- DISORT, ALBTRN
!   Calls- ASYMTX, ERRMSG
! +-------------------------------------------------------------------+


!     .. Scalar Arguments ..

      integer, intent(in) ::  MAZIM, NN, NSTR
!     ..
!     .. Array Arguments ..

      real(dk), intent(out) ::  EVAL(:), KK(:)
      real(dk), intent(out) ::  CC(:,:), EVECC(:,:), GC(:,:)
      real(dk), intent(out) ::  AMB(:,:), APB(:,:), ARRAY(:,:)
      real(dk), intent(in)  ::  CMU(:), CWT(:), GL(0:), YLMC(0:,:)
      real(dk), intent(inout)  :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)

!     ..
!     .. Local Scalars ..

      integer :: IER, IQ, JQ, KQ, L
      real(dk)    :: ALPHA, BETA, GPMIGM, GPPLGM, SUM

!                             ** Calculate quantities in Eqs. SS(5-6)
      DO IQ = 1, NN
         DO JQ = 1, NSTR
            SUM  = rZERO
            DO L = MAZIM, NSTR - 1
               SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
            ENDDO
            CC( IQ, JQ ) = ONEHALF*SUM*CWT( JQ )
         ENDDO

         DO JQ = 1, NN
!                             ** Fill remainder of array using symmetry
!                                relations  C(-mui,muj) = C(mui,-muj)
!                                and        C(-mui,-muj) = C(mui,muj)
            CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
            CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )

!                                       ** Get factors of coeff. matrix
!                                          of reduced eigenvalue problem

            ALPHA  = CC( IQ, JQ ) / CMU( IQ )
            BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
            AMB( IQ, JQ ) = ALPHA - BETA
            APB( IQ, JQ ) = ALPHA + BETA
         ENDDO

         AMB( IQ, IQ ) = AMB( IQ, IQ ) - rONE / CMU( IQ )
         APB( IQ, IQ ) = APB( IQ, IQ ) - rONE / CMU( IQ )
      ENDDO
!                      ** Finish calculation of coefficient matrix of
!                         reduced eigenvalue problem:  get matrix
!                         product (alfa+beta)*(alfa-beta); SS(12)
      DO IQ = 1, NN
         DO JQ = 1, NN
            SUM  = rZERO
            DO KQ = 1, NN
               SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
            ENDDO
            ARRAY( IQ, JQ ) = SUM
         ENDDO
      ENDDO
!                      ** Find (real) eigenvalues and eigenvectors

      CALL ASYMTX( ARRAY, EVECC, EVAL, IER, WKD, AAD, EVECCD, EVALD )

      IF( IER.GT.0 ) THEN
         WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',  &
            IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'
         CALL ERRMSG( 'ASYMTX--convergence problems',.True.)
      END IF

      DO IQ = 1, NN
         EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
         KK( IQ + NN ) = EVAL( IQ )
!                                      ** Add negative eigenvalue
         KK( NN + 1 - IQ ) = -EVAL( IQ )
      ENDDO

!                          ** Find eigenvectors (G+) + (G-) from SS(10)
!                             and store temporarily in APB array
      DO JQ = 1, NN
         DO IQ = 1, NN
            SUM  = rZERO
            DO KQ = 1, NN
               SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
            ENDDO
            APB( IQ, JQ ) = SUM / EVAL( JQ )
         ENDDO
      ENDDO

      DO JQ = 1, NN
         DO IQ = 1, NN
            GPPLGM = APB( IQ, JQ )
            GPMIGM = EVECC( IQ, JQ )
!                                ** Recover eigenvectors G+,G- from
!                                   their sum and difference; stack them
!                                   to get eigenvectors of full system
!                                   SS(7) (JQ = eigenvector number)
            EVECC( IQ,      JQ ) = ONEHALF*( GPPLGM + GPMIGM )
            EVECC( IQ + NN, JQ ) = ONEHALF*( GPPLGM - GPMIGM )
!                                ** Eigenvectors corresponding to
!                                   negative eigenvalues (corresp. to
!                                   reversing sign of 'k' in SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = ONEHALF * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = ONEHALF * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )
         ENDDO
      ENDDO

      END SUBROUTINE SOLEIG

      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,  &
                         FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM, &
                         NCOL, NCUT, NN, NSTR,                          &
                         PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

      use tuvx_linear_algebra_linpack, only : linear_algebra_linpack_t

!        Construct right-hand side vector B for general boundary
!        conditions STWJ(17) and solve system of equations obtained
!        from the boundary conditions and the continuity-of-
!        intensity-at-layer-interface equations.
!        Thermal emission contributes only in azimuthal independence.

!     I N P U T      V A R I A B L E S:

!       BDR      :  Surface bidirectional reflectivity
!       BEM      :  Surface bidirectional emissivity
!       BPLANK   :  Bottom boundary thermal emission
!       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
!                   scaled by Eq. SC(12); in banded form required
!                   by LINPACK solution routines
!       CMU      :  Abscissae for Gauss quadrature over angle cosine
!       CWT      :  Weights for Gauss quadrature over angle cosine
!       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
!       LYRCUT   :  Logical flag for truncation of comput. layer
!       MAZIM    :  Order of azimuthal component
!       ncol     :  Counts of columns in CBAND
!       NN       :  Order of double-Gauss quadrature (NSTR/2)
!       NCUT     :  Total number of computational layers considered
!       TPLANK   :  Top boundary thermal emission
!       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
!       ZZ       :  Beam source vectors in Eq. SS(19)
!       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
!       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
!       (remainder are DISORT input variables)

!   O U T P U T     V A R I A B L E S:

!       B        :  Right-hand side vector of Eq. SC(5) going into
!                   SGBSL; returns as solution vector of Eq. SC(12),
!                   constants of integration without exponential term
!
!      LL        :  Permanent storage for B, but re-ordered

!   I N T E R N A L    V A R I A B L E S:

!       IPVT     :  Integer vector of pivot indices
!       IT       :  Pointer for position in  B
!       NCD      :  Number of diagonals below or above main diagonal
!       RCOND    :  Indicator of singularity for CBAND
!       Z        :  Scratch array required by SGBCO

!   Called by- DISORT
!   Calls- SGBCO, ERRMSG, SGBSL
! +-------------------------------------------------------------------+


!     .. Scalar Arguments ..

      logical, intent(in) :: LAMBER, LYRCUT
      integer, intent(in) :: MAZIM, NCOL, NCUT, NN, NSTR
      real(dk),    intent(in) :: BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0
!     ..
!     .. Array Arguments ..

      integer, intent(inout) :: IPVT(:)
      real(dk), intent(out)      :: B(:), LL(:,:)
      real(dk), intent(inout)    :: Z(:)
      real(dk), intent(inout)    :: CBAND(:,:)
      real(dk), intent(in)       :: BDR(:,0:), BEM(:)
      real(dk), intent(in)       :: CMU(:), CWT(:)
      real(dk), intent(in)       :: EXPBEA(0:), TAUCPR(0:)
      real(dk), intent(in)       :: ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)
!     ..
!     .. Local Scalars ..

      integer ::   IPNT, IQ, IT, JQ, LC, NCD
      real(dk)    :: RCOND, SUM
      TYPE(linear_algebra_linpack_t) :: linpack

      B = rZERO
!                              ** Construct B,  STWJ(20a,c) for
!                                 parallel beam + bottom reflection +
!                                 thermal emission at top and/or bottom

      IF( MAZIM.GT.0 .AND. FBEAM.GT.0.0 ) THEN
!                                         ** Azimuth-dependent case
!                                            (never called if FBEAM = 0)
         IF( LYRCUT .OR. LAMBER ) THEN
!               ** No azimuthal-dependent intensity for Lambert surface;
!                  no intensity component for truncated bottom layer

            DO IQ = 1, NN
!                                                  ** Top boundary
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )
!                                                  ** Bottom boundary
               B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )
            ENDDO
         ELSE
            DO IQ = 1, NN
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )
               SUM  = rZERO
               DO JQ = 1, NN
                  SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ ) &
                               * ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
               ENDDO
               B( NCOL - NN + IQ ) = SUM
               IF( FBEAM.GT.0.0 ) B( NCOL - NN + IQ ) = SUM                &
                 + (BDR(IQ,0)*UMU0*FBEAM / PI - ZZ(IQ + NN,NCUT))*EXPBEA(NCUT)
            ENDDO
         END IF
!                             ** Continuity condition for layer
!                                interfaces of Eq. STWJ(20b)
         IT   = NN
         DO LC = 1, NCUT - 1

            DO IQ = 1, NSTR
               IT   = IT + 1
               B( IT ) = ( ZZ( IQ, LC+1 ) - ZZ( IQ, LC ) )*EXPBEA( LC )
            ENDDO
         ENDDO
      ELSE
!                                   ** Azimuth-independent case
         IF( FBEAM .EQ. rZERO ) THEN
            !  ** Top boundary
            DO IQ = 1, NN
               B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
            ENDDO

            IF( LYRCUT ) THEN
               !  ** No intensity component for truncated bottom layer
               DO IQ = 1, NN
                  B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) &
                                        - ZPLK1( IQ + NN, NCUT )*TAUCPR(NCUT)
               ENDDO
            ELSE
               DO IQ = 1, NN
                  SUM  = rZERO
                  DO JQ = 1, NN
                     SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ ) &
                                * (ZPLK0( NN + 1 - JQ,NCUT ) &
                                   + ZPLK1( NN + 1 - JQ,NCUT )*TAUCPR( NCUT ))
                  ENDDO
                  B( NCOL - NN + IQ ) = rTWO*SUM + BEM( IQ )*BPLANK  &
                                      - ZPLK0( IQ + NN, NCUT )     &
                                      - ZPLK1( IQ + NN, NCUT )*TAUCPR(NCUT)
               ENDDO
            END IF
            !  ** Continuity condition for layer interfaces, STWJ(20b)
            IT   = NN
            DO LC = 1, NCUT - 1
               DO IQ = 1, NSTR
                  IT   = IT + 1
                  B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) &
                          + ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )*TAUCPR(LC)
               ENDDO
            ENDDO
         ELSE
            DO IQ = 1, NN
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 ) &
                         - ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
            ENDDO

            IF( LYRCUT ) THEN
               DO IQ = 1, NN
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT)  &
                                  - ZPLK0(IQ+NN, NCUT) - ZPLK1(IQ+NN, NCUT)*TAUCPR(NCUT)
               ENDDO
            ELSE
               DO IQ = 1, NN
                  SUM  = rZERO
                  DO JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)       &
                                * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT)  &
                                  + ZPLK0(NN+1-JQ, NCUT)              &
                                  + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
                  ENDDO
                  B(NCOL-NN+IQ) = rTWO*SUM + (BDR(IQ,0) * UMU0*FBEAM/PI   &
                                      - ZZ(IQ+NN, NCUT)) * EXPBEA(NCUT) &
                                  + BEM(IQ) * BPLANK                     &
                                  - ZPLK0(IQ+NN, NCUT) - ZPLK1(IQ+NN, NCUT)*TAUCPR(NCUT)
               ENDDO
            END IF

            IT   = NN

            DO LC = 1, NCUT - 1
               DO IQ = 1, NSTR
                  IT   = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)  &
                          + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC)           &
                          + (ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC))*TAUCPR(LC)
               ENDDO
            ENDDO
         END IF
      END IF

!  ** Form L-U (lower/upper triangular) decomposition
!     of band matrix CBAND and test if it is nearly
!     singular (note: CBAND is destroyed)
!     (CBAND is in LINPACK packed format)

      RCOND  = rZERO
      NCD    = 3*NN - 1

      CALL linpack%SGBCO( CBAND, NCOL, NCD, NCD, IPVT, RCOND, Z )

      IF( rONE + RCOND == rONE ) &
        CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)

!     ** Solve linear system with coeff matrix CBAND
!     and R.H. side(s) B after CBAND has been L-U
!     decomposed.  Solution is returned in B.

      CALL linpack%SGBSL( CBAND, NCOL, NCD, NCD, IPVT, B, 0 )

!     ** Zero CBAND (it may contain 'foreign'
!        elements upon returning from LINPACK);
!        necessary to prevent errors

      CBAND = rZERO

      DO LC = 1, NCUT
         IPNT  = LC*NSTR - NN
         DO IQ = 1, NN
            LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
            LL( IQ + NN,     LC ) = B( IQ + IPNT )
         ENDDO
      ENDDO

      END SUBROUTINE SOLVE0

      SUBROUTINE SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER, MAZIM,  &
                         NN, NUMU, NSTR, ONLYFL, UMU,                &
                         USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM, RMU, &
                         solver_constants )

!       Specifies user's surface bidirectional properties, STWJ(21)

!   I N P U T     V A R I A B L E S:

!       DELM0  :  Kronecker delta, delta-sub-m0
!       HLPR   :  Legendre moments of surface bidirectional reflectivity
!                    (with 2K+1 factor included)
!       MAZIM  :  Order of azimuthal component
!       NN     :  Order of double-Gauss quadrature (NSTR/2)
!       YLM0   :  Normalized associated Legendre polynomial
!                 at the beam angle
!       YLMC   :  Normalized associated Legendre polynomials
!                 at the quadrature angles
!       YLMU   :  Normalized associated Legendre polynomials
!                 at the user angles
!       (remainder are DISORT input variables)

!    O U T P U T     V A R I A B L E S:

!       BDR :  Surface bidirectional reflectivity (computational angles)
!       RMU :  Surface bidirectional reflectivity (user angles)
!       BEM :  Surface directional emissivity (computational angles)
!       EMU :  Surface directional emissivity (user angles)

!    I N T E R N A L     V A R I A B L E S:

!       DREF      Directional reflectivity
!       NMUG   :  Number of angle cosine quadrature points on (0,1) for
!                   integrating bidirectional reflectivity to get
!                   directional emissivity (it is necessary to use a
!                   quadrature set distinct from the computational
!                   angles, because the computational angles may not be
!                   dense enough--NSTR may be too small--to give an
!                   accurate approximation for the integration).
!       GMU    :  The NMUG angle cosine quadrature points on (0,1)
!       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
!       YLMG   :  Normalized associated Legendre polynomials
!                   at the NMUG quadrature angles

!   Called by- DISORT
!   Calls- QGAUSN, LEPOLY, ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      logical, intent(in) :: LAMBER, ONLYFL, USRANG
      integer, intent(in) :: MAZIM, NN, NSTR, NUMU
      real(dk), intent(in)    :: ALBEDO, DELM0, FBEAM
!     ..
!     .. Array Arguments ..

      real(dk), intent(in)  ::  HLPR(0:), UMU(:)
      real(dk), intent(in)  ::  YLM0(0:), YLMC(0:,:), YLMU(0:,:)
      real(dk), intent(out) ::  BDR(:,0:), BEM(:), EMU(:), RMU(:,0:)
!     ..
!     .. Type Arguments ..
      type(solver_constants_t), intent(in) :: solver_constants
!     ..
!     .. Local Scalars ..

      integer  :: IQ, IU, JG, JQ, K
      real(dk) :: DREF, SGN, SUM
!     ..

      BDR = rZERO
      BEM = rZERO

      IF( LAMBER .AND. MAZIM == 0 ) THEN

         DO IQ = 1, NN
            BEM(IQ) = rONE - ALBEDO
            DO JQ = 0, NN
               BDR(IQ,JQ) = ALBEDO
            ENDDO
         ENDDO

      ELSE IF( .NOT. LAMBER ) THEN
!                                  ** Compute surface bidirectional
!                                     properties at computational angles
         DO IQ = 1, NN
            DO JQ = 1, NN
               SUM  = rZERO
               DO K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLMC(K,JQ+NN)
               ENDDO
               BDR(IQ,JQ) = (rTWO - DELM0)*SUM
            ENDDO

            IF( FBEAM > rZERO ) THEN
               SUM  = rZERO
               DO K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLM0(K)
               ENDDO
               BDR(IQ,0) = (rTWO - DELM0)*SUM
            END IF
         ENDDO

         IF( MAZIM == 0 ) THEN
            IF( NSTR > MAXSTR ) &
              CALL ERRMSG('SURFAC--parameter MAXSTR too small',.True.)
!                              ** Integrate bidirectional reflectivity
!                                 at reflection polar angles CMU and
!                                 incident angles solver_constants%GMU to get
!                                 directional emissivity at
!                                 computational angles CMU.
            DO IQ = 1, NN
               DREF  = rZERO
               DO JG = 1, NMUG
                  SUM  = rZERO
                  DO K = 0, NSTR - 1
                     SUM  = SUM + HLPR(K)*YLMC(K,IQ)*solver_constants%YLMG(K,JG)
                  ENDDO
                  DREF  = DREF + rTWO*solver_constants%GWT(JG)*solver_constants%GMU(JG)*SUM
               ENDDO
               BEM(IQ) = rONE - DREF
            ENDDO
         END IF
      END IF
!                                       ** Compute surface bidirectional
!                                          properties at user angles

      IF( .NOT. ONLYFL .AND. USRANG ) THEN
         EMU = rZERO
         RMU = rZERO
         DO IU = 1, NUMU
            IF( UMU( IU ) > rZERO ) THEN
               IF( LAMBER .AND. MAZIM == 0 ) THEN
                  DO IQ = 0, NN
                     RMU( IU, IQ ) = ALBEDO
                  ENDDO
                  EMU( IU ) = rONE - ALBEDO
               ELSE IF( .NOT.LAMBER ) THEN
                  DO IQ = 1, NN
                     SUM  = rZERO
                     DO K = MAZIM, NSTR - 1
                        SUM = SUM + HLPR(K)*YLMU(K,IU)*YLMC(K,IQ + NN)
                     ENDDO
                     RMU(IU,IQ) = (rTWO - DELM0)*SUM
                  ENDDO


                  IF( FBEAM > rZERO ) THEN
                     SUM  = rZERO
                     DO K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*YLM0( K )
                     ENDDO
                     RMU(IU,0) = (rTWO - DELM0)*SUM
                  END IF

                  IF( MAZIM == 0 ) THEN
!                               ** Integrate bidirectional reflectivity
!                                  at reflection angles UMU and
!                                  incident angles solver_constants%GMU to get
!                                  directional emissivity at
!                                  user angles UMU.
                     DREF  = rZERO
                     DO JG = 1, NMUG
                        SUM  = rZERO
                        DO K = 0, NSTR - 1
                           SUM = SUM + HLPR(K)*YLMU(K,IU)*solver_constants%YLMG(K,JG)
                        ENDDO
                        DREF  = DREF + rTWO*solver_constants%GWT( JG )*solver_constants%GMU( JG )*SUM
                     ENDDO
                     EMU( IU ) = rONE - DREF
                  END IF
               END IF
            END IF
         ENDDO
      END IF

      END SUBROUTINE SURFAC

      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL,   &
           MAZIM, NN, NSTR, YLM0, YLMC, CC,               &
           EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD,       &
           WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ, &
           OPRIM, LC, mu2, glsave, dgl, solver_constants)

!bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially 
!bm  unstable solution, the calculation is repeated with a slightly 
!bm  changed single scattering albedo; this process is iterates 
!bm  until a stable solution is found; as stable solutions may be 
!bm  reached either by increasing or by decreasing the single 
!bm  scattering albedo, both directions are explored ('upward' and
!bm  'downward' iteration); the solution which required the smaller 
!bm  change in the single scattering albedo is finally returned 
!bm  by SOLVEC.

!gy added glsave and dgl to call to allow adjustable dimensioning

!     .. Scalar Arguments ..

      integer, intent(in)  :: MAZIM, NN, NSTR, LC
      REAL(dk), intent(in) :: DELM0, FBEAM, PI, UMU0, OPRIM
      REAL(dk), intent(in) :: mu2

!     ..
!     .. Array Arguments ..

      integer, intent(inout)  ::   IPVT(:)
      
      REAL(dk), intent(inout) :: WK(:)
      REAL(dk), intent(in)    :: CMU(:), CWT(:), YLM0(0:), YLMC(0:,:)

      REAL(dk), intent(inout) :: AMB(:,:), APB(:,:), ARRAY(:,:)
      REAL(dk), intent(inout) :: EVAL(:), EVECC(:,:), KK(:), GC(:,:)
      REAL(dk), intent(inout) :: GLSAVE(0:), DGL(0:), CC(:,:), ZJ(:), ZZ(:)
      REAL(dk), intent(out)   :: GL(0:)

      integer  :: K
      REAL(dk) :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)

!     .. Type Arguments ..

      type(solver_constants_t), intent(in) :: solver_constants

!     ..
!bm   Variables for instability fix
      
      integer  :: UAGAIN, DAGAIN
      REAL(dk) :: MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR
      
      logical ::  DONE, NOUP, NODN, DEBUG, INSTAB
      logical ::  small_ssa, large_ssa
      
!bm   reset parameters

      DONE = .FALSE.
      NOUP = .FALSE.
      NODN = .FALSE.


!bm   flag for printing debugging output      
!      DEBUG  = .TRUE.
      DEBUG  = .FALSE.

!bm   instability parameter; the solution is considered 
!bm   unstable, if the RCOND reported by SGECO is smaller 
!bm   than MINRCOND
      MINRCOND = 1.0e-7_dk

!bm   if an instability is detected, the single scattering albedo
!bm   is iterated downwards in steps of DADD and upwards in steps 
!bm   of UADD; in practice, MINRCOND and -MINRCOND should 
!bm   be reasonable choices for these parameters
      DADD    = -MINRCOND
      UADD    = MINRCOND

      UAGAIN = 0
      DAGAIN = 0
      ADD   = DADD
      

!bm   save array GL( ) because it will be 
!bm   changed if an iteration should be neccessary
      GLSAVE(mazim:nstr-1) =  GL(mazim:nstr-1)
      
      SSA = OPRIM

!bm   in case of an instability reported by UPBEAM (INSTAB)
!bm   the single scattering albedo will be changed by a small 
!bm   amount (ADD); this is indicated by DAGAIN or UAGAIN 
!bm   being larger than 0; a change in the single scattering 
!bm   albedo is equivalent to scaling the array GL( )

      small_ssa = .false. ; large_ssa = .false.
CONV_LOOP: &
      DO WHILE( .true. )
         IF ( DAGAIN > 0 .OR. UAGAIN > 0)  THEN
            FACTOR = (SSA + ADD) / SSA
            GL(mazim:nstr-1) =  GL(mazim:nstr-1) * FACTOR
            SSA = SSA + ADD
         
!bm   if the single scattering albedo is now smaller than 0
!bm   the downward iteration is stopped and upward iteration 
!bm   is forced instead
            IF( SSA < solver_constants%DITHER) THEN
               NODN = .TRUE.
               DAGAIN = -1
               small_ssa = .true.
!bm   if the single scattering albedo is now larger than its maximum 
!bm   allowed value (1.0 - solver_constants%DITHER), the upward iteration is 
!bm   stopped and downward iteration is forced instead
            ELSEIF( SSA > (rONE - solver_constants%DITHER)) THEN
               NOUP = .TRUE.
               UAGAIN = -1
               large_ssa = .true.
            ENDIF
         ENDIF

!     ** Solve eigenfunction problem in Eq. STWJ(8B);
!        return eigenvalues and eigenvectors

         DO WHILE( .true. )
            IF( .not. (small_ssa .or. large_ssa) ) THEN
               CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL,          &
                            MAZIM, NN, NSTR, YLMC, CC, EVECC, EVAL, &
                            KK, GC, AAD, EVECCD, EVALD, WKD )

!     ** Calculate particular solutions of
!        q.SS(18) for incident beam source
               IF ( FBEAM > rZERO ) THEN
                  CALL  UPBEAM( mu2,                         &
                       ARRAY, CC, CMU, DELM0, FBEAM, GL,     &
                       IPVT, MAZIM, NN, NSTR, PI, UMU0, WK,  &
                       YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)
               ENDIF
      
!     ** Calculate particular solutions of
!        Eq. SS(15) for thermal emission source
!        (not available in psndo.f)
      
!bm   finished if the result is stable on the first try
               IF ( .NOT. INSTAB .AND. UAGAIN == 0 .AND. DAGAIN == 0 ) THEN
                  exit CONV_LOOP
               ENDIF

               IF( INSTAB )  THEN
!bm   downward iteration
                  IF( UAGAIN == 0 )  THEN
                     DAGAIN = DAGAIN + 1
!bm   upward iteration
                  ELSEIF( UAGAIN > 0 )  THEN
                     UAGAIN = UAGAIN + 1
                  ENDIF
                  cycle CONV_LOOP
               ENDIF
            ENDIF
      
!bm   ( DAGAIN .NE. 0 ) at this place means that the downward
!bm   iteration is finished 
            IF (.not. large_ssa ) THEN
               IF( small_ssa ) small_ssa = .false.
               IF (DAGAIN /= 0 .AND. UAGAIN == 0) THEN
!bm   save downward iteration data for later use and 
!bm   restore original input data
                  DGL(mazim:nstr-1) =  GL(mazim:nstr-1 )
                  GL(mazim:nstr-1 ) =  GLSAVE(mazim:nstr-1)
                  DSSA = SSA
                  SSA = OPRIM
!bm   start upward iteration
                  ADD = UADD
                  UAGAIN = UAGAIN + 1
                  cycle CONV_LOOP
               ENDIF
            ELSE
               large_ssa = .false.
            ENDIF

!bm   both iterations finished
            IF (DONE) THEN
               GL(mazim:nstr-1) = GLSAVE(mazim:nstr-1)
               exit CONV_LOOP
            ENDIF

!bm  if neither upward nor downward iteration converged, the 
!bm  original conditions are restored and SOLEIG/UPBEAM 
!bm  is called for the last time 
         
            IF (NOUP .AND. NODN) THEN
               GL(mazim:nstr-1) =  GLSAVE(mazim:nstr-1)
               SSA = OPRIM
               IF (DEBUG) THEN
                  write (*,*) '! *** Neither upward nor downward iteration'
                  write (*,*) '! *** converged; using original result.'
               ENDIF
               DONE = .TRUE.
               cycle CONV_LOOP
            ENDIF

!bm  if upward iteration did not converge, the stable downward conditions
!bm  are restored and SOLEIG/UPBEAM is called for the last time
            IF (NOUP) THEN
               GL(mazim:nstr-1) =  DGL(mazim:nstr-1)
               SSA = DSSA
               IF (DEBUG) THEN
                  write (*,*) '! *** The upward iteration did not converge.'
                  write (*,*) '! *** Had to iterate ', DAGAIN,' times in layer LC =', LC,';'
                  write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','
                  write (*,*) '! *** by a factor of ', SSA/OPRIM
               ENDIF
               DONE = .TRUE.
               cycle CONV_LOOP
            ENDIF

!bm  if downward iteration did not converge, we are done 
!bm  (the result of the upward iteration will be used)
            IF (NODN) THEN
               IF (DEBUG) THEN
                  write (*,*) '! *** The downward iteration did not converge.'
                  write (*,*) '! *** Had to iterate ', UAGAIN,' times in layer LC =', LC,';'
                  write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','
                  write (*,*) '! *** by a factor of ', SSA/OPRIM
               ENDIF
         
               GL(mazim:nstr-1) =  GLSAVE(mazim:nstr-1)
               exit CONV_LOOP
            ENDIF
      
!bm   if both iterations converged, and if the upward iteration 
!bm   required more steps than the downward iteration, the stable 
!bm   downward conditions are restored and SOLEIG/UPBEAM is 
!bm   called for the last time 
         
            IF (UAGAIN > DAGAIN) THEN
               GL(mazim:nstr-1) =  DGL(mazim:nstr-1)
               SSA = DSSA
               IF (DEBUG) THEN
                  write (*,*) '! *** Both iterations converged; using downward.'
                  write (*,*) '! *** Had to iterate ', DAGAIN,' times in layer LC =', LC,';'
                  write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','
                  write (*,*) '! *** by a factor of ', SSA/OPRIM
               ENDIF

               DONE = .TRUE.
               cycle CONV_LOOP
            ELSE
               IF (DEBUG) THEN
                  write (*,*) '! *** Both iterations converged; using upward.'
                  write (*,*) '! *** Had to iterate ', UAGAIN,' times in layer LC =', LC,';'
                  write (*,*) '! *** changed SSA from ',OPRIM, ' to ', SSA,','
                  write (*,*) '! *** by a factor of ', SSA/OPRIM
               ENDIF

               GL(mazim:nstr-1) =  GLSAVE(mazim:nstr-1)
               exit CONV_LOOP
            ENDIF
         ENDDO
      ENDDO CONV_LOOP
      
      END SUBROUTINE SOLVEC

      SUBROUTINE UPBEAM( mu2,                                            &
                         ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,  &
                         NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,         &
                         ZZ, MINRCOND, INSTAB )

      use tuvx_linear_algebra_linpack, only : linear_algebra_linpack_t

!         Finds the incident-beam particular solution of SS(18)

!   I N P U T    V A R I A B L E S:

!       CC     :  C-sub-ij in Eq. SS(5)
!       CMU    :  Abscissae for Gauss quadrature over angle cosine
!       DELM0  :  Kronecker delta, delta-sub-m0
!       GL     :  Delta-M scaled Legendre coefficients of phase function
!                    (including factors 2L+1 and single-scatter albedo)
!       MAZIM  :  Order of azimuthal component
!       YLM0   :  Normalized associated Legendre polynomial
!                    at the beam angle
!       YLMC   :  Normalized associated Legendre polynomial
!                    at the quadrature angles
!       (remainder are DISORT input variables)

!   O U T P U T    V A R I A B L E S:

!       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
!                 solution vector Z-sub-zero after solving that system

!       ZZ     :  Permanent storage for ZJ, but re-ordered

!   I N T E R N A L    V A R I A B L E S:

!       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
!       IPVT   :  Integer vector of pivot indices required by LINPACK
!       WK     :  Scratch array required by LINPACK

!   Called by- DISORT
!   Calls- SGECO, ERRMSG, SGESL
! +-------------------------------------------------------------------+


!     .. Scalar Arguments ..

      integer, intent(in)  :: MAZIM, NN, NSTR
      logical, intent(out) :: INSTAB
      real(dk), intent(in)     :: MINRCOND
      real(dk), intent(in)     :: DELM0, FBEAM, PI, UMU0
      real(dk), intent(in)     :: mu2
!     ..
!     .. Array Arguments ..

      integer, intent(inout) :: IPVT(:)

      real(dk), intent(in) :: CC(:,:), CMU(:)
      real(dk), intent(in) :: GL(0:), YLM0(0:), YLMC(0:,:)
      real(dk), intent(inout) :: WK(:)
      real(dk), intent(out) :: ARRAY(:,:)
      real(dk), intent(out) :: ZJ(:), ZZ(:)
!     ..
!     .. Local Scalars ..

      integer   :: IQ, JOB, JQ, K
      real(dk)      :: RCOND, SUM

      TYPE(linear_algebra_linpack_t) :: linpack
!     ..
!     .. External Subroutines ..

!     EXTERNAL  ERRMSG, SGECO, SGESL
!     ..


      DO IQ = 1, NSTR
         DO JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
         ENDDO

         ARRAY( IQ, IQ ) = rONE + CMU( IQ ) / mu2 + ARRAY( IQ, IQ )

         SUM  = rZERO
         DO K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
         ENDDO

         ZJ( IQ ) = ( rTWO - DELM0 )*FBEAM*SUM / ( rFOUR*PI )
      ENDDO

!                  ** Find L-U (lower/upper triangular) decomposition
!                     of ARRAY and see if it is nearly singular
!                     (NOTE:  ARRAY is destroyed)
      RCOND  = rZERO

      CALL linpack%SGECO( ARRAY, NSTR, IPVT, RCOND, WK )

!bm      IF( 1.0 + RCOND.EQ.1.0 )
!bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
!bm
!bm   replaced original check of RCOND by the following:

      INSTAB = .FALSE.
      IF( ABS(RCOND) .LT. MINRCOND )  THEN
         INSTAB = .TRUE.
         RETURN
      ENDIF

!                ** Solve linear system with coeff matrix ARRAY
!                   (assumed already L-U decomposed) and R.H. side(s)
!                   ZJ;  return solution(s) in ZJ
      JOB  = 0

      CALL linpack%SGESL( ARRAY, NSTR, IPVT, ZJ, JOB )

      DO IQ = 1, NN
         ZZ( IQ + NN )     = ZJ( IQ )
         ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
      ENDDO

      END SUBROUTINE UPBEAM

      SUBROUTINE ZEROAL( EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,  &
                         CMU, CWT, PSI, WK, Z0, Z1, ZJ,          &
                         HLPR, YLM0,                             &
                         ARRAY, CC, EVECC,                       &
                         GL, YLMC, YLMU,                         &
                         KK, LL, ZZ, ZPLK0, ZPLK1,               &
                         GC, LAYRU, UTAUPR,                      &
                         GU, Z0U, Z1U, ZBEAM,                    &
                         EVAL, AMB, APB, IPVT, Z,                &
                         RFLDIR, RFLDN, FLUP, UAVG, DFDT,        &
                         TRNMED, U0U, UU )

!         rZERO ARRAYS

!   Called by- DISORT
! --------------------------------------------------------------------

!     .. Array Arguments ..

      integer, intent(out) ::   IPVT(:), LAYRU(:)
      real(dk), intent(out)    ::                                      &
                AMB(:,:), APB(:,:), ARRAY(:,:), CC(:,:),           &
                CMU(:), CWT(:), DFDT(:), EVAL(:), EVECC(:,:),      &
                EXPBEA(:), FLUP(:), FLYR(:), GC(:,:,:), GL(:,:),   &
                GU(:,:,:), HLPR(:), KK(:,:), LL(:,:), OPRIM(:),    &
                PSI(:), RFLDIR(:), RFLDN(:), TAUCPR(:),            &
                TRNMED(:), U0U(:,:), UAVG(:), UTAUPR(:),           &
                UU(:,:,:),                                         &
                WK(:), XR0(:), XR1(:), YLM0(:), YLMC(:,:),         &
                YLMU(:,:), Z(:), Z0(:), Z0U(:,:), Z1(:), Z1U(:,:), &
                ZBEAM(:,:), ZJ(:), ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)

!     ..

      EXPBEA = rZERO
      FLYR   = rZERO
      OPRIM  = rZERO
      TAUCPR = rZERO
      XR0    = rZERO
      XR1    = rZERO

         CMU = rZERO
         CWT = rZERO
         PSI = rZERO
         WK  = rZERO
         Z0  = rZERO
         Z1  = rZERO
         ZJ  = rZERO

         HLPR = rZERO
         YLM0 = rZERO

         ARRAY = rZERO
         CC    = rZERO
         EVECC = rZERO

         GL = rZERO

         YLMC = rZERO

         YLMU = rZERO

         KK    = rZERO
         LL    = rZERO
         ZZ    = rZERO
         ZPLK0 = rZERO
         ZPLK1 = rZERO

         GC = rZERO

         LAYRU  = 0
         UTAUPR = rZERO

         GU = rZERO

         Z0U   = rZERO
         Z1U   = rZERO
         ZBEAM = rZERO

         EVAL = rZERO

         AMB = rZERO
         APB = rZERO

         IPVT = 0
         Z    = rZERO

         RFLDIR = rZERO
         RFLDN  = rZERO
         FLUP   = rZERO
         UAVG   = rZERO
         DFDT   = rZERO

         TRNMED = rZERO

         U0U = rZERO

         UU = rZERO

      END SUBROUTINE ZEROAL

      real(dk) FUNCTION DREF( MU, HL, NSTR, solver_constants )

!        Exact flux albedo for given angle of incidence, given
!        a bidirectional reflectivity characterized by its
!        Legendre coefficients ( NOTE** these will only agree
!        with bottom-boundary albedos calculated by DISORT in
!        the limit as number of streams go to infinity, because
!        DISORT evaluates the integral 'CL' only approximately,
!        by quadrature, while this routine calculates it exactly.)

!  INPUT :   MU     Cosine of incidence angle
!            HL     Legendre coefficients of bidirectional reflectivity
!          NSTR     Number of elements of HL to consider

!  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

!       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
!                   (vanishes for  L = 3, 5, 7, ... )
!       PL      P-sub-L
!       PLM1    P-sub-(L-1)
!       PLM2    P-sub-(L-2)

!   Called by- CHEKIN
!   Calls- ERRMSG
! +-------------------------------------------------------------------+

!     .. Scalar Arguments ..

      integer, intent(in)  :: NSTR
      real(dk), intent(in) :: MU
!     ..
!     .. Array Arguments ..

      real(dk), intent(in)   :: HL(0:NSTR)
!     ..
!     .. Type Arguments ..
      type(solver_constants_t), intent(in) :: solver_constants

!     ..
!     .. Local Scalars ..

      logical   :: PASS1
      integer   :: L
      real(dk)  :: PL, PLM1, PLM2
!     ..


      IF( NSTR < 2 .OR. ABS(MU) > rONE ) CALL ERRMSG( 'DREF--input argument error(s)',.True. )

      IF( NSTR > MAXTRM ) CALL ERRMSG( 'DREF--parameter MAXTRM too small',.True. )


      DREF  = HL(0) - rTWO*HL(1)*MU
      PLM2  = rONE
      PLM1  = - MU

      DO L = 2, NSTR - 1
!     ** Legendre polynomial recurrence
         PL = (real(2*L - 1,dk)*(-MU)*PLM1 - real(L-1,dk)*PLM2) / real(L,dk)
         IF( MOD( L,2 ) == 0 ) DREF = DREF + solver_constants%C(L)*HL(L)*PL
         PLM2  = PLM1
         PLM1  = PL
      ENDDO

      IF( DREF < rZERO .OR. DREF > rONE ) &
        CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )

      END FUNCTION DREF

      real(dk) FUNCTION RATIO( A, B )
! ---------------------------------------------------------------
!        Calculate ratio  A/B  with over- and under-flow protection
!        (thanks to Prof. Jeff Dozier for some suggestions here).
!        Since this routine takes two logs, it is no speed demon,
!        but it is invaluable for comparing results from two runs
!        of a program under development.
! ---------------------------------------------------------------

!     .. Scalar Arguments ..

      real(dk), intent(in) :: A, B
!     ..
!     .. Local Scalars ..

      logical  ::   PASS1
      real(dk) ::   ABSA, ABSB, LARGE, POWA, POWB, POWMAX, POWMIN, SMALL

!     ..
!     .. Intrinsic Functions ..

      INTRINSIC ABS, LOG10, SIGN
!     ..


      SMALL   = TINY( 1._dk )
      LARGE   = HUGE( 1._dk )
      POWMAX = LOG10( LARGE )
      POWMIN = LOG10( SMALL )

      IF( A == rZERO ) THEN
         IF( B == rZERO ) THEN
            RATIO  = rONE
         ELSE
            RATIO  = rZERO
         END IF
      ELSE IF( B == rZERO ) THEN
         RATIO  = SIGN( LARGE, A )
      ELSE
         ABSA   = ABS( A )
         ABSB   = ABS( B )
         POWA   = LOG10( ABSA )
         POWB   = LOG10( ABSB )

         IF( ABSA < SMALL .AND. ABSB < SMALL ) THEN
            RATIO  = rONE
         ELSE IF( POWA - POWB >= POWMAX ) THEN
            RATIO  = LARGE
         ELSE IF( POWA - POWB <= POWMIN ) THEN
            RATIO  = SMALL
         ELSE
            RATIO  = ABSA / ABSB
         END IF
!                      ** DONT use old trick of determining sign
!                      ** from A*B because A*B may (over/under)flow

         IF( ( A > rZERO .AND. B < rZERO ) .OR. &
             ( A < rZERO .AND. B > rZERO ) ) RATIO = -RATIO
      END IF

      END FUNCTION RATIO

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

!        Print out a warning or error message;  abort if error
!        after making symbolic dump (machine-specific)

      logical       FATAL, MsgLim, Cray
      character*(*) MESSAG
      integer       MaxMsg, NumMsg
      SAVE          MaxMsg, NumMsg, MsgLim
      DATA NumMsg / 0 /,  MaxMsg / 100 /,  MsgLim / .FALSE. /


      IF ( FATAL )  THEN
         WRITE ( *, '(//,2A,//)' )  ' ******* ERROR >>>>>>  ', MESSAG
         STOP
      END IF

      NumMsg = NumMsg + 1
      IF( MsgLim )  RETURN

      IF ( NumMsg.LE.MaxMsg )  THEN
         WRITE ( *, '(/,2A,/)' )  ' ******* WARNING >>>>>>  ', MESSAG
      ELSE
         WRITE ( *,99 )
         MsgLim = .True.
      ENDIF

      RETURN

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',  &
         'They will no longer be printed  <<<<<<<', // )

      END SUBROUTINE ErrMsg

      logical FUNCTION  WrtBad ( VarNam )

!          Write names of erroneous variables and return 'TRUE'

!      INPUT :   VarNam = Name of erroneous variable to be written
!                         ( character, any length )

      character*(*)  VarNam
      integer        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,'  in error  ****'
      IF ( NumMsg.EQ.MaxMsg ) CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      END FUNCTION WrtBad

      logical FUNCTION  WrtDim ( DimNam, MinVal )

!          Write name of too-small symbolic dimension and
!          the value it should be increased to;  return 'TRUE'

!      INPUT :  DimNam = Name of symbolic dimension which is too small
!                        ( character, any length )
!               Minval = Value to which that dimension should be
!                        increased (at least)

      character*(*)  DimNam
      integer        MinVal


      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,  &
                           '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      END FUNCTION WrtDim

      logical FUNCTION  TstBad( VarNam, RelErr )

!       Write name (VarNam) of variable failing self-test and its
!       percent error from the correct value;  return  'FALSE'.

      character*(*)  VarNam
      real(dk)           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )  &
             ' Output variable ', VarNam,' differed by ', 100.*RelErr,  &
             ' per cent from correct value.  Self-test failed.'

      END FUNCTION TstBad

end module tuvx_discrete_ordinate_util
