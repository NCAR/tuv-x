      module RTRANS

      IMPLICIT NONE

      private
      public :: rtlink, PSNDO

      INTEGER, PARAMETER :: MXCLY = 151
      INTEGER, PARAMETER :: MXULV = 151
      INTEGER, PARAMETER :: MXCMU = 32
      INTEGER, PARAMETER :: MXUMU = 32
      INTEGER, PARAMETER :: MXPHI = 3
      REAL, PARAMETER    :: rZERO = 0.0
      REAL, PARAMETER    :: rONE  = 1.0
* Discrete ordinate constants:
* For pseudo-spherical DISORT, PLANK, USRTAU and USRANG must be .FALSE.;
* ONLYFL must be .TRUE.; FBEAM = 1.; FISOT = 0.; IBCND = 0
      INTEGER, parameter :: NPHI  = 0
      INTEGER, parameter :: IBCND = 0
      REAL, parameter    :: ACCUR = 0.0001
      REAL, parameter    :: FBEAM = 1.
      REAL, parameter    :: FISOT = rZERO
      REAL, parameter    :: PHI0  = rZERO
      LOGICAL, parameter :: DELTAM = .true.
      LOGICAL, parameter :: LAMBER = .true.
      LOGICAL, parameter :: PLANK = .FALSE.
      LOGICAL, parameter :: USRANG = .FALSE.
      LOGICAL, parameter :: USRTAU = .FALSE.
      LOGICAL, parameter :: ONLYFL = .TRUE.
      LOGICAL, parameter :: PRNT(7) = .FALSE.

      REAL :: DITHER

      contains

      SUBROUTINE rtlink(nstr, nz, 
     $     albedo, zen,
     $     dsdh, nid,
     $     dtrl, 
     $     dto3, 
     $     dto2,
     $     dtso2,
     $     dtno2, 
     $     dtcld, omcld, gcld,
     $     dtaer, omaer, gaer,
     $     dtsnw, omsnw, gsnw,
     $     dt_any,om_any,g_any,
     $     edir, edn, eup, fdir, fdn, fup)

      use tuv_params, only : pi, largest,precis
      use abstract_radXfer, only : abstract_radXfer_t
      use delta_eddington, only : delta_eddington_t
      use disord,          only : disord_t

      IMPLICIT NONE

* input

      INTEGER, intent(in) :: nstr
      INTEGER, intent(in) :: nz
      INTEGER, intent(in) :: nid(0:)
      REAL, intent(in) :: albedo
      REAL, intent(in) :: zen
      REAL, intent(in) :: dsdh(0:,:)
      REAL, intent(in) :: dtrl(:)
      REAL, intent(in) :: dto3(:), dto2(:)
      REAL, intent(in) :: dtso2(:), dtno2(:)
      REAL, intent(in) :: dtcld(:), omcld(:), gcld(:)
      REAL, intent(in) :: dtaer(:), omaer(:), gaer(:)
      REAL, intent(in) :: dtsnw(:), omsnw(:), gsnw(:)
      REAL, intent(in) :: dt_any(:), om_any(:), g_any(:)

* output

      REAL, intent(out) :: edir(:), edn(:), eup(:)
      REAL, intent(out) :: fdir(:), fdn(:), fup(:)

* constants:

      REAL, parameter    :: dr = pi/180.
      REAL, parameter    :: fourPI = 4. * pi

* local:

      INTEGER :: i, ii
      REAL :: dtabs,dtsct,dscld,dsaer,dssnw,dsany,dacld,
     $        daaer,dasnw,daany
      REAL :: dt(nz-1), om(nz-1), g(nz-1)

* specific two ps2str

      REAL :: ediri(nz), edni(nz), eupi(nz)
      REAL :: fdiri(nz), fdni(nz), fupi(nz)

*  specific to psndo:

      INTEGER :: istr, iu
      REAL :: pmcld, pmray, pmaer, pmsnw, pmany
      REAL :: om1

      INTEGER ::  NLYR
      REAL    ::  UMU0
      REAL    ::  DTAUC(nz-1),
     $            PMOM(0:nstr,nz-1),
     $            SSALB(nz-1)

      REAL :: RFLDIR(nz), RFLDN(nz), FLUP(nz)
      REAL :: uavgso(nz), uavgup(nz), uavgdn(nz)

      TYPE(delta_eddington_t) :: delta_eddington_rt
      TYPE(disord_t)          :: disord_rt
      class(abstract_radXfer_t), allocatable :: radiative_xfer_obj

*_______________________________________________________________________

* initialize:

      fdir = 0.
      fup = 0.
      fdn = 0.
      edir = 0.
      eup = 0.
      edn = 0.


      NLYR = nz - 1

      if( nstr < 2 ) then
        allocate( delta_eddington_t :: radiative_xfer_obj )
      else
        allocate( disord_t :: radiative_xfer_obj )
      endif
      call radiative_xfer_obj%initialize( 
     $       nlyr, nstr, zen, nid, dsdh,
     $       dtrl, dto3, dto2, dtso2, dtno2,
     $       dtcld, omcld, gcld,
     $       dtaer, omaer, gaer,
     $       dtsnw, omsnw, gsnw,
     $       dt_any, om_any, g_any )

      layer_loop: DO i = 1, NLYR

         dscld = dtcld(i)*omcld(i)
         dacld = dtcld(i)*(rONE - omcld(i))

         dsaer = dtaer(i)*omaer(i)
         daaer = dtaer(i)*(rONE - omaer(i))

         dssnw = dtsnw(i)*omsnw(i)
         dasnw = dtsnw(i)*(rONE - omsnw(i))

         dsany = dt_any(i)*om_any(i)
         daany = dt_any(i)*(rONE - om_any(i))

         dtsct = dtrl(i) + dscld + dsaer + dssnw + dsany
         dtabs = dtso2(i) + dto2(i) + dto3 (i)
     $           + dtno2(i) + dacld + daaer + dasnw + daany

* put in a floor
         dtabs = max(dtabs,rONE/largest)
         dtsct = max(dtsct,rONE/largest)

* from bottum-up -> top-down

         ii = nz - i
         dt(ii) = dtsct + dtabs
         om(ii) = dtsct/(dtsct + dtabs)
         IF( dtsct == 1./largest ) then
           om(ii) = 1./largest
         ENDIF
         g(ii) = (gcld(i)*dscld + 
     $            gsnw(i)*dssnw +
     $            gaer(i)*dsaer +
     $            g_any(i)*dsany)/dtsct

         IF(nstr > 1) THEN

* DISORD parameters

           om1       = MIN( om(ii),rONE - PRECIS )
           ssalb(ii) = MAX( om1,PRECIS )
           dtauc(ii) = MAX( dt(ii),PRECIS )
*  phase function - assume Henyey-Greenstein for cloud and aerosol
*  and Rayleigh for molecular scattering
           pmom(0,ii) = rONE
           DO istr = 1, nstr
             pmcld = gcld(i)**ISTR
             pmaer = gaer(i)**ISTR
             pmsnw = gsnw(i)**ISTR
             pmany = g_any(i)**ISTR
             IF(istr == 2) THEN
                pmray = 0.1
             ELSE
                pmray = rZERO
             ENDIF
             pmom(istr,ii) = (pmcld*dscld + 
     $           pmaer*dsaer + 
     $           pmsnw*dssnw +
     $           pmany*dsany +
     $           pmray*dtrl(i)) / dtsct
           ENDDO
         ENDIF

      ENDDO layer_loop

* call rt routine:

      UMU0 = cos(zen*dr)
      IF( nstr < 2 ) THEN
         call delta_eddington_rt%calculate( 
     $        nlyr, nstr, albedo,
     $        fdir, fup, fdn, edir, eup, edn )

* output (top-down -> buttom-up)
         fdir(1:nz) = fdiri(nz:1:-1)
         fup(1:nz) = fupi(nz:1:-1)
         fdn(1:nz) = fdni(nz:1:-1)
         edir(1:nz) = ediri(nz:1:-1)
         eup(1:nz) = eupi(nz:1:-1)
         edn(1:nz) = edni(nz:1:-1)
      ELSE
         CALL  PSNDO( dsdh, nid,
     $        NLYR, DTAUC, SSALB, PMOM, 
     $        ALBEDO, NSTR, umu0,
     $        RFLDIR,RFLDN, FLUP, 
     $        uavgso, uavgup, uavgdn )

* output (top-down -> buttom-up)
        edir(1:nz) = rfldir(nz:1:-1)
        edn(1:nz)  = rfldn(nz:1:-1)
        eup(1:nz)  = flup(nz:1:-1)
        fdir(1:nz) =  fourPI * uavgso(nz:1:-1)
        fdn(1:nz)  =  fourPI * uavgdn(nz:1:-1)
        fup(1:nz)  =  fourPI * uavgup(nz:1:-1)
      ENDIF

      END SUBROUTINE rtlink

*=============================================================================*

      SUBROUTINE PSNDO( dsdh, nid,
     $                  NLYR, DTAUC, SSALB, PMOM,
     $                  ALBEDO, NSTR, UMU0,
     $                  RFLDIR, RFLDN, FLUP,
     $                  uavgso, uavgup, uavgdn )


      use tuv_params, only : PI

      INTEGER, intent(in) :: NLYR
      INTEGER, intent(in) :: NSTR
      INTEGER, intent(in) :: nid(0:)
      REAL, intent(in)    :: ALBEDO
      REAL, intent(in)    :: UMU0
      REAL, intent(in)    :: dsdh(0:,:)
      REAL, intent(in)    :: DTAUC(:)
      REAL, intent(inout) :: SSALB(:)
      REAL, intent(inout) :: PMOM(0:,:)
      REAL, intent(out)   :: RFLDIR(:), RFLDN(:), FLUP(:)
      REAL, intent(out)   :: uavgso(:)
      REAL, intent(out)   :: uavgup(:)
      REAL, intent(out)   :: uavgdn(:)

      REAL, PARAMETER     :: RPD  = PI / 180.0
c     spherical geometry
      REAL tausla(0:NLYR), tauslau(0:NLYR), mu2(0:NLYR)
c     ..

c     local variables
      REAL ::   DFDT(NLYR+1),
     &          HL(0:NSTR), PHI(MXPHI),
     &          TRNMED(NSTR), U0U(nstr,nlyr+1), UAVG(NLYR+1),
     &          UMU(NSTR), CWT(NSTR), UTAU(NLYR+1),
     &          UU(NSTR,NLYR,MXPHI)
c     ..
c     .. Local Scalars ..

      LOGICAL :: COMPAR, LYRCUT, PASS1
      INTEGER :: IQ, IU, J, KCONV, L, LC, LEV, LU, MAZIM, NAZ, NCOL,
     &           NCOS, NCUT, NN
      REAL    :: AZERR, AZTERM, BPLANK, COSPHI, DELM0,
     &           DUM, SGN, TPLANK
c     ..
c     .. Local Arrays ..

      INTEGER :: NTAU, NUMU
      INTEGER :: IPVT(NSTR*NLYR), LAYRU(NLYR+1)

      REAL ::   ANGCOS(1)
      REAL ::   AMB(NSTR/2,NSTR/2), APB(NSTR/2,NSTR/2),
     &          ARRAY(NSTR,NSTR),
     &          B(NSTR*NLYR), BDR(NSTR/2,0:NSTR/2), BEM(NSTR/2),
     &          CBAND(9*(NSTR/2)-2,NSTR*NLYR), CC(NSTR,NSTR),
     &          CMU(NSTR), DTAUCP(NLYR),
     &          EMU(NSTR), EVAL(NSTR/2), EVECC(NSTR, NSTR),
     &          EXPBEA(0:NLYR), FLDIR(NLYR+1), FLDN(NLYR+1),
     &          FLYR(NLYR), GC(NSTR,NSTR,NLYR),
     &          GL(0:NSTR,NLYR), GU(NSTR,NSTR,NLYR),
     &          HLPR(0:NSTR), KK(NSTR,NLYR), LL(NSTR,NLYR),
     &          OPRIM(NLYR), PHIRAD(MXPHI), PKAG(0:NLYR),
     &          PSI(NSTR), RMU(NSTR,0:NSTR/2), TAUC(0:NLYR),
     &          TAUCPR(0:NLYR), U0C(NSTR,NLYR+1), UTAUPR(NLYR+1),
     &          UUM(NSTR,NLYR), WK(NSTR), XR0(NLYR),
     &          XR1(NLYR), YLM0(0:NSTR,1), YLMC(0:NSTR,NSTR),
     &          YLMU(0:NSTR,NSTR), Z(NSTR*NLYR), Z0(NSTR),
     &          Z0U(NSTR,NLYR), Z1(NSTR), Z1U(NSTR,NLYR),
     &          ZBEAM(NSTR,NLYR), ZJ(NSTR),
     &          ZPLK0(NSTR,NLYR), ZPLK1(NSTR,NLYR), ZZ(NSTR,NLYR)

      REAL :: sindir(nlyr+1), sinup(nlyr+1), sindn(nlyr+1)

cgy added glsave and dgl to allow adjustable dimensioning in SOLVEC
      REAL GLSAVE(0:nstr), DGL(0:nstr)

      REAL(8) :: AAD(NSTR/2,NSTR/2), EVECCD(NSTR/2,NSTR/2),
     &           EVALD(NSTR/2), WKD(NSTR)

      REAL    :: PLKAVG

      SAVE      PASS1
      DATA      PASS1 / .TRUE. /

      IF( PASS1 ) THEN
         DITHER = 10.*R1MACH( 4 )
c                            ** Must dither more on Cray (14-digit prec)
         IF( DITHER < 1.E-10 ) DITHER = 10.*DITHER
         PASS1 = .FALSE.
      END IF
 
c     ** Calculate cumulative optical depth
c     and dither single-scatter albedo
c     to improve numerical behavior of
c     eigen{value/vector} computation

      TAUC = rZERO

      DO LC = 1, NLYR
         IF( SSALB(LC) == rONE ) THEN
           SSALB(LC) = rONE - DITHER
         ENDIF
         TAUC(LC) = TAUC(LC - 1) + DTAUC(LC)
      ENDDO
c                                ** Check input dimensions and variables

      CALL CHEKIN( NLYR, DTAUC, SSALB, PMOM,
     $             NTAU, UTAU, NSTR, NUMU, UMU, NPHI,
     $             PHI, UMU0, FISOT, ALBEDO,
     $             HL, TAUC )

c                                 ** Zero internal and output arrays

      CALL  ZEROAL( EXPBEA(1:), FLYR, OPRIM, TAUCPR(1:), XR0, XR1,
     $              CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     $              HLPR, YLM0(:,1), ARRAY, CC, EVECC,
     $              GL, YLMC, YLMU,
     $              KK, LL, ZZ, ZPLK0, ZPLK1,
     $              GC, LAYRU, UTAUPR,
     $              GU, Z0U, Z1U, ZBEAM,
     $              EVAL, AMB, APB, IPVT, Z,
     $              RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     $              TRNMED, U0U, UU )

c                                 ** Perform various setup operations

      CALL SETDIS( 
     $    dsdh, nid, tausla, tauslau, mu2,
     &    CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM, FLYR,
     &    GL, HL, HLPR, IBCND, LAMBER, LAYRU, LYRCUT,
     &    NCUT, NLYR, NTAU, NN, NSTR, PLANK,
     &    NUMU, ONLYFL, OPRIM, PMOM, SSALB, TAUC, TAUCPR, UTAU,
     &    UTAUPR, UMU, UMU0, USRTAU, USRANG )

c                                 ** Print input information
      IF ( PRNT(1) ) THEN
           CALL PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM,
     $                  NTAU, UTAU, NSTR, NUMU, UMU,
     $                  PHI, UMU0, FISOT, ALBEDO, HL,
     $                  FLYR, LYRCUT,
     $                  OPRIM, TAUC, TAUCPR, PRNT(7) )
      ENDIF

c                              ** Handle special case for getting albedo
c                                 and transmissivity of medium for many
c                                 beam angles at once
c                                   ** Calculate Planck functions

         BPLANK = rZERO
         TPLANK = rZERO
         PKAG   = rZERO

c ========  BEGIN LOOP TO SUM AZIMUTHAL COMPONENTS OF INTENSITY  =======
c           (EQ STWJ 5)

      KCONV  = 0
c                                    ** Azimuth-independent case

      IF( FBEAM == rZERO .OR. (1. - UMU0) < 1.E-5 .OR. ONLYFL .OR.
     $    (NUMU == 1 .AND. (1. - UMU(1)) < 1.E-5 ) ) THEN
        NAZ = 0
      ELSE
        NAZ  = NSTR - 1
      ENDIF

      AZIMUTH_LOOP: DO MAZIM = 0, NAZ

         IF( MAZIM.EQ.0 ) DELM0  = rONE
         IF( MAZIM.GT.0 ) DELM0  = rZERO

c                             ** Get normalized associated Legendre
c                                polynomials for
c                                (a) incident beam angle cosine
c                                (b) computational and user polar angle
c                                    cosines
         IF( FBEAM > rZERO ) THEN
            NCOS   = 1
            ANGCOS = -UMU0
            CALL LEPOLY( NCOS, MAZIM, NSTR - 1, ANGCOS, YLM0 )
         END IF


         IF( .NOT. ONLYFL .AND. USRANG ) THEN
            CALL LEPOLY( NUMU, MAZIM, NSTR-1, UMU, YLMU )
         ENDIF

         CALL LEPOLY( NN, MAZIM, NSTR-1, CMU, YLMC )

c                       ** Get normalized associated Legendre polys.
c                          with negative arguments from those with
c                          positive arguments; Dave/Armstrong Eq. (15)
         SGN  = -rONE

         DO L = MAZIM, NSTR - 1
            SGN  = - SGN
            DO IQ = NN + 1, NSTR
               YLMC(L,IQ) = SGN*YLMC(L,IQ - NN)
            ENDDO
         ENDDO
c     ** Specify users bottom reflectivity
c        and emissivity properties
         IF ( .NOT. LYRCUT ) THEN
           CALL  SURFAC( 
     $         ALBEDO, DELM0, FBEAM, HLPR, LAMBER,
     $         MAZIM, NN, NUMU, NSTR, ONLYFL,
     $         UMU, USRANG, YLM0(:,1), YLMC, YLMU, BDR, EMU, BEM, RMU )
         ENDIF


c ===================  BEGIN LOOP ON COMPUTATIONAL LAYERS  =============

         DO LC = 1, NCUT
            CALL SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL(0:,LC),
     &           MAZIM, NN, NSTR, YLM0(:,1), YLMC, CC, 
     &           EVECC, EVAL, KK(:,LC ), GC(:,:,LC), AAD, EVECCD,
     &           EVALD, WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0,
     &           ZJ, ZZ(:,LC), OPRIM(LC), LC, mu2(lc),
     &           glsave, dgl)
         ENDDO


c ===================  END LOOP ON COMPUTATIONAL LAYERS  ===============


c                      ** Set coefficient matrix of equations combining
c                         boundary and layer interface conditions

         CALL SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     $                LAMBER, LYRCUT, NCOL, NCUT,
     $                NN, NSTR, TAUCPR, WK )

c                      ** Solve for constants of integration in homo-
c                         geneous solution (general boundary conditions)

         CALL SOLVE0( 
     $       B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     &       FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,
     &       NCOL, NCUT, NN, NSTR, PI,
     &       TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

c                                  ** Compute upward and downward fluxes

         IF ( MAZIM == 0 ) THEN
            CALL FLUXES( tausla, tauslau,
     $                   CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                   NCUT, NN, NSTR, NTAU,
     $                   PI, PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR,
     $                   XR0, XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP,
     $                   FLDN, FLDIR, RFLDIR, RFLDN, UAVG, U0C,
     $                   uavgso, uavgup, uavgdn,
     $                   sindir, sinup, sindn)
         ENDIF

         IF( ONLYFL ) THEN
            IF( MXUMU >= NSTR ) THEN
c     ** Save azimuthal-avg intensities
c        at quadrature angles
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
c     ** Save azimuthally averaged intensities
            DO LU = 1, NTAU
               DO IU = 1, NUMU
                  U0U(IU,LU) = UUM(IU,LU)
                  DO J = 1, NPHI
                     UU(IU,LU,J) = UUM(IU,LU)
                  ENDDO
               ENDDO
            ENDDO
c                              ** Print azimuthally averaged intensities
c                                 at user angles

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
c                                ** Increment intensity by current
c                                   azimuthal component (Fourier
c                                   cosine series);  Eq SD(2)
            AZERR  = rZERO

            DO J = 1, NPHI
               COSPHI = COS( MAZIM*PHIRAD(J) )
               DO LU = 1, NTAU
                  DO IU = 1, NUMU
                     AZTERM = UUM(IU,LU)*COSPHI
                     UU(IU,LU,J) = UU(IU,LU,J) + AZTERM
                     AZERR = MAX( AZERR,
     $                       RATIO( ABS(AZTERM), ABS(UU(IU,LU,J)) ) )
                  ENDDO
               ENDDO
            ENDDO

            IF( AZERR <= ACCUR ) KCONV  = KCONV + 1

            IF( KCONV >= 2 ) THEN
               EXIT AZIMUTH_LOOP
            ENDIF

         ENDIF

      ENDDO AZIMUTH_LOOP

c ===================  END LOOP ON AZIMUTHAL COMPONENTS  ===============


c                                          ** Print intensities
      IF( PRNT( 5 ) .AND. .NOT. ONLYFL ) THEN
        CALL PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )
      ENDIF

      END SUBROUTINE PSNDO

      SUBROUTINE ASYMTX( AA, EVEC, EVAL, IER, WKD, AAD,
     $                   EVECD, EVALD )

c    =======  D O U B L E    P R E C I S I O N    V E R S I O N  ======

c       Solves eigenfunction problem for real asymmetric matrix
c       for which it is known a priori that the eigenvalues are real.

c       This is an adaptation of a subroutine EIGRF in the IMSL
c       library to use real instead of complex arithmetic, accounting
c       for the known fact that the eigenvalues and eigenvectors in
c       the discrete ordinate solution are real.  Other changes include
c       putting all the called subroutines in-line, deleting the
c       performance index calculation, updating many DO-loops
c       to Fortran77, and in calculating the machine precision
c       TOL instead of specifying it in a data statement.

c       EIGRF is based primarily on EISPACK routines.  The matrix is
c       first balanced using the Parlett-Reinsch algorithm.  Then
c       the Martin-Wilkinson algorithm is applied.

c       References:
c          Dongarra, J. and C. Moler, EISPACK -- A Package for Solving
c             Matrix Eigenvalue Problems, in Cowell, ed., 1984:
c             Sources and Development of Mathematical Software,
c             Prentice-Hall, Englewood Cliffs, NJ
c         Parlett and Reinsch, 1969: Balancing a Matrix for Calculation
c             of Eigenvalues and Eigenvectors, Num. Math. 13, 293-304
c         Wilkinson, J., 1965: The Algebraic Eigenvalue Problem,
c             Clarendon Press, Oxford

c   I N P U T    V A R I A B L E S:

c       AA    :  input asymmetric matrix, destroyed after solved
c        M    :  order of  AA
c       IA    :  first dimension of  AA
c    IEVEC    :  first dimension of  EVEC

c   O U T P U T    V A R I A B L E S:

c       EVEC  :  (unnormalized) eigenvectors of  AA
c                   ( column J corresponds to EVAL(J) )

c       EVAL  :  (unordered) eigenvalues of AA ( dimension at least M )

c       IER   :  if .NE. 0, signals that EVAL(IER) failed to converge;
c                   in that case eigenvalues IER+1,IER+2,...,M  are
c                   correct but eigenvalues 1,...,IER are set to zero.

c   S C R A T C H   V A R I A B L E S:

c       WKD   :  work area ( dimension at least 2*M )
c       AAD   :  double precision stand-in for AA
c       EVECD :  double precision stand-in for EVEC
c       EVALD :  double precision stand-in for EVAL

c   Called by- SOLEIG
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER, intent(out) :: IER
c     ..
c     .. Array Arguments ..

      REAL, intent(in)  :: AA(:,:)
      REAL, intent(out) :: EVAL(:), EVEC(:,:)
      REAL(8), intent(inout) :: WKD(:)
      REAL(8), intent(out)   :: AAD(:,:)
      REAL(8), intent(out)   :: EVALD(:), EVECD(:,:)
c     ..
c     .. Local Scalars ..

      INTEGER :: I, II, IN, J, K, KA, KKK, L, LB, LLL, N, N1, N2
      INTEGER :: M, IA, IEVEC
      LOGICAL :: NOCONV, NOTLAS

      REAL(8) :: COL, DISCRI, F, G, H,
     $           P, Q, R, REPL, RNORM, ROW, S, SCALE, SGN, T,
     $           TOL, UU, VV, W, X, Y, Z

      REAL(8), parameter :: C1 = 0.4375_8
      REAL(8), parameter :: C2 = 0.5_8
      REAL(8), parameter :: C3 = 0.75_8
      REAL(8), parameter :: C4 = 0.95_8
      REAL(8), parameter :: C5 = 16._8
      REAL(8), parameter :: C6 = 256._8
      REAL(8), parameter :: ZERO = 0._8
      REAL(8), parameter :: ONE = 1._8

      IER  = 0
      TOL  = D1MACH( 4 )

      M     = size(EVAL)
      IA    = size(AA,dim=1)
      IEVEC = size(EVECD,dim=1)

      IF( M < 1 .OR. IA < M .OR. IEVEC < M )
     &    CALL ERRMSG( 'ASYMTX--bad input variable(s)', .TRUE. )

c                           ** Handle 1x1 and 2x2 special cases

      IF( M == 1 ) THEN

         EVAL( 1 )    = AA( 1, 1 )
         EVEC( 1, 1 ) = rONE
         RETURN

      ELSE IF( M == 2 ) THEN

         DISCRI = ( AA( 1,1 ) - AA( 2,2 ) )**2
     &              + 4.*AA( 1, 2 )*AA( 2, 1 )

         IF( DISCRI.LT.rZERO )
     &       CALL ERRMSG( 'ASYMTX--complex evals in 2x2 case',.TRUE. )

         SGN  = rONE

         IF( AA( 1,1 ) < AA( 2,2 ) ) SGN  = -rONE

         EVAL( 1 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) + SGN*SQRT( DISCRI ) )
         EVAL( 2 ) = 0.5*( AA( 1,1 ) + AA( 2,2 ) - SGN*SQRT( DISCRI ) )
         EVEC( 1, 1 ) = rONE
         EVEC( 2, 2 ) = rONE

         IF( AA( 1,1 ) == AA( 2,2 ) .AND.
     $       ( AA( 2,1 ) == rZERO .OR. AA( 1,2 ) == rZERO ) ) THEN

            RNORM  = ABS( AA( 1,1 ) ) + ABS( AA( 1,2 ) ) +
     $               ABS( AA( 2,1 ) ) + ABS( AA( 2,2 ) )
            W  = TOL*RNORM
            EVEC( 2, 1 ) =   AA( 2, 1 ) / W
            EVEC( 1, 2 ) = - AA( 1, 2 ) / W
         ELSE
            EVEC( 2, 1 ) = AA( 2, 1 ) / ( EVAL( 1 ) - AA( 2,2 ) )
            EVEC( 1, 2 ) = AA( 1, 2 ) / ( EVAL( 2 ) - AA( 1,1 ) )
         END IF

         RETURN

      END IF
c                               ** Put s.p. matrix into d.p. matrix
      DO J = 1, M
         DO K = 1, M
            AAD(J,K) = REAL(AA(J,K),kind=8)
         ENDDO
      ENDDO

c                                ** Initialize output variables
      IER  = 0

      DO 40 I = 1, M
         EVALD(I) = ZERO

         DO 30 J = 1, M
            EVECD( I, J ) = ZERO
   30    CONTINUE

         EVECD( I, I ) = ONE
   40 CONTINUE

c                  ** Balance the input matrix and reduce its norm by
c                     diagonal similarity transformation stored in WK;
c                     then search for rows isolating an eigenvalue
c                     and push them down
      RNORM  = ZERO
      L  = 1
      K  = M

   50 CONTINUE
      KKK  = K

      DO 90 J = KKK, 1, -1

         ROW  = ZERO

         DO 60 I = 1, K

            IF( I.NE.J ) ROW  = ROW + ABS( AAD( J,I ) )

   60    CONTINUE

         IF( ROW.EQ.ZERO ) THEN

            WKD( K ) = J

            IF( J.NE.K ) THEN

               DO 70 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, K )
                  AAD( I, K ) = REPL
   70          CONTINUE

               DO 80 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( K, I )
                  AAD( K, I ) = REPL
   80          CONTINUE

            END IF

            K  = K - 1
            GO TO  50

         END IF

   90 CONTINUE
c                                ** Search for columns isolating an
c                                   eigenvalue and push them left
  100 CONTINUE
      LLL  = L

      DO 140 J = LLL, K

         COL  = ZERO

         DO 110 I = L, K

            IF( I.NE.J ) COL  = COL + ABS( AAD( I,J ) )

  110    CONTINUE

         IF( COL.EQ.ZERO ) THEN

            WKD( L ) = J

            IF( J.NE.L ) THEN

               DO 120 I = 1, K
                  REPL        = AAD( I, J )
                  AAD( I, J ) = AAD( I, L )
                  AAD( I, L ) = REPL
  120          CONTINUE

               DO 130 I = L, M
                  REPL        = AAD( J, I )
                  AAD( J, I ) = AAD( L, I )
                  AAD( L, I ) = REPL
  130          CONTINUE

            END IF

            L  = L + 1
            GO TO  100

         END IF

  140 CONTINUE

c                           ** Balance the submatrix in rows L through K
      DO 150 I = L, K
         WKD( I ) = ONE
  150 CONTINUE

  160 CONTINUE
      NOCONV = .FALSE.

      DO 220 I = L, K

         COL  = ZERO
         ROW  = ZERO

         DO 170 J = L, K

            IF( J.NE.I ) THEN

               COL  = COL + ABS( AAD( J,I ) )
               ROW  = ROW + ABS( AAD( I,J ) )

            END IF

  170    CONTINUE

         F  = ONE
         G  = ROW / C5
         H  = COL + ROW

  180    CONTINUE
         IF( COL.LT.G ) THEN

            F    = F*C5
            COL  = COL*C6
            GO TO  180

         END IF

         G  = ROW*C5

  190    CONTINUE
         IF( COL.GE.G ) THEN

            F    = F / C5
            COL  = COL / C6
            GO TO  190

         END IF
c                                                ** Now balance
         IF( ( COL + ROW ) / F.LT.C4*H ) THEN

            WKD( I ) = WKD( I )*F
            NOCONV = .TRUE.

            DO 200 J = L, M
               AAD( I, J ) = AAD( I, J ) / F
  200       CONTINUE

            DO 210 J = 1, K
               AAD( J, I ) = AAD( J, I )*F
  210       CONTINUE

         END IF

  220 CONTINUE


      IF( NOCONV ) GO TO  160
c                                   ** Is A already in Hessenberg form?
      IF( K-1 .LT. L+1 ) GO TO  370

c                                   ** Transfer A to a Hessenberg form
      DO 310 N = L + 1, K - 1

         H  = ZERO
         WKD( N + M ) = ZERO
         SCALE  = ZERO
c                                                 ** Scale column
         DO 230 I = N, K
            SCALE  = SCALE + ABS( AAD( I,N - 1 ) )
  230    CONTINUE

         IF( SCALE.NE.ZERO ) THEN

            DO 240 I = K, N, -1
               WKD( I + M ) = AAD( I, N - 1 ) / SCALE
               H  = H + WKD( I + M )**2
  240       CONTINUE

            G    = - SIGN( SQRT( H ), WKD( N + M ) )
            H    = H - WKD( N + M )*G
            WKD( N + M ) = WKD( N + M ) - G
c                                            ** Form (I-(U*UT)/H)*A
            DO 270 J = N, M

               F  = ZERO

               DO 250 I = K, N, -1
                  F  = F + WKD( I + M )*AAD( I, J )
  250          CONTINUE

               DO 260 I = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( I + M )*F / H
  260          CONTINUE

  270       CONTINUE
c                                    ** Form (I-(U*UT)/H)*A*(I-(U*UT)/H)
            DO 300 I = 1, K

               F  = ZERO

               DO 280 J = K, N, -1
                  F  = F + WKD( J + M )*AAD( I, J )
  280          CONTINUE

               DO 290 J = N, K
                  AAD( I, J ) = AAD( I, J ) - WKD( J + M )*F / H
  290          CONTINUE

  300       CONTINUE

            WKD( N + M ) = SCALE*WKD( N + M )
            AAD( N, N - 1 ) = SCALE*G

         END IF

  310 CONTINUE


      DO 360 N = K - 2, L, -1

         N1   = N + 1
         N2   = N + 2
         F  = AAD( N + 1, N )

         IF( F.NE.ZERO ) THEN

            F  = F*WKD( N + 1 + M )

            DO 320 I = N + 2, K
               WKD( I + M ) = AAD( I, N )
  320       CONTINUE

            IF( N + 1.LE.K ) THEN

               DO 350 J = 1, M

                  G  = ZERO

                  DO 330 I = N + 1, K
                     G  = G + WKD( I + M )*EVECD( I, J )
  330             CONTINUE

                  G  = G / F

                  DO 340 I = N + 1, K
                     EVECD( I, J ) = EVECD( I, J ) + G*WKD( I + M )
  340             CONTINUE

  350          CONTINUE

            END IF

         END IF

  360 CONTINUE


  370 CONTINUE

      N  = 1

      DO 390 I = 1, M

         DO 380 J = N, M
            RNORM  = RNORM + ABS( AAD( I,J ) )
  380    CONTINUE

         N  = I

         IF( I.LT.L .OR. I.GT.K ) EVALD( I ) = AAD( I, I )

  390 CONTINUE

      N  = K
      T  = ZERO

c                                      ** Search for next eigenvalues
  400 CONTINUE
      IF( N.LT.L ) GO TO  550

      IN  = 0
      N1  = N - 1
      N2  = N - 2
c                          ** Look for single small sub-diagonal element
  410 CONTINUE

      DO 420 I = L, N
         LB  = N + L - I

         IF( LB.EQ.L ) GO TO  430

         S  = ABS( AAD( LB - 1,LB - 1 ) ) + ABS( AAD( LB,LB ) )

         IF( S.EQ.ZERO ) S  = RNORM

         IF( ABS( AAD( LB, LB-1 ) ).LE. TOL*S ) GO TO  430

  420 CONTINUE


  430 CONTINUE
      X  = AAD( N, N )

      IF( LB.EQ.N ) THEN
c                                        ** One eigenvalue found
         AAD( N, N ) = X + T
         EVALD( N ) = AAD( N, N )
         N  = N1
         GO TO  400

      END IF

C next line has been included to avoid run time error caused by xlf

      IF ( ( N1.LE.0 ).OR.( N.LE.0 ) ) THEN
        WRITE(0,*) 'Subscript out of bounds in ASYMTX'
        STOP 9999
      ENDIF

      Y  = AAD( N1, N1 )
      W  = AAD( N, N1 )*AAD( N1, N )

      IF( LB.EQ.N1 ) THEN
c                                        ** Two eigenvalues found
         P  = ( Y - X )*C2
         Q  = P**2 + W
         Z  = SQRT( ABS( Q ) )
         AAD( N, N ) = X + T
         X  = AAD( N, N )
         AAD( N1, N1 ) = Y + T
c                                        ** Real pair
         Z  = P + SIGN( Z, P )
         EVALD( N1 ) = X + Z
         EVALD( N ) = EVALD( N1 )

         IF( Z.NE.ZERO ) EVALD( N ) = X - W / Z

         X  = AAD( N, N1 )
c                                  ** Employ scale factor in case
c                                     X and Z are very small
         R  = SQRT( X*X + Z*Z )
         P  = X / R
         Q  = Z / R
c                                             ** Row modification
         DO 440 J = N1, M
            Z  = AAD( N1, J )
            AAD( N1, J ) = Q*Z + P*AAD( N, J )
            AAD( N, J ) = Q*AAD( N, J ) - P*Z
  440    CONTINUE
c                                             ** Column modification
         DO 450 I = 1, N
            Z  = AAD( I, N1 )
            AAD( I, N1 ) = Q*Z + P*AAD( I, N )
            AAD( I, N ) = Q*AAD( I, N ) - P*Z
  450    CONTINUE
c                                          ** Accumulate transformations
         DO 460 I = L, K
            Z  = EVECD( I, N1 )
            EVECD( I, N1 ) = Q*Z + P*EVECD( I, N )
            EVECD( I, N ) = Q*EVECD( I, N ) - P*Z
  460    CONTINUE

         N  = N2
         GO TO  400

      END IF


      IF( IN.EQ.30 ) THEN

c                    ** No convergence after 30 iterations; set error
c                       indicator to the index of the current eigenvalue
         IER  = N
         GO TO  700

      END IF
c                                                  ** Form shift
      IF( IN.EQ.10 .OR. IN.EQ.20 ) THEN

         T  = T + X

         DO 470 I = L, N
            AAD( I, I ) = AAD( I, I ) - X
  470    CONTINUE

         S  = ABS( AAD( N,N1 ) ) + ABS( AAD( N1,N2 ) )
         X  = C3*S
         Y  = X
         W  = -C1*S**2

      END IF


      IN  = IN + 1

c                ** Look for two consecutive small sub-diagonal elements

C inhibit vectorization by CF77, as this will cause a run time error

      DO 480 J = LB, N2
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

         IF( I.EQ.LB ) GO TO  490

         UU   = ABS( AAD( I, I-1 ) )*( ABS( Q ) + ABS( R ) )
         VV   = ABS( P ) * ( ABS( AAD( I-1, I-1 ) ) + ABS( Z ) +
     &                       ABS( AAD( I+1, I+1 ) ) )

         IF( UU .LE. TOL*VV ) GO TO  490

  480 CONTINUE

  490 CONTINUE
      AAD( I+2, I ) = ZERO

c                      ** fpp vectorization of this loop triggers
c                         array bounds errors, so inhibit
CFPP$ NOVECTOR L
      DO 500 J = I + 3, N
         AAD( J, J - 2 ) = ZERO
         AAD( J, J - 3 ) = ZERO
  500 CONTINUE

c             ** Double QR step involving rows K to N and columns M to N

      DO 540 KA = I, N1

         NOTLAS = KA.NE.N1

         IF( KA.EQ.I ) THEN

            S  = SIGN( SQRT( P*P + Q*Q + R*R ), P )

            IF( LB.NE.I ) AAD( KA, KA - 1 ) = -AAD( KA, KA - 1 )

         ELSE

            P  = AAD( KA, KA - 1 )
            Q  = AAD( KA + 1, KA - 1 )
            R  = ZERO

            IF( NOTLAS ) R  = AAD( KA + 2, KA - 1 )

            X  = ABS( P ) + ABS( Q ) + ABS( R )

            IF( X.EQ.ZERO ) GO TO  540

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
c                                              ** Row modification
         DO 510 J = KA, M

            P  = AAD( KA, J ) + Q*AAD( KA + 1, J )

            IF( NOTLAS ) THEN

               P  = P + R*AAD( KA + 2, J )
               AAD( KA + 2, J ) = AAD( KA + 2, J ) - P*Z

            END IF

            AAD( KA + 1, J ) = AAD( KA + 1, J ) - P*Y
            AAD( KA, J ) = AAD( KA, J ) - P*X
  510    CONTINUE
c                                                 ** Column modification
         DO 520 II = 1, MIN( N, KA + 3 )

            P  = X*AAD( II, KA ) + Y*AAD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*AAD( II, KA + 2 )
               AAD( II, KA + 2 ) = AAD( II, KA + 2 ) - P*R

            END IF

            AAD( II, KA + 1 ) = AAD( II, KA + 1 ) - P*Q
            AAD( II, KA ) = AAD( II, KA ) - P
  520    CONTINUE
c                                          ** Accumulate transformations
         DO 530 II = L, K

            P  = X*EVECD( II, KA ) + Y*EVECD( II, KA + 1 )

            IF( NOTLAS ) THEN

               P  = P + Z*EVECD( II, KA + 2 )
               EVECD( II, KA + 2 ) = EVECD( II, KA + 2 ) - P*R

            END IF

            EVECD( II, KA + 1 ) = EVECD( II, KA + 1 ) - P*Q
            EVECD( II, KA ) = EVECD( II, KA ) - P
  530    CONTINUE

  540 CONTINUE

      GO TO  410
c                     ** All evals found, now backsubstitute real vector
  550 CONTINUE

      IF( RNORM.NE.ZERO ) THEN

         DO 580 N = M, 1, -1
            N2   = N
            AAD( N, N ) = ONE

            DO 570 I = N - 1, 1, -1
               W  = AAD( I, I ) - EVALD( N )

               IF( W.EQ.ZERO ) W  = TOL*RNORM

               R  = AAD( I, N )

               DO 560 J = N2, N - 1
                  R  = R + AAD( I, J )*AAD( J, N )
  560          CONTINUE

               AAD( I, N ) = -R / W
               N2   = I
  570       CONTINUE

  580    CONTINUE
c                      ** End backsubstitution vectors of isolated evals
         DO 600 I = 1, M

            IF( I.LT.L .OR. I.GT.K ) THEN

               DO 590 J = I, M
                  EVECD( I, J ) = AAD( I, J )
  590          CONTINUE

            END IF

  600    CONTINUE
c                                   ** Multiply by transformation matrix
         IF( K.NE.0 ) THEN

            DO 630 J = M, L, -1

               DO 620 I = L, K
                  Z  = ZERO

                  DO 610 N = L, MIN( J, K )
                     Z  = Z + EVECD( I, N )*AAD( N, J )
  610             CONTINUE

                  EVECD( I, J ) = Z
  620          CONTINUE

  630       CONTINUE

         END IF

      END IF


      DO 650 I = L, K

         DO 640 J = 1, M
            EVECD( I, J ) = EVECD( I, J )*WKD( I )
  640    CONTINUE
  650 CONTINUE

c                           ** Interchange rows if permutations occurred
      DO 670 I = L-1, 1, -1

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 660 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  660       CONTINUE

         END IF

  670 CONTINUE


      DO 690 I = K + 1, M

         J  = WKD( I )

         IF( I.NE.J ) THEN

            DO 680 N = 1, M
               REPL   = EVECD( I, N )
               EVECD( I, N ) = EVECD( J, N )
               EVECD( J, N ) = REPL
  680       CONTINUE

         END IF

  690 CONTINUE

c                         ** Put results into output arrays
  700 CONTINUE

      DO 720 J = 1, M

         EVAL( J ) = EVALD( J )

         DO 710 K = 1, M
            EVEC( J, K ) = EVECD( J, K )
  710    CONTINUE

  720 CONTINUE

      END SUBROUTINE ASYMTX

      SUBROUTINE CHEKIN( NLYR, DTAUC, SSALB, PMOM,
     $                   NTAU, UTAU, NSTR, NUMU,
     $                   UMU, NPHI, PHI, UMU0,
     $                   FISOT, ALBEDO, HL, TAUC )

c           Checks the input dimensions and variables

c   Calls- WRTBAD, WRTDIM, DREF, ERRMSG
c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER, intent(in) ::   NLYR, NPHI, NSTR, NTAU, NUMU
      REAL, intent(in)    ::   ALBEDO, FISOT, UMU0
c     ..
c     .. Array Arguments ..

      REAL, intent(in) ::
     $          DTAUC(:), HL(0:), PHI(:),
     $          PMOM(0:,:), SSALB(:),
     $          TAUC(0:), UMU(:)
      REAL, intent(inout) :: UTAU(:)
c     ..
c     .. Local Scalars ..

      LOGICAL  :: INPERR
      INTEGER  :: IRMU, IU, J, K, LC, LU
      REAL     :: FLXALB, RMU

c     ..

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
         IF( ANY( PMOM(:,LC) < -rONE ) .or. 
     $       ANY( PMOM(:,LC) > rONE ) ) THEN
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

         IF( IBCND == 1 .AND. 2*NUMU > MXUMU )
     $       INPERR = WRTBAD( 'MXUMU' )

         DO IU = 1, NUMU
            IF( UMU(IU) < -rONE .OR. UMU(IU) > rONE .OR.
     $          UMU(IU) == rZERO ) INPERR = WRTBAD( 'UMU' )

            IF( IBCND == 1 .AND. UMU(IU) < rZERO )
     $          INPERR = WRTBAD( 'UMU' )

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

         IF( ANY( PHI(:) < rZERO ) .OR. ANY( PHI(:) > 360.0 ) )
     $          INPERR = WRTBAD( 'PHI' )

      END IF


      IF( IBCND.LT.0 .OR. IBCND.GT.1 ) INPERR = WRTBAD( 'IBCND' )

      IF( IBCND.EQ.0 ) THEN

         IF( FBEAM.LT.rZERO ) INPERR = WRTBAD( 'FBEAM' )

         IF( FBEAM.GT.rZERO .AND. abs(UMU0).GT.rONE )
     &       INPERR = WRTBAD( 'UMU0' )

         IF( FBEAM.GT.rZERO .AND. ( PHI0.LT.rZERO .OR.PHI0.GT.360.0 ) )
     &       INPERR = WRTBAD( 'PHI0' )

         IF( FISOT.LT.rZERO ) INPERR = WRTBAD( 'FISOT' )

c                    ** Make sure flux albedo at dense mesh of incident
c                       angles does not assume unphysical values
         DO IRMU = 0, 100
            RMU  = REAL(IRMU)*0.01
            FLXALB = DREF(RMU,HL,NSTR)
            IF( FLXALB < rZERO .OR. FLXALB > rONE )
     $          INPERR = WRTBAD( 'HL' )
            ENDDO



      ELSE IF( IBCND == 1 ) THEN

         IF( ALBEDO < rZERO .OR. ALBEDO > rONE )
     &       INPERR = WRTBAD( 'ALBEDO' )

      END IF

      IF( ACCUR < rZERO .OR. ACCUR > 1.E-2 ) 
     $      INPERR = WRTBAD( 'ACCUR' )

      IF( MXCLY.LT.NLYR ) INPERR = WRTDIM( 'MXCLY', NLYR )

      IF( IBCND /= 1 ) THEN

         IF( USRTAU .AND. MXULV.LT.NTAU )
     $       INPERR = WRTDIM( 'MXULV',NTAU )

         IF( .NOT.USRTAU .AND. MXULV .LT. NLYR + 1 )
     $       INPERR = WRTDIM( 'MXULV', NLYR + 1 )

      ELSE

         IF( MXULV.LT.2 ) INPERR = WRTDIM( 'MXULV', 2 )

      END IF

      IF( MXCMU.LT.NSTR ) INPERR = WRTDIM( 'MXCMU', NSTR )

      IF( USRANG .AND. MXUMU.LT.NUMU ) INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( USRANG .AND. IBCND.EQ.1 .AND.MXUMU.LT.2*NUMU )
     &    INPERR = WRTDIM( 'MXUMU', NUMU )

      IF( .NOT.USRANG .AND. MXUMU.LT.NSTR )
     &    INPERR = WRTDIM( 'MXUMU', NSTR )

      IF( .NOT.ONLYFL .AND. IBCND.NE.1 .AND. MXPHI.LT.NPHI )
     &    INPERR = WRTDIM( 'MXPHI', NPHI )

      IF( INPERR )
     $    CALL ERRMSG( 'DISORT--input and/or dimension errors',.True.)

      END SUBROUTINE CHEKIN

      SUBROUTINE FLUXES( tausla, tauslau,
     $                   CMU, CWT, FBEAM, GC, KK, LAYRU, LL, LYRCUT,
     $                   NCUT, NN, NSTR, NTAU, PI,
     $                   PRNT, SSALB, TAUCPR, UMU0, UTAU, UTAUPR, XR0,
     $                   XR1, ZZ, ZPLK0, ZPLK1, DFDT, FLUP, FLDN, FLDIR,
     $                   RFLDIR, RFLDN, UAVG, U0C,
     $                   uavgso, uavgup, uavgdn,
     $                   sindir, sinup, sindn)

c       Calculates the radiative fluxes, mean intensity, and flux
c       derivative with respect to optical depth from the m=0 intensity
c       components (the azimuthally-averaged intensity)

c    I N P U T     V A R I A B L E S:

c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LAYRU    :  Layer number of user level UTAU
c       LL       :  Constants of integration in Eq. SC(1), obtained
c                     by solving scaled version of Eq. SC(5);
c                     exponential term of Eq. SC(12) not included
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Number of computational layer where absorption
c                     optical depth exceeds ABSCUT
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       UTAUPR   :  Optical depths of user output levels in delta-M
c                     coordinates;  equal to UTAU if no delta-M
c       XR0      :  Expansion of thermal source function in Eq. SS(14)
c       XR1      :  Expansion of thermal source function Eqs. SS(16)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)


c                   O U T P U T     V A R I A B L E S:

c       U0C      :  Azimuthally averaged intensities
c                   ( at polar quadrature angles )
c       (RFLDIR, RFLDN, FLUP, DFDT, UAVG are DISORT output variables)


c                   I N T E R N A L       V A R I A B L E S:

c       DIRINT   :  Direct intensity attenuated
c       FDNTOT   :  Total downward flux (direct + diffuse)
c       FLDIR    :  Direct-beam flux (delta-M scaled)
c       FLDN     :  Diffuse down-flux (delta-M scaled)
c       FNET     :  Net flux (total-down - diffuse-up)
c       FACT     :  EXP( - UTAUPR / UMU0 )
c       PLSORC   :  Planck source function (thermal)
c       ZINT     :  Intensity of m = 0 case, in Eq. SC(1)

c   Called by- DISORT
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER, intent(in) :: NCUT, NN, NSTR, NTAU
      REAL, intent(in)    :: FBEAM, PI, UMU0
      LOGICAL, intent(in) :: LYRCUT
c     ..
c     .. Array Arguments ..

      LOGICAL, intent(in) :: PRNT(:)
      INTEGER, intent(in) :: LAYRU(:)
      REAL,    intent(in) :: CMU(:), CWT(:),
     &          GC(:,:,:), KK(:,:), LL(:,:),
     &          SSALB(:),  TAUCPR(0:),
     &          UTAU(:), UTAUPR(:), XR0(:), XR1(:),
     &          ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)
      REAL, intent(in)  :: tausla(0:), tauslau(0:)
      REAL, intent(out) :: U0C(:,:)
      REAL, intent(out) :: RFLDIR(:), RFLDN(:), FLUP(:)
      REAL, intent(out) :: DFDT(:), UAVG(:)
      REAL, intent(out) :: uavgso(:), uavgup(:), uavgdn(:)
      REAL, intent(out) :: sindir(:), sinup(:), sindn(:)
      REAL, intent(out) :: FLDIR(:), FLDN(:)
c     ..
c     .. Local Scalars ..

      INTEGER :: IQ, JQ, LU, LYU
      REAL    :: ANG1, ANG2, DIRINT, FACT, FDNTOT, FNET, PLSORC, ZINT
c     ..

      IF( PRNT( 2 ) ) WRITE( *, 9000 )
c                                          ** Zero DISORT output arrays
      U0C   = 0.
      FLDIR = 0.
      FLDN  = 0.
      uavgso = 0.
      uavgup = 0.
      uavgdn = 0.
      sindir = 0.
      sinup  = 0.
      sindn  = 0.

c    ** Loop over user levels
      LEVEL_LOOP: DO LU = 1, NTAU

         LYU  = LAYRU( LU )

         IF( LYRCUT .AND. LYU > NCUT ) THEN
c                                                ** No radiation reaches
c                                                ** this level
            FDNTOT = rZERO
            FNET   = rZERO
            PLSORC = rZERO
            GO TO  70

         END IF

         IF( FBEAM > rZERO ) THEN
 
            FACT  = EXP( - tausla(LU-1) )
            DIRINT       = FBEAM*FACT
            FLDIR( LU )  = UMU0*( FBEAM*FACT )
            RFLDIR( LU ) = UMU0*FBEAM * EXP( -tauslau(lu-1) )
            sindir( LU ) = SQRT(1.-UMU0*UMU0)*FBEAM * 
     $                     EXP( -tauslau(lu-1) )

         ELSE

            DIRINT       = rZERO
            FLDIR( LU )  = rZERO
            RFLDIR( LU ) = rZERO
            sindir( LU ) = rZERO

         END IF


         DO IQ = 1, NN

            ZINT   = rZERO

            DO JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     $                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     $                  TAUCPR( LYU ) ) )
            ENDDO

            DO JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
            ENDDO

            U0C( IQ, LU ) = ZINT

            IF( FBEAM > rZERO ) THEN
              U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
            ENDIF

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( NN + 1 - IQ )*U0C( IQ, LU )
            uavgdn(lu) = uavgdn(lu) + cwt(nn+1-iq) * u0c( iq,lu )
            sindn(lu)  = sindn(lu)  + cwt(nn+1-iq) * 
     &                   SQRT(1.-CMU(NN+1-IQ)*CMU(NN+1-IQ))*
     &                   U0C( IQ, LU )
            FLDN( LU ) = FLDN( LU ) + CWT( NN + 1 - IQ )*
     &                   CMU( NN + 1 - IQ )*U0C( IQ, LU )
         ENDDO


         DO IQ = NN + 1, NSTR

            ZINT   = rZERO

            DO JQ = 1, NN
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU ) ) )
            ENDDO

            DO JQ = NN + 1, NSTR
               ZINT   = ZINT + GC( IQ, JQ, LYU )*LL( JQ, LYU )*
     &                  EXP( -KK( JQ,LYU )*( UTAUPR( LU ) -
     &                  TAUCPR( LYU - 1 ) ) )
            ENDDO

            U0C( IQ, LU ) = ZINT

            IF( FBEAM > rZERO ) THEN
              U0C( IQ, LU ) = ZINT + ZZ( IQ, LYU )*FACT
            ENDIF

            U0C( IQ, LU ) = U0C( IQ, LU ) + ZPLK0( IQ, LYU ) +
     &                      ZPLK1( IQ, LYU )*UTAUPR( LU )
            UAVG( LU ) = UAVG( LU ) + CWT( IQ - NN )*U0C( IQ, LU )
            uavgup(lu) = uavgup(lu) + cwt(iq-nn) * u0c( iq,lu )
            sinup (lu) = sinup(lu)  + cwt(iq-nn) * 
     &                   SQRT(1.-CMU(IQ-NN)*CMU(IQ-NN))*
     &                   U0C( IQ, LU )
            FLUP( LU ) = FLUP( LU ) + CWT( IQ - NN )*CMU( IQ - NN )*
     &                   U0C( IQ, LU )
         ENDDO


         FLUP( LU )  = 2.*PI*FLUP( LU )
         FLDN( LU )  = 2.*PI*FLDN( LU )
         FDNTOT      = FLDN( LU ) + FLDIR( LU )
         FNET        = FDNTOT - FLUP( LU )
         RFLDN( LU ) = FDNTOT - RFLDIR( LU )
         UAVG( LU )  = ( 2.*PI*UAVG( LU ) + DIRINT ) / ( 4.*PI )
         uavgso( lu ) = dirint / (4.*pi)
         uavgup( lu ) = (2.0 * pi * uavgup(lu) )/ (4.*pi)
         uavgdn( lu)  = (2.0 * pi * uavgdn(lu) )/ (4.*pi)
         sindn ( lu ) = 2.*PI*sindn ( LU )
         sinup ( lu ) = 2.*PI*sinup ( LU )

         PLSORC      = XR0( LYU ) + XR1( LYU )*UTAUPR( LU )
         DFDT( LU )  = ( 1.- SSALB( LYU ) ) * 4.*PI *
     &                 ( UAVG( LU ) - PLSORC )

   70    CONTINUE
         IF( PRNT( 2 ) ) WRITE( *, FMT = 9010 ) UTAU( LU ), LYU,
     &       RFLDIR( LU ), RFLDN( LU ), FDNTOT, FLUP( LU ), FNET,
     &       UAVG( LU ), PLSORC, DFDT( LU )

      ENDDO LEVEL_LOOP


      IF( PRNT( 3 ) ) THEN

         WRITE( *, FMT = 9020 )

         DO LU = 1, NTAU

            WRITE( *, FMT = 9030 ) UTAU( LU )

            DO IQ = 1, NN
               ANG1   = 180./ PI* ACOS( CMU( 2*NN - IQ + 1 ) )
               ANG2   = 180./ PI* ACOS( CMU( IQ ) )
               WRITE( *, 9040 ) ANG1, CMU(2*NN-IQ+1), U0C(IQ,LU),
     $                          ANG2, CMU(IQ),        U0C(IQ+NN,LU)
            ENDDO
         ENDDO

      END IF


 9000 FORMAT( //, 21X,
     $ '<----------------------- FLUXES ----------------------->', /,
     $ '   Optical  Compu    Downward    Downward    Downward     ',
     $ ' Upward                    Mean      Planck   d(Net Flux)', /,
     $ '     Depth  Layer      Direct     Diffuse       Total     ',
     $ 'Diffuse         Net   Intensity      Source   / d(Op Dep)', / )
 9010 FORMAT( F10.4, I7, 1P, 7E12.3, E14.3 )
 9020 FORMAT( / , / , ' ******** AZIMUTHALLY AVERAGED INTENSITIES',
     &      ' ( at polar quadrature angles ) *******' )
 9030 FORMAT( /, ' Optical depth =', F10.4, //,
     $  '     Angle (deg)   cos(Angle)     Intensity',
     $  '     Angle (deg)   cos(Angle)     Intensity' )
 9040 FORMAT( 2( 0P,F16.4,F13.5,1P,E14.3 ) )

      END SUBROUTINE FLUXES

      SUBROUTINE LEPOLY( NMU, M, TWONM1, MU, YLM )

c       Computes the normalized associated Legendre polynomial,
c       defined in terms of the associated Legendre polynomial
c       Plm = P-sub-l-super-m as

c             Ylm(MU) = sqrt( (l-m)!/(l+m)! ) * Plm(MU)

c       for fixed order m and all degrees from l = m to TWONM1.
c       When m.GT.0, assumes that Y-sub(m-1)-super(m-1) is available
c       from a prior call to the routine.

c       REFERENCE: Dave, J.V. and B.H. Armstrong, Computations of
c                  High-Order Associated Legendre Polynomials,
c                  J. Quant. Spectrosc. Radiat. Transfer 10,
c                  557-562, 1970.  (hereafter D/A)

c       METHOD: Varying degree recurrence relationship.

c       NOTE 1: The D/A formulas are transformed by
c               setting  M = n-1; L = k-1.
c       NOTE 2: Assumes that routine is called first with  M = 0,
c               then with  M = 1, etc. up to  M = TWONM1.
c       NOTE 3: Loops are written in such a way as to vectorize.

c  I N P U T     V A R I A B L E S:

c       NMU    :  Number of arguments of YLM
c       M      :  Order of YLM
c       MAXMU  :  First dimension of YLM
c       TWONM1 :  Max degree of YLM
c       MU(i)  :  Arguments of YLM (i = 1 to NMU)

c       If M.GT.0, YLM(M-1,i) for i = 1 to NMU is assumed to exist
c       from a prior call.

c  O U T P U T     V A R I A B L E:

c       YLM(l,i) :  l = M to TWONM1, normalized associated Legendre
c                   polynomials evaluated at argument MU(i)

c   Called by- DISORT, ALBTRN, SURFAC
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER, PARAMETER ::  MAXSQT = 1000
c     ..
c     .. Scalar Arguments ..

      INTEGER, intent(in) :: M, NMU, TWONM1
c     ..
c     .. Array Arguments ..

      REAL, intent(in)  :: MU(:)
      REAL, intent(out) :: YLM(0:,:)
c     ..
c     .. Local Scalars ..

      LOGICAL   :: PASS1
      INTEGER   :: I, L, NS
      REAL      :: TMP1, TMP2
c     ..
c     .. Local Arrays ..

      REAL ::  SQT( MAXSQT )
c     ..
c     .. External Subroutines ..

      SAVE      SQT, PASS1
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN
         PASS1  = .FALSE.
         DO NS = 1, MAXSQT
            SQT( NS ) = SQRT( REAL( NS ) )
         ENDDO
      END IF

      IF( 2*TWONM1 > MAXSQT )
     &    CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)


      IF( M == 0 ) THEN
c                             ** Upward recurrence for ordinary
c                                Legendre polynomials
         DO I = 1, NMU
            YLM(0,I) = rONE
            YLM(1,I) = MU(I)
         ENDDO

         DO L = 2, TWONM1
            DO I = 1, NMU
               YLM(L,I) = ((2*L - 1)*MU(I)*YLM(L - 1,I)
     $                           - (L - 1)*YLM(L - 2,I)) / L
            ENDDO
         ENDDO
      ELSE
         DO I = 1, NMU
c                               ** Y-sub-m-super-m; derived from
c                               ** D/A Eqs. (11,12)
            YLM(M,I) = -SQT(2*M - 1) / SQT(2*M )
     $                  *SQRT(rONE - MU(I)**2)*YLM(M - 1,I)

c                              ** Y-sub-(m+1)-super-m; derived from
c                              ** D/A Eqs.(13,14) using Eqs.(11,12)
            YLM(M + 1,I) = SQT(2*M + 1)*MU(I)*YLM(M,I)
         ENDDO
c                                   ** Upward recurrence; D/A EQ.(10)
         DO L = M + 2, TWONM1
            TMP1 = SQT(L - M )*SQT(L + M)
            TMP2 = SQT(L - M - 1)*SQT(L + M - 1)
            DO I = 1, NMU
               YLM(L,I) = 
     $          ((2*L - 1 )*MU(I)*YLM(L-1,I) - TMP2*YLM(L-2,I)) / TMP1
            ENDDO
         ENDDO
      END IF

      END SUBROUTINE LEPOLY

      SUBROUTINE PRAVIN( UMU, NUMU, UTAU, NTAU, U0U )

c        Print azimuthally averaged intensities at user angles

c   Called by- DISORT

c     LENFMT   Max number of polar angle cosines UMU that can be
c                printed on one line, as set in FORMAT statement
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      INTEGER   NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      U0U( :, : ), UMU( : ), UTAU( : )
c     ..
c     .. Local Scalars ..

      INTEGER   IU, IUMAX, IUMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NUMU.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *******  AZIMUTHALLY AVERAGED INTENSITIES ' //
     &   '(at user polar angles)  ********'

      LENFMT = 8
      NPASS  = 1 + (NUMU-1) / LENFMT

      WRITE( *,'(/,A,/,A)') '   Optical   Polar Angle Cosines',
     &                      '     Depth'

      DO 20 NP = 1, NPASS

         IUMIN  = 1 + LENFMT * ( NP - 1 )
         IUMAX  = MIN( LENFMT*NP, NUMU )
         WRITE( *,'(/,10X,8F14.5)') ( UMU(IU), IU = IUMIN, IUMAX )

         DO 10 LU = 1, NTAU
            WRITE( *, '(0P,F10.4,1P,8E14.4)' ) UTAU( LU ),
     &           ( U0U( IU,LU ), IU = IUMIN, IUMAX )
   10    CONTINUE

   20 CONTINUE


      END

      SUBROUTINE PRTINP( NLYR, DTAUC, DTAUCP, SSALB, PMOM,
     $                   NTAU, UTAU, NSTR, NUMU, UMU,
     $                   PHI, UMU0, FISOT,
     $                   ALBEDO, HL, FLYR, LYRCUT,
     $                   OPRIM, TAUC, TAUCPR, PRTMOM )

c        Print values of input variables

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Scalar Arguments ..

      LOGICAL, intent(in) :: LYRCUT, PRTMOM
      INTEGER, intent(in) :: NLYR, NSTR, NTAU, NUMU
      REAL, intent(in)    :: ALBEDO, FISOT, UMU0
c     ..
c     .. Array Arguments ..

      REAL, intent(in) ::  DTAUC(:), DTAUCP(:), FLYR(:), HL(0:),
     $          OPRIM(:), PHI(:), PMOM(0:,:), SSALB(:),
     $          TAUC(0:), TAUCPR(0:), UMU(:), UTAU(:)
c     ..
c     .. Local Scalars ..

      INTEGER :: IU, J, K, LC, LU
      REAL    :: YESSCT
c     ..


      WRITE( *, '(/,A,I4,A,I4)' ) ' No. streams =', NSTR,
     &       '     No. computational layers =', NLYR

      IF( IBCND /= 1 ) WRITE( *, '(I4,A,10F10.4,/,(26X,10F10.4))' )
     &    NTAU,' User optical depths :', ( UTAU(LU), LU = 1, NTAU )

      IF( .NOT. ONLYFL ) WRITE( *, '(I4,A,10F9.5,/,(31X,10F9.5))' )
     &    NUMU,' User polar angle cosines :',( UMU(IU), IU = 1, NUMU )

      IF( .NOT. ONLYFL .AND. IBCND /= 1 )
     &    WRITE( *, '(I4,A,10F9.2,/,(28X,10F9.2))' )
     &           NPHI,' User azimuthal angles :',( PHI(J), J = 1, NPHI )

      IF( .NOT. PLANK .OR. IBCND == 1 )
     &    WRITE( *, '(A)' ) ' No thermal emission'


      WRITE( *, '(A,I2)' ) ' Boundary condition flag: IBCND =', IBCND

      IF( IBCND == 0 ) THEN

         WRITE( *, '(A,1P,E11.3,A,0P,F8.5,A,F7.2,/,A,1P,E11.3)' )
     &          '    Incident beam with intensity =', FBEAM,
     &          ' and polar angle cosine = ', UMU0,
     &          '  and azimuth angle =', PHI0,
     &          '    plus isotropic incident intensity =', FISOT

         IF( LAMBER ) WRITE( *, '(A,0P,F8.4)' )
     &                '    Bottom albedo (Lambertian) =', ALBEDO

         IF( .NOT. LAMBER ) WRITE( *, '(A,/,(10X,10F9.5))' )
     &     '    Legendre coeffs of bottom bidirectional reflectivity :',
     &         ( HL( K ), K = 0, NSTR )

      ELSE IF( IBCND == 1 ) THEN

         WRITE(*,'(A)') '    Isotropic illumination from top and bottom'
         WRITE( *, '(A,0P,F8.4)' )
     &          '    Bottom albedo (Lambertian) =', ALBEDO
      END IF


      IF( DELTAM ) WRITE( *, '(A)' ) ' Uses delta-M method'
      IF( .NOT.DELTAM ) WRITE( *, '(A)' ) ' Does not use delta-M method'


      IF( IBCND == 1 ) THEN

         WRITE( *, '(A)' ) ' Calculate albedo and transmissivity of'//
     &                     ' medium vs. incident beam angle'

      ELSE IF( ONLYFL ) THEN

         WRITE( *, '(A)' )
     &          ' Calculate fluxes and azim-averaged intensities only'

      ELSE

         WRITE( *, '(A)' ) ' Calculate fluxes and intensities'

      END IF


      WRITE( *, '(A,1P,E11.2)' )
     &       ' Relative convergence criterion for azimuth series =',
     &       ACCUR

      IF( LYRCUT ) WRITE( *, '(A)' )
     &    ' Sets radiation = 0 below absorption optical depth 10'


c                                        ** Print layer variables
      IF( PLANK ) WRITE( *, FMT = 9180 )
      IF( .NOT. PLANK ) WRITE( *, FMT = 9190 )

      YESSCT = rZERO

      DO LC = 1, NLYR
         YESSCT = YESSCT + SSALB( LC )

         IF( PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4,F14.3)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM(1,LC)

         IF( .NOT.PLANK )
     &       WRITE(*,'(I4,2F10.4,F10.5,F12.5,2F10.4,F10.5,F9.4)')
     &             LC, DTAUC( LC ), TAUC( LC ), SSALB( LC ), FLYR( LC ),
     &             DTAUCP( LC ), TAUCPR( LC ), OPRIM( LC ), PMOM( 1,LC )
      ENDDO


      IF( PRTMOM .AND. YESSCT > rZERO ) THEN

         WRITE( *, '(/,A)' ) ' Layer   Phase Function Moments'

         DO LC = 1, NLYR
            IF( SSALB( LC ).GT.rZERO )
     &          WRITE( *, '(I6,10F11.6,/,(6X,10F11.6))' )
     &                 LC, ( PMOM( K, LC ), K = 0, NSTR )
         ENDDO

      END IF

c                ** (Read every other line in these formats)

 9180 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor   Temperature' )
 9190 FORMAT( /, 37X, '<------------- Delta-M --------------->', /,
     &'                   Total    Single                           ',
     &               'Total    Single', /,
     &'       Optical   Optical   Scatter   Truncated   ',
     &   'Optical   Optical   Scatter    Asymm', /,
     &'         Depth     Depth    Albedo    Fraction     ',
     &     'Depth     Depth    Albedo   Factor' )

      END

      SUBROUTINE PRTINT( UU, UTAU, NTAU, UMU, NUMU, PHI )

c         Prints the intensity at user polar and azimuthal angles

c     All arguments are DISORT input or output variables

c   Called by- DISORT

c     LENFMT   Max number of azimuth angles PHI that can be printed
c                on one line, as set in FORMAT statement
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER   NPHI, NTAU, NUMU
c     ..
c     .. Array Arguments ..

      REAL      PHI(:), UMU(:), UTAU(:), UU(:,:,:)
c     ..
c     .. Local Scalars ..

      INTEGER   IU, J, JMAX, JMIN, LENFMT, LU, NP, NPASS
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC MIN
c     ..


      IF( NPHI.LT.1 )  RETURN

      WRITE( *, '(//,A)' )
     &   ' *********  I N T E N S I T I E S  *********'

      LENFMT = 10
      NPASS  = 1 + (NPHI-1) / LENFMT

      WRITE( *, '(/,A,/,A,/,A)' )
     &   '             Polar   Azimuth angles (degrees)',
     &   '   Optical   Angle',
     &   '    Depth   Cosine'

      DO 30 LU = 1, NTAU

         DO 20 NP = 1, NPASS

            JMIN   = 1 + LENFMT * ( NP - 1 )
            JMAX   = MIN( LENFMT*NP, NPHI )

            WRITE( *, '(/,18X,10F11.2)' ) ( PHI(J), J = JMIN, JMAX )

            IF( NP.EQ.1 ) WRITE( *, '(F10.4,F8.4,1P,10E11.3)' )
     &             UTAU(LU), UMU(1), (UU(1, LU, J), J = JMIN, JMAX)
            IF( NP.GT.1 ) WRITE( *, '(10X,F8.4,1P,10E11.3)' )
     &                       UMU(1), (UU(1, LU, J), J = JMIN, JMAX)

            DO 10 IU = 2, NUMU
               WRITE( *, '(10X,F8.4,1P,10E11.3)' ) 
     &                 UMU( IU ), ( UU( IU, LU, J ), J = JMIN, JMAX )
   10       CONTINUE

   20    CONTINUE

   30 CONTINUE

      END SUBROUTINE PRTINT

      SUBROUTINE QGAUSN( M, GMU, GWT )

c       Compute weights and abscissae for ordinary Gaussian quadrature
c       on the interval (0,1);  that is, such that

c           sum(i=1 to M) ( GWT(i) f(GMU(i)) )

c       is a good approximation to

c           integral(0 to 1) ( f(x) dx )

c   INPUT :    M       order of quadrature rule

c   OUTPUT :  GMU(I)   array of abscissae (I = 1 TO M)
c             GWT(I)   array of weights (I = 1 TO M)

c   REFERENCE:  Davis, P.J. and P. Rabinowitz, Methods of Numerical
c                   Integration, Academic Press, New York, pp. 87, 1975

c   METHOD:  Compute the abscissae as roots of the Legendre
c            polynomial P-sub-M using a cubically convergent
c            refinement of Newton's method.  Compute the
c            weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note
c            that Newton's method can very easily diverge; only a
c            very good initial guess can guarantee convergence.
c            The initial guess used here has never led to divergence
c            even for M up to 1000.

c   ACCURACY:  relative error no better than TOL or computer
c              precision (machine epsilon), whichever is larger

c   INTERNAL VARIABLES:

c    ITER      : number of Newton Method iterations
c    MAXIT     : maximum allowed iterations of Newton Method
c    PM2,PM1,P : 3 successive Legendre polynomials
c    PPR       : derivative of Legendre polynomial
c    P2PRI     : 2nd derivative of Legendre polynomial
c    TOL       : convergence criterion for Legendre poly root iteration
c    X,XI      : successive iterates in cubically-convergent version
c                of Newtons Method (seeking roots of Legendre poly.)

c   Called by- SETDIS, SURFAC
c   Calls- D1MACH, ERRMSG
c +-------------------------------------------------------------------+

c     .. Scalar Arguments ..

      INTEGER, intent(in) :: M
c     ..
c     .. Array Arguments ..

      REAL, intent(out)  :: GMU(:), GWT(:)
c     ..
c     .. Local Scalars ..

      INTEGER, PARAMETER ::  MAXIT = 1000
      REAL(8), PARAMETER ::  ONE = 1.D0, TWO = 2.D0

      INTEGER ::  ITER, K, LIM, NN, NP1
      REAL    ::  CONA, PI, T
      REAL(8) ::  EN, NNP1, P, P2PRI, PM1, PM2, PPR, PROD,
     $            TMP, TOL, X, XI

      SAVE      PI, TOL

      DATA      PI / rZERO / 


      IF( PI == rZERO ) THEN
         PI   = 2.*ASIN( rONE )
         TOL  = 10.*D1MACH( 4 )
      END IF


      IF( M < 1 ) CALL ERRMSG( 'QGAUSN--Bad value of M',.True.)

      IF( M == 1 ) THEN
         GMU( 1 ) = 0.5
         GWT( 1 ) = rONE
         RETURN
      END IF

      EN   = M
      NP1  = M + 1
      NNP1 = M*NP1
      CONA = REAL( M - 1 ) / REAL( 8*M**3 )

      LIM  = M / 2

      DO 30 K = 1, LIM
c                                        ** Initial guess for k-th root
c                                           of Legendre polynomial, from
c                                           Davis/Rabinowitz (2.7.3.3a)
         T  = ( 4*K - 1 )*PI / ( 4*M + 2 )
         X  = COS( T + CONA / TAN( T ) )
         ITER = 0
c                                        ** Upward recurrence for
c                                           Legendre polynomials
   10    CONTINUE
         ITER   = ITER + 1
         PM2    = ONE
         PM1    = X

         DO 20 NN = 2, M
            P    = ( ( 2*NN - 1 )*X*PM1 - ( NN - 1 )*PM2 ) / NN
            PM2  = PM1
            PM1  = P
   20    CONTINUE
c                                              ** Newton Method
         TMP    = ONE / ( ONE - X**2 )
         PPR    = EN*( PM2 - X*P )*TMP
         P2PRI  = ( TWO*X*PPR - NNP1*P )*TMP
         XI     = X - ( P / PPR )*( ONE +
     &            ( P / PPR )*P2PRI / ( TWO*PPR ) )

c                                              ** Check for convergence
         IF( ABS( XI - X ) > TOL ) THEN

            IF( ITER.GT.MAXIT )
     &          CALL ERRMSG( 'QGAUSN--max iteration count',.True.)

            X  = XI
            GO TO  10

         END IF
c                             ** Iteration finished--calculate weights,
c                                abscissae for (-1,1)
         GMU( K ) = -X
         GWT( K ) = TWO / ( TMP*( EN*PM2 )**2 )
         GMU( NP1 - K ) = -GMU( K )
         GWT( NP1 - K ) = GWT( K )
   30 CONTINUE
c                                    ** Set middle abscissa and weight
c                                       for rules of odd order
      IF( MOD( M,2 ).NE.0 ) THEN

         GMU( LIM + 1 ) = 0.0
         PROD   = ONE

         DO 40 K = 3, M, 2
            PROD   = PROD * K / ( K - 1 )
   40    CONTINUE

         GWT( LIM + 1 ) = TWO / PROD**2
      END IF

c                                        ** Convert from (-1,1) to (0,1)
      DO 50 K = 1, M
         GMU( K ) = 0.5*GMU( K ) + 0.5
         GWT( K ) = 0.5*GWT( K )
   50 CONTINUE

      END SUBROUTINE QGAUSN

      SUBROUTINE SETDIS( dsdh, nid, tausla, tauslau, mu2,
     $                   CMU, CWT, DELTAM, DTAUC, DTAUCP, EXPBEA, FBEAM,
     $                   FLYR, GL, HL, HLPR, IBCND, LAMBER, LAYRU,
     $                   LYRCUT, NCUT, NLYR,
     $                   NTAU, NN, NSTR, PLANK, NUMU, ONLYFL, OPRIM,
     $                   PMOM, SSALB, TAUC, TAUCPR, UTAU, UTAUPR, UMU,
     $                   UMU0, USRTAU, USRANG )

c          Perform miscellaneous setting-up operations

c       INPUT :  all are DISORT input variables (see DOC file)

c       OUTPUT:  NTAU,UTAU   if USRTAU = FALSE
c                NUMU,UMU    if USRANG = FALSE
c                CMU,CWT     computational polar angles and
c                               corresponding quadrature weights
c                EXPBEA      transmission of direct beam
c                FLYR        truncated fraction in delta-M method
c                GL          phase function Legendre coefficients multi-
c                              plied by (2L+1) and single-scatter albedo
c                HLPR        Legendre moments of surface bidirectional
c                              reflectivity, times 2K+1
c                LAYRU       Computational layer in which UTAU falls
c                LYRCUT      flag as to whether radiation will be zeroed
c                              below layer NCUT
c                NCUT        computational layer where absorption
c                              optical depth first exceeds  ABSCUT
c                NN          NSTR / 2
c                OPRIM       delta-M-scaled single-scatter albedo
c                TAUCPR      delta-M-scaled optical depth
c                UTAUPR      delta-M-scaled version of  UTAU

c   Called by- DISORT
c   Calls- QGAUSN, ERRMSG
c ----------------------------------------------------------------------

      use tuv_params, only : largest

c     .. Scalar Arguments ..

      INTEGER, intent(in)  ::  IBCND, NLYR, NSTR
      INTEGER, intent(out) ::  NCUT, NN, NTAU, NUMU
      LOGICAL, intent(in)  ::  DELTAM, LAMBER, ONLYFL
      LOGICAL, intent(in)  ::  PLANK, USRANG, USRTAU
      LOGICAL, intent(out) ::  LYRCUT
      REAL, intent(in)     ::  FBEAM, UMU0

c geometry
      INTEGER, intent(in) :: nid(0:)
      REAL, intent(in)    :: dsdh(0:,:)
      REAL, intent(out)   :: tausla(0:), tauslau(0:), mu2(0:)

      REAL :: sum, sumu
c     ..
c     .. Array Arguments ..

      INTEGER, intent(out) :: LAYRU(:)
      REAL, intent(in)   :: DTAUC(:), HL(0:), SSALB(:), TAUC(0:)
      REAL, intent(out)  :: UTAU(:), UTAUPR(:)
      REAL, intent(out)  :: CMU(:), CWT(:), DTAUCP(:), EXPBEA(0:)
      REAL, intent(out)  :: FLYR(:), GL(0:,:), HLPR(0:), UMU(:)
      REAL, intent(out)  :: OPRIM(:), PMOM(0:,:), TAUCPR(0:) 

c     ..
c     .. Local Scalars ..

      REAL, PARAMETER :: ABSCUT = 10000.

      INTEGER   :: IQ, IU, K, LC, LU, I
      REAL      :: ABSTAU, F

      IF( .NOT. USRTAU ) THEN
c  ** Set output levels at computational layer boundaries
         NTAU  = NLYR + 1
         DO LC = 0, NTAU - 1
            UTAU(LC + 1) = TAUC(LC)
         ENDDO
      END IF
c                        ** Apply delta-M scaling and move description
c                           of computational layers to local variables
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
               GL(K,LC) = REAL(2*K + 1)*OPRIM(LC)*PMOM(K,LC)
            ENDDO

            F  = rZERO
         ELSE
c                                    ** Do delta-M transformation
            F  = PMOM(NSTR,LC)
            OPRIM(LC) = SSALB(LC) * (rONE - F) / ( rONE - F*SSALB(LC))
            DTAUCP(LC) = ( rONE - F*SSALB(LC))*DTAUC(LC)
            TAUCPR(LC) = TAUCPR(LC - 1) + DTAUCP(LC)

            DO K = 0, NSTR - 1
               GL(K,LC) = REAL(2*K + 1) * OPRIM(LC)
     $                    * (PMOM(K,LC) - F) / (rONE - F)
            ENDDO
         END IF

         FLYR(LC) = F
         EXPBEA(LC) = rZERO
      ENDDO
c 
* calculate slant optical depth
*              
         IF(umu0 <  rZERO) THEN
           IF(nid(0) < 0) THEN
             tausla(0) = largest
             tauslau(0) = largest
           ELSE
             sum = rZERO
             sumu = rZERO
             DO lc = 1, nid(0)
               sum = sum + 2.*dtaucp(lc)*dsdh(0,lc)
               sumu = sumu + 2.*dtauc(lc)*dsdh(0,lc)
             END DO
             tausla(0) = sum 
             tauslau(0) = sumu 
           END IF
         END IF

         expbea(0) = EXP( -tausla(0) )

*
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
               sum = sum + 2.*dtaucp(lu)*dsdh(lc,lu)
               sumu = sumu + 2.*dtauc(lu)*dsdh(lc,lu)
            ENDDO
            tausla(lc) = sum 
            tauslau(lc) = sumu 
            IF(tausla(lc) == tausla(lc-1)) THEN
              mu2(lc) = largest
            ELSE
              mu2(lc) = (taucpr(lc)-taucpr(lc-1))
     $                  /(tausla(lc)-tausla(lc-1))
              mu2(lc) = SIGN( MAX(ABS(mu2(lc)),1./largest),mu2(lc) )
            END IF
          END IF
          expbea(lc) = EXP( -tausla(lc) )
         ENDDO

c                      ** If no thermal emission, cut off medium below
c                         absorption optical depth = ABSCUT ( note that
c                         delta-M transformation leaves absorption
c                         optical depth invariant ).  Not worth the
c                         trouble for one-layer problems, though.

      LYRCUT = ABSTAU >= ABSCUT .AND. .NOT. PLANK .AND.
     $         IBCND /= 1 .AND. NLYR > 1

      IF( .NOT.LYRCUT ) NCUT = NLYR

c                             ** Set arrays defining location of user
c                             ** output levels within delta-M-scaled
c                             ** computational mesh
      DO LU = 1, NTAU
         DO LC = 1, NLYR
            IF( UTAU(LU) >= TAUC(LC - 1 ) .AND.
     $          UTAU(LU) <= TAUC(LC) ) EXIT
         ENDDO
         LC = MIN( NLYR,LC )

         UTAUPR(LU) = UTAU(LU)
         IF( DELTAM ) THEN
           UTAUPR( LU ) = TAUCPR(LC - 1)
     $                  + (rONE - SSALB(LC)*FLYR(LC))
     $                    * (UTAU(LU) - TAUC(LC-1))
         ENDIF
         LAYRU(LU) = LC
      ENDDO
c                      ** Calculate computational polar angle cosines
c                         and associated quadrature weights for Gaussian
c                         quadrature on the interval (0,1) (upward)
      NN   = NSTR / 2

      CALL QGAUSN( NN, CMU, CWT )
c                                  ** Downward (neg) angles and weights
      DO IQ = 1, NN
         CMU(IQ + NN) = -CMU(IQ)
         CWT(IQ + NN) = CWT(IQ)
      ENDDO


         DO IQ = 1, NN
C                      ** Dither mu2 if it is close to one of the 
C                         quadrature angles.
           DO  lc = 1, nlyr
             IF (  ABS(mu2(lc)) < 1.E5 ) THEN
               IF( ABS(rONE - ABS(mu2(lc))/CMU(IQ)) < 0.05 ) 
     $           mu2(lc) = mu2(lc)*0.999
             ENDIF
           END DO
         END DO

      IF( .NOT. USRANG .OR. ( ONLYFL .AND. MXUMU >= NSTR ) ) THEN
c                                   ** Set output polar angles to
c                                      computational polar angles
         NUMU = NSTR
         DO IU = 1, NN
            UMU(IU) = -CMU(NN + 1 - IU)
         ENDDO

         DO IU = NN + 1, NSTR
            UMU(IU) = CMU(IU - NN)
         ENDDO
      END IF


      IF( USRANG .AND. IBCND == 1 ) THEN
c                               ** Shift positive user angle cosines to
c                                  upper locations and put negatives
c                                  in lower locations
         DO IU = 1, NUMU
            UMU(IU+NUMU) = UMU(IU)
         ENDDO

         DO IU = 1, NUMU
            UMU(IU) = -UMU(2*NUMU + 1 - IU)
         ENDDO

         NUMU   = 2*NUMU
      END IF


      IF( .NOT. LYRCUT .AND. .NOT.LAMBER ) THEN
         DO K = 0, NSTR
            HLPR(K) = REAL((2*K + 1))*HL(K)
         ENDDO
      END IF

      END SUBROUTINE SETDIS

      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                   LAMBER, LYRCUT, NCOL, NCUT,
     &                   NN, NSTR, TAUCPR, WK )

c        Calculate coefficient matrix for the set of equations
c        obtained from the boundary conditions and the continuity-
c        of-intensity-at-layer-interface equations;  store in the
c        special banded-matrix format required by LINPACK routines

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       DELM0    :  Kronecker delta, delta-sub-m0
c       GC       :  Eigenvectors at polar quadrature angles, SC(1)
c       KK       :  Eigenvalues of coeff. matrix in Eq. SS(7)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       NN       :  Number of streams in a hemisphere (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                      scaled by Eq. SC(12); in banded form required
c                      by LINPACK solution routines
c       NCOL     :  Counts of columns in CBAND

c   I N T E R N A L    V A R I A B L E S:

c       IROW     :  Points to row in CBAND
c       JCOL     :  Points to position in layer block
c       LDA      :  Row dimension of CBAND
c       NCD      :  Number of diagonals below or above main diagonal
c       NSHIFT   :  For positioning number of rows in band storage
c       WK       :  Temporary storage for EXP evaluations

c   Called by- DISORT, ALBTRN
c +--------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL, intent(in)  :: LAMBER, LYRCUT
      INTEGER, intent(in)  :: NCUT, NN, NSTR
      INTEGER, intent(out) :: NCOL
      REAL, intent(in)     :: DELM0
c     ..
c     .. Array Arguments ..

      REAL, intent(in)  ::  BDR(:,0:), CMU(:),
     $          CWT(:), DTAUCP(:), GC(:,:,:),
     $          KK(:,:), TAUCPR(0:)
      REAL, intent(inout) :: WK(:)
      REAL, intent(out)   :: CBAND(:,:)
c     ..
c     .. Local Scalars ..

      INTEGER   :: IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      REAL      :: EXPA, SUM
c     ..

      CBAND = rZERO

      NCD    = 3*NN - 1
      LDA    = 3*NCD + 1
      NSHIFT = LDA - 2*NSTR + 1
      NCOL   = 0
c                         ** Use continuity conditions of Eq. STWJ(17)
c                            to form coefficient matrix in STWJ(20);
c                            employ scaling transformation STWJ(22)
      DO 60 LC = 1, NCUT

         DO 10 IQ = 1, NN
            WK( IQ ) = EXP( KK( IQ,LC )*DTAUCP( LC ) )
   10    CONTINUE

         JCOL  = 0

         DO 30 IQ = 1, NN

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 20 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )*WK( IQ )
               IROW  = IROW + 1
   20       CONTINUE

            JCOL  = JCOL + 1

   30    CONTINUE


         DO 50 IQ = NN + 1, NSTR

            NCOL  = NCOL + 1
            IROW  = NSHIFT - JCOL

            DO 40 JQ = 1, NSTR
               CBAND( IROW + NSTR, NCOL ) =   GC( JQ, IQ, LC )*
     &                                          WK( NSTR + 1 - IQ )
               CBAND( IROW, NCOL )        = - GC( JQ, IQ, LC )
               IROW  = IROW + 1
   40       CONTINUE

            JCOL  = JCOL + 1

   50    CONTINUE

   60 CONTINUE
c                  ** Use top boundary condition of STWJ(20a) for
c                     first layer

      JCOL  = 0

      DO 80 IQ = 1, NN

         EXPA  = EXP( KK( IQ,1 )*TAUCPR( 1 ) )
         IROW  = NSHIFT - JCOL + NN

         DO 70 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )*EXPA
            IROW  = IROW + 1
   70    CONTINUE

         JCOL  = JCOL + 1

   80 CONTINUE


      DO 100 IQ = NN + 1, NSTR

         IROW  = NSHIFT - JCOL + NN

         DO 90 JQ = NN, 1, -1
            CBAND( IROW, JCOL + 1 ) = GC( JQ, IQ, 1 )
            IROW  = IROW + 1
   90    CONTINUE

         JCOL  = JCOL + 1

  100 CONTINUE
c                           ** Use bottom boundary condition of
c                              STWJ(20c) for last layer

      NNCOL = NCOL - NSTR
      JCOL  = 0

      DO 130 IQ = 1, NN

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR

         DO 120 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

c                          ** No azimuthal-dependent intensity if Lam-
c                             bert surface; no intensity component if
c                             truncated bottom layer

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )

            ELSE

               SUM  = 0.0

               DO 110 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                     GC( NN + 1 - K, IQ, NCUT )
  110          CONTINUE

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT ) -
     &                                ( 1.+ DELM0 )*SUM
            END IF

            IROW  = IROW + 1

  120    CONTINUE

         JCOL  = JCOL + 1

  130 CONTINUE


      DO 160 IQ = NN + 1, NSTR

         NNCOL  = NNCOL + 1
         IROW   = NSHIFT - JCOL + NSTR
         EXPA   = WK( NSTR + 1 - IQ )

         DO 150 JQ = NN + 1, NSTR

            IF( LYRCUT .OR. ( LAMBER .AND. DELM0.EQ.0 ) ) THEN

               CBAND( IROW, NNCOL ) = GC( JQ, IQ, NCUT )*EXPA

            ELSE

               SUM  = 0.0

               DO 140 K = 1, NN
                  SUM  = SUM + CWT( K )*CMU( K )*BDR( JQ - NN, K )*
     &                         GC( NN + 1 - K, IQ, NCUT )
  140          CONTINUE

               CBAND( IROW, NNCOL ) = ( GC( JQ,IQ,NCUT ) -
     &                                ( 1.+ DELM0 )*SUM )*EXPA
            END IF

            IROW  = IROW + 1

  150    CONTINUE

         JCOL  = JCOL + 1

  160 CONTINUE

      END SUBROUTINE SETMTX

      SUBROUTINE SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL, MAZIM,
     $                   NN, NSTR, YLMC, CC, EVECC, EVAL, KK, GC,
     $                   AAD, EVECCD, EVALD, WKD )

c         Solves eigenvalue/vector problem necessary to construct
c         homogeneous part of discrete ordinate solution; STWJ(8b)
c         ** NOTE ** Eigenvalue problem is degenerate when single
c                    scattering albedo = 1;  present way of doing it
c                    seems numerically more stable than alternative
c                    methods that we tried

c   I N P U T     V A R I A B L E S:

c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2l+1 and single-scatter albedo)
c       CMU    :  Computational polar angle cosines
c       CWT    :  Weights for quadrature over polar angle cosine
c       MAZIM  :  Order of azimuthal component
c       NN     :  Half the total number of streams
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles CMU
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5); needed in SS(15&18)
c       EVAL   :  NN eigenvalues of Eq. SS(12) on return from ASYMTX
c                    but then square roots taken
c       EVECC  :  NN eigenvectors  (G+) - (G-)  on return
c                    from ASYMTX ( column j corresponds to EVAL(j) )
c                    but then  (G+) + (G-)  is calculated from SS(10),
c                    G+  and  G-  are separated, and  G+  is stacked on
c                    top of  G-  to form NSTR eigenvectors of SS(7)
c       GC     :  Permanent storage for all NSTR eigenvectors, but
c                    in an order corresponding to KK
c       KK     :  Permanent storage for all NSTR eigenvalues of SS(7),
c                    but re-ordered with negative values first ( square
c                    roots of EVAL taken and negatives added )

c   I N T E R N A L   V A R I A B L E S:

c       AMB,APB :  Matrices (alpha-beta), (alpha+beta) in reduced
c                    eigenvalue problem
c       ARRAY   :  Complete coefficient matrix of reduced eigenvalue
c                    problem: (alfa+beta)*(alfa-beta)
c       GPPLGM  :  (G+) + (G-) (cf. Eqs. SS(10-11))
c       GPMIGM  :  (G+) - (G-) (cf. Eqs. SS(10-11))
c       WKD     :  Scratch array required by ASYMTX

c   Called by- DISORT, ALBTRN
c   Calls- ASYMTX, ERRMSG
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER, intent(in) ::  MAZIM, NN, NSTR
c     ..
c     .. Array Arguments ..

      REAL, intent(out) ::  EVAL(:), KK(:)
      REAL, intent(out) ::  CC(:,:), EVECC(:,:), GC(:,:)
      REAL, intent(out) ::  AMB(:,:), APB(:,:), ARRAY(:,:)
      REAL, intent(in)  ::  CMU(:), CWT(:), GL(0:), YLMC(0:,:)
      REAL(8), intent(inout)  :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)

c     ..
c     .. Local Scalars ..

      INTEGER :: IER, IQ, JQ, KQ, L
      REAL    :: ALPHA, BETA, GPMIGM, GPPLGM, SUM

c                             ** Calculate quantities in Eqs. SS(5-6)
      DO 40 IQ = 1, NN

         DO 20 JQ = 1, NSTR

            SUM  = 0.0
            DO 10 L = MAZIM, NSTR - 1
               SUM  = SUM + GL( L )*YLMC( L, IQ )*YLMC( L, JQ )
   10       CONTINUE

            CC( IQ, JQ ) = 0.5*SUM*CWT( JQ )

   20    CONTINUE

         DO 30 JQ = 1, NN
c                             ** Fill remainder of array using symmetry
c                                relations  C(-mui,muj) = C(mui,-muj)
c                                and        C(-mui,-muj) = C(mui,muj)

            CC( IQ + NN, JQ ) = CC( IQ, JQ + NN )
            CC( IQ + NN, JQ + NN ) = CC( IQ, JQ )

c                                       ** Get factors of coeff. matrix
c                                          of reduced eigenvalue problem

            ALPHA  = CC( IQ, JQ ) / CMU( IQ )
            BETA   = CC( IQ, JQ + NN ) / CMU( IQ )
            AMB( IQ, JQ ) = ALPHA - BETA
            APB( IQ, JQ ) = ALPHA + BETA

   30    CONTINUE

         AMB( IQ, IQ ) = AMB( IQ, IQ ) - rONE / CMU( IQ )
         APB( IQ, IQ ) = APB( IQ, IQ ) - rONE / CMU( IQ )

   40 CONTINUE
c                      ** Finish calculation of coefficient matrix of
c                         reduced eigenvalue problem:  get matrix
c                         product (alfa+beta)*(alfa-beta); SS(12)
      DO 70 IQ = 1, NN

         DO 60 JQ = 1, NN

            SUM  = 0.
            DO 50 KQ = 1, NN
               SUM  = SUM + APB( IQ, KQ )*AMB( KQ, JQ )
   50       CONTINUE

            ARRAY( IQ, JQ ) = SUM

   60    CONTINUE

   70 CONTINUE
c                      ** Find (real) eigenvalues and eigenvectors

      CALL ASYMTX( ARRAY, EVECC, EVAL, IER, WKD, AAD,
     $             EVECCD, EVALD )

      IF( IER.GT.0 ) THEN

         WRITE( *, FMT = '(//,A,I4,A)' ) ' ASYMTX--eigenvalue no. ',
     &      IER, '  didnt converge.  Lower-numbered eigenvalues wrong.'

         CALL ERRMSG( 'ASYMTX--convergence problems',.True.)

      END IF

      DO 80 IQ = 1, NN
         EVAL( IQ )    = SQRT( ABS( EVAL( IQ ) ) )
         KK( IQ + NN ) = EVAL( IQ )
c                                      ** Add negative eigenvalue
         KK( NN + 1 - IQ ) = -EVAL( IQ )
   80 CONTINUE

c                          ** Find eigenvectors (G+) + (G-) from SS(10)
c                             and store temporarily in APB array
      DO 110 JQ = 1, NN

         DO 100 IQ = 1, NN

            SUM  = 0.
            DO 90 KQ = 1, NN
               SUM  = SUM + AMB( IQ, KQ )*EVECC( KQ, JQ )
   90       CONTINUE

            APB( IQ, JQ ) = SUM / EVAL( JQ )

  100    CONTINUE

  110 CONTINUE


      DO 130 JQ = 1, NN
CDIR$ IVDEP
         DO 120 IQ = 1, NN

            GPPLGM = APB( IQ, JQ )
            GPMIGM = EVECC( IQ, JQ )
c                                ** Recover eigenvectors G+,G- from
c                                   their sum and difference; stack them
c                                   to get eigenvectors of full system
c                                   SS(7) (JQ = eigenvector number)

            EVECC( IQ,      JQ ) = 0.5*( GPPLGM + GPMIGM )
            EVECC( IQ + NN, JQ ) = 0.5*( GPPLGM - GPMIGM )

c                                ** Eigenvectors corresponding to
c                                   negative eigenvalues (corresp. to
c                                   reversing sign of 'k' in SS(10) )
            GPPLGM = - GPPLGM
            EVECC(IQ,   JQ+NN) = 0.5 * ( GPPLGM + GPMIGM )
            EVECC(IQ+NN,JQ+NN) = 0.5 * ( GPPLGM - GPMIGM )
            GC( IQ+NN,   JQ+NN )   = EVECC( IQ,    JQ )
            GC( NN+1-IQ, JQ+NN )   = EVECC( IQ+NN, JQ )
            GC( IQ+NN,   NN+1-JQ ) = EVECC( IQ,    JQ+NN )
            GC( NN+1-IQ, NN+1-JQ ) = EVECC( IQ+NN, JQ+NN )

  120    CONTINUE

  130 CONTINUE

      END SUBROUTINE SOLEIG

      SUBROUTINE SOLVE0( B, BDR, BEM, BPLANK, CBAND, CMU, CWT, EXPBEA,
     $                   FBEAM, FISOT, IPVT, LAMBER, LL, LYRCUT, MAZIM,
     $                   NCOL, NCUT, NN, NSTR,
     $                   PI, TPLANK, TAUCPR, UMU0, Z, ZZ, ZPLK0, ZPLK1 )

      use linalgebra, only : linalgebra_t

c        Construct right-hand side vector B for general boundary
c        conditions STWJ(17) and solve system of equations obtained
c        from the boundary conditions and the continuity-of-
c        intensity-at-layer-interface equations.
c        Thermal emission contributes only in azimuthal independence.

c     I N P U T      V A R I A B L E S:

c       BDR      :  Surface bidirectional reflectivity
c       BEM      :  Surface bidirectional emissivity
c       BPLANK   :  Bottom boundary thermal emission
c       CBAND    :  Left-hand side matrix of linear system Eq. SC(5),
c                   scaled by Eq. SC(12); in banded form required
c                   by LINPACK solution routines
c       CMU      :  Abscissae for Gauss quadrature over angle cosine
c       CWT      :  Weights for Gauss quadrature over angle cosine
c       EXPBEA   :  Transmission of incident beam, EXP(-TAUCPR/UMU0)
c       LYRCUT   :  Logical flag for truncation of comput. layer
c       MAZIM    :  Order of azimuthal component
c       ncol     :  Counts of columns in CBAND
c       NN       :  Order of double-Gauss quadrature (NSTR/2)
c       NCUT     :  Total number of computational layers considered
c       TPLANK   :  Top boundary thermal emission
c       TAUCPR   :  Cumulative optical depth (delta-M-scaled)
c       ZZ       :  Beam source vectors in Eq. SS(19)
c       ZPLK0    :  Thermal source vectors Z0, by solving Eq. SS(16)
c       ZPLK1    :  Thermal source vectors Z1, by solving Eq. SS(16)
c       (remainder are DISORT input variables)

c   O U T P U T     V A R I A B L E S:

c       B        :  Right-hand side vector of Eq. SC(5) going into
c                   SGBSL; returns as solution vector of Eq. SC(12),
c                   constants of integration without exponential term
c
c      LL        :  Permanent storage for B, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       IPVT     :  Integer vector of pivot indices
c       IT       :  Pointer for position in  B
c       NCD      :  Number of diagonals below or above main diagonal
c       RCOND    :  Indicator of singularity for CBAND
c       Z        :  Scratch array required by SGBCO

c   Called by- DISORT
c   Calls- SGBCO, ERRMSG, SGBSL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL, intent(in) :: LAMBER, LYRCUT
      INTEGER, intent(in) :: MAZIM, NCOL, NCUT, NN, NSTR
      REAL,    intent(in) :: BPLANK, FBEAM, FISOT, PI, TPLANK, UMU0
c     ..
c     .. Array Arguments ..

      INTEGER, intent(inout) :: IPVT(:)
      REAL, intent(out)      :: B(:), LL(:,:)
      REAL, intent(inout)    :: Z(:)
      REAL, intent(inout)    :: CBAND(:,:)
      REAL, intent(in)       :: BDR(:,0:), BEM(:),
     $          CMU(:), CWT(:),
     $          EXPBEA(0:), TAUCPR(0:),
     $          ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)
c     ..
c     .. Local Scalars ..

      INTEGER ::   IPNT, IQ, IT, JQ, LC, NCD
      REAL    :: RCOND, SUM
      TYPE(linalgebra_t) :: linpack

      B = rZERO
c                              ** Construct B,  STWJ(20a,c) for
c                                 parallel beam + bottom reflection +
c                                 thermal emission at top and/or bottom

      IF( MAZIM.GT.0 .AND. FBEAM.GT.0.0 ) THEN

c                                         ** Azimuth-dependent case
c                                            (never called if FBEAM = 0)
         IF( LYRCUT .OR. LAMBER ) THEN

c               ** No azimuthal-dependent intensity for Lambert surface;
c                  no intensity component for truncated bottom layer

            DO 10 IQ = 1, NN
c                                                  ** Top boundary
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )
c                                                  ** Bottom boundary

               B( NCOL - NN + IQ ) = -ZZ( IQ + NN, NCUT )*EXPBEA( NCUT )

   10       CONTINUE


         ELSE

            DO 30 IQ = 1, NN

               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 )

               SUM  = 0.
               DO 20 JQ = 1, NN
                  SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                         ZZ( NN + 1 - JQ, NCUT )*EXPBEA( NCUT )
   20          CONTINUE

               B( NCOL - NN + IQ ) = SUM
               IF( FBEAM.GT.0.0 ) B( NCOL - NN + IQ ) = SUM +
     &             ( BDR( IQ,0 )*UMU0*FBEAM / PI - ZZ( IQ + NN,NCUT ) )*
     &             EXPBEA( NCUT )

   30       CONTINUE

         END IF
c                             ** Continuity condition for layer
c                                interfaces of Eq. STWJ(20b)
         IT   = NN

         DO 50 LC = 1, NCUT - 1

            DO 40 IQ = 1, NSTR
               IT   = IT + 1
               B( IT ) = ( ZZ( IQ, LC+1 ) - ZZ( IQ, LC ) )*EXPBEA( LC )
   40       CONTINUE

   50    CONTINUE


      ELSE
c                                   ** Azimuth-independent case

         IF( FBEAM.EQ.0.0 ) THEN

            DO 60 IQ = 1, NN
c                                      ** Top boundary

               B( IQ ) = -ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK

   60       CONTINUE


            IF( LYRCUT ) THEN
c                               ** No intensity component for truncated
c                                  bottom layer
               DO 70 IQ = 1, NN
c                                      ** Bottom boundary

                  B( NCOL - NN + IQ ) = - ZPLK0( IQ + NN, NCUT ) -
     &                                    ZPLK1( IQ + NN, NCUT )*
     &                                    TAUCPR( NCUT )
   70          CONTINUE


            ELSE

               DO 90 IQ = 1, NN

                  SUM  = 0.
                  DO 80 JQ = 1, NN
                     SUM  = SUM + CWT( JQ )*CMU( JQ )*BDR( IQ, JQ )*
     &                            ( ZPLK0( NN + 1 - JQ,NCUT ) +
     &                        ZPLK1( NN + 1 - JQ,NCUT )*TAUCPR( NCUT ) )
   80             CONTINUE

                  B( NCOL - NN + IQ ) = 2.*SUM + BEM( IQ )*BPLANK -
     &                                  ZPLK0( IQ + NN, NCUT ) -
     &                                  ZPLK1( IQ + NN, NCUT )*
     &                                  TAUCPR( NCUT )
   90          CONTINUE

            END IF
c                             ** Continuity condition for layer
c                                interfaces, STWJ(20b)
            IT   = NN
            DO 110 LC = 1, NCUT - 1

               DO 100 IQ = 1, NSTR
                  IT   = IT + 1
                  B( IT ) =   ZPLK0( IQ, LC + 1 ) - ZPLK0( IQ, LC ) +
     &                      ( ZPLK1( IQ, LC + 1 ) - ZPLK1( IQ, LC ) )*
     &                      TAUCPR( LC )
  100          CONTINUE

  110       CONTINUE


         ELSE

            DO 120 IQ = 1, NN
               B( IQ ) = - ZZ( NN + 1 - IQ, 1 )*EXPBEA( 0 ) -
     &                   ZPLK0( NN + 1 - IQ, 1 ) + FISOT + TPLANK
  120       CONTINUE

            IF( LYRCUT ) THEN

               DO 130 IQ = 1, NN
                  B(NCOL-NN+IQ) = - ZZ(IQ+NN, NCUT) * EXPBEA(NCUT)
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  130          CONTINUE


            ELSE

               DO 150 IQ = 1, NN

                  SUM  = 0.
                  DO 140 JQ = 1, NN
                     SUM = SUM + CWT(JQ) * CMU(JQ) * BDR(IQ,JQ)
     &                          * ( ZZ(NN+1-JQ, NCUT) * EXPBEA(NCUT)
     &                            + ZPLK0(NN+1-JQ, NCUT)
     &                            + ZPLK1(NN+1-JQ, NCUT) * TAUCPR(NCUT))
  140             CONTINUE

                  B(NCOL-NN+IQ) = 2.*SUM + ( BDR(IQ,0) * UMU0*FBEAM/PI
     &                                - ZZ(IQ+NN, NCUT) ) * EXPBEA(NCUT)
     &                            + BEM(IQ) * BPLANK
     &                            - ZPLK0(IQ+NN, NCUT)
     &                            - ZPLK1(IQ+NN, NCUT) * TAUCPR(NCUT)
  150          CONTINUE

            END IF


            IT   = NN

            DO 170 LC = 1, NCUT - 1

               DO 160 IQ = 1, NSTR

                  IT   = IT + 1
                  B(IT) = ( ZZ(IQ,LC+1) - ZZ(IQ,LC) ) * EXPBEA(LC)
     &                    + ZPLK0(IQ,LC+1) - ZPLK0(IQ,LC) +
     &                    ( ZPLK1(IQ,LC+1) - ZPLK1(IQ,LC) ) * TAUCPR(LC)
  160          CONTINUE

  170       CONTINUE

         END IF

      END IF
c                     ** Find L-U (lower/upper triangular) decomposition
c                        of band matrix CBAND and test if it is nearly
c                        singular (note: CBAND is destroyed)
c                        (CBAND is in LINPACK packed format)
      RCOND  = rZERO
      NCD    = 3*NN - 1

      CALL linpack%SGBCO( CBAND, NCOL, NCD, NCD, IPVT, RCOND, Z )

      IF( rONE + RCOND == rONE )
     $    CALL ERRMSG('SOLVE0--SGBCO says matrix near singular',.FALSE.)

c                   ** Solve linear system with coeff matrix CBAND
c                      and R.H. side(s) B after CBAND has been L-U
c                      decomposed.  Solution is returned in B.

      CALL linpack%SGBSL( CBAND, NCOL, NCD, NCD, IPVT, B, 0 )

c                   ** Zero CBAND (it may contain 'foreign'
c                      elements upon returning from LINPACK);
c                      necessary to prevent errors

      CBAND = rZERO

      DO 190 LC = 1, NCUT

         IPNT  = LC*NSTR - NN

         DO 180 IQ = 1, NN
            LL( NN + 1 - IQ, LC ) = B( IPNT + 1 - IQ )
            LL( IQ + NN,     LC ) = B( IQ + IPNT )
  180    CONTINUE

  190 CONTINUE

      END SUBROUTINE SOLVE0

      SUBROUTINE SURFAC( ALBEDO, DELM0, FBEAM, HLPR, LAMBER, MAZIM,
     $                   NN, NUMU, NSTR, ONLYFL, UMU,
     $                   USRANG, YLM0, YLMC, YLMU, BDR, EMU, BEM, RMU )

c       Specifies user's surface bidirectional properties, STWJ(21)

c   I N P U T     V A R I A B L E S:

c       DELM0  :  Kronecker delta, delta-sub-m0
c       HLPR   :  Legendre moments of surface bidirectional reflectivity
c                    (with 2K+1 factor included)
c       MAZIM  :  Order of azimuthal component
c       NN     :  Order of double-Gauss quadrature (NSTR/2)
c       YLM0   :  Normalized associated Legendre polynomial
c                 at the beam angle
c       YLMC   :  Normalized associated Legendre polynomials
c                 at the quadrature angles
c       YLMU   :  Normalized associated Legendre polynomials
c                 at the user angles
c       (remainder are DISORT input variables)

c    O U T P U T     V A R I A B L E S:

c       BDR :  Surface bidirectional reflectivity (computational angles)
c       RMU :  Surface bidirectional reflectivity (user angles)
c       BEM :  Surface directional emissivity (computational angles)
c       EMU :  Surface directional emissivity (user angles)

c    I N T E R N A L     V A R I A B L E S:

c       DREF      Directional reflectivity
c       NMUG   :  Number of angle cosine quadrature points on (0,1) for
c                   integrating bidirectional reflectivity to get
c                   directional emissivity (it is necessary to use a
c                   quadrature set distinct from the computational
c                   angles, because the computational angles may not be
c                   dense enough--NSTR may be too small--to give an
c                   accurate approximation for the integration).
c       GMU    :  The NMUG angle cosine quadrature points on (0,1)
c       GWT    :  The NMUG angle cosine quadrature weights on (0,1)
c       YLMG   :  Normalized associated Legendre polynomials
c                   at the NMUG quadrature angles

c   Called by- DISORT
c   Calls- QGAUSN, LEPOLY, ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER, PARAMETER ::   NMUG = 10, MAXSTR = 100
c     ..
c     .. Scalar Arguments ..

      LOGICAL, intent(in) :: LAMBER, ONLYFL, USRANG
      INTEGER, intent(in) :: MAZIM, NN, NSTR, NUMU
      REAL, intent(in)    :: ALBEDO, DELM0, FBEAM
c     ..
c     .. Array Arguments ..

      REAL, intent(in)  ::  HLPR(0:), UMU(:),
     $                      YLM0(0:), YLMC(0:,:), YLMU(0:,:)
      REAL, intent(out) ::  BDR(:,0:), BEM(:), EMU(:), RMU(:,0:)
c     ..
c     .. Local Scalars ..

      LOGICAL :: PASS1
      INTEGER :: IQ, IU, JG, JQ, K
      REAL    :: DREF, SGN, SUM
c     ..
c     .. Local Arrays ..

      REAL    :: GMU(NMUG), GWT(NMUG), YLMG(0:MAXSTR,NMUG)

      SAVE      PASS1, GMU, GWT, YLMG
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN
         PASS1  = .FALSE.

         CALL QGAUSN( NMUG, GMU, GWT )
         CALL LEPOLY( NMUG, 0, MAXSTR, GMU, YLMG )

c                       ** Convert Legendre polys. to negative GMU
         SGN  = - 1.0

         DO K = 0, MAXSTR
            SGN  = - SGN
            DO JG = 1, NMUG
               YLMG(K,JG) = SGN*YLMG(K,JG)
            ENDDO
         ENDDO
      END IF


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
c                                  ** Compute surface bidirectional
c                                     properties at computational angles
         DO IQ = 1, NN
            DO JQ = 1, NN
               SUM  = rZERO
               DO K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLMC(K,JQ+NN)
               ENDDO
               BDR(IQ,JQ) = (2.- DELM0)*SUM
            ENDDO

            IF( FBEAM > rZERO ) THEN
               SUM  = rZERO
               DO K = MAZIM, NSTR - 1
                  SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLM0(K)
               ENDDO
               BDR(IQ,0) = (2.- DELM0)*SUM
            END IF
         ENDDO

         IF( MAZIM == 0 ) THEN
            IF( NSTR > MAXSTR )
     $          CALL ERRMSG('SURFAC--parameter MAXSTR too small',.True.)
c                              ** Integrate bidirectional reflectivity
c                                 at reflection polar angles CMU and
c                                 incident angles GMU to get
c                                 directional emissivity at
c                                 computational angles CMU.
            DO IQ = 1, NN
               DREF  = rZERO
               DO JG = 1, NMUG
                  SUM  = rZERO
                  DO K = 0, NSTR - 1
                     SUM  = SUM + HLPR(K)*YLMC(K,IQ)*YLMG(K,JG)
                  ENDDO
                  DREF  = DREF + 2.*GWT(JG)*GMU(JG)*SUM
               ENDDO
               BEM(IQ) = rONE - DREF
            ENDDO
         END IF
      END IF
c                                       ** Compute surface bidirectional
c                                          properties at user angles

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
                     RMU(IU,IQ) = (2. - DELM0)*SUM
                  ENDDO


                  IF( FBEAM > rZERO ) THEN
                     SUM  = rZERO
                     DO K = MAZIM, NSTR - 1
                        SUM  = SUM + HLPR( K )*YLMU( K, IU )*YLM0( K )
                     ENDDO
                     RMU(IU,0) = (2. - DELM0)*SUM
                  END IF

                  IF( MAZIM == 0 ) THEN
c                               ** Integrate bidirectional reflectivity
c                                  at reflection angles UMU and
c                                  incident angles GMU to get
c                                  directional emissivity at
c                                  user angles UMU.
                     DREF  = rZERO
                     DO JG = 1, NMUG
                        SUM  = rZERO
                        DO K = 0, NSTR - 1
                           SUM = SUM + HLPR(K)*YLMU(K,IU)*YLMG(K,JG)
                        ENDDO
                        DREF  = DREF + 2.*GWT( JG )*GMU( JG )*SUM
                     ENDDO
                     EMU( IU ) = rONE - DREF
                  END IF
               END IF
            END IF
         ENDDO
      END IF

      END SUBROUTINE SURFAC

      SUBROUTINE SOLVEC( AMB, APB, ARRAY, CMU, CWT, GL,
     &     MAZIM, NN, NSTR, YLM0, YLMC, CC, 
     &     EVECC, EVAL, KK, GC, AAD, EVECCD, EVALD,
     &     WK, WKD, DELM0, FBEAM, IPVT, PI, UMU0, ZJ, ZZ,
     &     OPRIM, LC, mu2, glsave, dgl)

*bm  SOLVEC calls SOLEIG and UPBEAM; if UPBEAM reports a potenially 
*bm  unstable solution, the calculation is repeated with a slightly 
*bm  changed single scattering albedo; this process is iterates 
*bm  until a stable solution is found; as stable solutions may be 
*bm  reached either by increasing or by decreasing the single 
*bm  scattering albedo, both directions are explored ('upward' and
*bm  'downward' iteration); the solution which required the smaller 
*bm  change in the single scattering albedo is finally returned 
*bm  by SOLVEC.

cgy added glsave and dgl to call to allow adjustable dimensioning

c     .. Scalar Arguments ..

      INTEGER, intent(in) :: MAZIM, NN, NSTR, LC
      REAL, intent(in)    :: DELM0, FBEAM, PI, UMU0, OPRIM
      REAL, intent(in)    :: mu2

c     ..
c     .. Array Arguments ..

      INTEGER, intent(inout) ::   IPVT(:)
      
      REAL, intent(inout) :: WK(:)
      REAL, intent(in) ::  
     $     CMU(:), CWT(:), YLM0(0:), YLMC(0:,:)

      REAL, intent(inout) :: AMB(:,:), APB(:,:), ARRAY(:,:)
      REAL, intent(inout) :: EVAL(:), EVECC(:,:), KK(:), GC(:,:)
      REAL, intent(inout) :: GLSAVE(0:), DGL(0:), CC(:,:), ZJ(:), ZZ(:)
      REAL, intent(out)   :: GL(0:)

      INTEGER :: K
      REAL(8) :: AAD(:,:), EVALD(:), EVECCD(:,:), WKD(:)

*bm   Variables for instability fix
      
      INTEGER :: UAGAIN, DAGAIN
      REAL    :: MINRCOND, ADD, UADD, DADD, SSA, DSSA, FACTOR
      
      LOGICAL ::  DONE, NOUP, NODN, DEBUG, INSTAB
      
*bm   reset parameters

      DONE = .FALSE.
      NOUP = .FALSE.
      NODN = .FALSE.


*bm   flag for printing debugging output      
*      DEBUG  = .TRUE.
      DEBUG  = .FALSE.

*bm   instability parameter; the solution is considered 
*bm   unstable, if the RCOND reported by SGECO is smaller 
*bm   than MINRCOND
      MINRCOND = 5000. * R1MACH(4)

*bm   if an instability is detected, the single scattering albedo
*bm   is iterated downwards in steps of DADD and upwards in steps 
*bm   of UADD; in practice, MINRCOND and -MINRCOND should 
*bm   be reasonable choices for these parameters
      DADD    = -MINRCOND
      UADD    = MINRCOND

      UAGAIN = 0
      DAGAIN = 0
      ADD   = DADD
      

*bm   save array GL( ) because it will be 
*bm   changed if an iteration should be neccessary
      DO K = MAZIM, NSTR - 1
         GLSAVE( K ) =  GL( K )
      ENDDO
      
      SSA = OPRIM


*bm   in case of an instability reported by UPBEAM (INSTAB)
*bm   the single scattering albedo will be changed by a small 
*bm   amount (ADD); this is indicated by DAGAIN or UAGAIN 
*bm   being larger than 0; a change in the single scattering 
*bm   albedo is equivalent to scaling the array GL( )

 666  IF ( DAGAIN > 0 .OR. UAGAIN > 0)  THEN
         FACTOR = (SSA + ADD) / SSA
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GL( K ) * FACTOR
         ENDDO

         SSA = SSA + ADD
         
*bm   if the single scattering albedo is now smaller than 0
*bm   the downward iteration is stopped and upward iteration 
*bm   is forced instead

         IF( SSA < DITHER) THEN
            NODN = .TRUE.
            DAGAIN = -1
            goto 778
         ENDIF

*bm   if the single scattering albedo is now larger than its maximum 
*bm   allowed value (1.0 - DITHER), the upward iteration is 
*bm   stopped and downward iteration is forced instead

         IF( SSA > rONE - DITHER) THEN
            NOUP = .TRUE.
            UAGAIN = -1
            goto 888
         ENDIF
      ENDIF


c     ** Solve eigenfunction problem in Eq. STWJ(8B);
c        return eigenvalues and eigenvectors

 777     CALL SOLEIG( AMB, APB, ARRAY, CMU, CWT, GL,
     &     MAZIM, NN, NSTR, YLMC, CC, EVECC, EVAL,
     &     KK, GC, AAD, EVECCD, EVALD, WKD )

c     ** Calculate particular solutions of
c        q.SS(18) for incident beam source

      IF ( FBEAM > rZERO ) THEN
         CALL  UPBEAM( mu2,
     $        ARRAY, CC, CMU, DELM0, FBEAM, GL,
     $        IPVT, MAZIM, NN, NSTR, PI, UMU0, WK,
     $        YLM0, YLMC, ZJ, ZZ, MINRCOND, INSTAB)
      ENDIF
      
c     ** Calculate particular solutions of
c        Eq. SS(15) for thermal emission source
c        (not available in psndo.f)
      
*bm   finished if the result is stable on the first try
      IF ( (.NOT. INSTAB) .AND. 
     $     (UAGAIN == 0) .AND. (DAGAIN == 0)) THEN
         goto 999
      ENDIF

*bm   downward iteration
      IF( INSTAB .AND. UAGAIN == 0 )  THEN
         DAGAIN = DAGAIN + 1
         GOTO 666
      ENDIF
      
*bm   upward iteration
      IF( INSTAB .AND. UAGAIN > 0 )  THEN
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF


*bm   ( DAGAIN .NE. 0 ) at this place means that the downward
*bm   iteration is finished 

 778  IF (DAGAIN /= 0 .AND. UAGAIN == 0) THEN
         
*bm   save downward iteration data for later use and 
*bm   restore original input data
         DO K = MAZIM, NSTR - 1
            DGL( K ) =  GL( K )
            GL( K ) =  GLSAVE( K )
         ENDDO

         DSSA = SSA
         SSA = OPRIM

*bm   start upward iteration
         ADD = UADD
         UAGAIN = UAGAIN + 1
         GOTO 666
      ENDIF

*bm   both iterations finished
 888  IF (DONE) THEN
         goto 998
      ENDIF


*bm  if neither upward nor downward iteration converged, the 
*bm  original conditions are restored and SOLEIG/UPBEAM 
*bm  is called for the last time 
         
      IF (NOUP .AND. NODN) THEN
         
         DO K = MAZIM, NSTR - 1
            GL( K ) =  GLSAVE( K )
         ENDDO
         
         SSA = OPRIM
         
         IF (DEBUG) THEN
            write (*,*) '! *** Neither upward nor downward iteration'
            write (*,*) '! *** converged; using original result.'
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if upward iteration did not converge, the stable downward conditions
*bm  are restored and SOLEIG/UPBEAM is called for the last time
      IF (NOUP) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** The upward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ENDIF

*bm  if downward iteration did not converge, we are done 
*bm  (the result of the upward iteration will be used)
      IF (NODN) THEN
         IF (DEBUG) THEN
            write (*,*) '! *** The downward iteration did not converge.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $           ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF
         
         DONE = .TRUE.
         GOTO 998
      ENDIF

      
*bm   if both iterations converged, and if the upward iteration 
*bm   required more steps than the downward iteration, the stable 
*bm   downward conditions are restored and SOLEIG/UPBEAM is 
*bm   called for the last time 
         
      IF (UAGAIN > DAGAIN) THEN
         DO K = MAZIM, NSTR - 1
            GL( K ) =  DGL( K )
         ENDDO
         
         SSA = DSSA
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using downward.'
            write (*,*) '! *** Had to iterate ', DAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         GOTO 777
      ELSE
         
         IF (DEBUG) THEN
            write (*,*) '! *** Both iterations converged;',
     $           ' using upward.'
            write (*,*) '! *** Had to iterate ', UAGAIN,
     $        ' times in layer LC =', LC,';'
            write (*,*) '! *** changed SSA from ',
     $           OPRIM, ' to ', SSA,','
            write (*,*) '! *** by a factor of ', SSA/OPRIM
         ENDIF

         DONE = .TRUE.
         goto 998
      ENDIF
      
*bm   finally restore original input data
 998  DO K = MAZIM, NSTR - 1
         GL( K ) =  GLSAVE( K )
      ENDDO
      
 999  CONTINUE

      END SUBROUTINE SOLVEC

      SUBROUTINE UPBEAM( mu2,
     $                   ARRAY, CC, CMU, DELM0, FBEAM, GL, IPVT, MAZIM,
     $                   NN, NSTR, PI, UMU0, WK, YLM0, YLMC, ZJ,
     $                   ZZ, MINRCOND, INSTAB )

      use linalgebra, only : linalgebra_t

c         Finds the incident-beam particular solution of SS(18)

c   I N P U T    V A R I A B L E S:

c       CC     :  C-sub-ij in Eq. SS(5)
c       CMU    :  Abscissae for Gauss quadrature over angle cosine
c       DELM0  :  Kronecker delta, delta-sub-m0
c       GL     :  Delta-M scaled Legendre coefficients of phase function
c                    (including factors 2L+1 and single-scatter albedo)
c       MAZIM  :  Order of azimuthal component
c       YLM0   :  Normalized associated Legendre polynomial
c                    at the beam angle
c       YLMC   :  Normalized associated Legendre polynomial
c                    at the quadrature angles
c       (remainder are DISORT input variables)

c   O U T P U T    V A R I A B L E S:

c       ZJ     :  Right-hand side vector X-sub-zero in SS(19); also the
c                 solution vector Z-sub-zero after solving that system

c       ZZ     :  Permanent storage for ZJ, but re-ordered

c   I N T E R N A L    V A R I A B L E S:

c       ARRAY  :  Coefficient matrix in left-hand side of Eq. SS(19)
c       IPVT   :  Integer vector of pivot indices required by LINPACK
c       WK     :  Scratch array required by LINPACK

c   Called by- DISORT
c   Calls- SGECO, ERRMSG, SGESL
c +-------------------------------------------------------------------+


c     .. Scalar Arguments ..

      INTEGER, intent(in)  :: MAZIM, NN, NSTR
      LOGICAL, intent(out) :: INSTAB
      REAL, intent(in)     :: MINRCOND
      REAL, intent(in)     :: DELM0, FBEAM, PI, UMU0
      REAL, intent(in)     :: mu2
c     ..
c     .. Array Arguments ..

      INTEGER, intent(inout) :: IPVT(:)

      REAL, intent(in) :: CC(:,:), CMU(:),
     $                    GL(0:), YLM0(0:), YLMC(0:,:)
      REAL, intent(inout) :: WK(:)
      REAL, intent(out) :: ARRAY(:,:)
      REAL, intent(out) :: ZJ(:), ZZ(:)
c     ..
c     .. Local Scalars ..

      INTEGER   :: IQ, JOB, JQ, K
      REAL      :: RCOND, SUM

      TYPE(linalgebra_t) :: linpack
c     ..
c     .. External Subroutines ..

c     EXTERNAL  ERRMSG, SGECO, SGESL
c     ..


      DO IQ = 1, NSTR
         DO JQ = 1, NSTR
            ARRAY( IQ, JQ ) = -CC( IQ, JQ )
         ENDDO

         ARRAY( IQ, IQ ) = 1.+ CMU( IQ ) / mu2 + ARRAY( IQ, IQ )

         SUM  = rZERO
         DO K = MAZIM, NSTR - 1
            SUM  = SUM + GL( K )*YLMC( K, IQ )*YLM0( K )
         ENDDO

         ZJ( IQ ) = ( 2.- DELM0 )*FBEAM*SUM / ( 4.*PI )
      ENDDO

c                  ** Find L-U (lower/upper triangular) decomposition
c                     of ARRAY and see if it is nearly singular
c                     (NOTE:  ARRAY is destroyed)
      RCOND  = rZERO

      CALL linpack%SGECO( ARRAY, NSTR, IPVT, RCOND, WK )

*bm      IF( 1.0 + RCOND.EQ.1.0 )
*bm     &    CALL ERRMSG('UPBEAM--SGECO says matrix near singular',.FALSE.)
*bm
*bm   replaced original check of RCOND by the following:

      INSTAB = .FALSE.
      IF( ABS(RCOND) .LT. MINRCOND )  THEN
         INSTAB = .TRUE.
         RETURN
      ENDIF

c                ** Solve linear system with coeff matrix ARRAY
c                   (assumed already L-U decomposed) and R.H. side(s)
c                   ZJ;  return solution(s) in ZJ
      JOB  = 0

      CALL linpack%SGESL( ARRAY, NSTR, IPVT, ZJ, JOB )

CDIR$ IVDEP
      DO IQ = 1, NN
         ZZ( IQ + NN )     = ZJ( IQ )
         ZZ( NN + 1 - IQ ) = ZJ( IQ + NN )
      ENDDO

      END SUBROUTINE UPBEAM

      SUBROUTINE ZEROAL( EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,
     &                   CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     &                   HLPR, YLM0,
     &                   ARRAY, CC, EVECC,
     &                   GL, YLMC, YLMU,
     &                   KK, LL, ZZ, ZPLK0, ZPLK1,
     &                   GC, LAYRU, UTAUPR,
     &                   GU, Z0U, Z1U, ZBEAM,
     &                   EVAL, AMB, APB, IPVT, Z,
     &                   RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     &                   TRNMED, U0U, UU )

c         ZERO ARRAYS

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Array Arguments ..

      INTEGER, intent(out) ::   IPVT(:), LAYRU(:)
      REAL, intent(out)    ::
     $          AMB(:,:), APB(:,:), ARRAY(:,:), CC(:,:),
     &          CMU(:), CWT(:), DFDT(:), EVAL(:), EVECC(:,:),
     &          EXPBEA(:), FLUP(:), FLYR(:), GC(:,:,:), GL(:,:),
     &          GU(:,:,:), HLPR(:), KK(:,:), LL(:,:), OPRIM(:),
     &          PSI(:), RFLDIR(:), RFLDN(:), TAUCPR(:),
     &          TRNMED(:), U0U(:,:), UAVG(:), UTAUPR(:),
     $          UU(:,:,:),
     &          WK(:), XR0(:), XR1(:), YLM0(:), YLMC(:,:),
     &          YLMU(:,:), Z(:), Z0(:), Z0U(:,:), Z1(:), Z1U(:,:),
     &          ZBEAM(:,:), ZJ(:), ZPLK0(:,:), ZPLK1(:,:), ZZ(:,:)

c     ..

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

      REAL FUNCTION DREF( MU, HL, NSTR )

c        Exact flux albedo for given angle of incidence, given
c        a bidirectional reflectivity characterized by its
c        Legendre coefficients ( NOTE** these will only agree
c        with bottom-boundary albedos calculated by DISORT in
c        the limit as number of streams go to infinity, because
c        DISORT evaluates the integral 'CL' only approximately,
c        by quadrature, while this routine calculates it exactly.)

c  INPUT :   MU     Cosine of incidence angle
c            HL     Legendre coefficients of bidirectional reflectivity
c          NSTR     Number of elements of HL to consider

c  INTERNAL VARIABLES (P-sub-L is the L-th Legendre polynomial) :

c       CL      Integral from 0 to 1 of  MU * P-sub-L(MU)
c                   (vanishes for  L = 3, 5, 7, ... )
c       PL      P-sub-L
c       PLM1    P-sub-(L-1)
c       PLM2    P-sub-(L-2)

c   Called by- CHEKIN
c   Calls- ERRMSG
c +-------------------------------------------------------------------+

c     .. Parameters ..

      INTEGER, PARAMETER ::  MAXTRM = 100
c     ..
c     .. Scalar Arguments ..

      INTEGER, intent(in) :: NSTR
      REAL, intent(in)    :: MU
c     ..
c     .. Array Arguments ..

      REAL, intent(in)   :: HL(0:NSTR)
c     ..
c     .. Local Scalars ..

      LOGICAL   :: PASS1
      INTEGER   :: L
      REAL      :: CL, PL, PLM1, PLM2
c     ..
c     .. Local Arrays ..

      REAL      :: C( MAXTRM )
c     ..

      SAVE      PASS1, C
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN
         PASS1  = .FALSE.
         CL     = 0.125
         C(2)   = 10.*CL

         DO L = 4, MAXTRM, 2
            CL   = - CL*(L - 3) / (L + 2)
            C(L) = 2.*(2*L + 1)*CL
         ENDDO
      END IF


      IF( NSTR < 2 .OR. ABS(MU) > rONE )
     &    CALL ERRMSG( 'DREF--input argument error(s)',.True. )

      IF( NSTR > MAXTRM )
     &    CALL ERRMSG( 'DREF--parameter MAXTRM too small',.True. )


      DREF  = HL(0) - 2.*HL(1)*MU
      PLM2  = rONE
      PLM1  = - MU

      DO L = 2, NSTR - 1
c                                ** Legendre polynomial recurrence
         PL = ((2*L - 1)*(-MU)*PLM1 - (L-1)*PLM2) / L
         IF( MOD( L,2 ) == 0 ) DREF = DREF + C(L)*HL(L)*PL
         PLM2  = PLM1
         PLM1  = PL
      ENDDO

      IF( DREF < rZERO .OR. DREF > rONE )
     $    CALL ERRMSG( 'DREF--albedo value not in (0,1)',.False. )

      END FUNCTION DREF

      REAL FUNCTION RATIO( A, B )

c        Calculate ratio  A/B  with over- and under-flow protection
c        (thanks to Prof. Jeff Dozier for some suggestions here).
c        Since this routine takes two logs, it is no speed demon,
c        but it is invaluable for comparing results from two runs
c        of a program under development.

c        NOTE:  In Fortran90, built-in functions TINY and HUGE
c               can replace the R1MACH calls.
c ---------------------------------------------------------------

c     .. Scalar Arguments ..

      REAL, intent(in) :: A, B
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      REAL      ABSA, ABSB, HUGE, POWA, POWB, POWMAX, POWMIN, TINY
c     ..
c     .. External Functions ..

c     EXTERNAL  R1MACH
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC ABS, LOG10, SIGN
c     ..
      SAVE      PASS1, TINY, HUGE, POWMAX, POWMIN
      DATA      PASS1 / .TRUE. /
c     ..


      IF( PASS1 ) THEN
         TINY   = R1MACH( 1 )
         HUGE   = R1MACH( 2 )
         POWMAX = LOG10( HUGE )
         POWMIN = LOG10( TINY )
         PASS1  = .FALSE.
      END IF


      IF( A == rZERO ) THEN
         IF( B == rZERO ) THEN
            RATIO  = rONE
         ELSE
            RATIO  = rZERO
         END IF
      ELSE IF( B == rZERO ) THEN
         RATIO  = SIGN( HUGE, A )
      ELSE
         ABSA   = ABS( A )
         ABSB   = ABS( B )
         POWA   = LOG10( ABSA )
         POWB   = LOG10( ABSB )

         IF( ABSA < TINY .AND. ABSB < TINY ) THEN
            RATIO  = rONE
         ELSE IF( POWA - POWB >= POWMAX ) THEN
            RATIO  = HUGE
         ELSE IF( POWA - POWB <= POWMIN ) THEN
            RATIO  = TINY
         ELSE
            RATIO  = ABSA / ABSB
         END IF
c                      ** DONT use old trick of determining sign
c                      ** from A*B because A*B may (over/under)flow

         IF( ( A > rZERO .AND. B < rZERO ) .OR.
     $       ( A < rZERO .AND. B > rZERO ) ) RATIO = -RATIO
      END IF

      END FUNCTION RATIO

      SUBROUTINE  ErrMsg( MESSAG, FATAL )

c        Print out a warning or error message;  abort if error
c        after making symbolic dump (machine-specific)

      LOGICAL       FATAL, MsgLim, Cray
      CHARACTER*(*) MESSAG
      INTEGER       MaxMsg, NumMsg
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

   99 FORMAT( //,' >>>>>>  TOO MANY WARNING MESSAGES --  ',
     $   'They will no longer be printed  <<<<<<<', // )

      END SUBROUTINE ErrMsg

      LOGICAL FUNCTION  WrtBad ( VarNam )

c          Write names of erroneous variables and return 'TRUE'

c      INPUT :   VarNam = Name of erroneous variable to be written
c                         ( CHARACTER, any length )

      CHARACTER*(*)  VarNam
      INTEGER        MaxMsg, NumMsg
      SAVE  NumMsg, MaxMsg
      DATA  NumMsg / 0 /,  MaxMsg / 50 /


      WrtBad = .TRUE.
      NumMsg = NumMsg + 1
      WRITE ( *, '(3A)' )  ' ****  Input variable  ', VarNam,
     $                     '  in error  ****'
      IF ( NumMsg.EQ.MaxMsg )
     $   CALL  ErrMsg ( 'Too many input errors.  Aborting...', .TRUE. )

      END FUNCTION WrtBad

      LOGICAL FUNCTION  WrtDim ( DimNam, MinVal )

c          Write name of too-small symbolic dimension and
c          the value it should be increased to;  return 'TRUE'

c      INPUT :  DimNam = Name of symbolic dimension which is too small
c                        ( CHARACTER, any length )
c               Minval = Value to which that dimension should be
c                        increased (at least)

      CHARACTER*(*)  DimNam
      INTEGER        MinVal


      WRITE ( *, '(3A,I7)' )  ' ****  Symbolic dimension  ', DimNam,
     $                     '  should be increased to at least ', MinVal
      WrtDim = .TRUE.

      END FUNCTION WrtDim

      LOGICAL FUNCTION  TstBad( VarNam, RelErr )

c       Write name (VarNam) of variable failing self-test and its
c       percent error from the correct value;  return  'FALSE'.

      CHARACTER*(*)  VarNam
      REAL           RelErr


      TstBad = .FALSE.
      WRITE( *, '(/,3A,1P,E11.2,A)' )
     $       ' Output variable ', VarNam,' differed by ', 100.*RelErr,
     $       ' per cent from correct value.  Self-test failed.'

      END FUNCTION TstBad

      FUNCTION D1MACH(i)
*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= D1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   D1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., D1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular D1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines D1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-D1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, D1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  D1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), D1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+D1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, D1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, D1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665D for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      REAL(8) :: d1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      REAL(8) :: dmach(4) 
      SAVE dmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665d(dmach)
           doinit = .FALSE.
        ENDIF
        d1mach = dmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (D1MACH) <<<  invalid argument'
        STOP
      ENDIF

*!csm
*!!! over-ride by sm on 5/26/03.  For some compilers than don't allow
* calculation of d1mach(4).  Use value found on ACD server.

c      if( i .eq. 4) d1mach = 2.22e-15

      END FUNCTION D1MACH

C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665D(DMACH)
C-----------------------------------------------------------------------
C This subroutine is a double precision version of subroutine T665R.
C See code of T665R for detailed comments and explanation
C-----------------------------------------------------------------------
      REAL(8) :: DMACH(4)
      INTEGER :: I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     $           MINEXP,MX,NEGEP,NGRD,NXRES
CS    REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
CS   1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
      REAL(8) :: A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
     $           T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
CS    CONV(I) = REAL(I)
      CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      DMACH(1) = XMIN
      DMACH(2) = XMAX
      DMACH(3) = EPSNEG
      DMACH(4) = EPS
  520 RETURN

      END SUBROUTINE T665D

      FUNCTION R1MACH(i)

*-----------------------------------------------------------------------------*
*= PURPOSE:                                                                  =*
*= R1MACH calculates various machine constants in single precision.          =*
*-----------------------------------------------------------------------------*
*= PARAMETERS:                                                               =*
*=   I       -  INTEGER, identifies the machine constant (0<I<5)         (I) =*
*=   R1MACH  -  REAL, machine constant in single precision               (O) =*
*=      I=1     - the smallest non-vanishing normalized floating-point       =*
*=                power of the radix, i.e., R1MACH=FLOAT(IBETA)**MINEXP      =*
*=      I=2     - the largest finite floating-point number.  In              =*
*=                particular R1MACH=(1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP        =*
*=                Note - on some machines R1MACH will be only the            =*
*=                second, or perhaps third, largest number, being            =*
*=                too small by 1 or 2 units in the last digit of             =*
*=                the significand.                                           =*
*=      I=3     - A small positive floating-point number such that           =*
*=                1.0-R1MACH .NE. 1.0. In particular, if IBETA = 2           =*
*=                or  IRND = 0, R1MACH = FLOAT(IBETA)**NEGEPS.               =*
*=                Otherwise,  R1MACH = (IBETA**NEGEPS)/2.  Because           =*
*=                NEGEPS is bounded below by -(IT+3), R1MACH may not         =*
*=                be the smallest number that can alter 1.0 by               =*
*=                subtraction.                                               =*
*=      I=4     - the smallest positive floating-point number such           =*
*=                that  1.0+R1MACH .NE. 1.0. In particular, if either        =*
*=                IBETA = 2  or  IRND = 0, R1MACH=FLOAT(IBETA)**MACHEP.      =*
*=                Otherwise, R1MACH=(FLOAT(IBETA)**MACHEP)/2                 =*
*=  (see routine T665R for more information on different constants)          =*
*-----------------------------------------------------------------------------*

      REAL r1mach
      INTEGER i
   
      LOGICAL doinit
      DATA doinit/.TRUE./
      SAVE doinit

      REAL rmach(4) 
      SAVE rmach

      IF (( i .GE. 1 ) .AND. ( i .LE. 4 )) THEN
* compute constants at first call only
        IF (doinit) THEN
           CALL t665r(rmach)
           doinit = .FALSE.
        ENDIF
        r1mach = rmach(i)
      ELSE
        WRITE(0,*) '>>> ERROR (R1MACH) <<<  invalid argument'
        STOP
      ENDIF

      END FUNCTION


C      ALGORITHM 665, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 14, NO. 4, PP. 303-311.
      SUBROUTINE T665R(RMACH)
C-----------------------------------------------------------------------
C  This Fortran 77 subroutine is intended to determine the parameters
C   of the floating-point arithmetic system specified below.  The
C   determination of the first three uses an extension of an algorithm
C   due to M. Malcolm, CACM 15 (1972), pp. 949-951, incorporating some,
C   but not all, of the improvements suggested by M. Gentleman and S.
C   Marovich, CACM 17 (1974), pp. 276-277.  An earlier version of this
C   program was published in the book Software Manual for the
C   Elementary Functions by W. J. Cody and W. Waite, Prentice-Hall,
C   Englewood Cliffs, NJ, 1980.
C
C  The program as given here must be modified before compiling.  If
C   a single (double) precision version is desired, change all
C   occurrences of CS (CD) in columns 1 and 2 to blanks.
C
C  Parameter values reported are as follows:
C
C       IBETA   - the radix for the floating-point representation
C       IT      - the number of base IBETA digits in the floating-point
C                 significand
C       IRND    - 0 if floating-point addition chops
C                 1 if floating-point addition rounds, but not in the
C                   IEEE style
C                 2 if floating-point addition rounds in the IEEE style
C                 3 if floating-point addition chops, and there is
C                   partial underflow
C                 4 if floating-point addition rounds, but not in the
C                   IEEE style, and there is partial underflow
C                 5 if floating-point addition rounds in the IEEE style,
C                   and there is partial underflow
C       NGRD    - the number of guard digits for multiplication with
C                 truncating arithmetic.  It is
C                 0 if floating-point arithmetic rounds, or if it
C                   truncates and only  IT  base  IBETA digits
C                   participate in the post-normalization shift of the
C                   floating-point significand in multiplication;
C                 1 if floating-point arithmetic truncates and more
C                   than  IT  base  IBETA  digits participate in the
C                   post-normalization shift of the floating-point
C                   significand in multiplication.
C       MACHEP  - the largest negative integer such that
C                 1.0+FLOAT(IBETA)**MACHEP .NE. 1.0, except that
C                 MACHEP is bounded below by  -(IT+3)
C       NEGEPS  - the largest negative integer such that
C                 1.0-FLOAT(IBETA)**NEGEPS .NE. 1.0, except that
C                 NEGEPS is bounded below by  -(IT+3)
C       IEXP    - the number of bits (decimal places if IBETA = 10)
C                 reserved for the representation of the exponent
C                 (including the bias or sign) of a floating-point
C                 number
C       MINEXP  - the largest in magnitude negative integer such that
C                 FLOAT(IBETA)**MINEXP is positive and normalized
C       MAXEXP  - the smallest positive power of  BETA  that overflows
C       EPS     - the smallest positive floating-point number such
C                 that  1.0+EPS .NE. 1.0. In particular, if either
C                 IBETA = 2  or  IRND = 0, EPS = FLOAT(IBETA)**MACHEP.
C                 Otherwise,  EPS = (FLOAT(IBETA)**MACHEP)/2
C       EPSNEG  - A small positive floating-point number such that
C                 1.0-EPSNEG .NE. 1.0. In particular, if IBETA = 2
C                 or  IRND = 0, EPSNEG = FLOAT(IBETA)**NEGEPS.
C                 Otherwise,  EPSNEG = (IBETA**NEGEPS)/2.  Because
C                 NEGEPS is bounded below by -(IT+3), EPSNEG may not
C                 be the smallest number that can alter 1.0 by
C                 subtraction.
C       XMIN    - the smallest non-vanishing normalized floating-point
C                 power of the radix, i.e.,  XMIN = FLOAT(IBETA)**MINEXP
C       XMAX    - the largest finite floating-point number.  In
C                 particular  XMAX = (1.0-EPSNEG)*FLOAT(IBETA)**MAXEXP
C                 Note - on some machines  XMAX  will be only the
C                 second, or perhaps third, largest number, being
C                 too small by 1 or 2 units in the last digit of
C                 the significand.
C
C     Latest revision - April 20, 1987
C
C     Author - W. J. Cody
C              Argonne National Laboratory
C
C-----------------------------------------------------------------------
      REAL rmach(4)
      INTEGER I,IBETA,IEXP,IRND,IT,ITEMP,IZ,J,K,MACHEP,MAXEXP,
     1        MINEXP,MX,NEGEP,NGRD,NXRES
      REAL A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,T,TEMP,TEMPA,
     1     TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
CD    REAL(8) A,B,BETA,BETAIN,BETAH,CONV,EPS,EPSNEG,ONE,
CD   1                 T,TEMP,TEMPA,TEMP1,TWO,XMAX,XMIN,Y,Z,ZERO
C-----------------------------------------------------------------------
      CONV(I) = REAL(I)
CD    CONV(I) = DBLE(I)
      ONE = CONV(1)
      TWO = ONE + ONE
      ZERO = ONE - ONE
C-----------------------------------------------------------------------
C  Determine IBETA, BETA ala Malcolm.
C-----------------------------------------------------------------------
      A = ONE
   10 A = A + A
         TEMP = A+ONE
         TEMP1 = TEMP-A
         IF (TEMP1-ONE .EQ. ZERO) GO TO 10
      B = ONE
   20 B = B + B
         TEMP = A+B
         ITEMP = INT(TEMP-A)
         IF (ITEMP .EQ. 0) GO TO 20
      IBETA = ITEMP
      BETA = CONV(IBETA)
C-----------------------------------------------------------------------
C  Determine IT, IRND.
C-----------------------------------------------------------------------
      IT = 0
      B = ONE
  100 IT = IT + 1
         B = B * BETA
         TEMP = B+ONE
         TEMP1 = TEMP-B
         IF (TEMP1-ONE .EQ. ZERO) GO TO 100
      IRND = 0
      BETAH = BETA / TWO
      TEMP = A+BETAH
      IF (TEMP-A .NE. ZERO) IRND = 1
      TEMPA = A + BETA
      TEMP = TEMPA+BETAH
      IF ((IRND .EQ. 0) .AND. (TEMP-TEMPA .NE. ZERO)) IRND = 2
C-----------------------------------------------------------------------
C  Determine NEGEP, EPSNEG.
C-----------------------------------------------------------------------
      NEGEP = IT + 3
      BETAIN = ONE / BETA
      A = ONE
      DO 200 I = 1, NEGEP
         A = A * BETAIN
  200 CONTINUE
      B = A
  210 TEMP = ONE-A
         IF (TEMP-ONE .NE. ZERO) GO TO 220
         A = A * BETA
         NEGEP = NEGEP - 1
      GO TO 210
  220 NEGEP = -NEGEP
      EPSNEG = A
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 300
      A = (A*(ONE+A)) / TWO
      TEMP = ONE-A
      IF (TEMP-ONE .NE. ZERO) EPSNEG = A
C-----------------------------------------------------------------------
C  Determine MACHEP, EPS.
C-----------------------------------------------------------------------
  300 MACHEP = -IT - 3
      A = B
  310 TEMP = ONE+A
         IF (TEMP-ONE .NE. ZERO) GO TO 320
         A = A * BETA
         MACHEP = MACHEP + 1
      GO TO 310
  320 EPS = A
      TEMP = TEMPA+BETA*(ONE+EPS)
      IF ((IBETA .EQ. 2) .OR. (IRND .EQ. 0)) GO TO 350
      A = (A*(ONE+A)) / TWO
      TEMP = ONE+A
      IF (TEMP-ONE .NE. ZERO) EPS = A
C-----------------------------------------------------------------------
C  Determine NGRD.
C-----------------------------------------------------------------------
  350 NGRD = 0
      TEMP = ONE+EPS
      IF ((IRND .EQ. 0) .AND. (TEMP*ONE-ONE .NE. ZERO)) NGRD = 1
C-----------------------------------------------------------------------
C  Determine IEXP, MINEXP, XMIN.
C
C  Loop to determine largest I and K = 2**I such that
C         (1/BETA) ** (2**(I))
C  does not underflow.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
      I = 0
      K = 1
      Z = BETAIN
      T = ONE + EPS
      NXRES = 0
  400 Y = Z
         Z = Y * Y
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Z * ONE
         TEMP = Z * T
         IF ((A+A .EQ. ZERO) .OR. (ABS(Z) .GE. Y)) GO TO 410
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .EQ. Z) GO TO 410
         I = I + 1
         K = K + K
      GO TO 400
  410 IF (IBETA .EQ. 10) GO TO 420
      IEXP = I + 1
      MX = K + K
      GO TO 450
C-----------------------------------------------------------------------
C  This segment is for decimal machines only.
C-----------------------------------------------------------------------
  420 IEXP = 2
      IZ = IBETA
  430 IF (K .LT. IZ) GO TO 440
         IZ = IZ * IBETA
         IEXP = IEXP + 1
      GO TO 430
  440 MX = IZ + IZ - 1
C-----------------------------------------------------------------------
C  Loop to determine MINEXP, XMIN.
C  Exit from loop is signaled by an underflow.
C-----------------------------------------------------------------------
  450 XMIN = Y
         Y = Y * BETAIN
C-----------------------------------------------------------------------
C  Check for underflow here.
C-----------------------------------------------------------------------
         A = Y * ONE
         TEMP = Y * T
         IF (((A+A) .EQ. ZERO) .OR. (ABS(Y) .GE. XMIN)) GO TO 460
         K = K + 1
         TEMP1 = TEMP * BETAIN
         IF (TEMP1*BETA .NE. Y) GO TO 450
      NXRES = 3
      XMIN = Y
  460 MINEXP = -K
C-----------------------------------------------------------------------
C  Determine MAXEXP, XMAX.
C-----------------------------------------------------------------------
      IF ((MX .GT. K+K-3) .OR. (IBETA .EQ. 10)) GO TO 500
      MX = MX + MX
      IEXP = IEXP + 1
  500 MAXEXP = MX + MINEXP
C-----------------------------------------------------------------
C  Adjust IRND to reflect partial underflow.
C-----------------------------------------------------------------
      IRND = IRND + NXRES
C-----------------------------------------------------------------
C  Adjust for IEEE-style machines.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 2) .OR. (IRND .EQ. 5)) MAXEXP = MAXEXP - 2
C-----------------------------------------------------------------
C  Adjust for non-IEEE machines with partial underflow.
C-----------------------------------------------------------------
      IF ((IRND .EQ. 3) .OR. (IRND .EQ. 4)) MAXEXP = MAXEXP - IT
C-----------------------------------------------------------------
C  Adjust for machines with implicit leading bit in binary
C  significand, and machines with radix point at extreme
C  right of significand.
C-----------------------------------------------------------------
      I = MAXEXP + MINEXP
      IF ((IBETA .EQ. 2) .AND. (I .EQ. 0)) MAXEXP = MAXEXP - 1
      IF (I .GT. 20) MAXEXP = MAXEXP - 1
      IF (A .NE. Y) MAXEXP = MAXEXP - 2
      XMAX = ONE - EPSNEG
      IF (XMAX*ONE .NE. XMAX) XMAX = ONE - BETA * EPSNEG
      XMAX = XMAX / (BETA * BETA * BETA * XMIN)
      I = MAXEXP + MINEXP + 3
      IF (I .LE. 0) GO TO 520
      DO 510 J = 1, I
          IF (IBETA .EQ. 2) XMAX = XMAX + XMAX
          IF (IBETA .NE. 2) XMAX = XMAX * BETA
  510 CONTINUE
      RMACH(1) = XMIN
      RMACH(2) = XMAX
      RMACH(3) = EPSNEG
      RMACH(4) = EPS
  520 RETURN

      END SUBROUTINE T665R

      end module RTRANS
