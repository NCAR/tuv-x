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

c     spherical geometry
      REAL tausla(0:NLYR), tauslau(0:NLYR), mu2(0:NLYR)

c     .. \todo Figure out where these came from in original code
      INTEGER, PARAMETER :: MXPHI = 3
      LOGICAL, PARAMETER :: PRNT(7) = .FALSE.
      LOGICAL, PARAMETER :: ONLYFL = .TRUE.
      LOGICAL, PARAMETER :: USRANG = .FALSE.

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
      REAL    :: AZERR, AZTERM, BPLANK, COSPHI, DELM0, DITHER,
     &           DUM, RPD, SGN, TPLANK
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

      SAVE      PASS1, DITHER, RPD
      DATA      PASS1 / .TRUE. /

      IF( PASS1 ) THEN

         DITHER = 10.*R1MACH( 4 )

c                            ** Must dither more on Cray (14-digit prec)

         IF( DITHER.LT.1.E-10 ) DITHER = 10.*DITHER

         RPD  = PI / 180.0
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
     &           ZJ, ZZ(:,LC), OPRIM(LC), LC, DITHER, mu2(lc),
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
