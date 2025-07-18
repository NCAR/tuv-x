      MODULE LA_SRB_MOD

      IMPLICIT NONE

      private
      public :: init_la_srb, la_srb, sjo2

      INTEGER :: ila
      INTEGER, PARAMETER :: kla = 2
      INTEGER, PARAMETER :: nla = kla - 1
      REAL, PARAMETER :: wlla(kla) = (/ 121.4, 121.9 /)

      INTEGER :: isrb
      INTEGER, PARAMETER :: ksrb = 18
      INTEGER, PARAMETER :: nsrb = ksrb - 1
      REAL, PARAMETER :: wlsrb(ksrb) =
     $   (/175.4, 177.0, 178.6, 180.2, 181.8, 183.5, 185.2, 186.9,
     $     188.7, 190.5, 192.3, 194.2, 196.1, 198.0, 200.0, 202.0, 
     $     204.1, 206.2/)

      REAL(8) :: AC(20,nsrb)
      REAL(8) :: BC(20,nsrb) ! Chebyshev polynomial coeffs

      contains

      SUBROUTINE init_la_srb(wl)
*----------------------------------------------------------------------
* check that the user wavelength grid, WL(IW), is compatible 
* with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
* Also compute and save corresponding grid indices (ILA, ISRB)
*----------------------------------------------------------------------

      use tuv_params, only : precis

      real, intent(in) :: wl(:)

      integer :: iw, nw

      IF(wl(1) <= wlsrb(nsrb)) THEN
  
      nw = size(wl)
** locate Lyman-alpha wavelengths on grid
      ila = 0
      DO iw = 1, nw
        IF(ABS(wl(iw) - wlla(1)) < 10.*precis) THEN
          ila = iw
          EXIT
        ENDIF
      ENDDO

* check 
      IF(ila == 0) THEN
         WRITE(*,*) 'For wavelengths below 205.8 nm, only the'
         WRITE(*,*) 'pre-specified wavelength grid is permitted'
         WRITE(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
         STOP ' Lyman alpha grid mis-match - 1'
      ENDIF
      DO iw = 2, nla + 1
         IF(ABS(wl(ila + iw - 1) - wlla(iw)) > 10.*precis) THEN
            WRITE(*,*) 'Lyman alpha grid mis-match - 2'
            STOP
         ENDIF
      ENDDO

** locate Schumann-Runge wavelengths on grid
      isrb = 0
      DO iw = 1, nw
         IF(ABS(wl(iw) - wlsrb(1)) < 10.*precis) THEN
            isrb = iw
            EXIT
         ENDIF
      ENDDO

* check
      IF(isrb == 0) THEN
         WRITE(*,*) 'For wavelengths below 205.8 nm, only the'
         WRITE(*,*) 'pre-specified wavelength grid is permitted'
         WRITE(*,*) 'Use nwint=-156, or edit subroutine gridw.f'
         STOP ' SRB grid mis-match - 1'
      ENDIF
      DO iw = 2, nsrb + 1
         IF(ABS(wl(isrb + iw - 1) - wlsrb(iw)) > 10.* precis) THEN
            WRITE(*,*) ' SRB grid mismatch - w'
            STOP
         ENDIF
      ENDDO

c------------------------------------------
c      Loads Chebyshev polynomial Coeff.
c------------------------------------------
      CALL INIT_XS

      ENDIF

      END SUBROUTINE init_la_srb

      SUBROUTINE la_srb(z,tlev,wl,vcol,scol,o2xs1,dto2,o2xs)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute equivalent optical depths for O2 absorption, and O2 effective    =*
*=  absorption cross sections, parameterized in the Lyman-alpha and SR bands =*
*-----------------------------------------------------------------------------* 
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)=*
*=            working wavelength grid                                        =*
*=  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)=*
*=            altitude layer                                                 =*
*=  ZEN     - REAL, solar zenith angle                                    (I)=*
*=                                                                           =*
*=  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)=*
*=                                                                           =*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified altitude and wavelength.  Includes Herzberg     =*
*=            continuum.                                                     =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : largest, precis


      REAL, intent(in) :: wl(:)
      REAL, intent(in) :: z(:)

      REAL, intent(in)  :: tlev(:)
      REAL, intent(in)  :: vcol(:), scol(:)
      REAL, intent(in)  :: o2xs1(:)
      REAL, intent(inout) :: dto2(:,:), o2xs(:,:)

      INTEGER :: nz, nw, i, iz, iw
      REAL    :: secchi(size(z))
      REAL    :: o2col(size(scol))

* Lyman-alpha variables
* O2 optical depth and equivalent cross section in the Lyman-alpha region
* Wavelengths for Lyman alpha and SRB parameterizations:

      REAL    :: dto2la(size(vcol),kla-1), o2xsla(size(z),kla-1)

* grid on which Koppers' parameterization is defined
* O2 optical depth and equivalent cross section on Koppers' grid

      REAL    :: dto2k(size(vcol),ksrb-1), o2xsk(size(z),ksrb-1)

      nw = size(wl)
*----------------------------------------------------------------------
* initalize O2 cross sections 
*----------------------------------------------------------------------
      DO iw = 1, nw - 1   
         o2xs(:,iw) = o2xs1(iw)
      ENDDO

      IF(wl(1) <= wlsrb(nsrb)) THEN
        nz = size(z)
*----------------------------------------------------------------------
* Slant O2 column and x-sections.
*----------------------------------------------------------------------
        o2col(:) = 0.2095 * scol(:)
*----------------------------------------------------------------------
* Effective secant of solar zenith angle.  
* Use 2.0 if no direct sun (value for isotropic radiation)
* For nz, use value at nz-1
*----------------------------------------------------------------------
        WHERE( scol(1:nz-1) <= .1*largest )
          secchi(1:nz-1) = scol(1:nz-1)/vcol(1:nz-1)
        ELSEWHERE
          secchi(1:nz-1) = 2.
        ENDWHERE
        secchi(nz) = secchi(nz-1)
*---------------------------------------------------------------------
* Lyman-Alpha parameterization, output values of O2 optical depth
* and O2 effective (equivalent) cross section
*----------------------------------------------------------------------
        CALL lymana(o2col,secchi,dto2la,o2xsla)
        DO iw = ila, ila + nla - 1
          dto2(:,iw) = dto2la(:,iw - ila + 1)
          o2xs(:,iw) = o2xsla(:,iw - ila + 1)
        ENDDO
*------------------------------------------------------------------------------
* Koppers' parameterization of the SR bands, output values of O2
* optical depth and O2 equivalent cross section 
*------------------------------------------------------------------------------
        CALL schum(o2col,tlev,secchi,dto2k,o2xsk)
        DO iw = isrb, isrb + nsrb - 1
          dto2(:,iw) = dto2k(:,iw - isrb + 1)
          o2xs(:,iw) = o2xsk(:,iw - isrb + 1)
        ENDDO
      ENDIF

      END SUBROUTINE la_srb

*=============================================================================*

      SUBROUTINE lymana(o2col,secchi,dto2la,o2xsla)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
*=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
*=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
*=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
*=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer                                                 =*
*=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
*-----------------------------------------------------------------------------*

* input:

      REAL, intent(in) :: o2col(:)
      REAL, intent(in) :: secchi(:)

* output

      REAL, intent(out) :: dto2la(:,:), o2xsla(:,:)

* local variables

      REAL, parameter :: xsmin = 1.e-20
      real(8), parameter :: b(3) = 
     $      (/ 6.8431D-01, 2.29841D-01,  8.65412D-02 /)
      real(8), parameter :: c(3) = 
     $      (/8.22114D-21, 1.77556D-20,  8.22112D-21 /)
      real(8), parameter :: d(3) = 
     $      (/ 6.0073D-21, 4.28569D-21,  1.28059D-20 /)
      real(8), parameter :: e(3) = 
     $      (/ 8.21666D-21, 1.63296D-20,  4.85121D-17 /)

      INTEGER :: nz
      INTEGER :: iz, i
      REAL(8) :: coldens
      REAL(8) :: sigma(3), tau(3)
      REAL(8) :: rm(size(o2col)), ro2(size(o2col))

* calculate reduction factors at every altitude

      rm(:)  = 0.0_8
      ro2(:) = 0.0_8
      nz = size(o2col)
      DO iz = 1, nz
        coldens = REAL(o2col(iz),kind=8)
        sigma   = c * coldens
        WHERE( sigma < 100.e8_8 )
          tau = EXP( -sigma )
        ELSEWHERE
          tau = 0.0_8
        ENDWHERE
        rm(iz) = dot_product( tau,b )

        sigma   = e * coldens
        WHERE( sigma < 100.e8_8 )
          tau = EXP( -sigma ) 
        ELSEWHERE
          tau = 0.0_8
        ENDWHERE
        ro2(iz) = dot_product( tau,d )
      ENDDO

* calculate effective O2 optical depths and effective O2 cross sections

      DO iz = 1, nz-1
         IF (rm(iz) > 1.0D-100) THEN
            IF (ro2(iz) > 1.D-100) THEN
               o2xsla(iz,1) = ro2(iz)/rm(iz)
            ELSE
               o2xsla(iz,1) = xsmin
            ENDIF

            IF (rm(iz+1) > 0._8) THEN
               dto2la(iz,1) = LOG(rm(iz+1)) / secchi(iz+1) 
     $                      - LOG(rm(iz))   / secchi(iz)
            ELSE
               dto2la(iz,1) = 1000.
            ENDIF
         ELSE
            dto2la(iz,1) = 1000.
            o2xsla(iz,1) = xsmin
         ENDIF
      ENDDO

* do top layer separately

      IF(rm(nz) > 1.D-100) THEN
         o2xsla(nz,1) = ro2(nz)/rm(nz)
      ELSE
         o2xsla(nz,1) = xsmin
      ENDIF

      END SUBROUTINE LYMANA

*=============================================================================*

      SUBROUTINE SCHUM(o2col, tlev, secchi, dto2, o2xsk)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
*=  The algorithm is based on parameterization of G.A. Koppers, and          =*
*=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
*=  Final values do include effects from the Herzberg continuum.             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  TLEV    - tmeperature at each level                                   (I)=*
*=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  O2XSK  - REAL, molecular absorption cross section in SR bands at     (O)=*
*=            each specified wavelength.  Includes Herzberg continuum        =*
*-----------------------------------------------------------------------------*
    
      use tuv_params, only : precis

      REAL, intent(in) :: o2col(:)
      REAL, intent(in) :: tlev(:), secchi(:)

      REAL, intent(out) :: dto2(:,:), o2xsk(:,:)

      INTEGER :: nz
      INTEGER :: i, k, ktop, ktop1, kbot

      REAL, parameter :: colmin = exp( 38. )
      REAL :: X, NORM
      REAL :: XS(nsrb)
      REAL :: o2col1(size(o2col))

      real, parameter :: xslod(nsrb) =
     $          (/ 6.2180730E-21, 5.8473627E-22, 5.6996334E-22,
     $             4.5627094E-22, 1.7668250E-22, 1.1178808E-22,
     $             1.2040544E-22, 4.0994668E-23, 1.8450616E-23,
     $             1.5639540E-23, 8.7961075E-24, 7.6475608E-24,
     $             7.6260556E-24, 7.5565696E-24, 7.6334338E-24,
     $             7.4371992E-24, 7.3642966E-24 /)

c------------------------------------------
*sm	 Initialize cross sections to values
*sm	 at large optical depth
c------------------------------------------
      DO i = 1, nsrb
         o2xsk(:,i) = xslod(i)
      ENDDO

c------------------------------------------
c     Calculate cross sections
*sm:  Set smallest O2col = exp(38.) molec cm-2
*sm     to stay in range of parameterization
*sm     given by Koppers et al. at top of atm.
c------------------------------------------
      nz = size(o2col)
      ktop = nz
      kbot = 0
      NORM = 1./REAL(nz-1)
      DO k = 1,nz
         o2col1(k) = MAX(o2col(k),colmin)
         x  = LOG(o2col1(k))
         
         IF (x < 38.0) THEN
            ktop1 = k-1
            ktop  = MIN(ktop1,ktop)
         ELSE IF (x > 56.0) THEN
            kbot = k
         ELSE
            CALL EFFXS( x, tlev(k), xs )
            o2xsk(k,1:nsrb) = xs(1:nsrb)
         ENDIF
      END DO

c------------------------------------------
c  fill in cross section where X is out of range 
c  by repeating edge table values
c------------------------------------------
*sm do not allow kbot = nz to avoid division by zero in
*   no light case.
       
      IF(kbot == nz) kbot = nz - 1

      DO i=1,nsrb
         o2xsk(1:kbot,i) = o2xsk(kbot+1,i)
      END DO
      
      DO i=1,nsrb
         o2xsk(ktop+1:nz,i) = o2xsk(ktop,i)
      END DO

c------------------------------------------
c  Calculate incremental optical depth
c------------------------------------------
      DO i = 1,nsrb
         DO k = 1,nz-1
c... calculate an optical depth weighted by density
*sm:  put in mean value estimate, if in shade

            IF (ABS(1. - o2col1(k+1)/o2col1(k)) <= 2.*precis) THEN
               dto2(k,i) = o2xsk(k+1,i)*o2col1(k+1)*NORM
            ELSE
               dto2(k,i) = ABS(
     $           ( o2xsk(k+1,i)*o2col1(k+1) - o2xsk(k,i)*o2col1(k) )
     $           / ( 1. + LOG(o2xsk(k+1,i)/o2xsk(k,i)) 
     $           / LOG(o2col1(k+1)/o2col1(k)) ) )
c... change to vertical optical depth
               dto2(k,i) = 2. * dto2(k,i) / (secchi(k)+secchi(k+1))
            ENDIF
         END DO
      END DO 

      END SUBROUTINE SCHUM

      SUBROUTINE INIT_XS
C-------------------------------------------------------------
C	polynomial coeffs necessary to calculate O2 effective
C       cross-sections
C-------------------------------------------------------------

      INTEGER :: IN_LUN	! file unit number
      INTEGER :: IOST		! i/o status
      INTEGER :: I, J

      IN_LUN = 11

      OPEN (UNIT=IN_LUN, FILE=
     $  'odat/DATAE1/O2/effxstex.txt',FORM='FORMATTED')

      READ( IN_LUN, 901 )
      DO I = 1,20
        READ( IN_LUN, 903 ) ( AC(I,J), J=1,nsrb )
      ENDDO
      READ( IN_LUN, 901 )
      DO I = 1,20
        READ( IN_LUN, 903 ) ( BC(I,J), J=1,nsrb )
      ENDDO

 901  FORMAT( / )
 903  FORMAT( 17(E23.14,1x))

      CLOSE (IN_LUN)

      END SUBROUTINE INIT_XS

*=============================================================================*

      FUNCTION CHEBEV(a,b,c,m,x) result(poly)
C-------------------------------------------------------------
C
C     Chebyshev evaluation algorithm
C     See Numerical recipes p193
C
C-------------------------------------------------------------
      
      INTEGER, intent(in) :: M
      REAL, intent(in)    :: A, B, X
      REAL(8), intent(in) :: C(M)

      REAL :: poly

      INTEGER :: J
      REAL    :: D,DD,SV,Y,Y2


      IF ((X-A)*(X-B) > 0.) THEN
        WRITE(6,*) 'X NOT IN RANGE IN CHEBEV', X
	poly = 0.0
      ELSE
         D = 0.
         DD = 0.
         Y=(2.*X-A-B)/(B-A)
         Y2=2.*Y
         DO J=M,2,-1
           SV=D
           D=Y2*D-DD+C(J)
           DD=SV
         ENDDO
         poly = Y*D - DD + 0.5*C(1)
      ENDIF
      
      END FUNCTION CHEBEV

      SUBROUTINE EFFXS( X, T, XS )
C-------------------------------------------------------------
C     Subroutine for evaluating the effective cross section
C     of O2 in the Schumann-Runge bands using parameterization
C     of G.A. Koppers, and D.P. Murtagh [ref. Ann.Geophys., 14
C     68-79, 1996]
C      
C     method:
C     ln(xs) = A(X)[T-220]+B(X)
C     X = log of slant column of O2
C     A,B calculated from Chebyshev polynomial coeffs
C     AC and BC using NR routine chebev.  Assume interval
C     is 38<ln(NO2)<56.
C
C     Revision History:
C
C     drm 2/97  initial coding
C-------------------------------------------------------------

	REAL, intent(in)  :: T, X
	REAL, intent(out) :: XS(:)

	INTEGER :: I
	REAL    :: A(nsrb), B(nsrb) 

	CALL CALC_PARAMS( X, A, B )

	XS(:) = EXP( A(:)*( T - 220.) + B(:) )

	END SUBROUTINE EFFXS

*=============================================================================*

	SUBROUTINE CALC_PARAMS( X, A, B )
C-------------------------------------------------------------
C
C       calculates coefficients (A,B), used in calculating the
C	effective cross section, for 17 wavelength intervals
C       as a function of log O2 column density (X)
C       Wavelength intervals are defined in WMO1985
C
C-------------------------------------------------------------

	REAL, intent(in)  :: X
	REAL, intent(out) :: A(:), B(:)

	INTEGER :: I

C       call Chebyshev Evaluation routine to calc A and B from
C	set of 20 coeficients for each wavelength

	DO I=1,size(A)
	  A(I) = CHEBEV(38.0 , 56.0, AC(1,I), 20, X)
	  B(I) = CHEBEV(38.0 , 56.0, BC(1,I), 20, X)
	ENDDO

	END SUBROUTINE CALC_PARAMS

*=============================================================================*

       SUBROUTINE sjo2(xso2,sq)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Update the weighting function (cross section x quantum yield) for O2     =*
*=  photolysis.  The strong spectral variations in the O2 cross sections are =*
*=  parameterized into a few bands for Lyman-alpha (121.4-121.9 nm, one band)=*
*=  and Schumann-Runge (174.4-205.8, 17 bands) regions. The parameterizations=*
*=  depend on the overhead O2 column, and therefore on altitude and solar    =*
*=  zenith angle, so they need to be updated at each time/zenith step.       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  XSO2   - REAL, molecular absorption cross section in SR bands at      (I)=*
*=           each specified altitude and wavelength.  Includes Herzberg      =*
*=            continuum.                                                     =*
*=  NJ     - INTEGER, index of O2 photolysis in array SQ                  (I)=*
*=  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)=*
*=           photolysis reaction, at each wavelength and each altitude level =*
*-----------------------------------------------------------------------------*

* calling parameters

      REAL, intent(in)    :: xso2(:,:)
      REAL, intent(inout) :: sq(:,:)

* local

      INTEGER :: iw, nz

* O2 + hv -> O + O
* quantum yield assumed to be unity
* assign cross section values at all wavelengths and at all altitudes
*      qy = 1.

      nz = size(xso2,dim=1)
      DO iw = 1, size(xso2,dim=2)
         sq(1:nz,iw) = xso2(1:nz,iw)
      ENDDO

      END SUBROUTINE sjo2

      END MODULE LA_SRB_MOD
