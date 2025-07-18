      SUBROUTINE schum(nz, o2col,tlev,secchi,dto2,o2xsk,pathname)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the equivalent absorption cross section of O2 in the SR bands. =*
*=  The algorithm is based on parameterization of G.A. Koppers, and          =*
*=  D.P. Murtagh [ref. Ann.Geophys., 14 68-79, 1996]                         =*
*=  Final values do include effects from the Herzberg continuum.             =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  TLEV    - tmeperature at each level                                   (I)=*
*=  SECCHI  - ratio of slant to vertical o2 columns                       (I)=*
*=  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer at each specified wavelength                    =*
*=  O2XSK  - REAL, molecular absorption cross section in SR bands at      (O)=*
*=            each specified wavelength.  Includes Herzberg continuum        =*
*-----------------------------------------------------------------------------*


      IMPLICIT NONE
      INCLUDE 'params'

      INTEGER nz
      REAL o2col(kz), o2col1(kz)
      REAL tlev(kz), secchi(kz)

      REAL dto2(kz,17), o2xsk(kz,17)

      CHARACTER*80 pathname
      INTEGER i, k, ktop, ktop1, kbot

      REAL XS(17), X
      REAL xslod(17)
      LOGICAL firstcall
      SAVE firstcall
      DATA firstcall /.TRUE./

      DATA xslod  /6.2180730E-21, 5.8473627E-22, 5.6996334E-22,
     $             4.5627094E-22, 1.7668250E-22, 1.1178808E-22,
     $             1.2040544E-22, 4.0994668E-23, 1.8450616E-23,
     $             1.5639540E-23, 8.7961075E-24, 7.6475608E-24,
     $             7.6260556E-24, 7.5565696E-24, 7.6334338E-24,
     $             7.4371992E-24, 7.3642966E-24/
c------------------------------------------
C	 Initialize values
c------------------------------------------   
      dto2(:,:) = 0.0

c------------------------------------------
c sm	 Initialize cross sections to values
c sm	 at large optical depth
c------------------------------------------

      DO k = 1, nz
         DO i = 1, 17
            o2xsk(k,i) = xslod(i)
         ENDDO	
      ENDDO

c------------------------------------------
c      Loads Chebyshev polynomial Coeff.
c------------------------------------------

      if (firstcall) then 
        call INIT_XS(pathname)
	firstcall = .FALSE.
      endif

c------------------------------------------
c     Calculate cross sections
c sm:  Set smallest O2col = exp(38.) molec cm-2
c sm     to stay in range of parameterization
c sm     given by Koppers et al. at top of atm.
c------------------------------------------

      ktop = nz
      kbot = 0

c     EXP(38.) = 3.185e16
c     EXP(56.) = 2.091e24
      DO k=1,nz    !! loop for alt
         o2col1(k) = MAX(o2col(k),EXP(38.))

         x  = ALOG(o2col1(k))
         
         IF (x .LT. 38.0) THEN
            ktop1 = k-1
            write(*,*) ktop1
            ktop  = MIN(ktop1,ktop)
         ELSE IF (x .GT. 56.0) THEN
            kbot = k
         ELSE
            CALL effxs( x, tlev(k), xs )
            DO i=1,17
               o2xsk(k,i) = xs(i)
            END DO
         ENDIF

      END DO                    !! finish loop for alt

c------------------------------------------
c  fill in cross section where X is out of range 
c  by repeating edge table values
c------------------------------------------
       
c sm do not allow kbot = nz to avoid division by zero in
c   no light case.
       
      IF(kbot .EQ. nz) kbot = nz - 1

      DO k=1,kbot
         DO i=1,17
            o2xsk(k,i) = o2xsk(kbot+1,i)
         END DO
      END DO
      
      DO k=ktop+1,nz
         DO i=1,17
            o2xsk(k,i) = o2xsk(ktop,i)
         END DO
      END DO

c------------------------------------------
c  Calculate incremental optical depths 
c------------------------------------------

      DO i=1,17                   ! loop over wavelength

         DO k=1,nz-1            ! loop for alt

c... calculate an optical depth weighted by density
c sm:  put in mean value estimate, if in shade

            IF (ABS(1. - o2col1(k+1)/o2col1(k)) .LE. 2.*precis) THEN

               dto2(k,i) = o2xsk(k+1,i)*o2col1(k+1)/(nz-1)

            ELSE

            dto2(k,i) = ABS(
     $           ( o2xsk(k+1,i)*o2col1(k+1) - o2xsk(k,i)*o2col1(k) )
     $           / ( 1. + ALOG(o2xsk(k+1,i)/o2xsk(k,i)) 
     $           / ALOG(o2col1(k+1)/o2col1(k)) ) )

c... change to vertical optical depth

            dto2(k,i) = 2. * dto2(k,i) / (secchi(k)+secchi(k+1))

            ENDIF

         END DO
         dto2(nz,i) = 0.0       ! set optical depth to zero at top


      END DO 

      return
      end

C-------------------------------------------------------------
      SUBROUTINE EFFXS( X, T, XS )
C-------------------------------------------------------------
C
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
C
C-------------------------------------------------------------

	IMPLICIT NONE

	REAL*4 NO2, T, X
	REAL*4 XS(17)
	REAL*4 A(17), B(17) 
	INTEGER I

	CALL CALC_PARAMS( X, A, B )

	DO I = 1,17
	  XS(I) = EXP( A(I)*( T - 220.) + B(I) )
	ENDDO

        RETURN

	END

C-------------------------------------------------------------
	SUBROUTINE CALC_PARAMS( X, A, B )
C-------------------------------------------------------------
C
C       calculates coefficients (A,B), used in calculating the
C	effective cross section, for 17 wavelength intervals
C       as a function of log O2 column density (X)
C       Wavelength intervals are defined in WMO1985
C
C-------------------------------------------------------------

	IMPLICIT NONE

	REAL*4 X
	REAL*4 A(17), B(17)

	REAL*4   CHEBEV

	REAL*8 AC(20,17)
        REAL*8 BC(20,17) ! Chebyshev polynomial coeffs
	REAL*4 WAVE_NUM(17)
	COMMON /XS_COEFFS/ AC, BC, WAVE_NUM

	INTEGER I

C       call Chebyshev Evaluation routine to calc A and B from
C	set of 20 coeficients for each wavelength

	DO I=1,17
	  A(I) = CHEBEV(38.0 , 56.0, AC(1,I), 20, X)
	  B(I) = CHEBEV(38.0 , 56.0, BC(1,I), 20, X)
	ENDDO

	RETURN

	END

C-------------------------------------------------------------
	SUBROUTINE INIT_XS(pathname)
C-------------------------------------------------------------
C
C       loads COMMON block XS_COEFFS containing the Chebyshev
C	polynomial coeffs necessary to calculate O2 effective
C       cross-sections
C
C-------------------------------------------------------------
	REAL*8 AC(20,17)
	REAL*8 BC(20,17) ! Chebyshev polynomial coeffs
	REAL*4 WAVE_NUM(17)
	COMMON /XS_COEFFS/ AC, BC, WAVE_NUM
	
        Character*80 pathname

C       locals
	INTEGER*4 IN_LUN	! file unit number
	INTEGER*4 IOST		! i/o status
	INTEGER*4 I, J

        IN_LUN = 11

	OPEN (UNIT=IN_LUN, FILE=
     $       TRIM(pathname)//'effxstex.txt',FORM='FORMATTED')

	READ( IN_LUN, 901 )
	DO I = 1,20
	  READ( IN_LUN, 903 ) ( AC(I,J), J=1,17 )
	ENDDO
	READ( IN_LUN, 901 )
	DO I = 1,20
	  READ( IN_LUN, 903 ) ( BC(I,J), J=1,17 )
	ENDDO

 901    FORMAT( / )
 903    FORMAT( 17(E23.14,1x))

 998	CLOSE (IN_LUN)
	
	DO I=1,17
	  WAVE_NUM(18-I) = 48250. + (500.*I)
	ENDDO

        END

C-------------------------------------------------------------
	FUNCTION chebev(a,b,c,m,x)
C-------------------------------------------------------------
C
C     Chebyshev evaluation algorithm
C     See Numerical recipes p193
C
C-------------------------------------------------------------
      
	INTEGER M
        REAL*4 CHEBEV,A,B,X
	REAL*8 C(M)
        INTEGER J
        REAL D,DD,SV,Y,Y2

        IF ((X-A)*(X-B).GT.0.) THEN
	  WRITE(6,*) 'X NOT IN RANGE IN CHEBEV', X
	  CHEBEV = 0.0
	  RETURN
        ENDIF

	D=0.
        DD=0.
        Y=(2.*X-A-B)/(B-A)
        Y2=2.*Y
        DO 11 J=M,2,-1
          SV=D
          D=Y2*D-DD+C(J)
          DD=SV
 11     CONTINUE
        CHEBEV=Y*D-DD+0.5*C(1)
      
	RETURN
        END

