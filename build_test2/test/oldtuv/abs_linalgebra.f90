
   module abs_linalgebra

   implicit none

   private
   public :: abs_linalgebra_t, sscal, saxpy, sasum, sdot, isamax

   type, abstract :: abs_linalgebra_t
     contains
     procedure(SGBCO), deferred :: SGBCO
     procedure(SGBFA), deferred :: SGBFA
     procedure(SGBSL), deferred :: SGBSL
     procedure(SGECO), deferred :: SGECO
     procedure(SGEFA), deferred :: SGEFA
     procedure(SGESL), deferred :: SGESL
     procedure(tridiag), deferred :: tridiag
   end type abs_linalgebra_t

   interface
      SUBROUTINE SGBCO( this, ABD, N, ML, MU, IPVT, RCOND, Z )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in)  :: N, ML, MU
        INTEGER, intent(out) :: IPVT(:)
	REAL, intent(out)    :: RCOND
	REAL, intent(inout)  :: ABD(:,:)
	REAL, intent(out)    :: Z(:)
      END SUBROUTINE SGBCO

      SUBROUTINE SGBFA( this, ABD, N, ML, MU, IPVT, INFO )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in)  :: N, ML, MU
        INTEGER, intent(out) ::  INFO
        INTEGER, intent(out) :: IPVT(:)
	REAL, intent(inout)  :: ABD(:,:)
      END SUBROUTINE SGBFA

      SUBROUTINE SGBSL( this, ABD, N, ML, MU, IPVT, B, JOB )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in) ::  N, ML, MU, JOB
        INTEGER, intent(in) ::  IPVT(:)
        REAL, intent(in)    ::  ABD(:,:)
        REAL, intent(inout) ::  B(:)
      END SUBROUTINE SGBSL

      SUBROUTINE SGECO( this, A, N,IPVT, RCOND, Z )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in)  :: N
        INTEGER, intent(out) :: IPVT(:)
	REAL, intent(out)    :: RCOND
        REAL, intent(inout)  :: A(:,:)
	REAL, intent(out)    :: Z(:)
      END SUBROUTINE SGECO

      SUBROUTINE SGEFA( this, A, N, IPVT, INFO )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in)  :: N
        INTEGER, intent(out) :: INFO
        INTEGER, intent(out) :: IPVT(:)
        REAL, intent(inout)  :: A(:,:)
      END SUBROUTINE SGEFA

      SUBROUTINE SGESL( this, A, N, IPVT, B, JOB )
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(inout) :: this
        INTEGER, intent(in) :: N, JOB
        INTEGER, intent(in) :: IPVT(:)
        REAL, intent(in)    :: A(:,:)
        REAL, intent(inout) :: B(:)
      END SUBROUTINE SGESL

      FUNCTION tridiag(this, a, b, c, r) result(u)
        import abs_linalgebra_t
        class(abs_linalgebra_t), intent(in) :: this
        REAL, intent(in) :: a(:), b(:), c(:), r(:)
        REAL :: u(size(b))
      END FUNCTION tridiag
   end interface

   contains

   INTEGER FUNCTION ISAMAX( N, SX, INCX )

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR OF INTEREST
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- ISAMAX   FIRST I, I = 1 TO N, TO MAXIMIZE
!                         ABS(SX(1+(I-1)*INCX))

        INTEGER, intent(in) :: N, INCX
	REAL, intent(in)    :: SX(:)

        INTEGER :: I,II
	REAL :: SMAX, XMAG


	IF( N.LE.0 ) THEN
	   ISAMAX = 0
	ELSE IF( N.EQ.1 ) THEN
	   ISAMAX = 1
	ELSE
	   SMAX = 0.0
	   II = 1
	   DO I = 1, 1+(N-1)*INCX, INCX
	      XMAG = ABS(SX(I))
	      IF( SMAX.LT.XMAG ) THEN
	         SMAX = XMAG
	         ISAMAX = II
	      ENDIF
	      II = II + 1
           ENDDO
	ENDIF

   END FUNCTION  ISAMAX

   REAL FUNCTION SASUM( N, SX, INCX )

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR TO BE SUMMED
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'

! --OUTPUT-- SASUM   SUM FROM 0 TO N-1 OF  ABS(SX(1+I*INCX))

        INTEGER, intent(in) :: N, INCX
	REAL, intent(in)    :: SX(:)

        INTEGER :: I, M

	SASUM = 0.0
	IF( N.LE.0 )  RETURN
	IF( INCX.NE.1 ) THEN
!                                          ** NON-UNIT INCREMENTS
	    DO 10 I = 1, 1+(N-1)*INCX, INCX
	       SASUM = SASUM + ABS(SX(I))
   10     CONTINUE
	ELSE
!                                          ** UNIT INCREMENTS
	   M = MOD(N,6)
	   IF( M.NE.0 ) THEN
!                             ** CLEAN-UP LOOP SO REMAINING VECTOR 
!                             ** LENGTH IS A MULTIPLE OF 6.
	      DO 30  I = 1, M
	        SASUM = SASUM + ABS(SX(I))
   30       CONTINUE
	   ENDIF
!                              ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 6
	     SASUM = SASUM + ABS(SX(I))   + ABS(SX(I+1)) + ABS(SX(I+2)) &
                           + ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50    CONTINUE
	ENDIF

   END FUNCTION SASUM

   REAL FUNCTION SDOT( N, SX, INCX, SY, INCY )

!          S.P. DOT PRODUCT OF VECTORS  'X'  AND  'Y'

!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'

! --OUTPUT--
!     SDOT   SUM FOR I = 0 TO N-1 OF  SX(LX+I*INCX) * SY(LY+I*INCY),
!            WHERE  LX = 1          IF INCX .GE. 0, 
!                      = (-INCX)*N  IF INCX .LT. 0,
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

        INTEGER, intent(in) :: N, INCX, INCY
	REAL, intent(in)    :: SX(:), SY(:)

        INTEGER :: I, M, IX, IY

	SDOT = 0.0
	IF( N.LE.0 )  RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SDOT = SDOT + SX(I) * SY(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

!                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,5)
	   IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	         SDOT = SDOT + SX(I) * SY(I)
   20       CONTINUE
	   ENDIF
!                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 5
	      SDOT = SDOT + SX(I)*SY(I) + SX(I+1)*SY(I+1) &
                          + SX(I+2)*SY(I+2) + SX(I+3)*SY(I+3) &
                          + SX(I+4)*SY(I+4)
   30    CONTINUE

	ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SDOT = SDOT + SX(IX) * SY(IY)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

   END FUNCTION SDOT

   SUBROUTINE SSCAL( N, SA, SX, INCX )

!         CALCULATE  X = A*X  (X = VECTOR, A = SCALAR)

!  --INPUT--  N  NUMBER OF ELEMENTS IN VECTOR
!            SA  SINGLE PRECISION SCALE FACTOR
!            SX  SING-PREC ARRAY, LENGTH 1+(N-1)*INCX, CONTAINING VECTOR
!          INCX  SPACING OF VECTOR ELEMENTS IN 'SX'
! --OUTPUT-- SX  REPLACE  SX(1+I*INCX)  WITH  SA * SX(1+I*INCX) 
!                FOR I = 0 TO N-1

        INTEGER, intent(in) :: N, INCX
	REAL, intent(in)    :: SA
	REAL, intent(inout) :: SX(:)

        INTEGER :: I, M

	IF( N.LE.0 ) RETURN

	IF( INCX.NE.1 ) THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SX(I) = SA * SX(I)
   10     CONTINUE

	ELSE

	   M = MOD(N,5)
	   IF( M.NE.0 ) THEN
!                           ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                           ** IS A MULTIPLE OF 5.
	      DO 30  I = 1, M
	         SX(I) = SA * SX(I)
   30       CONTINUE
	   ENDIF
!                             ** UNROLL LOOP FOR SPEED
	   DO 50  I = M+1, N, 5
	      SX(I)   = SA * SX(I)
	      SX(I+1) = SA * SX(I+1)
	      SX(I+2) = SA * SX(I+2)
	      SX(I+3) = SA * SX(I+3)
	      SX(I+4) = SA * SX(I+4)
   50    CONTINUE

	ENDIF

        END SUBROUTINE SSCAL

   SUBROUTINE SAXPY( N, SA, SX, INCX, SY, INCY )

!          Y = A*X + Y  (X, Y = VECTORS, A = SCALAR)
!  --INPUT--
!        N  NUMBER OF ELEMENTS IN INPUT VECTORS 'X' AND 'Y'
!       SA  SINGLE PRECISION SCALAR MULTIPLIER 'A'
!       SX  SING-PREC ARRAY CONTAINING VECTOR 'X'
!     INCX  SPACING OF ELEMENTS OF VECTOR 'X' IN 'SX'
!       SY  SING-PREC ARRAY CONTAINING VECTOR 'Y'
!     INCY  SPACING OF ELEMENTS OF VECTOR 'Y' IN 'SY'
! --OUTPUT--
!       SY   FOR I = 0 TO N-1, OVERWRITE  SY(LY+I*INCY) WITH 
!                 SA*SX(LX+I*INCX) + SY(LY+I*INCY), 
!            WHERE LX = 1          IF INCX .GE. 0,
!                     = (-INCX)*N  IF INCX .LT. 0
!            AND LY IS DEFINED IN A SIMILAR WAY USING INCY.

        INTEGER, intent(in) :: N, INCX, INCY
	REAL, intent(in)    :: SA
	REAL, intent(in)    :: SX(:)
	REAL, intent(inout) :: SY(:)

        INTEGER :: I, M, IX, IY

	IF( N.LE.0 .OR. SA.EQ.0.0 ) RETURN

	IF ( INCX.EQ.INCY .AND. INCX.GT.1 )  THEN

	    DO 10  I = 1, 1+(N-1)*INCX, INCX
	       SY(I) = SY(I) + SA * SX(I)
   10     CONTINUE

	ELSE IF ( INCX.EQ.INCY .AND. INCX.EQ.1 )  THEN

!                                        ** EQUAL, UNIT INCREMENTS
	   M = MOD(N,4)
	   IF( M .NE. 0 ) THEN
!                            ** CLEAN-UP LOOP SO REMAINING VECTOR LENGTH
!                            ** IS A MULTIPLE OF 4.
	      DO 20  I = 1, M
	        SY(I) = SY(I) + SA * SX(I)
   20       CONTINUE
	   ENDIF
!                              ** UNROLL LOOP FOR SPEED
	   DO 30  I = M+1, N, 4
	      SY(I)   = SY(I)   + SA * SX(I)
	      SY(I+1) = SY(I+1) + SA * SX(I+1)
	      SY(I+2) = SY(I+2) + SA * SX(I+2)
	      SY(I+3) = SY(I+3) + SA * SX(I+3)
   30    CONTINUE

	ELSE
!               ** NONEQUAL OR NONPOSITIVE INCREMENTS.
	   IX = 1
	   IY = 1
	   IF( INCX.LT.0 )  IX = 1 + (N-1)*(-INCX)
	   IF( INCY.LT.0 )  IY = 1 + (N-1)*(-INCY)
	   DO 40  I = 1, N
	      SY(IY) = SY(IY) + SA*SX(IX)
	      IX = IX + INCX
	      IY = IY + INCY
   40    CONTINUE

	ENDIF

        END SUBROUTINE SAXPY

   end module abs_linalgebra
