      SUBROUTINE LEPOLY( NMU, M, MAXMU, TWONM1, MU, YLM )

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

      INTEGER   MAXSQT
      PARAMETER ( MAXSQT = 1000 )
c     ..
c     .. Scalar Arguments ..

      INTEGER   M, MAXMU, NMU, TWONM1
c     ..
c     .. Array Arguments ..

      REAL      MU( NMU ), YLM( 0:MAXMU, NMU )
c     ..
c     .. Local Scalars ..

      LOGICAL   PASS1
      INTEGER   I, L, NS
      REAL      TMP1, TMP2
c     ..
c     .. Local Arrays ..

      REAL      SQT( MAXSQT )
c     ..
c     .. External Subroutines ..

c     EXTERNAL  ERRMSG
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC FLOAT, SQRT
c     ..
      SAVE      SQT, PASS1
      DATA      PASS1 / .TRUE. /


      IF( PASS1 ) THEN

         PASS1  = .FALSE.

         DO 10 NS = 1, MAXSQT
            SQT( NS ) = SQRT( FLOAT( NS ) )
   10    CONTINUE

      END IF

      IF( 2*TWONM1.GT.MAXSQT )
     &    CALL ERRMSG('LEPOLY--need to increase param MAXSQT',.True.)


      IF( M.EQ.0 ) THEN
c                             ** Upward recurrence for ordinary
c                                Legendre polynomials
         DO 20 I = 1, NMU
            YLM( 0, I ) = rONE
            YLM( 1, I ) = MU( I )
   20    CONTINUE


         DO 40 L = 2, TWONM1

            DO 30 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L - 1, I ) -
     &                         ( L - 1 )*YLM( L - 2, I ) ) / L
   30       CONTINUE

   40    CONTINUE


      ELSE

         DO 50 I = 1, NMU
c                               ** Y-sub-m-super-m; derived from
c                               ** D/A Eqs. (11,12)

            YLM( M, I ) = - SQT( 2*M - 1 ) / SQT( 2*M )*
     &                      SQRT( 1.- MU(I)**2 )*YLM( M - 1, I )

c                              ** Y-sub-(m+1)-super-m; derived from
c                              ** D/A Eqs.(13,14) using Eqs.(11,12)

            YLM( M + 1, I ) = SQT( 2*M + 1 )*MU( I )*YLM( M, I )

   50    CONTINUE

c                                   ** Upward recurrence; D/A EQ.(10)
         DO 70 L = M + 2, TWONM1

            TMP1  = SQT( L - M )*SQT( L + M )
            TMP2  = SQT( L - M - 1 )*SQT( L + M - 1 )

            DO 60 I = 1, NMU
               YLM( L, I ) = ( ( 2*L - 1 )*MU( I )*YLM( L-1, I ) -
     &                         TMP2*YLM( L-2, I ) ) / TMP1
   60       CONTINUE

   70    CONTINUE

      END IF

      END SUBROUTINE LEPOLY
