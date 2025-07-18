
   module linalgebra

   use abs_linalgebra, only : abs_linalgebra_t

   implicit none

   public :: linalgebra_t

   type, extends(abs_linalgebra_t) :: linalgebra_t
     contains
     procedure :: SGBCO
     procedure :: SGBFA
     procedure :: SGBSL
     procedure :: SGECO
     procedure :: SGEFA
     procedure :: SGESL
     procedure :: tridiag
   end type linalgebra_t

   REAL, PARAMETER :: rZERO = 0.0
   REAL, PARAMETER :: rONE  = 1.0

   contains

   SUBROUTINE SGBCO( this, ABD, N, ML, MU, IPVT, RCOND, Z )
!------------------------------------------------------------
!         FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION 
!         AND ESTIMATES THE CONDITION OF THE MATRIX.
!
!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
!
!     IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
!     TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.
!------------------------------------------------------------

   use abs_linalgebra, only : SASUM, SDOT, SASUM, SAXPY, SSCAL

   class(linalgebra_t), intent(inout) :: this

   INTEGER, intent(in)  :: N, ML, MU
   INTEGER, intent(out) :: IPVT(:)
   REAL, intent(out)    :: RCOND
   REAL, intent(inout)  :: ABD(:,:)
   REAL, intent(out)    :: Z(:)

   REAL    :: EK, T, WK, WKM
   REAL    :: ANORM, S, SM, YNORM
   INTEGER :: IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM


!                       ** COMPUTE 1-NORM OF A
   ANORM = rZERO
   L = ML + 1
   IS = L + MU
   DO J = 1, N
     ANORM = MAX(ANORM, SASUM(L,ABD(IS:,J), 1))
     IF (IS .GT. ML + 1) IS = IS - 1
     IF (J .LE. MU) L = L + 1
     IF (J .GE. N - ML) L = L - 1
   ENDDO

   CALL this%SGBFA(ABD, N, ML, MU, IPVT, INFO)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                     ** SOLVE TRANS(U)*W = E
	EK = rONE
	Z(1:N) = rZERO

	M = ML + MU + 1
	JU = 0
	DO K = 1, N
	   IF (Z(K) /= rZERO) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) > ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K))/ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (ABD(M,K) /= rZERO) THEN
	      WK  = WK /ABD(M,K)
	      WKM = WKM/ABD(M,K)
	   ELSE
	      WK  = rONE
	      WKM = rONE
	   ENDIF
	   KP1 = K + 1
	   JU = MIN(MAX(JU, MU+IPVT(K)), N)
	   MM = M
	   IF (KP1 <= JU) THEN
	      DO J = KP1, JU
	         MM = MM - 1
	         SM = SM + ABS(Z(J)+WKM*ABD(MM,J))
	         Z(J) = Z(J) + WK*ABD(MM,J)
	         S = S + ABS(Z(J))
              ENDDO
	      IF (S < SM) THEN
	         T = WKM - WK
	         WK = WKM
	         MM = M
	         DO J = KP1, JU
	            MM = MM - 1
	            Z(J) = Z(J) + T*ABD(MM,J)
                 ENDDO
	      ENDIF
	   ENDIF
	   Z(K) = WK
        ENDDO

	S = rONE / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

!                         ** SOLVE TRANS(L)*Y = W
	DO KB = 1, N
	   K = N + 1 - KB
	   LM = MIN(ML, N-K)
	   IF (K < N) THEN
             Z(K) = Z(K) + SDOT(LM, ABD(M+1:,K), 1, Z(K+1:), 1)
	   ENDIF
	   IF (ABS(Z(K)) > rONE) THEN
	      S = rONE / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	ENDDO

	S = rONE / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)

	YNORM = rONE
!                         ** SOLVE L*V = Y
	DO K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   LM = MIN(ML, N-K)
	   IF (K < N) THEN
             CALL SAXPY(LM, T, ABD(M+1:,K), 1, Z(K+1:), 1)
	   ENDIF
	   IF (ABS(Z(K)) > rONE) THEN
	      S = rONE / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	ENDDO

	S = rONE/SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM
!                           ** SOLVE  U*Z = W
	DO KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) > ABS(ABD(M,K))) THEN
	      S = ABS(ABD(M,K)) / ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (ABD(M,K) /= rZERO) THEN
              Z(K) = Z(K)/ABD(M,K)
           ELSE
	      Z(K) = rONE
           ENDIF
	   LM = MIN(K, M) - 1
	   LA = M - LM
	   LZ = K - LM
	   T = -Z(K)
	   CALL SAXPY(LM, T, ABD(LA:,K), 1, Z(LZ:), 1)
	ENDDO
!                              ** MAKE ZNORM = 1.0
	S = rONE / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM /= rZERO) THEN
           RCOND = YNORM/ANORM
        ELSE
	   RCOND = rZERO
        ENDIF

   END SUBROUTINE SGBCO

   SUBROUTINE SGBFA( this, ABD, N, ML, MU, IPVT, INFO )

   use abs_linalgebra, only : SAXPY, SSCAL, ISAMAX

!         FACTORS A REAL BAND MATRIX BY ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.

!     INPUT:  SAME AS 'SGBCO'

!     ON RETURN:

!        ABD,IPVT    SAME AS 'SGBCO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
!                     CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     (SEE 'SGBCO' FOR DESCRIPTION OF BAND STORAGE MODE)

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX
!                       FROM FORTRAN: MAX, MIN

   class(linalgebra_t), intent(inout) :: this
   INTEGER, intent(in)   ::  N, ML, MU
   INTEGER, intent(out)  ::  INFO
   INTEGER, intent(out)  ::  IPVT(:)
   REAL, intent(inout)   ::  ABD(:,:)

   REAL     :: T
   INTEGER  :: I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1


	M = ML + MU + 1
	INFO = 0
!                        ** ZERO INITIAL FILL-IN COLUMNS
	J0 = MU + 2
	J1 = MIN(N, M) - 1
	DO JZ = J0, J1
	   I0 = M + 1 - JZ
	   DO I = I0, ML
	      ABD(I,JZ) = rZERO
           ENDDO
        ENDDO
	JZ = J1
	JU = 0

!                       ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	NM1 = N - 1
        DO K = 1, NM1
           KP1 = K + 1
!                                  ** ZERO NEXT FILL-IN COLUMN
	   JZ = JZ + 1
	   IF (JZ <= N) THEN
	      DO I = 1, ML
	         ABD(I,JZ) = rZERO
              ENDDO
	   ENDIF
!                                  ** FIND L = PIVOT INDEX
	   LM = MIN(ML, N-K)
	   L = ISAMAX(LM+1, ABD(M:,K), 1) + M - 1
	   IPVT(K) = L + K - M

	   IF (ABD(L,K) == rZERO) THEN
!          ** ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
!                                ** INTERCHANGE IF NECESSARY
	      IF (L /= M) THEN
	         T = ABD(L,K)
	         ABD(L,K) = ABD(M,K)
	         ABD(M,K) = T
	      ENDIF
!                                   ** COMPUTE MULTIPLIERS
	      T = -rONE / ABD(M,K)
	      CALL SSCAL(LM, T, ABD(M+1:,K), 1)

!                               ** ROW ELIMINATION WITH COLUMN INDEXING

	      JU = MIN(MAX(JU, MU+IPVT(K)), N)
	      MM = M
	      DO J = KP1, JU
	         L = L - 1
	         MM = MM - 1
	         T = ABD(L,J)
	         IF (L /= MM) THEN
	            ABD(L,J) = ABD(MM,J)
	            ABD(MM,J) = T
	         ENDIF
	         CALL SAXPY(LM, T, ABD(M+1:,K), 1, ABD(MM+1:,J), 1)
              ENDDO
           ENDIF
        ENDDO

	IPVT(N) = N
	IF (ABD(M,N) == rZERO) INFO = N

   END SUBROUTINE SGBFA

   SUBROUTINE SGBSL( this, ABD, N, ML, MU, IPVT, B, JOB )

   use abs_linalgebra, only : SAXPY, SDOT

!         SOLVES THE REAL BAND SYSTEM
!            A * X = B  OR  TRANSPOSE(A) * X = B
!         USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     INPUT:

!        ABD     REAL(LDA, N)
!                THE OUTPUT FROM SBGCO OR SGBFA.

!        N       INTEGER
!                THE ORDER OF THE ORIGINAL MATRIX.

!        ML      INTEGER
!                NUMBER OF DIAGONALS BELOW THE MAIN DIAGONAL.

!        MU      INTEGER
!                NUMBER OF DIAGONALS ABOVE THE MAIN DIAGONAL.

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SBGCO OR SGBFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B , WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
!        OR SGBFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
!        10 CONTINUE

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
!                       FROM FORTRAN: MIN

   class(linalgebra_t), intent(inout) :: this

   INTEGER, intent(in) ::  N, ML, MU, JOB
   INTEGER, intent(in) ::  IPVT(:)
   REAL, intent(in)    ::  ABD(:,:)
   REAL, intent(inout) ::  B(:)

	REAL     :: T
	INTEGER  :: K,KB,L,LA,LB,LM,M,NM1


        M = MU + ML + 1
        NM1 = N - 1
        IF (JOB == 0) THEN
!                               ** JOB = 0 , SOLVE  A * X = B
!                               ** FIRST SOLVE L*Y = B
	   IF (ML /= 0) THEN
	      DO K = 1, NM1
	         LM = MIN(ML, N-K)
	         L = IPVT(K)
	         T = B(L)
	         IF (L /= K) THEN
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
	         CALL SAXPY( LM, T, ABD(M+1:,K), 1, B(K+1:), 1 )
	      ENDDO
	   ENDIF
!                           ** NOW SOLVE  U*X = Y
	   DO KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / ABD(M,K)
	      LM = MIN(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = -B(K)
	      CALL SAXPY(LM, T, ABD(LA:,K), 1, B(LB:), 1)
	   ENDDO
        ELSE
!                          ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                  ** FIRST SOLVE  TRANS(U)*Y = B
	   DO K = 1, N
	      LM = MIN(K, M) - 1
	      LA = M - LM
	      LB = K - LM
	      T = SDOT(LM, ABD(LA:,K), 1, B(LB:), 1)
	      B(K) = (B(K) - T)/ABD(M,K)
	   ENDDO
!                                  ** NOW SOLVE TRANS(L)*X = Y
	   IF (ML /= 0) THEN
	      DO KB = 1, NM1
	         K = N - KB
	         LM = MIN(ML, N-K)
	         B(K) = B(K) + SDOT(LM, ABD(M+1:,K), 1, B(K+1:), 1)
	         L = IPVT(K)
	         IF (L /= K) THEN
	            T = B(L)
	            B(L) = B(K)
	            B(K) = T
	         ENDIF
	      ENDDO
	   ENDIF
        ENDIF

   END SUBROUTINE SGBSL

   SUBROUTINE SGESL( this, A, N, IPVT, B, JOB )

   use abs_linalgebra, only : SAXPY, SDOT

!         SOLVES THE REAL SYSTEM
!            A * X = B  OR  TRANS(A) * X = B
!         USING THE FACTORS COMPUTED BY SGECO OR SGEFA.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE OUTPUT FROM SGECO OR SGEFA.

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!        IPVT    INTEGER(N)
!                THE PIVOT VECTOR FROM SGECO OR SGEFA.

!        B       REAL(N)
!                THE RIGHT HAND SIDE VECTOR.

!        JOB     INTEGER
!                = 0         TO SOLVE  A*X = B ,
!                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
!                            TRANS(A)  IS THE TRANSPOSE.

!     ON RETURN

!        B       THE SOLUTION VECTOR  X .

!     ERROR CONDITION

!        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
!        ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
!        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
!        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
!        CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
!        OR SGEFA HAS SET INFO .EQ. 0 .

!     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
!     WITH  P  COLUMNS
!           CALL SGECO(A,LDA,N,IPVT,RCOND,Z)
!           IF (RCOND IS TOO SMALL) GO TO ...
!           DO 10 J = 1, P
!              CALL SGESL(A,LDA,N,IPVT,C(1,J),0)
!        10 CONTINUE

   class(linalgebra_t), intent(inout) :: this

        INTEGER, intent(in) :: N, JOB
        INTEGER, intent(in) :: IPVT(:)
	REAL, intent(in)    :: A(:,:)
	REAL, intent(inout) :: B(:)

	REAL     :: T
	INTEGER  :: K,KB,L,NM1


	NM1 = N - 1
	IF (JOB == 0) THEN
!                                 ** JOB = 0 , SOLVE  A * X = B
!                                     ** FIRST SOLVE  L*Y = B
	   DO K = 1, NM1
	      L = IPVT(K)
	      T = B(L)
	      IF (L /= K) THEN
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
	      CALL SAXPY( N-K, T, A(K+1:,K), 1, B(K+1:), 1 )
           ENDDO
!                                    ** NOW SOLVE  U*X = Y
	   DO KB = 1, N
	      K = N + 1 - KB
	      B(K) = B(K) / A(K,K)
	      T = -B(K)
	      CALL SAXPY( K-1, T, A(1:,K), 1, B, 1 )
           ENDDO

	ELSE
!                         ** JOB = NONZERO, SOLVE  TRANS(A) * X = B
!                                    ** FIRST SOLVE  TRANS(U)*Y = B
	   DO K = 1, N
	      T = SDOT( K-1, A(1:,K), 1, B, 1 )
	      B(K) = (B(K) - T) / A(K,K)
           ENDDO
!                                    ** NOW SOLVE  TRANS(L)*X = Y
	   DO KB = 1, NM1
	      K = N - KB
	      B(K) = B(K) + SDOT( N-K, A(K+1:,K), 1, B(K+1:), 1 )
	      L = IPVT(K)
	      IF (L /= K) THEN
	         T = B(L)
	         B(L) = B(K)
	         B(K) = T
	      ENDIF
           ENDDO

	ENDIF

   END SUBROUTINE SGESL

   SUBROUTINE SGECO( this, A, N, IPVT, RCOND, Z )

   use abs_linalgebra, only : SAXPY, SSCAL, SDOT, SASUM

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
!         AND ESTIMATES THE CONDITION OF THE MATRIX.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!         IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
!         TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.

!     ON ENTRY

!        A       REAL(LDA, N)
!                THE MATRIX TO BE FACTORED.

!        N       INTEGER
!                THE ORDER OF THE MATRIX  A .

!     ON RETURN

!        A       AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
!                WHICH WERE USED TO OBTAIN IT.
!                THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
!                L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
!                TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.

!        IPVT    INTEGER(N)
!                AN INTEGER VECTOR OF PIVOT INDICES.

!        RCOND   REAL
!                AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
!                FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
!                IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
!                RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
!                IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
!                           1.0 + RCOND .EQ. 1.0
!                IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
!                PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
!                EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
!                UNDERFLOWS.

!        Z       REAL(N)
!                A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
!                IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
!                AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

!     ROUTINES CALLED:  FROM LINPACK: SGEFA
!                       FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
!                       FROM FORTRAN: ABS, AMAX1, SIGN

   class(linalgebra_t), intent(inout) :: this

   INTEGER, intent(in)    :: N
   INTEGER, intent(out)   :: IPVT(:)
   REAL, intent(out)      :: RCOND
   REAL, intent(inout)    :: A(:,:)
   REAL, intent(out)      :: Z(:)

   REAL ::  EK,T,WK,WKM
   REAL ::  ANORM,S,SM,YNORM
   INTEGER  :: INFO,J,K,KB,KP1,L

!                        ** COMPUTE 1-NORM OF A
	ANORM = 0.0E0
	DO 10 J = 1, N
	   ANORM = MAX( ANORM, SASUM(N,A(1:,J),1) )
   10 CONTINUE
!                                      ** FACTOR
      CALL this%SGEFA(A,N,IPVT,INFO)

!     RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
!     TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
!     CHOSEN TO CAUSE MAXIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
!     TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
!     OVERFLOW.

!                        ** SOLVE TRANS(U)*W = E
	EK = 1.0E0
	DO 20 J = 1, N
	   Z(J) = 0.0E0
   20 CONTINUE

	DO 100 K = 1, N
	   IF (Z(K) .NE. 0.0E0) EK = SIGN(EK, -Z(K))
	   IF (ABS(EK-Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K)) / ABS(EK-Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      EK = S*EK
	   ENDIF
	   WK = EK - Z(K)
	   WKM = -EK - Z(K)
	   S = ABS(WK)
	   SM = ABS(WKM)
	   IF (A(K,K) .NE. 0.0E0) THEN
	      WK  = WK  / A(K,K)
	      WKM = WKM / A(K,K)
	   ELSE
	      WK  = 1.0E0
	      WKM = 1.0E0
	   ENDIF
	   KP1 = K + 1
	   IF (KP1 .LE. N) THEN
	      DO 60 J = KP1, N
	         SM = SM + ABS(Z(J)+WKM*A(K,J))
	         Z(J) = Z(J) + WK*A(K,J)
	         S = S + ABS(Z(J))
   60       CONTINUE
	      IF (S .LT. SM) THEN
	         T = WKM - WK
	         WK = WKM
	         DO 70 J = KP1, N
	            Z(J) = Z(J) + T*A(K,J)
   70          CONTINUE
	      ENDIF
	   ENDIF
	   Z(K) = WK
  100 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
!                                ** SOLVE TRANS(L)*Y = W
	DO 120 KB = 1, N
	   K = N + 1 - KB
	   IF (K .LT. N) THEN
             Z(K) = Z(K) + SDOT(N-K, A(K+1:,K), 1, Z(K+1:), 1)
	   ENDIF
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	   ENDIF
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
  120 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
!                                 ** SOLVE L*V = Y
	YNORM = 1.0E0
	DO 140 K = 1, N
	   L = IPVT(K)
	   T = Z(L)
	   Z(L) = Z(K)
	   Z(K) = T
	   IF (K .LT. N) CALL SAXPY(N-K, T, A(K+1:,K), 1, Z(K+1:), 1)
	   IF (ABS(Z(K)) .GT. 1.0E0) THEN
	      S = 1.0E0/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
  140 CONTINUE

	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
!                                  ** SOLVE  U*Z = V
	YNORM = S*YNORM
	DO 160 KB = 1, N
	   K = N + 1 - KB
	   IF (ABS(Z(K)) .GT. ABS(A(K,K))) THEN
	      S = ABS(A(K,K))/ABS(Z(K))
	      CALL SSCAL(N, S, Z, 1)
	      YNORM = S*YNORM
	   ENDIF
	   IF (A(K,K) .NE. 0.0E0) Z(K) = Z(K)/A(K,K)
	   IF (A(K,K) .EQ. 0.0E0) Z(K) = 1.0E0
	   T = -Z(K)
	   CALL SAXPY(K-1, T, A(1:,K), 1, Z, 1)
  160 CONTINUE
!                                   ** MAKE ZNORM = 1.0
	S = 1.0E0 / SASUM(N, Z, 1)
	CALL SSCAL(N, S, Z, 1)
	YNORM = S*YNORM

	IF (ANORM .NE. 0.0E0) RCOND = YNORM/ANORM
	IF (ANORM .EQ. 0.0E0) RCOND = 0.0E0
	
   END SUBROUTINE SGECO

   SUBROUTINE SGEFA( this, A, N, IPVT, INFO )

   use abs_linalgebra, only : SAXPY, SSCAL, ISAMAX

!         FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.

!         REVISION DATE:  8/1/82
!         AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)

!     SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
!     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
!     (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .

!     INPUT:  SAME AS 'SGECO'

!     ON RETURN:

!        A,IPVT  SAME AS 'SGECO'

!        INFO    INTEGER
!                = 0  NORMAL VALUE.
!                = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
!                     CONDITION FOR THIS SUBROUTINE, BUT IT DOES
!                     INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
!                     IF CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
!                     INDICATION OF SINGULARITY.

!     ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, ISAMAX

   class(linalgebra_t), intent(inout) :: this

   INTEGER, intent(in)  ::  N
   INTEGER, intent(out) ::  INFO
   INTEGER, intent(out) ::  IPVT(:)
   REAL, intent(inout)  ::  A(:,:)

	REAL     :: T
	INTEGER  :: J,K,KP1,L,NM1


!                      ** GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
	INFO = 0
	NM1 = N - 1
	DO 60 K = 1, NM1
	   KP1 = K + 1
!                                            ** FIND L = PIVOT INDEX
	   L = ISAMAX( N-K+1, A(K:,K), 1) + K-1
	   IPVT(K) = L

	   IF (A(L,K) .EQ. 0.0E0) THEN
!                                     ** ZERO PIVOT IMPLIES THIS COLUMN 
!                                     ** ALREADY TRIANGULARIZED
	      INFO = K
	   ELSE
!                                     ** INTERCHANGE IF NECESSARY
	      IF (L .NE. K) THEN
	         T = A(L,K)
	         A(L,K) = A(K,K)
	         A(K,K) = T
	      ENDIF
!                                     ** COMPUTE MULTIPLIERS
	      T = -1.0E0 / A(K,K)
	      CALL SSCAL( N-K, T, A(K+1:,K), 1 )

!                              ** ROW ELIMINATION WITH COLUMN INDEXING
	      DO 30 J = KP1, N
	         T = A(L,J)
	         IF (L .NE. K) THEN
	            A(L,J) = A(K,J)
	            A(K,J) = T
	         ENDIF
	         CALL SAXPY( N-K, T, A(K+1:,K), 1, A(K+1:,J), 1 )
   30       CONTINUE

	   ENDIF

   60 CONTINUE

	IPVT(N) = N
	IF (A(N,N) == rZERO) INFO = N

   END SUBROUTINE SGEFA

   FUNCTION tridiag(this, a, b, c, r) result(u)
!_______________________________________________________________________
! solves tridiagonal system.  From Numerical Recipies, p. 40
!_______________________________________________________________________

! input:
      REAL, intent(in) :: a(:), b(:), c(:), r(:)
      class(linalgebra_t), intent(in) :: this

! output:
      REAL :: u(size(b))

! local:
      INTEGER :: j
      INTEGER :: n

      REAL :: bet
      REAL :: gam(2*size(b))
!_______________________________________________________________________

      IF (b(1) == rZERO) THEN
        write(*,*) 'tridiag: pivot 1 is zero; halting'
        STOP 'tridiag: zero pivot'
      ENDIF
      n = size(b)
      bet   = b(1)
      u(1) = r(1)/bet
      DO j = 2, n   
         gam(j) = c(j - 1)/bet
         bet = b(j) - a(j)*gam(j)
         IF (bet == rZERO) THEN
           write(*,'('' tridiag: pivot '',i4,''is zero; halting'')') j
           STOP 'tridiag: zero pivot'
         ENDIF
         u(j) = (r(j) - a(j)*u(j - 1))/bet
      ENDDO

      DO j = n - 1, 1, -1  
         u(j) = u(j) - gam(j + 1)*u(j + 1)
      ENDDO

      END FUNCTION tridiag

   end module linalgebra
