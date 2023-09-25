! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_linear_algebra_linpack

   use musica_constants,               only : dk => musica_dk
   use tuvx_linear_algebra,            only : linear_algebra_t

   implicit none

   public :: linear_algebra_linpack_t

   type, extends(linear_algebra_t) :: linear_algebra_linpack_t
     ! Linear algebra functions using LINPACK
   contains
     procedure :: SGBCO
     procedure :: SGBFA
     procedure :: SGBSL
     procedure :: SGECO
     procedure :: SGEFA
     procedure :: SGESL
     procedure :: tridiag
   end type linear_algebra_linpack_t

   real(dk), parameter ::    rZERO = 0.0_dk
   real(dk), parameter ::    rONE  = 1.0_dk

   contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBCO( this, ABD, N, ML, MU, IPVT, RCOND, Z )
    ! FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
    ! AND ESTIMATES THE CONDITION OF THE MATRIX.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! IF  RCOND  IS NOT NEEDED, SGBFA IS SLIGHTLY FASTER.
    ! TO SOLVE  A*X = B , FOLLOW SBGCO BY SGBSL.

    use tuvx_linear_algebra,           only : SASUM, SDOT, SASUM, SAXPY, SSCAL

    class(linear_algebra_linpack_t), intent(inout) :: this

    real(dk), intent(inout) :: ABD(:,:) ! Banded matrix A(Ax = b) in a packed format
    integer,  intent(in)    :: N        ! Order of the original matrix
    integer,  intent(in)    :: ML       ! Number of diagonals below the main diagonal
    integer,  intent(in)    :: MU       ! Number of diagonals above the main diagonal
    integer,  intent(out)   :: IPVT(:)  ! Pivot vector
    real(dk), intent(out)   :: RCOND    ! Estimate of the condition of the inverse of ABD
    real(dk), intent(out)   :: Z(:)     ! Working vector

    real(dk) :: EK, T, WK, WKM
    real(dk) :: ANORM, S, SM, YNORM
    integer  :: IS, INFO, J, JU, K, KB, KP1, L, LA, LM, LZ, M, MM

    ! COMPUTE 1-NORM OF A
    ANORM = rZERO
    L = ML + 1
    IS = L + MU
    do J = 1, N
      ANORM = max( ANORM, SASUM( L, ABD( IS :, J ), 1 ) )
      if (IS .GT. ML + 1) IS = IS - 1
      if (J .LE. MU) L = L + 1
      if (J .GE. N - ML) L = L - 1
    enddo

    call this%SGBFA( ABD, N, ML, MU, IPVT, INFO )

    ! RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
    ! ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
    ! TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
    ! CHOSEN TO CAUSE maxIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
    ! TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
    ! OVERFLOW.
    ! SOLVE TRANS(U)*W = E
    EK = rONE
    Z( : N ) = rZERO

    M = ML + MU + 1
    JU = 0
    do K = 1, N
      if( Z( K ) /= rZERO ) EK = SIGN( EK, -Z( K ) )
      if( abs( EK - Z( K )) > abs( ABD( M, K ) ) ) then
        S = abs( ABD( M, K ) ) / abs( EK - Z( K ) )
        call SSCAL( N, S, Z, 1 )
        EK = S * EK
      endif
      WK = EK - Z( K )
      WKM = -EK - Z( K )
      S = abs( WK )
      SM = abs( WKM )
      if( ABD( M, K ) /= rZERO ) then
        WK  = WK / ABD( M, K )
        WKM = WKM / ABD( M, K )
      else
        WK  = rONE
        WKM = rONE
      endif
      KP1 = K + 1
      JU = min( max( JU, MU + IPVT( K ) ), N )
      MM = M
      if( KP1 <= JU ) then
        do J = KP1, JU
          MM = MM - 1
          SM = SM + abs( Z( J ) + WKM * ABD( MM, J ) )
          Z( J ) = Z( J ) + WK * ABD( MM, J )
          S = S + abs( Z( J ) )
        enddo
        if( S < SM ) then
          T = WKM - WK
          WK = WKM
          MM = M
          do J = KP1, JU
            MM = MM - 1
            Z( J ) = Z( J ) + T * ABD( MM, J )
          enddo
        endif
      endif
      Z(K) = WK
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )

    ! SOLVE TRANS(L)*Y = W
    do KB = 1, N
      K = N + 1 - KB
      LM = min( ML, N - K )
      if( K < N ) then
        Z( K ) = Z( K ) + SDOT( LM, ABD( M + 1 :, K ), 1, Z( K + 1 : ), 1 )
      endif
      if( abs( Z( K ) ) > rONE ) then
        S = rONE / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
      endif
      L = IPVT( K )
      T = Z( L )
      Z( L ) = Z( K )
      Z( K ) = T
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )

    YNORM = rONE
    ! SOLVE L*V = Y
    do K = 1, N
      L = IPVT( K )
      T = Z( L )
      Z( L ) = Z( K )
      Z( K ) = T
      LM = min( ML, N - K )
      if( K < N ) then
        call SAXPY( LM, T, ABD( M + 1 :, K ), 1, Z( K + 1 : ), 1 )
      endif
      if( abs( Z( K ) ) > rONE ) then
        S = rONE / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
        YNORM = S * YNORM
      endif
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    YNORM = S * YNORM

    ! SOLVE  U*Z = W
    do KB = 1, N
      K = N + 1 - KB
      if( abs( Z( K ) ) > abs( ABD( M, K ) ) ) then
        S = abs( ABD( M, K ) ) / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
        YNORM = S * YNORM
      endif
      if( ABD( M, K ) /= rZERO ) then
        Z( K ) = Z( K ) / ABD( M, K )
      else
        Z( K ) = rONE
      endif
      LM = min( K, M ) - 1
      LA = M - LM
      LZ = K - LM
      T = -Z( K )
      call SAXPY( LM, T, ABD( LA :, K ), 1, Z( LZ : ), 1 )
    enddo
    ! MAKE ZNORM = 1.0
    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    YNORM = S * YNORM

    if( ANORM /= rZERO ) then
      RCOND = YNORM / ANORM
    else
      RCOND = rZERO
    endif

  end subroutine SGBCO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBFA( this, ABD, N, ML, MU, IPVT, INFO )
    ! FACTORS A REAL BAND MATRIX BY ELIMINATION.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! SGBFA IS USUALLY CALLED BY SBGCO, BUT IT CAN BE CALLED
    ! DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
    !
    !
    ! ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, isamax
    ! FROM FORTRAN: max, min

    use tuvx_linear_algebra, only : SAXPY, SSCAL, isamax

    class(linear_algebra_linpack_t), intent(inout) :: this
    real(dk), intent(inout) ::  ABD(:,:) ! Banded matrix A(Ax = b) in a packed format
    integer,  intent(in)    ::  N        ! Order of the original matrix
    integer,  intent(in)    ::  ML       ! Number of diagonals below the main diagonal
    integer,  intent(in)    ::  MU       ! Number of diagonals above the main diagonal
    integer,  intent(out)   ::  INFO     ! = 0  NORMAL VALUE.
    ! = K  IF  U(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
    ! CONDITION FOR THIS SUBROUTINE, BUT IT DOES
    ! INDICATE THAT SGBSL WILL DIVIDE BY ZERO IF
    ! CALLED.  USE  RCOND  IN SBGCO FOR A RELIABLE
    ! INDICATION OF SINGULARITY.
    integer,  intent(out)   ::  IPVT(:)  ! Pivot vector

    real(dk) :: T
    integer  :: I,I0,J,JU,JZ,J0,J1,K,KP1,L,LM,M,MM,NM1

    M = ML + MU + 1
    INFO = 0
    ! ** ZERO INITIAL FILL-IN COLUMNS
    J0 = MU + 2
    J1 = min( N, M ) - 1
    do JZ = J0, J1
      I0 = M + 1 - JZ
      do I = I0, ML
        ABD( I, JZ ) = rZERO
      enddo
    enddo
    JZ = J1
    JU = 0

    ! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
    NM1 = N - 1
    do K = 1, NM1
      KP1 = K + 1
      ! ZERO NEXT FILL-IN COLUMN
      JZ = JZ + 1
      if( JZ <= N ) then
        do I = 1, ML
          ABD( I, JZ ) = rZERO
        enddo
      endif
      ! FIND L = PIVOT INDEX
      LM = min( ML, N - K )
      L = isamax( LM + 1, ABD( M :, K ), 1 ) + M - 1
      IPVT( K ) = L + K - M

      if (ABD(L,K) == rZERO) then
        ! ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
        INFO = K
      else
        ! INTERCHANGE if NECESSARY
        if( L /= M ) then
          T = ABD( L, K )
          ABD( L, K ) = ABD( M, K )
          ABD( M, K ) = T
        endif
        ! COMPUTE MULTIPLIERS
        T = -rONE / ABD( M, K )
        call SSCAL( LM, T, ABD( M + 1 :, K ), 1 )
        ! ROW ELIMINATION WITH COLUMN INDEXING
        JU = min( max( JU, MU + IPVT( K ) ), N )
        MM = M
        do J = KP1, JU
          L = L - 1
          MM = MM - 1
          T = ABD( L, J )
          if( L /= MM ) then
            ABD( L, J ) = ABD( MM, J )
            ABD( MM, J ) = T
          endif
          call SAXPY( LM, T, ABD( M + 1 :, K ), 1, ABD( MM + 1 :,J ), 1 )
        enddo
      endif
    enddo

    IPVT( N ) = N
    if( ABD( M, N ) == rZERO ) INFO = N

  end subroutine SGBFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBSL( this, ABD, N, ML, MU, IPVT, B, JOB )
    ! SOLVES THE REAL BAND SYSTEM
    ! A * X = B  OR  TRANSPOSE(A) * X = B
    ! USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! ERROR CONDITION
    !
    ! A DIVISION BY ZERO WILL OCCUR if THE INPUT FACTOR CONTAINS A
    ! ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
    ! BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
    ! SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
    ! CALLED CORRECTLY AND IF SBGCO HAS SET RCOND .GT. 0.0
    ! OR SGBFA HAS SET INFO .EQ. 0 .
    !
    ! TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
    ! WITH  P  COLUMNS
    ! call SGBCO(ABD,LDA,N,ML,MU,IPVT,RCOND,Z)
    ! if (RCOND IS TOO SMALL) GO TO ...
    ! do 10 J = 1, P
    ! call SGBSL(ABD,LDA,N,ML,MU,IPVT,C(1,J),0)
    ! 10 CONTINUE
    !
    ! ROUTINES CALLED:  FROM BLAS:    SAXPY, SDOT
    ! FROM FORTRAN: min

    use tuvx_linear_algebra, only : SAXPY, SDOT

    class(linear_algebra_linpack_t), intent(inout) :: this
    real(dk), intent(in)    ::  ABD(:,:) ! Output from SGBCO or SGBFA
    integer,  intent(in)    ::  N        ! Order of the original matrix
    integer,  intent(in)    ::  ML       ! Number of diagonals below the main diagonal
    integer,  intent(in)    ::  MU       ! Number of diagonals above the main diagonal
    integer,  intent(in)    ::  IPVT(:)  ! Pivot vector from SGBCO or SGBFA
    real(dk), intent(inout) ::  B(:)     ! Right-hand side vector
    integer,  intent(in)    ::  JOB      ! = 0 TO SOLVE  A*X = B ,
    ! = NONZERO  TO SOLVE  TRANS(A)*X = B , WHERE TRANS(A)  IS THE TRANSPOSE.

    real(dk)     :: T
    integer  :: K,KB,L,LA,LB,LM,M,NM1

    M = MU + ML + 1
    NM1 = N - 1
    if( JOB == 0 ) then
      ! JOB = 0 , SOLVE  A * X = B
      ! FIRST SOLVE L*Y = B
      if( ML /= 0 ) then
        do K = 1, NM1
          LM = min( ML, N - K )
          L = IPVT( K )
          T = B( L )
          if( L /= K ) then
            B( L ) = B( K )
            B( K ) = T
          endif
          call SAXPY( LM, T, ABD( M + 1 :, K ), 1, B( K + 1 : ), 1 )
        enddo
      endif
      ! NOW SOLVE  U*X = Y
      do KB = 1, N
        K = N + 1 - KB
        B( K ) = B( K ) / ABD( M, K )
        LM = min( K, M ) - 1
        LA = M - LM
        LB = K - LM
        T = -B( K )
        call SAXPY( LM, T, ABD( LA :, K ), 1, B( LB : ), 1 )
      enddo
    else
      ! JOB = NONZERO, SOLVE  TRANS(A) * X = B
      ! FIRST SOLVE  TRANS(U)*Y = B
      do K = 1, N
        LM = min( K, M ) - 1
        LA = M - LM
        LB = K - LM
        T = SDOT( LM, ABD( LA :, K ), 1, B( LB : ), 1 )
        B( K ) = ( B( K ) - T ) / ABD( M, K )
      enddo
      ! NOW SOLVE TRANS(L)*X = Y
      if( ML /= 0 ) then
        do KB = 1, NM1
          K = N - KB
          LM = min( ML, N - K )
          B( K ) = B( K ) + SDOT( LM, ABD( M + 1 :,K ), 1, B( K + 1 : ), 1 )
          L = IPVT( K )
          if( L /= K ) then
            T = B( L )
            B( L ) = B( K )
            B( K ) = T
          endif
        enddo
      endif
    endif

  end subroutine SGBSL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGESL( this, A, N, IPVT, B, JOB )
    ! SOLVES THE real SYSTEM
    ! A * X = B  OR  TRANS(A) * X = B
    ! USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! ERROR CONDITION
    !
    ! A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
    ! ZERO ON THE DIAGONAL.  TECHNICALLY, THIS INDICATES SINGULARITY,
    ! BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
    ! SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
    ! CALLED CORRECTLY AND IF SGECO HAS SET RCOND .GT. 0.0
    ! OR SGEFA HAS SET INFO .EQ. 0 .
    !
    ! TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
    ! WITH  P  COLUMNS
    ! call SGECO(A,LDA,N,IPVT,RCOND,Z)
    ! if (RCOND IS TOO SMALL) GO TO ...
    ! do 10 J = 1, P
    ! call SGESL(A,LDA,N,IPVT,C(1,J),0)
    ! 10 CONTINUE

    use tuvx_linear_algebra, only : SAXPY, SDOT

    class(linear_algebra_linpack_t), intent(inout) :: this

    integer,  intent(in)    :: N       ! Order of matrix A
    integer,  intent(in)    :: JOB     ! = 0  TO SOLVE  A*X = B ,
    ! = NONZERO  TO SOLVE TRANS(A)*X = B  WHERE TRANS(A)  IS THE TRANSPOSE.
    integer,  intent(in)    :: IPVT(:) ! Pivot vector from SGECO or SGERFA
    real(dk), intent(in)    :: A(:,:)  ! Output from SGECO or SGEFA
    real(dk), intent(inout) :: B(:)    ! Right-hand side vector

    real(dk)     :: T
    integer  :: K,KB,L,NM1

    NM1 = N - 1
    if( JOB == 0 ) then
      ! JOB = 0 , SOLVE  A * X = B
      ! FIRST SOLVE  L*Y = B
      do K = 1, NM1
        L = IPVT( K )
        T = B( L )
        if( L /= K ) then
           B( L ) = B( K )
           B( K ) = T
        endif
        call SAXPY( N-K, T, A(K+1:,K), 1, B(K+1:), 1 )
      enddo
      ! NOW SOLVE  U*X = Y
      do KB = 1, N
        K = N + 1 - KB
        B( K ) = B( K ) / A( K, K )
        T = -B( K )
        call SAXPY( K - 1, T, A( 1 :, K ), 1, B, 1 )
      enddo
    else
      ! JOB = NONZERO, SOLVE  TRANS(A) * X = B
      ! FIRST SOLVE  TRANS(U)*Y = B
      do K = 1, N
        T = SDOT( K - 1, A( 1 :, K ), 1, B, 1 )
        B( K ) = ( B( K ) - T ) / A( K, K )
      enddo
      ! NOW SOLVE  TRANS(L)*X = Y
      do KB = 1, NM1
        K = N - KB
        B( K ) = B( K ) + SDOT( N - K, A( K + 1 :, K ), 1, B( K + 1 : ), 1 )
        L = IPVT( K )
        if( L /= K ) then
          T = B( L )
          B( L ) = B( K )
          B( K ) = T
        endif
      enddo
    endif

  end subroutine SGESL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGECO( this, A, N, IPVT, RCOND, Z )
    ! FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
    ! AND ESTIMATES THE CONDITION OF THE MATRIX.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! IF  RCOND  IS NOT NEEDED, SGEFA IS SLIGHTLY FASTER.
    ! TO SOLVE  A*X = B , FOLLOW SGECO BY SGESL.
    !
    ! ROUTINES CALLED:  FROM LINPACK: SGEFA
    ! FROM BLAS:    SAXPY, SDOT, SSCAL, SASUM
    ! FROM FORTRAN: abs, Amax1, SIGN

    use tuvx_linear_algebra, only : SAXPY, SSCAL, SDOT, SASUM

    class(linear_algebra_linpack_t), intent(inout) :: this
    real(dk), intent(inout) :: A(:,:) ! AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
    ! WHICH WERE USED TO OBTAIN IT.
    ! THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
    ! L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
    ! TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
    integer,  intent(in)    :: N       ! Order of matrix A
    integer,  intent(out)   :: IPVT(:) ! Pivot vector
    real(dk), intent(out)   :: RCOND   ! AN ESTIMATE OF THE RECIPROCAL CONDITION OF  A .
    ! FOR THE SYSTEM  A*X = B , RELATIVE PERTURBATIONS
    ! IN  A  AND  B  OF SIZE  EPSILON  MAY CAUSE
    ! RELATIVE PERTURBATIONS IN  X  OF SIZE  EPSILON/RCOND .
    ! IF  RCOND  IS SO SMALL THAT THE LOGICAL EXPRESSION
    ! 1.0 + RCOND .EQ. 1.0
    ! IS TRUE, THEN  A  MAY BE SINGULAR TO WORKING
    ! PRECISION.  IN PARTICULAR,  RCOND  IS ZERO  IF
    ! EXACT SINGULARITY IS DETECTED OR THE ESTIMATE
    ! UNDERFLOWS.
    real(dk), intent(out)   :: Z(:) ! A WORK VECTOR WHOSE CONTENTS ARE USUALLY UNIMPORTANT.
    ! IF  A  IS CLOSE TO A SINGULAR MATRIX, THEN  Z  IS
    ! AN APPROXIMATE NULL VECTOR IN THE SENSE THAT
    ! NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .

    real(dk) ::  EK,T,WK,WKM
    real(dk) ::  ANORM,S,SM,YNORM
    integer  :: INFO,J,K,KB,KP1,L

    ! COMPUTE 1-NORM OF A
    ANORM = rZERO
    do J = 1, N
      ANORM = max( ANORM, SASUM(N,A(1:,J),1) )
    enddo
    ! FACTOR
    call this%SGEFA(A,N,IPVT,INFO)

    ! RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
    ! ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  TRANS(A)*Y = E .
    ! TRANS(A)  IS THE TRANSPOSE OF A .  THE COMPONENTS OF  E  ARE
    ! CHOSEN TO CAUSE maxIMUM LOCAL GROWTH IN THE ELEMENTS OF W  WHERE
    ! TRANS(U)*W = E .  THE VECTORS ARE FREQUENTLY RESCALED TO AVOID
    ! OVERFLOW.

    ! ** SOLVE TRANS(U)*W = E
    EK = rONE
    do J = 1, N
      Z( J ) = rZERO
    enddo

    do K = 1, N
      if( Z( K ) .NE. rZERO ) EK = SIGN( EK, -Z( K ) )
      if( abs( EK - Z( K ) ) .GT. abs( A( K, K ) ) ) then
        S = abs( A( K, K ) ) / abs( EK - Z( K ) )
        call SSCAL( N, S, Z, 1 )
        EK = S * EK
      endif
      WK = EK - Z( K )
      WKM = -EK - Z( K )
      S = abs( WK )
      SM = abs( WKM )
      if( A( K, K ) .NE. rZERO ) then
        WK  = WK  / A( K, K )
        WKM = WKM / A( K, K )
      else
        WK  = rONE
        WKM = rONE
      endif
      KP1 = K + 1
      if( KP1 .LE. N ) then
        do J = KP1, N
          SM = SM + abs( Z( J ) + WKM * A( K, J ) )
          Z( J ) = Z( J ) + WK * A( K, J )
          S = S + abs( Z( J ) )
        enddo
        if( S .LT. SM ) then
          T = WKM - WK
          WK = WKM
          do J = KP1, N
            Z( J ) = Z( J ) + T * A( K, J )
          enddo
        endif
      endif
      Z( K ) = WK
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    ! SOLVE TRANS(L)*Y = W
    do KB = 1, N
      K = N + 1 - KB
      if( K .LT. N ) then
        Z( K ) = Z( K ) + SDOT( N - K, A( K + 1 :, K ), 1, Z( K + 1 : ), 1 )
      endif
      if( abs( Z( K ) ) .GT. rONE ) then
        S = rONE / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
      endif
      L = IPVT( K )
      T = Z( L )
      Z( L ) = Z( K )
      Z( K ) = T
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    ! SOLVE L*V = Y
    YNORM = rONE
    do K = 1, N
      L = IPVT( K )
      T = Z( L )
      Z( L ) = Z( K )
      Z( K ) = T
      if( K .LT. N ) call SAXPY( N - K, T, A( K + 1 :, K ), 1, Z( K + 1 : ), 1 )
      if( abs( Z( K ) ) .GT. rONE ) then
        S = rONE / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
        YNORM = S * YNORM
      endif
    enddo

    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    ! SOLVE  U*Z = V
    YNORM = S * YNORM
    do KB = 1, N
      K = N + 1 - KB
      if( abs( Z( K ) ) .GT. abs( A( K, K ) ) ) then
        S = abs( A( K, K ) ) / abs( Z( K ) )
        call SSCAL( N, S, Z, 1 )
        YNORM = S * YNORM
      endif
      if( A( K, K ) .NE. rZERO ) Z( K ) = Z( K ) / A( K, K )
      if( A( K, K ) .EQ. rZERO ) Z( K ) = 1.0E0
      T = -Z( K )
      call SAXPY( K - 1, T, A( 1 :, K ), 1, Z, 1 )
    enddo
    ! MAKE ZNORM = 1.0
    S = rONE / SASUM( N, Z, 1 )
    call SSCAL( N, S, Z, 1 )
    YNORM = S * YNORM

    if( ANORM .NE. rZERO ) RCOND = YNORM / ANORM
    if( ANORM .EQ. rZERO ) RCOND = rZERO

  end subroutine SGECO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGEFA( this, A, N, IPVT, INFO )
    ! FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
    !
    ! REVISION DATE:  8/1/82
    ! AUTHOR:  MOLER, C. B., (U. OF NEW MEXICO)
    !
    ! SGEFA IS USUALLY CALLED BY SGECO, BUT IT CAN BE CALLED
    ! DIRECTLY WITH A SAVING IN TIME IF RCOND IS NOT NEEDED.
    ! (TIME FOR SGECO) = (1 + 9/N)*(TIME FOR SGEFA) .
    !
    ! ROUTINES CALLED:  FROM BLAS:    SAXPY, SSCAL, isamax

    use tuvx_linear_algebra, only : SAXPY, SSCAL, isamax

    class(linear_algebra_linpack_t), intent(inout) :: this
    real(dk), intent(inout) :: A(:,:) ! AN UPPER TRIANGULAR MATRIX AND THE MULTIPLIERS
    ! WHICH WERE USED TO OBTAIN IT.
    ! THE FACTORIZATION CAN BE WRITTEN  A = L*U , WHERE
    ! L  IS A PRODUCT OF PERMUTATION AND UNIT LOWER
    ! TRIANGULAR MATRICES AND  U  IS UPPER TRIANGULAR.
    integer,  intent(in)    :: N       ! Order of matrix A
    integer,  intent(out)   :: IPVT(:) ! Pivot vector
    integer,  intent(out)   ::  INFO   ! = 0  NORMAL VALUE.
    ! = K  REALU(K,K) .EQ. 0.0 .  THIS IS NOT AN ERROR
    ! CONDITION FOR THIS SUBROUTINE, BUT IT doES
    ! INDICATE THAT SGESL OR SGEDI WILL DIVIDE BY ZERO
    ! if CALLED.  USE  RCOND  IN SGECO FOR A RELIABLE
    ! INDICATION OF SINGULARITY.

    real(dk)     :: T
    integer  :: J,K,KP1,L,NM1

    ! GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
    INFO = 0
    NM1 = N - 1
    do K = 1, NM1
      KP1 = K + 1
      ! FIND L = PIVOT INDEX
      L = isamax( N - K + 1, A( K :, K ), 1 ) + K - 1
      IPVT( K ) = L
      if( A( L, K ) .EQ. rZERO ) then
        ! ZERO PIVOT IMPLIES THIS COLUMN ALREADY TRIANGULARIZED
        INFO = K
      else
        ! INTERCHANGE if NECESSARY
        if( L .NE. K ) then
          T = A( L, K )
          A( L, K ) = A( K, K )
          A( K, K ) = T
        endif
        ! COMPUTE MULTIPLIERS
        T = -rONE / A( K, K )
        call SSCAL( N - K, T, A( K + 1 :, K ), 1 )

        ! ** ROW ELIMINATION WITH COLUMN INDEXING
        do J = KP1, N
          T = A( L, J )
          if( L .NE. K ) then
            A( L, J ) = A( K, J )
            A( K, J ) = T
          endif
          call SAXPY( N - K, T, A( K + 1 :, K ), 1, A( K + 1 :, J ), 1 )
        enddo
      endif
    enddo

    IPVT( N ) = N
    if( A( N, N ) == rZERO ) INFO = N

  end subroutine SGEFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function tridiag( this, a, b, c, r ) result( u )
    ! The system to be solved is:
    !
    ! .. math::
    !
    !    \begin{pmatrix}
    !      b_1 & c_1 & 0   & \cdots &    0     &    0    &    0    \\
    !      a_2 & b_2 & c_2 & \cdots &    0     &    0    &    0    \\
    !          &     &     & \cdots &          &         &         \\
    !       0  &  0  &  0  & \cdots & a_{N-1}  & b_{N-1} & c_{N-1} \\
    !       0  &  0  &  0  & \cdots &    0     &   a_N   &   b_N
    !    \end{pmatrix}
    !    \times
    !    \begin{pmatrix}
    !      u_1 \\
    !      u_2 \\
    !      \cdots \\
    !      u_{N-1} \\
    !      u_N
    !    \end{pmatrix}
    !    =
    !    \begin{pmatrix}
    !      r_1 \\
    !      r_2 \\
    !      \cdots \\
    !      r_{N-1} \\
    !      r_N
    !    \end{pmatrix}
    !
    ! ..
    !   | b1 c1  0 ...                |   |  u1  |   |  r1  |
    !   | a2 b2 c2 ...                |   |  u2  |   |  r2  |
    !   |          ...                | * | ...  | = | ...  |
    !   |          ... aN-1 bN-1 cN-1 |   | uN-1 |   | rN-1 |
    !   |          ...   0   aN   bN  |   |  uN  |   |  rN  |
    !
    use musica_assert,    only : assert_msg
    use musica_constants, only : dk => musica_dk

    class(linear_algebra_linpack_t), intent(in) :: this
    real(dk),            intent(in) :: a(:) ! lower diagonal
    real(dk),            intent(in) :: b(:) ! primary diagonal
    real(dk),            intent(in) :: c(:) ! upper diagonal
    real(dk),            intent(in) :: r(:) ! right-hand side vector
    real(dk)                        :: u( size( b ) ) ! result vector

    integer :: i
    real(dk) :: denom
    real(dk) :: cp( size( b ) )

    cp(1) = c(1) / b(1)
    u(1) = r(1) / b(1)
    do i = 2, size( b )
      denom = 1.0 / ( b(i) - a(i) * cp(i-1) )
      cp(i) = c(i) * denom
      u(i) = ( r(i) - a(i) * u(i-1) ) * denom
    end do
    do i = size( b ) - 1, 1, -1
      u(i) = u(i) - cp(i) * u(i+1)
    end do

  end function tridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_linear_algebra_linpack
