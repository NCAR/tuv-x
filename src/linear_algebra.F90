! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_linear_algebra
  ! General interface and functions for linear algebra types

   use musica_constants,               only : dk => musica_dk

   implicit none

   private
   public :: linear_algebra_t, sscal, saxpy, sasum, sdot, isamax

   real(dk), parameter :: rZERO = 0.0_dk

   type, abstract :: linear_algebra_t
    ! Linear algebra functions
   contains
     procedure(SGBCO), deferred :: SGBCO
     procedure(SGBFA), deferred :: SGBFA
     procedure(SGBSL), deferred :: SGBSL
     procedure(SGECO), deferred :: SGECO
     procedure(SGEFA), deferred :: SGEFA
     procedure(SGESL), deferred :: SGESL
     procedure(tridiag), deferred :: tridiag
   end type linear_algebra_t

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBCO( this, ABD, N, ML, MU, IPVT, RCOND, Z )
    ! Factors a real band matrix by Gaussian elimination and estimates the
    ! condition of the matrix
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(inout) :: ABD(:,:) ! Banded matrix A(Ax = b) in a packed format
    integer,                 intent(in)    :: N        ! Order of the original matrix
    integer,                 intent(in)    :: ML       ! Number of diagonals below the main diagonal
    integer,                 intent(in)    :: MU       ! Number of diagonals above the main diagonal
    integer,                 intent(out)   :: IPVT(:)  ! Pivot vector
    real(dk),                intent(out)   :: RCOND    ! Estimate of the condition of the inverse of ABD
    real(dk),                intent(out)   :: Z(:)     ! Working vector
  end subroutine SGBCO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBFA( this, ABD, N, ML, MU, IPVT, INFO )
    ! Factors a real band matrix by elimination
    !
    ! SGBFA is usually called by SBGCO, but it can be called directly if
    ! RCOND is not needed.
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(inout) :: ABD(:,:) ! Banded matrix A(Ax = b) in a packed format
    integer,                 intent(in)    :: N        ! Order of the original matrix
    integer,                 intent(in)    :: ML       ! Number of diagonals below the main diagonal
    integer,                 intent(in)    :: MU       ! Number of diagonals above the main diagonal
    integer,                 intent(out)   :: IPVT(:)  ! Pivot vector
    integer,                 intent(out)   :: INFO     ! 1 = ABD is singular; 0 = otherwise
  end subroutine SGBFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBSL( this, ABD, N, ML, MU, IPVT, B, JOB )
    ! Solves the real band system
    ! A * X = B or transpose(A) * X = B
    ! using the factors computed by SGBCO or SGBFA
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(in)    :: ABD(:,:) ! Output of SGBCO or SGBFA
    integer,                 intent(in)    :: N        ! Order of the original matrix
    integer,                 intent(in)    :: ML       ! Number of diagonals below the main diagonal
    integer,                 intent(in)    :: MU       ! Number of diagonals above the main diagonal
    integer,                 intent(in)    :: IPVT(:)  ! Pivot vector
    real(dk),                intent(inout) :: B(:)     ! Right-hand side vector
    integer,                 intent(in)    :: JOB      ! 0 = solve A*X=B, otherwise = solve transpose(A)*X=B
  end subroutine SGBSL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGECO( this, A, N, IPVT, RCOND, Z )
    ! Factors a real matrix by Gaussian Elimination and estimates the
    ! condition of the matrix
    !
    ! If RCOND is not needed, SGEFA Is slightly faster
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(inout) :: A(:,:) ! INPUT: The matrix to be factored
    ! OUTPUT: An Upper Triangular matrix and the multipliers
    ! which were used to obtain it.
    ! The factorization can be written  A = L*U , where
    ! L  is a product of perumtation and unit lower
    ! triangular matrices and  U  is upper triangular.
    integer,                 intent(in)    :: N       ! The order of matrix A
    integer,                 intent(out)   :: IPVT(:) ! Pivot vector
    real(dk),                intent(out)   :: RCOND   ! An estimate of the
    ! reciprocal condition of  A .
    ! For the system  A*X = B , relative perturbations
    ! in  A  and  B  Oof size EPSILON  may cause
    ! relative perturbations in  X  of size  EPSILON/RCOND .
    ! If  RCOND  is so small that the logical expression
    ! 1.0 + RCOND .EQ. 1.0
    ! is true, then  A  may be singular to working
    ! precision.  In particular,  RCOND  is zero  if
    ! exact singularity is detected or the estimate
    ! underflows.
    real(dk),                intent(out)   :: Z(:)    ! Working vector
  end subroutine SGECO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGEFA( this, A, N, IPVT, INFO )
    ! Factors a real matrix by Gaussian Elimination
    !
    ! SGEFA is usually called by SGECO, but it can be called directly when
    ! RCOND is not needed.
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(inout) :: A(:,:) ! INPUT: The matrix to be factored
    ! OUTPUT: An Upper Triangular matrix and the multipliers
    ! which were used to obtain it.
    ! The factorization can be written  A = L*U , where
    ! L  is a product of perumtation and unit lower
    ! triangular matrices and  U  is upper triangular.
    integer,                 intent(in)    :: N       ! The order of matrix A
    integer,                 intent(out)   :: IPVT(:) ! Pivot vector
    integer,                 intent(out)   :: INFO    ! 0 = normal value
    ! K = If  U(K,K) .EQ. 0.0 .  This is not an error
    ! condition for this subroutine, but it does
    ! indicate that SGESL or SGEDI will divide by zero
    ! if called.  Use  RCOND  in SGECO for a reliable
    ! indication of singularity.
  end subroutine SGEFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGESL( this, A, N, IPVT, B, JOB )
    ! Solves the real system
    ! A * X = B or transpose(A) * X = B
    ! using the factors computed by SGECO or SGEFA
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(inout) :: this
    real(dk),                intent(in)    :: A(:,:)  ! Output of SGECO or SGEFA
    integer,                 intent(in)    :: N       ! Order of the A matrix
    integer,                 intent(in)    :: IPVT(:) ! Pivot vector
    real(dk),                intent(inout) :: B(:)    ! Right-hand side vector
    integer,                 intent(in)    :: JOB     ! 0 = solve A*X=B, otherwise = solve transpose(A)*X=B
  end subroutine SGESL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function tridiag( this, a, b, c, r ) result( u )
    ! Solves the tridiagonal system
    !
    ! The system to be solved is:
    !
    ! | b1 c1  0 ...                |   |  u1  |   |  r1  |
    ! | a2 b2 c2 ...                |   |  u2  |   |  r2  |
    ! |          ...                | * | ...  | = | ...  |
    ! |          ... aN-1 bN-1 cN-1 |   | uN-1 |   | rN-1 |
    ! |          ...  aN   bN   cN  |   |  uN  |   |  rN  |
    !
    use musica_constants,              only : dk => musica_dk
    import linear_algebra_t
    class(linear_algebra_t), intent(in) :: this
    real(dk),                intent(in) :: a(:) ! lower diagonal
    real(dk),                intent(in) :: b(:) ! primary diagonal
    real(dk),                intent(in) :: c(:) ! upper diagonal
    real(dk),                intent(in) :: r(:) ! right-hand side vector
    real(dk)                            :: u( size( b ) ) ! result vector
  end function tridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function isamax( N, SX, INCX )
    ! Returns the first I from 1 to N to maximize
    ! abs( SX( 1 +( I - 1 ) * INCX ) )

    integer,  intent(in) :: N     ! Number of elements in the vector of interest
    integer,  intent(in) :: INCX  ! Spacing of vector elements in SX
    real(dk), intent(in) :: SX(:) ! Vector of interest

    integer :: I,II
    real(dk) :: SMAX, XMAG

    if( N <= 0 ) then
      isamax = 0
    elseif( N == 1 ) then
      isamax = 1
    else
      SMAX = rZERO
      II = 1
      do I = 1, 1 + ( N - 1 ) * INCX, INCX
        XMAG = abs( SX( I ) )
        if( SMAX .lt. XMAG ) then
          SMAX = XMAG
          isamax = II
        endif
        II = II + 1
      enddo
    endif

   end function isamax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dk) function sasum( N, SX, INCX )
    ! Returns the sum from 1 to N-1 of abs( SX( 1 + I * INCX ) )

    integer,  intent(in) :: N     ! Number of elements in the vector of interest
    integer,  intent(in) :: INCX  ! Spacing of vector elements in SX
    real(dk), intent(in) :: SX(:) ! Vector of interest

    integer :: I, M

    sasum = rZERO
    if( N > 0 ) then
      if( INCX.ne. 1 ) then
        ! non-unit increments
        do I = 1, 1 + ( N - 1 ) * INCX, INCX
          sasum = sasum + abs( SX( I ) )
        enddo
      else
        ! unit increments
        M = mod( N, 6 )
        if( M .ne. 0 ) then
          ! clean-up loop so remaining vector
          ! length is a multiple of 6.
          do I = 1, M
            sasum = sasum + abs( SX( I ) )
          enddo
        endif
        ! unroll loop for speed
        do I = M + 1, N, 6
          sasum = sasum + abs( SX( I ) ) + abs( SX( I + 1 ) )                 &
                        + abs( SX( I + 2 ) ) + abs( SX( I + 3 ) )             &
                        + abs( SX( I + 4 ) ) + abs( SX( I + 5 ) )
        enddo
      endif
    endif

  end function sasum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dk) function SDOT( N, SX, INCX, SY, INCY )
    ! Dot product of vectors X and Y

    integer,  intent(in) :: N     ! Number of elements in vectors X and Y
    real(dk), intent(in) :: SX(:) ! Vector X
    integer,  intent(in) :: INCX  ! Spacing of elements in vector X
    real(dk), intent(in) :: SY(:) ! Vector Y
    integer,  intent(in) :: INCY  ! Spacing of elements in vector Y

    integer :: I, M, IX, IY

    SDOT = rZERO
    if( N > 0 ) then
      if ( INCX .eq. INCY .and. INCX .gt. 1 )  then
        do I = 1, 1 + ( N - 1 ) * INCX, INCX
           SDOT = SDOT + SX( I ) * SY( I )
        enddo
      else if ( INCX .eq. INCY .and. INCX .eq. 1 )  then
        ! equal, unit increments
        M = mod( N, 5 )
        if( M .ne. 0 ) then
          ! clean-up loop so remaining vector length
          ! is a multiple of 4.
          do I = 1, M
            SDOT = SDOT + SX( I ) * SY( I )
          enddo
        endif
        ! unroll loop for speed
        do I = M + 1, N, 5
          SDOT = SDOT + SX( I ) * SY( I ) + SX( I + 1 ) * SY( I + 1 )         &
                      + SX( I + 2 ) * SY( I + 2 ) + SX( I + 3 ) * SY( I + 3 ) &
                      + SX( I + 4 ) * SY( I + 4 )
        enddo
      else
        ! nonequal or nonpositive increments.
        IX = 1
        IY = 1
        if( INCX .lt. 0 ) IX = 1 + ( N - 1 ) * ( -INCX )
        if( INCY .lt. 0 ) IY = 1 + ( N - 1 ) * ( -INCY )
        do I = 1, N
          SDOT = SDOT + SX( IX ) * SY( IY )
          IX = IX + INCX
          IY = IY + INCY
        enddo
      endif
    endif

  end function SDOT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sscal( N, SA, SX, INCX )
    ! Calculates  X = A*X  (X = vector, A = scalar)

    integer,  intent(in)    :: N     ! Number of elements in the vector of interest
    real(dk), intent(in)    :: SA    ! Scale factor
    real(dk), intent(inout) :: SX(:) ! Vector of interest
    integer,  intent(in)    :: INCX  ! Spacing of the vector elements in SX

    integer :: I, M

    if( N > 0 ) then
      if( INCX.ne. 1 ) then
        do I = 1, 1 + ( N - 1 ) * INCX, INCX
          SX( I ) = SA * SX( I )
        enddo
      else
        M = mod( N, 5 )
        if( M .ne. 0 ) then
          ! clean-up loop so remaining vector length
          ! is a multiple of 5.
          do I = 1, M
            SX( I ) = SA * SX( I )
          enddo
        endif
        ! unroll loop for speed
        do I = M+1, N, 5
          SX( I )     = SA * SX( I )
          SX( I + 1 ) = SA * SX( I + 1 )
          SX( I + 2 ) = SA * SX( I + 2 )
          SX( I + 3 ) = SA * SX( I + 3 )
          SX( I + 4 ) = SA * SX( I + 4 )
        enddo
      endif
    endif

  end subroutine sscal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine saxpy( N, SA, SX, INCX, SY, INCY )
    ! Calculates Y = A*X + Y  (X, Y = vectors, A = scalar)

    integer,  intent(in)    :: N     ! Number of elements in the vector of interest
    real(dk), intent(in)    :: SA    ! Scalar multiplier
    real(dk), intent(in)    :: SX(:) ! Vector X
    integer,  intent(in)    :: INCX  ! Spacing of elements in X
    real(dk), intent(inout) :: SY(:) ! Vector Y
    integer,  intent(in)    :: INCY  ! Spacing of elements in Y

    integer :: I, M, IX, IY

    if( N > 0 .and. SA /= rZERO ) then
      if ( INCX .eq. INCY .and. INCX .gt. 1 ) then
        do I = 1, 1 + ( N - 1 ) * INCX, INCX
          SY( I ) = SY( I ) + SA * SX( I )
        enddo
      else if ( INCX .eq. INCY .and. INCX .eq. 1 )  then
        ! equal, unit increments
        M = mod( N, 4 )
        if( M .ne. 0 ) then
          ! clean-up loop so remaining vector length
          ! is a multiple of 4.
          do I = 1, M
            SY( I ) = SY( I ) + SA * SX( I )
          enddo
        endif
        ! unroll loop for speed
        do I = M + 1, N, 4
          SY( I )     = SY( I )     + SA * SX( I )
          SY( I + 1 ) = SY( I + 1 ) + SA * SX( I + 1 )
          SY( I + 2 ) = SY( I + 2 ) + SA * SX( I + 2 )
          SY( I + 3 ) = SY( I + 3 ) + SA * SX( I + 3 )
        enddo
      else
        ! nonequal or nonpositive increments.
        IX = 1
        IY = 1
        if( INCX .lt. 0 )  IX = 1 + ( N - 1 ) * ( -INCX )
        if( INCY .lt. 0 )  IY = 1 + ( N - 1 ) * ( -INCY )
        do I = 1, N
          SY( IY ) = SY( IY ) + SA * SX( IX )
          IX = IX + INCX
          IY = IY + INCY
        enddo
      endif
    endif

  end subroutine saxpy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_linear_algebra
