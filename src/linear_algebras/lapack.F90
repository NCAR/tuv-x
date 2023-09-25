! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_linear_algebra_lapack

  use musica_constants,                only : dk => musica_dk
  use tuvx_linear_algebra,             only : linear_algebra_t

  implicit none

  private
  public :: linear_algebra_lapack_t

  type, extends(linear_algebra_t) :: linear_algebra_lapack_t
    ! Linear algebra functions using LAPACK
  contains
    procedure :: SGBCO
    procedure :: SGBFA
    procedure :: SGBSL
    procedure :: SGECO
    procedure :: SGEFA
    procedure :: SGESL
    procedure :: tridiag
  end type linear_algebra_lapack_t

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBCO( this, ABD, N, ML, MU, IPVT, RCOND, Z )
    ! FACTORS A REAL BAND MATRIX BY GAUSSIAN ELIMINATION
    ! AND ESTIMATES THE CONDITION OF THE MATRIX.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
    real(dk), intent(inout) :: ABD(:,:) ! Banded matrix A(Ax = b) in a packed format
    integer,  intent(in)    :: N        ! Order of the original matrix
    integer,  intent(in)    :: ML       ! Number of diagonals below the main diagonal
    integer,  intent(in)    :: MU       ! Number of diagonals above the main diagonal
    integer,  intent(out)   :: IPVT(:)  ! Pivot vector
    real(dk), intent(out)   :: RCOND    ! Estimate of the condition of the inverse of ABD
    real(dk), intent(out)   :: Z(:)     ! Working vector

    call die_msg(590789017, "LAPACK function not implemented")

  end subroutine SGBCO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBFA( this, ABD, N, ML, MU, IPVT, INFO )
    ! FACTORS A REAL BAND MATRIX BY ELIMINATION.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
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

    call die_msg(252380743, "LAPACK function not implemented")

  end subroutine SGBFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGBSL( this, ABD, N, ML, MU, IPVT, B, JOB )
    ! SOLVES THE REAL BAND SYSTEM
    ! A * X = B  OR  TRANSPOSE(A) * X = B
    ! USING THE FACTORS COMPUTED BY SBGCO OR SGBFA.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
    real(dk), intent(in)    ::  ABD(:,:) ! Output from SGBCO or SGBFA
    integer,  intent(in)    ::  N        ! Order of the original matrix
    integer,  intent(in)    ::  ML       ! Number of diagonals below the main diagonal
    integer,  intent(in)    ::  MU       ! Number of diagonals above the main diagonal
    integer,  intent(in)    ::  IPVT(:)  ! Pivot vector from SGBCO or SGBFA
    real(dk), intent(inout) ::  B(:)     ! Right-hand side vector
    integer,  intent(in)    ::  JOB      ! = 0 TO SOLVE  A*X = B ,
    ! = NONZERO  TO SOLVE  TRANS(A)*X = B , WHERE TRANS(A)  IS THE TRANSPOSE.    

    call die_msg(759492682, "LAPACK function not implemented")

  end subroutine SGBSL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGESL( this, A, N, IPVT, B, JOB )
    ! SOLVES THE real SYSTEM
    ! A * X = B  OR  TRANS(A) * X = B
    ! USING THE FACTORS COMPUTED BY SGECO OR SGEFA.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
    integer,  intent(in)    :: N       ! Order of matrix A
    integer,  intent(in)    :: JOB     ! = 0  TO SOLVE  A*X = B ,
    ! = NONZERO  TO SOLVE TRANS(A)*X = B  WHERE TRANS(A)  IS THE TRANSPOSE.
    integer,  intent(in)    :: IPVT(:) ! Pivot vector from SGECO or SGERFA
    real(dk), intent(in)    :: A(:,:)  ! Output from SGECO or SGEFA
    real(dk), intent(inout) :: B(:)    ! Right-hand side vector

    call die_msg(366604622, "LAPACK function not implemented")

  end subroutine SGESL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGECO( this, A, N, IPVT, RCOND, Z )
    ! FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION
    ! AND ESTIMATES THE CONDITION OF THE MATRIX.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
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

    call die_msg(473658660, "LAPACK function not implemented")

  end subroutine SGECO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine SGEFA( this, A, N, IPVT, INFO )
    ! FACTORS A REAL MATRIX BY GAUSSIAN ELIMINATION.
    use musica_assert,    only : die_msg

    class(linear_algebra_lapack_t), intent(inout) :: this
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

    call die_msg(415820101, "LAPACK function not implemented")

  end subroutine SGEFA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function tridiag( this, a, b, c, r ) result( u )
    ! Solves a tridiagonal system.
    !
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

    use musica_assert,                 only : assert_msg
    external :: dgtsv

    class(linear_algebra_lapack_t), intent(in) :: this
    real(dk),            intent(in) :: a(:) ! lower diagonal
    real(dk),            intent(in) :: b(:) ! primary diagonal
    real(dk),            intent(in) :: c(:) ! upper diagonal
    real(dk),            intent(in) :: r(:) ! right-hand side vector
    real(dk)                        :: u( size( b ) ) ! result vector

    integer :: info

    u(:) = r(:)
    call dgtsv( size( b ), 1, a( 2 : size( a ) ), b, c( 1 : size( c ) - 1 ),  &
                u, size( b ), info )
    call assert_msg(236877362, info == 0, "Tridiagonal solver failure")

  end function tridiag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_linear_algebra_lapack    