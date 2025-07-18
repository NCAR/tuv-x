      SUBROUTINE SETMTX( BDR, CBAND, CMU, CWT, DELM0, DTAUCP, GC, KK,
     &                   LAMBER, LYRCUT, MI, MI9M2, MXCMU, NCOL, NCUT,
     &                   NNLYRI, NN, NSTR, TAUCPR, WK )

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
c   Calls- ZEROIT
c +--------------------------------------------------------------------+


c     .. Scalar Arguments ..

      LOGICAL   LAMBER, LYRCUT
      INTEGER   MI, MI9M2, MXCMU, NCOL, NCUT, NN, NNLYRI, NSTR
      REAL      DELM0
c     ..
c     .. Array Arguments ..

      REAL      BDR( MI, 0:MI ), CBAND( MI9M2, NNLYRI ), CMU( MXCMU ),
     &          CWT( MXCMU ), DTAUCP( * ), GC( MXCMU, MXCMU, * ),
     &          KK( MXCMU, * ), TAUCPR( 0:* ), WK( MXCMU )
c     ..
c     .. Local Scalars ..

      INTEGER   IQ, IROW, JCOL, JQ, K, LC, LDA, NCD, NNCOL, NSHIFT
      REAL      EXPA, SUM
c     ..
c     .. External Subroutines ..

      EXTERNAL  ZEROIT
c     ..
c     .. Intrinsic Functions ..

      INTRINSIC EXP
c     ..


      CALL ZEROIT( CBAND, MI9M2*NNLYRI )

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

      END
