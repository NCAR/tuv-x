      SUBROUTINE ZEROAL( EXPBEA, FLYR, OPRIM, TAUCPR, XR0, XR1,
     &                   CMU, CWT, PSI, WK, Z0, Z1, ZJ,
     &                   HLPR, YLM0,
     &                   ARRAY, CC, EVECC,
     &                   GL,
     &                   YLMC,
     &                   YLMU,
     &                   KK, LL, ZZ, ZPLK0, ZPLK1,
     &                   GC,
     &                   LAYRU, UTAUPR,
     &                   GU,
     &                   Z0U, Z1U, ZBEAM,
     &                   EVAL,
     &                   AMB, APB,
     &                   IPVT, Z,
     &                   RFLDIR, RFLDN, FLUP, UAVG, DFDT,
     &                   ALBMED, TRNMED,
     &                   U0U,
     &                   UU )

c         ZERO ARRAYS; NDn is dimension of all arrays following
c         it in the argument list

c   Called by- DISORT
c --------------------------------------------------------------------

c     .. Array Arguments ..

      INTEGER, intent(out) ::   IPVT( : ), LAYRU( : )
      REAL, intent(out)    ::
     $          ALBMED( : ), AMB( : ), APB( : ), ARRAY(:,:), CC(:,:),
     &          CMU( : ), CWT( : ), DFDT( : ), EVAL( : ), EVECC(:,:),
     &          EXPBEA( : ), FLUP( : ), FLYR( : ), GC( : ), GL( : ),
     &          GU( : ), HLPR( : ), KK( : ), LL( : ), OPRIM( : ),
     &          PSI( : ), RFLDIR( : ), RFLDN( : ), TAUCPR( : ),
     &          TRNMED( : ), U0U( : ), UAVG( : ), UTAUPR( : ), UU( : ),
     &          WK( : ), XR0( : ), XR1( : ), YLM0( : ), YLMC( : ),
     &          YLMU( : ), Z( : ), Z0( : ), Z0U( : ), Z1( : ), Z1U( : ),
     &          ZBEAM( : ), ZJ( : ), ZPLK0( : ), ZPLK1( : ), ZZ( : )

c     ..

      EXPBEA = 0.0
      FLYR   = 0.0
      OPRIM  = 0.0
      TAUCPR = 0.0
      XR0    = 0.0
      XR1    = 0.0

         CMU( : ) = 0.0
         CWT( : ) = 0.0
         PSI( : ) = 0.0
         WK( : )  = 0.0
         Z0( : )  = 0.0
         Z1( : )  = 0.0
         ZJ( : )  = 0.0

         HLPR( : ) = 0.0
         YLM0( : ) = 0.0

         ARRAY( :, : ) = 0.0
         CC( :, : )    = 0.0
         EVECC( :, : ) = 0.0

         GL( : ) = 0.0

         YLMC( : ) = 0.0

         YLMU( : ) = 0.0

         KK( : )    = 0.0
         LL( : )    = 0.0
         ZZ( : )    = 0.0
         ZPLK0( : ) = 0.0
         ZPLK1( : ) = 0.0

         GC( : ) = 0.0

         LAYRU( : )  = 0
         UTAUPR( : ) = 0.0

         GU( : ) = 0.0

         Z0U( : )   = 0.0
         Z1U( : )   = 0.0
         ZBEAM( : ) = 0.0

         EVAL( : ) = 0.0

         AMB( : ) = 0.0
         APB( : ) = 0.0

         IPVT( : ) = 0
         Z( : )    = 0.0

         RFLDIR( : ) = 0.
         RFLDN( : )  = 0.
         FLUP( : )   = 0.
         UAVG( : )   = 0.
         DFDT( : )   = 0.

         ALBMED( : ) = 0.
         TRNMED( : ) = 0.

         U0U( : ) = 0.

         UU( : ) = 0.

      END SUBROUTINE ZEROAL
