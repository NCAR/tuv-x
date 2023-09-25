      MODULE SPHERICAL_GEOM

      IMPLICIT NONE

      private
      public :: sphers, airmas

      contains

*=============================================================================*

      SUBROUTINE sphers(z, zen, dsdh, nid)

      use debug, only : diagout

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate slant path over vertical depth ds/dh in spherical geometry.    =*
*=  Calculation is based on:  A.Dahlback, and K.Stamnes, A new spheric model =*
*=  for computing the radiation field available for photolysis and heating   =*
*=  at twilight, Planet.Space Sci., v39, n5, pp. 671-683, 1991 (Appendix B)  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  ZEN     - REAL, solar zenith angle (degrees)                          (I)=*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  double precision fix for shallow layers - Julia Lee-Taylor Dec 2000      =*
*-----------------------------------------------------------------------------*
*= This program is free software;  you can redistribute it and/or modify     =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= To contact the authors, please mail to:                                   =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu                                            =*
*-----------------------------------------------------------------------------*
 
      use tuv_params, only : radius, pi

* input
      REAL, intent(in) :: zen
      REAL, intent(in) ::  z(:)

* output
      INTEGER, intent(out) :: nid(0:)
      REAL, intent(out)    :: dsdh(0:,:)

* more program constants
      REAL :: re
      REAL :: ze(size(z))

      REAL, parameter ::  dr = pi/180._8

* local 

      real(8), parameter :: rZERO = 0._8
      real(8), parameter :: rONE  = 1._8
      real(8), parameter :: rNINETY  = 90._8

      INTEGER :: nz
      INTEGER :: i, j, k
      INTEGER :: id
      INTEGER :: nlayer
      REAL(8) :: zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
      REAL    :: zd(0:size(z)-1)

*-----------------------------------------------------------------------------

      zenrad = zen*dr

      nz = size(z)
* number of layers:
      nlayer = nz - 1

* include the elevation above sea level to the radius of the earth:
      re = radius + z(1)
* correspondingly z changed to the elevation above earth surface:
      ze = z - z(1)

* inverse coordinate of z
      zd(0) = ze(nz)
      DO k = 1, nlayer
        zd(k) = ze(nz - k)
      END DO

* initialize dsdh(i,j), nid(i)
      nid  = 0
      dsdh = rZERO

* calculate ds/dh of every layer
      layer_loop: DO i = 0, nlayer
        rpsinz = (re + zd(i)) * SIN(zenrad)
        IF ( (zen > rNINETY) .AND. (rpsinz < re) ) THEN
           nid(i) = -1
        ELSE
*
* Find index of layer in which the screening height lies
*
           IF( zen <= rNINETY ) THEN
              id = i
           ELSE
              DO j = 1, nlayer
                 IF( (rpsinz < ( zd(j-1) + re ) ) .AND.
     $               (rpsinz >= ( zd(j) + re )) ) id = j
              ENDDO
           END IF
 
           DO j = 1, id
             sm = rONE
             IF(j == id .AND. id == i .AND. zen > rNINETY)
     $          sm = -rONE
             rj = re + zd(j-1)
             rjp1 = re + zd(j)
             dhj = zd(j-1) - zd(j)
             ga = rj*rj - rpsinz*rpsinz
             gb = rjp1*rjp1 - rpsinz*rpsinz
             ga = MAX( rZERO,ga )
             gb = MAX( rZERO,gb )
 
             IF(id > i .AND. j == id) THEN
                dsj = SQRT( ga )
             ELSE
                dsj = SQRT( ga ) - sm*SQRT( gb )
             END IF
             dsdh(i,j) = dsj / dhj
           ENDDO
 
           nid(i) = id
        END IF
      ENDDO layer_loop

*-----------------------------------------------------------------------------

      END SUBROUTINE sphers

*=============================================================================*

      SUBROUTINE airmas(dsdh, nid, aircol, vcol, scol)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate vertical and slant air columns, in spherical geometry, as a    =*
*=  function of altitude.                                                    =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  DSDH    - REAL, slant path of direct beam through each layer crossed  (O)=*
*=            when travelling from the top of the atmosphere to layer i;     =*
*=            DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                            =*
*=  NID     - INTEGER, number of layers crossed by the direct beam when   (O)=*
*=            travelling from the top of the atmosphere to layer i;          =*
*=            NID(i), i = 0..NZ-1                                            =*
*=  VCOL    - REAL, output, vertical air column, molec cm-2, above level iz  =*
*=  SCOL    - REAL, output, slant air column in direction of sun, above iz   =*
*=            also in molec cm-2                                             =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : largest

* Input:
      INTEGER, intent(in) :: nid(0:)
      REAL, intent(in)    :: dsdh(0:,:)
      REAL, intent(in)    :: aircol(:)

* output: 

      REAL, intent(out)   :: vcol(:), scol(:)

* internal:

      INTEGER :: nz
      INTEGER :: id, j
      REAL    :: accum

* calculate vertical and slant column from each level:
* work downward

      nz = size(aircol)
      accum = aircol(nz)
      DO id = 1, nz - 1
         accum = accum + aircol(nz-id)
         vcol(nz-id) = accum
      ENDDO

      scol(nz) = dsdh(1,1)*aircol(nz)
      DO id = 1, nz - 1
         accum = scol(nz)
         IF(nid(id) < 0) THEN
            accum = largest
         ELSE
* single pass layers:
            DO j = 1, MIN(nid(id), id)
               accum = accum + aircol(nz-j)*dsdh(id,j)
            ENDDO
* double pass layers:
            DO j = MIN(nid(id),id)+1, nid(id)
               accum = accum + 2.*aircol(nz-j)*dsdh(id,j)
            ENDDO
         ENDIF
         scol(nz - id) = accum
      ENDDO
      
      END SUBROUTINE airmas

      END MODULE SPHERICAL_GEOM
