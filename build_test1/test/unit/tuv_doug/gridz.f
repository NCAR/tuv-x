      SUBROUTINE gridz(z1,zout,izout,nz,z)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the altitude grid for all interpolations and radiative transfer   =*
*=  calculations.  Grid may be irregularly spaced.  All altitudes are in     =*
*=  kilometers (km).  The altitude at index 1 specifies the surface elevation=*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  nz  - INTEGER, number of altitude points (levels)                     (O)=*
*=  z   - REAL, vector of altitude levels (in km)                         (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  Original                                                                 =*
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
*= Copyright (C) 1994,95,96,97,98,99,2000,01,02  University Corporation for  =*
*= Atmospheric Research                                                      =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input surface elevation, altitude for output

      REAL z1, zout

* output: altitude working grid, index for output

      REAL z(kz)
      INTEGER nz
      INTEGER izout

* local:

      REAL zincr
      INTEGER i
      LOGICAL ok
*_______________________________________________________________________

* set vertical grid of the atmosphere.  All values should be in km.
* User specifies upright grid (surface at lowest km value, increasing
* upwards:
*     -  nz = total number of user levels
*     -  z(I) = altitude in km for each level.
* Note "levels" are vertical points
*      "layers" are vertical distances between levels

* set atmospheric level altitudes (in real km), including 
* top-most level.
* non-uniform spacing is possible 
* z(1) is the elevation of the surface (km asl), and can be specified either
* here or in the main progarm.
 
      z(1) = z1
      nz = 151
      zincr = 1.
      DO 10, i = 2, nz
         z(i) =  z(1) + FLOAT(i-1) * zincr
*         print*, z(i)
   10 CONTINUE
      
* Insert additional altitude for selected outputs.
* if within 1 meter of existing grid, do not insert (to avoid
*     numerical problems with very thin spherical shells)
* Set on 11/26/17?? (skip logic and goes to 24 continue
      zout = 0.0
      DO i = 1, nz
         IF(ABS(z(i) - zout) .LT. 0.001) THEN
            izout = i
            GO TO 24
         ENDIF
      ENDDO

* locate index for new altitude

      izout = 0
      DO i = 1, nz
         IF(z(i) .GT. zout) THEN
            izout = i
            GO TO 22
         ENDIF
      ENDDO
 22   CONTINUE
*      IF(izout .LE. 1) STOP 'zout not in range - '

* shift overlying levels and insert new point

      nz = nz + 1
      DO i = nz, izout + 1, -1
         z(i) = z(i-1)
      ENDDO
      z(izout) = zout

 24   CONTINUE

*     write to record:
      print*, nz, z(1), z(nz)
      WRITE(kout,*)'z-grid:',nz,z(1),z(nz)

* check grid for assorted improprieties:

      CALL gridck(kz,nz,z,ok)

      IF (.NOT. ok) THEN
         WRITE(kout,*)'STOP in GRIDZ:  The z-grid does not make sense'
!         STOP
      ENDIF
*_______________________________________________________________________

      RETURN
      END
