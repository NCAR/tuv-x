      SUBROUTINE inter1(ng,xg,yg, n,x,y)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Map input data given on single, discrete points, onto a discrete target  =*
*=  grid.                                                                    =*
*=  The original input data are given on single, discrete points of an       =*
*=  arbitrary grid and are being linearly interpolated onto a specified      =*
*=  discrete target grid.  A typical example would be the re-gridding of a   =*
*=  given data set for the vertical temperature profile to match the speci-  =*
*=  fied altitude grid.                                                      =*
*=  Some caution should be used near the end points of the grids.  If the    =*
*=  input data set does not span the range of the target grid, the remaining =*
*=  points will be set to zero, as extrapolation is not permitted.           =*
*=  If the input data does not encompass the target grid, use ADDPNT to      =*
*=  expand the input array.                                                  =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NG  - INTEGER, number of points in the target grid                    (I)=*
*=  XG  - REAL, target grid (e.g. altitude grid)                          (I)=*
*=  YG  - REAL, y-data re-gridded onto XG                                 (O)=*
*=  N   - INTEGER, number of points in the input data set                 (I)=*
*=  X   - REAL, grid on which input data are defined                      (I)=*
*=  Y   - REAL, input y-data                                              (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/95  Loop 10 restructured                                              =*
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
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

* output:
      REAL yg(ng)

* local:
      REAL slope
      INTEGER jsave, i, j
*_______________________________________________________________________

      jsave = 1
      DO 20, i = 1, ng
         yg(i) = 0.
         j = jsave
   10    CONTINUE
            IF ((x(j) .GT. xg(i)) .OR. (xg(i) .GE. x(j+1))) THEN
               j = j+1
               IF (j .LE. n-1) GOTO 10
*        ---- end of loop 10 ----
            ELSE
               slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
               yg(i) = y(j) + slope * (xg(i) - x(j))
               jsave = j
             ENDIF
   20 CONTINUE
*_______________________________________________________________________

      RETURN
      END
