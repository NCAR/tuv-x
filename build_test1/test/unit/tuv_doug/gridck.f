      SUBROUTINE gridck(k,n,x,ok)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Check a grid X for various improperties.  The values in X have to comply =*
*=  with the following rules:                                                =*
*=  1) Number of actual points cannot exceed declared length of X            =*
*=  2) Number of actual points has to be greater than or equal to 2          =*
*=  3) X-values must be non-negative                                         =*
*=  4) X-values must be unique                                               =*
*=  5) X-values must be in ascending order                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  K  - INTEGER, length of X as declared in the calling program          (I)=*
*=  N  - INTEGER, number of actual points in X                            (I)=*
*=  X  - REAL, vector (grid) to be checked                                (I)=*
*=  OK - LOGICAL, .TRUE. -> X agrees with rules 1)-5)                     (O)=*
*=                .FALSE.-> X violates at least one of 1)-5)                 =*
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
*= Copyright (C) 1994,95,96  University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

      INCLUDE 'params'

* input:
      INTEGER k, n
      REAL x(k)

* output:
      LOGICAL ok

* local:
      INTEGER i
*_______________________________________________________________________

      ok = .TRUE.

* check if dimension meaningful and within bounds

      IF (n .GT. k) THEN
         ok = .false.
         WRITE(kout,100)
         RETURN
      ENDIF         
  100 FORMAT('Number of data exceeds dimension')

      IF (n .LT. 2) THEN
         ok = .FALSE.
         WRITE(kout,101)
         RETURN
      ENDIF
  101 FORMAT('Too few data, number of data points must be >= 2')

* disallow negative grid values

      IF(x(1) .LT. 0.) THEN
         ok = .FALSE.
         WRITE(kout,105)
         RETURN
      ENDIF
  105 FORMAT('Grid cannot start below zero')

* check sorting

      DO 10, i = 2, n
         IF( x(i) .LE. x(i-1)) THEN
            ok = .FALSE.
            WRITE(kout,110)
            RETURN
         ENDIF
   10 CONTINUE
  110 FORMAT('Grid is not sorted or contains multiple values')
*_______________________________________________________________________

      RETURN
      END
