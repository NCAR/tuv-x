      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
*=  ascending order                                                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
*=  Y    - REAL vector of length LD, y-values                            (IO)=*
*=  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
*=         program                                                           =*
*=  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
*=         N < LD.  On exit, N is incremented by 1.                          =*
*=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
*=  YNEW - REAL, y-value of point to be added                             (I)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  02/95  Original                                                          =*
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

* calling parameters

      INTEGER ld, n
      REAL x(ld), y(ld)
      REAL xnew, ynew
      INTEGER ierr

* local variables

      INTEGER insert
      INTEGER i

*-----------------------------------------------------------------------

* initialize error flag

      ierr = 0

* check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         STOP
      ENDIF

      insert = 1
      i = 2

* check, whether x is already sorted.
* also, use this loop to find the point at which xnew needs to be inserted
* into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//
     >                'in ascending order!'
           STOP
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

* if <xnew,ynew> needs to be appended at the end, just do so,
* otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN
 
         x(n+1) = xnew
         y(n+1) = ynew
  
      ELSE

* shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

* insert new point

         x(insert) = xnew
         y(insert) = ynew
  
      ENDIF

* increase total number of elements in x, y

      n = n+1

      END
