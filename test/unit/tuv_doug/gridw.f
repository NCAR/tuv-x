      SUBROUTINE gridw(nw,wl,wc,wu,pathname)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the altitude grid for all interpolations and radiative transfer   =*
*=  calculations.  Grid may be irregularly spaced.  Wavelengths are in nm.   =*
*=  No gaps are allowed within the wavelength grid.                          =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW  - INTEGER, number of wavelength grid _points_                     (O)=*
*=  WL  - REAL, vector carrying the lower limit of each wavel. interval   (O)=*
*=  WC  - REAL, vector carrying the center wavel of each wavel. interval  (O)=*
*=              (wc(i) = 0.5*(wl(i)+wu(i), i = 1..NW-1)                      =*
*=  WU  - REAL, vector carrying the upper limit of each wavel. interval   (O)=*
*=
*=  MOPT- INTEGER OPTION for wave-length IF 3 good for JO2                (O)=*
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

* output:

      REAL wl(kw), wc(kw), wu(kw)
      INTEGER nw

* local:

      INTEGER mopt
      REAL wincr
      INTEGER iw

      CHARACTER*80 fi
      CHARACTER*20 wlabel
      CHARACTER*80 pathname

      LOGICAL ok
*_______________________________________________________________________

**** chose wavelength grid

* some pre-set options
*     mopt = 1    equal spacing
*     mopt = 2    grid defined in data table
*     mopt = 3    user-defined


      mopt = 2
      IF (mopt .EQ. 1) GO TO 1
      IF (mopt .EQ. 2) GO TO 2
      IF (mopt .EQ. 3) GO TO 3

*_______________________________________________________________________

 1    CONTINUE
      wlabel = 'equal spacing'
      nw = 140 + 1
      wincr = 1.0
      DO 10, iw = 1, nw-1
         wl(iw) = 280. + wincr*FLOAT(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
 10   CONTINUE
      wl(nw) = wu(nw-1)
      GO TO 9

*_______________________________________________________________________

 2    CONTINUE

* Input from table.  In this example:
* Wavelength grid will be read from a file.
* First line of table is:  nw = number of wavelengths (no. of intervals + 1)
* Then, nw wavelengths are read in, and assigned to wl(iw)
* Finally, wu(iw) and wc(iw) are computed from wl(iw)
      wlabel = 'from file'
      OPEN(unit=kin,file=pathname,status='old')
      READ(kin,*) nw
      DO iw = 1, nw
         READ(kin,*) wl(iw)
      ENDDO
      CLOSE(kin)
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO
      GO TO 9

*_______________________________________________________________________

 3    CONTINUE

* user-defined grid.  In this example, a single calculation is used to 
* obtain results for two 1 nm wide intervals centered at 310 and 400 nm:
* interval 1 : 1 nm wide, centered at 310 nm
* interval 3 : 2 nm wide, centered at 400 nm
* (inteval 2 : 310.5 - 399.5 nm, required to connect intervals 1 & 3)

      nw = 4
      wl(1) = 309.5
      wl(2) = 310.5
      wl(3) = 399.5
      wl(4) = 400.5
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO

*_______________________________________________________________________

 9    CONTINUE

* write to record

      WRITE(kout,*) 'w-grid:',wlabel,nw,wl(1),wl(nw)

* check grid for assorted improprieties:

      CALL gridck(kw,nw,wl,ok)

      IF (.NOT. ok) THEN
         WRITE(kout,*)'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

      do iw = 1,nw-1
       write(88,*) iw, wl(iw),wl(iw+1),wc(iw)
      enddo
*_______________________________________________________________________

      RETURN
      END





