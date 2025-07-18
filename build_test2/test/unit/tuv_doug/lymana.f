      SUBROUTINE lymana(nz,o2col,secchi,dto2la,o2xsla)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Calculate the effective absorption cross section of O2 in the Lyman-Alpha=*
*=  bands and an effective O2 optical depth at all altitudes.  Parameterized =*
*=  after:  Chabrillat, S., and G. Kockarts, Simple parameterization of the  =*
*=  absorption of the solar Lyman-Alpha line, Geophysical Research Letters,  =*
*=  Vol.24, No.21, pp 2659-2662, 1997.                                       =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  O2COL   - REAL, slant overhead O2 column (molec/cc) at each specified (I)=*
*=            altitude                                                       =*
*=  DTO2LA  - REAL, optical depth due to O2 absorption at each specified  (O)=*
*=            vertical layer                                                 =*
*=  O2XSLA  - REAL, molecular absorption cross section in LA bands        (O)=*
*-----------------------------------------------------------------------------*
*=  EDIT HISTORY:                                                            =*
*=  01/98  Original                                                          =*
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
*= Copyright (C) 1994 - 1998 University Corporation for Atmospheric Research =*
*-----------------------------------------------------------------------------*

      IMPLICIT NONE

* input:

      INCLUDE 'params'
      INTEGER nz
      REAL o2col(kz)
      REAL secchi(kz)

* output

      REAL dto2la(kz,*), o2xsla(kz,*)

* local variables

      DOUBLE PRECISION rm(kz), ro2(kz)
      DOUBLE PRECISION b(3), c(3), d(3), e(3)
      DATA b/ 6.8431D-01, 2.29841D-01,  8.65412D-02/,
     >     c/8.22114D-21, 1.77556D-20,  8.22112D-21/,
     >     d/ 6.0073D-21, 4.28569D-21,  1.28059D-20/,
     >     e/8.21666D-21, 1.63296D-20,  4.85121D-17/

      INTEGER iz, i
      REAL xsmin
*------------------------------------------------------------------------------*
!     sm:  set minimum cross section
      xsmin = 1.D-20

      DO iz = 1, nz
        rm(iz) = 0.D+00
        ro2(iz) = 0.D+00                
        DO i = 1, 3
           rm(iz) = rm(iz) + b(i) * DEXP(-c(i) * DBLE(o2col(iz)))
        END DO
        ! TUV-x logic difference
        ! DO i = 1, 2
        DO i = 1, 3
           ro2(iz) = ro2(iz) + d(i) * DEXP(-e(i) * DBLE(o2col(iz)))
        ENDDO
       
      ENDDO
      
* calculate effective O2 optical depths and effective O2 cross sections
      DO iz = 1, nz-1

         IF (rm(iz) .GT. 1.0D-100) THEN
            IF (ro2(iz) .GT. 1.D-100) THEN
               o2xsla(iz,1) = ro2(iz)/rm(iz)
            ELSE
               ! TUV-x logic difference
               ! o2xsla(iz,1) = 0.               
               o2xsla(iz,1) = xsmin
            ENDIF

            IF (rm(iz+1) .GT. 0.) THEN

               dto2la(iz,1) = LOG(rm(iz+1)) / secchi(iz+1) 
     $                      - LOG(rm(iz))   / secchi(iz)

            ELSE
               dto2la(iz,1) = 1000.
            ENDIF
         ELSE
            dto2la(iz,1) = 1000.
            o2xsla(iz,1) = xsmin
         ENDIF

      ENDDO

* do top layer separately

      dto2la(nz,1) = 0.
      IF(rm(nz) .GT. 1.D-100) THEN
         o2xsla(nz,1) = ro2(nz)/rm(nz)
      ELSE
         o2xsla(nz,1) = xsmin
      ENDIF

*------------------------------------------------------------------------------*
      END
