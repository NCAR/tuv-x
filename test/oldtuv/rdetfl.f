      MODULE ETFL

      use debug,      only : diagout

      IMPLICIT NONE

      private
      public :: rdetfl

      contains

      SUBROUTINE rdetfl(nw,wl,f)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read and re-grid extra-terrestrial flux data.                            =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*
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

      use tuv_params, only : kin, kout, hc, deltax

      integer, parameter :: kdata = 100000

* input: (wavelength grid)
      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(:)

* output: (extra terrestrial solar flux)
      REAL, intent(out) :: f(:)

* INTERNAL:

      INTEGER :: iw

* work arrays for input data files:

      INTEGER :: nhead, n, i, ierr
      REAL :: dum
      REAL :: x1(kdata)
      REAL :: y1(kdata)
      CHARACTER(len=40) :: fil

* data gridded onto wl(kw) grid:

      REAL :: yg1(size(wl))
      REAL :: yg2(size(wl))
      REAL :: yg3(size(wl))
      REAL :: yg4(size(wl))
      REAL :: wrk(size(wl)-1)

      INTEGER :: msun

      msun = 15

* simple files are read and interpolated here in-line. Reading of 
* more complex files may be done with longer code in a read#.f subroutine.


         WRITE(kout,*) 'odat/DATAE1/SUN/susim_hi.flx'
         CALL read1(nw,wl,yg1)

         nhead = 5
         fil = 'odat/DATAE1/SUN/atlas3_1994_317_a.dat'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         n = 5160
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), y1(i)
            y1(i) = y1(i) * 1.E-3
         ENDDO
         CLOSE (kin)
         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)

      write(*,*) ' '
      write(*,*) 'Diagnostics for atlas3_1994_317_a'
      write(*,*) 'read1: size model lambdaGrid = ',nw
      write(*,*) 'read1: lambdaGrid'
      write(*,'(1p10g15.7)') wl(:nw)
      write(*,*) ' '
      write(*,*) 'read1: size inputGrid = ',n
      write(*,*) 'read1: inputGrid'
      write(*,'(1p10g15.7)') x1(:n)
      write(*,*) ' '
      write(*,*) 'read1: size inputData = ',n
      write(*,*) 'read1: inputData'
      write(*,'(1p10g15.7)') y1(:n)

      call diagout( 'atlas.inputGrid.old', x1(:n) )
      call diagout( 'atlas.inputData.old', y1(:n) )

         CALL inter2(nw,wl,yg2,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         

      call diagout( 'atlas.interpolated.old', yg2(:nw-1) )

      write(*,*) ' '
      write(*,*) 'read1: size yg2 = ',size(yg2)
      write(*,*) 'read1: interpolated Etfl'
      write(*,'(1p10g15.7)') yg2(:nw-1)

         fil = 'odat/DATAE1/SUN/neckel.flx'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
         nhead = 11
         n = 496
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) dum, y1(i)
            if (dum < 630.0) then
              x1(i) = dum - 0.5
            elseif (dum >= 630.0 .and. dum < 870.0) then
              x1(i) = dum - 1.0
            else
              x1(i) = dum - 2.5
            endif
            y1(i) = y1(i) * 1.E4 * hc / (dum * 1.E-9)
         ENDDO
         CLOSE (kin)

         x1(n+1) = x1(n) + 2.5
         y1(n+1) = 0.0

      call diagout( 'neckel.inputGrid.old', x1(:n+1) )
      call diagout( 'neckel.inputData.old', y1(:n+1) )

      write(*,*) ' '
      write(*,*) 'Diagnostics for neckel.flx'
      write(*,*) 'read1: size model lambdaGrid = ',nw
      write(*,*) 'read1: lambdaGrid'
      write(*,'(1p10g15.7)') wl(:nw)
      write(*,*) ' '
      write(*,*) 'read1: size inputGrid = ',n+1
      write(*,*) 'read1: inputGrid'
      write(*,'(1p10g15.7)') x1(:n+1)
      write(*,*) ' '
      write(*,*) 'read1: size inputData = ',n+1
      write(*,*) 'read1: inputData'
      write(*,'(1p10g15.7)') y1(:n+1)
         call inter4(nw,wl,yg3,n+1,x1,y1,0)
      write(*,*) ' '
      write(*,*) 'read1: size yg3 = ',size(yg3)
      write(*,*) 'read1: interpolated Etfl'
      write(*,'(1p10g15.7)') yg3(:nw-1)
      call diagout( 'neckel.interpolated.old', yg3(:nw-1) )

         nhead = 8
         fil = 'odat/DATAE1/SUN/sao2010.solref.converted'
         write(kout,*) fil
         OPEN(UNIT=kin,FILE=fil,STATUS='old')
!        n = 80099 - nhead
         n = 80101 - nhead
         DO i = 1, nhead
            READ(kin,*)
         ENDDO
         DO i = 1, n
            READ(kin,*) x1(i), dum, y1(i), dum
         ENDDO
         CLOSE (kin)

         CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
         CALL addpnt(x1,y1,kdata,n,          0.,0.)
         CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
         CALL addpnt(x1,y1,kdata,n,      1.e+38,0.)
      call diagout( 'sao2010.inputGrid.old', x1(:n) )
      call diagout( 'sao2010.inputData.old', y1(:n) )
         CALL inter2(nw,wl,yg4,n,x1,y1,ierr)
         IF (ierr .NE. 0) THEN
            WRITE(*,*) ierr, fil
            STOP
         ENDIF         
      call diagout( 'sao2010.interpolated.old', yg4(:nw-1) )

*    for wl(iw) .lt. 150.01                                 susim_hi.flx
*    for wl(iw) .ge. 150.01 .and. wl(iw) .lt. 200.07        atlas3.flx
*    for wl(iw) .ge. 200.07 .and. wl(iw) .le. 1000.99       Chance and Kurucz 2010
*    for wl(iw) .gt. 1000.99                                Neckel & Labs 

         wrk = 0.0
         DO iw = 1, nw-1
            IF (wl(iw) < 150.01) THEN
               f(iw) = yg1(iw)
               wrk(iw) = yg1(iw)
            ELSE IF ((wl(iw) >= 150.01) .AND. wl(iw) < 200.07) THEN
               f(iw) = yg2(iw)
            ELSE IF ((wl(iw) >= 200.07) .AND. wl(iw) < 1000.99)THEN
               f(iw) = yg4(iw)
            ELSE IF (wl(iw) > 1000.99) THEN
               f(iw) = yg3(iw)
            ENDIF
         ENDDO

      call diagout( 'susim.etfl.old', wrk )
      call diagout( 'etfl.old', f(:nw-1) )

      END SUBROUTINE rdetfl

*=============================================================================*

      SUBROUTINE read1(nw,wl,f)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kin, deltax

* input: (wavelength grid)
      INTEGER, intent(in) :: nw
      REAL, intent(in) :: wl(:)

* output: (extra terrestrial solar flux)
      REAL, intent(inout) :: f(:)

* local:

      INTEGER :: ierr
      INTEGER :: i, j, n
      REAL :: lambda
      REAL :: lambda_hi(10000),irrad_hi(10000)
      CHARACTER(len=40) :: FIL

*_______________________________________________________________________

******* SUSIM irradiance 
*_______________________________________________________________________
* VanHoosier, M. E., J.-D. F. Bartoe, G. E. Brueckner, and
* D. K. Prinz, Absolute solar spectral irradiance 120 nm -
* 400 nm (Results from the Solar Ultraviolet Spectral Irradiance
* Monitor - SUSIM- Experiment on board Spacelab 2), 
* Astro. Lett. and Communications, 1988, vol. 27, pp. 163-168.
*     SUSIM SL2 high resolution (0.15nm) Solar Irridance data.
*     Irradiance values are given in milliwatts/m^2/nanomenters
*     and are listed at 0.05nm intervals.  The wavelength given is
*     the center wavelength of the 0.15nm triangular bandpass.
*     Normalized to 1 astronomical unit.
*  DATA for wavelengths > 350 nm are unreliable
* (Van Hoosier, personal communication, 1994).
*_______________________________________________________________________

** high resolution

      fil = 'odat/DATAE1/SUN/susim_hi.flx'
      OPEN(UNIT=kin,FILE=fil,STATUS='old')
      DO i = 1, 7
         READ(kin,*)
      ENDDO

      DO i = 1, 559
         READ(kin,*)lambda,(irrad_hi(10*(i-1)+j), j=1, 10)
      ENDDO

      CLOSE (kin)

* compute wavelengths, convert from mW to W

      n = 559*10
      DO i = 1, n
         lambda_hi(i) = 120.5 + REAL(i-1)*.05
         irrad_hi(i)  = 1.e-3 * irrad_hi(i)
      ENDDO
*_______________________________________________________________________

      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(1)*(1.-deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,                 0.,0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,
     >            lambda_hi(n)*(1.+deltax),0.)
      CALL addpnt(lambda_hi,irrad_hi,10000,n,              1.e38,0.)

      write(*,*) ' '
      write(*,*) 'Diagnostics for susim_hi.flx'
      write(*,*) 'read1: size model lambdaGrid = ',nw
      write(*,*) 'read1: lambdaGrid'
      write(*,'(1p10g15.7)') wl(:nw)
      write(*,*) ' '
      write(*,*) 'read1: size inputGrid = ',n
      write(*,*) 'read1: inputGrid'
      write(*,'(1p10g15.7)') lambda_hi(:n)
      write(*,*) ' '
      write(*,*) 'read1: size inputData = ',n
      write(*,*) 'read1: inputData'
      write(*,'(1p10g15.7)') irrad_hi(:n)

      call diagout( 'susim.inputGrid.old', lambda_hi(:n) )
      call diagout( 'susim.inputData.old', irrad_hi(:n) )

      CALL inter2(nw,wl,f,n,lambda_hi,irrad_hi,ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) ierr, fil
         STOP
      ENDIF

      call diagout( 'susim.interpolated.old', f(:nw-1) )

      write(*,*) ' '
      write(*,*) 'read1: size f = ',size(f)
      write(*,*) 'read1: interpolated Etfl'
      write(*,'(1p10g15.7)') f(:nw-1)

      END SUBROUTINE read1

*=============================================================================*

      SUBROUTINE read2(nw,wl,f)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Read extra-terrestrial flux data.  Re-grid data to match specified       =*
*=  working wavelength grid.                                                 =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  F      - REAL, spectral irradiance at the top of the atmosphere at    (O)=*
*=           each specified wavelength                                       =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kin

* input: (wavelength grid)
      INTEGER, intent(in) :: nw
      REAL, intent(in)    :: wl(:)
* output: (extra terrestrial solar flux)
      REAL, intent(out)   :: f(:)

*
      INTEGER :: iw
      REAL    :: yg(size(wl))

* local:

      INTEGER :: i, n
      INTEGER :: IDUM
      REAL :: DUM
      REAL :: x1(1000), y1(1000) 
      REAL :: x2(1000)
      REAL :: x3(1000)

*_______________________________________________________________________

*********WMO 85 irradiance

      OPEN(UNIT=kin,FILE=
     $ 'odat/DATAE1/SUN/wmo85.flx',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO

      n = 158
      DO i = 1, n
         READ(kin,*) idum, x1(i),x2(i),y1(i), dum, dum, dum
         x3(i) = 0.5 * (x1(i) + x2(i))
      ENDDO

      CLOSE (kin)

      x1(n+1) = x2(n)

C inter2: INPUT : average value in each bin 
C         OUTPUT: average value in each bin
C inter3: INPUT : total area in each bin
C         OUTPUT: total area in each bin

      CALL inter3(nw,wl,yg, n+1,x1,y1,0)

      DO iw = 1, nw-1
* from quanta s-1 cm-2 bin-1 to  watts m-2 nm-1
* 1.e4 * ([hc =] 6.62E-34 * 2.998E8)/(wc*1e-9) 
         
C the scaling by bin width needs to be done only if
C inter3 is used for interpolation
         yg(iw) = yg(iw) / (wl(iw+1)-wl(iw))
         f(iw) = yg(iw) * 1.e4 * (6.62E-34 * 2.998E8) / 
     $        ( 0.5 * (wl(iw+1)+wl(iw)) * 1.e-9)

      ENDDO

      END SUBROUTINE read2

      END MODULE ETFL
