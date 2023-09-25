      MODULE GRIDS

      IMPLICIT NONE

      private
      public :: gridw, gridz, gridt

      contains

      SUBROUTINE gridw(wstart, wstop, nwint, wl,wc,wu)

      use tuv_params, only : kw, kin
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Create the wavelength grid for all interpolations and radiative transfer =*
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

* input:

      INTEGER, intent(inout) :: nwint
      REAL, intent(in)       :: wstart, wstop

* output:

      REAL, allocatable, intent(out) :: wl(:), wc(:), wu(:)

* local:

      integer, allocatable :: wn(:)
      INTEGER :: mopt, nw
      INTEGER :: iw, i
      REAL    :: wincr

      CHARACTER(len=40) :: fi
      CHARACTER(len=20) :: wlabel

      INTEGER :: mrefr
      REAL    :: airout
      REAL    :: dum
      LOGICAL :: ok
*_______________________________________________________________________

**** chose wavelength grid

* some pre-set options
*     mopt = 1    equal spacing
*     mopt = 2    grid defined in data table
*     mopt = 3    user-defined
*     mopt = 4    fast-TUV, troposheric wavelengths only
*     mopt = 5    high resolution grid for O3 isotopologue study
*     mopt = 6    create uniform grid in air-wavelength scale

*     mopt = 10  Landgraf and Crutzen, 1998
*     mopt = 11  fastJ, Wild et al. 2000
*     mopt = 12  fastJ2, Bian and Prather, 2002
*     mopt = 13  UV-b, UV-a, Visible

      select case( nwint )
        case( -7 )
          mopt = 4
        case( -10 )
          mopt = 10
        case( -11 )
          mopt = 11
        case( -12 )
          mopt = 12
        case( -13 )
          mopt = 13
        case( -156 )
          mopt = 2
        case default
          mopt = 1
      end select

      SELECT CASE( mopt )

*_______________________________________________________________________

        CASE( 1 )

      wlabel = 'equal spacing'
      nw = nwint + 1
      allocate( wl(nw), wc(nw-1), wu(nw-1) )

      wincr = (wstop - wstart) / REAL (nwint)
      DO iw = 1, nw-1
         wl(iw) = wstart + wincr*real(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = .5*(wl(iw) + wu(iw))
      ENDDO
      wl(nw) = wu(nw-1)

*_______________________________________________________________________

        CASE( 2 )

* Input from table.  In this example:
* Wavelength grid will be read from a file.
* First line of table is:  nw = number of wavelengths (no. of intervals + 1)
* Then, nw wavelengths are read in, and assigned to wl(iw)
* Finally, wu(iw) and wc(iw) are computed from wl(iw)

c      wlabel = 'isaksen.grid'
      wlabel = 'combined.grid'

      fi = 'odat/DATAE1/GRIDS/'//wlabel
      OPEN(unit=kin,file=fi,status='old')
      READ(kin,*) nw
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      DO iw = 1, nw
         READ(kin,*) wl(iw)
      ENDDO
      CLOSE(kin)
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO

*_______________________________________________________________________

        CASE( 3 )

* user-defined grid.  In this example, a single calculation is used to 
* obtain results for two 1 nm wide intervals centered at 310 and 400 nm:
* interval 1 : 1 nm wide, centered at 310 nm
* interval 3 : 2 nm wide, centered at 400 nm
* (inteval 2 : 310.5 - 399.5 nm, required to connect intervals 1 & 3)

      nw = 4
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      wl(1) = 309.5
      wl(2) = 310.5
      wl(3) = 399.5
      wl(4) = 400.5
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO

*_______________________________________________________________________

        CASE( 4 )
      wlabel = 'fast-TUV tropospheric grid'
      
      fi = 'odat/DATAE1/GRIDS/fast_tuv.grid'
      OPEN(UNIT=kin,FILE=fi,STATUS='old')
      DO iw = 1, 4
         READ(kin,*)
      ENDDO

* skip wavelength shorter than 289.9 nm

      DO iw = 1, 10
         READ(kin,*)
      ENDDO
      nw = 8
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      DO iw = 1, nw-1
         READ(kin,*) dum, wl(iw), dum, dum
      ENDDO
      wl(nw) = dum
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO

*_______________________________________________________________________

        CASE( 5 )

* use standard grid up to 205.8 nm
* elsewhere, use 10 cm-1 grid to 1000 nm

      nw = 3859 + 38
      allocate( wl(nw), wc(nw-1), wu(nw-1), wn(nw) )

      wlabel = 'combined.grid'
      fi = 'odat/DATAE1/GRIDS/'//wlabel
      OPEN(unit=kin,file=fi,status='old')
      READ(kin,*) 
      DO iw = 1, 38
         READ(kin,*) wl(iw)
      ENDDO
      CLOSE(kin)


      DO i = 1, 3859
         iw = 3859 - i + 39
         wn(iw) = 10000 + 10*(i-1)
         wl(iw) = 1.E7/real(wn(iw))
      ENDDO

      nwint = nw - 1
      DO iw = 1, nwint
         wu(iw) = wl(iw+1)
         wc(iw) = (wl(iw) + wu(iw))/2.
      ENDDO

*_______________________________________________________________________

        CASE( 6 )

***** Correction for air-vacuum wavelength shift:
* The TUV code assumes that all working wavelengths are strictly IN-VACUUM. This is for ALL
* spectral data including extraterrestrial fluxes, ozone (and other) absorption cross sections,
* and various weighting functons (action spectra, photolysis cross sections, instrument spectral
* response functions).  If the original data are specified in-air, conversion to in-vacuum must be
* made when reading those data.

*  Occasionally, users may want their results to be given for wavelengths measured IN-AIR.
*   The shift between IN-VACUUM and IN-AIR wavelengths depends on the index of refraction
*   of air, which in turn depends on the local density of air, which in turn depends on
*   altitude, temperature, etc.
*  Here, we provide users with the option to use a wavelength grid IN-AIR, at the air density
*   corresponding to the selected altitude, airout.
*   The actual radiative transfer calculations will be done strictly with IN-VACUUM values.  

* create grid that will be nicely spaced in air wavelengths.

      wlabel = 'grid in air wavelengths'
      nw = nwint + 1
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      wincr = (wstop - wstart) / real (nwint)
      DO iw = 1, nw-1
         wl(iw) = wstart + wincr*real(iw-1)
         wu(iw) = wl(iw) + wincr
         wc(iw) = ( wl(iw) + wu(iw) )/2.
      ENDDO
      wl(nw) = wu(nw-1)

* shift by refractive index to vacuum wavelengths, for use in tuv

      airout = 2.45e19
      mrefr = 1
      CALL wshift(mrefr, nw,    wl, airout)
      CALL wshift(mrefr, nwint, wc, airout)
      CALL wshift(mrefr, nwint, wu, airout)

c      do iw = 1, nw-1
c         write(33,333) wl(iw), wc(iw), wu(iw)
c      enddo
c 333  format(3(0pf11.4))

*_______________________________________________________________________
* Landgraf and Crutzen 1998
        CASE( 10 )
      nw = 6
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      wl(1) = 289.0
      wl(2) = 305.5
      wl(3) = 313.5
      wl(4) = 337.5
      wl(5) = 422.5
      wl(6) = 752.5
      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO
      
*_______________________________________________________________________
* Wild 2000
        CASE( 11 )
      nw = 8
      allocate( wl(nw), wc(nw-1), wu(nw-1) )
      wl(1) = 289.00
      wl(2) = 298.25
      wl(3) = 307.45
      wl(4) = 312.45
      wl(5) = 320.30
      wl(6) = 345.0
      wl(7) = 412.5
      wl(8) = 850.0

      DO iw = 1, nw-1
         wu(iw) = wl(iw+1)
         wc(iw) = 0.5*(wl(iw) + wu(iw))
      ENDDO

*_______________________________________________________________________
* Bian and Prather 2002
        CASE( 12 )
          nw = 8
          allocate( wl(nw), wc(nw-1), wu(nw-1) )
          wl = (/ 291., 298.3, 307.5, 312.5, 320.3, 345., 412.5, 850. /)

          DO iw = 1, nw-1
            wu(iw) = wl(iw+1)
            wc(iw) = 0.5*(wl(iw) + wu(iw))
          ENDDO

*_______________________________________________________________________
* UV-b, UV-A, Vis

        CASE( 13 )
          nw = 4
          allocate( wl(nw), wc(nw-1), wu(nw-1) )
          wl = (/ 280., 315., 400., 700. /)

          DO iw = 1, nw-1
            wu(iw) = wl(iw+1)
            wc(iw) = 0.5*(wl(iw) + wu(iw))
          ENDDO

      CASE DEFAULT
        write(*,*) ' '
        write(*,*) 'gridw: option ',mopt,' is invalid'
        stop 'Bad Option'
      END SELECT
*_______________________________________________________________________

* check grid for assorted improprieties:

      CALL gridck(kw,wl,ok)

      IF (.NOT. ok) THEN
         WRITE(*,*)'STOP in GRIDW:  The w-grid does not make sense'
         STOP
      ENDIF

      END SUBROUTINE gridw

*=============================================================================*

      SUBROUTINE gridz(zstart, zstop, nz, z, zout, izout)
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

      use tuv_params, only : kz

* input surface elevation, altitude for output

      REAL, intent(in)    :: zstart, zstop
      REAL, intent(inout) :: zout

* output: altitude working grid, index for output

      INTEGER, intent(inout)  :: nz
      INTEGER, intent(out)    :: izout
      REAL, allocatable, intent(out) :: z(:)

* local:

      INTEGER :: i, n, nlev
      INTEGER :: grid_opt
      LOGICAL :: ok
      REAL    :: zincr
*_______________________________________________________________________

* Set vertical grid of the atmosphere: atmospheric level altitudes 
* (in real km), including top-most level.
* User specifies grid (surface at lowest km value), increasing
* upwards:
*     -  nz = total number of user levels
*     -  z(I) = altitude in km for each level.
* z(1) is the elevation of the surface (km asl), and can be specified either
* here or in the main program.
* Non-uniform spacing is possible:
*     - zincr = altitude increment between current and previous level (km)
*     - nlev = number of levels in current equally-spaced section
*     - n = index of top level of equally-spaced section
* Note "levels" are vertical points
*      "layers" are vertical distances between levels


* Grid selection options:
* 1 = standard equally spaced grid, manual
* 2 = standard equally spaced grid, auto generated
* 3 = variable spacing grid, example for snow
* 4 = mirage z-grid for Mexico City
* 5 = arbitrary user-defined grid

      grid_opt = 1
      SELECT CASE( grid_opt )
        CASE( 1 )
*-----grid option 1: manual -----------------
* entire grid (nz levels) in increments zincr 
      WRITE(*,*) 'equally spaced z-grid'
      zincr = (zstop - zstart) / REAL(nz - 1)
      allocate( z(nz) )
      z(1) = zstart
      DO i = 2, nz
         z(i) = z(1) + zincr*REAL(i-1)
      ENDDO

        CASE( 2 )
*-----grid option 2: automatic -----------------
* entire grid (nz levels) in increments zincr 

      WRITE(*,*) 'equally spaced z-grid'
      zincr = (zstop - zstart) / real(nz - 1)
      nlev = nz-1
      n = 1
      CALL buildz(zincr, nlev, n, z)

        CASE( 3 )
*-----grid option 3: variable grid example----------------------
*-----copy & edit this section for non-uniform grid----
* the example provided below is high vertical resolution in 
*   snow, with atmosphere above it.
      WRITE(*,*) 'snow-atmosphere grid'
* 0.-10. cm from ground, in 1 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-5
      nlev = 10
      n = 1
      CALL buildz(zincr,nlev,n,z)

* 10-90 cm from ground, in 10 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-4
      nlev = 8
      CALL buildz(zincr,nlev,n,z)

* 90-95 cm from ground, in 1x 5 cm increment ( 1 cm = 1e-5 km):
      zincr = 5.e-5
      nlev = 1
      CALL buildz(zincr,nlev,n,z)

* 95-99 cm from ground, in 4x 1 cm increments ( 1 cm = 1e-5 km):
      zincr = 1.e-5
      nlev = 4
      CALL buildz(zincr,nlev,n,z)

* 99-99.5 cm from ground, in 1x 0.5 cm increment ( 1 cm = 1e-5 km):
      zincr = 5.e-6
      nlev = 1
      CALL buildz(zincr,nlev,n,z)

* 99.5 centimeters - 1m, in 0.1 cm increments (1 cm = 1e-5 km):
      zincr = 1.e-6
      nlev = 5
      CALL buildz(zincr,nlev,n,z)

*atmosphere
* 1.-10. m in 1 m increments
      zincr = 1.e-3
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 10.-100 m in 10 m increments
      zincr = 1.e-2
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 100.- 1000. meters, in 100 m increments
      zincr = 1.e-1
      nlev = 9
      CALL buildz(zincr,nlev,n,z)

* 1.-2. km in 1 km increments
      zincr = 1.
      nlev =  1
      CALL buildz(zincr,nlev,n,z)

* 2.-80. km in 2 km increments
      zincr = 2.
      nlev =  39
      CALL buildz(zincr,nlev,n,z)

        CASE( 4 )
*-----grid option 4:  grid for Mexico City

      WRITE(*,*) 'mirage z-grid'

* grid for mirage km: incr(range)i 
* 0.1(0-4)   2-41
* 0.2(4-8)   42-61
*   1(8-30)  62-83
*   2(30-50) 84-93 
*   5(50-80) 94-99

      nz = 99
      z(1) = zstart
      DO i = 2, 41
         z(i) = z(1) + 0.1*real(i-1)
      ENDDO
      DO i = 42, 61
         z(i) = z(41) + 0.2*real(i-41)
      ENDDO
      DO i = 62, 83
         z(i) = z(61) + 1.*real(i-61)
      ENDDO
      DO i = 84, 93
         z(i) = z(83) + 2.*real(i-83)
      ENDDO
      DO i = 94, 99
         z(i) = z(93) + 5.*real(i-93)
      ENDDO
		
        CASE( 5 )
*-----grid option 5: user defined

* insert your grid values here:
* specify:
*  nz = total number of altitudes
* Table:  z(iz), where iz goes from 1 to nz
* trivial example of 2-layer (3-altitudes) shown below, user should modify
      WRITE(*,*) 'user-defined grid, named...'
      nz = 3
      z(1) = 0.
      z(2) = 10.
      z(3) = 80.

        CASE( 6 )
*-----end of user options.
*-----grid option 6: high resolution window

* insert your grid values here:
* specify:
*  nz = total number of altitudes
* Table:  z(iz), where iz goes from 1 to nz

      WRITE(*,*) 'user-defined grid, named...'

      END SELECT

* Insert additional altitude for selected outputs.

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
            EXIT
         ENDIF
      ENDDO
      IF(izout .LE. 1) STOP 'zout not in range - '

* shift overlying levels and insert new point

      z        = [z(1:izout-1),zout,z(izout:)]

 24   CONTINUE

* check grid for assorted improprieties:

      CALL gridck(kz,z,ok)

      IF (.NOT. ok) THEN
         WRITE(*,*)'STOP in GRIDZ:  The z-grid does not make sense'
         STOP
      ENDIF

      END SUBROUTINE gridz

*=============================================================================*
      SUBROUTINE buildz(zincr,nlev,n,z)
*-----------------------------------------------------------------------------*
*= Purpose: to construct the altitude grid from parameters in gridz          =*
*-----------------------------------------------------------------------------*


      INTEGER, intent(in)    :: nlev
      INTEGER, intent(inout) :: n
      REAL, intent(in)       :: zincr
      REAL, intent(inout)    :: z(:)

      INTEGER :: i, j

      j = 0
      DO i = n + 1, n + nlev
        j = j + 1
        z(i) = z(n) + real(j)*zincr
      ENDDO
      n = n + nlev

      END SUBROUTINE buildz

*=============================================================================*

      SUBROUTINE gridt(lat, lon, tmzone,
     $     iyear, imonth, iday,
     $     lzenit, tstart, tstop,
     $     t, sza, sznoon, esrm2)

      use orbit, only        : calend, sunae

*-----------------------------------------------------------------------------*
*=  Subroutine to create time (or solar zenith angle) grid                   =*
*=  Also computes earth-sun distance (1/R**2) correction.                    =*
*-----------------------------------------------------------------------------*

* INPUTS
      INTEGER, intent(in) :: iyear, imonth, iday
      REAL, intent(in)    :: lat, lon, tmzone
      REAL, intent(in)    :: tstart, tstop
      LOGICAL, intent(in) :: lzenit

* OUTPUTS
      REAL, intent(out) :: sznoon
      REAL, intent(out) :: sza(:)
      REAL, intent(out) :: t(:), esrm2(:)

* INTERNAL
      INTEGER :: it
      INTEGER :: nt
      INTEGER :: jday, nday
      REAL    :: ut, dt
      REAL    :: az, el, elnoon, soldia, soldst
      LOGICAL :: oky, okm, okd

*  switch for refraction correction to solar zenith angle. Because
* this is only for the observed sza at the surface, do not use.

      LOGICAL, parameter :: lrefr = .FALSE.

***************
      nt = size(sza)
      IF(nt == 1) THEN
         dt = 0.
      ELSE
         dt = (tstop - tstart) / real(nt - 1)
      ENDIF

      DO it = 1, nt
         t(it) = tstart + dt * real(it - 1)
* solar zenith angle calculation:
*  If lzenit = .TRUE., use selected solar zenith angles, also
*  set Earth-Sun distance to 1 AU.
         IF (lzenit) THEN
            sza(it) = t(it)
            esrm2(it) = 1.
*  If lzenit = .FALSE., compute solar zenith angle for specified
* location, date, time of day.  Assume no refraction (lrefr = .FALSE.)
*  Also calculate corresponding
* Earth-Sun correcton factor. 
         ELSE
            CALL calend(iyear, imonth, iday,
     $           jday, nday, oky, okm, okd)
            IF( oky .AND. okm .AND. okd) THEN
               ut = t(it) - tmzone
               CALL sunae(iyear, jday, ut, lat, lon, lrefr,
     $              elnoon, az, el, soldia, soldst )
               sza(it) = 90. - el
               sznoon = 90. - elnoon
               esrm2(it) = 1./(soldst*soldst)
            ELSE
               WRITE(*,*) '**** incorrect date specification'
               STOP ' in gridt '
            ENDIF
         ENDIF
      ENDDO

      END SUBROUTINE gridt

*=============================================================================*

      SUBROUTINE gridck(k,x,ok)
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

      use tuv_params, only : kout

* input:
      INTEGER, intent(in) :: k
      REAL, intent(in)    :: x(:)

* output:
      LOGICAL, intent(out) :: ok

      INTEGER :: n

      ok = .false.
      n  = size(x)

* check if dimension meaningful and within bounds

      IF (n > k) THEN
         WRITE(kout,100)
      ELSEIF (n < 2) THEN
         WRITE(kout,101)

* disallow negative grid values
      ELSEIF(x(1) < 0.) THEN
         WRITE(kout,105)

* check sorting
      ELSEIF( any( x(1:n-1) >= x(2:n) ) ) THEN
         WRITE(kout,110)
      ELSE
         ok = .TRUE.
      ENDIF

  100 FORMAT('Number of data exceeds dimension')
  101 FORMAT('Too few data, number of data points must be >= 2')
  105 FORMAT('Grid cannot start below zero')
  110 FORMAT('Grid is not sorted or contains multiple values')
*_______________________________________________________________________

      END SUBROUTINE gridck

      END MODULE GRIDS
