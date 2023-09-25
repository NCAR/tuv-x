      MODULE SET_SO2_OPTICAL_DTAU

      IMPLICIT NONE

      private
      public :: setso2

      contains

      FUNCTION setso2(ipbl, zpbl, xpbl,
     $     so2new, z, nbins, so2xs, tlay, dcol) result(dtso2)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of SO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead SO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  SO2NEW - REAL, overhead SO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If SO2NEW < 0, no scaling is done    =*
*=  NZ     - INTEGER, number of specified altitude levels in the working  (I)=*
*=           grid                                                            =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  NW     - INTEGER, number of specified intervals + 1 in working        (I)=*
*=           wavelength grid                                                 =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  SO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : largest

********
* input:
********
      INTEGER, intent(in) :: ipbl
      INTEGER, intent(in) :: nbins
      REAL, intent(in)    :: zpbl, xpbl

* grids:
      REAL, intent(in) :: z(:)

      REAL, intent(in) :: so2new
* absorption cross sections 
      REAL, intent(in) :: so2xs(:)

* mid-layer temperature and layer air column
      REAL, intent(in) :: tlay(:), dcol(:)

********
* output:
********

      REAL :: dtso2(size(z)-1,nbins)

********
* local:
********

      REAL, parameter :: km2cm = 1.e5
      REAL, parameter :: hscale = 4.5*km2cm
      REAL, parameter :: ppb = 1.e-9
      REAL, parameter :: ppt = 1.e-12

      INTEGER :: nz
      REAL    :: cz(size(z)-1)

* sulfur dioxide profile data:
      INTEGER, parameter :: nd = 3

      REAL :: zd(nd), so2(nd), cd(nd-1)
      REAL :: colold, scale

* other:

      INTEGER :: i, l

      nz = size(z)
*_______________________________________________________________________
* Data input:

* Example:  set to 1 ppb in lowest 1 km, set to zero above that.
* - do by specifying concentration at 3 altitudes.
      zd(1:2) = (/ 0.,1. /)
      zd(3)   = zd(2)* 1.000001
      so2(:) = (/ 2.69e10, 2.69e10, 10./largest /)

* compute column increments (alternatively, can specify these directly)
      DO i = 1, nd - 1
         cd(i) = (so2(i+1)+so2(i)) * km2cm * .5 * (zd(i+1)-zd(i))
      ENDDO

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      cd(nd-1) = cd(nd-1) + hscale * so2(nd)

***********
*********** end data input.

* Compute column increments on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

**** Scaling of vertical profile by ratio of new to old column:
* If old column is near zero (less than 1 molec cm-2), 
* use constant mixing ratio profile (nominal 1 ppt before scaling) 
* to avoid numerical problems when scaling.

      IF(sum(cz(:)) < 1.) THEN
         cz(:) = ppt * dcol(:)
      ENDIF
      colold = sum(cz(1:nz-1))
      scale =  2.687e16 * so2new / colold
      cz(:) = cz(:) * scale

*! overwrite for specified pbl height, set concentration here

      IF(ipbl > 0) THEN
         write(*,*) 'pbl SO2 = ', xpbl, ' ppb'
         DO i = 1, nz-1
            IF (i <= ipbl) THEN
               cz(i) = xpbl * ppb * dcol(i)
            ELSE
               cz(i) = 0.
            ENDIF
         ENDDO
      ENDIF

************************************
* calculate sulfur optical depth for each layer, with optional temperature 
* correction.  Output, dtso2(kz,kw)

      DO l = 1, nbins
         dtso2(:,l) = cz(:)*so2xs(l)
      ENDDO

      END FUNCTION setso2

      END MODULE SET_SO2_OPTICAL_DTAU
