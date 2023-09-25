      MODULE SET_NO2_OPTICAL_DTAU

      IMPLICIT NONE

      private
      public :: setno2

      contains

      FUNCTION setno2(ipbl, zpbl, xpbl, 
     $                no2new, z, nbins, no2xs, 
     $                tlay, dcol) result(dtno2)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of NO2 molecules, and corresponding absorption=*
*=  optical depths.  Subroutine includes a shape-conserving scaling method   =*
*=  that allows scaling of the entire profile to a given overhead NO2        =*
*=  column amount.                                                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NO2NEW - REAL, overhead NO2 column amount (molec/cm^2) to which       (I)=*
*=           profile should be scaled.  If NO2NEW < 0, no scaling is done    =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  NO2XS  - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)=*
*=           each specified wavelength                                       =*
*=  TLAY   - REAL, temperature (K) at each specified altitude layer       (I)=*
*=  DTNO2  - REAL, optical depth due to NO2 absorption at each            (O)=*
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
      REAL, intent(in) :: no2new

* absorption cross sections 
      REAL, intent(in) :: no2xs(:,:)

* mid-layer temperature, layer air column
      REAL, intent(in) :: tlay(:), dcol(:)

********
* output:
********

      REAL :: dtno2(size(z)-1,nbins)

********
* local:
********

* nitrogen dioxide profile data:

      INTEGER, parameter :: nd = 3
      REAL, parameter    :: km2cm = 1.e5
      REAL, parameter    :: hscale = 4.5*km2cm
      REAL, parameter    :: ppb = 1.e-9
      REAL, parameter    :: ppt = 1.e-12

      REAL :: cz(size(z)-1)
      REAL :: zd(nd), no2(nd), cd(nd-1)
      REAL :: colold, scale
      REAL :: sno2

* other:

      INTEGER :: nz
      INTEGER :: i, l

*_______________________________________________________________________
* Data input:

* Example:  set to 1 ppb in lowest 1 km, set to zero above that.
* - do by specifying concentration at 3 altitudes.
      zd(1:2) = (/ 0.,1. /)
      zd(3)   = zd(2)* 1.000001
      no2(:) = (/ 2.69e10, 2.69e10, 10./largest /)

* compute column increments (alternatively, can specify these directly)

      DO i = 1, nd - 1
         cd(i) = (no2(i+1)+no2(i)) * km2cm * .5 * (zd(i+1)-zd(i))
      ENDDO

* Include exponential tail integral from top level to infinity.
* fold tail integral into top layer
* specify scale height near top of data (use ozone value)

      cd(nd-1) = cd(nd-1) + hscale * no2(nd)

***********
*********** end data input.

      nz = size(z)
* Compute column increments and total column on standard z-grid.  

      CALL inter3(nz,z,cz, nd,zd,cd, 1)

**** Scaling of vertical profile by ratio of new to old column:
* If old column is near zero (less than 1 molec cm-2), 
* use constant mixing ratio profile (nominal 1 ppt before scaling) 
* to avoid numerical problems when scaling.

      IF(sum(cz(:)) < 1.) THEN
         DO i = 1, nz-1
            cz(i) = ppt * dcol(i)
         ENDDO
      ENDIF
      colold = sum(cz(:))
      scale =  2.687e16 * no2new / colold

      cz(:) = cz(:) * scale

*! overwrite for specified pbl height

      IF(ipbl > 0) THEN
         write(*,*) 'pbl NO2 = ', xpbl, ' ppb'
         DO i = 1, nz-1
            IF (i .LE. ipbl) THEN
               cz(i) = xpbl * ppb * dcol(i)
            ELSE
               cz(i) = 0.
            ENDIF
         ENDDO
      ENDIF

************************************
* calculate optical depth for each layer.  Output: dtno2(kz,kw)

      DO l = 1, nbins
         dtno2(:,l) = cz(:)*no2xs(:,l)
      ENDDO

      END FUNCTION setno2

      END MODULE SET_NO2_OPTICAL_DTAU
