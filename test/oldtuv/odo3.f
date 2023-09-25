      MODULE O3_OPTICAL_DEPTH

      IMPLICIT NONE

      private
      public :: odo3

      contains

*=============================================================================*

      FUNCTION odo3(z,wl,o3xs,c) result(dto3)

*-----------------------------------------------------------------------------*
*=  NAME:  Optical Depths of O3
*=  PURPOSE:                                                                 =*
*=  Compute ozone optical depths as a function of altitude and wavelength    +*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)=*
*=           each specified wavelength and altitude                          =*
*=  C      - REAL, ozone vertical column increments, molec cm-2, for each (I)=*
*=           layer                                                           =*
*=  DTO3   - REAL, optical depth due to ozone absorption at each          (O)=*
*=           specified altitude at each specified wavelength                 =*
*-----------------------------------------------------------------------------*

********
* input:
********

* grids:

      REAL, intent(in) :: wl(:)
      REAL, intent(in) :: z(:)

* ozone absorption cross section, functions of wavelength and altitude

      REAL, intent(in) :: o3xs(:,:)

* ozone vertical column increments

      REAL, intent(in) :: c(:)

********
* output:
********

      REAL              :: dto3(size(z)-1,size(wl)-1)

********
* internal:
********

      INTEGER :: iw

*_______________________________________________________________________

* calculate ozone optical depth for each layer, with temperature 
* correction.  Output, dto3(kz,kw)

      DO iw = 1, size(wl)-1
         dto3(:,iw) = c(:) * o3xs(:,iw)
      ENDDO

      END FUNCTION odo3

      END MODULE O3_OPTICAL_DEPTH
