      MODULE ALBEDO_MOD

      IMPLICIT NONE

      private
      public :: setalb

      contains
*=============================================================================*

      FUNCTION setalb(albnew,wl) result(albedo)

*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set the albedo of the surface.  The albedo is assumed to be Lambertian,  =*
*=  i.e., the reflected light is isotropic, and independent of direction     =*
*=  of incidence of light.  Albedo can be chosen to be wavelength dependent. =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  ALBEDO  - REAL, surface albedo at each specified wavelength           (O)=*
*-----------------------------------------------------------------------------*


* input: (wavelength working grid data)
      REAL, intent(in) :: albnew
      REAL, intent(in) :: wl(:)

* output:
      REAL :: albedo(size(wl)-1)

      albedo(:) = albnew

* alternatively, can input wavelenght-dependent values if avaialble.
*_______________________________________________________________________

      END FUNCTION setalb

      END MODULE ALBEDO_MOD
