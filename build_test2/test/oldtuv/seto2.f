      MODULE SET_O2_OPTICAL_TAU

      IMPLICIT NONE

      private
      public :: seto2

      contains

      FUNCTION seto2(z,wl, cz, o2xs1) result(dto2)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of air molecules.  Subroutine includes a      =*
*=  shape-conserving scaling method that allows scaling of the entire        =*
*=  profile to a given sea-level pressure.                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  WL      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=            and each specified wavelength                                  =*
*=  CZ      - REAL, number of air molecules per cm^2 at each specified    (O)=*
*=            altitude layer                                                 =*
*-----------------------------------------------------------------------------*

* input: (grids)

      REAL, intent(in) :: wl(:)
      REAL, intent(in) :: z(:)
      REAL, intent(in) :: cz(:)
      REAL, intent(in) :: o2xs1(:)

* output:
*  O2 absorption optical depth per layer at each wavelength

      REAL :: dto2(size(z)-1,size(wl)-1)

      REAL, parameter :: o2vmr = .2095

      INTEGER :: iw

*_______________________________________________________________________
*  Assumes that O2 = 20.95 % of air density.  If desire different O2 
*    profile (e.g. for upper atmosphere) then can load it here.

      DO iw = 1, size(wl)-1
         dto2(:,iw) = o2vmr * cz(:) * o2xs1(iw)
      ENDDO  

      END FUNCTION seto2

      END MODULE SET_O2_OPTICAL_TAU
