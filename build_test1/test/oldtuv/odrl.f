      module rayleigh

      IMPLICIT NONE

      private
      public :: odrl

      contains

      FUNCTION odrl(wc, c) result(dtrl)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  WC      - REAL, vector of lower limits of wavelength intervals in     (I)=*
*=            working wavelength grid                                        =*
*=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
*=            altitude layer                                                 =*
*=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
*=            and each specified wavelength                                  =*
*-----------------------------------------------------------------------------*

* input: 
      REAL, intent(in) ::  wc(:), c(:)

* output:
* Rayleigh optical depths

      REAL              :: dtrl(size(c),size(wc))

* other:

      REAL    :: srayl, wmicrn, xx 
      INTEGER :: iw
      
* compute Rayleigh cross sections and depths:

      DO iw = 1,size(wc)
* Rayleigh scattering cross section from WMO 1985 (originally from
* Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
* An empirical formula for its calculation in the homoshpere, Planet.
* Space Sci., 32, 1467-1468, 1984.
         wmicrn =  1.e-3*wc(iw)
         IF( wmicrn <= 0.55) THEN
            xx = 3.6772 + 0.389*wmicrn + 0.09426/wmicrn
         ELSE
            xx = 4. + 0.04
         ENDIF
         srayl = 4.02e-28/(wmicrn)**xx
         dtrl(:,iw) = c(:)*srayl
      ENDDO

      END FUNCTION odrl

      end module rayleigh
