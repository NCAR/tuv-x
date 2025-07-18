
      module CLOUD

      IMPLICIT NONE

      private
      public :: setcld

      contains

      SUBROUTINE setcld(taucld,zbase,ztop,
     $                  z,wl,dtcld,omcld,gcld)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set cloud properties for each specified altitude layer.  Properties      =*
*=  may be wavelength dependent.                                             =*
*=  Assumes horizontally infinite homogeneous cloud layers.
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  NW      - INTEGER, number of specified intervals + 1 in working       (I)=*
*=            wavelength grid                                                =*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTCLD   - REAL, optical depth due to absorption by clouds at each     (O)=*
*=            altitude and wavelength                                        =*
*=  OMCLD   - REAL, single scattering albedo due to clouds at each        (O)=*
*=            defined altitude and wavelength                                =*
*=  GCLD    - REAL, cloud asymmetry factor at each defined altitude and   (O)=*
*=            wavelength                                                     =*
*-----------------------------------------------------------------------------*

***** input

* new total cloud optical depth:
      REAL, intent(in) :: taucld
      REAL, intent(in) :: zbase, ztop

* (grids)
      REAL, intent(in) :: wl(:)
      REAL, intent(in) :: z(:)

***** Output: 
      REAL, intent(out) :: dtcld(:,:), omcld(:,:), gcld(:,:)

      INTEGER :: nz, nlyr
      INTEGER :: nbins

***** specified default data:
      REAL, allocatable :: zd(:), data_dtcld(:), omd(:), gd(:)
      REAL, allocatable :: womd(:), wgd(:)

* other:

      REAL :: mdl_dtcld(size(z)-1)
      REAL :: omz(size(z)-1)
      REAL :: gz(size(z)-1)
      INTEGER :: i, iw, n
*_______________________________________________________________________

* Set up clouds:
* All clouds are assumed to be infinite homogeneous layers
* Can have different clouds at different altitudes.
*   If multiple cloud layers are specified, non-cloudy layers
*   between them (if any) must be assigned zero optical depth.
* Set cloud optical properties:
*   data_dtcld(i) = optical depth of i_th cloudy layer
*   omd(i) = singel scattering albedo of i_th  cloudy layer
*   gd(i) = asymmetry factorof i_th  cloudy layer
* Cloud top and bottom can be set to any height zd(i), but if they don't
* match the z-grid (see subroutine gridz.f), they will be interpolated to
* the z-grid.

* Example:  set two separate cloudy layers:
*  cloud 1:  
*     base = 4 km
*     top  = 7 km
*     optical depth = 20.  (6.67 per km)
*     single scattering albedo = 0.9999
*     asymmetry factor = 0.85
*  cloud 2:
*     base = 9 km
*     top  = 11 km
*     optical depth = 5.  (2.50 per km)
*     single scattering albedo = 0.99999
*     asymmetry factor = 0.85

      nbins = size(wl) - 1
      nz = size(z) ;  nlyr = nz - 1

* cloud 1

      zd = (/ zbase,ztop /)
      n = size(zd)
      allocate( data_dtcld(n-1), omd(n-1), gd(n-1) )
      data_dtcld(1:n-1) = taucld
      omd(1:n-1) = .9999
      gd(1:n-1) = .85
      allocate( womd(n-1), wgd(n-1) )

******************
* compute integrals and averages over grid layers:
* for g and omega, use averages weighted by optical depth

      womd = omd * data_dtcld
      wgd  = gd  * data_dtcld

      CALL inter3(nz,z,mdl_dtcld,  n, zd,data_dtcld, 0)
      CALL inter3(nz,z,omz, n, zd,womd, 0)
      CALL inter3(nz,z,gz , n, zd,wgd, 0)

      WHERE (mdl_dtcld(1:nlyr) > 0.)
        omz(1:nlyr) = omz(1:nlyr)/mdl_dtcld(1:nlyr)
        gz(1:nlyr)  = gz(1:nlyr) /mdl_dtcld(1:nlyr)
      ELSEWHERE
        omz(1:nlyr) = 1.
        gz(1:nlyr) = 0.
      ENDWHERE
      
* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO iw = 1, nbins
        dtcld(1:nlyr,iw) = mdl_dtcld(1:nlyr)
        omcld(1:nlyr,iw) = omz(1:nlyr)
        gcld (1:nlyr,iw) = gz(1:nlyr)
      ENDDO

      END SUBROUTINE setcld

      END module CLOUD
