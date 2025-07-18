      MODULE COLUMN_ATM

      IMPLICIT NONE

      private
      public :: vpair

      contains

      SUBROUTINE vpair(psurf, z, con, col, col_lasrb)
*-----------------------------------------------------------------------------*
*=  NAME:  Vertial Profile of AIR
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of air molecules.  Subroutine includes a      =*
*=  shape-conserving scaling method that allows scaling of the entire        =*
*=  profile to a given sea-level pressure.                                   =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  PSURF   - REAL, surface pressure (mb) to which profile should be      (I)=*
*=            scaled.  If PSURF < 0, no scaling is done                      =*
*=  NZ      - INTEGER, number of specified altitude levels in the working (I)=*
*=            grid                                                           =*
*= outputs are on z-grid:
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  CON     - REAL, air density (molec/cc) at each specified altitude     (O)=* 
*=  COL     - REAL, number of air molecules per cm^2 in each specified    (O)=*
*=            altitude layer (column vertical increment                      =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kin, kout

********* input:
    
      REAL, intent(in) ::  psurf
      REAL, intent(in) :: z(:)

********* output:
* con(iz) = air density (molec cm-3) at level iz
* col(iz) =  column amount (molec cm-2) for layer iz      
      REAL, intent(out) :: con(:)
      REAL, intent(out) :: col(:)
      REAL, intent(out) :: col_lasrb(:)

* specified air profile data:
      INTEGER :: nd
      REAL    :: hscale
      REAL, allocatable :: zd(:), air(:), cd(:)


* local:
      REAL, PARAMETER :: pconv = 980.665 * 1.E-3 * 28.9644 / 6.022169E23
      REAL, PARAMETER :: km2cm = 1.e5

      INTEGER :: nz, nlyr
      REAL :: scale
      REAL :: pold
      REAL :: altitude, number_density_air

* other:
      INTEGER :: i, istat
      REAL, allocatable :: airlog(:), conlog(:)

      interface
        FUNCTION inter1(xtarget, xsrc,ysrc) result( ytarget )
          REAL, intent(in) :: xtarget(:)
          REAL, intent(in) :: xsrc(:), ysrc(:)
          REAL :: ytarget(size(xtarget))
        END FUNCTION inter1
      end interface

*_______________________________________________________________________

* The objective of this subroutine is to take the input air profile and interpolate
*   it to the working grid, z(i = 1, nz).  The desired outputs are con(i) and col(i).
* Input vertical profiles can be specified in various ways. For example: 
*   altitude vs. concentration (molecules cm-3)
*   altitude vs. pressure (mbar) and temperature (K)
*   altitude vs. column increments (molec cm-2)
* The interpolation scheme will depend on the specific type of input data.
* Here, the US Standard Atmosphere is given as altitude vs. concentration (also called
*   number density)

* _________SECTION 1:  Read in vertical profile of concentration

      WRITE(kout,*) 'air concentrations: USSA, 1976'

      allocate( zd(0), air(0) )

      OPEN(kin,FILE='odat/DATAE1/ATM/ussa.dens',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO

      nd = 1
      DO
        READ(kin,*,iostat=istat) altitude, number_density_air
        if( istat /= 0 ) then
          exit
        endif
        zd = [zd,altitude] ; air = [air,number_density_air]
      ENDDO

      CLOSE(kin)
* add 1 meter to top, to avoid interpolation end-problem if z-grid is up to 120 km
      nd = size(zd)
      zd(nd) = zd(nd) + 0.001

* scale height, km, at top of atmosphere:
      hscale = 8.01

********************** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  
* For air, this is best done using logarithms of concentration.
*   take logs
*   if z-grid extends beyond available data, stop (no extrapolation allowed)
*   interpolate log of air(nd) onto z grid 
*   re-exponentiate to get gridded concentrations

      airlog = LOG(air)

      nz = size(z)
      nlyr = nz - 1
      IF(z(nz) > zd(nd)) THEN
        STOP 'in vpair: ztop < zdata'
      ENDIF

      allocate( conlog(nz) )
      conlog = inter1(z, zd,airlog)

      con = EXP(conlog)

* Find gridded column increments in z-grid:
*   use log intergration

*  replaced logarithmic integral by geometric average, because of precision
*   problems with very thin layers (e.g. 1 cm in snow)         
* - sm, jan 2018
      col(1:nlyr) = km2cm * (z(2:nz) - z(1:nlyr))
     $                    * sqrt(con(2:nz))*sqrt(con(1:nlyr))

* Add exponential tail integral at top of atmosphere:
*   this is folded into layer nz-1, making this layer "heavy'.  
*   The layer nz is not used. The radiative transfer 
*   calculation is based on nz-1 layers (not nz).

      col_lasrb(1:nlyr) = col(1:nlyr)
      col_lasrb(nz) = km2cm * hscale * con(nz)
      col(nlyr) = col(nlyr) + col_lasrb(nz)
  
* Scale by input surface pressure:
* min value = 1 molec cm-2

      IF(psurf > 0.) THEN
        pold  = pconv * MAX( sum(col(1:nlyr)),1. )
        scale = psurf/pold
        col(:) = col(:) * scale
        con(:) = con(:) * scale
        col_lasrb(:) = col_lasrb(:) * scale
      ENDIF

      END SUBROUTINE vpair

      END MODULE COLUMN_ATM
