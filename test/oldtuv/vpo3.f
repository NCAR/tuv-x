      MODULE COLUMN_O3
*=============================================================================*

      IMPLICIT NONE

      private
      public :: vpo3

      contains

      FUNCTION vpo3(ipbl, zpbl, mr_pbl, to3new, z, aircol) result(col)
*-----------------------------------------------------------------------------*
*=  NAME:  Vertical Profiles of Ozone = vpo3                                 =*
*=  PURPOSE:                                                                 =*
*=  Computes O3 column increments, col(i), molec cm-2 for each layer i of    =* 
*=  the working grid z(i = 1, nz).                                           =*
*=  Normally, col(i) values are computed from input vertical profiles of     =*
*=  concentrations (molec cm-3), that are then interpolated and integrated   =*
*=  over each layer.  The default example here uses the US Standard          =*
*=  Atmosphere (1976) mid-latitude concentration profile.                    =*
*=  Users can substitute their own concentration profile, as long as         =*
*=  appropriate adjustments are made in Section 1 to input the data, and in  =*
*=  section 2 to interpolate to the working grid.                            =*
*=  A scale factor is provided to allow changing the total column amount,    =*
*=  but conserving the shape of the profile.                                 =*
*=  An option to insert PBL pollutants is provided                           =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  TO3NEW - REAL,  Dobson Units, new total column amount from surface       =* 
*=            to space, to be scaled.  If TO3NEW < 0, no scaling is done.    =*
*=  Z      - REAL, specified altitude working grid (km)                   (I)=*
*=  AIRCOL(KZ) = REAL, air column increment (molec cm-2), provided here in   =*
*=  case it is necessary to convert from mixing ratio units (e.g. ppb).      =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kin, kout, nzero, pzero

********
* inputs:
********
      INTEGER, intent(in) :: ipbl
      REAL, intent(in) :: to3new
      REAL, intent(in) :: zpbl, mr_pbl
      REAL, intent(in) :: z(:)
      REAL, intent(in) :: aircol(:)
********
* output:
********

      REAL    :: col(size(z)-1)

********
* internal
********
      REAL, parameter :: km2cm = 1.e5

      INTEGER :: nlyr, nz
      INTEGER :: i, nd
      REAL    :: hscale
      REAL    :: to3inp, scale
      REAL    :: rfact
      REAL    :: altitude, o3_density
      
      REAL    :: con(size(z))
      REAL, allocatable :: zd(:), xd(:)

      interface
        FUNCTION inter1(xtarget, xsrc,ysrc) result( ytarget )
          REAL, intent(in) :: xtarget(:)
          REAL, intent(in) :: xsrc(:), ysrc(:)
          REAL :: ytarget(size(xtarget))
        END FUNCTION inter1
      end interface

*-----------------------------------------------------------------------------*
* The objective of this subroutine is to calculate the vertical increments 
*   in the O3 column, for each layer of the working grid z(i = 1, nz).
* The input O3 profiles can be specified in different ways, and each case
*   will require careful consideration of the interpolation scheme.  Some
*   examples of possible input data are:
*    altitude vs. O3 concentration (number density), molec cm-3
*    altitude vs. O3 mixing ratio (e.g. parts per billion) relative to air
*    altitude vs. O3 column increments, molec cm-2, between specific altitudes.
* Special caution is required with mixed inputs, e.g. ppb in boundary layer and 
*   molec cm-3 above the boundary layer.
*-----------------------------------------------------------------------------*

*-----------------------------------------------------------------------------*
* _________SECTION 1:  Read in vertical profile of concentration
* Default is US Standard Atmosphere (1976)
* If a different vertical concentration profile is specified, the code
* in this section (Section 1) should be replaced accordingly
*-----------------------------------------------------------------------------*
      allocate( zd(0),xd(0) )
      WRITE(kout,*) 'ozone profile: USSA, 1976'
      OPEN(kin,FILE='odat/DATAE1/ATM/ussa.ozone',STATUS='old')
      DO i = 1, 7
        READ(kin,*)
      ENDDO
      nd = 39
      DO i = 1, nd
         READ(kin,*) altitude, o3_density
         zd = [zd,altitude] ; xd = [xd,o3_density]
      ENDDO
      CLOSE(kin)

*-----------------------------------------------------------------------------*
* Ussa data stop at 74 km.  Add values up to 121 km, 
* assuming exponential decay from 74 km up, with scale height of
*  4.5 km.
*-----------------------------------------------------------------------------*

      hscale = 4.5
      rfact  = EXP(-1./hscale)

      DO
        nd = nd + 1
        zd = [zd,zd(nd-1)+ 1.]
        xd = [xd,xd(nd-1) * rfact]
        IF(zd(nd) >= 121.) EXIT
      ENDDO

*********** end data input.

* ________SECTION 2:  Compute column increments on standard z-grid.  

* linear interpolation

      nz   = size(z)
      nlyr = nz - 1
      con = inter1(z, zd,xd )

      write(*,*) 'vpo3: data z grid'
      write(*,'(1p10g15.7)') zd
      write(*,*) ' '
      write(*,*) 'vpo3: o3 on data z grid'
      write(*,'(1p10g15.7)') xd
     
* compute column increments

      DO i = 1, nlyr
         col(i) = 0.5 * (con(i) + con(i+1)) * (z(i+1) - z(i)) * km2cm
      ENDDO

*-----------------------------------------------------------------------------*
* Add exponential tail integral at top of atmosphere:
*   this is folded into layer nz-1, making this layer "heavy'.  
*   The layer nz is not used. The radiative transfer 
*   calculation is based on nz-1 layers (not nz).
*-----------------------------------------------------------------------------*
      col(nlyr) = col(nlyr) + km2cm * hscale * con(nz)

*-----------------------------------------------------------------------------*
***** Scaling to new total ozone
* to3old = total o3 column, in Dobson Units, old value
* to3new = total o3 column, in Dobson Units, new value
*    (1 DU = 2.687e16)
* If to3new is not negative, scale to new total column value, to3new :
* (to3new = 0. is a possible input, to see effect of zero ozone)
*-----------------------------------------------------------------------------*
      IF (to3new > nzero) THEN
         to3inp = sum(col(1:nlyr))/2.687e16
         IF(to3inp < pzero) THEN
           STOP 'in vpo3: to3old is too small'
         ENDIF
         scale = to3new/to3inp
         DO i = 1, nlyr
            col(i) = col(i) * scale
            con(i) = con(i) * scale
         ENDDO
         con(nz) = con(nz) * scale
      ENDIF

      write(*,*) ' '
      write(*,*) 'vpo3: o3 on mdl z grid edges'
      write(*,'(1p10g15.7)') con
      write(*,*) ' '

*-----------------------------------------------------------------------------*
*! overwrite column increments for specified pbl height
* use mixing ratio in pbl
*-----------------------------------------------------------------------------*
      IF(ipbl > 0) THEN
         write(*,*) 'pbl O3 = ', mr_pbl, ' ppb'
         DO i = 1, nlyr
            IF (i <= ipbl) THEN
               col(i) = mr_pbl*1.E-9 * aircol(i)
            ENDIF
         ENDDO
      ENDIF

      END FUNCTION vpo3

      END MODULE COLUMN_O3
