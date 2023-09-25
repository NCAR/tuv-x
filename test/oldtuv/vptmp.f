      MODULE COLUMN_TEMP

      IMPLICIT NONE

      private
      public :: vptmp

      contains

      SUBROUTINE vptmp(z,tlev,tlay)
*-----------------------------------------------------------------------------*
*   NAME: Vertical Profile of TeMPerature
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of temperatures.  Temperature values are      =*
*=  needed to compute some cross sections and quantum yields.  Distinguish   =*
*=  between temperature at levels and layers.                                =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  TLEV    - REAL, temperature (K) at each specified altitude level      (O)=*
*=  TLAY    - REAL, temperature (K) at each specified altitude layer      (O)=*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kin, kout

* input: (altitude working grid)
      REAL, intent(in) :: z(:)

* output:
      REAL, intent(out) :: tlev(:), tlay(:)

* local:
      INTEGER :: nz, nztmp
      INTEGER :: i, nd, istat
      REAL    :: altitude, temperature
      REAL, allocatable :: zd(:), td(:)

      interface
        FUNCTION inter1(xtarget, xsrc,ysrc) result( ytarget )
          REAL, intent(in) :: xtarget(:)
          REAL, intent(in) :: xsrc(:), ysrc(:)
          REAL :: ytarget(size(xtarget))
        END FUNCTION inter1
      end interface
*_______________________________________________________________________


* read in temperature profile

      WRITE(kout,*) 'air temperature: USSA, 1976'
      allocate( zd(0),td(0) )
      OPEN(kin,FILE='odat/DATAE1/ATM/ussa.temp',STATUS='old')
      DO i = 1, 3
         READ(kin,*)
      ENDDO

      DO
         READ(kin,*,iostat=istat) altitude, temperature
         if( istat /= 0 ) then
           exit
         endif
         zd = [zd,altitude] ; td = [td,temperature]
      ENDDO

      CLOSE(kin)

* use constant temperature to infinity:  

      nd = size(zd)
      zd(nd) = 1.E10

* alternative input temperature data could include, e.g., a read file here:

***********
*********** end data input.

* interpolate onto z-grid
      nz = size(z)
      tlev = inter1(z, zd,td)

      write(*,*) 'vptmp: data z grid'
      write(*,'(1p10g15.7)') zd
      write(*,*) ' '
      write(*,*) 'vptmp: Temp on data z grid'
      write(*,'(1p10g15.7)') td
      write(*,*) ' '
      write(*,*) 'vptmp: Temp on mdl z grid edges'
      write(*,'(1p10g15.7)') tlev
      write(*,*) ' '

* compute layer-averages
      tlay(1:nz-1) = .5*(tlev(2:nz) + tlev(1:nz-1))
c     tlay(nz) = tlay(nz-1)
      
      END SUBROUTINE vptmp

      END MODULE COLUMN_TEMP
