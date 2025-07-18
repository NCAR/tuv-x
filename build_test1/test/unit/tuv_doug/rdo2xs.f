      SUBROUTINE rdo2xs(nw,wl,wc,o2xs1,pn)
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Read O2 absorption cross section.  Except the SR bands and L-alpha line  !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!---------------------------------------------------------------------------!
!  EDIT HISTORY:                                                            !
!  02/02  By Xuexi                                                          !
!---------------------------------------------------------------------------!
      IMPLICIT NONE
      INCLUDE 'params'

!-----------------------------------------------------------------------------!
!     ... input                                                               !
!-----------------------------------------------------------------------------!
      real, intent(in)    :: wl(kw)
      real, intent(in)    :: wc(kw)
      integer, intent(in) :: nw

!-----------------------------------------------------------------------------!
!     ... output                                                              !
!-----------------------------------------------------------------------------!
      real, intent(out)   ::  o2xs1(kw)

!... Internal

      integer    i, iw, n, kdata, ierr
      parameter (kdata = 200)
      real       x1(kdata), y1(kdata)
      real       x, y
      character*80 pn

!------------------------------------------------------------------------
!   NOTE: Output O2 xsect, is temporary and will be over-written in 
!          Lyman-alpha and Schumann-Runge wavelength bands.
!------------------------------------------------------------------------
!     ... data                                                                
!------------------------------------------------------------------------
!    Read O2 absorption cross section data:
!     116.65 to 203.05 nm = from Brasseur and Solomon 1986
!     205 to 240 nm = Yoshino et al. 1988 (same as JPL06)
!
!    Note that subroutine seto2.f will over-write values in the 
!       spectral regions corresponding to:
!       Lyman-alpha (LA: 121.4-121.9 nm, Chabrillat and Kockarts
!                     parameterization 
!       Schumann-Runge bands (SRB: 174.4-205.8 nm, Koppers 
!                     parameteriaztion)
!-----------------------------------------------------------------------
      n = 0

      OPEN(UNIT=kin,FILE=Trim(pn)//'XS_O2_brasseur.txt')
      DO i = 1, 7
         READ(kin,*)
      ENDDO
      DO i = 1, 78
         READ(kin,*) x, y
         IF (x .LE. 204.) THEN
            n = n + 1
            x1(n) = x
            y1(n) = y
!            print*, x1(n), y1(n)
         ENDIF
      ENDDO
      CLOSE(kin)

      OPEN(UNIT=kin,
     $     FILE=Trim(pn)//'XS_O2_yoshino.txt',STATUS='old')
      DO i = 1, 8
         READ(kin,*)
      ENDDO
      DO i = 1, 36
         n = n + 1
         READ(kin,*) x, y
         y1(n) = y*1.E-24
         x1(n) = x
!         print*, x1(n), y1(n)
      END DO
      CLOSE (kin)
      
!-----------------------------------------------------------------------------
!        Add termination points and interpolate onto the 
!        user grid (set in subroutine gridw):
!-----------------------------------------------------------------------------
      CALL addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),y1(1))
      CALL addpnt(x1,y1,kdata,n,0.               ,y1(1))
      CALL addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      CALL addpnt(x1,y1,kdata,n,              1.E+38,0.)

!      CALL inter2(nw,wl,o2xs1, n,x1,y1, ierr)
      print*, "* interp4 used in rdo2xs.f"
      ierr = 0
      CALL inter4(nw,wl,o2xs1, n+1,x1,y1, ierr)
!---------------------------------------------------------------
!     ... Check routine 
!      do iw = 1,51
!        print*, iw, wc(iw), o2xs1(iw)
!      enddo
!      stop
!---------------------------------------------------------------  
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, 'O2 -> O + O'
         STOP
      ENDIF

      end subroutine rdo2xs






