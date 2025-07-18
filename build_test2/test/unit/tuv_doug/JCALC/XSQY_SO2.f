      SUBROUTINE XSQY_SO2(nw,wl,wc,nz,tlev,airden,j,sq,jlabel,pn)
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide the product (cross section) x (quantum yield) for photolysis:    !
!           SO2 + hv -> Products                                            !
!                                                                           !
!  Cross section from Mike Mills, CU/LASP, Base on:                         !
!  1. Yung, Y.L., and W.B. Demore (1982) Photochemistry of the Stratosphere !
!  of Venus: Implications for Atmospheric Evolution, Icarus, 51, 199-247.   !
!  2. Okabe, H. In Photochemistry of Small Molecules; John Wiley and Sons   !
!     Inc.: New York, 1978; pp 248-249                                      !
!                                                                           !
!  Quantum yield = 1.0                                                      !
!---------------------------------------------------------------------------!
!  PARAMETERS:                                                              !
!  NW     - INTEGER, number of specified intervals + 1 in working        (I)!
!           wavelength grid                                                 !
!  WL     - REAL, vector of lower limits of wavelength intervals in      (I)!
!           working wavelength grid                                         !
!  WC     - REAL, vector of center points of wavelength intervals in     (I)!
!           working wavelength grid                                         !
!  NZ     - INTEGER, number of altitude levels in working altitude grid  (I)!
!  TLEV   - REAL, temperature (K) at each specified altitude level       (I)!
!  AIRDEN - REAL, air density (molec/cc) at each altitude level          (I)!
!  J      - INTEGER, counter for number of weighting functions defined  (IO)!
!  SQ     - REAL, cross section x quantum yield (cm^2) for each          (O)!
!           photolysis reaction defined, at each defined wavelength and     !
!           at each defined altitude level                                  !
!  JLABEL - CHARACTER*50, string identifier for each photolysis reaction (O)!
!           defined                                                         !
!---------------------------------------------------------------------------!
      IMPLICIT NONE
      INCLUDE 'params'

!---------------------------------------------------------------------------!
!     ... input                                                             !
!---------------------------------------------------------------------------!
      real, intent(in) :: wl(kw)
      real, intent(in) :: wc(kw)
      real, intent(in) :: tlev(kz)
      real, intent(in) :: airden(kz)

      integer, intent(in) :: nz
      integer, intent(in) :: nw

      character*80, intent(in) ::  pn
      character*60, intent(out) :: jlabel(kj)
      real, intent(out) :: sq(kj,kz,kw)

!---------------------------------------------------------------------------!
!     ... input/output                                                      !
!---------------------------------------------------------------------------!
      integer, intent(inout) :: j

!---------------------------------------------------------------------------!
!     ... local                                                             !
!---------------------------------------------------------------------------!
      integer kdata
      parameter (kdata=300)
      integer  i, n, ierr, iw
      real x_min(kdata), x_max(kdata), x(kdata), y(kdata)
      real yg(kw)
      real qy

!-----------------------------------------------
!     ... SO2 photolysis 
!-----------------------------------------------
      j = j+1
      jlabel(j) = 'SO2 + hv -> SO + O'

!-----------------------------------------------
!     ... SO2 cross sections
!----------------------------------------------
      OPEN(UNIT=kin,FILE=TRIM(pn)//'XS_SO2_mills.txt',
     $     STATUS='old')
      DO i = 1, 13
         READ(kin,*)
      ENDDO
      n = 125
      DO i = 1, n
         READ(kin,*) x_min(i), x_max(i), y(i)
         x(i) = (x_min(i)+x_max(i)) / 2.0
      ENDDO
      
      CLOSE(kin)

      CALL addpnt(x,y,kdata,n,x(1)*(1.-deltax),0.)
      CALL addpnt(x,y,kdata,n,               0.,0.)
      CALL addpnt(x,y,kdata,n,x(n)*(1.+deltax),0.)
      CALL addpnt(x,y,kdata,n,           1.e+38,0.)
      CALL inter2(nw,wl,yg,n,x,y,ierr)
      IF (ierr .NE. 0) THEN
         WRITE(*,*) ierr, jlabel(j)
         STOP
      ENDIF

!-----------------------------------------------
!     ... combine
!-----------------------------------------------
      qy = 1.0

      DO iw = 1, nw - 1
         DO i = 1, nz
            sq(j,i,iw) = yg(iw) * qy
         ENDDO
      ENDDO

      end subroutine XSQY_SO2
