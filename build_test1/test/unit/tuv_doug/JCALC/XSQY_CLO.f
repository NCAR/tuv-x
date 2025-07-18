      subroutine XSQY_CLO(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield):                        !
!           ClO + hv -> Cl + O                                                !
!   cross section: JPL06                                                      !
!   quantum yield: is unity.                                                  !
!-----------------------------------------------------------------------------!
!   parameters:                                                               !
!   nw     - integer, number of specified intervals + 1 in working        (i) !
!            wavelength grid                                                  !
!   wl     - real, vector of lower limits of wavelength intervals in      (i) !
!            working wavelength grid                                          !
!   wc     - real, vector of center points of wavelength intervals in     (i) !
!            working wavelength grid                                          !
!   nz     - integer, number of altitude levels in working altitude grid  (i) !
!   tlev   - real, temperature (k) at each specified altitude level       (i) !
!   airlev - real, air density (molec/cc) at each altitude level          (i) !
!   j      - integer, counter for number of weighting functions defined  (io) !
!   sq     - real, cross section x quantum yield (cm^2) for each          (o) !
!            photolysis reaction defined, at each defined wavelength and      !
!            at each defined altitude level                                   !
!   jlabel - character*60, string identifier for each photolysis reaction (o) !
!            defined                                                          !
!-----------------------------------------------------------------------------!
!   edit history:                                                             !
!   07/27/07  Doug Kinnison                                                   !
!-----------------------------------------------------------------------------!
      implicit none
      include 'params'

!-----------------------------------------------------------------------------!
!     ... input                                                               !
!-----------------------------------------------------------------------------!
      real, intent(in) :: wl(kw)
      real, intent(in) :: wc(kw)
      real, intent(in) :: tlev(kz)
      real, intent(in) :: airlev(kz)

      integer, intent(in) :: nz
      integer, intent(in) :: nw

      character*80, intent(in)  ::  pn
      character*60, intent(out) :: jlabel(kj)
      real, intent(out)         :: sq(kj,kz,kw)

!-----------------------------------------------------------------------------!
!     ... input/output                                                        !
!-----------------------------------------------------------------------------!
      integer, intent(inout) :: j

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=300)
      integer i, iw, n, idum, ierr, iz
      real x1(kdata)
      real y1(kdata)
      real yg(kw)
      real qy
 
!----------------------------------------------
!     ... jlabel(j) = 'ClO + hv -> Cl + O'
!----------------------------------------------
      j = j+1
      jlabel(j) = 'ClO + hv -> Cl + O'

!----------------------------------------------------
!     ... cross sections from JPL06 recommendation
!----------------------------------------------------
      open(kin,file=TRIM(pn)//'XS_CLO_JPL06.txt',status='old')

      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n
        read(kin,*) x1(i), y1(i)
      enddo
      close(kin)
      
      call addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n,          0.,0.)
      call addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n,        1e38,0.)

      call inter2(nw,wl,yg,n,x1,y1,ierr)

      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif
!-------------------------------------------------------
!     ... quantum yield (assumed) to be unity (JPL06)
!-------------------------------------------------------
      qy = 1.0

      do iw = 1, nw-1
        do iz = 1, nz
          sq(j,iz,iw) = qy * yg(iw)
 
        enddo
      enddo
!-------------------------------------------------------
!     ... Check routine (no temperature dependence
!      print*,'jclo'
!      do iw = 30, 72
!        print*, iw, wc(iw), (qy * yg(iw))
!      enddo
!      stop
!-------------------------------------------------------
      
      end subroutine XSQY_CLO
