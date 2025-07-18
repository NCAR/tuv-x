      subroutine XSQY_HNO3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product of (cross section) x (quantum yield) for photolysis       !
!         hno3 + hv -> oh + no2                                               !
!   cross section: burkholder et al., 1993 (and JPL06)                        !
!   quantum yield: assumed to be unity                                        !
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
!   05/98  original, adapted from former jspec1 subroutine                    !
!   01/15/08 minor update,dek                                                 !
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

      character*80, intent(in) ::  pn
      character*60, intent(out) :: jlabel(kj)
      real, intent(out) :: sq(kj,kz,kw)

!-----------------------------------------------------------------------------!
!     ... input/output                                                        !
!-----------------------------------------------------------------------------!
      integer, intent(inout) :: j

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=100)
      integer n1, n2
      integer i, iw, n, idum, iz
      integer ierr
      real x1 (kdata), x2 (kdata)
      real y1 (kdata), y2 (kdata)
      real yg1(kw),    yg2(kw)
      real yg( kw)
      real tin(nz)
 
!----------------------------------------------
!     ... tin set to tlev
!----------------------------------------------
      tin(:) = tlev(:)

!----------------------------------------------
!     ... jlabel(j) = 'HNO3 -> OH + NO2
!----------------------------------------------
       j = j + 1
       jlabel(j) = 'HNO3 + hv -> OH + NO2'

!-----------------------------------------------------------------------
!     ... hno3 cross section parameters from burkholder et al. 1993
!-----------------------------------------------------------------------
      open(kin,file=TRIM(pn)//'XS_HNO3_JPL06.txt',status='old')

!...  read lambda and cross sections
      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n
        read(kin,*) x1(i), y1(i)
      enddo

!...  read lambda and T-dep coeff.      
      read(kin,*)
      do i = 1, n
        read(kin,*) x2(i), y2(i)
      enddo
      close(kin)
 
      call addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n,          0.,0.)
      call addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n,        1e38,0.)
      call inter2(nw,wl,yg1,n,x1,y1,ierr)
      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif

      n= 80
      call addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
      call addpnt(x2,y2,kdata,n,               0.,0.)
      call addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),0.)
      call addpnt(x2,y2,kdata,n,            1.e+38,0.)
      call inter2(nw,wl,yg2,n,x2,y2,ierr)
      if (ierr .ne. 0) then
         write(*,*) ierr, jlabel(j)
         stop
      endif

!--------------------------------------------------
!     ... quantum yield = 1
!         correct for temperature dependence
!--------------------------------------------------
      do iw = 1, nw - 1
         do iz = 1, nz
            sq(j,iz,iw) = yg1(iw)
     $           * exp( yg2(iw)/1.e3*(tin(iz)-298.) )
         enddo
      enddo

!-------------------------------------------------------
!     ... Check routine (no temperature dependence
!      iz = 1
!      do iw = 29, 79
!        print*, iw, wc(iw), sq(j,iz,iw)
!      enddo
!      stop
!-------------------------------------------------------
      
      end subroutine XSQY_HNO3
