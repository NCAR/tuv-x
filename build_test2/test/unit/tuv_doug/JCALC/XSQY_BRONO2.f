      subroutine XSQY_BRONO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for brono2 photolysis:  !
!         BrONO2 + hv -> products                                             !
!                                                                             !
!   cross section: jpl 06 recommendation                                      !
!   quantum yield: jpl 06 recommendation                                      !
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
!   07/27/07: Doug Kinnison                                                   !
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

      real, intent(out) :: sq(kj,kz,kw)

!-----------------------------------------------------------------------------!
!     ... input/output                                                        !
!-----------------------------------------------------------------------------!
      integer, intent(inout) :: j

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer kdata
      parameter(kdata=200)
      integer i, iw, n, n1, idum, ierr, iz
      real x1  (kdata), y1(kdata)
      real xin (kdata)
      real a1  (kdata), a2(kdata)
      real ytmp(nz,kdata)
      real ytd (nz,kw)
      real yg1 (kw)
      real tin (nz)
      real qy1
      real qy2
!-----------------------------------------------
!     ... tin set to tlev 
!-----------------------------------------------
      tin(:) = tlev(:)

!-----------------------------------------------
!     ... jlabel(j) = 'BrONO2 + hv -> Br + NO3'
!     ... jlabel(j) = 'BrONO2 + hv -> BrO + NO2'
!-----------------------------------------------
      j = j+1
      jlabel(j) = 'BrONO2 + hv -> Br + NO3'

!-----------------------------------------------
!     ... cross sections from JPL06 
!-----------------------------------------------
      open(kin,
     &     file=TRIM(pn)//'XS_BRONO2_JPL06.txt',status='old')
      
      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do iw = 1, n
        read(kin,*) xin(iw), y1(iw)
      enddo
      close(kin)

!-----------------------------------------------
!     ... Read in temperature dep coeff
!-----------------------------------------------
      open(kin,
     &     file=TRIM(pn)//'XS_BRONO2_td_JPL06.txt',status='old')
      
      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do iw = 1, n
        read(kin,*) xin(iw), a1(iw), a2(iw)
      enddo
      close(kin)

!-----------------------------------------------
!     ... derive T-dep (200-296K)
!-----------------------------------------------
      do  iz = 1, nz
        do iw = 1 , n
          if ((tin(iz) .GE. 200.) .AND. (tin(iz) .LE. 296.)) Then
            ytmp(iz,iw) = y1(iw) * 
     &                   ( 1. + 
     &                     A1(iw)* (tin(iz)-296.) +  
     &                     A2(iw)*((tin(iz)-296.)**2) 
     &                    )      
          endif
          if (tin(iz) .LT. 200.) then
             ytmp(iz,iw) = y1(iw) * 
     &                   ( 1. + 
     &                     A1(iw)* (200.-296.) +  
     &                     A2(iw)*((200.-296.)**2) 
     &                    )      
          endif
          if (tin(iz) .GT. 296.) then
             ytmp(iz,iw) = y1(iw) 
          endif
        enddo
      enddo
      
!-----------------------------------------------
!     ... Interpolate
!-----------------------------------------------
      do iz = 1, nz
        n1 = n
        y1 = ytmp(iz,:)
        x1 = xin

        call addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
        call addpnt(x1,y1,kdata,n1,          0.,0.)
        call addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
        call addpnt(x1,y1,kdata,n1,        1e38,0.)
        call inter2(nw,wl,yg1,n1,x1,y1,ierr)
        ytd(iz,:) = yg1(:)
 
        if (ierr .ne. 0) then
          write(*,*) ierr, jlabel(j)
          stop
        endif
      enddo

!----------------------------------------------
!     ...Quantum yields JPL06
!     ...This recommendation is only for >300nm
!        However, it is used at all wavelengths
!----------------------------------------------
      qy1 = 0.85
      qy2 = 0.15
      do iw = 1, nw-1
         do iz = 1, nz
            sq(j,iz,iw)   = qy1 * ytd(iz,iw)
            sq(j+1,iz,iw) = qy2 * ytd(iz,iw)
         enddo
      enddo
	
      j = j+1
      jlabel(j) = 'BrONO2 + hv -> BrO + NO2'
      
      end subroutine XSQY_BRONO2
