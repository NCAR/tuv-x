      subroutine XSQY_HO2NO2(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product of (cross section) x (quantum yield) for hno4 photolysis  !
!     1)   HO2NO2 + hv -> HO2 + NO2                                           !
!     2)   HO2NO2 + hv -> OH +  NO3                                           !
!   cross sections and QY from JPL06                                          !
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
!   06/01  modified by doug kinnison                                          !
!   01/08  modified by Doug Kinnison                                          !
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
      integer i, iw, iz, n, n1, idum, ierr, icnt
      real x1   (kdata), x2(kdata), wcb(kdata)
      real y1   (kdata), aa(kdata), bb (kdata)
      real ytmp (nz,kdata), ycomb(nz,kdata)
      real ytd  (nz,kw), yg(kw)
      real Q(nz), tin(nz), t

!----------------------------------------------
!     ... tin set to tlev
!----------------------------------------------
      tin(:) = tlev(:) 

!----------------------------------------------
!     ... jlabel(j) = 'HO2NO2 -> HO2 + NO2
!         jlabel(j) = 'HO2NO2 -> OH + NO3
!----------------------------------------------
      j = j + 1
      jlabel(j) = 'HO2NO2 + hv -> OH + NO3'

!----------------------------------------------
!     ...ho2no2 cross sections plus T-dep. 
!        (Burkholder et al., 2002.)
!----------------------------------------------
      open(kin,file=TRIM(pn)//'XS_HO2NO2_JPL06.txt',status='old')

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
      read(kin,*) idum, n1
      do i = 1, n1
        read(kin,*) x2(i), aa(i), bb(i)
      enddo
      close(kin)

!----------------------------------------------
!     ...Derive T-dep Burkholder et al., 2002.)
!----------------------------------------------
      do iz = 1, nz
        do iw = 1, n1
          t           = MAX(280.,MIN(tin(iz),350.))
          Q(iz)       = 1 + exp(-988./(0.69*t))
          ytmp(iz,iw) = ( aa(iw)/Q(iz) + bb(iw)*(1-1/Q(iz)))*1e-20
        enddo
      enddo
!----------------------------------------------
!     ... Check routine
!      iz = 1
!      do iw = 1, n1
!        print*, iw, x2(iw), ytmp(iz,iw)
!      enddo
!      stop
!----------------------------------------------
!     ... Combine cross sections
      do iz = 1, nz
        icnt = 1

!     ... < 280 nm
!     ... x1(iw) goes from 190-350nm
        do iw = 1, n
          IF (x1(iw) .LT. 280.) THEN
            ycomb(iz,icnt) = y1(iw)
            wcb  (icnt)    = x1(iw)
            icnt = icnt + 1
          ENDIF     
        enddo 
!     ... 280-350 nm
        do iw = 1, n1
            ycomb(iz,icnt) = ytmp(iz,iw)
            wcb  (icnt)    = x2  (iw)
            icnt = icnt+1
        enddo
       enddo

!...   Test No TD   
!       do iz = 1, nz
!        icnt = 1
!        do iw = 1, n
!         ycomb(iz,icnt) = y1(iw)
!         wcb  (icnt)    = x1(iw)
!         icnt = icnt + 1
!        enddo
!       enddo
!----------------------------------------------
!     ... Check routine
!      iz = 1
!      print*,"tin=", tin(iz)
!      do iw = 1, icnt-1
!        print*, iw, wcb(iw), ycomb(iz,iw)
!      enddo
!      stop
!----------------------------------------------
!     ... Interpolate Combine cross sections
      do iz = 1, nz
        n  = icnt-1
        y1 = ycomb(iz,:)
        x1 = wcb

        call addpnt(x1,y1,kdata,n,x1(1)*(1.-deltax),0.)
        call addpnt(x1,y1,kdata,n,               0.,0.)
        call addpnt(x1,y1,kdata,n,x1(n)*(1.+deltax),0.)
        call addpnt(x1,y1,kdata,n,           1.e+38,0.)
        call inter2(nw,wl,yg,n,x1,y1,ierr)
        ytd(iz,:) = yg(:)
        
        if (ierr .ne. 0) then
           write(*,*) ierr, jlabel(j)
           stop
        endif
      enddo

!-------------------------------------------------    
!       iz = 1 
!       do iw = 24, 80
!         print*, iw, wc(iw), ytd(iz,iw), tin(iz)
!       enddo
!       stop
!-------------------------------------------------
      do iw = 1, nw - 1
         IF (wc(iw) .LT. 200.0) THEN
             do iz = 1, nz
               sq(j,  iz,iw) = 0.30 * ytd(iz,iw)
	       sq(j+1,iz,iw) = 0.70 * ytd(iz,iw)
             enddo
           ENDIF
           IF (wc(iw) .GE. 200.0) THEN
             do iz = 1, nz
               sq(j,  iz,iw) = 0.20 * ytd(iz,iw)
	       sq(j+1,iz,iw) = 0.80 * ytd(iz,iw)
             enddo
           ENDIF
        enddo

!-------------------------------------------------- 
!       iz = 1 
!       do iw = 24, 80
!         print*, wc(iw), sq(j,iz,iw), sq(j+1,iz,iw) 
!         print*, sq(j,iz,iw)+sq(j+1,iz,iw) 
!       enddo
!       stop
!-------------------------------------------------
      j = j + 1
      jlabel(j) = 'HO2NO2 + hv -> HO2 + NO2'

      end subroutine XSQY_HO2NO2
