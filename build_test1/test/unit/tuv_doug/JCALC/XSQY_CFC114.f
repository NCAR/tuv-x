      subroutine XSQY_CFC114(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for cfc113 photolysis:  !
!          CFC114 + hv -> 2Cl                                                 !
!   cross section: from JPL06 recommendation                                  !
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
!   01/06/12  Doug Kinnison                                                   !
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
      parameter(kdata=300)
      integer i, iw, n, idum, nloop, n1
      integer ierr, iz, iwc, icnt
      real x1   (kdata),    y1   (kdata)
      real xin  (kdata),    yin  (kdata)
      real wctmp(kdata),    wcb  (kdata)
      real ytmp (nz,kdata), ycomb(nz,kdata)
      real yg1(kw),         tin  (nz)
      real ytd (nz,kw)
      real AA(5), BB(5), lp(5)
      real qy, ysave

      AA(1) = -160.50
      AA(2) =  2.4807
      AA(3) = -1.5202e-2
      AA(4) =  3.8412e-5
      AA(5) = -3.4373e-8

      BB(1) =  -1.5296
      BB(2) =   3.5248e-2
      BB(3) =  -2.9951e-4
      BB(4) =   1.1129e-6
      BB(5) =  -1.5259e-9

      lp(1) = 0.0
      lp(2) = 1.0
      lp(3) = 2.0
      lp(4) = 3.0
      lp(5) = 4.0

!---------------------------------------------------
!     ... tin set to tlev
!---------------------------------------------------
      tin(:)   = tlev(:)
!---------------------------------------------------
!     ... jlabel(j) = 'cfc114 -> 2cl'
!---------------------------------------------------
      j = j+1
      jlabel(j) = 'CFC114 + hv -> 2Cl'

!---------------------------------------------------
!    Derive temperature dependence 
!---------------------------------------------------
! Temperature dependence good between 
!     210-300K and 172 nm-220 nm
!---------------------------------------------------
      iwc      = 1
      ytmp(:,:)= 0.0 

      do iw = 1, nw-1
      
        IF ((wc(iw) .GE. 172.) .AND. (wc(iw) .LE.220.)) THEN
      
          do iz = 1, nz
            
            IF (tin(iz) .LT. 210.) THEN
              do nloop = 1, 5
                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                    +  AA(nloop)* (wc(iw)**lp(nloop))
     &                    + (210.0-273.0)*BB(nloop)*wc(iw)**lp(nloop)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF

            IF ((tin(iz) .GE. 210.).AND.(tin(iz) .LE. 300.)) THEN
              do nloop = 1,5

                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                    +  AA(nloop)* (wc(iw)**lp(nloop)) 
     &                    + (tin(iz)-273.0)*BB(nloop)*wc(iw)**lp(nloop)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF

            IF (tin(iz) .GT. 300.) THEN
              do nloop = 1, 5
                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                    +  AA(nloop)* (wc(iw)**lp(nloop)) 
     &                    + (300.0-273.0)*BB(nloop)*wc(iw)**lp(nloop)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF
          enddo
         iwc = iwc+ 1
         
        ENDIF
      
      enddo
!---------------------------------------------------
!     ... For wavelengths >220nm 
!---------------------------------------------------
      open(kin,file=TRIM(pn)//'XS_CFC114_JPL10.txt',status='old')

      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n
        read(kin,*) xin(i), yin(i)
      enddo
      close(kin)

!---------------------------------------------------
!     ... Combine cross sections
!---------------------------------------------------
      do iz = 1, nz
        icnt = 1

!     ... 172-220 nm
        do i = 1, iwc-1
          ycomb(iz,icnt) = 10**(ytmp(iz,i))
          wcb  (icnt)    = wctmp(i)
          icnt = icnt+1
        enddo
!     ... >20nm
        do i = 1, n
          IF (xin(i) .GT. 220.) THEN
            ycomb(iz,icnt) = yin(i)
            wcb  (icnt)    = xin(i)
            icnt = icnt+1
          ENDIF
        enddo
      enddo
!---------------------------------------------------
!     ... interpolate
!---------------------------------------------------
      do iz = 1, nz
        n1 = icnt-1
        y1 = ycomb(iz,:)
        x1 = wcb
!---------------------------------------------------    
!       iz = 1       
!       do iw = 1, icnt-1
!         print*, iw, wcb(iw), ycomb(iz,iw), tin(iz)
!       enddo
!       stop
!---------------------------------------------------  
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

!----------------------------------------------------     
!       iz = 1 
!       do iw = 10, 64
!         print*, iw, wc(iw), ytd(iz,iw), tin(iz)
!       enddo
!       stop
!----------------------------------------------------
!----------------------------------------------------
!     ...quantum yield assumed to be unity
!----------------------------------------------------
      qy = 1.
      do iw = 1, nw-1
        do iz = 1, nz
          sq(j,iz,iw) = qy * ytd(iz,iw)
        enddo
      enddo
      
      end subroutine XSQY_CFC114
