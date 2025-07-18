      subroutine XSQY_CHBR3(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for chbr3 photolysis:   !
!          CHBr3 + hv -> 3Br                                                  !
!   cross section: from Papanastasiou et al, ACP, 2014                        !
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
!   07/02/14  Doug Kinnison                                                   !
!   09/17/14  added <260
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
      integer i, iw, n, idum, nloopAA, nloopBB, n1
      integer ierr, iz, iwc, icnt

      real x1   (kdata),   y1   (kdata)
      real xin  (kdata),   yin  (kdata)
      real wctmp(kdata),   wcb  (kdata)
      real ytmp (nz,kdata),ycomb(nz,kdata)
      real yg1  (kw),      ytd  (nz,kw)
      real qy
      real AA(6), BB(5), lp(6)
      real tin(nz)

      AA(1) =  -32.6067
      AA(2) =   0.10308
      AA(3) =   6.39e-5
      AA(4) =  -7.7392e-7
      AA(5) =  -2.2513e-9
      AA(6) =   6.1376e-12

      BB(1) =   0.1582
      BB(2) =  -0.0014758
      BB(3) =   3.8058e-6
      BB(4) =   9.187e-10
      BB(5) =  -1.0772e-11
 
      lp(1) = 0.0
      lp(2) = 1.0
      lp(3) = 2.0
      lp(4) = 3.0
      lp(5) = 4.0
      lp(6) = 5.0
 
!----------------------------------------------
!     ... set tin to tlev
!----------------------------------------------
      tin(:)   = tlev(:)
!----------------------------------------------
!     ... jlabel(j) = 'CHBr3 + hv -> 3Br'
!----------------------------------------------
      j = j+1
      jlabel(j) = 'CHBr3 + hv -> 3Br'

!----------------------------------------------
!    Derive temperature dependence 
!----------------------------------------------
!    Temperature dependence good between 
!        260-330K and 260 nm to 345 nm
!        99% of the loss in this region
!----------------------------------------------
      iwc      = 1
      ytmp(:,:)= 0.0 

      do iw = 1, nw-1
      

!    Extrapolate to 357.5nm with TUV grid
        IF ((wc(iw) .GE. 260.) .AND. (wc(iw) .LE. 362.)) THEN      
          do iz = 1, nz
            
            IF (tin(iz) .LT. 260.) THEN
              do nloopAA = 1, 6
                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                  +  AA(nloopAA)* (wc(iw)**lp(nloopAA))                 
              enddo
              do nloopBB = 1, 5
                ytmp(iz,iwc) = ytmp(iz,iwc) + 
     &            (296.0-260.0)*BB(nloopBB)*wc(iw)**lp(nloopBB)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF

            IF ((tin(iz) .GE. 260.).AND.(tin(iz) .LE. 330.)) THEN
               do nloopAA = 1, 6
                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                  +  AA(nloopAA)* (wc(iw)**lp(nloopAA))                 
              enddo
              do nloopBB = 1, 5
                ytmp(iz,iwc) = ytmp(iz,iwc) + 
     &            (296.0-tin(iz))*BB(nloopBB)*wc(iw)**lp(nloopBB)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF

            IF (tin(iz) .GT. 330.) THEN
              do nloopAA = 1, 6
                ytmp(iz,iwc) = ytmp(iz,iwc)
     &                  +  AA(nloopAA)* (wc(iw)**lp(nloopAA))                 
              enddo
              do nloopBB = 1, 5
                ytmp(iz,iwc) = ytmp(iz,iwc) + 
     &            (296.0-tin(iz))*BB(nloopBB)*wc(iw)**lp(nloopBB)
              enddo
              wctmp(iwc) = wc(iw)
            ENDIF

          enddo
         iwc = iwc+ 1
         
        ENDIF
      
      enddo

!----------------------------------------------
!     ... For wavelengths >310 nm and <240 nm
!----------------------------------------------
      open(kin,file=TRIM(pn)//'XS_CHBR3_JPL10.txt',status='old')

      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n
        read(kin,*) xin(i), yin(i)
      enddo
      close(kin)

!----------------------------------------------
!     ... Combine cross sections
!----------------------------------------------
      do iz = 1, nz
        icnt = 1

!     ... < 260nm
        do i = 1, n
          IF (xin(i) .LT. 260.) THEN
            ycomb(iz,icnt) = yin(i)
            wcb  (icnt)    = xin(i)
            icnt = icnt + 1
          ENDIF
        enddo
!     ... 260-362 nm
        do i = 1, iwc-1
          ycomb(iz,icnt) = 10**(ytmp(iz,i))
          wcb  (icnt)    = wctmp(i)
          icnt = icnt+1
        enddo

      enddo
!----------------------------------------------
!     ... interpolate
!----------------------------------------------
      do iz = 1, nz
        n1 = icnt-1
        y1 = ycomb(iz,:)
        x1 = wcb
!----------------------------------------------     
!         
!       do iw = 1, icnt-1
!         print*, iw, wcb(iw), ycomb(iz,iw), tin(iz)
!       enddo
!       stop
!----------------------------------------------    
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
!       iz = 1 
!       do iw = 1, nw-1
!         print*, iw, wc(iw), ytd(iz,iw), tin(iz)
!       enddo
!       stop
!----------------------------------------------     
!----------------------------------------------
!     ...quantum yield assumed to be unity
!----------------------------------------------
      qy = 1.

      do iw = 1, nw-1
        do iz = 1, nz
          sq(j,iz,iw) = qy * ytd(iz,iw)
        enddo
      enddo
      
      end subroutine XSQY_CHBR3
