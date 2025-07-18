      subroutine XSQY_N2O5(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)

!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield):                        !
!           N2O5 + hv -> NO3 + NO2                                            !
!           N2O5 + hv -. NO3 + NO + O                                         !
!   cross section: JPL06                                                      !
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
!   01/17/08  Doug Kinnison                                                   !
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
      integer   kdata
      parameter (kdata=300)
      real x1   (kdata)
      real y1   (kdata)
      real Xin  (kdata)
      real Yin  (kdata)
      real wctmp(kdata)
      real wcb  (kdata)
      real ytmp (nz,kdata)
      real ycomb(nz,kdata)
      real ytd  (nz,kw)
      real yg   (kw)
      real yg1  (kw)
      real qy_O3p
      real ww (16), AA(16), BB(16)
      real tin(nz)
      real XS_harwood (nz,16)
      real xstest (nz,kw), qytest(nz,kw), qytest1(nz,kw)
 
      integer i, iw, n, idum, ierr
      integer iz, icnt, iwc, n1

      data ww / 260.0,   270.0,   280.0,    290.0,   300.0,  
     $          310.0,   320.0,   330.0,    340.0,   350.0,
     $          360.0,   370.0,   380.0,    390.0,   400.0,
     $          410.0 /
  
      data aa / -18.27,  -18.42,  -18.59,   -18.72,  -18.84,
     $          -18.90,  -18.93,  -18.87,   -18.77,  -18.71,
     $          -18.31,  -18.14,  -18.01,   -18.42,  -18.59,
     $          -18.13 /
      
      data bb / -0.091,  -0.104,  -0.112,   -0.135,  -0.170, 
     $          -0.226,  -0.294,  -0.388,   -0.492,  -0.583, 
     $          -0.770,  -0.885,  -0.992,   -0.949,  -0.966,
     $          -1.160 /

!----------------------------------------------------
!     ... tin set to tlev
!----------------------------------------------------
      tin(:) = tlev(:)

!----------------------------------------------------
!     ... Calculate the T-dep XS (233-295K)
!         and 260-410nm
!----------------------------------------------------
      do iw = 1, 16
         do iz = 1, nz
           IF (tin(iz) .LT. 200.0) THEN
             XS_harwood(iz,iw) = 10**(aa(iw) + (1000.*bb(iw)/200.0))
           ENDIF
           IF ((tin(iz) .GE. 200.0) .AND. (tin(iz) .LE. 295.)) THEN
             XS_harwood(iz,iw) = 10**(aa(iw) + (1000.*bb(iw)/tin(iz)))
           ENDIF
           IF (tin(iz) .GT. 295.0) THEN
             XS_harwood(iz,iw) = 10**(aa(iw) + (1000.*bb(iw)/295.0))
           ENDIF
         enddo
      enddo    
!---------------------------------------
!      iz = 1
!      print*, 'tin=', tin(iz)
!      do iw = 1, 16 
!        print*, ww(iw), XS_harwood(iz,iw)
!      enddo    
!----------------------------------------
!----------------------------------------
!     ... cross sections from JPL06 
!----------------------------------------
      open(kin,file=TRIM(pn)//'XS_N2O5_JPL06.txt',status='old')

      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n
        read(kin,*) Xin(i), Yin(i)
      enddo
      close(kin)
!---------------------------------------
!      do iw = 1, n 
!        print*, iw, xin(iw), yin(iw)
!      enddo
!      stop
!---------------------------------------
!     ... Combine cross sections
      do iz = 1, nz
        icnt = 1

!     ... < 260 nm
        do i = 1, n
          IF (xin(i) .LT. 260.) THEN
            ycomb(iz,icnt) = yin(i)
            wcb  (icnt)    = xin(i)
            icnt = icnt + 1
          ENDIF
        enddo
!     ... 260-410 nm
        do i = 1, 16
          ycomb(iz,icnt) = (xs_harwood(iz,i))
          wcb  (icnt)    =  ww(i)
          icnt = icnt+1
        enddo
!     ... >410 nm
        do i = 1, n
          IF (xin(i) .GT. 410.) THEN
            ycomb(iz,icnt) = yin(i)
            wcb  (icnt)    = xin(i)
            icnt = icnt+1
          ENDIF
        enddo
      enddo

!     ... Interpolate to TUV grid 
      do iz = 1, nz
         n1 = icnt-1
         y1 = ycomb(iz,:)
         x1 = wcb
!----------------------------------------------------------         
!       do iw = 1, icnt-1
!       print*, iw, wcb(iw), ycomb(iz,iw), tin(iz)
!       enddo       
!       stop
!---------------------------------------------------------- 
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
!-------------------------------------------------------     
!       iz = 1 
!       do iw = 33, 95
!         print*, iw, wc(iw), ytd(iz,iw), tin(iz)
!       enddo
!       stop
!-------------------------------------------------------
!-------------------------------------------------------
!     ... quantum yield (JPL06)
!-------------------------------------------------------
!     Branch #1 label
      j = j+1
      jlabel(j) = 'N2O5 + hv -> NO3 + NO2'

      do iw = 1, nw-1

        if (wc(iw) .GE. 300.0) THEN 
           qy_O3p = 0.0
          do iz = 1, nz
            sq(j,  iz,iw) = 1.0 * ytd(iz,iw)
            sq(j+1,iz,iw) = qy_O3p
          enddo
        endif

        if (wc(iw) .LT. 300.0) THEN
          qy_O3p = min( 1., 3.832441 - 0.012809638 * wc(iw) )
          qy_O3p = max( 0., qy_O3p )
          do iz = 1, nz
            sq(j,  iz,iw) = (1.0-qy_O3p)*ytd(iz,iw)
            sq(j+1,iz,iw) =      qy_O3p *ytd(iz,iw)

          enddo
        endif   
       
       enddo     

!-------------------------------------------------------
!     Branch #2 label
      j = j+1
      jlabel(j) = 'N2O5 + hv -> NO3 + NO + O'

!-------------------------------------------------------
!     ... Check routine 
!      iz=1
!      do iw = 30,100
!        print*, wc(iw), sq(j-1,iz,iw), 
!     $                  sq(j,iz,iw),sq(j-1,iz,iw)+sq(j,iz,iw) 
!      enddo
!      stop
!-------------------------------------------------------
!      print*, '* N2O5 TD allowed down to 200K'
      end subroutine  XSQY_N2O5
