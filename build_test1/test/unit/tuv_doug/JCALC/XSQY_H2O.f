      subroutine XSQY_H2O(nw,wl,wc,nz,tlev,airlev,j,sq,jlabel,pn)
!-----------------------------------------------------------------------------!
!   purpose:                                                                  !
!   provide product (cross section) x (quantum yield) for hcl photolysis:     !
!           H2O + hv -> products                                              !
!   cross section: taken from three sources                                   !
!     1) JPL06 (jpl97-4), 175.5 - 189.3                                       !
!     2) Cantrell et al., grl, 24, 17, 2195-2198, 1997,  183.0 - 193.0 nm     !
!     3) Yoshino et al.,  chemical physics, 211 (1996) 387-391, 120.38-188.03 !
!                                                                             !
!   quantum yield: is unity between 175.5 and 189.3                           !
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
!   06/11/01   original, dek addition                                         !
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
      parameter(kdata=7000)
      integer i, iw, n1, n2, n3, n4, idum, ierr, iz
      real x1(kdata), y1(kdata)
      real x2(kdata), y2(kdata)
      real x3(kdata), y3(kdata)
      real x4(kdata), y4(kdata)
      real yg(kw), yg1(kw), yg2(kw), yg3(kw), yg4(kw)
      real qy

!----------------------------------------------
!     ... jlabel(j) = 'h2o -> prod'
!----------------------------------------------
      j = j+1
      jlabel(j) = 'H2O + hv -> H + OH'
      j = j+1
      jlabel(j) = 'H2O + hv -> H2 + O(1D)'
	j = j+1
      jlabel(j) = 'H2O + hv -> 2H + O(3P)'

!----------------------------------------------------
!     ... cross sections from JPL06 recommendation
!----------------------------------------------------
      open(kin,
     &   file=TRIM(pn)//'XS_H2O_JPL06.txt',status='old')
      read(kin,*) idum, n1
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n1
        read(kin,*) x1(i), y1(i)
        y1(i) = y1(i) * 1e-20
      enddo
      close(kin)

      call addpnt(x1,y1,kdata,n1,x1(1)*(1.-deltax),0.)
      call addpnt(x1,y1,kdata,n1,          0.,0.)
      call addpnt(x1,y1,kdata,n1,x1(n1)*(1.+deltax),0.)
      call addpnt(x1,y1,kdata,n1,        1e38,0.)

      call inter2(nw,wl,yg1,n1,x1,y1,ierr)
 
      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif

!----------------------------------------------------
!     ...cross sections from Cantrell et al., 1997
!----------------------------------------------------
      open(kin,
     &  file=TRIM(pn)//'XS_H2O_cantrell_1996.txt',status='old')
      read(kin,*) idum, n2
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n2
        read(kin,*) x2(i), y2(i)
        y2(i) = y2(i) * 1e-20
      enddo
      close(kin)

      call addpnt(x2,y2,kdata,n2,x2(1)*(1.-deltax),0.)
      call addpnt(x2,y2,kdata,n2,          0.,0.)
      call addpnt(x2,y2,kdata,n2,x2(n2)*(1.+deltax),0.)
      call addpnt(x2,y2,kdata,n2,        1e38,0.)

      call inter2(nw,wl,yg2,n2,x2,y2,ierr)
 
      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif
!----------------------------------------------------
!     ... cross sections from Yoshino et al., 1996
!----------------------------------------------------
      open(kin,
     &  file=TRIM(pn)//'XS_H2O_yoshino_1996.txt',status='old')
      read(kin,*) idum, n3
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n3
        read(kin,*) x3(i), y3(i)
      enddo
      close(kin)

      call addpnt(x3,y3,kdata,n3,x3(1)*(1.-deltax),0.)
      call addpnt(x3,y3,kdata,n3,          0.,0.)
      call addpnt(x3,y3,kdata,n3,x3(n3)*(1.+deltax),0.)
      call addpnt(x3,y3,kdata,n3,        1e38,0.)
      
      call inter2(nw,wl,yg3,n3,x3,y3,ierr)
 
      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif
!----------------------------------------------------
!     ... cross sections from Ranjan et al. 2020
!         189.059nm to 2019.316nm
!----------------------------------------------------
      open(kin,
     &  file=TRIM(pn)//'XS_H2O_Ranjan.txt',status='old')
      read(kin,*) idum, n4
      do i = 1, idum-2
        read(kin,*)
      enddo
      do i = 1, n4
        read(kin,*) x4(i), y4(i)
      enddo
      close(kin)

      call addpnt(x4,y4,kdata,n4,x4(1)*(1.-deltax),0.)
      call addpnt(x4,y4,kdata,n4,          0.,0.)
      call addpnt(x4,y4,kdata,n4,x4(n4)*(1.+deltax),0.)
      call addpnt(x4,y4,kdata,n4,        1e38,0.)
      
      call inter2(nw,wl,yg4,n4,x4,y4,ierr)
 
      if (ierr .ne. 0) then
        write(*,*) ierr, jlabel(j)
        stop
      endif
!--------------------------------------------------------
!     ...combine data sets (i.e., Yoshino et al., 1996 
!        and Cantrell et al., 1997)
!--------------------------------------------------------    
      do i = 1, nw-1
        if (wc(i) .lt. 183.0) then
	    yg(i) = yg3(i)
	  elseif (wc(i) .le. 191.0) then
             yg(i) = yg2(i)
          elseif (wc(i) .le. 220.0) then
             yg(i) = yg4(i)
	  else
	    yg(i) = 0.
	  endif
       enddo

!-------------------------------------------------------
!     ... Check routine (no temperature dependence)
!      do iw = 1, 45
!        print*, iw, wc(iw), yg(iw)
!      enddo
!      stop
!-------------------------------------------------------

!------------------------------------------------------
!     ... quantum yield assumed to be unity (jpl97-4)
!------------------------------------------------------
!     ... 105 to 145 nm
!         (JPL 1997 which references Stief, L.J., W.A. 
!         Payne, and R. B. Klemm, A flash
!         photolysis-resonance fluoresence study of the 
!         formation of O(1D) in the photolysis of water 
!         and the reaction of O(1D) with H2, Ar, and He, 
!         J. Chem. Phys., 62, 4000, 1975.)

      do iw = 1, nw-1

        if (wc(iw) .le. 145.0) then

	     do iz = 1, nz
             sq(j-2,iz,iw) = yg(iw) * 0.890
	       sq(j-1,iz,iw) = yg(iw) * 0.110
	       sq(j,  iz,iw) = yg(iw) * 0.0
           end do

	  end if

!     ... > 145nm
!         JPL97
        if (wc(iw) .gt. 145.0) then

	     do iz = 1, nz
               sq(j-2,iz,iw) = yg(iw) * 1.0
	         sq(j-1,iz,iw) = yg(iw) * 0.0
	         sq(j,  iz,iw) = yg(iw) * 0.0
           end do

	  end if

	end do	! end wavelength loop

!     ... Overwrite Lyamn Alpha
!         Slanger, T.G., and G. Black, Photodissociative 
!         channels at 1216A for H2O, NH3 and CH4,
!         J. Chem. Phys., 77, 2432, 1982.)
!         **** Caution, Lyman Alpha is always wavelength postion #2 ****
      iw = 2
      do iz = 1, nz
        sq(j-2,iz,iw)   = yg(iw) * 0.780
	  sq(j-1,iz,iw) = yg(iw) * 0.100
	  sq(j,  iz,iw) = yg(iw) * 0.120
      end do

      end subroutine XSQY_H2O     
