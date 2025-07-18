      SUBROUTINE XSQY_ACETONE(nw,wl,wc,nz,tlev,airden,j,sq,jlabel,pn)
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Provide product (cross section) x (quantum yield) for CH3COCH3 photolysis!
!          CH3COCH3 + hv -> Products                                        !
!                                                                           !
!  Cross section:  Choice between                                           !
!                   (1) Calvert and Pitts                                   !
!                   (2) Martinez et al., 1991, alson in IUPAC 97            !
!                   (3) NOAA, 1998, unpublished as of 01/98                 !
!  Quantum yield:  Choice between                                           !
!                   (1) Gardiner et al, 1984                                !
!                   (2) IUPAC 97                                            !
!                   (3) McKeen et al., 1997                                 !
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

!-----------------------------------------------------------------------------!
!     ... input/output                                                        !
!-----------------------------------------------------------------------------!
      integer, intent(inout) :: j

!-----------------------------------------------------------------------------!
!     ... local                                                               !
!-----------------------------------------------------------------------------!
      integer  kdata
      parameter (kdata=150)

      integer  i, n, n1, n2, n3, iw, ierr, iz, idum
      real x1(kdata), x2(kdata), A(kdata), B(kdata), C(kdata)
      real y1(kdata), y2(kdata), y3(kdata)
      real xs(nz,kdata), sig(nz,kw)
      real yg(kw), yg1(kw), yg2(kw), yg3(kw)
      real tin(nz), AD(nz)
      real qytot(kw), qyCO(kw), qyCH3CO(kw)
      real AA0, a0, b0
      real AA1, a1, b1, t, qy
      real AA2, AA3, AA4, a2, b2, a3, b3, c3, a4, b4

      !!! TUV-x MOD - initializing qyCO, qyCH3CO !!!
      qyCO(:) = 0.0
      qyCH3CO(:) = 0.0
      !!! end TUV-x mod !!!

!---------------------------------------------
!     ... tin set to tlev
!---------------------------------------------
      tin(:) = tlev(:)
      AD (:) = airden(:)

!---------------------------------------------
!     ... CH3COCH3 photodissociation
!---------------------------------------------
      j = j + 1
      jlabel(j) = 'CH3COCH3 + hv -> CH3CO3 + CH3O2'

!---------------------------------------------
!     ... cross sections JPL06
!---------------------------------------------
      open(kin,file=TRIM(pn)//'XS_ACETONE_JPL06.txt',status='old')

      read(kin,*) idum, n
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n
        read(kin,*) x1(i), y1(i)
!        print*, x1(i), y1(i)
      enddo
      close(kin)

!---------------------------------------------
!     ... cross sections TD coeff JPL06
!---------------------------------------------
      open(kin,file=TRIM(pn)//'XS_ACETONE_TD_JPL06.txt',status='old')

      read(kin,*) idum, n1
      do i = 1, idum-2
        read(kin,*)
      enddo

      do i = 1, n1
        read(kin,*) x1(i), A(i), B(i), C(i)
         A(i) = A(i)*1e-3
         B(i) = B(i)*1e-5
         C(i) = C(i)*1e-8
!        print*, x1(i), y1(i), A(i), B(i), C(i)
      enddo
      close(kin)
!      stop
!---------------------------------------------
!     ... Derive XS at given temperature
!---------------------------------------------
    
      do iz = 1, nz

        do iw = 1, n1

          if ((tin(iz) .GE. 235.) .AND. (tin(iz) .LE. 298.)) Then
            xs(iz,iw) = y1(iw) *( 1 + (A(iw)*tin(iz)) +
     &                       (B(iw)*tin(iz)**2)  +
     &                       (C(iw)*tin(iz)**3) )
           
          endif

          if (tin(iz) .LT. 235.) then
            xs(iz,iw) = y1(iw) *( 1 + (A(iw)*235.) +
     &                       (B(iw)*(235.)**2)  +
     &                       (C(iw)*(235.)**3) )
   
          endif

          if (tin(iz) .GT. 298.) then
            xs(iz,iw) = y1(iw) *( 1 + (A(iw)*298.) +
     &                       (B(iw)*(298.)**2)  +
     &                       (C(iw)*(298.)**3) )
   
          endif

        enddo

        n     = n1
        x2(:) = x1(:)
        y2(:) = xs(iz,:)

!---------------------------------------------
!     ... Interpolate
!---------------------------------------------
        call addpnt(x2,y2,kdata,n,x2(1)*(1.-deltax),0.)
        call addpnt(x2,y2,kdata,n,          0.,0.)
        call addpnt(x2,y2,kdata,n,x2(n)*(1.+deltax),0.)
        call addpnt(x2,y2,kdata,n,        1e38,0.)
        call inter2(nw,wl,yg,   n,x2,     y2,ierr)
 
        sig(iz,:) = yg(:)

        if (ierr .ne. 0) then
          write(*,*) ierr, jlabel(j)
          stop
        endif

      enddo

!---------------------------------------------
!     ... Check Routine
!      iz = 10
!      print*, 'tin=', tin(iz)
!      do iw = 40, 80
!        print*, iw, wc(iw), sig(iz,iw)
!      enddo
!      stop
!---------------------------------------------
!---------------------------------------------
!     ... quantum yield JPL06
!---------------------------------------------
      DO iz = 1, nz

         T = min(tin(iz), 295.)
         T = max(T, 218.)
             
         DO iw = 1, nw-1

               IF ((wc(iw) .GE. 279.).AND.(wc(iw) .LT. 327.) ) THEN

                  a0 = 0.350* (T/295.)**(-1.28)
                  b0 = 0.068* (T/295.)**(-2.65)
                 AA0 = (a0 / (1-a0))* exp(b0*(wc(iw)-248.))
                 qyCO(iw) = 1. / (1. + AA0)
 !                print*, 'qyCO', qyCO(iw)

               ENDIF

               IF ((wc(iw) .GE. 279.).AND.(wc(iw) .LT. 302.)) THEN

                  a1 = 1.6e-19* (T/295.)**(-2.38) 
                  b1 = 0.55e-3* (T/295.)**(-3.19)
                 AA1 = a1* exp(-b1*((1e7/wc(iw)) - 33113.))
                 qyCH3CO(iw) = (1-qyCO(iw)) / (1 + AA1*AD(iz))

  !               print*, 'qyCO', qyCO(iw), 'qyCH3CO', qyCH3CO(iw)

               ELSEIF ((wc(iw) .GE. 302.).AND.(wc(iw) .LE. 327.5)) THEN

                  a2= 1.62e-17* (T/295.)**(-10.03)
                  b2= 1.79e-3 * (T/295.)**(-1.364)
                 AA2= a2* exp(-b2*((1e7/wc(iw))-30488.))

                  a3= 26.29*   (T/295.)**(-6.59)
                  b3= 5.72e-7* (T/295.)**(-2.93)
                  c3= 30006.*  (T/295.)**(-0.064)
                 AA3= a3* exp(-b3*((1e7/wc(iw))-c3)**2)

                  a4= 1.67e-15* (T/295.)**(-7.25)
                  b4= 2.08e-3*  (T/295.)**(-1.16)
                 AA4= a4* exp(-b4*((1e7/wc(iw)) - 30488.))

                 qyCH3CO(iw) = ((1 + AA4*AD(iz) + AA3) /
     &                         ((1 + AA2*AD(iz) + AA3)*
     &                          (1 + AA4*AD(iz))))*(1-qyCO(iw))
      
!                  print*, 'qyCH3CO', qyCH3CO(iw)          
   
               ELSEIF (wc(iw) .GT. 327.5) THEN
                  qytot(iw)  = 0.
               ENDIF

               qytot(iw) = qyCO(iw) + qyCH3CO(iw)
               
               if (wc(iw) .LT. 279.) then
                 qytot(iw) = 1.0
               endif

               qytot(iw) = max(0., qytot(iw))
               qytot(iw) = min(1., qytot(iw))


               sq(j,iz,iw) = sig(iz,iw)*qytot(iw)

         ENDDO
      ENDDO
!---------------------------------------------
!     ... Check Routine
!      iz = 10
!      print*, 'tin=', tin(iz)
!      do iw = 40, 80
!        print*, iw, wc(iw), sig(iz,iw), qytot(iw)
!        print*, iw, wc(iw), sq(j,iz,iw)
!      enddo
!      stop
!---------------------------------------------
!     c210417 there are issues with Actone being less than zero.
!      do iz = 1, nz-1
!         do iw = 1, nw-1
!            IF (sq(j,iz,iw) .LE. 0.0) THEN
!               sq(j,iz,iw) = 0.0
!            ENDIF
!         ENDDO
!      ENDDO
      
      end subroutine XSQY_ACETONE
