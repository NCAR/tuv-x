      SUBROUTINE la_srb(nz,z,tlev,nw,wl,o2col,vcol,scol,
     $                  o2xs1,dto2,o2xs,pathname)
!---------------------------------------------------------------------------!
!  PURPOSE:                                                                 !
!  Compute equivalent optical depths for O2 absorption, and O2 effective    !
!  absorption cross sections, parameterized in the Lyman-alpha and SR bands !
!---------------------------------------------------------------------------! 
!  PARAMETERS:                                                              !
!  NZ      - INTEGER, number of specified altitude levels in the working (I)!
!            grid                                                           !
!  Z       - REAL, specified altitude working grid (km)                  (I)!
!  NW      - INTEGER, number of specified intervals + 1 in working       (I)!
!            wavelength grid                                                !
!  WL      - REAL, vector of lxower limits of wavelength intervals in    (I)!
!            working wavelength grid                                        !
!  CZ      - REAL, number of air molecules per cm^2 at each specified    (I)!
!            altitude layer                                                 !
!  ZEN     - REAL, solar zenith angle                                    (I)!
!                                                                           !
!  O2XS1   - REAL, O2 cross section from rdo2xs                          (I)!
!                                                                           !
!  DTO2    - REAL, optical depth due to O2 absorption at each specified  (O)!
!            vertical layer at each specified wavelength                    !
!  O2XS    - REAL, molecular absorption cross section in SR bands at     (O)!
!            each specified altitude and wavelength.  Includes Herzberg     !
!            continuum.                                                     !
!---------------------------------------------------------------------------!
!  EDIT HISTORY:                                                            !
!  02/02  Major revision only over-write LA and SRB                         !
!  02/02  add Koppers and delete Korcarts                                   !
!  02/98  Included Lyman-alpha parameterization                             !
!  03/97  Fix dto2 problem at top level (nz)                                !
!  02/97  Changed offset for grid-end interpolation to relative number      !
!         (x * (1 +- deltax))                                               !
!  08/96  Modified for early exit, no redundant read of data and smaller    !
!         internal grid if possible;  internal grid uses user grid points   !
!         whenever possible                                                 !
!  07/96  Modified to work on internal grid and interpolate final values    !
!         onto the user-defined grid                                        !
!---------------------------------------------------------------------------!
      implicit none
      include 'params'

      integer  nz, nw, iz, iw
      real wl(kw), z(kz)
      real vcol (kz),  scol (kz)
      real o2col(kz),  o2xs1(kw)
      real dto2(kz,kw), o2xs(kz,kw)
      real secchi(kz)
      real tlev(kz)
      character*80 pathname

!----------------------------------------------------------------------
!     Lyman-alpha variables
!     O2 optical depth and equivalent cross section in the 
!        Lyman-alpha region
!----------------------------------------------------------------------
      integer ila, nla, kla
      parameter (kla = 2)
      real wlla(kla)
      real dto2la(kz, kla-1), o2xsla(kz, kla-1)
      save ila

!---------------------------------------------------------------------- 
!     Grid on which Koppers' parameterization is defined
!     O2 optical depth and equivalent cross section on Koppers' grid
!----------------------------------------------------------------------
      integer isrb, nsrb, ksrb
      parameter(ksrb = 18)
      real  wlsrb(ksrb)
      real  dto2k(kz, ksrb-1), o2xsk(kz, ksrb-1)
      save  isrb

      integer i

      logical call1
      data call1/.TRUE./
      save call1

!----------------------------------------------------------------------
!      Wavelengths for Lyman alpha and SRB parameterizations:
!----------------------------------------------------------------------
      data nla /1/
      data wlla/ 121.0, 122.0/

      data nsrb /17/
      data wlsrb/174.4, 177.0, 178.6, 180.2, 181.8, 183.5, 185.2, 186.9,
     $     188.7, 190.5, 192.3, 194.2, 196.1, 198.0, 200.0, 202.0, 
     $     204.1, 205.8/

!----------------------------------------------------------------------
!      initalize O2 cross sections 
!----------------------------------------------------------------------
      DO iz = 1, nz
         DO iw =1, nw - 1   
            o2xs(iz,iw) = o2xs1(iw)
         ENDDO  
      ENDDO

      IF(wl(1) .GT. wlsrb(nsrb)) RETURN

!----------------------------------------------------------------------
! On first call, check that the user wavelength grid, WL(IW), is compatible 
! with the wavelengths for the parameterizations of the Lyman-alpha and SRB.
! Also compute and save corresponding grid indices (ILA, ISRB)
!----------------------------------------------------------------------
      IF (call1) THEN

!     locate Lyman-alpha wavelengths on grid

         ila = 0
         DO iw = 1, nw
            IF(ABS(wl(iw) - wlla(1)) .LT. 10.*precis) THEN
               ila = iw
               GO TO 5
            ENDIF
         ENDDO
 5       CONTINUE
        
!     check 
         IF(ila .EQ. 0) STOP ' Lyman alpha grid mis-match - 1'
         DO i = 2, nla + 1
            IF(ABS(wl(ila + i - 1) - wlla(i)) .GT. 10.*precis) THEN
               WRITE(*,*) 'Lyman alpha grid mis-match - 2'
               STOP
            ENDIF
         ENDDO

!     locate Schumann-Runge wavelengths on grid
         isrb = 0
         DO iw = 1, nw
            IF(ABS(wl(iw) - wlsrb(1)) .LT. 10.*precis) THEN
               isrb = iw
               GO TO 6
            ENDIF
         ENDDO
 6       CONTINUE
         

!     check
         IF(isrb .EQ. 0) STOP ' SRB grid mis-match - 1'
         DO i = 2, nsrb + 1
            IF(ABS(wl(isrb + i - 1) - wlsrb(i)) .GT. 10.* precis) THEN
               WRITE(*,*) ' SRB grid mismatch - w'
               STOP
            ENDIF
         ENDDO

         IF (call1) call1 = .FALSE.
      ENDIF

!----------------------------------------------------------------------
! Effective secant of solar zenith angle.  
! Use 2.0 if no direct sun (value for isotropic radiation)
! For nz, use value at nz-1
!----------------------------------------------------------------------
      DO i = 1, nz - 1
         secchi(i) = scol(i)/vcol(i)
         IF(scol(i) .GT. largest/10.) secchi(i) = 2.
      ENDDO
      secchi(nz) = secchi(nz-1)

!---------------------------------------------------------------------
! Lyman-Alpha parameterization, output values of O2 optical depth
! and O2 effective (equivalent) cross section
!---------------------------------------------------------------------
      CALL lymana(nz,o2col,secchi,dto2la,o2xsla)

      DO iw = ila, ila + nla - 1
         DO iz = 1, nz
            dto2(iz,iw) = dto2la(iz, iw - ila + 1)
            o2xs(iz,iw) = o2xsla(iz, iw - ila + 1)
         ENDDO
      ENDDO

!----------------------------------------------------------------------
! Koppers' parameterization of the SR bands, output values of O2
! optical depth and O2 equivalent cross section 
!----------------------------------------------------------------------

      CALL schum(nz,o2col,tlev,secchi,dto2k,o2xsk,pathname)
      DO iw = isrb, isrb + nsrb - 1
         DO iz = 1, nz
            dto2(iz,iw) = dto2k(iz, iw - isrb + 1)
            o2xs(iz,iw) = o2xsk(iz, iw - isrb + 1)
         ENDDO
      ENDDO


      RETURN
      END





