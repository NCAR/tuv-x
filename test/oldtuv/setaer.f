      MODULE SET_AEROSOL_OPTICAL_PROPERTIES

      IMPLICIT NONE

      private
      public :: setaer

      contains

      SUBROUTINE setaer(
     $     ipbl, zpbl, aod330,
     $     tau550, ssaaer, alpha,
     $     z, wc, dtaer, omaer, gaer)
*-----------------------------------------------------------------------------*
*=  PURPOSE:                                                                 =*
*=  Set up an altitude profile of aerosols, and corresponding absorption     =*
*=  optical depths, single scattering albedo, and asymmetry factor.          =*
*=  Single scattering albedo and asymmetry factor can be selected for each   =*
*=  input aerosol layer (do not have to correspond to working altitude       =*
*=  grid).  See loop 27.                                                     =*
*-----------------------------------------------------------------------------*
*=  PARAMETERS:                                                              =*
*=  Z       - REAL, specified altitude working grid (km)                  (I)=*
*=  WL     - REAL, vector of lower limits of wavelength intervals in      (I)=*
*=           working wavelength grid                                         =*
*=  DTAER   - REAL, optical depth due to absorption by aerosols at each   (O)=*
*=            altitude and wavelength                                        =*
*=  OMAER   - REAL, single scattering albedo due to aerosols at each      (O)=*
*=            defined altitude and wavelength                                =*
*=  GAER    - REAL, aerosol asymmetry factor at each defined altitude and (O)=*
*=            wavelength                                                     =*
*-----------------------------------------------------------------------------*

      use debug,      only : diagout
      use tuv_params, only : kout, nzero, pzero

      INTEGER, parameter :: kdata = 51

* input:

      INTEGER, intent(in) :: ipbl
      REAL, intent(in) :: zpbl
      REAL, intent(in) :: aod330
      REAL, intent(in) :: tau550
      REAL, intent(in) :: ssaaer, alpha
      REAL, intent(in) :: wc(:)
      REAL, intent(in) :: z(:)

* output: (on converted grid)

      REAL, intent(out) :: dtaer(:,:), omaer(:,:), gaer(:,:)

* local:

      INTEGER :: nz
      INTEGER :: nbins
      INTEGER :: i, iw, nd
      REAL    :: colold
      REAL    :: wscale
*     REAL    :: zd(kdata), aer(kdata)
      REAL, allocatable :: zd(:), aer(:)
      REAL    :: cd(kdata-1), omd(kdata-1), gd(kdata-1)
      REAL    :: womd(kdata-1), wgd(kdata-1)

      REAL    :: cz(size(z)-1)
      REAL    :: omz(size(z)-1)
      REAL    :: gz(size(z)-1)
      REAL    :: aodw(size(wc)), ssaw(size(wc))
      REAL    :: fract(ipbl)

* Aerosol data from Elterman (1968)
* These are vertical optical depths per km, in 1 km
* intervals from 0 km to 50 km, at 340 nm.
* This is one option.  User can specify different data set.

      aer = (/
     1     2.40E-01,1.06E-01,4.56E-02,1.91E-02,1.01E-02,7.63E-03,
     2     5.38E-03,5.00E-03,5.15E-03,4.94E-03,4.82E-03,4.51E-03,
     3     4.74E-03,4.37E-03,4.28E-03,4.03E-03,3.83E-03,3.78E-03,
     4     3.88E-03,3.08E-03,2.26E-03,1.64E-03,1.23E-03,9.45E-04,
     5     7.49E-04,6.30E-04,5.50E-04,4.21E-04,3.22E-04,2.48E-04,
     6     1.90E-04,1.45E-04,1.11E-04,8.51E-05,6.52E-05,5.00E-05,
     7     3.83E-05,2.93E-05,2.25E-05,1.72E-05,1.32E-05,1.01E-05,
     8     7.72E-06,5.91E-06,4.53E-06,3.46E-06,2.66E-06,2.04E-06,
     9     1.56E-06,1.19E-06,9.14E-07 /)
*_______________________________________________________________________

      nz    = size(z)
      nbins = size(wc)

* Altitudes corresponding to Elterman profile, from bottom to top:

      WRITE(kout,*)'aerosols:  Elterman (1968) continental profile'
      nd = 51
      zd = (/ (REAL(i-1),i=1,kdata) /)

* assume these are point values (at each level), so find column
* increments

      omd = ssaaer
      gd  = .61
      DO i = 1, nd - 1
         cd(i) = .5 * (aer(i+1) + aer(i))
      ENDDO

      call diagout( 'rawOD.old',aer )
      call diagout( 'inpaerOD.old',cd )
      write(*,*) 'setaer: hardwired OD'
      write(*,'(1p10g15.7)') aer
      write(*,'(1p10g15.7)') cd

*********** end data input.

* Compute integrals and averages over grid layers:
* for g and omega, use averages weighted by optical depth

      DO i = 1, nd-1
         womd(i) = omd(i) * cd(i)
         wgd(i)  = gd(i) * cd(i)
      ENDDO

      CALL inter3(nz,z,cz, nd,zd,cd, 1)
      CALL inter3(nz,z,omz, nd, zd,womd, 1)
      CALL inter3(nz,z,gz , nd, zd,wgd, 1)

      call diagout( 'cz.aer.old',cz )
      call diagout( 'omz.aer.old',omz )
      call diagout( 'gz.aer.old',gz )

      WHERE( cz(:) > 0. )
         omz(:) = omz(:)/cz(:)
         gz(:)  = gz(:) /cz(:)
      ELSEWHERE
         omz(:) = 1.
         gz(:) = 0.
      ENDWHERE

* old column at 340 nm
*  (minimum value is pzero = 10./largest)

      colold = MAX(sum(cz(1:nz-1)),pzero)

* scale with new column tau at 550 nm

      IF(tau550 > nzero) THEN
         DO i = 1, nz-1
            cz(i) = cz(i) * (tau550/colold) * (550./340.)**alpha 
         ENDDO
      ENDIF

* assign at all wavelengths
* (can move wavelength loop outside if want to vary with wavelength)

      DO iw = 1, nbins
* Elterman's data are for 340 nm, so assume optical depth scales 
* inversely with first power of wavelength.
         wscale = (340./wc(iw))**alpha
* optical depths:
         dtaer(:,iw) = cz(:)  * wscale
         omaer(:,iw) = omz(:)
         gaer(:,iw)  = gz(:)
      ENDDO

*! overwrite for pbl:

      IF(ipbl > 0) THEN	
         write (*,*) 'pbl aerosols, aod330 = ', aod330
* create wavelength-dependent optical depth and single scattering albedo:
         DO iw = 1, nbins
            aodw(iw) = aod330*(wc(iw)/330.)**(-1.0)
            IF(wc(iw) < 400.) THEN
               ssaw(iw) = 0.6
            ELSE
               ssaw(iw) = 0.9
            ENDIF
         ENDDO

* divide aod among pbl layers, overwrite Elterman profile in pbl

         DO i = 1, ipbl
            fract(i) = (z(i+1) - z(i))/zpbl
         ENDDO

         DO iw = 1, nbins
            dtaer(1:ipbl,iw) = aodw(iw) * fract(1:ipbl)
            omaer(1:ipbl,iw) = ssaw(iw)
         ENDDO
      ENDIF

      END SUBROUTINE setaer

      END MODULE SET_AEROSOL_OPTICAL_PROPERTIES
