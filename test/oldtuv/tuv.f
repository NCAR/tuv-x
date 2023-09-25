      PROGRAM tuv
*-----------------------------------------------------------------------------*
*=    Tropospheric Ultraviolet-Visible (TUV) radiation model                 =*
*=    Version 5.4                                                            =*
*=    November 2018                                                          =*
*-----------------------------------------------------------------------------*
*= Developed by Sasha Madronich with important contributions from:           =*
*= Chris Fischer, Siri Flocke, Julia Lee-Taylor, Bernhard Meyer,             =*
*= Irina Petropavlovskikh,  Xuexi Tie, and Jun Zen.                          =*
*= Special thanks to Knut Stamnes and co-workers for the development of the  =*
*= Discrete Ordinates code, and to Warren Wiscombe and co-workers for the    =*
*= development of the solar zenith angle subroutine. Citations for the many  =*
*= data bases (e.g. extraterrestrial irradiances, molecular spectra) may be  =*
*= found in the data files headers and/or in the subroutines that read them. =*
*=              To contact the author, write to:                             =*
*= Sasha Madronich, NCAR/ACD, P.O.Box 3000, Boulder, CO, 80307-3000, USA  or =*
*= send email to:  sasha@ucar.edu  or tuv@acd.ucar.edu                       =*
*-----------------------------------------------------------------------------*
*= This program is free software; you can redistribute it and/or modify      =*
*= it under the terms of the GNU General Public License as published by the  =*
*= Free Software Foundation;  either version 2 of the license, or (at your   =*
*= option) any later version.                                                =*
*= The TUV package is distributed in the hope that it will be useful, but    =*
*= WITHOUT ANY WARRANTY;  without even the implied warranty of MERCHANTIBI-  =*
*= LITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public     =*
*= License for more details.                                                 =*
*= To obtain a copy of the GNU General Public License, write to:             =*
*= Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.   =*
*-----------------------------------------------------------------------------*
*= Copyright (C) 1994-2018 by the University Corporation for Atmospheric     =*
*= Research, extending to all called subroutines, functions, and data unless =*
*= another source is specified.                                              =*
*-----------------------------------------------------------------------------*

      use tuv_params, only : kz, kw, kj, ks, kout, kt, hc, kdom, nzero
      use rtlink_mod,  only : rtlink
      use rayleigh, only : odrl
      use CLOUD, only  : setcld
      use ALBEDO_MOD, only : setalb
      use COLUMN_O3, only   : vpo3
      use COLUMN_TEMP, only : vptmp
      use COLUMN_ATM, only  : vpair
      use O3_OPTICAL_DEPTH, only  : odo3
      use SET_O2_OPTICAL_TAU, only  : seto2
      use SET_SO2_OPTICAL_DTAU, only : setso2
      use SET_NO2_OPTICAL_DTAU, only : setno2
      use SET_AEROSOL_OPTICAL_PROPERTIES, only : setaer
      use SET_SNOW_OPTICAL_PROPERTIES, only : setsnw
      use ARCHIVE, only : saver1, saver2, outpt1
      use LA_SRB_MOD, only : init_la_srb, la_srb, sjo2
      use GRIDS, only : gridz, gridw, gridt
      use ETFL, only : rdetfl
      use RDXS, only : rdo2xs, rdo3xs, rdno2xs, rdso2xs
      use SPHERICAL_GEOM, only : sphers, airmas
      use MO_SWCHEM, only : swchem
      use micm_photo_kinetics, only : photo_kinetics_t
      use micm_spectral_wght_warehouse, only : spectral_wght_warehouse_t
      use micm_environment, only : environment_t
      use musica_constants, only : musica_dk, musica_rk
      use musica_config,    only : config_t
      use musica_string,    only : string_t
      use debug,            only : diagout

      IMPLICIT NONE

      real, parameter :: rZERO = 0.0

* Wavelength grid:

      INTEGER :: nw, iw, nwint, nbins
      REAL, allocatable :: wl(:), wc(:), wu(:), wdelta(:)
      REAL    :: wstart, wstop

* Altitude grid

      INTEGER :: nz, nlyr, iz, izout
      REAL    :: zstart, zstop, zout
      REAL, allocatable :: z(:)

* Solar zenith angle and azimuth
* slant pathlengths in spherical geometry

      INTEGER, allocatable :: nid(:)
      REAL                 :: zen, sznoon
      REAL, allocatable    :: sza(:)
      REAL, allocatable    :: dsdh(:,:)

* Extra terrestrial solar flux and earth-Sun distance ^-2

      REAL, allocatable :: f(:), etf(:), esfact(:)

* Ozone absorption cross section

      INTEGER :: mabs
      INTEGER :: strpos
      REAL, allocatable :: o3xs(:,:)

* O2 absorption cross section

      REAL, allocatable :: o2xs(:,:), o2xs1(:)

* SO2 absorption cross section
     
      REAL, allocatable :: so2xs(:)

* NO2 absorption cross section
     
      REAL, allocatable :: no2xs(:,:)

* Atmospheric optical parameters

      REAL, allocatable :: co3(:)
      REAL, allocatable :: albedo(:)
      REAL, allocatable :: dtrl(:,:)
      REAL, allocatable :: dto3(:,:), dto2(:,:)
      REAL, allocatable :: dtso2(:,:), dtno2(:,:)
      REAL, allocatable :: dtcld(:,:), omcld(:,:), gcld(:,:)
      REAL, allocatable :: dtaer(:,:), omaer(:,:), gaer(:,:)
      REAL, allocatable :: dtsnw(:,:), omsnw(:,:), gsnw(:,:)
      REAL, allocatable :: dt_any(:,:), om_any(:,:), g_any(:,:)
      REAL, allocatable :: aircon(:), aircol(:), aircol_lasrb(:)
      REAL, allocatable :: vcolo2(:), scolo2(:)
      REAL, allocatable :: tlev(:)
      REAL, allocatable :: tlay(:)

* Spectral irradiance and actinic flux (scalar irradiance)

      REAL, allocatable :: edir(:), edn(:), eup(:)
      REAL, allocatable :: sirrad(:,:)
      REAL, allocatable :: fdir(:), fdn(:), fup(:)
      REAL, allocatable :: saflux(:,:)
      REAL, allocatable :: radField(:,:)

* Spectral weighting functions and weighted radiation

      INTEGER :: ns, is
      REAL :: drdw
      REAL :: sw(ks,kw), rate(ks,kz), dose(ks)
      CHARACTER(len=50) :: slabel(ks)
      CHARACTER(len=:), allocatable :: annotatedRate

* Photolysis coefficients (j-values)

      INTEGER :: nj, ij, nprate, nsw
      INTEGER :: tpflag(kj)
      REAL    :: djdw
      REAL    :: valj(kj,kz)
      REAL    :: sj(kj,kz,kw)
      REAL, allocatable :: xsqy(:,:,:)
      REAL, allocatable :: spcwght(:,:)
      REAL, allocatable :: jval_(:,:)
      CHARACTER(len=50) :: jlabel(kj)
      CHARACTER(len=50), allocatable :: annotatedjlabel(:)

*-----------------------------------------------------------------------------*
**** Re-scaling factors (can be read from input file)
* New surface albedo and surface pressure (milli bar)
* Total columns of O3, SO2, NO2 (Dobson Units)
* Cloud optical depth, altitude of base and top
* Aerosol optical depth at 550 nm, single scattering albedo, Angstrom alpha
*-----------------------------------------------------------------------------*
      REAL :: alsurf, psurf
      REAL :: o3_tc, so2_tc, no2_tc
      REAL :: taucld, zbase, ztop
      REAL :: tauaer, ssaaer, alpha

*-----------------------------------------------------------------------------*
* Location: Lat and Lon (deg.), surface elev (km)
* Altitude, temperature for specific outputs
*-----------------------------------------------------------------------------*
      REAL :: lat, lon
      REAL :: zaird, ztemp

*-----------------------------------------------------------------------------*
* Time and/or solar zenith angle
*-----------------------------------------------------------------------------*
      INTEGER :: iyear, imonth, iday
      INTEGER :: it, nt
      REAL    :: tstart, tstop
      REAL    :: tmzone
      REAL, allocatable :: t(:)
      LOGICAL :: lzenit

* number of radiation streams
      INTEGER :: nstr

* input/output control
      LOGICAL :: intrct
      CHARACTER(len=6) :: inpfil, outfil

      INTEGER :: iout

      REAL :: dirsun, difdn, difup

      CHARACTER(len=1) :: again

* Save arrays for output:

      LOGICAL :: lirrad, laflux, lrates, ljvals, lmmech
      INTEGER :: isfix, ijfix, itfix, izfix, iwfix
      INTEGER :: nms, nmj
      INTEGER :: ims(ks), imj(kj)

      REAL svj_zj(kz,kj), svj_tj(kt,kj), svj_zt(kz,kt)
      REAL svr_zs(kz,ks), svr_ts(kt,ks), svr_zt(kz,kt)
      REAL, allocatable :: svf_zw(:,:), svf_tw(:,:), svf_zt(:,:)
      REAL, allocatable :: svi_zw(:,:), svi_tw(:,:), svi_zt(:,:)

* Planetary boundary layer height and pollutant concentrations

      INTEGER :: ipbl
      REAL    :: zpbl
      REAL    :: o3pbl, so2pbl, no2pbl, aod330

* WRF-Chem output control

      LOGICAL :: wrfchm

* radiative transfer cross section configuration file

      LOGICAL            :: Obj_radXfer_xsects
      LOGICAL            :: Obj_photo_rates
      LOGICAL            :: Obj_spectral_wghts
      CHARACTER(len=256) :: radXfer_config_filespec
      CHARACTER(len=256) :: photo_rate_config_filespec
      CHARACTER(len=256) :: spectral_wght_config_filespec
      CHARACTER(len=256) :: command_option
      CHARACTER(len=2)   :: number
      type(string_t)     :: command_string
      type(string_t)     :: delim
      type(string_t), allocatable :: command_tokens(:)

***** Surface waters (lakes, ocean)
*   sdom = spectral absorption by Dissolved Organic Matter (DOM) 
*          in lakes and ocean
*   h2oabs = sdom at specific wavenght

      INTEGER :: jdom, jd
      CHARACTER(len=50) :: dlabel(kdom)
      REAL    :: sdom(kdom,kw)
      REAL    :: h2oabs

*     ydepth = depth in meters
*   irradiances normalized to unity incidence:
*     se_0 = solar beam irradiance just below surface
*     de_0 = diffuse irradiance just below surface (sum over all angles)
*     se_y = solar beam irradiance at ydepth
*     de_y = diffuse irradiance at ydepth  (sum over all angles)
*     se_int = integral of solar irradiance from surface to ydepth 
*     de_int = integral of diffuse irradiance from surface to ydepth
*   spectral irradiances (total = direct + diffuse-dn) in units of W m-2 nm-1:
*     we_0 = just below air-water interface
*     we_int = from 0 to ydepth

      REAL :: ydepth
      REAL :: se_0, de_0, se_y, de_y, se_int, de_int
      REAL :: we_0, we_int

*  in-water dose rates, doses - DNA weighted unless otherwise specified
*     dose in air just above surface, computed by time integration of rate(13,1)
*     wrate0, wdose0 = dose rate, dose, just below surface
*     wratei, wdosei = dose rate, dose, integrated from surface to ydepth

      REAL :: adose
      REAL :: wratei, wrate0
      REAL :: wdosei, wdose0

****** Other user-defined variables here:

* spectrum locator indices:
      integer :: js_dna, js_uvi, jd_dom
      integer :: id

* radiative transfer cross section object, config, environment
      type(environment_t)            :: environment
      type(config_t)                 :: radXfer_config
* photo rate config
      type(config_t)                  :: photo_rate_config
      type(photo_kinetics_t), pointer :: photo_kinetics
      CHARACTER(len=50)               :: header
* spectral wght object, config
      type(config_t)                 :: spectral_wght_config
      type(spectral_wght_warehouse_t), pointer :: 
     $       spectral_wght_warehouse

* radiators to include in calculations
      logical :: do_rayleigh, do_o2, do_o3, do_aerosols, do_clouds

* --- END OF DECLARATIONS ---------------------------------------------

* get radiative transfer cross section configuration filespec
      radXfer_config_filespec    = ' '
      photo_rate_config_filespec = ' '
      spectral_wght_config_filespec = ' '
      Obj_radXfer_xsects = .false.
      Obj_photo_rates    = .false.
      Obj_spectral_wghts = .false.
      do_rayleigh        = .false.
      do_o2              = .false.
      do_o3              = .false.
      do_aerosols        = .false.
      do_clouds          = .false.
      delim = '='
      do is = 1, COMMAND_ARGUMENT_COUNT( )
        CALL GET_COMMAND_ARGUMENT( is, command_option )
        command_string = trim( command_option )
        command_tokens = 
     $    command_string%split( delim,compress=.true. )
        command_tokens(1) = command_tokens(1)%to_upper()
        select case( command_tokens(1)%to_char() )
          case( 'RADXFER_CONFIG_FILESPEC' )
            write(*,*) 'Processing radXfer json config file'
            radXfer_config_filespec = command_tokens(2)%to_char()
            CALL radXfer_config%from_file( radXfer_config_filespec )
            Obj_radXfer_xsects = .true.
          case( 'PHOTO_RATE_CONFIG_FILESPEC' )
            write(*,*) 'Processing photo_rate json config file'
            photo_rate_config_filespec = command_tokens(2)%to_char()
            CALL 
     $       photo_rate_config%from_file( photo_rate_config_filespec )
            Obj_photo_rates = .true.
          case( 'SPECTRAL_WGHT_CONFIG_FILESPEC' )
            write(*,*) 'Processing spectral_wght json config file'
            spectral_wght_config_filespec = command_tokens(2)%to_char()
            CALL 
     $       spectral_wght_config%from_file( 
     $                     spectral_wght_config_filespec )
            Obj_spectral_wghts = .true.
          case( 'DO_RAYLEIGH' )
            do_rayleigh = .true.
          case( 'DO_O2' )
            do_o2       = .true.
          case( 'DO_O3' )
            do_o3       = .true.
          case( 'DO_AEROSOLS' )
            do_aerosols = .true.
          case( 'DO_CLOUDS' )
            do_clouds   = .true.
          case default
            write(*,*) 'tuv: ',trim(command_tokens(1)%to_char())
            write(*,*) '       is not a valid keyword'
            Stop 'Arg Error'
        end select
      enddo

      write(*,*) 
     $  'TUV: uses xsect objects in radXfer = ',Obj_radXfer_xsects
      write(*,*) 
     $  'TUV: uses photo_rate objects = ',Obj_photo_rates
      write(*,*) 
     $  'TUV: uses spectral wght objects = ',Obj_spectral_wghts

* re-entry point

 1000 CONTINUE

* Open log file:

      OPEN(UNIT=kout,FILE='odat/OUTPUTS/'//'tuvlog'//'.txt',
     &     STATUS='UNKNOWN')

* ___ SECTION 1: SIMPLE INPUT VARIABLES --------------------------------
******* Read simple input variables from a file:

* can read interactively (intrct = .TRUE.) 
* or in batch mode (intrct = .FALSE.)

      intrct = .TRUE.
      IF ( .NOT. intrct) inpfil = 'usrinp'

      CALL rdinp(intrct, 
     $     inpfil, outfil, nstr,   lat,    lon,    tmzone,
     $     iyear,  imonth, iday,   zstart, zstop,  nz,
     $     wstart, wstop,  nwint,  tstart, tstop,  nt,
     $     lzenit, alsurf, psurf,  o3_tc,  so2_tc, no2_tc,
     $     taucld, zbase,  ztop,   tauaer, ssaaer, alpha,
     $     dirsun, difdn,  difup,  zout,   zaird,  ztemp,
     $     lirrad, laflux, lmmech, lrates, isfix,  nms,
     $     ljvals, ijfix,  nmj,    iwfix,  itfix,  izfix,
     $     ims,    slabel, imj,    jlabel)

      IF(outfil .EQ. 'screen') THEN
         iout = 6
      ELSE
         iout = 30
      ENDIF         

      nlyr = nz - 1                                    ! number of layers
      allocate( nid(0:nlyr) )
      allocate( dsdh(0:nlyr,nlyr) )

************* Can overwrite basic inputs here manually:
* Input and output files:
*   inpfil = input file name
*   outfil = output file name
* Radiative transfer scheme:
*   nstr = number of streams
*          If nstr < 2, will use 2-stream Delta Eddington
*          If nstr > 1, will use nstr-stream discrete ordinates
* Location (geographic):
*   lat = LATITUDE (degrees, North = positive)
*   lon = LONGITUDE (degrees, East = positive)
*   tmzone = Local time zone difference (hrs) from Universal Time (ut):  
*            ut = timloc - tmzone
* Date:
*   iyear = year (1950 to 2050)
*   imonth = month (1 to 12)
*   iday = day of month
* Time of day grid:
*   tstart = starting time, local hours
*   tstop = stopping time, local hours
*   nt = number of time steps
*   lzenit = switch for solar zenith angle (sza) grid rather than time 
*             grid. If lzenit = .TRUE. then 
*                tstart = first sza in deg., 
*                tstop = last sza in deg., 
*                nt = number of sza steps. 
*                esfact = 1. (Earth-sun distance = 1.000 AU)
* Vertical grid:
*   zstart = surface elevation above sea level, km
*   zstop = top of the atmosphere (exospheric), km
*   nz = number of vertical levels, equally spaced
*        (nz will increase by +1 if zout does not match altitude grid)
* Wavlength grid:
*   wstart = starting wavelength, nm
*   wstop  = final wavelength, nm
*   nwint = number of wavelength intervals, equally spaced
*           if nwint < 0, the standard atmospheric wavelength grid, not
*           equally spaced, from 120 to 735 nm, will be used. In this
*           case, wstart and wstop values are ignored.
* Surface condition:
*   alsurf = surface albedo, wavelength independent
*   psurf = surface pressure, mbar.  Set to negative value to use
*           US Standard Atmosphere, 1976 (USSA76)
* Column amounts of absorbers (in Dobson Units, from surface to space):
*          Vertical profile for O3 from USSA76.  For SO2 and NO2, vertical
*          concentration profile is 2.69e10 molec cm-3 between 0 and 
*          1 km above sea level, very small residual (10/largest) above 1 km.
*   o3_tc = ozone (O3)
*   so2_tc = sulfur dioxide (SO2)
*   no2_tc = nitrogen dioxide (NO2)
* Cloud, assumed horizontally uniform, total coverage, single scattering
*         albedo = 0.9999, asymmetry factor = 0.85, indep. of wavelength,
*         and also uniform vertically between zbase and ztop:
*   taucld = vertical optical depth, independent of wavelength
*   zbase = altitude of base, km above sea level
*   ztop = altitude of top, km above sea level
* Aerosols, assumed vertical provile typical of continental regions from
*         Elterman (1968):
*   tauaer = aerosol vertical optical depth at 550 nm, from surface to space. 
*           If negative, will default to Elterman's values (ca. 0.235 
*           at 550 nm).
*   ssaaer = single scattering albedo of aerosols, wavelength-independent.
*   alpha = Angstrom coefficient = exponent for wavelength dependence of 
*           tauaer, so that  tauaer1/tauaer2  = (w2/w1)**alpha.
* Directional components of radiation, weighting factors:
*   dirsun = direct sun
*   difdn = down-welling diffuse
*   difup = up-welling diffuse
*        e.g. use:
*        dirsun = difdn = 1.0, difup = 0 for total down-welling irradiance
*        dirsun = difdn = difup = 1.0 for actinic flux from all directions
*        dirsun = difdn = 1.0, difup = -1 for net irradiance
* Output altitude:
*   zout = altitude, km, for desired output.
*        If not within 1 m of altitude grid, an additional
*        level will be inserted and nz will be increased by +1.
*   zaird = air density (molec. cm-3) at zout.  Set to negative value for
*        default USSA76 value interpolated to zout.
*   ztemp = air temperature (K) at zout.  Set to negative value for
*        default USSA76 value interpolated to zout.
* Output options, logical switches:
*   lirrad = output spectral irradiance
*   laflux = output spectral actinic flux
*   lmmech = output for NCAR Master Mechanism use
*   lrates = output dose rates (UVB, UVA, CIE/erythema, etc.)
* Output options, integer selections:
*   isfix:  if > 0, output dose rate for action spectrum is=isfix, tabulated
*           for different times and altitudes.
*   ijfix:  if > 0, output j-values for reaction ij=ijfix, tabulated
*           for different times and altitudes.
*   iwfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at wavelength iw=iwfix, tabulated for different times
*           and altitudes.
*   itfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at time it=itfix, tabulated for different altitudes
*           and wavelengths.
*   izfix:  if > 0, output spectral irradiance and/or spectral actinic
*           flux at altitude iz=izfix, tabulated for different times
*           and wavelengths.
*   nms:    number of dose rates that will be reported. Selections must be 
*           made interactively, or by editing input file.
*   nmj:    number of j-values that will be reported. Selections must be 
*           made interactively, or by editing input file.

      IF(nstr < 2) THEN
         WRITE(kout,*) 'Delta-Eddington 2-stream radiative transfer' 
      ELSE
         WRITE(kout,*) 'Discrete ordinates ', 
     $        nstr, '-stream radiative transfer' 
      ENDIF

      WRITE(*,*) 'calculating....'

* ___ SECTION 2: SET GRIDS _________________________________________________

* altitudes (creates altitude grid, locates index for selected output, izout)

      CALL gridz(zstart, zstop, nz, z, zout, izout)
      IF(izfix > 0) izout = izfix

* time/zenith (creates time/zenith angle grid, starting at tstart)

      allocate( sza(nt), esfact(nt), t(nt) )
      CALL gridt(lat, lon, tmzone,
     $     iyear, imonth, iday,
     $     lzenit, tstart, tstop,
     $     t, sza, sznoon, esfact)


* wavelength grid, user-set range and spacing. 
* NOTE:  Wavelengths are in vacuum, and therefore independent of altitude.
* To use wavelengths in air, see options in subroutine gridw

      CALL gridw(wstart, wstop, nwint, wl, wc, wu)
      nbins = size(wc)
      nw    = nbins + 1

      allocate( etf(nbins) )
      allocate( wdelta(nbins) )
      allocate( svf_zw(nz,nbins), svf_tw(nt,nbins), svf_zt(nz,nt) )
      allocate( svi_zw(nz,nbins), svi_tw(nt,nbins), svi_zt(nz,nt) )

      wdelta = wl(2:) - wl(:nbins)

* ___ SECTION 3: SET UP VERTICAL PROFILES OF TEMPERATURE, AIR DENSITY, and OZONE

***** Temperature vertical profile, Kelvin 
*   can overwrite temperature at altitude z(izout)

      allocate( tlev(nz),tlay(nlyr) )
      CALL vptmp(z, tlev,tlay)
      IF(ztemp > nzero) tlev(izout) = ztemp
      call diagout( 'vptmp.old', tlev )

*****  Air density (molec cm-3) vertical profile 
*   can overwrite air density at altitude z(izout)

      allocate( aircon(nz), aircol(nlyr), aircol_lasrb(nz) )
      CALL vpair(psurf, z, aircon, aircol, aircol_lasrb)
      IF(zaird > nzero) aircon(izout) = zaird
      call diagout( 'vpair.old', aircon )

*****
*! PBL pollutants will be added if zpbl > 0.
* CAUTIONS:  
* 1. The top of the PBL, zpbl in km, should be on one of the z-grid altitudes.
* 2. Concentrations, column increments, and optical depths
*       will be overwritten between surface and zpbl.
* 3. Inserting PBL constituents may change their total column amount.
* 4. Above pbl, the following are used:
*       for O3:  USSA or other profile
*       for NO2 and SO2: set to zero.
*       for aerosols: Elterman
* Turning on pbl will affect subroutines:
* vpo3, setno2, setso2, and setaer. See there for details

      zpbl = -999.

* locate z-index for top of pbl

      ipbl = 0
      IF(zpbl > rZERO) THEN
         DO iz = 1, nlyr
            IF(z(iz+1) > (z(1) + zpbl*1.00001)) THEN
               EXIT
            ENDIF
         ENDDO

         ipbl = iz - 1
         write(*,*) 'top of PBL index, height (km) ', ipbl, z(ipbl)

* specify pbl concetrations, in parts per billion

         o3pbl = 100.
         so2pbl = 10.
         no2pbl = 50.

* PBL aerosol optical depth at 330 nm
* (to change ssa and g of pbl aerosols, go to subroutine setair.f)

         aod330 = 0.8

      ENDIF

***** Ozone vertical profile

      allocate( co3(nlyr) )
      co3 = vpo3(ipbl, zpbl, o3pbl, o3_tc, z, aircol)
      call diagout( 'vpco3.old', co3 )

* ___ SECTION 4: READ SPECTRAL DATA ____________________________

* read (and grid) extra terrestrial flux data:
      
      allocate( f(nbins) )
      CALL rdetfl(nw,wl, f)

* read cross section data for 
*    O3 (temperature-dependent)
*    O2 (will overwrite at Lyman-alpha and SRB wavelengths
*            see subroutine la_srb.f)
*    NO2 (temperature dependent)
*    SO2 

      allocate( o2xs1(nbins) )
      allocate( o3xs(nlyr,nbins) )
      allocate( so2xs(nbins) )
      allocate( no2xs(nlyr,nbins) )
      no2xs = rZERO
*    standard TUV cross section handling
*     if( .not. Obj_radXfer_xsects ) then
        mabs = 1
        CALL rdo3xs(mabs,nlyr,tlay,nw,wl, o3xs)
        CALL rdo2xs(nw,wl, o2xs1)
*       CALL rdno2xs(nz,tlay,nw,wl, no2xs)
        CALL rdso2xs(nw,wl, so2xs)
        call diagout( 'o3xs.old',o3xs )
        OPEN(unit=44,file='odat/OUTPUTS/o3xs_old',form='unformatted')
        WRITE(unit=44) o3xs
        CLOSE(unit=44)
*    new cross section objects
*     else
*    intialize rad transfer cross section object
*       radXfer_xsect_warehouse => radXfer_xsect_warehouse_t( 
*    $                   radXfer_config,real(wl,kind=musica_dk) )
*       do iz = 1,nlyr
*         environment%temperature = real(tlay(iz),kind=musica_dk)
*         environment%number_density_air = 
*    $          real(aircon(iz),kind=musica_dk)
*         CALL 
*    $     radXfer_xsect_warehouse%update_for_new_environmental_state
*    $              ( environment )
*         if( iz == 1 ) then
*           o2xs1 = 
*    $       real( radXfer_xsect_warehouse%cross_section_values_(:,2),
*    $             kind=musica_rk )
*           so2xs = 
*    $       real( radXfer_xsect_warehouse%cross_section_values_(:,4),
*    $             kind=musica_rk )
*         endif
*         o3xs(iz,:) = 
*    $       real( radXfer_xsect_warehouse%cross_section_values_(:,1),
*    $             kind=musica_rk )
*         no2xs(iz,:) = 
*    $       real( radXfer_xsect_warehouse%cross_section_values_(:,3),
*    $             kind=musica_rk )
*       enddo
*       OPEN(unit=44,file='odat/OUTPUTS/o3xs_new',form='unformatted')
*       WRITE(unit=44) o3xs
*       CLOSE(unit=44)
*     endif

****** Spectral weighting functions 
* (Some of these depend on temperature T and pressure P, and therefore
*  on altitude z.  Therefore they are computed only after the T and P profiles
*  are set above with subroutines settmp and setair.)
* Photo-physical   set in swphys.f (transmission functions)
* Photo-biological set in swbiol.f (action spectra)
* Photo-chemical   set in swchem.f (cross sections x quantum yields)* Physical 
*   and biological weigthing functions are assumed to depend
*   only on wavelength.
* Chemical weighting functions (product of cross-section x quantum yield)
*   for many photolysis reactions are known to depend on temperature
*   and/or pressure, and therefore are functions of wavelength and altitude.
* Output:
* from swphys & swbiol:  sw(ks,kw) - for each weighting function slabel(ks)
* from swchem:  sj(kj,kz,kw) - for each reaction jlabel(kj)
* For swchem, need to know temperature and pressure profiles.

      CALL swphys(nw,wl,wc, ns,sw,slabel)
      CALL swbiol(nw,wl,wc, ns,sw,slabel)
      OPEN(unit=33,file='odat/OUTPUTS/annotatedslabels.old',
     &form='formatted')
      DO IS = 1,NS
        WRITE(33,'(a)') trim(slabel(is))
      ENDDO
      CLOSE(unit=33)

      IF( .not. Obj_photo_rates ) THEN
        CALL swchem(nw,wl,nz,tlev,aircon, nj,sj,jlabel,tpflag)
      ENDIF

* output spectral weights
      OPEN(unit=44,file='odat/OUTPUTS/sw_org',form='unformatted')
      do is = 1,ns
        header = slabel(is)
        WRITE(44) header
        WRITE(44) sw(is,1:nw-1)
      enddo
      CLOSE(unit=44)

** Read other spectral data
* absorption coefficients for Dissolved Organic Matter (DOM) in surface waters

      CALL swdom(nw,wl,wc, jdom,dlabel,sdom)

* locate indices for some specific spectra:

      js_dna = 0
      js_uvi = 0
      DO is = 1, ns
         IF(slabel(is) .EQ. 
     $        'DNA damage, in vitro (Setlow, 1974)               ') 
     $        js_dna = is

         IF(slabel(is) .EQ. 
     $        'UV index (WMO, 1994; Webb et al., 2011)')           
     $        js_uvi = is
      ENDDO

      jd_dom = 0
      DO jd = 1, jdom
         if(dlabel(jd) .eq. 
     $        'Generic DOM absorption')
     $        jd_dom = jd
      ENDDO

c      write(*,*) js_dna, js_uvi, jd_dom

**** The following CALL is normally commented out.
* Subroutine newlst regenerates the list of weighting functions 
* (molecular and biological spectra) when new ones are added, to 
* update the default input files (defin1, defin2. etc.).  User
* input files, e.g. usrinp, should be similarly updated. 
* The program STOPS at the completion of newlst.
* If not in use, newlst.o can be safely removed from Makefile.

c      CALL newlst(ns,slabel,nj,jlabel)

**** Option for writing look-up tables of 
* (molecular cross sections x quantum yields) 
* for WRF-Chem, at selected temperatures and pressures. 
* STOPs after tables are written.

      wrfchm = .FALSE.
      IF (inpfil .EQ. 'defin5') wrfchm = .TRUE.
      IF (wrfchm) CALL wrflut(nw, wl, nz, tlev, aircon)

* ___ SECTION 5: SET ATMOSPHERIC OPTICAL DEPTH INCREMENTS _____________________

* Rayleigh optical depth increments:

      allocate( dtrl(nlyr,nbins) )
      dtrl = odrl( wc, aircol )
      call diagout( 'dtrl.old', dtrl )
      
* O2 vertical profile and O2 absorption optical depths
*   For now, O2 densitiy assumed as 20.95% of air density, can be changed
*   in subroutine.
*   Optical depths in Lyman-alpha and SRB will be over-written
*   in subroutine la_srb.f

      allocate( dto2(nlyr,nbins) )
      dto2 = seto2(z,wl,aircol,o2xs1)

* Ozone optical depths

      allocate( dto3(nlyr,nbins) )
      dto3 = odo3(z,wl,o3xs,co3)
      deallocate( co3 )
      call diagout( 'dto3.old', dto3 )

* SO2 vertical profile and optical depths

      allocate( dtso2(nlyr,nbins) )
      dtso2 = setso2(ipbl, zpbl, so2pbl,
     $     so2_tc, z, nbins, so2xs, tlay, aircol)

* NO2 vertical profile and optical depths

      allocate( dtno2(nlyr,nbins) )
      dtno2 = setno2(ipbl, zpbl, no2pbl, 
     $     no2_tc, z, nbins, no2xs, tlay, aircol)

* Cloud vertical profile, optical depths, single scattering albedo, asymmetry factor
      allocate( dtcld(nlyr,nbins),omcld(nlyr,nbins),gcld(nlyr,nbins) )
      CALL setcld(taucld,zbase,ztop,
     $            z,wl, dtcld,omcld,gcld)

* Aerosol vertical profile, optical depths, single scattering albedo, asymmetry factor

      allocate( dtaer(nlyr,nbins),omaer(nlyr,nbins),gaer(nlyr,nbins) )
      CALL setaer(ipbl, zpbl, aod330,
     $     tauaer, ssaaer, alpha,
     $     z, wc, dtaer, omaer, gaer)
      call diagout( 'dtaer.old', dtaer )

* Snowpack physical and optical depths, single scattering albedo, asymmetry factor

      allocate( dtsnw(nlyr,nbins),omsnw(nlyr,nbins),gsnw(nlyr,nbins) )
      CALL setsnw(z, wl, dtsnw, omsnw, gsnw)

* Surface albedo

      allocate( albedo(nbins) )
      albedo = setalb(alsurf,wl)

* Set any additional absorber or scatterer:
* Must populate dt_any(kz,kw), om_any(kz,kw), g_any(kz,kw) manually
* This allows user to put in arbitrary absorber or scatterer
* could write a subroutine, e.g.:
C      CALL setany(nz,z,nw,wl,aircol, dt_any,om_any, g_any)
* or write manually here.

      allocate(
     $   dt_any(nlyr,nbins),om_any(nlyr,nbins),g_any(nlyr,nbins) )
      dt_any = rZERO
      om_any = rZERO
      g_any  = rZERO

* ___ SECTION 6: TIME/SZA LOOP  _____________________________________

* Initialize any time-integrated quantities here

      adose = rZERO
      wdose0 = rZERO
      wdosei = rZERO
      dose(1:ks) = rZERO

      write(*,*) 'Date, Lat, Lon, Min_SZA'
      write(*,222) iyear,imonth,iday,lat,lon,sznoon
 222  format(i4,'/',i2,'/',i2,3(1x,F7.3))

* Initialize lymana-alpha, schumann-runge bands
      call init_la_srb(wl)

* Initialize photo rate constant objects
      if( Obj_photo_rates ) then
        photo_kinetics => 
     $    photo_kinetics_t( photo_rate_config,real(wl,kind=musica_dk) )
        nprate = size(photo_kinetics%cross_section_objs_)
        if( .not. allocated( xsqy ) ) then
          allocate( xsqy(nz,nbins,nprate) )
        endif
        if( .not. allocated( jval_ ) ) then
          allocate( jval_(nz,nprate) )
        endif
* get cross section, quantum yield; form product
        do iz = 1,nz
          environment%temperature = real(tlev(iz),kind=musica_dk)
          environment%number_density_air = 
     $          real(aircon(iz),kind=musica_dk)
          CALL photo_kinetics%update_for_new_environmental_state( 
     $              environment, nbins )
          do ij = 1,nprate
            xsqy(iz,:,ij) = real( 
     $                      photo_kinetics%cross_section_values_(:,ij)
     $                    * photo_kinetics%quantum_yield_values_(:,ij),
     $                      kind=musica_rk )
          enddo
        enddo
      endif
* Initialize spectral wght objects
      if( Obj_spectral_wghts ) then
        spectral_wght_warehouse => 
     $    spectral_wght_warehouse_t( 
     $        spectral_wght_config,real(wl,kind=musica_dk) )
        nsw = size(spectral_wght_warehouse%spectral_wght_objs_)
        if( .not. allocated( spcwght ) ) then
          allocate( spcwght(nbins,nsw) )
        endif
        environment%temperature = real(tlev(1),kind=musica_dk)
        environment%number_density_air = 
     $          real(aircon(1),kind=musica_dk)
* get spectral weights
        CALL 
     $  spectral_wght_warehouse%update_for_new_environmental_state
     $              (environment, nbins )
        do is = 1,nsw
            spcwght(:,is) = real( 
     $        spectral_wght_warehouse%spectral_wght_values_(:,is)
     $        ,kind=musica_rk )
        enddo
* output spectral weights
        OPEN(unit=44,file='odat/OUTPUTS/sw_new',form='unformatted')
        do is = 1,nsw
          header = spectral_wght_warehouse%spectral_wght_key(is)
          WRITE(44) header
          WRITE(44) spcwght(:,is)
        enddo
        CLOSE(unit=44)
      endif

      allocate( edir(nz), edn(nz), eup(nz) )
      allocate( sirrad(nz,nbins) )
      allocate( fdir(nz), fdn(nz), fup(nz) )
      allocate( saflux(nz,nbins) )
      allocate( radField(nz,nbins) )


* Loop over time or solar zenith angle (zen):
      sza_loop: DO it = 1, nt

         zen = sza(it)

         WRITE(*,200) it, zen, esfact(it)
         WRITE(kout,200) it, zen, esfact(it)
 200     FORMAT('step = ', I4,' sza = ', F9.3, 
     $        ' Earth-sun factor = ', F10.7)


* correction for earth-sun distance
         etf(1:nbins) = f(1:nbins) * esfact(it)

* ____ SECTION 7: CALCULATE ZENITH ANGLE-DEPENDENT QUANTITIES __________

* slant path lengths for spherical geometry

         CALL sphers(z,zen, dsdh,nid)
         if( .not. allocated(vcolo2) ) then
           allocate( vcolo2(nlyr),scolo2(nz) )
         endif
         CALL airmas(dsdh,nid, aircol_lasrb,vcolo2,scolo2)

* Recalculate effective O2 optical depth and cross sections for Lyman-alpha
* and Schumann-Runge bands, must know zenith angle
* Then assign O2 cross section to sj(1,*,*)

         if( .not. allocated(o2xs) ) then
           allocate( o2xs(nz,nbins) )
         endif
         CALL la_srb(z,tlev,wl,vcolo2,scolo2,o2xs1,dto2,o2xs)
         CALL sjo2(o2xs,sj(1,:,:))
         call diagout( 'dto2.old',dto2 )
* Output cross section, quantum yield product for regression testing
         write(number,'(i2.2)') it
         OPEN(unit=33,
     $        file='odat/OUTPUTS/xsqy.'//number//'.old',
     $        form='unformatted')
         write(unit=33) 
     $     reshape( sj(1:nj,1:nz,1:nbins),(/nz,nbins,nj/),
     $              order=(/3,1,2/) )
         CLOSE(unit=33)
* Output spectral wgths for regression testing
         OPEN(unit=33,
     $        file='odat/OUTPUTS/sw.'//number//'.old',
     $        form='unformatted')
         write(unit=33) transpose( sw(1:ns,1:nbins) )
         CLOSE(unit=33)

         if( .not. do_rayleigh ) dtrl = rZERO
         if( .not. do_o2       ) dto2 = rZERO
         if( .not. do_o3       ) dto3 = rZERO
         if( .not. do_aerosols ) then
           dtaer = rZERO ; gaer = rZERO ; omaer = rZERO
         endif
         if( .not. do_clouds   ) then
           omcld = rZERO ; omsnw = rZERO
         endif
         if( all( dtrl == 0. ) ) write(*,*) 'TUV: dtrl = 0'
         if( all( dto3 == 0. ) ) write(*,*) 'TUV: dto3 = 0'
         if( all( dto2 == 0. ) ) write(*,*) 'TUV: dto2 = 0'
         if( all( dtso2 == 0. ) ) write(*,*) 'TUV: dtso2 = 0'
         if( all( dtno2 == 0. ) ) write(*,*) 'TUV: dtno2 = 0'
         if( all( dtcld == 0. ) ) write(*,*) 'TUV: dtcld = 0'
         if( all( omcld == 0. ) ) write(*,*) 'TUV: omcld = 0'
         if( all( gcld == 0. ) ) write(*,*) 'TUV: gcld = 0'
         if( all( dtaer == 0. ) ) write(*,*) 'TUV: dtaer = 0'
         if( all( omaer == 0. ) ) write(*,*) 'TUV: omaer = 0'
         if( all( gaer == 0. ) ) write(*,*) 'TUV: gaer = 0'
         if( all( dtsnw == 0. ) ) write(*,*) 'TUV: dtsnw = 0'
         if( all( omsnw == 0. ) ) write(*,*) 'TUV: omsnw = 0'
         if( all( gsnw == 0. ) ) write(*,*) 'TUV: gsnw = 0'

* ____ SECTION 8: WAVELENGTH LOOP ______________________________________

* initialize for wavelength integration

         rate = rZERO
         valj = rZERO
         
         wrate0 = rZERO
         wratei = rZERO

***** Main wavelength loop:
         wave_loop: DO iw = 1, nbins
*-----------------------------------------------------------------------------*
** monochromatic radiative transfer. Outputs are:
*  normalized irradiances     edir(iz), edn(iz), eup(iz) 
*  normalized actinic fluxes  fdir(iz), fdn(zi), fup(iz)
*  where 
*  dir = direct beam, dn = down-welling diffuse, up = up-welling diffuse
*-----------------------------------------------------------------------------*
            CALL rtlink(nstr,
     $           albedo(iw), zen,
     $           dsdh,nid,
     $           dtrl(:,iw),
     $           dto3(:,iw),
     $           dto2(:,iw),
     $           dtso2(:,iw),
     $           dtno2(:,iw),
     $           dtcld(:,iw), omcld(:,iw), gcld(:,iw),
     $           dtaer(:,iw),omaer(:,iw),gaer(:,iw),
     $           dtsnw(:,iw),omsnw(:,iw),gsnw(:,iw),
     $           dt_any(:,iw),om_any(:,iw),g_any(:,iw),
     $           edir, edn, eup, fdir, fdn, fup)

             radField(:,iw) = fdir(:) + fup(:) + fdn(:)
* Spectral irradiance, W m-2 nm-1
* for downwelling only, use difup = 0.
            sirrad(:,iw) = etf(iw) * 
     $        (dirsun*edir(:) + difdn*edn(:) + difup*eup(:))

* Spectral actinic flux, quanta s-1 nm-1 cm-2, all directions:
*    units conversion:  1.e-4 * (wc*1e-9) / hc
            saflux(:,iw) = etf(iw) * (1.e-13 * wc(iw) / hc) *
     $           (dirsun*fdir(:) + difdn*fdn(:) + difup*fup(:))

**** Save irradiances and actinic fluxes for output
            CALL saver1(it, itfix, iw, iwfix,  izout,
     $           sirrad, saflux,
     $           svi_zw, svf_zw, svi_zt, svf_zt, svi_tw, svf_tw)

            CYCLE wave_loop

*** Accumulate weighted integrals over wavelength, at all altitudes:
            DO iz = 1, nz
* Weighted irradiances (dose rates) W m-2
               DO is = 1, ns
                  drdw = sirrad(iz,iw) * sw(is,iw) 
                  rate(is,iz) = rate(is,iz) + drdw * (wu(iw) - wl(iw))
               ENDDO
* Photolysis rate coefficients (J-values) s-1
               DO ij = 1, nj
                  djdw = saflux(iz,iw) * sj(ij,iz,iw)
                  valj(ij,iz) = valj(ij,iz) + djdw * (wu(iw) - wl(iw))
               ENDDO
            ENDDO

************ In-water radiation:
*   Input:  
*     ydepth, in meters, for which radiation field is desired
*     h2oabs = absorption coeff of DOM in H2O at this wavelength, 1/meter
*     zen = solar zenith angle
*   Output from subroutine waters:
*     se_0 = solar beam irradiance just below surface
*     de_0 = diffuse irradiance just below surface (sum over all angles)
*     se_y = solar beam irradiance at ydepth
*     de_y = diffuse irradiance at ydepth  (sum over all angles)
*     se_i = integral of solar beam irradiance from surface to ydepth 
*     de_i = integral of diffuse irradiance from surface to ydepth
            ydepth = 1.
            h2oabs = sdom(jd_dom,iw)
            CALL waters(zen,h2oabs,ydepth, 
     $           se_0,de_0,se_y,de_y,se_int,de_int)
            
* calculate spectral irradiances in water:
* irradiance just below air-water interface:
            we_0 = etf(iw) * (se_0*edir(1)*dirsun + 
     $           de_0*edn(1)*difdn)

* average spectral irradiance in water, from 0 to ydepth = integral / ydepth
            we_int = etf(iw) * (se_int*edir(1)*dirsun + 
     $           de_int*edn(1)*difdn) / ydepth

* calculate DNA-weighted irradiance, just below surface and
* averaged from 0 to ydepth
            drdw = we_0 * sw(js_dna,iw)
            wrate0 = wrate0 + drdw * (wu(iw)-wl(iw))

            drdw = we_int * sw(js_dna,iw)
            wratei = wratei + drdw * (wu(iw)-wl(iw))
         ENDDO wave_loop

* compute photo rate constants
      if( it == 1 .and. Obj_photo_rates ) then
        OPEN(unit=44,file='odat/OUTPUTS/jval_new',form='unformatted')
        jval_ = rZERO
        do ij = 1,nprate
          do iw = 1,nbins
            jval_(:,ij) = jval_(:,ij) 
     $                  + saflux(:,iw)*xsqy(:,iw,ij)*wdelta(iw)
          enddo
          header = photo_kinetics%reaction_key(ij)
          WRITE(44) header
          WRITE(44) jval_(:,ij)
        enddo
        CLOSE(unit=44)
      endif

         write(number,'(i2.2)') it
         call diagout( 'radField.' // number // '.old',radField )
         CYCLE sza_loop

**** integrate doses over time: 
* adose = dose in air just above surface
* wdose0 = dose in water just below surface
* wdosei = dose averaged over ydepth
         adose = adose + 
     $        rate(13,1) * 3600.* (tstop - tstart)/float(nt-1)
         wdose0 = wdose0 + 
     $        wrate0 * 3600.* (tstop - tstart)/float(nt-1)
         wdosei = wdosei + 
     $        wratei * 3600.* (tstop - tstart)/float(nt-1)

* Save dose rates and j-values for output
         CALL saver2(it,itfix, nz,izout, ns,isfix,ims, nj,ijfix,imj,
     $        rate, valj,
     $        svr_zs, svj_zj, svr_zt, svj_zt, svr_ts, svj_tj)

      ENDDO sza_loop

* output in-water doses

c      write(*,222) adose, wdose0, wdosei, dlabel(jd_dom)
c 222  format(3(0pf10.4,1x),a50)

**output all Js at zout

c      do iz = 1, nz
c         if(z(iz) .eq. zout) then
c            do ij = 1, nj
c               write(44,444) valj(ij,iz), jlabel(ij)
c            enddo
c         endif
c      enddo
c 444  format(1pe11.4,1x,a50)
*^^^^^^^^^^^^^^^^ end time/zenith loop

* ____ SECTION 9: OUTPUT ______________________________________________

      call outpt1( outfil, iout, 
     $     lirrad, laflux, lrates, ljvals, lmmech, lzenit,
     $     nms, ims, nmj, imj,
     $     z, tlev, aircon, izout,
     $     wl, etf, iwfix,
     $     t, sza, itfix,
     $     ns, slabel, isfix, nj, jlabel, ijfix,
     $     svj_zj, svj_tj, svj_zt,
     $     svr_zs, svr_ts, svr_zt,
     $     svf_zw, svf_tw, svf_zt,
     $     svi_zw, svi_tw, svi_zt )

*_______________________________________________________________________

      IF(intrct) THEN
         WRITE(*,*) 'do you want to do another calculation?'
         WRITE(*,*) 'y = yes'
         WRITE(*,*) 'any other key = no'
         READ(*,1001) again
 1001    FORMAT(A1)
         IF(again .EQ. 'y' .OR. again .EQ. 'Y') GO TO 1000
      ENDIF

      CLOSE(iout)
      CLOSE(kout)

      CONTAINS

      function get_photorate_index( photorateNames, match ) 
     $                      result( index )

      character(len=*), intent(in) :: match
      character(len=*), intent(in) :: photorateNames(:)
      integer                      :: index

      integer :: n

      index = -1
      do n = 1,size(photorateNames)
        if( trim(match) == trim(photorateNames(n)) ) then
          index = n
          exit
        endif
      enddo

      end function get_photorate_index

      function rmspace( inString ) result( outString  )

      character(len=*), intent(in) :: inString
      character(len=:), allocatable  :: outString

      integer   :: n
      character :: chr

      do n = 1,len_trim(inString)
        chr = inString(n:n)
        if( chr /= ' ' ) then
          outString = outString // chr
        endif
      enddo

      end function rmspace

      END PROGRAM tuv
