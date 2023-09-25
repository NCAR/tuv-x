! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!solar zenith angle from time type
module micm_earth_sun_distance

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_Profile,     only : base_profile_t
  use musica_assert,    only : die_msg

  implicit none

  private
  public :: earth_sun_distance_t

  type, extends(base_profile_t) :: earth_sun_distance_t
  contains
    !> Initialize grid
    procedure :: initialize
    final     :: finalize
  end type earth_sun_distance_t

contains
  !> Initialize distance between sun, earth in AU
  subroutine initialize( this, profile_config, gridWareHouse )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use micm_1d_grid,  only : base_grid_t
    use micm_grid_warehouse,  only : grid_warehouse_t

    !> Arguments
    class(earth_sun_distance_t), intent(inout) :: this
    type(config_t), intent(inout)         :: profile_config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> Local variables
    real(dk), parameter :: NINETY  = 90._dk
    integer(ik) :: n, tNdx
    integer(ik) :: Year, Month, Day
    integer(ik) :: Jday
    real(dk)    :: tmzone, ut, soldst
    real(dk)    :: Lon, Lat
    character(len=*), parameter :: Iam = 'earth sun distance initialize: '
    type(string_t) :: Handle
    class(base_grid_t), pointer :: timeGrid

    Handle = 'Time, hrs'
    timeGrid => gridWareHouse%get_grid( Handle )
    this%ncells_ = timeGrid%ncells_

    call profile_config%get( 'Handle', this%handle_, Iam, default='None' )

    allocate( this%edge_val_(0) )

    !> Map solar zenith angle as function of time
    call profile_config%get( 'Year', Year, Iam, default=2002 )
    call profile_config%get( 'Month', Month, Iam, default=3 )
    call profile_config%get( 'Day', Day, Iam, default=21 )
    call profile_config%get( 'TimeZone', tmzone, Iam, default=0.0_dk )
    call profile_config%get( 'Longitude', Lon, Iam, default=0.0_dk )
    call profile_config%get( 'Latitude', Lat, Iam, default=0.0_dk )
    Jday = JulianDayofYear(Year, Month, Day )
    do tNdx = 1_ik,this%ncells_+1_ik
      ut = timeGrid%edge_(tNdx) - tmzone
      soldst = earthSunDistance(Year, Jday, ut, Lat, Lon )
      this%edge_val_  = [this%edge_val_,soldst]
    enddo
    this%mid_val_ = .5_dk &
                  *(this%edge_val_(1_ik:this%ncells_) + this%edge_val_(2_ik:this%ncells_+1_ik))
    this%delta_val_ = (this%edge_val_(2_ik:this%ncells_+1_ik) - this%edge_val_(1_ik:this%ncells_))

  end subroutine initialize

  function JulianDayofYear( year, month, day ) result( julianday )

!-----------------------------------------------------------------------------*
!= calculates julian day corresponding to specified year, month, day         =*
!= also checks validity of date                                              =*
!-----------------------------------------------------------------------------*

      !> Arguments
      integer(ik), intent(in)  :: year, month, day

      integer(ik)              :: julianday

      !> Local variables
      integer(ik) :: m
      integer(ik) :: dayofyear
      integer(ik), allocatable :: daysinmonth(:)
      LOGICAL(lk) :: leapyr

      daysinmonth = (/31_ik,28_ik,31_ik,30_ik,31_ik,30_ik,31_ik,31_ik,30_ik,31_ik,30_ik,31_ik/)

      !> check date inputs
      IF(year < 1950_ik .OR. year > 2050_ik) THEN
        call die_msg( 460768716, "1950 <= Year <= 2050" )
      ENDIF

      IF(month < 1_ik .OR. month > 12_ik) THEN
        call die_msg( 460768717, "1 <= Month <= 12" )
      ENDIF

      !> leap year
      leapyr = MOD(year,4_ik) == 0_ik
      IF ( leapyr ) THEN
         daysinmonth(2) = 29_ik
      ELSE
         daysinmonth(2) = 28_ik
      ENDIF

      IF (day > daysinmonth(month)) THEN
        call die_msg( 460768718, "Day > days in month" )
      ENDIF

      dayofyear = sum( daysinmonth(1:month-1) )
      julianday = dayofyear + day

      end function JulianDayofYear

      function earthSunDistance( YEAR, DAY, HOUR, LAT, LONG ) result( soldst )

!     Calculates the local solar azimuth and elevation angles, and
!     the distance to and angle subtended by the Sun, at a specific 
!     location and time using approximate formulas in The Astronomical 
!     Almanac.  Accuracy of angles is 0.01 deg or better (the angular 
!     width of the Sun is about 0.5 deg, so 0.01 deg is more than
!     sufficient for most applications).

!     Unlike many GCM (and other) sun angle routines, this
!     one gives slightly different sun angles depending on
!     the year.  The difference is usually down in the 4th
!     significant digit but can slowly creep up to the 3rd
!     significant digit after several decades to a century.

!     A refraction correction appropriate for the "US Standard
!     Atmosphere" is added, so that the returned sun position is
!     the APPARENT one.  The correction is below 0.1 deg for solar
!     elevations above 9 deg.  To remove it, comment out the code
!     block where variable REFRAC is set, and replace it with
!     REFRAC = 0.0.

!   References:

!     Michalsky, J., 1988: The Astronomical Almanac's algorithm for
!        approximate solar position (1950-2050), Solar Energy 40,
!        227-235 (but the version of this program in the Appendix
!        contains errors and should not be used)

!     The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
!        D.C. (published every year): the formulas used from the 1995
!        version are as follows:
!        p. A12: approximation to sunrise/set times
!        p. B61: solar elevation ("altitude") and azimuth
!        p. B62: refraction correction
!        p. C24: mean longitude, mean anomaly, ecliptic longitude,
!                obliquity of ecliptic, right ascension, declination,
!                Earth-Sun distance, angular diameter of Sun
!        p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)


!   Authors:  Dr. Joe Michalsky (joe@asrc.albany.edu)
!             Dr. Lee Harrison (lee@asrc.albany.edu)
!             Atmospheric Sciences Research Center
!             State University of New York
!             Albany, New York

!   Modified by:  Dr. Warren Wiscombe (wiscombe@climate.gsfc.nasa.gov)
!                 NASA Goddard Space Flight Center
!                 Code 913
!                 Greenbelt, MD 20771


!   WARNING:  Do not use this routine outside the year range
!             1950-2050.  The approximations become increasingly
!             worse, and the calculation of Julian date becomes
!             more involved.

!   Input:

!      YEAR     year (integer; range 1950 to 2050)

!      DAY      day of year at LAT-LONG location (integer; range 1-366)

!      HOUR     hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0)
!               = (local hour) + (time zone number)
!                 + (Daylight Savings Time correction; -1 or 0)
!               where (local hour) range is 0 to 24,
!               (time zone number) range is -12 to +12, and
!               (Daylight Time correction) is -1 if on Daylight Time
!               (summer half of year), 0 otherwise;  
!               Example: 8:30 am Eastern Daylight Time would be
!
!                           HOUR = 8.5 + 5 - 1 = 12.5

!      LAT      latitude [degrees]
!               (REAL; range -90.0 to 90.0; north is positive)

!      LONG     longitude [degrees]
!               (REAL; range -180.0 to 180.0; east is positive)


!   Output:

!      EL       solar elevation angle [-90 to 90 degs]; 
!               solar zenith angle = 90 - EL

!      SOLDST   distance to sun [Astronomical Units, AU]
!               (1 AU = mean Earth-sun distance = 1.49597871E+11 m
!                in IAU 1976 System of Astronomical Constants)


!   Local Variables:

!     DEC       Declination (radians)

!     ECLONG    Ecliptic longitude (radians)

!     GMST      Greenwich mean sidereal time (hours)

!     HA        Hour angle (radians, -pi to pi)

!     JD        Modified Julian date (number of days, including 
!               fractions thereof, from Julian year J2000.0);
!               actual Julian date is JD + 2451545.0

!     LMST      Local mean sidereal time (radians)

!     MNANOM    Mean anomaly (radians, normalized to 0 to 2*pi)

!     MNLONG    Mean longitude of Sun, corrected for aberration 
!               (deg; normalized to 0-360)

!     OBLQEC    Obliquity of the ecliptic (radians)

!     RA        Right ascension  (radians)

!     REFRAC    Refraction correction for US Standard Atmosphere (degs)

! --------------------------------------------------------------------
!   Uses double precision for safety and because Julian dates can
!   have a large number of digits in their full form (but in practice
!   this version seems to work fine in single precision if you only
!   need about 3 significant digits and aren't doing precise climate
!   change or solar tracking work).
! --------------------------------------------------------------------

!   Why does this routine require time input as Greenwich Mean Time 
!   (GMT; also called Universal Time, UT) rather than "local time"?
!   Because "local time" (e.g. Eastern Standard Time) can be off by
!   up to half an hour from the actual local time (called "local mean
!   solar time").  For society's convenience, "local time" is held 
!   constant across each of 24 time zones (each technically 15 longitude
!   degrees wide although some are distorted, again for societal 
!   convenience).  Local mean solar time varies continuously around a 
!   longitude belt;  it is not a step function with 24 steps.  
!   Thus it is far simpler to calculate local mean solar time from GMT,
!   by adding 4 min for each degree of longitude the location is
!   east of the Greenwich meridian or subtracting 4 min for each degree
!   west of it.  

! --------------------------------------------------------------------

!   TIME
!   
!   The measurement of time has become a complicated topic.  A few
!   basic facts are:
!   
!   (1) The Gregorian calendar was introduced in 1582 to replace 
!   Julian calendar; in it, every year divisible by four is a leap 
!   year just as in the Julian calendar except for centurial years
!   which must be exactly divisible by 400 to be leap years.  Thus 
!   2000 is a leap year, but not 1900 or 2100.

!   (2) The Julian day begins at Greenwich noon whereas the calendar 
!   day begins at the preceding midnight;  and Julian years begin on
!   "Jan 0" which is really Greenwich noon on Dec 31.  True Julian 
!   dates are a continous count of day numbers beginning with JD 0 on 
!   1 Jan 4713 B.C.  The term "Julian date" is widely misused and few
!   people understand it; it is best avoided unless you want to study
!   the Astronomical Almanac and learn to use it correctly.

!   (3) Universal Time (UT), the basis of civil timekeeping, is
!   defined by a formula relating UT to GMST (Greenwich mean sidereal
!   time).  UTC (Coordinated Universal Time) is the time scale 
!   distributed by most broadcast time services.  UT, UTC, and other
!   related time measures are within a few sec of each other and are
!   frequently used interchangeably.

!   (4) Beginning in 1984, the "standard epoch" of the astronomical
!   coordinate system is Jan 1, 2000, 12 hr TDB (Julian date 
!   2,451,545.0, denoted J2000.0).  The fact that this routine uses
!   1949 as a point of reference is merely for numerical convenience.
! --------------------------------------------------------------------

      !> Arguments
      integer(ik), intent(in) ::  YEAR, DAY
      real(dk), intent(in)    ::  HOUR, LAT, LONG
      real(dk)                ::  soldst

      !> Local variables
      integer(ik), parameter :: iZERO = 0_ik
      real(dk), parameter :: Day2Hrs = 24._dk
      real(dk), parameter :: NINETY  = 90._dk
      real(dk), parameter :: THREE60 = 360._dk
      real(dk), parameter :: rZERO   = 0._dk
      real(dk), parameter :: rONE    = 1._dk
      real(dk), parameter :: rTWO    = 2._dk
      real(dk), parameter :: PI    = rTWO*ASIN( rONE )
      real(dk), parameter :: TWOPI = rTWO*PI
      real(dk), parameter :: RPD   = PI/180._dk

      integer(ik) ::  DELTA, LEAP
      real(dk) :: DEC, DEN, ECLONG, GMST, HA, JD, LMST, &
                  MNANOM, MNLONG, NUM, OBLQEC, RA, &
                  REFRAC, TIME

!                    ** current Julian date (actually add 2,400,000 
!                    ** for true JD);  LEAP = leap days since 1949;
!                    ** 32916.5 is midnite 0 jan 1949 minus 2.4e6
      DELTA  = YEAR - 1949_ik
      LEAP   = DELTA / 4_ik
      JD     = 32916.5_dk + real(DELTA*365_ik + LEAP + DAY,dk) + HOUR / Day2Hrs

!                    ** last yr of century not leap yr unless divisible
!                    ** by 400 (not executed for the allowed YEAR range,
!                    ** but left in so our successors can adapt this for 
!                    ** the following 100 years)
      if( MOD( YEAR, 100_ik ) == iZERO .AND. MOD( YEAR, 400_ik ) /= iZERO ) then
        JD = JD - rONE
      endif

!                     ** ecliptic coordinates
!                     ** 51545.0 + 2.4e6 = noon 1 jan 2000
      TIME  = JD - 51545.0_dk

!                    ** force mean longitude between 0 and 360 degs
      MNLONG = 280.460_dk + 0.9856474_dk*TIME
      MNLONG = MOD( MNLONG, THREE60 )
      IF( MNLONG < rZERO ) MNLONG = MNLONG + THREE60

!                    ** mean anomaly in radians between 0 and 2*pi
      MNANOM = 357.528_dk + 0.9856003_dk*TIME
      MNANOM = MOD( MNANOM, THREE60 )
      IF( MNANOM < rZERO ) MNANOM = MNANOM + THREE60

      MNANOM = MNANOM*RPD

!   ** distance to sun in A.U.
      SOLDST = 1.00014_dk - 0.01671_dk*COS(MNANOM) - 0.00014_dk*COS(rTWO*MNANOM)

      end function earthSunDistance

  subroutine finalize( this )

  type(earth_sun_distance_t), intent(inout) :: this

  if( allocated( this%edge_val_ ) ) then
    deallocate( this%edge_val_ )
  endif
  if( allocated( this%mid_val_ ) ) then
    deallocate( this%mid_val_ )
  endif
  if( allocated( this%delta_val_ ) ) then
    deallocate( this%delta_val_ )
  endif

  end subroutine finalize

end module micm_earth_sun_distance
