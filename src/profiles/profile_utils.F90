! Copyright (C) 2022 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_utils
  ! Utility function for use in profile types
  ! :f:func:`tuvx_profile_utils/earth_sun_distance` and
  ! :f:func:`tuvx_profile_utils/solar_zenith_angle` are based on 
  ! 
  ! Michalsky, J., 1988: The Astronomical Almanac's algorithm for
  ! approximate solar position (1950-2050), Solar Energy 40,
  ! 227-235 (but the version of this program in the Appendix
  ! contains errors and should not be used)
  ! `doi:10.1016/0038-092X(88)90045-X <https://doi.org/10.1016/0038-092X(88)90045-X>`_
  !
  ! The Astronomical Almanac, U.S. Gov't Printing Office, Washington,
  ! D.C. (published every year): the formulas used from the 1995
  ! version are as follows
  ! 
  ! - p. A12: approximation to sunrise/set times
  ! - p. B61: solar elevation ("altitude") and azimuth
  ! - p. B62: refraction correction
  ! - p. C24: mean longitude, mean anomaly, ecliptic longitude,
  !           obliquity of ecliptic, right ascension, declination,
  !           Earth-Sun distance, angular diameter of Sun
  ! - p. L2:  Greenwich mean sidereal time (ignoring T^2, T^3 terms)
  ! 
  ! These two functions calculate the local solar azimuth and
  ! elevation angles, and
  ! the distance to and angle subtended by the Sun, at a specific
  ! location and time using approximate formulas in The Astronomical
  ! Almanac.  Accuracy of angles is 0.01 deg or better (the angular
  ! width of the Sun is about 0.5 deg, so 0.01 deg is more than
  ! sufficient for most applications).
  !
  ! Unlike many GCM (and other) sun angle routines, this
  ! one gives slightly different sun angles depending on
  ! the year.  The difference is usually down in the 4th
  ! significant digit but can slowly creep up to the 3rd
  ! significant digit after several decades to a century.


  use musica_constants, only : &
    dk => musica_dk, ik => musica_ik, lk => musica_lk
  use musica_assert, only : die_msg

  implicit none

  private
  public :: julian_day_of_year, earth_sun_distance, solar_zenith_angle

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function julian_day_of_year( year, month, day ) result( julianday )
    ! Calculates julian day corresponding to specified year, month, day
    ! also checks validity of date     

    !> Arguments
    integer(ik), intent(in)  :: year, &  ! integer; range 1950 to 2050
                                month, & ! The integer month, e.g., 11
                                day      ! day of year at LAT-LONG location integer; range 1-366

    integer(ik)              :: julianday ! The `Julian day <https://en.wikipedia.org/wiki/Julian_day>`_

    !> Local variables
    integer(ik) :: dayofyear
    integer(ik), allocatable :: daysinmonth(:)

    daysinmonth = (/ &
      31_ik, 28_ik, 31_ik, & ! Jan, Feb, Mar
      30_ik, 31_ik, 30_ik, &
      31_ik, 31_ik, 30_ik, &
      31_ik, 30_ik, 31_ik  & ! ... Dec
    /)

    !> check date inputs
    IF(year < 1950_ik .OR. year > 2050_ik) THEN
      call die_msg( 460768716, "1950 <= Year <= 2050" )
    ENDIF

    IF(month < 1_ik .OR. month > 12_ik) THEN
      call die_msg( 460768717, "1 <= Month <= 12" )
    ENDIF

    !> leap year
    IF ( is_leap(year) ) THEN
        daysinmonth(2) = 29_ik
    ENDIF

    IF (day > daysinmonth(month)) THEN
      call die_msg( 460768718, "Day > days in month" )
    ENDIF

    dayofyear = sum( daysinmonth(1:month-1) )
    julianday = dayofyear + day

  end function julian_day_of_year

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function earth_sun_distance( year, day, hour ) result( solar_distance )
    ! Calculate the solar distance to the sun in AU

    !> arguments
    integer(ik), intent(in) ::  year, & ! integer; range 1950 to 2050
                                day     ! day of year at LAT-LONG location integer; range 1-366
    real(dk), intent(in)    ::  hour    ! hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0) (local hour) + (time zone number) + (Daylight Savings Time correction; -1 or 0) where (local hour) range is 0 to 24, (time zone number) range is -12 to +12, and (Daylight Time correction) is -1 if on Daylight Time (summer half of year), 0 otherwise; Example: 8:30 am Eastern Daylight Time would be HOUR = 8.5 + 5 - 1 = 12.5

    !> local variables
    real(dk) ::  solar_distance !  distance to sun [Astronomical Units, AU] (1 AU = mean Earth-sun distance = 1.49597871E+11 m in IAU 1976 System of Astronomical Constants)
    real(dk) :: mean_anomaly
    real(dk) :: time

    time  = calculate_time( year, day, hour )
    mean_anomaly = calculate_mean_anomaly( time )

    !   ** distance to sun in a.u.
    solar_distance = 1.00014_dk - 0.01671_dk*cos(mean_anomaly) &
      - 0.00014_dk*cos(2._dk*mean_anomaly)

  end function earth_sun_distance

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function solar_zenith_angle( year, day, hour, lat, long ) &
      result( solar_elevation )
      ! Calculates the solar zenith angle

    ! arguments
    integer(ik), intent(in) ::  year, & ! integer; range 1950 to 2050
                                day     ! day of year at LAT-LONG location integer; range 1-366
    real(dk), intent(in)    ::  hour, & !  hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0) (local hour) + (time zone number) + (Daylight Savings Time correction; -1 or 0) where (local hour) range is 0 to 24, (time zone number) range is -12 to +12, and (Daylight Time correction) is -1 if on Daylight Time (summer half of year), 0 otherwise; Example: 8:30 am Eastern Daylight Time would be HOUR = 8.5 + 5 - 1 = 12.5
                                lat,  & ! latitude [degrees], range -90.0 to 90.0; north is positive
                                long    ! longitude [degrees], range -180.0 to 180.0; east is positive

    ! local variables
    real(dk)            ::  solar_elevation ! angle [-90 to 90 degs]; solar zenith angle = 90 - EL
    real(dk), parameter :: day2hrs = 24._dk
    real(dk), parameter :: ninety  = 90._dk
    real(dk), parameter :: three60 = 360._dk
    real(dk), parameter :: rzero   = 0._dk
    real(dk), parameter :: rone    = 1._dk
    real(dk), parameter :: rtwo    = 2._dk
    real(dk), parameter :: pi    = 2._dk*asin( 1._dk )
    real(dk), parameter :: twopi = 2._dk*pi
    real(dk), parameter :: d2r   = pi/180._dk

    real(dk) :: mean_anomaly, mean_long
    real(dk) :: dec, eclong, ha, oblqec, ra, time

    time  = calculate_time( year, day, hour )

    mean_anomaly = calculate_mean_anomaly( time )
    mean_long = calculate_mean_long( time )

    oblqec = calculate_obliquity( time )
    eclong = calculate_ecliptic_longitude( mean_anomaly, mean_long )

    ra = calculate_right_ascension( eclong, oblqec )
    dec  = calculate_declination( eclong, oblqec )

    ha = calculate_hour_angle( time, hour, long, ra )

    !     noon when ha = 0
    solar_elevation  = asin( sin( dec )*sin( lat*d2r ) + &
      cos( dec )*cos( lat*d2r )*cos( ha ) )

    !                     ** convert elevation to degrees
    solar_elevation = solar_elevation / d2r

    !  ======== refraction correction for u.s. standard atmos. ==========
    if( solar_elevation < -ninety .or. solar_elevation > ninety ) then
      call die_msg( 460768815, "-90 =< solar elevation <= 90" )
    endif

  end function solar_zenith_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function calculate_time( year, day, hour ) result( time )
    ! Calculate the time.
    ! The measurement of time has become a complicated topic.  A few
    ! basic facts are:
    !
    ! (1) The Gregorian calendar was introduced in 1582 to replace
    ! Julian calendar; in it, every year divisible by four is a leap
    ! year just as in the Julian calendar except for centurial years
    ! which must be exactly divisible by 400 to be leap years.  Thus
    ! 2000 is a leap year, but not 1900 or 2100.
    !
    ! (2) The Julian day begins at Greenwich noon whereas the calendar
    ! day begins at the preceding midnight;  and Julian years begin on
    ! "Jan 0" which is really Greenwich noon on Dec 31.  True Julian
    ! dates are a continous count of day numbers beginning with JD 0 on
    ! 1 Jan 4713 B.C.  The term "Julian date" is widely misused and few
    ! people understand it; it is best avoided unless you want to study
    ! the Astronomical Almanac and learn to use it correctly.
    !
    ! (3) Universal Time (UT), the basis of civil timekeeping, is
    ! defined by a formula relating UT to GMST (Greenwich mean sidereal
    ! time).  UTC (Coordinated Universal Time) is the time scale
    ! distributed by most broadcast time services.  UT, UTC, and other
    ! related time measures are within a few sec of each other and are
    ! frequently used interchangeably.
    !
    ! (4) Beginning in 1984, the "standard epoch" of the astronomical
    ! coordinate system is Jan 1, 2000, 12 hr TDB (Julian date
    ! 2,451,545.0, denoted J2000.0).  The fact that this routine uses
    ! 1949 as a point of reference is merely for numerical convenience.

    integer(ik), intent(in) ::  year, & ! integer; range 1950 to 2050
                                day     ! day of year at LAT-LONG location integer; range 1-366
    real(dk), intent(in)    ::  hour    !  hour of DAY [GMT or UT] (REAL; range -13.0 to 36.0) (local hour) + (time zone number) + (Daylight Savings Time correction; -1 or 0) where (local hour) range is 0 to 24, (time zone number) range is -12 to +12, and (Daylight Time correction) is -1 if on Daylight Time (summer half of year), 0 otherwise; Example: 8:30 am Eastern Daylight Time would be HOUR = 8.5 + 5 - 1 = 12.5

    !> Local variables
    real(dk), parameter :: day_to_hours = 24._dk
    real(dk)    ::  time
    integer(ik) ::  delta, leap
    real(dk)    :: JD
    !     JD        Modified Julian date (number of days, including
    !               fractions thereof, from Julian year J2000.0);
    !               actual Julian date is JD + 2451545.0

    !                    ** current Julian date (actually add 2,400,000
    !                    ** for true JD);  leap = leap days since 1949;
    !                    ** 32916.5 is midnite 0 jan 1949 minus 2.4e6
    delta  = year - 1949_ik
    leap   = delta / 4_ik
    JD     = 32916.5_dk + real(delta*365_ik + leap + day,dk) + &
      hour / day_to_hours

    !                    ** last yr of century not leap yr unless divisible
    !                    ** by 400 (not executed for the allowed year range,
    !                    ** but left in so our successors can adapt this for
    !                    ** the following 100 years)
    if( mod( year, 100_ik ) == 0_ik .and. mod( year, 400_ik ) /= 0_ik ) then
      JD = JD - 1._dk
    endif

    !                     ** ecliptic coordinates
    !                     ** 51545.0 + 2.4e6 = noon 1 jan 2000
    time  = JD - 51545.0_dk
  end function calculate_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_mean_long( time ) result( mean_long )
    ! Calculate the mean longitude of the sun at a given time

    !> Arguments
    real(dk), intent(in) :: time

    !> Local variables
    real(dk), parameter :: THREE60 = 360._dk
    real(dk), parameter :: PI    = 2._dk*ASIN( 1._dk )
    real(dk), parameter :: degrees_to_radians   = PI/180._dk
    real(dk) :: mean_long ! Mean longitude of Sun, corrected for aberration (deg; normalized to 0-360)

    !                    ** force mean longitude between 0 and 360 degs
    mean_long = 280.460_dk + 0.9856474_dk*time
    mean_long = MOD( mean_long, THREE60 )
    IF( mean_long < 0._dk ) mean_long = mean_long + THREE60

  end function calculate_mean_long

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_mean_anomaly( time ) result( mean_anomaly )
    ! Calculate the mean anamoly at a given time

    !> Arguments
    real(dk), intent(in) :: time

    !> Local variables
    real(dk), parameter :: THREE60 = 360._dk
    real(dk), parameter :: PI    = 2._dk*ASIN( 1._dk )
    real(dk), parameter :: degrees_to_radians   = PI/180._dk
    real(dk)  :: mean_anomaly ! radians, normalized to 0 to 2*pi

    !                    ** mean anomaly in radians between 0 and 2*pi
    mean_anomaly = 357.528_dk + 0.9856003_dk*time
    mean_anomaly = MOD( mean_anomaly, THREE60 )
    IF( mean_anomaly < 0._dk ) mean_anomaly = mean_anomaly + THREE60

    mean_anomaly = mean_anomaly * degrees_to_radians

  end function calculate_mean_anomaly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_declination( ecliptic_longitude, obliquity ) result (dec)
    ! Calculate the declimation angle

    !> arguments
    real(dk), intent(in) :: ecliptic_longitude, & ! the ecliptic longitude
                            obliquity             ! the obliquity

    !> local variables
    real(dk) :: dec ! declination (radians)
    dec  = asin( sin( obliquity )*sin( ecliptic_longitude ) )

  end function calculate_declination

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_obliquity( time ) result ( oblqec )
    ! Calculate the obliquity

    !> arguments
    real(dk), intent(in) :: time

    !> local variables
    real(dk) :: oblqec ! Obliquity of the ecliptic (radians)
    real(dk), parameter :: pi    = 2._dk*asin( 1._dk )
    real(dk), parameter :: degrees_to_radians   = pi/180._dk

    oblqec = 23.439_dk - 0.0000004_dk*time
    oblqec = oblqec * degrees_to_radians
  end function calculate_obliquity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_ecliptic_longitude( mean_anomaly, mean_long ) &
      result( eclong )
      ! Calculate the ecliptic longitude

    !> arguments
    real(dk), intent(in)    ::  mean_anomaly, mean_long

    !> local variables
    real(dk), parameter :: three60 = 360._dk
    real(dk), parameter :: rzero   = 0._dk
    real(dk), parameter :: rtwo    = 2._dk
    real(dk), parameter :: pi    = 2._dk*asin( 1._dk )
    real(dk), parameter :: degrees_to_radians   = pi/180._dk

    real(dk) :: eclong ! Ecliptic longitude (radians)

    eclong = mean_long + 1.915_dk*sin( mean_anomaly ) &
                    + 0.020_dk*sin( rtwo*mean_anomaly )
    eclong = mod( eclong, three60 )
    if( eclong < rzero ) eclong = eclong + three60

    eclong = eclong * degrees_to_radians

  end function calculate_ecliptic_longitude

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_right_ascension( ecliptic_longitude, obliquity ) &
      result ( ra )
      ! Calculate the right ascension

    !> Arguments
    real(dk), intent(in) :: ecliptic_longitude, obliquity

    !> Local variables
    real(dk), parameter :: rzero   = 0._dk
    real(dk), parameter :: pi    = 2._dk*asin( 1._dk )
    real(dk), parameter :: twopi = 2._dk*pi
    real(dk) :: den, num
    real(dk) :: ra ! Right ascension  (radians)

    !                    ** right ascension
    num  = cos( obliquity )*sin( ecliptic_longitude )
    den  = cos( ecliptic_longitude )
    ra   = atan( num / den )

    !                    ** force right ascension between 0 and 2*pi
    if( den < rzero ) then
        ra  = ra + pi
    else if( num < rzero ) then
        ra  = ra + twopi
    end if
  end function calculate_right_ascension

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate_hour_angle( time, hour, long, right_ascension ) &
      result( ha )
      ! Calculate the hour angle

    !> arguments
    real(dk), intent(in)    ::  time, hour, long, right_ascension

    !> local variables
    real(dk), parameter :: day2hrs = 24._dk
    real(dk), parameter :: rzero   = 0._dk
    real(dk), parameter :: pi    = 2._dk*asin( 1._dk )
    real(dk), parameter :: twopi = 2._dk*pi
    real(dk), parameter :: d2r   = pi/180._dk

    ! lmst      Local mean sidereal time (radians)
    real(dk) :: gmst, lmst
    real(dk) :: ha ! Hour angle (radians, -pi to pi)

    !                    ** greenwich mean sidereal time in hours
    gmst = 6.697375_dk + 0.0657098242_dk*time + hour
    !                    ** hour not changed to sidereal time since
    !                    ** 'time' includes the fractional day

    gmst  = mod( gmst, day2hrs)
    if( gmst < rzero ) gmst   = gmst + day2hrs

    !                    ** local mean sidereal time in radians
    lmst  = gmst + long / 15._dk
    lmst  = mod( lmst, day2hrs )
    if( lmst < rzero ) lmst   = lmst + day2hrs

    lmst   = lmst*15._dk*d2r

    !                    ** hour angle in radians between -pi and pi
    ha  = lmst - right_ascension

    if( ha < -pi ) then
      ha  = ha + twopi
    elseif( ha > pi )  then
      ha  = ha - twopi
    endif

  end function calculate_hour_angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function is_leap( year ) result( leap )
    ! Determines if a function is a `leap year <https://en.wikipedia.org/wiki/Leap_year#Algorithm>`_

    !> Arguments
    integer(ik), intent(in) ::  year ! integer year

    !> Local variables
    logical(lk) :: leap

    leap = .false.

    if( mod( year, 4_ik ) == 0 .and. &
      ( mod( year, 100 ) /= 0_ik .or. mod( year, 400_ik ) == 0) &
    ) then
      leap = .true.
    endif

  end function is_leap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_utils
