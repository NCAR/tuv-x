! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_constants
! General usage constants

  use musica_constants,            only : dk => musica_dk

  implicit none

  ! Numerical constants
  real(dk), parameter :: deltax = 1.e-5_dk         ! delta for adding points at beginning or end of data grids
  real(dk), parameter :: largest = 1.E+36_dk       ! largest number in the machine

  real(dk), parameter :: pzero =  &                ! small positive number
                           10._dk / largest
  real(dk), parameter :: nzero =  &                ! small negative number
                           -10._dk / largest
  real(dk), parameter :: precis = 1.e-7_dk         ! machine precision

  ! Physical constants
  real(dk), parameter :: pi = 3.1415926535898_dk   ! Pi
  real(dk), parameter :: radius = 6.371E+3_dk      ! Radius of the Earth [km]
  real(dk), parameter :: hc = 6.626068e-34_dk * 2.99792458e8_dk ! Plank's constants x speed of light [J m]

end module tuvx_constants
