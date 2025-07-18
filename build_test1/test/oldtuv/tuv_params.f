      MODULE TUV_PARAMS

      IMPLICIT NONE

* BROADLY USED PARAMETERS:
*_________________________________________________
* i/o file unit numbers
      INTEGER, PARAMETER :: kout = 53, kin = 12
*_________________________________________________
* altitude, wavelength, time (or solar zenith angle) grids
      INTEGER, PARAMETER :: kz = 125, kw = 1000, kt = 100
*_________________________________________________
* number of weighting functions
      INTEGER, PARAMETER :: ks = 60, kj = 150, kdom = 200
* delta for adding points at beginning or end of data grids
      REAL, PARAMETER :: deltax = 1.e-5

* some constants...

* pi:
      REAL, PARAMETER :: pi = 3.1415926535898

* radius of the earth, km:
      REAL, PARAMETER :: radius = 6.371E+3

* Planck constant x speed of light, J m

      REAL, PARAMETER :: hc = 6.626068E-34 * 2.99792458E8

* largest number of the machine:
      REAL, PARAMETER :: largest = 1.E+36

* small numbers (positive and negative)
      REAL, PARAMETER :: pzero = +10./largest
      REAL, PARAMETER :: nzero = -10./largest

* machine precision
      REAL, PARAMETER :: precis = 1.e-7

* More physical constants:
*_________________________________________________________________
* Na = 6.022142E23  mol-1       = Avogadro constant
* kb = 1.38065E-23  J K-1       = Boltzmann constant 
* R  = 8.31447      J mol-1 K-1 = molar gas constant
* h  = 6.626068E-34 J s         = Planck constant 
* c  = 2.99792458E8 m s-1       = speed of light in vacuum 
* G  = 6.673E-11    m3 kg-1 s-2 = Netwonian constant of gravitation
* sb = 5.67040E-8   W m-2 K-4   = Stefan-Boltzmann constant
*_________________________________________________________________
* (1) From NIST Reference on Constants, Units, and Uncertainty
* http://physics.nist.gov/cuu/index.html Oct. 2001.
* (2) These constants are not assigned to variable names;  in other 
* words this is not Fortran code, but only a text table for quick 
* reference.  To use, you must declare a variable name/type and
* assign the value to that variable. Or assign as parameter (see
* example for pi above).

      END MODULE TUV_PARAMS
