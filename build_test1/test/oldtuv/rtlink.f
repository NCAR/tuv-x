      module rtlink_mod

      IMPLICIT NONE

      private
      public :: rtlink

      contains

      SUBROUTINE rtlink(nstr,
     $     albedo, zen,
     $     dsdh, nid,
     $     dtrl, 
     $     dto3, 
     $     dto2,
     $     dtso2,
     $     dtno2, 
     $     dtcld, omcld, gcld,
     $     dtaer, omaer, gaer,
     $     dtsnw, omsnw, gsnw,
     $     dt_any,om_any,g_any,
     $     edir, edn, eup, fdir, fdn, fup)

      use tuv_params,       only : pi
      use abstract_radXfer, only : abstract_radXfer_t
      use delta_eddington,  only : delta_eddington_t
      use disord,           only : disord_t

      IMPLICIT NONE

* input

      INTEGER, intent(in) :: nstr
      INTEGER, intent(in) :: nid(0:)
      REAL, intent(in) :: albedo
      REAL, intent(in) :: zen
      REAL, intent(in) :: dsdh(0:,:)
      REAL, intent(in) :: dtrl(:)
      REAL, intent(in) :: dto3(:), dto2(:)
      REAL, intent(in) :: dtso2(:), dtno2(:)
      REAL, intent(in) :: dtcld(:), omcld(:), gcld(:)
      REAL, intent(in) :: dtaer(:), omaer(:), gaer(:)
      REAL, intent(in) :: dtsnw(:), omsnw(:), gsnw(:)
      REAL, intent(in) :: dt_any(:), om_any(:), g_any(:)

* output

      REAL, intent(out) :: edir(:), edn(:), eup(:)
      REAL, intent(out) :: fdir(:), fdn(:), fup(:)

* constants:

      REAL, parameter    :: rZERO = 0.0
      REAL, parameter    :: dr = pi/180.
      REAL, parameter    :: fourPI = 4. * pi

* local:

      INTEGER :: nlyr, nz
      class(abstract_radXfer_t), allocatable :: radiative_xfer_obj

*_______________________________________________________________________

* initialize:

      fdir = rZERO
      fup  = rZERO
      fdn  = rZERO
      edir = rZERO
      eup  = rZERO
      edn  = rZERO

      nz   = size(fdir)
      nlyr = nz - 1
      if( nstr < 2 ) then
        allocate( delta_eddington_t :: radiative_xfer_obj )
      else
        allocate( disord_t :: radiative_xfer_obj )
      endif
* call rt routines
* initialize
      call radiative_xfer_obj%initialize( 
     $       nlyr, nstr, zen, nid, dsdh,
     $       dtrl, dto3, dto2, dtso2, dtno2,
     $       dtcld, omcld, gcld,
     $       dtaer, omaer, gaer,
     $       dtsnw, omsnw, gsnw,
     $       dt_any, om_any, g_any )
* calculate radiative field
      call radiative_xfer_obj%calculate( 
     $        nlyr, nstr, albedo,
     $        fdir, fup, fdn, edir, eup, edn)

* output (top-down -> bottom-up)
      fdir(1:nz) = fdir(nz:1:-1)
      fup(1:nz)  = fup(nz:1:-1)
      fdn(1:nz)  = fdn(nz:1:-1)
      edir(1:nz) = edir(nz:1:-1)
      eup(1:nz)  = eup(nz:1:-1)
      edn(1:nz)  = edn(nz:1:-1)

      if( nstr > 1 ) then
        fdir = fourPI * fdir
        fup  = fourPI * fup
        fdn  = fourPI * fdn
      endif

      END SUBROUTINE rtlink

      end module rtlink_mod
