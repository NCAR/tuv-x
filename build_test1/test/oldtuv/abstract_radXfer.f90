
   module abstract_radXfer

   implicit none

   private
   public :: abstract_radXfer_t

   type, abstract :: abstract_radXfer_t
     contains
     procedure(initial),   deferred :: initialize
     procedure(run),       deferred :: calculate
   end type abstract_radXfer_t

   interface
      subroutine initial( this, nlyr, nstr, zen, nid, dsdh, &
                          dtrl, dto3, dto2, dtso2, dtno2, &
                          dtcld, omcld, gcld, &
                          dtaer, omaer, gaer, &
                          dtsnw, omsnw, gsnw, &
                          dtany, omany, gany )

      import abstract_radXfer_t
      class(abstract_radXfer_t), intent(inout) :: this
      integer, intent(in)  :: nlyr, nstr
      integer, intent(in)  :: nid(0:)
      real, intent(in)     :: zen
      real, intent(in)     :: dsdh(0:,:)
      real, intent(in)     :: dtrl(:), dto3(:), dto2(:), dtso2(:), dtno2(:)
      real, intent(in)     :: dtcld(:), omcld(:), gcld(:)
      real, intent(in)     :: dtaer(:), omaer(:), gaer(:)
      real, intent(in)     :: dtsnw(:), omsnw(:), gsnw(:)
      real, intent(in)     :: dtany(:), omany(:), gany(:)

      end subroutine initial

      subroutine run( this, nlyr, nstr, albedo, &
                            fdr, fup, fdn, edr, eup, edn )
        import abstract_radXfer_t
        class(abstract_radXfer_t), intent(inout) :: this
        integer, intent(in)  :: nlyr, nstr
	real, intent(in)     :: albedo
        real, intent(out)    :: fdr(:), fup(:), fdn(:)
        real, intent(out)    :: edr(:), eup(:), edn(:)
      end subroutine run

   end interface

   end module abstract_radXfer
