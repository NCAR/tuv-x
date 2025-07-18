
   module disord

   use abstract_radXfer, only : abstract_radXfer_t

   implicit none

   private
   public :: disord_t

   type, extends(abstract_radXfer_t) :: disord_t
     real                 :: umu0
     integer, allocatable :: nid(:)
     real, allocatable    :: dsdh(:,:)
     real, allocatable    :: dtauc(:), ssalb(:)
     real, allocatable    :: pmom(:,:)
     contains
     procedure :: initialize
     procedure :: calculate
     final     :: finalize
   end type disord_t

   REAL, PARAMETER :: rZERO = 0.0
   REAL, PARAMETER :: rONE  = 1.0

   contains

   subroutine initialize( this, nlyr, nstr, zen, nid, dsdh, &
                          dtrl, dto3, dto2, dtso2, dtno2, &
                          dtcld, omcld, gcld, &
                          dtaer, omaer, gaer, &
                          dtsnw, omsnw, gsnw, &
                          dtany, omany, gany )

   use tuv_params, only : pi, largest, precis

   class(disord_t), intent(inout) :: this
   integer, intent(in)  :: nlyr, nstr
   integer, intent(in)  :: nid(0:)
   real, intent(in)     :: zen
   real, intent(in)     :: dsdh(0:,:)
   real, intent(in)     :: dtrl(:), dto3(:), dto2(:), dtso2(:), dtno2(:)
   real, intent(in)     :: dtcld(:), omcld(:), gcld(:)
   real, intent(in)     :: dtaer(:), omaer(:), gaer(:)
   real, intent(in)     :: dtsnw(:), omsnw(:), gsnw(:)
   real, intent(in)     :: dtany(:), omany(:), gany(:)

   real, parameter :: dr = pi/180.
   real, parameter :: floor = rONE/largest

   integer :: strNdx
   real :: dscld(nlyr), dacld(nlyr)
   real :: dsaer(nlyr), daaer(nlyr)
   real :: dssnw(nlyr), dasnw(nlyr)
   real :: dsany(nlyr), daany(nlyr)
   real :: dtsct(nlyr), dtabs(nlyr)
   real :: om(nlyr)
   real :: pmcld(nlyr), pmaer(nlyr), pmsnw(nlyr), pmany(nlyr), pmray(nlyr)

   allocate( this%dtauc(nlyr), this%ssalb(nlyr) )
   allocate( this%pmom(0:nstr,nlyr) )

   this%umu0 = cos( zen*dr )
   this%nid  = nid
   this%dsdh = dsdh

   dscld = dtcld * omcld
   dacld = dtcld * (rONE - omcld)

   dsaer = dtaer * omaer
   daaer = dtaer * (rONE - omaer)

   dssnw = dtsnw * omsnw
   dasnw = dtsnw * (rONE - omsnw)

   dsany = dtany * omany
   daany = dtany * (rONE - omany)

   dtsct = max( dtrl + dscld + dsaer + dssnw + dsany,floor )
   dtabs = max( dto2 + dto3 + dtso2 + dtno2 + dacld + daaer + dasnw + daany,floor )

   this%dtauc(nlyr:1:-1) = max( dtsct + dtabs,precis )
   om(nlyr:1:-1) = dtsct/(dtsct + dtabs)
   where( dtsct == floor )
     om(nlyr:1:-1) = floor
   endwhere
   this%ssalb(nlyr:1:-1) = max( min( om(nlyr:1:-1),rONE - precis ),precis )

   this%pmom(0,:) = rONE
   do strNdx = 1,nstr
     pmcld = gcld**strNdx
     pmaer = gaer**strNdx
     pmsnw = gsnw**strNdx
     pmany = gany**strNdx
     if( strNdx == 2 ) then
       pmray = 0.1
     else
       pmray = rZERO
     endif
     this%pmom( strNdx, nlyr:1:-1 ) = &
       (pmcld*dscld + pmaer*dsaer + pmsnw*dssnw + pmany*dsany + pmray*dtrl)/dtsct
   enddo

   end subroutine initialize

   subroutine calculate( this, nlyr, nstr, albedo, &
                         fdr, fup, fdn, edr, eup, edn )

   use tuv_params,  only : pi
   use DISORD_SUBS, only : PSNDO

   class(disord_t), intent(inout) :: this

   real, parameter     :: fourPi = 4. * pi
   INTEGER, intent(in) :: nlyr, nstr
   REAL, intent(in)    :: albedo
   REAL, intent(out)   :: eup(:), edn(:), edr(:)
   REAL, intent(out)   :: fup(:), fdn(:), fdr(:)

   CALL PSNDO( this%dsdh, this%nid, &
               NLYR, this%DTAUC, this%SSALB, this%PMOM, &
               ALBEDO, NSTR, this%umu0, &
               edr, edn, eup, fdr, fup, fdn )

   end subroutine calculate

   subroutine finalize( this )

   type(disord_t), intent(inout) :: this

   if( allocated( this%nid ) ) then
     deallocate( this%nid )
   endif
   if( allocated( this%dsdh ) ) then
     deallocate( this%dsdh )
   endif
   if( allocated( this%dtauc ) ) then
     deallocate( this%dtauc )
   endif
   if( allocated( this%ssalb ) ) then
     deallocate( this%ssalb )
   endif
   if( allocated( this%pmom ) ) then
     deallocate( this%pmom )
   endif

   end subroutine finalize

   end module disord
