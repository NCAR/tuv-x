
   module delta_eddington

   use abstract_radXfer, only : abstract_radXfer_t

   implicit none

   private
   public :: delta_eddington_t

   logical :: initialized = .true.

   type, extends(abstract_radXfer_t) :: delta_eddington_t
     real                 :: umu0
     integer, allocatable :: nid(:)
     real, allocatable    :: dsdh(:,:)
     real, allocatable    :: dt(:), om(:), g(:)
   contains
     procedure :: initialize
     procedure :: calculate
     final     :: finalize
   end type delta_eddington_t

   REAL(8), PARAMETER :: largest = 1.e36_8
   REAL(8), PARAMETER :: rZERO  = 0.0_8
   REAL(8), PARAMETER :: rONE   = 1.0_8
   REAL(8), PARAMETER :: rTWO   = 2.0_8
   REAL(8), PARAMETER :: rTHREE = 3.0_8
   REAL(8), PARAMETER :: rFOUR  = 4.0_8
   REAL(8), PARAMETER :: rSEVEN = 7.0_8

   contains

   subroutine initialize( this, nlyr, nstr, zen, nid, dsdh, &
                          dtrl, dto3, dto2, dtso2, dtno2, &
                          dtcld, omcld, gcld, &
                          dtaer, omaer, gaer, &
                          dtsnw, omsnw, gsnw, &
                          dtany, omany, gany )

   use tuv_params, only : pi

   class(delta_eddington_t), intent(inout) :: this
   integer, intent(in)  :: nlyr, nstr
   integer, intent(in)  :: nid(0:)
   real, intent(in)     :: zen
   real, intent(in)     :: dsdh(0:,:)
   real, intent(in)     :: dtrl(:), dto3(:), dto2(:), dtso2(:), dtno2(:)
   real, intent(in)     :: dtcld(:), omcld(:), gcld(:)
   real, intent(in)     :: dtaer(:), omaer(:), gaer(:)
   real, intent(in)     :: dtsnw(:), omsnw(:), gsnw(:)
   real, intent(in)     :: dtany(:), omany(:), gany(:)

   real(8), parameter :: dr = pi/180._8
   real(8), parameter :: floor = rONE/largest

   real :: dscld(nlyr), dacld(nlyr)
   real :: dsaer(nlyr), daaer(nlyr)
   real :: dssnw(nlyr), dasnw(nlyr)
   real :: dsany(nlyr), daany(nlyr)
   real :: dtsct(nlyr), dtabs(nlyr)

   allocate( this%nid(0:nlyr) )
   allocate( this%dsdh(0:nlyr,nlyr) )
   allocate( this%dt(nlyr), this%om(nlyr), this%g(nlyr) )

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

   !> reorder from bottum-up to top-down
   this%dt(nlyr:1:-1) = dtsct + dtabs
   this%om(nlyr:1:-1) = dtsct/(dtsct + dtabs)
   where( dtsct == floor )
     this%om(nlyr:1:-1) = floor
   endwhere

   this%g(nlyr:1:-1) = (gcld*dscld + gaer*dsaer + gsnw*dssnw + gany*dsany)/dtsct

   end subroutine initialize

   subroutine calculate( this, nlyr, nstr, albedo, &
                         fdr, fup, fdn, edr, eup, edn )

   use tuv_params, only : precis
   use linalgebra, only : linalgebra_t
   use debug,      only : diagout

   class(delta_eddington_t), intent(inout) :: this

   INTEGER, intent(in) :: nlyr, nstr
   REAL, intent(in)    :: albedo
   REAL, intent(out)   :: fup(:), fdn(:), fdr(:)
   REAL, intent(out)   :: eup(:), edn(:),edr(:)


!******
! local:
!******
      REAL :: sum
      REAL :: tausla(0:nlyr), tauc(0:nlyr)
      REAL :: mu2(0:nlyr)

! internal coefficients and matrix
      INTEGER :: row
      REAL    :: lam(nlyr),taun(nlyr),bgam(nlyr)
      REAL    :: e1(nlyr),e2(nlyr),e3(nlyr),e4(nlyr)
      REAL    :: cup(nlyr),cdn(nlyr)
      REAL    :: cuptn(nlyr),cdntn(nlyr)
      REAL    :: mu1(nlyr)
      REAL    :: a(2*nlyr),b(2*nlyr),d(2*nlyr)
      REAL    :: e(2*nlyr),y(2*nlyr)

!******
! other:
!******
      REAL :: pifs, fdn0, surfem, tempg
      REAL :: f, g, om
      REAL :: gam1, gam2, gam3, gam4
      REAL :: gi(nlyr), omi(nlyr)

      INTEGER :: mrows, lev
      INTEGER :: i, j
      REAL :: expon, expon0, expon1, divisr, temp, up, dn
      REAL :: ssfc

      TYPE(linalgebra_t) :: linpack

! Some additional program constants:

      REAL(8), parameter :: eps = 1.E-3_8
!_______________________________________________________________________

! MU = cosine of solar zenith angle
! RSFC = surface albedo
! TAUU =  unscaled optical depth of each layer
! OMU  =  unscaled single scattering albedo
! GU   =  unscaled asymmetry factor
! NLAYER = number of layers in the atmosphere
! NLEVEL = nlayer + 1 = number of levels

! initial conditions:  pi*solar flux = 1;  diffuse incidence = 0
      pifs = rONE
      fdn0 = rZERO

! emission at surface (for night light pollution, set pifs = 0, surfem = 1.)

      surfem = rZERO

!************* compute coefficients for each layer:

! GAM1 - GAM4 = 2-stream coefficients, different for different approximations
! EXPON0 = calculation of e when TAU is zero
! EXPON1 = calculation of e when TAU is TAUN
! CUP and CDN = calculation when TAU is zero
! CUPTN and CDNTN = calc. when TAU is TAUN
! DIVISR = prevents division by zero

      tauc   = rZERO
      tausla = rZERO
      mu2    = rONE/SQRT(largest)

! delta-scaling. Has to be done for delta-Eddington approximation, 
! delta discrete ordinate, Practical Improved Flux Method, delta function,
! and Hybrid modified Eddington-delta function methods approximations

   associate( mu => this%umu0, rsfc => albedo, nid => this%nid, dsdh => this%dsdh, &
              tauu => this%dt, gu => this%g, omu => this%om )

       DO i = 1, nlyr
         f     = gu(i)*gu(i)
         gi(i) = (gu(i) - f)/(rONE - f)
         omi(i) = (rONE - f)*omu(i)/(rONE - omu(i)*f)       
         taun(i) = (rONE - omu(i)*f)*tauu(i)
       ENDDO
       if( initialized ) then
         call diagout( 'tauu.old',tauu )
         call diagout( 'gu.old',gu )
         call diagout( 'omu.old',omu )
         call diagout( 'taun.old',taun )
       endif

! calculate slant optical depth at the top of the atmosphere when zen>90.
! in this case, higher altitude of the top layer is recommended which can 
! be easily changed in gridz.f.

         IF(mu < rZERO) THEN
           IF(nid(0) < 0) THEN
             tausla(0) = largest
           ELSE
             sum = rZERO
             DO j = 1, nid(0)
              sum = sum + rTWO*taun(j)*dsdh(0,j)
             END DO
             tausla(0) = sum 
           END IF
         END IF

        layer_loop: DO i = 1, nlyr

          g  = gi(i)
          om = omi(i)
          tauc(i) = tauc(i-1) + taun(i)

! stay away from 1 by precision.  For g, also stay away from -1

          tempg = MIN(abs(g),rONE - precis)
          g = SIGN(tempg,g)
          om = MIN(om,rONE-precis)


! calculate slant optical depth

          IF(nid(i) < 0) THEN
            tausla(i) = largest
          ELSE
            sum = rZERO
            DO j = 1, MIN(nid(i),i)
               sum = sum + taun(j)*dsdh(i,j)
            ENDDO
            DO j = MIN(nid(i),i)+1,nid(i)
               sum = sum + rTWO*taun(j)*dsdh(i,j)
            ENDDO
            tausla(i) = sum 
            IF(tausla(i) == tausla(i-1)) THEN
              mu2(i) = SQRT(largest)
            ELSE
              mu2(i) = (tauc(i)-tauc(i-1))/(tausla(i)-tausla(i-1))
              mu2(i) = SIGN( MAX(ABS(mu2(i)),rONE/SQRT(largest)),mu2(i) )
            END IF
    
            if( initialized .and. i == 2 ) then
              write(*,*)'TUV: dsdh diagnostic'
              write(*,*) dsdh(i,1:i)
            endif
          END IF

!** the following gamma equations are from pg 16,289, Table 1
!** save mu1 for each approx. for use in converting irradiance to actinic flux

! Eddington approximation(Joseph et al., 1976, JAS, 33, 2452):

        gam1 =  (rSEVEN - om*(rFOUR + rTHREE*g))/rFOUR
        gam2 = -(rONE - om*(rFOUR - rTHREE*g))/rFOUR
        gam3 = (rTWO - rTHREE*g*mu)/rFOUR
        gam4 = rONE - gam3
        mu1(i) = 0.5_8

         lam(i) = sqrt(gam1*gam1 - gam2*gam2)

         IF( gam2 /= rZERO) THEN
            bgam(i) = (gam1 - lam(i))/gam2
         ELSE
            bgam(i) = rZERO
         ENDIF

         expon = EXP(-lam(i)*taun(i))

! e1 - e4 = pg 16,292 equation 44
         
         e1(i) = rONE + bgam(i)*expon
         e2(i) = rONE - bgam(i)*expon
         e3(i) = bgam(i) + expon
         e4(i) = bgam(i) - expon

! the following sets up for the C equations 23, and 24
! found on page 16,290
! prevent division by zero (if LAMBDA=1/MU, shift 1/MU^2 by EPS = 1.E-3
! which is approx equiv to shifting MU by 0.5*EPS* (MU)**3

         expon0 = EXP(-tausla(i-1))
         expon1 = EXP(-tausla(i))
          
         divisr = lam(i)*lam(i) - rONE/(mu2(i)*mu2(i))
         temp   = MAX(eps,abs(divisr))
         divisr = SIGN(temp,divisr)

         up = om*pifs*((gam1 - rONE/mu2(i))*gam3 + gam4*gam2)/divisr
         dn = om*pifs*((gam1 + rONE/mu2(i))*gam4 + gam2*gam3)/divisr
         
! cup and cdn are when tau is equal to zero
! cuptn and cdntn are when tau is equal to taun

         cup(i) = up*expon0
         cdn(i) = dn*expon0
         cuptn(i) = up*expon1
         cdntn(i) = dn*expon1

         if( initialized .and. i == 3 ) then
           write(*,*) 'TUV: cup diagnostic'
           write(*,*) expon, expon0, expon1, divisr, temp, up, dn
           write(*,*) lam(i), mu2(i), gam1, gam2, gam3, gam4
           write(*,*) tauc(i-1:i), tausla(i-1:i)
         endif
 
      ENDDO layer_loop

       if( initialized ) then
         call diagout( 'e1.old',e1 )
         call diagout( 'e2.old',e2 )
         call diagout( 'e3.old',e3 )
         call diagout( 'e4.old',e4 )
         call diagout( 'cup.old',cup )
         call diagout( 'cdn.old',cdn )
         call diagout( 'cuptn.old',cuptn )
         call diagout( 'cdntn.old',cdntn )
         call diagout( 'lam.old',lam )
         call diagout( 'tausla.old',tausla )
         call diagout( 'mu2.old',mu2 )
       endif

!**************** set up matrix ******
! ssfc = pg 16,292 equation 37  where pi Fs is one (unity).

      ssfc = rsfc*mu*EXP(-tausla(nlyr))*pifs + surfem

! MROWS = the number of rows in the matrix

      mrows = 2*nlyr     
      
! the following are from pg 16,292  equations 39 - 43.
! set up first row of matrix:

      a(1) = rZERO
      b(1) = e1(1)
      d(1) = -e2(1)
      e(1) = fdn0 - cdn(1)

! set up odd rows 3 thru (MROWS - 1):

      i = 0
      DO row = 3, mrows - 1, 2
         i = i + 1
         a(row) = e2(i)*e3(i) - e4(i)*e1(i)
         b(row) = e1(i)*e1(i + 1) - e3(i)*e3(i + 1)
         d(row) = e3(i)*e4(i + 1) - e1(i)*e2(i + 1)
         e(row) = e3(i)*(cup(i + 1) - cuptn(i)) &
                + e1(i)*(cdntn(i) - cdn(i + 1))
      ENDDO

! set up even rows 2 thru (MROWS - 2): 

      i = 0
      DO row = 2, mrows - 2, 2
         i = i + 1
         a(row) = e2(i + 1)*e1(i) - e3(i)*e4(i + 1)
         b(row) = e2(i)*e2(i + 1) - e4(i)*e4(i + 1)
         d(row) = e1(i + 1)*e4(i + 1) - e2(i + 1)*e3(i + 1)
         e(row) = (cup(i + 1) - cuptn(i))*e2(i + 1) &
                - (cdn(i + 1) - cdntn(i))*e4(i + 1)
      ENDDO

! set up last row of matrix at MROWS:

      a(mrows) = e1(nlyr) - rsfc*e3(nlyr)
      b(mrows) = e2(nlyr) - rsfc*e4(nlyr)
      d(mrows) = rZERO
      e(mrows) = ssfc - cuptn(nlyr) + rsfc*cdntn(nlyr)

      if( initialized ) then
        call diagout( 'a.old',a )
        call diagout( 'b.old',b )
        call diagout( 'd.old',d )
        call diagout( 'e.old',e )
        write(*,*) 'e diagnostic'
        write(*,*) e(5), e1(2), e3(2), cup(3), cdn(3), cuptn(2), cdntn(2)
        initialized = .false.
      endif
! solve tri-diagonal system:

      y = linpack%tridiag(a, b, d, e)

!*** unfold solution of matrix, compute output fluxes:
      
! the following equations are from pg 16,291  equations 31 & 32

      fdr(1) = pifs * EXP( -tausla(0) )
      edr(1) = mu * fdr(1)
      edn(1) = fdn0
      eup(1) =  y(1)*e3(1) - y(2)*e4(1) + cup(1)
      fdn(1) = edn(1)/mu1(1)
      fup(1) = eup(1)/mu1(1)

      j   = 1
      row = 1 
      DO lev = 2, nlyr + 1
         fdr(lev) = pifs * EXP(-tausla(lev-1))
         edr(lev) =  mu *fdr(lev)
         edn(lev) =  y(row)*e3(j) + y(row + 1)*e4(j) + cdntn(j)
         eup(lev) =  y(row)*e1(j) + y(row + 1)*e2(j) + cuptn(j)
         fdn(lev) = edn(lev)/mu1(j)
         fup(lev) = eup(lev)/mu1(j)

         row = row + 2
         j = j + 1
      ENDDO

   end associate

   end subroutine calculate

   subroutine finalize( this )

   type(delta_eddington_t), intent(inout) :: this

   if( allocated( this%nid ) ) then
     deallocate( this%nid )
   endif
   if( allocated( this%dsdh ) ) then
     deallocate( this%dsdh )
   endif
   if( allocated( this%dt ) ) then
     deallocate( this%dt )
   endif
   if( allocated( this%om ) ) then
     deallocate( this%om )
   endif
   if( allocated( this%g ) ) then
     deallocate( this%g )
   endif

   end subroutine finalize

   end module delta_eddington
