! This file contains the following derived types,
! related to interpolation of input data
!     inter1
!     inter2
!     inter3
!     inter4
!=============================================================================*

   module interpolation

   use musica_constants, only : ik => musica_ik, dk => musica_dk, lk => musica_lk

   implicit none

   integer(ik), parameter :: iONE = 1_ik
   real(ik),    parameter :: rZERO = 0.0_dk

   type, abstract :: abs_interpolator_t
     character(len=10) :: handle_
   contains
     procedure(interpolate), deferred :: interpolate
   end type abs_interpolator_t

   interface
     function interpolate( this, xtarget, xsrc,ysrc, FoldIn ) result( ytarget )

     use musica_constants, only : ik => musica_ik, dk => musica_dk

     import abs_interpolator_t

     class(abs_interpolator_t), intent(inout) :: this
     integer(ik), intent(in), optional        :: FoldIn
     real(dk), intent(in) :: xtarget(:)
     real(dk), intent(in) :: xsrc(:), ysrc(:)
     real(dk)             :: ytarget(size(xtarget))

     end function interpolate
   end interface

   type, extends(abs_interpolator_t) :: interp1_t
   contains
     procedure :: interpolate => inter1
   end type interp1_t

   type, extends(abs_interpolator_t) :: interp2_t
   contains
     procedure :: interpolate => inter2
   end type interp2_t

   type, extends(abs_interpolator_t) :: interp3_t
   contains
     procedure :: interpolate => inter3
   end type interp3_t

   type, extends(abs_interpolator_t) :: interp4_t
   contains
     procedure :: interpolate => inter4
   end type interp4_t

   contains

      FUNCTION inter1(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points, onto a discrete target  =*
!=  grid.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified      =*
!=  discrete target grid.  A typical example would be the re-gridding of a   =*
!=  given data set for the vertical temperature profile to match the speci-  =*
!=  fied altitude grid.                                                      =*
!=  Some caution should be used near the end points of the grids.  If the    =*
!=  input data set does not span the range of the target grid, the remaining =*
!=  points will be set to zero, as extrapolation is not permitted.           =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - integer(ik), number of points in the target grid                    (I)=*
!=  XG  - real(dk), target grid (e.g. altitude grid)                          (I)=*
!=  YG  - real(dk), y-data re-gridded onto XG                                 (O)=*
!=  N   - integer(ik), number of points in the input data set                 (I)=*
!=  X   - real(dk), grid on which input data are defined                      (I)=*
!=  Y   - real(dk), input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

! input:
      class(interp1_t), intent(inout) :: this
      integer(ik), intent(in), optional :: FoldIn
      real(dk), intent(in) :: xtarget(:)
      real(dk), intent(in) :: xsrc(:), ysrc(:)

! output:
      real(dk) :: ytarget(size(xtarget))

! local:
      integer(ik) :: n
      integer(ik) :: jsave, i, j
      real(dk)    :: slope
!_______________________________________________________________________

      n = size(xsrc)
      jsave = iONE
      ytarget = rZERO
      DO i = 1_ik,size(xtarget)
        j = jsave
        DO
          IF( xtarget(i) < xsrc(j) .or. xtarget(i) >= xsrc(j+iONE) ) THEN
            j = j+iONE
            IF (j >= n) THEN
              EXIT
            ENDIF
          ELSE
            slope = (ysrc(j+iONE)-ysrc(j)) / (xsrc(j+iONE)-xsrc(j))
            ytarget(i) = ysrc(j) + slope * (xtarget(i) - xsrc(j))
            jsave = j
            EXIT
          ENDIF
        ENDDO
      ENDDO

      END FUNCTION inter1

!=============================================================================*

      FUNCTION inter2(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on single, discrete points onto a set of target     =*
!=  bins.                                                                    =*
!=  The original input data are given on single, discrete points of an       =*
!=  arbitrary grid and are being linearly interpolated onto a specified set  =*
!=  of target bins.  In general, this is the case for most of the weighting  =*
!=  functions (action spectra, molecular cross section, and quantum yield    =*
!=  data), which have to be matched onto the specified wavelength intervals. =*
!=  The average value in each target bin is found by averaging the trapezoi- =*
!=  dal area underneath the input data curve (constructed by linearly connec-=*
!=  ting the discrete input values).                                         =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data set does not span the range of the target grid, an error      =*
!=  message is printed and the execution is stopped, as extrapolation of the =*
!=  data is not permitted.                                                   =*
!=  If the input data does not encompass the target grid, use ADDPNT to      =*
!=  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG  - integer(ik), number of bins + 1 in the target grid                  (I)=*
!=  XG  - real(dk), target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - real(dk), y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - integer(ik), number of points in input grid                         (I)=*
!=  X   - real(dk), grid on which input data are defined                      (I)=*
!=  Y   - real(dk), input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

      use musica_assert, only : die_msg

! input:
      class(interp2_t), intent(inout) :: this
      integer(ik), intent(in), optional :: FoldIn
      real(dk), intent(in)    :: xtarget(:)
      real(dk), intent(in)    :: xsrc(:), ysrc(:)

! output:
      real(dk)                :: ytarget(size(xtarget))

! local:
      integer(ik) :: ng, n
      integer(ik) :: ngintv
      integer(ik) :: i, k, jstart
      real(dk) :: area, xgl, xgu
      real(dk) :: darea, slope
      real(dk) :: a1, a2, b1, b2

      associate( xg => xtarget, yg => ytarget, x => xsrc, y => ysrc )

!  test for correct ordering of data, by increasing value of x
      if( any( x(1:n-1) >= x(2:n) ) ) then
        call die_msg( 1000,'src grid must be monotonically increasing' )
      endif

      if( any( xg(1:ng-1) >= xg(2:ng) ) ) then
        call die_msg( 1000,'target grid must be monotonically increasing' )
      endif

! check for xg-values outside the x-range
      IF( (x(1) > xg(1)) .OR. (x(n) < xg(ng)) ) THEN
        call die_msg( 1001,'src and target grid do not overlap' )
      ENDIF

      yg = rZERO
      n  = size(x)
      ng = size(xg)
!  find the integral of each grid interval and use this to 
!  calculate the average y value for the interval      
!  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1_ik
      ngintv = ng - 1_ik
      DO i = 1_ik,ngintv
! initalize:
            area = 0.0_dk
            xgl = xg(i)
            xgu = xg(i+1_ik)

!  discard data before the first grid interval and after the 
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval

            k = jstart
            IF (k < n) THEN
!  if both points are before the first grid, go to the next point
              DO
                IF( x(k+1) <= xgl ) THEN
                  jstart = k - 1_ik
                  k = k+1_ik
                  IF( k >= n) THEN
                    EXIT
                  ENDIF
                ELSE
                  EXIT
                ENDIF
              ENDDO
!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
              DO WHILE( k < n .AND. x(k) < xgu )
                jstart = k-1_ik
! compute x-coordinates of increment
                a1 = MAX(x(k),xgl)
                a2 = MIN(x(k+1_ik),xgu)
!  if points coincide, contribution is zero
                IF (x(k+1) == x(k)) THEN
                   darea = 0.0_dk
                ELSE
                   slope = (y(k+1_ik) - y(k))/(x(k+1_ik) - x(k))
                   b1 = y(k) + slope*(a1 - x(k))
                   b2 = y(k) + slope*(a2 - x(k))
                   darea = .5_dk*(a2 - a1)*(b2 + b1)
                ENDIF
!  find the area under the trapezoid from a1 to a2
                area = area + darea
! go to next point
                k = k+1
              ENDDO
            ENDIF

!  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)
      ENDDO

      end associate

      END FUNCTION inter2

!=============================================================================*

      FUNCTION inter3(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on a set of bins onto a different set of target     =*
!=  bins.                                                                    =*
!=  The input data are given on a set of bins (representing the integral     =*
!=  of the input quantity over the range of each bin) and are being matched  =*
!=  onto another set of bins (target grid).  A typical example would be an   =*
!=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
!=  vals, that has to be matched onto the working wavelength grid.           =*
!=  The resulting area in a given bin of the target grid is calculated by    =*
!=  simply adding all fractional areas of the input data that cover that     =*
!=  particular target bin.                                                   =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data do not span the full range of the target grid, the area in    =*
!=  the "missing" bins will be assumed to be zero.  If the input data extend =*
!=  beyond the upper limit of the target grid, the user has the option to    =*
!=  integrate the "overhang" data and fold the remaining area back into the  =*
!=  last target bin.  Using this option is recommended when re-gridding      =*
!=  vertical profiles that directly affect the total optical depth of the    =*
!=  model atmosphere.                                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG     - integer(ik), number of bins + 1 in the target grid               (I)=*
!=  XG     - real(dk), target grid (e.g. working wavelength grid);  bin i     (I)=*
!=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!=  YG     - real(dk), y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!=           y-value for bin i (i = 1..NG-1)                                 =*
!=  N      - integer(ik), number of bins + 1 in the input grid                (I)=*
!=  X      - real(dk), input grid (e.g. data wavelength grid);  bin i is      (I)=*
!=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!=  Y      - real(dk), input y-data on grid X;  Y(i) specifies the            (I)=*
!=           y-value for bin i (i = 1..N-1)                                  =*
!=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!=           FoldIn = 0 -> No folding of "overhang" data                     =*
!=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!=                         last target bin                                   =*
!-----------------------------------------------------------------------------*

      use musica_assert, only : die_msg

! input:
      class(interp3_t), intent(inout) :: this
      integer(ik), intent(in), optional :: FoldIn
      real(dk), intent(in)    :: xtarget(:)
      real(dk), intent(in)    :: xsrc(:), ysrc(:)

! output:
      real(dk)              :: ytarget(size(xtarget))

! local:
      integer(ik) :: nfrom, nto
      integer(ik) :: jstart, i, j, k, ntobins
      real(dk)    :: a1, a2, sum
      real(dk)    :: tail
      logical(lk) :: do_Foldin

! check whether flag given is legal
      if( present(FoldIn) ) then
        IF( FoldIn /= 0_ik .and. FoldIn /= 1_ik ) THEN
          call die_msg( 2000,'Foldin must be 0 or 1' )
        ENDIF
        do_FoldIn = FoldIn == 1_ik
      else
        call die_msg( 2001,'Foldin argument not present' )
      endif

      associate( xto => xtarget, yto => ytarget, xfrom => xsrc, yfrom => ysrc )

      nto   = size(xto)
      nfrom = size(xfrom)
! do interpolation

      jstart = iONE

      yto = rZERO
      ntobins = nto - 1
      DO  i = 1, ntobins
        sum = rZERO
        j = jstart

        IF (j < nfrom) THEN
          DO
            IF( xfrom(j+1) < xto(i) ) THEN
              jstart = j
              j = j+1
              IF( j >= nfrom ) THEN
                EXIT
              ENDIF
            ELSE
              EXIT
            ENDIF
          ENDDO

          DO WHILE( (xfrom(j) <= xto(i+1)) .and. (j < nfrom) )
            a1 = MAX(xfrom(j),xto(i))
            a2 = MIN(xfrom(j+1),xto(i+1))

            sum = sum + yfrom(j) * (a2-a1)/(xfrom(j+1)-xfrom(j))
            j = j+1
          ENDDO
          yto(i) = sum 
        ENDIF
      ENDDO
      

! if wanted, integrate data "overhang" and fold back into last bin
      IF (do_FoldIn ) THEN
         j = j-1
         a1 = xto(nto)           ! upper limit of last interpolated bin
         a2 = xfrom(j+1)         ! upper limit of last input bin considered
         
!        do folding only if grids don't match up and there is more input 
         IF( a2 > a1 .OR. j+1 < nfrom ) THEN
            tail = yfrom(j) * (a2-a1)/(xfrom(j+1)-xfrom(j))
            DO k = j+1, nfrom-1
               tail = tail + yfrom(k) * (xfrom(k+1)-xfrom(k))
            ENDDO
            yto(ntobins) = yto(ntobins) + tail
         ENDIF
      ENDIF

      end associate

      END FUNCTION inter3

!=============================================================================*

      FUNCTION inter4(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Map input data given on a set of bins onto a different set of target     =*
!=  bins.                                                                    =*
!=  The input data are given on a set of bins (representing the integral     =*
!=  of the input quantity over the range of each bin) and are being matched  =*
!=  onto another set of bins (target grid).  A typical example would be an   =*
!=  input data set spcifying the extra-terrestrial flux on wavelength inter- =*
!=  vals, that has to be matched onto the working wavelength grid.           =*
!=  The resulting area in a given bin of the target grid is calculated by    =*
!=  simply adding all fractional areas of the input data that cover that     =*
!=  particular target bin.                                                   =*
!=  Some caution should be used near the endpoints of the grids.  If the     =*
!=  input data do not span the full range of the target grid, the area in    =*
!=  the "missing" bins will be assumed to be zero.  If the input data extend =*
!=  beyond the upper limit of the target grid, the user has the option to    =*
!=  integrate the "overhang" data and fold the remaining area back into the  =*
!=  last target bin.  Using this option is recommended when re-gridding      =*
!=  vertical profiles that directly affect the total optical depth of the    =*
!=  model atmosphere.                                                        =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  NG     - integer(ik), number of bins + 1 in the target grid               (I)=*
!=  XG     - real(dk), target grid (e.g. working wavelength grid);  bin i     (I)=*
!=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!=  YG     - real(dk), y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!=           y-value for bin i (i = 1..NG-1)                                 =*
!=  N      - integer(ik), number of bins + 1 in the input grid                (I)=*
!=  X      - real(dk), input grid (e.g. data wavelength grid);  bin i is      (I)=*
!=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!=  Y      - real(dk), input y-data on grid X;  Y(i) specifies the            (I)=*
!=           y-value for bin i (i = 1..N-1)                                  =*
!=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!=           FoldIn = 0 -> No folding of "overhang" data                     =*
!=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!=                         last target bin                                   =*
!-----------------------------------------------------------------------------*

      use musica_assert, only : die_msg
      
! input:
      class(interp4_t), intent(inout) :: this
      integer(ik), intent(in), optional :: FoldIn
      real(dk), intent(in) :: xtarget(:)
      real(dk), intent(in) :: xsrc(:), ysrc(:)

! output:
      real(dk)             :: ytarget(size(xtarget))

! local:
      integer(ik) :: n, ng
      integer(ik) :: jstart, i, j, k
      real(dk)    :: a1, a2, sum
      real(dk)    :: tail
      logical(lk) :: do_FoldIn

! check whether flag given is legal
      if( present(FoldIn) ) then
        IF( FoldIn /= 0_ik .and. FoldIn /= 1_ik ) THEN
          call die_msg( 2000,'Foldin must be 0 or 1' )
          do_FoldIn = FoldIn == 1_ik
        ENDIF
      else
        call die_msg( 2001,'Foldin argument not present' )
      endif

      associate( xg => xtarget, yg => ytarget, x => xsrc, y => ysrc )

      n = size(x)
      ng = size(xg)
! do interpolation

      jstart = 1
      yg = rZERO

      DO i = 1, ng - 1
         sum = 0.
         j = jstart
         IF (j < n) THEN
           DO
             IF (x(j+1) .LT. xg(i)) THEN
               jstart = j
               j = j+1
               IF (j >= n) THEN
                 EXIT
               ENDIF
             ELSE
               EXIT
             ENDIF               
           ENDDO

           DO WHILE( x(j) <= xg(i+1) .and. j < n )
              a1 = MAX(x(j),xg(i))
              a2 = MIN(x(j+1),xg(i+1))
              sum = sum + y(j) * (a2-a1)
              j = j+1
           ENDDO
           yg(i) = sum /(xg(i+1)-xg(i))
        ENDIF
      ENDDO


! if wanted, integrate data "overhang" and fold back into last bin

      IF (do_FoldIn) THEN

         j = j-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(j+1)     ! upper limit of last input bin considered

!        do folding only if grids don't match up and there is more input 
         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
           tail = y(j) * (a2-a1)/(x(j+1)-x(j))
           DO k = j+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF

      end associate

      END FUNCTION inter4

   end module interpolation
