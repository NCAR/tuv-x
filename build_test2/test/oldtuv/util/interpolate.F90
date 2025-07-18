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
     real(dk), intent(in)  :: xtarget(:)
     real(dk), intent(in)  :: xsrc(:), ysrc(:)
     real(dk), allocatable :: ytarget(:)

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

      function inter1(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
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
      real(dk), allocatable :: ytarget(:)

! local:
      integer(ik) :: n
      integer(ik) :: jsave, i, j
      real(dk)    :: slope

      allocate( ytarget(size(xtarget)) )

      n = size(xsrc)
      jsave = iONE
      ytarget = rZERO
      DO i = iONE,size(xtarget)
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

      END function inter1

!=============================================================================*

      function inter2(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
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
      real(dk), allocatable   :: ytarget(:)

! local:
      integer(ik) :: ng, n
      integer(ik) :: ngintv
      integer(ik) :: i, k, jstart
      real(dk) :: area, xgl, xgu
      real(dk) :: darea, slope
      real(dk) :: a1, a2, b1, b2

      allocate( ytarget(size(xtarget)-1) )

!  test for correct ordering of data, by increasing value of x
      if( any( xsrc(1:n-1) >= xsrc(2:n) ) ) then
        call die_msg( 1000,'src grid must be monotonically increasing' )
      endif

      if( any( xtarget(1:ng-1) >= xtarget(2:ng) ) ) then
        call die_msg( 1000,'target grid must be monotonically increasing' )
      endif

! check for xg-values outside the x-range
      IF( (xsrc(1) > xtarget(1)) .OR. (xsrc(n) < xtarget(ng)) ) THEN
        call die_msg( 1001,'src and target grid do not overlap' )
      ENDIF

      ytarget = rZERO
      n  = size(xsrc)
      ng = size(xtarget)
!  find the integral of each grid interval and use this to 
!  calculate the average y value for the interval      
!  xgl and xgu are the lower and upper limits of the grid interval

      jstart = iONE
      ngintv = ng - iONE
      DO i = iONE,ngintv
! initalize:
        area = rZERO
        xgl = xtarget(i)
        xgu = xtarget(i+iONE)
!  discard data before the first grid interval and after the 
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval
        k = jstart
        IF (k < n) THEN
!  if both points are before the first grid, go to the next point
          DO
            IF( xsrc(k+iONE) <= xgl ) THEN
              jstart = k - iONE
              k = k + iONE
              IF( k >= n) THEN
                EXIT
              ENDIF
            ELSE
              EXIT
            ENDIF
          ENDDO
!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
          DO WHILE( k < n .AND. xsrc(k) < xgu )
            jstart = k - iONE
! compute x-coordinates of increment
            a1 = MAX(xsrc(k),xgl)
            a2 = MIN(xsrc(k+iONE),xgu)
!  if points coincide, contribution is zero
            IF (xsrc(k+iONE) == xsrc(k)) THEN
              darea = rZERO
            ELSE
              slope = (ysrc(k+iONE) - ysrc(k))/(xsrc(k+iONE) - xsrc(k))
              b1 = ysrc(k) + slope*(a1 - xsrc(k))
              b2 = ysrc(k) + slope*(a2 - xsrc(k))
              darea = .5_dk*(a2 - a1)*(b2 + b1)
            ENDIF
!  find the area under the trapezoid from a1 to a2
            area = area + darea
! go to next point
            k = k + iONE
          ENDDO
        ENDIF
!  calculate the average y after summing the areas in the interval
        ytarget(i) = area/(xgu - xgl)
      ENDDO

      END function inter2

!=============================================================================*

      function inter3(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
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
      real(dk), allocatable   :: ytarget(:)

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

      allocate( ytarget(size(xtarget)-1) )

      associate( xto => xtarget, yto => ytarget, xfrom => xsrc, yfrom => ysrc )

      nto   = size(xto)
      nfrom = size(xfrom)
! do interpolation

      jstart = iONE

      yto     = rZERO
      ntobins = nto - iONE
      DO  i = 1, ntobins
        sum = rZERO
        j = jstart

        IF (j < nfrom) THEN
          DO
            IF( xfrom(j+iONE) < xto(i) ) THEN
              jstart = j
              j = j+iONE
              IF( j >= nfrom ) THEN
                EXIT
              ENDIF
            ELSE
              EXIT
            ENDIF
          ENDDO

          DO WHILE( (xfrom(j) <= xto(i+iONE)) .and. (j < nfrom) )
            a1 = MAX(xfrom(j),xto(i))
            a2 = MIN(xfrom(j+iONE),xto(i+iONE))
            sum = sum + yfrom(j) * (a2-a1)/(xfrom(j+iONE)-xfrom(j))
            j = j+iONE
          ENDDO
          yto(i) = sum 
        ENDIF
      ENDDO
      

! integrate data "overhang" and fold back into last target bin
      IF (do_FoldIn ) THEN
         j = j-iONE
         a1 = xto(nto)           ! upper limit of last interpolated bin
         a2 = xfrom(j+iONE)      ! upper limit of last input bin considered
!        do folding only if grids don't match up and there is more input 
         IF( a2 > a1 .OR. j+iONE < nfrom ) THEN
            tail = yfrom(j) * (a2-a1)/(xfrom(j+iONE)-xfrom(j))
            DO k = j+iONE, nfrom-iONE
               tail = tail + yfrom(k) * (xfrom(k+iONE)-xfrom(k))
            ENDDO
            yto(ntobins) = yto(ntobins) + tail
         ENDIF
      ENDIF

      end associate

      END function inter3

!=============================================================================*

      function inter4(this, xtarget, xsrc,ysrc, FoldIn) result( ytarget )
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
      real(dk), intent(in)  :: xtarget(:)
      real(dk), intent(in)  :: xsrc(:), ysrc(:)

! output:
      real(dk), allocatable :: ytarget(:)

! local:
      integer(ik) :: n, ng
      integer(ik) :: jstart, i, j, k
      real(dk)    :: a1, a2, sum
      real(dk)    :: tail
      logical(lk) :: do_FoldIn

! check whether flag given is legal
      if( present(FoldIn) ) then
        IF( FoldIn /= 0_ik .and. FoldIn /= iONE ) THEN
          call die_msg( 2000,'Foldin must be 0 or 1' )
          do_FoldIn = FoldIn == 1_ik
        ENDIF
      else
        call die_msg( 2001,'Foldin argument not present' )
      endif

      allocate( ytarget(size(xtarget)-1) )

      n  = size(xsrc)
      ng = size(xtarget)
! do interpolation

      jstart  = iONE
      ytarget = rZERO

      DO i = 1, ng - 1
         sum = rZERO
         j = jstart
         IF (j < n) THEN
           DO
             IF (xsrc(j+iONE) < xtarget(i)) THEN
               jstart = j
               j = j+iONE
               IF (j >= n) THEN
                 EXIT
               ENDIF
             ELSE
               EXIT
             ENDIF               
           ENDDO

           DO WHILE( xsrc(j) <= xtarget(i+iONE) .and. j < n )
              a1 = MAX(xsrc(j),xtarget(i))
              a2 = MIN(xsrc(j+iONE),xtarget(i+1))
              sum = sum + ytarget(j) * (a2-a1)
              j = j+iONE
           ENDDO
           ytarget(i) = sum /(xtarget(i+iONE) - xtarget(i))
        ENDIF
      ENDDO
! integrate data "overhang" and fold back into last bin
      IF (do_FoldIn) THEN
         j = j - iONE
         a1 = xtarget(ng)     ! upper limit of last interpolated bin
         a2 = xsrc(j+iONE)    ! upper limit of last input bin considered

!        do folding only if grids don't match up and there is more input 
         IF ((a2 > a1) .OR. (j+iONE < n)) THEN
           tail = ysrc(j) * (a2 - a1)/(xsrc(j+iONE) - xsrc(j))
           DO k = j+iONE, n-iONE
              tail = tail + ysrc(k) * (xsrc(k+iONE) - xsrc(k))
           ENDDO
           ytarget(ng-iONE) = ytarget(ng-iONE) + tail
         ENDIF
      ENDIF

      END function inter4

   end module interpolation
