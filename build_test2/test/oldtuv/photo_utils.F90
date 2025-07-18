
  module photo_utils

  use musica_constants, only : musica_ik, musica_dk
  use musica_assert,    only : die_msg

  implicit none

  private
  public :: inter2, inter4, addpnt

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter    :: rZERO = 0.0_musica_dk

  contains
  
    subroutine inter2(xto,yto,xfrom,yfrom,ierr)
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
!=  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!=  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!=        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!=  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!=        bin i (i = 1..NG-1)                                                =*
!=  N   - INTEGER, number of points in input grid                         (I)=*
!=  X   - REAL, grid on which input data are defined                      (I)=*
!=  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

      integer(musica_ik), intent(out) :: ierr
      real(musica_dk), intent(in)  :: xfrom(:), yfrom(:)
      real(musica_dk), intent(in)  :: xto(:)
      real(musica_dk), intent(out) :: yto(:)

! local:
      character(len=*), parameter :: Iam = 'inter2: '

      integer(musica_ik) :: nto, nfrom
      integer(musica_ik) :: ntom1, nfromm1
      integer(musica_ik) :: i, k
      real(musica_dk)    :: area, xtol, xtou
      real(musica_dk)    :: slope
      real(musica_dk)    :: a1, a2, b1, b2

      ierr   = 0_musica_ik
      nfrom  = size(xfrom)
      nto    = size(xto)
      nfromm1 = nfrom - iONE
      ntom1   = nto - iONE

!-----------------------------------------------------------------------------*
!  check data grid for monotonicity
!-----------------------------------------------------------------------------*
      if( any( xfrom(2:nfrom) <= xfrom(1:nfromm1) ) ) then
        call die_msg( 100000051, Iam//'data grid not monotonically increasing' )
      endif
!-----------------------------------------------------------------------------*
!  do model grid x values lie completley inside data grid x values?
!-----------------------------------------------------------------------------*
      IF ( (xfrom(1) > xto(1)) .or. (xfrom(nfrom) < xto(nto)) ) THEN
        call die_msg( 100000052, Iam//'Data do not span grid; Use ADDPNT to expand data and re-run.' )
      ENDIF
!-----------------------------------------------------------------------------*
!  find the integral of each grid interval and use this to 
!  calculate the average y value for the interval      
!  xtol and xtou are the lower and upper limits of the to grid interval
!-----------------------------------------------------------------------------*
to_interval_loop: &
      do i = iONE,ntom1
        xtol = xto(i)
        xtou = xto(i+1)
        a1 = MAX(xfrom(1),xtol)
        a2 = MIN(xfrom(nfrom),xtou)
        if( a2 > a1 ) then
          area = rZERO
from_interval_loop: &
          do k = 1,nfromm1
            a1 = MAX(xfrom(k),xtol)
            a2 = MIN(xfrom(k+1),xtou)
            if( a2 > a1 ) then
              slope = (yfrom(k+1) - yfrom(k))/(xfrom(k+1) - xfrom(k))
              b1 = yfrom(k) + slope*(a1 - xfrom(k))
              b2 = yfrom(k) + slope*(a2 - xfrom(k))
              area = area + .5_musica_dk*(a2 - a1)*(b2 + b1)
            endif
          enddo from_interval_loop
          yto(i) = area/(xtou - xtol)
        else
          yto(i) = rZERO
        endif
      enddo to_interval_loop

      end subroutine inter2

      subroutine inter4(xto,yto, xfrom,yfrom, FoldIn)
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
!=  NG     - INTEGER, number of bins + 1 in the target grid               (I)=*
!=  XG     - REAL, target grid (e.g. working wavelength grid);  bin i     (I)=*
!=           is defined as [XG(i),XG(i+1)] (i = 1..NG-1)                     =*
!=  YG     - REAL, y-data re-gridded onto XG;  YG(i) specifies the        (O)=*
!=           y-value for bin i (i = 1..NG-1)                                 =*
!=  N      - INTEGER, number of bins + 1 in the input grid                (I)=*
!=  X      - REAL, input grid (e.g. data wavelength grid);  bin i is      (I)=*
!=           defined as [X(i),X(i+1)] (i = 1..N-1)                           =*
!=  Y      - REAL, input y-data on grid X;  Y(i) specifies the            (I)=*
!=           y-value for bin i (i = 1..N-1)                                  =*
!=  FoldIn - Switch for folding option of "overhang" data                 (I)=*
!=           FoldIn = 0 -> No folding of "overhang" data                     =*
!=           FoldIn = 1 -> Integerate "overhang" data and fold back into     =*
!=                         last target bin                                   =*
!-----------------------------------------------------------------------------*

! input:
      REAL(musica_dk), intent(in)    :: xfrom(:), yfrom(:)
      REAL(musica_dk), intent(in)    :: xto(:)

      INTEGER(musica_ik), intent(in) :: FoldIn

! output:
      REAL(musica_dk), intent(out)   :: yto(:)

! local:
      character(len=*), parameter :: Iam = 'inter4: '

      INTEGER(musica_ik) :: nfrom, nto
      INTEGER(musica_ik) :: nfromm1, ntom1
      REAL(musica_dk) :: a1, a2, sum
      REAL(musica_dk) :: tail
      INTEGER(musica_ik) :: jstart, i, j, k

! check whether flag given is legal
      IF ( FoldIn /= 0_musica_ik .AND. FoldIn /= 1_musica_ik ) THEN
        call die_msg( 100000053, Iam//'Value for FOLDIN invalid. Must be 0 or 1' )
      ENDIF
      nfrom   = size(xfrom)
      nfromm1 = nfrom - iONE
      nto     = size(xto)
      ntom1   = nto - iONE
! do interpolation
      jstart = 1
      do i = 1, ntom1
        yto(i) = rZERO
        sum   = rZERO
        j = jstart
        if( j <= nfromm1 ) then
          do while( xfrom(j+1) < xto(i) )
            jstart = j
            j = j+1
            IF( j > nfromm1 ) then
              exit
            endif
          enddo

          do while( (xfrom(j) <= xto(i+1)) .and. (j <= nfromm1) )
            a1 = MAX(xfrom(j),xto(i))
            a2 = MIN(xfrom(j+1),xto(i+1))
            sum = sum + yfrom(j) * (a2-a1)
            j = j+1
          enddo
          yto(i) = sum /(xto(i+1) - xto(i))
        endif
      enddo

! if requested, integrate data "overhang" and fold back into last bin
      if( FoldIn == iONE ) then
        j = j-1
        a1 = xto(nto)       ! upper limit of last interpolated bin
        a2 = xfrom(j+1)     ! upper limit of last input bin considered
! do folding only if grids don't match up and there is more input 
        if( a2 > a1 .or. j+1 < nfrom ) then
          tail = yfrom(j) * (a2-a1)/(xfrom(j+1) - xfrom(j))
          do k = j+1, nfromm1
            tail = tail + yfrom(k) * (xfrom(k+1) - xfrom(k))
          enddo
          yto(ntom1) = yto(ntom1) + tail
        endif
      endif

      end subroutine inter4

      SUBROUTINE addpnt ( x, y, xnew, ynew )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!=  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!=  Y    - REAL vector of length LD, y-values                            (IO)=*
!=  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!=  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*

! arguments
      REAL(musica_dk), allocatable, intent(inout) :: x(:), y(:)
      REAL(musica_dk), intent(in)  :: xnew, ynew

! local variables
      character(len=*), parameter :: Iam = 'addpnt: '

      INTEGER(musica_ik) :: n
      INTEGER(musica_ik) :: insertNdx
      REAL(musica_dk), allocatable    :: wrk(:)
      logical            :: found
      

      n = size(x)
!-----------------------------------------------------------------------------*
!  check data grid for monotonicity
!-----------------------------------------------------------------------------*
      if( any( x(2:n) <= x(1:n-1) ) ) then
        write(*,*) Iam,'grid not monotonically increasing'
        stop 'GridErr'
      endif
!-----------------------------------------------------------------------------*
!  does xnew == any x value?
!-----------------------------------------------------------------------------*
      if( any( x(:) == xnew ) ) then
        write(*,*) Iam,'xnew exactly matches a grid x value'
        stop 'GridErr'
      endif
!-----------------------------------------------------------------------------*
! find the index at which xnew needs to be inserted into x
!-----------------------------------------------------------------------------*
      found = .true.
      if( xnew < x(1) ) then
        insertNdx = 1
      else if( xnew > x(n) ) then
        insertNdx = n
      else
        found = .false.
        do insertNdx = 2,n
          if (x(insertNdx) > xnew ) then
            found = .true.
            exit
          endif
        enddo
      endif
      if( .not. found ) then
        write(*,*) Iam,'something really wrong; all stop'
        stop 'codeErr'
      endif
!-----------------------------------------------------------------------------*
! increment x,y arrays, then insert xnew,ynew
!-----------------------------------------------------------------------------*
      if( insertNdx == 1 ) then
        x = [xnew,x]
        y = [ynew,y]
      elseif( insertNdx == n ) then
        x = [x,xnew]
        y = [y,ynew]
      else
        wrk = [x(:insertNdx-1),xnew]
        x   = [wrk,x(insertNdx:)]
        wrk = [y(:insertNdx-1),ynew]
        y   = [wrk,y(insertNdx:)]
      endif

!     x = [x,xnew]
!     y = [y,ynew]
!     if( xnew < x(n) ) then
!       x(insertNdx:) = eoshift( x(insertNdx:),shift=-1,boundary=xnew )
!       y(insertNdx:) = eoshift( y(insertNdx:),shift=-1,boundary=ynew )
!     endif

      end subroutine addpnt

  end module photo_utils
