      SUBROUTINE terint(ng,xg,yg,  n,x,y, c1,c2)

      IMPLICIT NONE

* SUBROUTINE to TERminate and INTerpolate a 1D input data array.
* INPUTS:
* working grid:

      INTEGER, intent(in) :: ng
      REAL, intent(in)    :: xg(ng)

* data array to be gridded:

      INTEGER, intent(in) :: n
      REAL, intent(in)    :: x(n), y(n)

* multipliers for 4 y-values to be added:

      INTEGER, intent(in) :: c1, c2

* OUTPUT:
* gridded y-values, with proper termination at both ends

      REAL, intent(out)  :: yg(ng)

* INTERNAL:
* input terminator points to be added:
      INTEGER :: i, ig

      REAL :: xt1, yt1
      REAL :: xt2, yt2
      REAL :: xt3, yt3
      REAL :: xt4, yt4

* delta for adding points at beginning or end of data grids

      REAL, PARAMETER :: deltax = 1.E-5
      REAL, PARAMETER :: large = 1.e36

* error flag:

      INTEGER :: ierr

* extended array:

      INTEGER :: n1, k1
      REAL    :: x1(n+4), y1(n+4)

* check for valid input values of c1, c2:

      ierr = 0
      if(c1 /= 0 .and. c1 /= 1) ierr = 1
      if(c2 /= 0 .and. c2 /= 1) ierr = 1
      IF (ierr /= 0) THEN
         WRITE(*,*) ierr, 'entering terint.f'
         STOP
      ENDIF

* assign x-coordinate to 4 new points

      xt1 = -large
      xt2 = x(1)*(1. - deltax)
      xt3 = x(n)*(1. + deltax)
      xt4 = large

* terminator multipliers  c1,c2 = e.g., = 1.,0.
* there are really only two practical options:  constant or zero

      yt1 = real(c1) * y(1)
      yt2 = real(c1) * y(1)

      yt3 = real(c2) * y(n)
      yt4 = real(c2) * y(n)

* transcribe input x,y array to avoid modifying it

      n1 = n
      x1(:n) = x(:n)
      y1(:n) = y(:n)

* extended data array, with up to 4 terminator points:
* note that n1 gets incremented by 1 for each call to ADDPNT      

      k1 = n1 + 4
      call addpnt(x1,y1,k1,n1, xt1,yt1)
      call addpnt(x1,y1,k1,n1, xt2,yt2)
      call addpnt(x1,y1,k1,n1, xt3,yt3)
      call addpnt(x1,y1,k1,n1, xt4,yt4)

* point to grid interpolation:

      CALL inter2(ng,xg,yg, n1,x1,y1,ierr)
      IF (ierr /= 0) THEN
         WRITE(*,*) ierr, 'exiting terint.f'
         STOP
      ENDIF

      END SUBROUTINE terint
