! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_util
  ! Utility functions used in TUV-x

  use musica_constants,                only : musica_dk

  implicit none

  private
  public :: add_point

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine add_point( x, y, xnew, ynew )
    ! Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in
    ! ascending order

    real(musica_dk), allocatable, intent(inout) :: x(:) ! Grid
    real(musica_dk), allocatable, intent(inout) :: y(:) ! Data
    real(musica_dk),              intent(in)    :: xnew ! Grid point to add
    real(musica_dk),              intent(in)    :: ynew ! Data point to add

    character(len=*), parameter :: Iam = 'Add point to gridded data'
    integer :: n
    integer :: insertNdx
    real(musica_dk), allocatable    :: wrk(:)
    logical            :: found

    n = size( x )

    !  check data grid for monotonicity
    if( any( x( 2 : n ) <= x( 1 : n - 1 ) ) ) then
      write(*,*) Iam, 'grid not monotonically increasing'
      stop 3
    endif

    !  does xnew == any x value?
    if( any( x(:) == xnew ) ) then
      write(*,*) Iam, 'xnew exactly matches a grid x value'
      stop 3
    endif

    ! find the index at which xnew needs to be inserted into x
    found = .true.
    if( xnew < x(1) ) then
      insertNdx = 1
    else if( xnew > x( n ) ) then
      insertNdx = n + 1
    else
      found = .false.
      do insertNdx = 2, n
        if ( x( insertNdx ) > xnew ) then
          found = .true.
          exit
        endif
      enddo
    endif
    if( .not. found ) then
      write(*,*) Iam, 'something really wrong; all stop'
      stop 3
    endif

    ! increment x,y arrays, then insert xnew,ynew
    if( insertNdx == 1 ) then
      x = [ xnew, x ]
      y = [ ynew, y ]
    elseif( insertNdx == n + 1 ) then
      x = [ x, xnew ]
      y = [ y, ynew ]
    else
      wrk = [ x( : insertNdx - 1 ), xnew ]
      x   = [ wrk, x( insertNdx : ) ]
      wrk = [ y( : insertNdx - 1 ), ynew ]
      y   = [ wrk, y( insertNdx : ) ]
    endif

  end subroutine add_point

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_util
