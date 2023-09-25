! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_spherical_geometry
  ! Calculates paths in a spherical geometry

      use musica_constants,            only : dk => musica_dk
      use tuvx_constants,              only : radius, pi
      use tuvx_grid_warehouse,         only : grid_warehouse_ptr

      implicit none

      private
      public :: spherical_geometry_t

      type :: spherical_geometry_t
        ! Vertical and slant column calculator
        integer, allocatable  :: nid_(:) ! number of layers crossed by the direct beam when travelling from the top of the atmosphere to layer i
        real(dk)              :: solar_zenith_angle_ ! the solar zenith angle in degrees
        real(dk), allocatable :: dsdh_(:,:) ! slant path of direct beam through each layer crossed when travelling from the top of the atmosphere to layer i
        type(grid_warehouse_ptr) :: height_grid_ ! pointer to the height grid in the grid warehouse
      contains
        procedure :: set_parameters
        procedure :: air_mass
        ! Returns the number of bytes needed to pack the calculator onto a
        ! buffer
        procedure :: pack_size
        ! Packs the calculator onto a character buffer
        procedure :: mpi_pack
        ! Unpacks a calculator from a character buffer
        procedure :: mpi_unpack
        final     :: finalize
      end type spherical_geometry_t

      real(dk), parameter ::  rZERO   = 0.0_dk
      real(dk), parameter ::  d2r     = pi / 180._dk
      real(dk), parameter ::  NINETY  = 90._dk

      interface spherical_geometry_t
        module procedure constructor
      end interface spherical_geometry_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( grid_warehouse ) result( this )
    ! Creates a spherical geometry calculator

    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_grid,                     only : grid_t

    type(spherical_geometry_t), pointer       :: this ! A :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`
    type(grid_warehouse_t),     intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    class(grid_t), pointer :: zGrid

    allocate( this )

    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    zGrid => grid_warehouse%get_grid( this%height_grid_ )

    allocate( this%nid_( 0 : zGrid%ncells_ ) )
    allocate( this%dsdh_( 0 : zGrid%ncells_, zGrid%ncells_ ) )

    this%nid_(:) = 0.0_dk
    this%solar_zenith_angle_ = 0.0_dk
    this%dsdh_(:,:) = 0.0_dk

    deallocate( zGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine set_parameters(this, zen, grid_warehouse)
    ! Calculates the slant path over vertical depth ds/dh in spherical geometry.
    ! calculation is based on:  a.dahlback, and k.stamnes, a new spheric model
    ! for computing the radiation field available for photolysis and heating
    ! at twilight, planet.space sci., v39, n5,
    ! pp. 671-683, 1991 (appendix b) `doi:10.1016/0032-0633(91)90061-E
    ! <https://doi.org/10.1016/0032-0633(91)90061-E>`_

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    real(dk),                    intent(in)    :: zen ! solar zenith angle (degrees)
    class(spherical_geometry_t), intent(inout) :: this ! A :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`
    type(grid_warehouse_t),      intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    real(dk) :: re
    real(dk), allocatable :: ze(:)

    integer :: nz ! number of specified altitude levels in the working grid
    integer :: i, j
    integer :: id
    integer :: nlayer
    real(dk)    :: sinrad, zenrad, rpsinz, rj, rjp1, dsj, dhj, ga, gb, sm
    real(dk), allocatable    :: zd(:)

    class(grid_t), pointer :: zGrid

    zenrad = zen * d2r

    zGrid => grid_warehouse%get_grid( this%height_grid_ )

    nlayer = zGrid%ncells_
    nz     = nlayer + 1

    ! include the elevation above sea level to the radius of the earth:
    re = radius + zGrid%edge_(1)
    ! correspondingly z changed to the elevation above earth surface:
    ze = zGrid%edge_ - zGrid%edge_(1)

    allocate( zd( 0 : nlayer ) )
    ! inverse coordinate of z
    zd( 0 : nlayer ) = ze( nz : 1 : -1 )

    ! initialize dsdh, nid
    this%nid_  = 0
    this%dsdh_ = rZERO

    sinrad = sin( zenrad )
    ! calculate ds/dh of every layer
    layer_loop: do i = 0, nlayer
      rpsinz = ( re + zd( i ) ) * sinrad
      if ( zen > NINETY .and. rpsinz < re ) then
        id = -1
      else
        ! find index of layer in which the screening height lies
        id = i
        if( zen > NINETY ) then
          id = -1
          do j = 1, nlayer
            if( rpsinz < ( zd( j - 1 ) + re ) .and.                           &
                rpsinz >= ( zd( j ) + re) ) id = j
          enddo
        end if

        do j = 1, id
          sm = 1.0_dk
          if( j == id .and. id == i .and. zen > NINETY) sm = -1.0_dk
          rj = re + zd( j - 1 )
          rjp1 = re + zd( j )
          dhj = zd( j - 1 ) - zd( j )
          ga = rj * rj - rpsinz * rpsinz
          gb = rjp1 * rjp1 - rpsinz * rpsinz
          ga = max( rZERO, ga )
          gb = max( rZERO, gb )

          if( id > i .and. j == id ) then
            dsj = sqrt( ga )
          else
            dsj = sqrt( ga ) - sm * sqrt( gb )
          end if
          this%dsdh_( i, j ) = dsj / dhj
        enddo
      end if
      this%nid_( i ) = id
    enddo layer_loop

    this%solar_zenith_angle_ = zen

    deallocate( zGrid )

  end subroutine set_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine air_mass( this, aircol, vcol, scol )
    !  calculate vertical and slant air columns, in spherical geometry, as a
    !  function of altitude.

    use tuvx_constants, only : largest

    class(spherical_geometry_t), intent(in)  :: this ! This :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`
    real(dk),                    intent(in)  :: aircol(:) ! layer density of air (height) [molec cm-2]
    real(dk),                    intent(out) :: vcol(:) ! vertical air column, molec cm-2, above level iz
    real(dk),                    intent(out) :: scol(:) ! slant air column in direction of sun, above iz also in molec cm-2

    integer :: nz ! number of specified altitude levels in the working (i) grid
    integer :: nlayer
    integer :: id, j
    real(dk)    :: accum

    ! calculate vertical and slant column from each level:
    ! work downward

    nz = size( aircol )
    nlayer = nz - 1
    accum = aircol( nz )
    do id = nlayer, 1, -1
      accum = accum + aircol( id )
      vcol( id ) = accum
    enddo

    scol( nz ) = this%dsdh_( 1, 1 ) * aircol( nz )
    do id = 1, nlayer
       accum = scol( nz )
       if( this%nid_( id ) < 0 ) then
          accum = largest
       else
          ! single pass layers:
          do j = 1, min( this%nid_( id ), id )
             accum = accum + aircol( nz - j ) * this%dsdh_( id, j )
          enddo
          ! double pass layers:
          do j = min( this%nid_( id ), id ) + 1, this%nid_( id )
             accum = accum + 2._dk * aircol( nz - j ) * this%dsdh_( id, j )
          enddo
       endif
       scol( nz - id ) = accum
    enddo

  end subroutine air_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the calculator onto a
    ! buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack_size

    class(spherical_geometry_t), intent(in) :: this ! Calculator to pack
    integer,                     intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    real(dk), allocatable :: temp_1D(:), temp_2D(:,:)

    call assert( 919114707, allocated( this%nid_ ) )
    call assert( 296325650, lbound( this%nid_, 1 ) == 0 )
    allocate( temp_1D( size( this%nid_ ) ) )
    temp_1D(:) = this%nid_(0:)
    call assert( 459765008, allocated( this%dsdh_ ) )
    call assert( 119451200, lbound( this%dsdh_, 1 ) == 0 )
    allocate( temp_2D( size( this%dsdh_, 1 ), size( this%dsdh_, 2 ) ) )
    temp_2D(:,:) = this%dsdh_(0:,:)
    pack_size = musica_mpi_pack_size( temp_1D,                  comm ) +      &
                musica_mpi_pack_size( this%solar_zenith_angle_, comm ) +      &
                musica_mpi_pack_size( temp_2D,                  comm ) +      &
                this%height_grid_%pack_size( comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the calculator onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(spherical_geometry_t), intent(in)    :: this      ! calculator to pack
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    real(dk), allocatable :: temp_1D(:), temp_2D(:,:)

    prev_pos = position
    call assert( 175836520, allocated( this%nid_ ) )
    call assert( 623204366, lbound( this%nid_, 1 ) == 0 )
    allocate( temp_1D( size( this%nid_ ) ) )
    temp_1D(:) = this%nid_(0:)
    call assert( 788096963, allocated( this%dsdh_ ) )
    call assert( 400473210, lbound( this%dsdh_, 1 ) == 0 )
    allocate( temp_2D( size( this%dsdh_, 1 ), size( this%dsdh_, 2 ) ) )
    temp_2D(:,:) = this%dsdh_(0:,:)
    call musica_mpi_pack( buffer, position, temp_1D,                  comm )
    call musica_mpi_pack( buffer, position, this%solar_zenith_angle_, comm )
    call musica_mpi_pack( buffer, position, temp_2D,                  comm )
    call this%height_grid_%mpi_pack( buffer, position, comm )
    call assert( 326119554, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a calculator from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(spherical_geometry_t), intent(out)   :: this      ! calculator to be unpacked
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos
    real(dk), allocatable :: temp_1D(:), temp_2D(:,:)

    prev_pos = position
    call musica_mpi_unpack( buffer, position, temp_1D,                  comm )
    allocate( this%nid_( 0 : size( temp_1D ) - 1 ) )
    this%nid_(0:size(this%nid_)-1) = temp_1D(1:size(temp_1D))
    call musica_mpi_unpack( buffer, position, this%solar_zenith_angle_, comm )
    call musica_mpi_unpack( buffer, position, temp_2D,                  comm )
    allocate( this%dsdh_( 0 : size( temp_2D, 1 ) - 1, size( temp_2D, 2 ) ) )
    this%dsdh_(0:,:) = temp_2D(1:,:)
    call this%height_grid_%mpi_unpack( buffer, position, comm )
    call assert( 758599069, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Clean up memory

    type(spherical_geometry_t), intent(inout) :: this ! A :f:type:`~tuvx_spherical_geometry/spherical_geometry_t`

    if( allocated( this%nid_ ) ) then
      deallocate( this%nid_ )
    endif
    if( allocated( this%dsdh_ ) ) then
      deallocate( this%dsdh_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_spherical_geometry
