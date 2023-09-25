! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_from_csv_file
! read a grid defined in a CSV file. See
! :ref:`configuration-grids` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid,                       only : grid_t

  implicit none

  public :: grid_from_csv_file_t

  type, extends(grid_t) :: grid_from_csv_file_t
  contains
  end type grid_from_csv_file_t

  !> Constructor
  interface grid_from_csv_file_t
    module procedure constructor
  end interface grid_from_csv_file_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result ( this )
    ! Initialize grid

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(config_t), intent(inout) :: config ! The grid config. See :ref:`configuration-grids` for more details

    ! local variables
    integer, parameter :: Ok = 0
    integer, parameter :: inUnit = 20
    character(len=*), parameter :: Iam = 'From_csv_file grid initialize: '
    type(grid_from_csv_file_t), pointer :: this

    integer  :: istat
    real(dk) :: Value
    logical  :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec
    type(string_t) :: required_keys(4), optional_keys(0)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "file path"
    required_keys(4) = "name"

    call assert_msg( 606960546,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "grid from CSV file." )

    allocate( this )

    call config%get( 'file path', Filespec, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )

    inquire( file = Filespec%to_char( ), exist = found )
    if( .not. found ) then
      call die_msg( 560768215, "File " // Filespec%to_char( ) // " not found" )
    endif

    open( unit = inUnit, file = Filespec%to_char( ), iostat = istat )
    if( istat /= Ok ) then
        call die_msg( 560768225, "Error reading " // Filespec%to_char( ) )
    endif

    ! The first line of the file contains the number of grid cells,
    ! throw it away
    read( inUnit, *, iostat = istat ) InputLine
    if( istat /= Ok ) then
      call die_msg( 560768226, "Error reading " // Filespec%to_char( ) )
    endif

    ! Now read the gride edge boundaries until the end of the file
    allocate( this%edge_(0) )
    do
      read( inUnit, *, iostat = istat ) Value
      if( istat /= Ok ) then
        exit
      endif
      this%edge_ = [ this%edge_, Value ]
    enddo

    this%ncells_ = size( this%edge_ ) - 1
    allocate( this%mid_( this%ncells_ ) )
    allocate( this%delta_( this%ncells_ ) )
    this%mid_(:) = .5_dk *                                                    &
      ( this%edge_( 1 : this%ncells_ ) + this%edge_( 2 : this%ncells_ + 1 ) )
    this%delta_(:) = &
      this%edge_( 2 : this%ncells_ + 1 ) - this%edge_( 1 : this%ncells_ )
    close( unit = inUnit )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_from_csv_file
