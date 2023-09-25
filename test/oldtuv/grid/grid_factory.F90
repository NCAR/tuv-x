! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_grid_factory module

!> Build grid objects
module micm_grid_factory

  use micm_1d_grid,                only : base_grid_t
  use micm_1d_equal_delta_grid,    only : equalDelta_t

  implicit none

  private
  public :: grid_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function grid_builder( config ) result( new_grid_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use micm_1d_equal_delta_grid,      only : equalDelta_t
    use micm_1d_from_csv_file_grid,    only : fromCsvFile_t
    use micm_1d_grid_from_config,      only : fromConfig_t

    !> Arguments
    !> Grid configuration data
    type(config_t), intent(inout) :: config

    !> New grid object
    class(base_grid_t), pointer :: new_grid_t

    !> Local variables
    type(string_t) :: grid_type
    character(len=*), parameter :: Iam = 'Grid builder: '

    write(*,*) Iam,'entering'
    new_grid_t => null()
    call config%get( 'Grid type', grid_type, Iam )

    select case( grid_type%to_char() )
      case( 'Equal interval' )
        allocate( equalDelta_t :: new_grid_t )
      case( 'From csv file' )
        allocate( fromCsvFile_t :: new_grid_t )
      case( 'From config file' )
        allocate( fromConfig_t :: new_grid_t )
      case default
        call die_msg( 460768215, "Invalid grid type: '" // grid_type%to_char()//"'" )
    end select

    call new_grid_t%initialize( config )
    write(*,*) Iam,'exiting'

  end function grid_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_grid_factory
