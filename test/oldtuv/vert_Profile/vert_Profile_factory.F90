! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_Profile_factory module

!> Build Profile objects
module micm_Profile_factory


  implicit none

  private
  public :: Profile_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function Profile_builder( config, gridWareHouse ) result( new_Profile_t )

    use musica_assert,                   only : die_msg
    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_from_csv_file_Profile,    only  : fromCsvFile_t
    use micm_air_from_csv_file_Profile, only : airfromCsvFile_t
    use micm_o2_from_csv_file_Profile, only  : o2fromCsvFile_t
    use micm_o3_from_csv_file_Profile, only  : o3fromCsvFile_t
    use micm_Profile,                   only : base_profile_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_grid_warehouse,             only : grid_warehouse_t

    !> Arguments
    !> Grid configuration data
    type(config_t), intent(inout)         :: config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> New Profile object
    class(base_profile_t), pointer :: new_Profile_t

    !> Local variables
    character(len=*), parameter :: Iam = 'Vert Profile builder: '
    type(string_t) :: Profile_type

    new_Profile_t => null()
    call config%get( 'Vert Profile type', Profile_type, Iam )

    select case( Profile_type%to_char() )
      case( 'From csv file' )
        allocate( fromCsvFile_t :: new_Profile_t )
      case( 'Air from csv file' )
        allocate( airfromCsvFile_t :: new_Profile_t )
      case( 'O2 from csv file' )
        allocate( o2fromCsvFile_t :: new_Profile_t )
      case( 'O3 from csv file' )
        allocate( o3fromCsvFile_t :: new_Profile_t )
      case default
        call die_msg( 460768215, "Invalid Profile type: '" // Profile_type%to_char()//"'" )
    end select

    call new_Profile_t%initialize( config, gridWareHouse )

  end function Profile_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_Profile_factory
