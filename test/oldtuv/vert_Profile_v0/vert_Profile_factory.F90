! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_vert_Profile_factory module

!> Build vert_Profile objects
module micm_vert_Profile_factory


  implicit none

  private
  public :: vert_Profile_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function vert_Profile_builder( config, gridWareHouse ) result( new_vert_Profile_t )

    use musica_assert,                   only : die_msg
    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use micm_from_csv_file_vert_Profile,    only  : fromCsvFile_t
    use micm_air_from_csv_file_vert_Profile, only : airfromCsvFile_t
    use micm_o2_from_csv_file_vert_Profile, only  : o2fromCsvFile_t
    use micm_o3_from_csv_file_vert_Profile, only  : o3fromCsvFile_t
    use micm_vert_Profile,               only : abs_vert_Profile_t
    use micm_1d_grid,                    only : base_grid_t
    use micm_grid_warehouse,             only : grid_warehouse_t

    !> Arguments
    !> Grid configuration data
    type(config_t), intent(inout)         :: config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> New vert Profile object
    class(abs_vert_Profile_t), pointer :: new_vert_Profile_t

    !> Local variables
    character(len=*), parameter :: Iam = 'Vert Profile builder: '
    type(string_t) :: vert_Profile_type

    write(*,*) Iam,'entering'

    new_vert_Profile_t => null()
    call config%get( 'Vert Profile type', vert_Profile_type, Iam )

    select case( vert_Profile_type%to_char() )
      case( 'From csv file' )
        allocate( fromCsvFile_t :: new_vert_Profile_t )
      case( 'Air from csv file' )
        allocate( airfromCsvFile_t :: new_vert_Profile_t )
      case( 'O2 from csv file' )
        allocate( o2fromCsvFile_t :: new_vert_Profile_t )
      case( 'O3 from csv file' )
        allocate( o3fromCsvFile_t :: new_vert_Profile_t )
      case default
        call die_msg( 460768215, "Invalid vert Profile type: '" // vert_Profile_type%to_char()//"'" )
    end select

    call new_vert_Profile_t%initialize( config, gridWareHouse )

    write(*,*) Iam,'exiting'

  end function vert_Profile_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_vert_Profile_factory
