! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_cross_section_factory module

!> Builder of cross section calculators
module micm_radXfer_cross_section_factory

  use micm_radXfer_abs_cross_section_type,  only : abs_cross_section_t
  use micm_radXfer_base_cross_section_type, only : base_cross_section_t
  use micm_radXfer_tint_cross_section_type, only : tint_cross_section_t
  use micm_radXfer_o3_tint_cross_section_type,  only : o3_tint_cross_section_t
  use micm_radXfer_rayliegh_cross_section_type, only : rayliegh_cross_section_t

  implicit none

  private
  public :: cross_section_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function cross_section_builder( config, gridWareHouse, ProfileWareHouse ) result( new_cross_section_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use musica_constants,              only : musica_dk
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_Profile_warehouse,        only : Profile_warehouse_t

    !> New rate constant calculator
    class(abs_cross_section_t), pointer :: new_cross_section_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config
    !> The warehouses
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

    type(string_t) :: cross_section_type
    character(len=*), parameter :: Iam = 'cross section builder: '

    write(*,*) Iam,'entering'
    new_cross_section_t => null( )
    call config%get( 'cross section type', cross_section_type, Iam )

    select case( cross_section_type%to_char() )
      case( 'base cross section' )
        allocate( base_cross_section_t :: new_cross_section_t )
      case( 'tint cross section' )
        allocate( tint_cross_section_t :: new_cross_section_t )
      case( 'O3 cross section' )
        allocate( o3_tint_cross_section_t :: new_cross_section_t )
      case( 'SO2 cross section' )
        allocate( base_cross_section_t :: new_cross_section_t )
      case( 'NO2 tint cross section' )
        allocate( tint_cross_section_t :: new_cross_section_t )
      case( 'Air cross section' )
        allocate( rayliegh_cross_section_t :: new_cross_section_t )
      case default
        call die_msg( 450768214, "Invalid cross section type: '"//              &
                                 cross_section_type%to_char( )//"'" )
    end select
    call new_cross_section_t%initialize( config, gridWareHouse, ProfileWareHouse )
    write(*,*) Iam,'exiting'

  end function cross_section_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_radXfer_cross_section_factory
