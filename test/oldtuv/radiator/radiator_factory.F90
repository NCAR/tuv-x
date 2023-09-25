! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_radiator_factory module

!> Build radiator objects
module micm_radiator_factory

  use micm_abs_radiator_type,       only : abs_radiator_t

  implicit none

  private
  public :: radiator_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function radiator_builder( config, gridWareHouse ) result( new_radiator_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_base_radiator_type,       only : base_radiator_t
    use micm_aerosol_radiator_type,    only : aerosol_radiator_t

    !> Arguments
    !> Radiator configuration data
    type(config_t), intent(inout)         :: config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> New radiator object
    class(abs_radiator_t), pointer :: new_radiator_t

    !> Local variables
    type(string_t) :: radiator_type
    character(len=*), parameter :: Iam = 'Radiator builder: '

    write(*,*) Iam,'entering'
    new_radiator_t => null()
    call config%get( 'radiator type', radiator_type, Iam )

    select case( radiator_type%to_char() )
      case( 'Base radiator' )
        allocate( base_radiator_t :: new_radiator_t )
      case( 'Aerosol radiator' )
        allocate( aerosol_radiator_t :: new_radiator_t )
      case default
        call die_msg( 460768245, "Invalid radiator type: '" // radiator_type%to_char()//"'" )
    end select

    call new_radiator_t%initialize( config, gridWareHouse )
    write(*,*) Iam,'exiting'

  end function radiator_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_radiator_factory
