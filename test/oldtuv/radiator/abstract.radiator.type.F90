! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The radiator module

!> The abstract radiator type and related functions
module micm_abs_radiator_type

  use musica_constants,                only : musica_dk
  use micm_environment,                only : environment_t
  use musica_string,                   only : string_t

  implicit none
  private

  public :: abs_radiator_t, radiator_state_t, radiator_ptr

  !> radiator state type
  type :: radiator_state_t
    !> layer optical depth
    real(kind=musica_dk), allocatable :: layer_OD_(:,:)
    !> layer single scattering albedo
    real(kind=musica_dk), allocatable :: layer_SSA_(:,:)
    !> layer asymmetry factor
    real(kind=musica_dk), allocatable :: layer_G_(:,:)
  contains
    final :: finalize
  end type radiator_state_t

  !> radiator abstract type
  type, abstract :: abs_radiator_t
    type(string_t)         :: handle_
    type(radiator_state_t) :: state_
  contains
    !> Calculate the photo rate cross section
    procedure(initial),     deferred :: initialize
    procedure(upDateState), deferred :: upDateState
  end type abs_radiator_t

  !> Pointer type for building sets of radiator objects
  type :: radiator_ptr
    class(abs_radiator_t), pointer :: val_ => null( )
  end type radiator_ptr

interface

  !> Initialize the base cross section type
  subroutine initial( this, radiator_config, gridWareHouse )
    use musica_config,        only : config_t
    use micm_grid_warehouse,  only : grid_warehouse_t

    import abs_radiator_t

    !> radiator obj
    class(abs_radiator_t), intent(inout)  :: this
    !> radiator configuration
    type(config_t),        intent(inout)  :: radiator_config
    !> grid warehouse
    type(grid_warehouse_t), intent(inout) :: gridWareHouse
 end subroutine initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate cross section and layer density
  subroutine upDateState( this, gridWareHouse, ProfileWareHouse, radXferXsectWareHouse )

    use micm_Profile_warehouse,        only : Profile_warehouse_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_radXfer_xsect_warehouse,  only : radXfer_xsect_warehouse_t
    use musica_constants,              only : musica_dk

    import abs_radiator_t

    !> Update radiator state
    class(abs_radiator_t), intent(inout)           :: this
    type(grid_warehouse_t), intent(inout)          :: gridWareHouse
    type(Profile_warehouse_t), intent(inout)       :: ProfileWareHouse
    type(radXfer_xsect_warehouse_t), intent(inout) :: radXferXsectWareHouse
    !> Radiator state
  end subroutine upDateState

end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radiator state obj
  subroutine finalize( this )

  !> object declaration
  type(radiator_state_t) :: this

  if( allocated( this%layer_OD_ ) ) then
    deallocate( this%layer_OD_ )
  endif
  if( allocated( this%layer_SSA_ ) ) then
    deallocate( this%layer_SSA_ )
  endif
  if( allocated( this%layer_G_ ) ) then
    deallocate( this%layer_G_ )
  endif

  end subroutine finalize

end module micm_abs_radiator_type
