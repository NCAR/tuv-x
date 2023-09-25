! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
! Profile specified in json config file
module micm_srfAlbedo_Profile_from_config

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_Profile,     only : base_profile_t

  implicit none

  private
  public :: srfAlbedofromConfig_t

  type, extends(base_profile_t) :: srfAlbedofromConfig_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type srfAlbedofromConfig_t

contains
  !> Initialize grid
  subroutine initialize( this, Profile_config, gridWareHouse )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : die_msg
    use micm_1d_grid,  only : base_grid_t
    use micm_grid_warehouse,  only : grid_warehouse_t

    !> Arguments
    class(srfAlbedofromConfig_t), intent(inout) :: this
    type(config_t), intent(inout)               :: Profile_config
    type(grid_warehouse_t), intent(inout)       :: gridWareHouse

    !> Local variables
    character(len=*), parameter :: Iam = 'From config Profile initialize: '
    integer(ik)                   :: ndx
    real(dk)                      :: uniformValue
    class(base_grid_t), pointer :: lambdaGrid
    type(string_t)                :: Handle
 
    !> Get the handle
    call Profile_config%get( 'Handle', this%handle_, Iam, default = 'None' )

    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )

    !> Get values from config file
    call Profile_config%get( "Uniform Value", uniformValue, Iam )

    this%ncells_ = lambdaGrid%ncells_
    this%edge_val_ = (/ (uniformValue,ndx=1,this%ncells_+1_ik) /)
    this%mid_val_ = .5_dk &
                   *(this%edge_val_(1_ik:this%ncells_) + this%edge_val_(2_ik:this%ncells_+1_ik))
    this%delta_val_ = (this%edge_val_(2_ik:this%ncells_+1_ik) - this%edge_val_(1_ik:this%ncells_))

  end subroutine initialize

end module micm_srfAlbedo_Profile_from_config
