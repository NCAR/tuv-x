! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
! Profile specified in json config file
module micm_Profile_from_config

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_Profile,     only : base_profile_t

  implicit none

  private
  public :: fromConfig_t

  type, extends(base_profile_t) :: fromConfig_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type fromConfig_t

contains
  !> Initialize grid
  subroutine initialize( this, Profile_config, gridWareHouse )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : die_msg
    use micm_1d_grid,  only : base_grid_t
    use micm_grid_warehouse,  only : grid_warehouse_t

    !> Arguments
    class(fromConfig_t), intent(inout) :: this
    type(config_t), intent(inout)      :: Profile_config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> Local variables
    character(len=*), parameter :: Iam = 'From config Profile initialize: '
    integer(ik)                 :: ndx
    real(dk)                    :: uniformValue
    logical(lk)                 :: found
    type(string_t)              :: gridHandle
    class(base_grid_t), pointer :: theGrid
 
    !> Get the handle
    call Profile_config%get( 'Handle', this%handle_, Iam, default = 'None' )

    !> Get values from config file
    call Profile_config%get( "Values", this%edge_val_, Iam, found=found )
    if( .not. found ) then
      call Profile_config%get( "Uniform value", uniformValue, Iam, found=found )
      if( found ) then
        call Profile_config%get( "Grid", gridHandle, Iam, found=found )
        if( found ) then
          theGrid => gridWareHouse%get_grid( gridHandle )
          this%edge_val_ = (/ (uniformValue,ndx=1,theGrid%ncells_+1_ik) /)
        else
          call die_msg( 123456,"Grid " // gridHandle%to_char() // " not in grid warehouse" )
        endif
      else
        call die_msg( 123457,"Neither 'Values' or 'Uniform value' keyword specified" )
      endif
    endif

    this%ncells_ = size(this%edge_val_) - 1_ik
    this%mid_val_ = .5_dk &
                   *(this%edge_val_(1_ik:this%ncells_) + this%edge_val_(2_ik:this%ncells_+1_ik))
    this%delta_val_ = (this%edge_val_(2_ik:this%ncells_+1_ik) - this%edge_val_(1_ik:this%ncells_))

  end subroutine initialize

end module micm_Profile_from_config
