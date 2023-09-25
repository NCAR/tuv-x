! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!1d grid specified in json config file
module micm_1d_grid_from_config

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_1d_grid,     only : base_grid_t

  implicit none

  public :: fromConfig_t

  type, extends(base_grid_t) :: fromConfig_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type fromConfig_t

contains
  !> Initialize grid
  subroutine initialize( this, grid_config )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : die_msg

    !> Arguments
    class(fromConfig_t), intent(inout) :: this
    type(config_t), intent(inout)      :: grid_config

    !> Local variables
    character(len=*), parameter :: Iam = 'From config grid initialize: '
 
    !> Get the handle
    call grid_config%get( 'Handle', this%handle_, Iam, default = 'None' )

    !> Get values from config file
    call grid_config%get( "Values", this%edge_, Iam )

    this%ncells_ = size(this%edge_) - 1_ik
    this%mid_ = .5_dk &
                   *(this%edge_(1_ik:this%ncells_) + this%edge_(2_ik:this%ncells_+1_ik))
    this%delta_ = (this%edge_(2_ik:this%ncells_+1_ik) - this%edge_(1_ik:this%ncells_))

  end subroutine initialize

end module micm_1d_grid_from_config
