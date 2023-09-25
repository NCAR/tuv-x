! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!one dimension, equally spaced  grid type
module micm_1d_equal_delta_grid

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_1d_grid,     only : base_grid_t

  implicit none

  public :: equalDelta_t

  type, extends(base_grid_t) :: equalDelta_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type equalDelta_t

contains
  !> Initialize grid
  subroutine initialize( this, grid_config )
      
    use musica_config, only : config_t
    use musica_string, only : string_t

    !> arguments
    class(equalDelta_t), intent(inout) :: this
    type(config_t), intent(inout)      :: grid_config
    !> local variables
    integer(ik) :: n
    real(dk)    :: Lower_val, Upper_val, Delta_val
    character(len=*), parameter :: Iam = 'EqualDelta grid initialize: '
    logical(lk) :: found

    call grid_config%get( 'Grid begins at', Lower_val, Iam )
    call grid_config%get( 'Grid ends at', Upper_val, Iam )
    call grid_config%get( 'Grid cell delta', Delta_val, Iam )
    call grid_config%get( 'Handle', this%handle_, Iam, found = found )
    if( .not. found ) then
      this%handle_ = "None"
    endif

    this%ncells_ = int( (Upper_val - Lower_val)/Delta_val,kind=ik )
    if( mod((Upper_val - Lower_val),Delta_val ) /= 0._dk ) then
      this%ncells_ = this%ncells_ + 1
    endif
    allocate( this%mid_(this%ncells_) )
    allocate( this%delta_(this%ncells_) )
    allocate( this%edge_(this%ncells_+1_ik) )
    do n = 1,this%ncells_+1_ik
      this%edge_(n) = min( real((n - 1_ik),kind=dk)*Delta_val + Lower_val,Upper_val )
    enddo
    this%mid_(:) = .5_dk &
                   *(this%edge_(1_ik:this%ncells_) + this%edge_(2_ik:this%ncells_+1_ik))
    this%delta_(:) = (this%edge_(2_ik:this%ncells_+1_ik) - this%edge_(1_ik:this%ncells_))

  end subroutine initialize

end module micm_1d_equal_delta_grid
