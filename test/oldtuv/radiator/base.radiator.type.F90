! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The base_micm_radiator module

!> The base_radiator_t type and related functions
!!
module micm_base_radiator_type

  use musica_constants,       only : dk => musica_dk, ik => musica_ik
  use musica_string,          only : string_t
  use micm_abs_radiator_type, only : abs_radiator_t

  implicit none

  private
  public :: base_radiator_t

  !> base radiator type
  type, extends(abs_radiator_t) :: base_radiator_t
  contains
    !> Initialize radiator
    procedure :: initialize
    !> Update radiator for new environmental conditions
    procedure :: upDateState
  end type base_radiator_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize radiator_t object
  subroutine initialize( this, radiator_config, gridWareHouse )

    use musica_config,        only : config_t
    use micm_grid_warehouse,  only : grid_warehouse_t
    use micm_1d_grid,         only : base_grid_t

    !> radiator object
    class(base_radiator_t), intent(inout) :: this
    !> radiator configuration object
    type(config_t), intent(inout)         :: radiator_config
    !> grid warehouse
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> local variables
    character(len=*), parameter   :: Iam = "Radiator initialize: "
    type(string_t)                :: Handle
    class(base_grid_t), pointer :: zGrid, lambdaGrid

    write(*,*) ' '
    write(*,*) Iam,'entering'

    Handle = 'Vertical Z' ; zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength' ; lambdaGrid => gridWareHouse%get_grid( Handle )

!-----------------------------------------------------------------------------
!> Get radiator "Handle"
!-----------------------------------------------------------------------------
    call radiator_config%get( 'Handle', this%handle_, Iam )
    write(*,*) Iam // 'handle = ',this%handle_%to_char()

!> allocate radiator state_ variables
    allocate( this%state_%layer_OD_(zGrid%ncells_,lambdaGrid%ncells_) )
    allocate( this%state_%layer_SSA_(zGrid%ncells_,lambdaGrid%ncells_) )
    allocate( this%state_%layer_G_(zGrid%ncells_,lambdaGrid%ncells_) )
    write(*,*) Iam // 'state_%layer_OD_ is allocated = ',allocated(this%state_%layer_OD_)
    write(*,*) Iam // 'state_%layer_SSA_ is allocated = ',allocated(this%state_%layer_SSA_)
    write(*,*) Iam // 'state_%layer_G_ is allocated = ',allocated(this%state_%layer_G_)
    write(*,*) Iam // 'state_%layer_OD_ is (',size(this%state_%layer_OD_,dim=1),' x ',size(this%state_%layer_OD_,dim=2),')'

    write(*,*) ' '
    write(*,*) Iam,'exiting'

  end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update radiator state
  subroutine upDateState( this, gridWareHouse, ProfileWareHouse, radXferXsectWareHouse )

    use musica_assert,                 only : die_msg
    use micm_Profile_warehouse,        only : Profile_warehouse_t
    use micm_Profile,                  only : base_profile_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_1d_grid,                  only : base_grid_t
    use micm_radXfer_xsect_warehouse,  only : radXfer_xsect_warehouse_t
    use micm_radXfer_abs_cross_section_type, only : abs_cross_section_t
    use debug,                         only : diagout

    !> Arguments
    !> radiator obj
    class(base_radiator_t), intent(inout) :: this
    !> Grid warehouse
    type(grid_warehouse_t), intent(inout) :: gridWareHouse
    !> Profile warehouse
    type(Profile_warehouse_t), intent(inout)       :: ProfileWareHouse
    !> RadXfer cross section warehouse
    type(radXfer_xsect_warehouse_t), intent(inout) :: radXferXsectWareHouse

    !> Local variables
    real(dk) , parameter  :: km2cm = 1.e5_dk

    integer(ik) :: wNdx
    real(dk), allocatable :: CrossSection(:,:)
    character(len=*), parameter :: Iam = 'base radiator upDateState: '
    type(string_t)      :: Handle
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid
    class(base_profile_t), pointer  :: radiatorProfile
    class(abs_cross_section_t), pointer :: radiatorCrossSection

    write(*,*) ' '
    write(*,*) Iam,'entering'

    write(*,*) Iam // 'handle = ',this%handle_%to_char()
!-----------------------------------------------------------------------------
!> get specific grids and profiles
!-----------------------------------------------------------------------------
    Handle = 'Vertical Z'
    zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength'
    lambdaGrid => gridWareHouse%get_grid( Handle )
    write(*,*) Iam // 'nlyr,nbins = ',zGrid%ncells_,lambdaGrid%ncells_

    !> Note: uses radiator handle for Profile handle
    radiatorProfile => ProfileWareHouse%get_Profile( this%handle_ )

    !> Note: uses radiator handle for cross section handle
    radiatorCrossSection => radXferXsectWareHouse%get_radXfer_cross_section( this%handle_ )

    !> check radiator state type allocation
    if( .not. allocated( this%state_%layer_OD_ ) ) then
      call die_msg( 2222222,"In radiator%upDateState radiator state not allocate" )
    else
      write(*,*) Iam // 'radiator state is allocated'
    endif
    write(*,*) Iam // 'size OD = ',size(this%state_%layer_OD_,dim=1),' x ', &
                                  size(this%state_%layer_OD_,dim=2)

    !> set radiator state members
    CrossSection = radiatorCrossSection%calculate( gridWareHouse, ProfileWareHouse )
    call diagout( 'o2xs.new',CrossSection )
    do wNdx = 1,lambdaGrid%ncells_
      this%state_%layer_OD_(:,wNdx) = radiatorProfile%layer_dens_ * CrossSection(:,wNdx)
    enddo

    !> Settings for a gas phase radiator
    if( this%handle_ == 'Air' ) then
      this%state_%layer_SSA_ = 1._dk
      this%state_%layer_G_   = 0._dk
    else
      this%state_%layer_SSA_ = 0._dk
      this%state_%layer_G_   = 0._dk
    endif

    write(*,*) ' '
    write(*,*) Iam,'exiting'

  end subroutine upDateState

end module micm_base_radiator_type
