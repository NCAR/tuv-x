! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp radiative transfer module
module radXfer_core

  use musica_config,    only : config_t
  use musica_string,    only : string_t
  use musica_assert,    only : assert
  use musica_constants, only : ik => musica_ik, dk => musica_dk
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_grid,    only : grid_t
  use tuvx_profile_warehouse, only : Profile_warehouse_t
  use tuvx_profile,           only : profile_t
  use tuvx_cross_section_warehouse, only : cross_section_warehouse_t
  use tuvx_cross_section, only : cross_section_t

  implicit none

  private
  public :: radXfer_core_t

  type :: radXfer_core_t
    type(grid_warehouse_t), pointer          :: theGridWarehouse_ => null( )
    type(Profile_warehouse_t), pointer       :: theProfileWarehouse_ => null( )
    type(cross_section_warehouse_t), pointer ::                               &
        theradXferXsectWarehouse_ => null( )
  contains
    procedure :: test => run
    final     :: finalize
  end type radXfer_core_t

  interface radXfer_core_t
    module procedure constructor
  end interface radXfer_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config_flsp ) result( radXfer_core_obj )

    !> Arguments
    type(string_t), intent(in) :: config_flsp

    class(radXfer_core_t), pointer :: radXfer_core_obj

    !> Local variables
    character(len=*), parameter :: Iam = 'radXfer core constructor: '
    type(config_t)              :: tst_config, child_config

    ! main configuration -> config type
    call tst_config%from_file( config_flsp%to_char() )

    ! Instantiate object
    allocate( radXfer_core_t :: radXfer_core_obj )

    ! Initialize grid warehouse
    call tst_config%get( "grids", child_config, Iam )
    radXfer_core_obj%theGridWarehouse_ => grid_warehouse_t( child_config )

    ! Initialize profile warehouse
    call tst_config%get( "profiles", child_config, Iam )
    radXfer_core_obj%theProfileWarehouse_ =>                                  &
        profile_warehouse_t( child_config, radXfer_core_obj%theGridWareHouse_ )

    ! Initialize radXfer xsect warehouse
    call tst_config%get( "cross sections", child_config, Iam )
    radXfer_core_obj%theradXferXsectWarehouse_ =>                             &
        cross_section_warehouse_t( child_config,                              &
                                   radXfer_core_obj%theGridWareHouse_,        &
                                   radXfer_core_obj%theProfileWarehouse_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run( this )

  use tuvx_profile_air, only : profile_air_t

  !> Arguments
  class(radXfer_core_t), intent(inout)  :: this

  !> Local variables
  character(len=*), parameter :: Iam = 'radXfer core run: '

  real(dk)                    :: tstCrossSection
  real(dk), allocatable       :: aCrossSection(:,:)

  class(grid_t), pointer       :: zGrid, lambdaGrid
  class(profile_t), pointer  :: AirProfile, TemperatureProfile
  class(cross_section_t), pointer :: RaylieghCrossSection
  type(string_t)                      :: Handle

    write(*,*) Iam // 'entering'

    !> Get copy of grid
    zGrid => this%theGridWarehouse_%get_grid( "height", "km" )
    call assert( 412238768, zGrid%ncells_ .eq. 120_ik )
    call assert( 412238769, all( zGrid%delta_ .eq. 1._dk ) )

    !> Get copy of wavelength grid
    lambdaGrid => this%theGridWarehouse_%get_grid( "wavelength", "nm" )
    call assert( 412238766, all( lambdaGrid%edge_ > 0._dk ) )
    call assert( 412238767, all( lambdaGrid%delta_ > 0._dk ) )

    !> Get copy of the Air Profile
    AirProfile => this%theProfileWarehouse_%get_profile( "air", "molecule cm-3" )
    call assert( 412238771, all( AirProfile%delta_val_ < 0._dk ) )
    call assert( 412238771, all( AirProfile%layer_dens_ > 0._dk ) )
    write(*,*) ' '
    write(*,*) Iam // 'Air layer density'
    write(*,'(1p10g15.7)') AirProfile%layer_dens_

    write(*,*) ' '
    write(*,*) Iam // 'Air burden density'
    write(*,'(1p10g15.7)') AirProfile%burden_dens_

    !> Get copy of the temperature Profile
    TemperatureProfile => this%theProfileWarehouse_%get_profile( "temperature", "K" )
    call assert( 412238772, all( TemperatureProfile%edge_val_ < 400._dk ) )
    call assert( 412238772, all( TemperatureProfile%edge_val_ > 150._dk ) )
    call assert( 412238773, all( abs(TemperatureProfile%delta_val_) < 20._dk ) )

    !> Get copy of the rayliegh cross section
    Handle = 'air'
    RaylieghCrossSection => this%theradXferXsectWareHouse_%get( Handle )
    aCrossSection = RaylieghCrossSection%calculate( this%theGridWareHouse_, this%theProfileWareHouse_ )
    call assert( 412238776, all( aCrossSection >= 0._dk ) )
    call assert( 412238776, all( aCrossSection < 1._dk ) )

    write(*,*) ' '
    write(*,*) Iam // 'aCrossSection is (',size(aCrossSection,dim=1),' x ',size(aCrossSection,dim=2),')'

    tstCrossSection = aCrossSection(1,1)
    call assert( 412238774, all( aCrossSection(:,1) == tstCrossSection ) )

    write(*,*) ' '
    write(*,*) Iam // 'Rayliegh cross section'
    write(*,'(1p10g15.7)') aCrossSection(1,:)
    call assert( 412238775, all( aCrossSection(1,:) == aCrossSection(zGrid%ncells_,:) ) )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( TemperatureProfile )
    deallocate( AirProfile )
    deallocate( RaylieghCrossSection )
    write(*,*) Iam // 'exiting'

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radXfer core
  subroutine finalize( this )

    !> radXfer core
    type(radXfer_core_t), intent(inout) :: this

  !> Local variables
    character(len=*), parameter :: Iam = 'radXfer core finalize: '

    if( associated( this%theGridWareHouse_ ) ) then
      deallocate( this%theGridWareHouse_ )
    end if

    if( associated( this%theProfileWareHouse_ ) ) then
      deallocate( this%theProfileWareHouse_ )
    end if

    if( associated( this%theradXferXsectWareHouse_ ) ) then
      deallocate( this%theradXferXsectWareHouse_ )
    end if

  end subroutine finalize

end module radXfer_core
