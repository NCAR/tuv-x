! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp radiative transfer module
module radiator_core

  use musica_config,    only : config_t
  use musica_string,    only : string_t
  use musica_assert,    only : assert
  use musica_constants, only : ik => musica_ik, dk => musica_dk
  use tuvx_grid_warehouse, only : grid_warehouse_t
  use tuvx_grid,        only : grid_t
  use tuvx_profile_warehouse, only : Profile_warehouse_t
  use tuvx_profile,           only : profile_t
  use tuvx_cross_section_warehouse,        only : cross_section_warehouse_t
  use tuvx_cross_section, only : cross_section_t
  use tuvx_radiator_warehouse,    only : radiator_warehouse_t
  use tuvx_radiator,     only : radiator_t
  use tuvx_radiator,     only : radiator_state_t
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif
  implicit none

  private
  public :: radiator_core_t

  type :: radiator_core_t
    type(grid_warehouse_t), pointer          :: theGridWarehouse_ => null( )
    type(Profile_warehouse_t), pointer       ::                               &
        theProfileWarehouse_ => null( )
    type(cross_section_warehouse_t), pointer ::                               &
        theradXferXsectWarehouse_ => null( )
    type(radiator_warehouse_t), pointer      ::                               &
        theRadiatorWarehouse_ => null( )
  contains
    procedure :: test => run
    final     :: finalize
  end type radiator_core_t

  interface radiator_core_t
    module procedure constructor
  end interface radiator_core_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config_flsp ) result( radiator_core_obj )

    !> Arguments
    type(string_t), intent(in)         :: config_flsp

    class(radiator_core_t), pointer :: radiator_core_obj

    !> Local variables
    character(len=*), parameter :: Iam = 'radiator core constructor: '
    type(config_t)              :: tst_config, child_config

    ! main configuration -> config type
    call tst_config%from_file( config_flsp%to_char() )

    ! Instantiate object
    allocate( radiator_core_t :: radiator_core_obj )

    ! Initialize grid warehouse
    call tst_config%get( "grids", child_config, Iam )
    radiator_core_obj%theGridWarehouse_ => grid_warehouse_t( child_config )

    ! Initialize profile warehouse
    call tst_config%get( "profiles", child_config, Iam )
    radiator_core_obj%theProfileWarehouse_ =>                                 &
        profile_warehouse_t( child_config,                                    &
                             radiator_core_obj%theGridWareHouse_ )

    ! Initialize radXfer xsect warehouse
    call tst_config%get( "cross sections", child_config, Iam )
    radiator_core_obj%theradXferXsectWarehouse_ =>                            &
        cross_section_warehouse_t( child_config,                              &
                                   radiator_core_obj%theGridWareHouse_,       &
                                   radiator_core_obj%theProfileWarehouse_ )

    ! Initialize radiator warehouse
    call tst_config%get( "radiators", child_config, Iam )
    radiator_core_obj%theRadiatorWarehouse_ =>                                &
        radiator_warehouse_t( child_config,                                   &
                              radiator_core_obj%theGridWareHouse_,            &
                              radiator_core_obj%theProfileWarehouse_,         &
                              radiator_core_obj%theradXferXsectWarehouse_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run( this )

    use musica_mpi
    use tuvx_radiator_warehouse,       only : warehouse_iterator_t

    class(radiator_core_t), intent(inout)  :: this

    ! Local variables
    character(len=*), parameter :: Iam = 'radiator core run: '

    real(dk)                    :: tstCrossSection
    real(dk), allocatable       :: aCrossSection(:,:)

    class(grid_t), pointer       :: zGrid, lambdaGrid
    class(profile_t), pointer       :: AirProfile, TemperatureProfile
    class(cross_section_t), pointer :: RaylieghCrossSection
    class(radiator_t), pointer     :: RaylieghRadiator
    class(radiator_t), pointer     :: aRadiator
    type(warehouse_iterator_t), pointer :: iter
    type(string_t)                      :: Handle
    type(radiator_state_t), allocatable :: RadiatorState
    type(radiator_warehouse_t), pointer :: radiators
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    logical :: found
    integer, parameter :: comm = MPI_COMM_WORLD

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
    deallocate( RaylieghCrossSection )

    write(*,*) ' '
    write(*,*) Iam // 'aCrossSection is (',size(aCrossSection,dim=1),' x ',size(aCrossSection,dim=2),')'

    tstCrossSection = aCrossSection(1,1)
    call assert( 412238774, all( aCrossSection(:,1) == tstCrossSection ) )

    write(*,*) ' '
    write(*,*) Iam // 'Rayliegh cross section'
    write(*,'(1p10g15.7)') aCrossSection(1,:)
    call assert( 412238775, all( aCrossSection(1,:) == aCrossSection(zGrid%ncells_,:) ) )

    ! Get copy of the rayliegh radiator and test MPI functions
    if( musica_mpi_rank( comm ) == 0 ) then
      Handle = 'air'
      RaylieghRadiator => this%theRadiatorWarehouse_%get_radiator( Handle )
      call RaylieghRadiator%update_state( this%theGridWareHouse_,             &
                                          this%theProfileWareHouse_,          &
                                          this%theradXferXsectWareHouse_ )
      radiators => this%theRadiatorWarehouse_
      pack_size = RaylieghRadiator%pack_size( comm ) +                             &
                  radiators%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call RaylieghRadiator%mpi_pack( buffer, pos , comm )
      call radiators%mpi_pack( buffer, pos , comm )
      call assert( 994121788, pos <= pack_size )
    end if

    call musica_mpi_bcast( pack_size , comm )
    if( musica_mpi_rank( comm ) .ne. 0 ) allocate( buffer( pack_size ) )
    call musica_mpi_bcast( buffer , comm )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      pos = 0
      allocate( RaylieghRadiator )
      allocate( radiators        )
      call RaylieghRadiator%mpi_unpack( buffer, pos , comm )
      call radiators%mpi_unpack( buffer, pos , comm )
      call assert( 761314313, pos <= pack_size )
    end if

     ! Evaluate radiator state
    call assert( 312238775, all( RaylieghRadiator%state_%layer_OD_ >= 0._dk ) )
    write(*,*) Iam // 'layer_OD_ is (',size(RaylieghRadiator%state_%layer_OD_,dim=1),' x ', &
                      size(RaylieghRadiator%state_%layer_OD_,dim=2),')'
    call assert( 312238776, all( RaylieghRadiator%state_%layer_SSA_ >= 0._dk ) )
    write(*,*) Iam // 'layer_SSA_ is (',size(RaylieghRadiator%state_%layer_SSA_,dim=1),' x ', &
                      size(RaylieghRadiator%state_%layer_SSA_,dim=2),')'
    call assert( 312238777, all( RaylieghRadiator%state_%layer_G_ >= 0._dk ) )
    call assert( 312238778, all( RaylieghRadiator%state_%layer_SSA_ == 1._dk ) )
    write(*,*) Iam // 'layer_G_ is (',size(RaylieghRadiator%state_%layer_G_,dim=1),' x ', &
                      size(RaylieghRadiator%state_%layer_G_,dim=2),')'
    write(*,*) ' '
    write(*,*) Iam // 'Air radiator OD @ top of model'
    write(*,'(1p10g15.7)') RaylieghRadiator%state_%layer_OD_(zGrid%ncells_,:)
    write(*,*) ' '
    write(*,*) Iam // 'Air radiator OD @ ground'
    write(*,'(1p10g15.7)') RaylieghRadiator%state_%layer_OD_(1,:)

    write(*,*) Iam // 'Before radiator iterator test'

    !> Test warehouse iterator and MPI passed warehouse
    found = .false.
    iter => radiators%get_iterator()
    do while( iter%next() )
      nullify( aRadiator )
      aRadiator => radiators%get_radiator( iter )
      if( radiators%name( iter ) == 'air' ) then
        found = .true.
        call assert( 356694336, all( aRadiator%state_%layer_OD_ >= 0._dk ) )
        call assert( 869070582, all( aRadiator%state_%layer_SSA_ >= 0._dk ) )
        call assert( 133963180, all( aRadiator%state_%layer_G_ >= 0._dk ) )
        call assert( 298855777, all( aRadiator%state_%layer_SSA_ == 1._dk ) )
      end if
    enddo
    deallocate(iter)
    call assert( 771640568, found )

    if( musica_mpi_rank( comm ) .ne. 0 ) then
      deallocate( RaylieghRadiator )
      deallocate( radiators        )
    end if
    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( AirProfile )
    deallocate( TemperatureProfile )
    write(*,*) Iam // 'exiting'

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radXfer core
  subroutine finalize( this )

    !> radXfer core
    type(radiator_core_t), intent(inout) :: this

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

    if( associated( this%theRadiatorWareHouse_ ) ) then
      deallocate( this%theradiatorWareHouse_ )
    end if

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module radiator_core
