! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_quantum_yield_no2_tint
  !> The no2 tint quantum yield type and related functions

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_quantum_yield,              only : quantum_yield_t

  implicit none

  private
  public :: quantum_yield_no2_tint_t

  type quantum_yield_data_t
    real(dk), allocatable :: temperature(:) ! Temperature grid [K]
    real(dk), allocatable :: deltaT(:)      ! Temperature difference between grid points [K]
    real(dk), allocatable :: array(:,:)     ! Quantum yield parameters (wavelength, temperature)
  contains
    ! Returns the number of bytes required to pack the data
    procedure :: pack_size => data_pack_size
    ! Packs the data onto a character buffer
    procedure :: mpi_pack => data_mpi_pack
    ! Unpacks data from a character buffer
    procedure :: mpi_unpack => data_mpi_unpack
  end type quantum_yield_data_t

  type, extends(quantum_yield_t) :: quantum_yield_no2_tint_t
    ! Calculator for tint quantum yield
    type(quantum_yield_data_t), allocatable :: parameters(:)
  contains
    procedure :: calculate => run
    ! Returns the number of bytes required to pack the quantum yield
    procedure :: pack_size
    ! Packs the quantum yield onto a character buffer
    procedure :: mpi_pack
    ! Unpacks a quantum yield from a character buffer
    procedure :: mpi_unpack
  end type quantum_yield_no2_tint_t

  interface quantum_yield_no2_tint_t
    module procedure constructor
  end interface quantum_yield_no2_tint_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result( this )
    ! Constructor

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_conserving_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(quantum_yield_no2_tint_t),  pointer :: this ! This :f:type:`~tuvx_quantum_yield_no2_tint/quantum_yield_no2_tint_t`
    type(config_t),            intent(inout) :: config ! Quantum yield configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = 'no2 tint quantum yield constructor'
    character(len=*), parameter :: Hdr = 'quantum_yield_'
    real(dk), parameter    :: rZERO = 0.0_dk
    real(dk), parameter    :: rONE  = 1.0_dk

    integer     :: nTemps, nParms
    integer     :: parmNdx, fileNdx, Ndxl, Ndxu
    real(dk)    :: tmp
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical     :: found, monopos
    type(netcdf_t),   allocatable :: netcdf_obj
    type(string_t),   allocatable :: netcdfFiles(:)
    class(grid_t),    pointer     :: lambdaGrid
    type(interpolator_conserving_t) :: interpolator

    allocate( this )

    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ =                                               &
        profile_warehouse%get_ptr( "temperature", "K" )
    this%air_profile_ = profile_warehouse%get_ptr( "air", "molecule cm-3" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! get quantum yield netcdf filespec
    call config%get( 'netcdf files', netcdfFiles, Iam, found = found )
    call assert_msg( 252847407, found,                                        &
                     Iam//'must have at least one netcdf input file' )
    allocate( this%parameters( size( netcdfFiles ) ) )
    file_loop: do fileNdx = 1, size( netcdfFiles )
      allocate( netcdf_obj )
      ! read netcdf file quantum yield data
      call netcdf_obj%read_netcdf_file(                                       &
                               file_path = netcdfFiles( fileNdx )%to_char( ), &
                               variable_name = Hdr )
      nParms = size( netcdf_obj%parameters, dim = 2 )
      call assert_msg( 235314124, nParms >= 2, Iam//'File: '//                &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' array must have 2 or more parameters' )
      associate( Qyield => this%parameters( fileNdx ) )
      ! interpolation temperatures must be in netcdf file
      call assert_msg( 264376965, allocated( netcdf_obj%temperature ),        &
                       Iam//'File: '//                                        &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' does not have interpolation temperatures' )
      Qyield%temperature = netcdf_obj%temperature
      nTemps = size( Qyield%temperature )
      ! must have two or more interpolation temperatures
      call assert_msg( 393489175, nTemps >= 2, Iam//'File: '//                &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' temperature array has < 2 entries' )
      call assert_msg( 167399246, nTemps >= nParms, Iam//'File: '//           &
                       trim( netcdfFiles( fileNdx )%to_char( ) )//            &
                       ' temperature array < number parameters' )
      Qyield%deltaT = Qyield%temperature( 2 : nParms ) -                      &
                        Qyield%temperature( 1 : nParms - 1 )
      monopos = all( Qyield%deltaT > rZERO )
      if( .not. monopos ) then
        call assert_msg( 606144012, .not. any( Qyield%deltaT > rZERO ),       &
                         Iam//'File: '//                                      &
                         trim( netcdfFiles( fileNdx )%to_char( ) )//          &
                         ' temperature array not monotonic' )
        do Ndxl = 1, nParms / 2
          Ndxu = nParms - Ndxl + 1
          tmp = Qyield%temperature( Ndxl )
          Qyield%temperature( Ndxl ) = Qyield%temperature( Ndxu )
          Qyield%temperature( Ndxu ) = tmp
          data_parameter = netcdf_obj%parameters(:,Ndxl)
          netcdf_obj%parameters( :, Ndxl ) =                                  &
              netcdf_obj%parameters( :, Ndxu )
          netcdf_obj%parameters( :, Ndxu ) = data_parameter
        enddo
        Qyield%deltaT = Qyield%temperature( 2 : nParms ) -                    &
                          Qyield%temperature( 1 : nParms - 1 )
      endif
      ! interpolate from data to model wavelength grid
      if( allocated( netcdf_obj%wavelength ) ) then
        if( .not. allocated( this%parameters( fileNdx )%array) ) then
          allocate( this%parameters( fileNdx )%array(                         &
                                              lambdaGrid%ncells_, nParms ) )
        endif
        do parmNdx = 1, nParms
          data_lambda    = netcdf_obj%wavelength
          data_parameter = netcdf_obj%parameters( :, parmNdx )
          call this%add_points( config, data_lambda, data_parameter )
          this%parameters( fileNdx )%array( :, parmNdx ) =                    &
                interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                          x_source = data_lambda,             &
                                          y_source = data_parameter )
        enddo
      else
        this%parameters( fileNdx )%array = netcdf_obj%parameters
      endif
      end associate
      deallocate( netcdf_obj )
    enddo file_loop

    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the quantum yield for the environmental conditions
  function run( this, grid_warehouse, profile_warehouse )                     &
      result( quantum_yield )

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(quantum_yield_no2_tint_t), intent(in)    :: this ! This :f:type:`~tuvx_quantum_yield_no2_tint/quantum_yield_no2_tint_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    real(kind=dk), allocatable               :: quantum_yield(:,:) ! Calculated quantum_yield

    !> Local variables
    character(len=*), parameter :: Iam = 'no2 tint quantum yield calculate'
    integer     :: nTemp
    integer     :: fileNdx, tNdx, vertNdx
    real(dk)    :: Tadj, Tstar
    real(dk),         allocatable :: WrkQuantumYield(:,:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    allocate( wrkQuantumYield( lambdaGrid%ncells_, zGrid%ncells_ + 1 ) )
    wrkQuantumYield = 0.0_dk

    do fileNdx = 1, size( this%parameters )
      associate( Temp => this%parameters( fileNdx )%temperature,              &
                 wrkQyield => this%parameters( fileNdx ) )
      nTemp = size( Temp )
      do vertNdx = 1, zGrid%ncells_ + 1
        Tadj  = mdlTemperature%edge_val_( vertNdx )
        do tNdx = 2, nTemp
          if( Tadj <= Temp( tNdx ) ) then
            exit
          endif
        enddo
        tndx = min( nTemp, tNdx ) - 1
        Tstar = ( Tadj - Temp( tNdx ) ) / wrkQyield%deltaT( tNdx )
        WrkQuantumYield( :, vertNdx ) = WrkQuantumYield( :, vertNdx ) +       &
                    wrkQyield%array( :, tNdx ) +                              &
                    Tstar * ( wrkQyield%array( :, tNdx + 1 ) -                &
                              wrkQyield%array( :, tNdx ) )
      enddo
      end associate
    enddo

    quantum_yield = transpose( max( WrkQuantumYield, 0.0_dk ) )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the quantum yield

    use musica_mpi,                    only : musica_mpi_pack_size

    class(quantum_yield_no2_tint_t), intent(in) :: this ! quantum yield to pack
    integer,                         intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: i_data

    pack_size = this%quantum_yield_t%pack_size( comm )                        &
                + musica_mpi_pack_size( allocated( this%parameters ), comm )
    if( allocated( this%parameters ) ) then
      pack_size = pack_size                                                   &
                  + musica_mpi_pack_size( size( this%parameters ), comm )
      do i_data = 1, size( this%parameters )
        pack_size = pack_size + this%parameters( i_data )%pack_size( comm )
      end do
    end if
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the quantum yield onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(quantum_yield_no2_tint_t), intent(in)    :: this      ! quantum yield to pack
    character,                       intent(inout) :: buffer(:) ! memory buffer
    integer,                         intent(inout) :: position  ! current buffer position
    integer,                         intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_data

    prev_pos = position
    call this%quantum_yield_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, allocated( this%parameters ),     &
                          comm )
    if( allocated( this%parameters ) ) then
      call musica_mpi_pack( buffer, position, size( this%parameters ), comm )
      do i_data = 1, size( this%parameters )
        call this%parameters( i_data )%mpi_pack( buffer, position, comm )
      end do
    end if
    call assert( 284574988, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a quantum yield from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(quantum_yield_no2_tint_t), intent(out)   :: this      ! quantum yield to be unpacked
    character,                       intent(inout) :: buffer(:) ! memory buffer
    integer,                         intent(inout) :: position  ! current buffer position
    integer,                         intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos, i_data, n_data
    logical :: alloced

    prev_pos = position
    call this%quantum_yield_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, alloced, comm )
    if( alloced ) then
      call musica_mpi_unpack( buffer, position, n_data, comm )
      allocate( this%parameters( n_data ) )
      do i_data = 1, n_data
        call this%parameters( i_data )%mpi_unpack( buffer, position, comm )
      end do
    end if
    call assert( 247628485, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function data_pack_size( this, comm ) result( pack_size )
    ! Returns the number of bytes required to pack the data

    use musica_mpi,                    only : musica_mpi_pack_size

    class(quantum_yield_data_t), intent(in) :: this ! data to be packed
    integer,                     intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = musica_mpi_pack_size( this%temperature, comm ) +              &
                musica_mpi_pack_size( this%deltaT,      comm ) +              &
                musica_mpi_pack_size( this%array,       comm )
#else
    pack_size = 0
#endif

  end function data_pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_mpi_pack( this, buffer, position, comm )
    ! Packs the data onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(quantum_yield_data_t), intent(in)    :: this      ! data to be packed
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_pack( buffer, position, this%temperature, comm )
    call musica_mpi_pack( buffer, position, this%deltaT,      comm )
    call musica_mpi_pack( buffer, position, this%array,       comm )
    call assert( 768462639, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine data_mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine data_mpi_unpack( this, buffer, position, comm )
    ! Unpacks data from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(quantum_yield_data_t), intent(out)   :: this      ! data to be unpacked
    character,                   intent(inout) :: buffer(:) ! memory buffer
    integer,                     intent(inout) :: position  ! current buffer position
    integer,                     intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call musica_mpi_unpack( buffer, position, this%temperature, comm )
    call musica_mpi_unpack( buffer, position, this%deltaT,      comm )
    call musica_mpi_unpack( buffer, position, this%array,       comm )
    call assert( 811412867, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine data_mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_quantum_yield_no2_tint
