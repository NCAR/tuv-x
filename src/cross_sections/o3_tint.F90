! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_o3_tint
! Calculates the cross section for ozone with temperature
! interpolation

  ! Including musica_config at the module level to avoid an ICE
  ! with Intel 2022.1 compiler
  use musica_config,                   only : config_t
  use musica_constants,                only : dk => musica_dk
  use tuvx_cross_section,              only : cross_section_t

  implicit none

  private
  public :: cross_section_o3_tint_t

  !> Calculator for o3_tint_cross_section
  type, extends(cross_section_t) :: cross_section_o3_tint_t
    real(dk) :: v185(1), v195(1), v345(1)
  contains
    !> Calculate the cross section
    procedure :: calculate
    !> Returns the number of bytes required to pack the cross section onto
    !> a buffer
    procedure :: pack_size
    !> Packs the cross section onto a character buffer
    procedure :: mpi_pack
    !> Unpacks a cross section from a character buffer into the object
    procedure :: mpi_unpack
    !> refraction
    procedure :: refraction
  end type cross_section_o3_tint_t

  !> Constructor
  interface cross_section_o3_tint_t
    module procedure constructor
  end interface cross_section_o3_tint_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize cross_section_o3_tint_t object

    use musica_assert,                 only : assert_msg, die_msg
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_interpolate,              only : interpolator_t
    use tuvx_netcdf,                   only : netcdf_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    type(cross_section_o3_tint_t), pointer       :: this ! A :f:type:`~tuvx_cross_section_o3_tint/cross_section_o3_tint_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    ! Local variables
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk
    real(dk), parameter  :: refracDensity = 2.45e19_dk
    real(dk), parameter  :: w185 = 185._dk
    real(dk), parameter  :: w195 = 195._dk
    real(dk), parameter  :: w345 = 345._dk

    character(len=*), parameter :: Iam = 'o3 tint cross section constructor'
    character(len=*), parameter :: Hdr = 'cross_section_'

    integer :: parmNdx, fileNdx, Ndxl, Ndxu
    integer :: nParms, nTemps
    real(dk)              :: tmp
    real(dk), allocatable :: refracNdx(:)
    real(dk), allocatable :: data_lambda(:)
    real(dk), allocatable :: data_parameter(:)
    logical :: found, monopos
    character(len=:), allocatable :: msg
    type(netcdf_t),   allocatable :: netcdf_obj
    type(config_t) :: netcdf_files, netcdf_file
    class(iterator_t), pointer :: iter
    type(string_t) :: file_path
    class(grid_t),    pointer     :: lambdaGrid
    type(config_t) :: interpolator_config
    class(interpolator_t), pointer :: interpolator
    type(string_t) :: required_keys(2), optional_keys(3)

    required_keys(1) = "type"
    required_keys(2) = "netcdf files"
    optional_keys(1) = "lower extrapolation"
    optional_keys(2) = "upper extrapolation"
    optional_keys(3) = "name"
    call assert_msg( 988762814,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "o3 tint cross section." )

    allocate( this )

    ! Get model wavelength grids
    this%wavelength_grid_ = grid_warehouse%get_ptr( "wavelength", "nm" )
    this%height_grid_ = grid_warehouse%get_ptr( "height", "km" )
    this%temperature_profile_ =                                               &
        profile_warehouse%get_ptr( "temperature", "K" )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    ! get cross section netcdf filespec
    call config%get( 'netcdf files', netcdf_files, Iam )
    iter => netcdf_files%get_iterator( )

    ! ensure the correct number of NetCDF Files are present
    call assert_msg( 430013090, netcdf_files%number_of_children( ) .ge. 2,    &
                     Iam//' must have at least two netcdf input files' )
    allocate( this%cross_section_parms( netcdf_files%number_of_children( ) ) )
    allocate( refracNdx(0) )

    fileNdx = 0
    do while( iter%next( ) )
      call netcdf_files%get( iter, netcdf_file, Iam )
      fileNdx = fileNdx + 1

      allocate( netcdf_obj )
      ! read netcdf cross section parameters
      call netcdf_file%get( "file path", file_path, Iam )
      call netcdf_obj%read_netcdf_file( file_path = file_path%to_char( ),     &
                                        variable_name = Hdr )
      nParms = size( netcdf_obj%parameters, dim = 2 )
      ! must have at least one parameter
      call assert_msg( 469152250, nParms >= 1,                                &
                       Iam//'File: '//file_path//                             &
                       '  array must have 1 or more parameters' )
      call assert_msg( 744909953, nParms >= 4 .or. fileNdx .ne. 2,            &
                       Iam//'File: '//file_path//                             &
                       '  array must have 4 or more parameters' )

      ! refraction index
      call assert_msg( 218168625, allocated( netcdf_obj%wavelength ),         &
                       "Missing wavelengths in O3 temperature integrated "//  &
                       "cross section data file '"//file_path//"'" )
      refracNdx = this%refraction( netcdf_obj%wavelength, refracDensity )
      netcdf_obj%wavelength = refracNdx * netcdf_obj%wavelength

      associate( Xsection => this%cross_section_parms( fileNdx ) )

      ! interpolation temperatures must be in netcdf file
      call assert_msg( 724757315, allocated( netcdf_obj%temperature ),        &
                       Iam//'File: '//file_path//                             &
                       ' must have interpolation temperatures' )
      Xsection%temperature = netcdf_obj%temperature
      nTemps = size( Xsection%temperature )
      if( nTemps > 1 ) then
        Xsection%deltaT = Xsection%temperature( 2 : nParms )                  &
                          - Xsection%temperature( 1 : nParms - 1 )
        monopos = all( Xsection%deltaT > rZERO )
        if( .not. monopos ) then
          if( any( Xsection%deltaT > rZERO ) ) then
            write(msg,*) Iam//'File: '//file_path//                           &
                         '  temperature array not monotonic'
            call die_msg( 175583000, msg )
          endif
          do Ndxl = 1, nParms / 2
            Ndxu = nParms - Ndxl + 1
            tmp = Xsection%temperature( Ndxl )
            Xsection%temperature( Ndxl ) = Xsection%temperature( Ndxu )
            Xsection%temperature( Ndxu ) = tmp
            data_parameter = netcdf_obj%parameters( :, Ndxl )
            netcdf_obj%parameters( :, Ndxl ) =                                &
                netcdf_obj%parameters( :, Ndxu )
            netcdf_obj%parameters( :, Ndxu ) = data_parameter
          enddo
          Xsection%deltaT = Xsection%temperature( 2 : nParms )                &
                            - Xsection%temperature( 1 : nParms - 1 )
        endif
      endif

      ! interpolate from data to model wavelength grid
      if( allocated( netcdf_obj%wavelength ) ) then
        if( .not. allocated( Xsection%array ) ) then
          allocate( Xsection%array( lambdaGrid%ncells_, nParms ) )
        endif
        call netcdf_file%get( "interpolator", interpolator_config, Iam,       &
                              found = found )
        if( .not. found ) then
          call interpolator_config%empty( )
          call interpolator_config%add( "type", "conserving", Iam )
        end if
        interpolator => interpolator_t( interpolator_config )
        do parmNdx = 1, nParms
          data_lambda    = netcdf_obj%wavelength
          data_parameter = netcdf_obj%parameters( :, parmNdx )
          call this%add_points( config, data_lambda, data_parameter )
          Xsection%array( :, parmNdx ) =                                      &
                interpolator%interpolate( x_target = lambdaGrid%edge_,        &
                                          x_source = data_lambda,             &
                                          y_source = data_parameter )
        enddo
        deallocate( interpolator )
      else
        Xsection%array = netcdf_obj%parameters
      endif
      end associate
      deallocate( netcdf_obj )
    enddo
    this%v185 = this%refraction( (/ w185 /), refracDensity ) * w185
    this%v195 = this%refraction( (/ w195 /), refracDensity ) * w195
    this%v345 = this%refraction( (/ w345 /), refracDensity ) * w345
    this%cross_section_parms(2)%array(:,4) =                                  &
        this%cross_section_parms(1)%array(:,1)

    deallocate( iter       )
    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile,                  only : profile_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=dk), allocatable                    :: cross_section(:,:) ! Calculated cross section
    class(cross_section_o3_tint_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_o3_tint/cross_section_o3_tint_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    character(len=*), parameter :: Iam = 'o3 tint cross section calculate'
    real(dk), parameter :: rZERO = 0.0_dk
    real(dk), parameter :: rONE  = 1.0_dk
    integer :: k
    integer :: nTemp
    integer :: fileNdx, tNdx, wNdx, nzdim
    real(dk)    :: Tadj, Tstar
    real(dk),         allocatable :: modelTemp(:)
    real(dk),         allocatable :: wrkCrossSection(:,:)
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    class(profile_t), pointer     :: mdlTemperature

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )
    mdlTemperature =>                                                         &
        profile_warehouse%get_profile( this%temperature_profile_ )

    ! temperature at model cell midpoint or edge
    nzdim     = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
        modelTemp = mdlTemperature%mid_val_
      else
        modelTemp = mdlTemperature%edge_val_
      endif
    else
      modelTemp = mdlTemperature%edge_val_
    endif

    allocate( wrkCrossSection( lambdaGrid%ncells_, nzdim ) )
    wrkCrossSection = rZERO

vert_loop:                                                                    &
    do k = 1, nzdim
lambda_loop:                                                                  &
    do wNdx = 1, lambdaGrid%ncells_
      associate( lambda => lambdaGrid%edge_( wNdx ) )
        if( lambda < this%v185(1) ) then
          fileNdx = 3
        elseif( this%v185(1) <= lambda .and. lambda < this%v195(1) ) then
          fileNdx = 4
        elseif( this%v195(1) <= lambda .and. lambda < this%v345(1) ) then
          fileNdx = 2
        else
          wrkCrossSection( wNdx, k ) = wrkCrossSection( wNdx, k )             &
                                + this%cross_section_parms(1)%array( wNdx, 1 )
          cycle lambda_loop
        endif
      end associate
      associate( dataTemp => this%cross_section_parms( fileNdx )%temperature, &
                 wrkXsect => this%cross_section_parms( fileNdx ) )
      nTemp = size( dataTemp )
      Tadj  = min( max( modelTemp( k ),dataTemp(1) ),dataTemp( nTemp ) )
      do tNdx = 2, nTemp
        if( Tadj <= dataTemp( tNdx ) ) then
          exit
        endif
      enddo
      tNdx = tNdx - 1
      Tstar = ( Tadj - dataTemp( tNdx ) ) / wrkXsect%deltaT( tNdx )
      wrkCrossSection( wNdx, k ) = wrkCrossSection( wNdx, k )                 &
                                   + wrkXsect%array( wNdx, tNdx )             &
                                 + Tstar * ( wrkXsect%array( wNdx, tNdx + 1 ) &
                                             - wrkXsect%array( wNdx, tNdx ) )
      end associate
    enddo lambda_loop
    enddo vert_loop

    cross_section = transpose( wrkCrossSection )

    deallocate( zGrid )
    deallocate( lambdaGrid )
    deallocate( mdlTemperature )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the size of a character buffer required to pack the cross
    ! section

    use musica_mpi,                    only : musica_mpi_pack_size

    class(cross_section_o3_tint_t), intent(in) :: this ! cross section to be packed
    integer,                        intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%cross_section_t%pack_size( comm ) +                      &
                musica_mpi_pack_size( this%v185(1), comm ) +                  &
                musica_mpi_pack_size( this%v195(1), comm ) +                  &
                musica_mpi_pack_size( this%v345(1), comm )
#else
    pack_size = this%cross_section_t%pack_size( comm )
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the cross section onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(cross_section_o3_tint_t), intent(in)    :: this      ! cross section to be packed
    character,                      intent(inout) :: buffer(:) ! memory buffer
    integer,                        intent(inout) :: position  ! current buffer position
    integer,                        intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%cross_section_t%mpi_pack( buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%v185(1), comm )
    call musica_mpi_pack( buffer, position, this%v195(1), comm )
    call musica_mpi_pack( buffer, position, this%v345(1), comm )
    call assert( 582324821, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks a cross section from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(cross_section_o3_tint_t), intent(out)   :: this      ! cross section to be unpacked
    character,                      intent(inout) :: buffer(:) ! memory buffer
    integer,                        intent(inout) :: position  ! current buffer position
    integer,                        intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%cross_section_t%mpi_unpack( buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%v185(1), comm )
    call musica_mpi_unpack( buffer, position, this%v195(1), comm )
    call musica_mpi_unpack( buffer, position, this%v345(1), comm )
    call assert( 560718944, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function refraction( this, wavelength, atmDensity ) result( refrac )
    ! Calculate the refractive index for standard air
    !
    ! (dry air at 15 deg. C, 101.325 kPa, 0.03% CO2)
    ! from CRC Handbook, originally from Edlen, B., Metrologia, 2, 71, 1966.
    ! valid from 200 nm to 2000 nm
    ! beyond this range, use constant value

    class(cross_section_o3_tint_t), intent(in) :: this ! A :f:type:`~tuvx_cross_section_o3_tint/cross_section_o3_tint_t`
    real(dk),                       intent(in) :: atmDensity ! air density [molecule cm-3]
    real(dk),                       intent(in) :: wavelength(:) ! Vacuum wavelength [nm]
    real(dk)                                   :: refrac( size( wavelength ) ) ! Refractive index of standard air

    real(dk), parameter :: rONE = 1.0_dk
    real(dk), parameter :: divisor = 2.69e19_dk * 273.15_dk / 288.15_dk
    integer :: wNdx
    real(dk) :: wrk, sig, sigsq


    do wNdx = 1, size( wavelength )
      wrk = min( max( wavelength( wNdx ), 200._dk ), 2000._dk )
      sig = 1.e3_dk / wrk
      sigsq = sig * sig

      wrk = 8342.13_dk + 2406030._dk / ( 130._dk - sigsq )                    &
          + 15997._dk / ( 38.9_dk - sigsq )

      ! adjust to local air density
      wrk = wrk * atmDensity / divisor

      ! index of refraction
      refrac(wNdx) = rONE + 1.e-8_dk * wrk
    enddo

    end function refraction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_o3_tint
