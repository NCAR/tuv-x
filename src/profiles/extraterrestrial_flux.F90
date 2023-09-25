! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_extraterrestrial_flux
  ! Extraterrestrial flux profile type

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t

  implicit none

  public :: profile_extraterrestrial_flux_t

  type, extends(profile_t) :: profile_extraterrestrial_flux_t
  contains
    final     :: finalize
  end type profile_extraterrestrial_flux_t

  !> Constructor
  interface profile_extraterrestrial_flux_t
    module procedure constructor
  end interface profile_extraterrestrial_flux_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result ( this )
    ! Initialize this profile

    use musica_assert,                 only : assert_msg, die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_constants,                only : hc, deltax
    use tuvx_diagnostic_util,          only : diagout
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_grid,                     only : grid_t
    use tuvx_interpolate
    use tuvx_util,                     only : add_point

    type(profile_extraterrestrial_flux_t), pointer :: this ! This f:type:`~tuvx_profile_extraterrestrial_flux/profile_extraterrestrial_flux_t`
    type(config_t), intent(inout)          :: config ! A profile config
    type(grid_warehouse_t), intent(inout)  :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! Local variables
    character(len=*), parameter :: Iam = &
      'extraterrestrial flux Profile initialize: '

    integer, parameter :: Ok = 0
    integer, parameter :: inUnit = 20
    real(dk),    parameter :: bin_edge(0:4) = (/ &
      0.0_dk,150.01_dk,200.07_dk,1000.99_dk,real(huge(0.0_dk),dk) &
    /)
    character(len=*), parameter :: comment = '#!$%*'
    class(grid_t), pointer :: lambdaGrid

    integer  :: istat
    integer  :: fileNdx, nFiles, nBins, nLines, Line
    real(dk), allocatable :: inputGrid(:), inputData(:), tmpinputGrid(:)
    real(dk), allocatable :: interpolatedEtfl(:)
    logical :: found
    character(len=132)             :: InputLine, trimInputLine
    character(len=512)             :: IoMsg
    type(string_t)                 :: defaultInterpolator
    type(string_t), allocatable    :: Filespec(:), Interpolator(:)
    class(interpolator_t), pointer :: theInterpolator
    type(string_t) :: required_keys(5), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "file path"
    required_keys(4) = "interpolator"
    required_keys(5) = "name"
    optional_keys(1) = "enable diagnostics"

    call assert_msg( 389428513,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "extraterrestrial flux profile." )

    allocate( this )

    call config%get( 'enable diagnostics', this%enable_diagnostics, Iam,       &
      default=.false. )

    defaultInterpolator = 'conserving'

    ! Get the configuration settings
    call config%get( 'file path', Filespec, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'interpolator', Interpolator, Iam, found=found )
    nFiles = size(Filespec)
    if( .not. found ) then
      allocate( Interpolator( nFiles ) )
      Interpolator = defaultInterpolator
    endif

    lambdaGrid => grid_warehouse%get_grid( "wavelength", "nm" )
    nBins = lambdaGrid%ncells_

    file_loop: do fileNdx = 1, nFiles
      if(Interpolator( fileNdx ) == "") then
        Interpolator( fileNdx ) = defaultInterpolator
      endif

      ! Does input grid file exist?
      inquire( file = Filespec( fileNdx )%to_char( ), exist = found )
      if( .not. found ) then
        call die_msg( 560768215, "File " &
          // Filespec( fileNdx )%to_char() // " not found" )
      endif

      open( unit = inUnit, file = Filespec( fileNdx )%to_char( ),             &
            iostat = istat, iomsg = IoMsg )

      if( istat /= Ok ) then
        call die_msg( 560768231, "Error opening " // &
          Filespec(fileNdx)%to_char() )
      endif

      ! Determine number of lines in file
      nLines = 0
      do
        read( inUnit, '(a)', iostat = istat ) InputLine
        if( istat == Ok ) then
          trimInputLine = adjustl( InputLine )
          if( verify( trimInputLine(1:1), comment ) /= 0 ) then
            nLines = nLines + 1
            cycle
          endif
        else
          rewind( unit = inUnit )
          exit
        endif
      enddo

      ! Skip the header
      do
        read( inUnit, '(a)', iostat = istat ) InputLine
        if( istat /= Ok ) then
          call die_msg( 560768227, "Error reading " // &
            Filespec( fileNdx )%to_char( ) )
        else
          trimInputLine = adjustl( InputLine )
          if( verify( trimInputLine(1:1), comment ) /= 0 ) then
            exit
          endif
        endif
      enddo

      allocate( inputGrid( nLines ) )
      allocate( inputData( nLines ) )

      ! Read the data
      Line = 1
      do
        trimInputLine = adjustl( InputLine )
        if( verify( trimInputLine(1:1), comment ) /= 0 ) then
          read( InputLine, *, iostat = istat, iomsg = IoMsg )                 &
                inputGrid( Line ), inputData( Line )

          if( istat /= Ok ) then
            call die_msg( 560768229, "Invalid data format in " //             &
              Filespec( fileNdx )%to_char( ) )
          endif
          Line = Line + 1
        endif
        read( inUnit, '(a)', iostat = istat ) InputLine
        if( istat /= Ok ) then
          exit
        endif
      enddo

      close( unit = inUnit )

      ! special handling for neckel.flx
      if( index( Filespec( fileNdx )%to_char( ), 'neckel.flx' ) /= 0 ) then
        allocate( tmpinputGrid, source = inputGrid )
        where( inputGrid < 630._dk )
          tmpinputGrid = inputGrid - 0.5_dk
        elsewhere( inputGrid >= 630._dk .and. inputGrid < 870._dk )
          tmpinputGrid = inputGrid - 1.0_dk
        elsewhere( inputGrid >= 870._dk)
          tmpinputGrid = inputGrid - 2.5_dk
        endwhere
        inputData = 1.e13_dk * hc * inputData / inputGrid
        inputGrid = tmpinputGrid
        inputGrid = [ inputGrid, inputGrid( size( inputGrid ) ) + 2.5_dk ]
        inputData = [ inputData, 0.0_dk ]
        deallocate( tmpinputGrid )
      else
        ! extend inputGrid,inputData to cover model photolysis grid
        call add_point( x = inputGrid, y = inputData,                         &
            xnew = ( 1.0_dk - deltax ) * inputGrid(1), ynew = 0.0_dk )
        call add_point( x = inputGrid, y = inputData,                         &
            xnew = 0.0_dk, ynew = 0.0_dk )
        call add_point( x = inputGrid, y = inputData,                         &
            xnew = ( 1.0_dk + deltax ) * inputGrid( size( inputGrid ) ),      &
            ynew = 0.0_dk )
        call add_point( x = inputGrid, y = inputData, xnew = 1.e38_dk,        &
            ynew = 0.0_dk )
      endif
      ! assign interpolator for this dataset
      select case( Interpolator( fileNdx )%to_char( ) )
        case( 'linear' )
          allocate( interpolator_linear_t :: theInterpolator )
        case( 'conserving' )
          allocate( interpolator_conserving_t :: theInterpolator )
        case( 'fractional source' )
          allocate( interpolator_fractional_source_t :: theInterpolator )
        case( 'fractional target' )
          allocate( interpolator_fractional_target_t :: theInterpolator )
        case default
          call die_msg( 560768275, "interpolator " // &
            Interpolator( fileNdx )%to_char() // " not a valid selection" )
      end select

      ! interpolate from source to model wavelength grid
      interpolatedEtfl = theInterpolator%interpolate( &
        lambdaGrid%edge_, inputGrid, inputData )
      if( .not. allocated( this%mid_val_ ) ) then
        allocate( this%mid_val_,mold=interpolatedEtfl )
        this%mid_val_ = 0.0_dk
      endif

      ! assign interpolated source to model etfl
      where( bin_edge( fileNdx - 1 ) <= lambdaGrid%edge_( : nBins ) .and.     &
             lambdaGrid%edge_( : nBins ) < bin_edge( fileNdx ) )
         this%mid_val_ = interpolatedEtfl
      endwhere

      ! test diagnostics
      associate( enable => this%enable_diagnostics )
      if( index( Filespec( fileNdx )%to_char( ), 'susim' ) /= 0 ) then
        call diagout( 'susim.inputGrid.new', inputGrid, enable )
        call diagout( 'susim.inputData.new', inputData, enable )
        call diagout( 'susim.interpolated.new', interpolatedEtfl, enable )
        call diagout( 'susim.etfl.new', this%mid_val_, enable )
      elseif( index( Filespec( fileNdx )%to_char( ), 'atlas' ) /= 0 ) then
        call diagout( 'atlas.inputGrid.new', inputGrid, enable )
        call diagout( 'atlas.inputData.new', inputData, enable )
        call diagout( 'atlas.interpolated.new', interpolatedEtfl, enable )
        call diagout( 'atlas.etfl.new', this%mid_val_, enable )
      elseif( index( Filespec( fileNdx )%to_char( ), 'neckel' ) /= 0 ) then
        call diagout( 'neckel.inputGrid.new', inputGrid, enable )
        call diagout( 'neckel.inputData.new', inputData, enable )
        call diagout( 'neckel.interpolated.new', interpolatedEtfl, enable )
        call diagout( 'neckel.etfl.new', this%mid_val_, enable )
      elseif( index( Filespec( fileNdx )%to_char( ), 'sao2010' ) /= 0 ) then
        call diagout( 'sao2010.inputGrid.new', inputGrid, enable )
        call diagout( 'sao2010.inputData.new', inputData, enable )
        call diagout( 'sao2010.interpolated.new', interpolatedEtfl, enable )
        call diagout( 'sao2010.etfl.new', this%mid_val_, enable )
      endif
      end associate

      deallocate( inputGrid, inputData )
      deallocate( theInterpolator )

    enddo file_loop

    ! test diagnostics (in original units of W m-2 nm-1)
    call diagout( 'etfl.new', this%mid_val_, this%enable_diagnostics )

    ! convert ET flux values from units of W m-2 nm-1 to photons cm-2 s-1
    this%mid_val_(:) = 1.0e-13 * this%mid_val_ * lambdaGrid%mid_ *            &
                       lambdaGrid%delta_ / hc

    deallocate( lambdaGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_extraterrestrial_flux_t), intent(inout) :: this ! This f:type:`~tuvx_profile_extraterrestrial_flux/profile_extraterrestrial_flux_t`

    if( allocated( this%edge_val_ ) ) deallocate( this%edge_val_ )
    if( allocated( this%mid_val_ ) )  deallocate( this%mid_val_ )
    if( allocated( this%delta_val_ ) ) deallocate( this%delta_val_ )
    if( allocated( this%layer_dens_ ) ) deallocate( this%layer_dens_ )

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_extraterrestrial_flux
