! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The aerosol_micm_radiator module

!> The aerosol_radiator_t type and related functions
!!
module micm_aerosol_radiator_type

  use musica_constants,       only : dk => musica_dk, ik => musica_ik
  use musica_string,          only : string_t
  use micm_abs_radiator_type, only : abs_radiator_t

  implicit none

  private
  public :: aerosol_radiator_t

  !> aerosol radiator type
  type, extends(abs_radiator_t) :: aerosol_radiator_t
  contains
    !> Initialize radiator
    procedure :: initialize
    !> Update radiator for new environmental conditions
    procedure :: upDateState
  end type aerosol_radiator_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Initialize radiator_t object
  subroutine initialize( this, radiator_config, gridWareHouse )

    use musica_assert,        only : die_msg
    use musica_config,        only : config_t
    use micm_grid_warehouse,  only : grid_warehouse_t
    use micm_1d_grid,         only : base_grid_t
    use interpolation
    use tuv_params,           only : nzero, pzero
    use debug,                only : diagout

    !> Arguments
    !> Radiator object
    class(aerosol_radiator_t), intent(inout) :: this
    !> Radiator configuration object
    type(config_t), intent(inout)         :: radiator_config
    !> Grid warehouse
    type(grid_warehouse_t), intent(inout) :: gridWareHouse

    !> Local variables
    character(len=*), parameter   :: Iam = "Aerosol radiator initialize: "
    real(dk), parameter           :: scaling_factor = 550._dk/340._dk

    integer                       :: k, nInputBins, binNdx
    real(dk)                      :: tau550, alpha, wscaling, ODscaling
    real(dk)                      :: coldens
    real(dk), allocatable         :: input_OD(:), rad_OD(:)
    real(dk), allocatable         :: input_SSA(:)
    real(dk), allocatable         :: input_G(:)
    real(dk), allocatable         :: input_zgrid(:)
    real(dk), allocatable         :: winput_SSA(:), winput_G(:)
    type(string_t)                :: Handle
    type(config_t)                :: Aerosol_config
    class(base_grid_t), pointer :: zGrid, lambdaGrid
    class(abs_interpolator_t), pointer :: theInterpolator

    write(*,*) ' '
    write(*,*) Iam,'entering'

    Handle = 'Vertical Z' ; zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength' ; lambdaGrid => gridWareHouse%get_grid( Handle )

!-----------------------------------------------------------------------------
!> Get radiator "Handle"
!-----------------------------------------------------------------------------
    call radiator_config%get( 'Handle', this%handle_, Iam )
    write(*,*) Iam // 'handle = ',this%handle_%to_char()

!> allocate radiator state variables
    allocate( this%state_%layer_OD_(zGrid%ncells_,lambdaGrid%ncells_) )
    allocate( this%state_%layer_SSA_(zGrid%ncells_,lambdaGrid%ncells_) )
    allocate( this%state_%layer_G_(zGrid%ncells_,lambdaGrid%ncells_) )

!> read json config
    call radiator_config%get( "Optical depth", Aerosol_config, Iam )
    call Aerosol_config%get( "Values", input_OD, Iam )
    nInputBins = size(input_OD)
    if( nInputBins > 1 ) then
!> interpolate input OD to state variable
      write(*,*) Iam // 'OD from config'
      write(*,*) Iam // 'size input_OD = ',nInputBins
      write(*,*) Iam // 'input_OD'
      write(*,'(1p10g15.7)') input_OD
      call diagout( 'rawOD.new',input_OD )
      input_OD(:nInputBins-1) = .5_dk*(input_OD(:nInputBins-1)+input_OD(2:))
      write(*,'(1p10g15.7)') input_OD(:nInputBins-1)
      call diagout( 'inpaerOD.new',input_OD(:nInputBins-1) )

      allocate( interp3_t :: theInterpolator )
      input_zgrid = (/ (real(k,dk),k=0,nInputBins-1) /)
      write(*,*) Iam // 'input zgrid'
      write(*,'(1p10g15.7)') input_zgrid
      rad_OD = theInterpolator%interpolate( zGrid%edge_, input_zgrid,input_OD, 1 )
      call diagout( 'cz.aer.new',rad_OD )
      write(*,*) 'size interpolated_OD = ',size(rad_OD)
      write(*,*) 'size interpolated_OD = ',sizeof(rad_OD)
      write(*,*) Iam // 'interpolated OD'
      write(*,'(1p10g15.7)') rad_OD
      do binNdx = 1,lambdaGrid%ncells_
        this%state_%layer_OD_(:,binNdx) = rad_OD
      enddo
    else
      this%state_%layer_OD_ = input_OD(1)
    endif
    
    call radiator_config%get( "Single scattering albedo", Aerosol_config, Iam )
    call Aerosol_config%get( "Values", input_SSA, Iam )
    if( size(input_SSA) > 1 ) then
      call die_msg( 234999,"Aerosol single scattering albedo must be a scalar" )
!> interpolate input SSA to state variable
    else
      winput_SSA = input_OD(:nInputBins-1) * input_SSA(1)
      this%state_%layer_SSA_(:,1) = theInterpolator%interpolate( zGrid%edge_, input_zgrid,winput_SSA, 1 )
      call diagout( 'omz.aer.new',this%state_%layer_SSA_(:,1) )
!     this%state_%layer_SSA_ = input_SSA(1)
      do binNdx = 2,lambdaGrid%ncells_
        this%state_%layer_SSA_(:,binNdx) = this%state_%layer_SSA_(:,1) 
      enddo
      write(*,*) Iam // 'SSA from config'
      write(*,*) Iam // 'size SSA = ',size(input_SSA)
      write(*,*) input_SSA
    endif

    call radiator_config%get( "Asymmetry factor", Aerosol_config, Iam )
    call Aerosol_config%get( "Values", input_G, Iam )
    if( size(input_G) > 1 ) then
      call die_msg( 234999,"Aerosol asymmetry factor must be a scalar" )
!> interpolate input G to state variable
    else
      winput_G = input_OD(:nInputBins-1) * input_G(1)
      this%state_%layer_G_(:,1) = theInterpolator%interpolate( zGrid%edge_, input_zgrid,winput_G, 1 )
      call diagout( 'gz.aer.new',this%state_%layer_G_(:,1) )
!     this%state_%layer_G_ = input_G(1)
      do binNdx = 2,lambdaGrid%ncells_
        this%state_%layer_G_(:,binNdx) = this%state_%layer_G_(:,1) 
      enddo
      write(*,*) Iam // 'G from config'
      write(*,*) Iam // 'size G = ',size(input_G)
      write(*,*) input_G
    endif

    call radiator_config%get( "550 optical depth", tau550, Iam, default=0._dk )
    call radiator_config%get( "Alpha", alpha, Iam, default=1._dk )
    write(*,*) Iam // 'tau550, alpha from config'
    write(*,*) tau550, alpha

    if( tau550 > 0.235_dk ) then
      coldens = max( sum( this%state_%layer_OD_(:,1) ),pzero )
      ODscaling = (tau550/coldens) * scaling_factor**alpha
      do binNdx = 1,lambdaGrid%ncells_
        this%state_%layer_OD_(:,binNdx) = this%state_%layer_OD_(:,binNdx) * ODscaling
      enddo
    endif

    do binNdx = 1,lambdaGrid%ncells_
      wscaling = (340._dk/lambdaGrid%mid_(binNdx))**alpha
      this%state_%layer_OD_(:,binNdx) = this%state_%layer_OD_(:,binNdx) * wscaling
      where( rad_OD > 0._dk )
        this%state_%layer_SSA_(:,binNdx) = this%state_%layer_SSA_(:,binNdx)/rad_OD
        this%state_%layer_G_(:,binNdx)   = this%state_%layer_G_(:,binNdx)/rad_OD
      elsewhere
        this%state_%layer_SSA_(:,binNdx) = 1._dk
        this%state_%layer_G_(:,binNdx)   = 0._dk
      endwhere
    enddo

    write(*,*) Iam // 'layer OD @ lambda = ',lambdaGrid%mid_(1)
    write(*,'(1p10g15.7)') this%state_%layer_OD_(:,1)
    write(*,*) Iam // 'layer OD @ lambda = ',lambdaGrid%mid_(lambdaGrid%ncells_)
    write(*,'(1p10g15.7)') this%state_%layer_OD_(:,lambdaGrid%ncells_)
    write(*,*) ' '
    write(*,*) Iam // 'layer SSA @ lambda = ',lambdaGrid%mid_(1)
    write(*,'(1p10g15.7)') this%state_%layer_SSA_(:,1)
    write(*,*) Iam // 'layer SSA @ lambda = ',lambdaGrid%mid_(lambdaGrid%ncells_)
    write(*,'(1p10g15.7)') this%state_%layer_SSA_(:,lambdaGrid%ncells_)
    write(*,*) ' '
    write(*,*) Iam // 'layer G @ lambda = ',lambdaGrid%mid_(1)
    write(*,'(1p10g15.7)') this%state_%layer_G_(:,1)
    write(*,*) Iam // 'layer G @ lambda = ',lambdaGrid%mid_(lambdaGrid%ncells_)
    write(*,'(1p10g15.7)') this%state_%layer_G_(:,lambdaGrid%ncells_)

    write(*,*) ' '
    write(*,*) Iam,'exiting'

!   stop 'Debugging'

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

    !> Arguments
    !> radiator obj
    class(aerosol_radiator_t), intent(inout) :: this
    !> Grid warehouse
    type(grid_warehouse_t), intent(inout)    :: gridWareHouse
    !> Profile warehouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse
    !> RadXfer cross section warehouse
    type(radXfer_xsect_warehouse_t), intent(inout) :: radXferXsectWareHouse

    !> Local variables
    integer(ik) :: wNdx
    character(len=*), parameter :: Iam = 'Aerosol radiator upDateState: '
    type(string_t)                :: Handle
    class(base_grid_t), pointer :: zGrid
    class(base_grid_t), pointer :: lambdaGrid

    write(*,*) ' '
    write(*,*) Iam,'entering'

    write(*,*) Iam // 'handle = ',this%handle_%to_char()
!-----------------------------------------------------------------------------
!> get specific grids and profiles
!-----------------------------------------------------------------------------
    Handle = 'Vertical Z' ; zGrid => gridWareHouse%get_grid( Handle )
    Handle = 'Photolysis, wavelength' ; lambdaGrid => gridWareHouse%get_grid( Handle )
    write(*,*) Iam // 'nlyr,nbins = ',zGrid%ncells_,lambdaGrid%ncells_

    !> check that radiator state is allocated
    if( .not. allocated( this%state_%layer_OD_ ) ) then
      call die_msg( 2222222,"In radiator%upDateState radiator state not allocate" )
    endif

    write(*,*) ' '
    write(*,*) Iam,'exiting'

  end subroutine upDateState

end module micm_aerosol_radiator_type
