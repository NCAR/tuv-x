! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_cross_section_rayliegh
! Rayleigh scattering cross section from WMO 1985 (originally from
! Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
! An empirical formula for its calculation in the homoshpere, Planet.
! Space Sci., 32, 1467-1468, 1984. 
! `doi:10.1016/10.1016/0032-0633(84)90089-8 
! <https://doi.org/10.1016/0032-0633(84)90089-8>`_

  use tuvx_cross_section,              only : cross_section_t,                &
                                              base_constructor

  implicit none

  private
  public :: cross_section_rayliegh_t

  !> Calculator for rayliegh_cross_section
  type, extends(cross_section_t) :: cross_section_rayliegh_t
  contains
    !> Calculate the cross section
    procedure :: calculate
  end type cross_section_rayliegh_t

  !> Constructor
  interface cross_section_rayliegh_t
    module procedure constructor
  end interface cross_section_rayliegh_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse, profile_warehouse )           &
      result ( this )
    ! Initialize the cross section

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    class(cross_section_t),    pointer       :: this ! This :f:type:`~tuvx_cross_section/cross_section_t`
    type(config_t),            intent(inout) :: config ! Cross section configuration data
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`

    type(string_t) :: required_keys(1), optional_keys(1)

    required_keys(1) = "type"
    optional_keys(1) = "name"
    call assert_msg( 587251674,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "rayliegh cross section." )
    allocate( cross_section_rayliegh_t :: this )
    call base_constructor( this, config, grid_warehouse, profile_warehouse )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function calculate( this, grid_warehouse, profile_warehouse, at_mid_point ) &
      result( cross_section )
    ! Calculate the cross section for a given set of environmental conditions

    use musica_constants,              only : musica_dk
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t
    use tuvx_profile_warehouse,        only : profile_warehouse_t

    real(kind=musica_dk), allocatable              :: cross_section(:,:) ! Calculated cross section
    class(cross_section_rayliegh_t), intent(in)    :: this ! A :f:type:`~tuvx_cross_section_rayliegh/cross_section_rayliegh_t`
    type(grid_warehouse_t),    intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`
    type(profile_warehouse_t), intent(inout) :: profile_warehouse ! A :f:type:`~tuvx_profile_warehouse/profile_warehouse_t`
    logical, optional,           intent(in)    :: at_mid_point ! Flag indicating whether cross-section data should be at mid-points on the wavelength grid.  If this is false or omitted, cross-section data are calculated at interfaces on the wavelength grid.

    ! Local variables
    integer :: colndx, nzdim
    class(grid_t),    pointer     :: zGrid
    class(grid_t),    pointer     :: lambdaGrid
    character(len=*), parameter   :: Iam = 'rayliegh cross section calculate'
    real(musica_dk), allocatable  :: pwr(:), wrk(:)
    real(musica_dk), allocatable  :: wrkCrossSection(:,:)

    zGrid => grid_warehouse%get_grid( this%height_grid_ )
    lambdaGrid => grid_warehouse%get_grid( this%wavelength_grid_ )

    nzdim = zGrid%ncells_ + 1
    if( present( at_mid_point ) ) then
      if( at_mid_point ) then
        nzdim = nzdim - 1
      endif
    endif

    allocate( wrkCrossSection( lambdaGrid%ncells_,nzdim ) )

    allocate( pwr( lambdaGrid%ncells_ ) )
    wrk = 1.e-3_musica_dk * lambdaGrid%mid_
    where( wrk <= 0.55_musica_dk )
      pwr = 3.6772_musica_dk + 0.389_musica_dk * wrk + 0.09426_musica_dk / wrk
    elsewhere
      pwr = 4.04_musica_dk
    endwhere

    wrkCrossSection(:,1) = 4.02e-28_musica_dk / ( wrk )**pwr

    do colndx = 2, nzdim
      wrkCrossSection( :, colndx ) = wrkCrossSection(:,1)
    enddo

    cross_section = transpose( wrkCrossSection )

    deallocate( zGrid )
    deallocate( lambdaGrid )

  end function calculate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_cross_section_rayliegh
