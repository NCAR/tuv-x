! Copyright (C) 2024 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

program test_la_srb
  !> Tests against Doug's LA and SR band calculations

  use musica_assert,                   only : assert, almost_equal
  use musica_constants,                only : dk => musica_dk
  use musica_config,                   only : config_t
  use tuvx_cross_section,              only : cross_section_t
  use tuvx_grid_warehouse,             only : grid_warehouse_t
  use tuvx_profile_warehouse,          only : profile_warehouse_t
  use tuvx_la_sr_bands,                only : la_sr_bands_t
  use tuvx_test_utils,                 only : check_values

  implicit none

  character(len=*), parameter :: my_name = "LUT LA/SRB test"
  character(len=*), parameter :: conf_l = 'test/data/la_srb_bands.config.json'
  type(config_t) :: grid_config, la_config, profile_config, o2_config
  class(grid_warehouse_t), pointer :: grids => null( )
  class(profile_warehouse_t), pointer :: profiles => null( )
  class(la_sr_bands_t), pointer :: la_sr_bands => null( )
  class(cross_section_t), pointer :: o2_cross_section => null( )
  character, allocatable :: buffer(:)
  
  call la_config%from_file( conf_l )
  call la_config%get( "grids", grid_config, my_name )
  call la_config%get( "profiles", profile_config, my_name )
  call la_config%get( "O2 cross section", o2_config, my_name )

  grids => grid_warehouse_t( grid_config )
  profiles => profile_warehouse_t( profile_config, grids )
  la_sr_bands => la_sr_bands_t( la_config, grids, profiles )
  o2_cross_section => cross_section_t( o2_config, grids, profiles )

  call compare_o2_cross_sections( la_sr_bands, grids, profiles )

  deallocate( grids )
  deallocate( profiles )
  deallocate( la_sr_bands )
  deallocate( o2_cross_section )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine compare_o2_cross_sections( la_sr_bands, grids, profiles )

    use tuvx_grid,                     only : grid_t
    use tuvx_profile,                  only : profile_t
    use tuvx_spherical_geometry,       only : spherical_geometry_t

    class(la_sr_bands_t),      intent(inout) :: la_sr_bands
    class(grid_warehouse_t),    intent(inout) :: grids
    class(profile_warehouse_t), intent(inout) :: profiles

    character(len=80) :: file_path
    real(dk), allocatable :: air_vertical_column(:), air_slant_column(:)
    real(dk), allocatable :: o2_vertical_column(:), o2_slant_column(:)
    real(dk), allocatable :: tuvx_o2_optical_depth(:,:),                      &
                             tuvx_o2_cross_section(:,:)
    real, dimension(151) :: lut_heights, lut_temperature,                     &
                            lut_air_vertical_column, lut_air_slant_column,    &
                            lut_o2_column
    real, dimension(700) :: lut_wavelength_edges, lut_wavelength_centers,     &
                            lut_o2_base_cross_section
    real, dimension(151,700) :: lut_o2_cross_section, lut_o2_optical_depth
    class(grid_t), pointer :: heights, wavelengths, lut_wavelengths
    class(spherical_geometry_t), pointer :: geometry
    class(profile_t), pointer :: air, o2, temperature
    real(dk), allocatable :: solar_zenith_angles(:)
    integer :: i_sza, i_height, i_wl, n_heights, n_wavelengths, i_output_height
    integer :: output_heights(4) = (/ 1, 50, 100, 150 /)
    real(dk) :: rel_tol

    heights         => grids%get_grid( "height", "km" )
    wavelengths     => grids%get_grid( "wavelength", "nm" )
    lut_wavelengths => grids%get_grid( "LUT wavelength", "nm" )
    air             => profiles%get_profile( "air", "molecule cm-3" )
    o2              => profiles%get_profile( "O2", "molecule cm-3" )
    temperature     => profiles%get_profile( "temperature", "K" )
    geometry        => spherical_geometry_t( grids )
    solar_zenith_angles = (/ 0.0_dk, 13.2_dk, 45.0_dk, 87.3_dk, 90.0_dk /)
    allocate( air_vertical_column( air%ncells_ ),                             &
              air_slant_column( air%ncells_ + 1 ) )
    allocate( o2_vertical_column( o2%ncells_ ),                               &
              o2_slant_column( o2%ncells_ + 1) )
    allocate( tuvx_o2_optical_depth( heights%ncells_, wavelengths%ncells_ ) )
    lut_heights(:)               = huge(1.0)
    lut_temperature(:)           = huge(1.0)
    lut_air_vertical_column(:)   = huge(1.0)
    lut_air_slant_column(:)      = huge(1.0)
    lut_o2_column(:)             = huge(1.0)
    lut_wavelength_edges(:)      = huge(1.0)
    lut_wavelength_centers(:)    = huge(1.0)
    lut_o2_base_cross_section(:) = huge(1.0)
    lut_o2_cross_section(:,:)    = huge(1.0)
    lut_o2_optical_depth(:,:)    = huge(1.0)

    do i_sza = 1, size( solar_zenith_angles )

      ! calculate slant O2 column
      call geometry%set_parameters( solar_zenith_angles( i_sza ), grids )
      call geometry%air_mass( air%exo_layer_dens_,                            &
                              air_vertical_column,                            &
                              air_slant_column )
      call geometry%air_mass( o2%exo_layer_dens_,                             &
                              o2_vertical_column,                             &
                              o2_slant_column )

      tuvx_o2_optical_depth(:,:) = 0.0_dk
      lut_o2_cross_section(:,:)  = 0.0
      lut_o2_optical_depth(:,:)  = 0.0

      ! get TUV-x O2 optical depths and cross sections
      call la_sr_bands%optical_depth( grids, profiles, air_vertical_column,   &
          air_slant_column, tuvx_o2_optical_depth, geometry )
      tuvx_o2_cross_section = o2_cross_section%calculate( grids, profiles )
      call la_sr_bands%cross_section( grids, profiles, air_vertical_column,   &
          air_slant_column, tuvx_o2_cross_section, geometry )

      ! get LUT O2 optical depths and cross sections
      n_heights                              = heights%ncells_ + 1
      n_wavelengths                          = lut_wavelengths%ncells_ + 1
      lut_heights(1:n_heights)               = real( heights%edge_(:) )
      lut_temperature(1:n_heights)           = real( temperature%edge_val_(:) )
      lut_wavelength_edges(1:n_wavelengths)  = real( lut_wavelengths%edge_(:) )
      lut_wavelength_centers(1:n_wavelengths-1) =                             &
                                               real( lut_wavelengths%mid_(:) )
      lut_o2_column(1:n_heights)             = real( o2_slant_column(:) )
      lut_air_vertical_column(1:air%ncells_) = real( air_vertical_column(:) )
      lut_air_slant_column(1:air%ncells_+1)  = real( air_slant_column(:) )

      file_path = "test/unit/tuv_doug/INPUT/XSQY/"
      call rdo2xs( n_wavelengths, lut_wavelength_edges,                       &
                   lut_wavelength_centers, lut_o2_base_cross_section,         &
                   file_path )

      call la_srb( n_heights, lut_heights, lut_temperature,                   &
                   n_wavelengths, lut_wavelength_edges, lut_o2_column,        &
                   lut_air_vertical_column, lut_air_slant_column,             &
                   lut_o2_base_cross_section, lut_o2_optical_depth,           &
                   lut_o2_cross_section, file_path )

      do i_height = 1, n_heights - 1
        do i_wl = 1, n_wavelengths - 1
          rel_tol = 1.0e-4
          if ( i_wl == 1 .or. i_wl == 3 ) cycle
          if ( i_wl == 20 .or. i_wl == 38 ) rel_tol = 0.5_dk
          if ( i_wl == 2 .and. i_height >= 112 ) rel_tol = 0.05_dk
          call assert( 624510149,                                             &
                       almost_equal( tuvx_o2_cross_section( i_height, i_wl ), &
                      real( lut_o2_cross_section( i_height, i_wl ), kind=dk ),&
                                     relative_tolerance = rel_tol ) )
          call assert( 746904813,                                             &
                       almost_equal( tuvx_o2_optical_depth( i_height, i_wl ), &
                      real( lut_o2_optical_depth( i_height, i_wl ), kind=dk ),&
                                     relative_tolerance = rel_tol ) )
        end do
      end do
      deallocate( tuvx_o2_cross_section )
    end do

    deallocate( heights )
    deallocate( wavelengths )
    deallocate( lut_wavelengths )
    deallocate( air )
    deallocate( o2 )
    deallocate( temperature )
    deallocate( geometry )
    deallocate( tuvx_o2_optical_depth )
    deallocate( air_vertical_column )
    deallocate( air_slant_column )
    deallocate( o2_vertical_column )
    deallocate( o2_slant_column )

  end subroutine compare_o2_cross_sections

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_la_srb