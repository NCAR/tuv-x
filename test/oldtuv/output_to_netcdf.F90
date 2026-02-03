! Copyright (C) 2020-2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
! Outputs dose and photolysis rates to NetCDF files
module output_to_netcdf

  implicit none

  private
  public :: output_dose_rates, output_photolysis_rates

contains

  function clean_string( str )
    ! Cleans a string by removing leading and trailing spaces
    character(len=*), intent(in) :: str
    character(len=len(str)) :: clean_string
    integer :: i
    clean_string = trim(adjustl(str))
    do i = 1, len(clean_string)
      select case (clean_string(i:i))
        case ('/')
          clean_string(i:i) = '_'
      end select
    end do
  end function clean_string

  subroutine output_dose_rates( file_name, dose_rate_labels, dose_rates, &
      number_of_szas, number_of_rates, number_of_levels )

    ! Outputs dose rates to a NetCDF file
    use musica_assert, only : assert_msg
    use musica_io_netcdf, only : io_netcdf_t
    use musica_string, only : string_t
    use musica_constants, only : musica_dk
    use musica_io, only : io_t
    
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: dose_rate_labels(:)
    real(musica_dk), intent(in) :: dose_rates(:,:,:) ! 3D array: sza, rate, level
    integer, intent(in) :: number_of_szas
    integer, intent(in) :: number_of_rates
    integer, intent(in) :: number_of_levels
    integer :: i_rate
    class(io_t), pointer :: nc_file
    type(string_t) :: file_name_str, var_name, n_level_name, n_sza_name
    
    ! Open the NetCDF file for writing
    file_name_str = file_name
    nc_file => io_netcdf_t(file_name_str, .false.)
    call assert_msg( 118377217, associated(nc_file), "Failed to open file: "//file_name )
    n_sza_name = 'number_of_solar_zenith_angles'
    n_level_name = 'number_of_levels'
    do i_rate = 1, number_of_rates
      var_name = clean_string(dose_rate_labels(i_rate))
      call nc_file%write( var_name, (/ n_sza_name, n_level_name /), &
                          dose_rates(1:number_of_szas, i_rate, 1:number_of_levels), &
                          'old TUV dose rate output' )
    end do
    deallocate( nc_file )
    write(*,*) "Dose rates successfully written to: ", file_name
  end subroutine output_dose_rates

  subroutine output_photolysis_rates( file_name, photolysis_rate_labels, &
      photolysis_rates, number_of_szas, number_of_rates, number_of_levels )

    ! Outputs photolysis rates to a NetCDF file
    use musica_assert, only : assert_msg
    use musica_io_netcdf, only : io_netcdf_t
    use musica_string, only : string_t
    use musica_constants, only : musica_dk
    use musica_io, only : io_t
    
    character(len=*), intent(in) :: file_name
    character(len=*), intent(in) :: photolysis_rate_labels(:)
    real(musica_dk), intent(in) :: photolysis_rates(:,:,:) ! 3D array: sza, rate, level
    integer, intent(in) :: number_of_szas
    integer, intent(in) :: number_of_rates
    integer, intent(in) :: number_of_levels
    integer :: i_rate
    class(io_t), pointer :: nc_file
    type(string_t) :: file_name_str, var_name, n_level_name, n_sza_name
    
    ! Open the NetCDF file for writing
    file_name_str = file_name
    nc_file => io_netcdf_t(file_name_str, .false.)
    call assert_msg( 118377219, associated(nc_file), "Failed to open file: "//file_name )
    n_sza_name = 'number_of_solar_zenith_angles'
    n_level_name = 'number_of_levels'
    do i_rate = 1, number_of_rates
      var_name = clean_string(photolysis_rate_labels(i_rate))
      call nc_file%write( var_name, (/ n_sza_name, n_level_name /), &
                          photolysis_rates(1:number_of_szas, i_rate, 1:number_of_levels), &
                          'old TUV photolysis rate output' )
    end do
    deallocate( nc_file )
    write(*,*) "Photolysis rates successfully written to: ", file_name
  end subroutine output_photolysis_rates
end module output_to_netcdf