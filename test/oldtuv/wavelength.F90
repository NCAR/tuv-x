! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
module micm_photolysis_wavelength_grid

  use musica_constants, only : musica_dk, musica_ik

  implicit none

  type photolysis_wavelength_grid_t
    !> number of wavelength grid cells
    integer(musica_ik) :: nwave
    !> extra terrestial flux (photons/cm^2)
    real(musica_dk), allocatable :: etf(:)
    !> wavelength grid cell centers (nm)
    real(musica_dk), allocatable :: wcenter(:)
    !> wavelength grid cell edges (nm)
    real(musica_dk), allocatable :: wedge(:)
  end type photolysis_wavelength_grid_t

  type(photolysis_wavelength_grid_t) :: wavelength_grid

contains

  function wavelength_grid_initialize(filespec) result(retcode)

    use musica_io,                     only : io_t
    use musica_io_netcdf,              only : io_netcdf_t
    use musica_string,                 only : string_t

    integer(musica_ik)              :: retcode
    character(len=*), intent(in)    :: filespec

    integer(musica_ik), parameter :: noErr = 0_musica_ik
    character(len=*), parameter   :: Iam = 'wavelength_grid_initialize: '
    integer(musica_ik), allocatable :: dims(:)
    class(io_t), pointer :: nc_file
    type(string_t) :: file_path

    file_path = filespec
    nc_file => io_netcdf_t( file_path )

    var_name = 'etf'
    call nc_file%read( var_name, wavelength_grid%etf, Iam )
    var_name = 'wc'
    call nc_file%read( var_name, wavelength_grid%wcenter, Iam )
    var_name = 'wl'
    call nc_file%read( var_name, wavelength_grid%wedge, Iam )
    wavelength_grid%nwave = size( wavelength_grid%etf, 1 )

    deallocate( nc_file )

  end function wavelength_grid_initialize

end module micm_photolysis_wavelength_grid
