! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_netcdf
  ! NetCDF I/O class

  use musica_constants,                only : musica_dk

  implicit none

  private
  public :: netcdf_t, clean_string

  type netcdf_t
    ! NetCDF I/O
    real(musica_dk), allocatable :: wavelength(:)
    real(musica_dk), allocatable :: temperature(:)
    real(musica_dk), allocatable :: parameters(:,:)
  contains
    procedure :: read_netcdf_file => run
    final     :: finalize
  end type netcdf_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine run( this, file_path, variable_name )
    ! Reads data from a NetCDF file

    use musica_assert,                    only : assert_msg
    use musica_io,                        only : io_t
    use musica_io_netcdf,                 only : io_netcdf_t
    use musica_string,                    only : string_t

    class(netcdf_t), intent(inout) :: this
    character(len=*), intent(in)   :: file_path
    character(len=*), intent(in)   :: variable_name

    integer, parameter :: noErr = 0
    character(len=*), parameter :: Iam = "read_netcdf_file: "

    class(io_t), pointer :: nc_file
    integer :: nLambda
    type(string_t) :: var_name
    type(string_t) :: file_path_str

    !  open the netcdf file
    file_path_str = file_path
    nc_file => io_netcdf_t( file_path_str, read_only = .true. )

    !  parameter array must be in netcdf file, read it
    var_name = trim( variable_name ) // 'parameters'
    call assert_msg( 118377216, nc_file%exists( var_name, Iam ),              &
                     "NetCDF file '"//file_path//"' must have a "//           &
                     "'parameters' variable." )
    call nc_file%read( var_name, this%parameters, Iam )
    nLambda = size( this%parameters, 1 )

    ! if it exists, read wavelength array
    var_name = 'wavelength'
    if( nc_file%exists( var_name, Iam ) ) then
      call nc_file%read( var_name, this%wavelength, Iam )
      call assert_msg( 944197086, size( this%wavelength, 1 ) == nLambda,      &
                       "Wavelength and parameters array size mismatch in '"// &
                       file_path//"'" )
    endif

    ! if it exists, read temperature array
    var_name = 'temperature'
    if( nc_file%exists( var_name, Iam ) ) then
      call nc_file%read( var_name, this%temperature, Iam )
    endif

    !  close the netcdf file
    deallocate( nc_file )

  end subroutine run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleans up memory of the NetCDF class

    type(netcdf_t), intent(inout) :: this

    if( allocated( this%wavelength ) ) then
      deallocate( this%wavelength )
    endif
    if( allocated( this%temperature ) ) then
      deallocate( this%temperature )
    endif
    if( allocated( this%parameters ) ) then
      deallocate( this%parameters )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function clean_string( from )

    use musica_string,                 only : string_t

    type(string_t), intent(in) :: from

    clean_string = from%replace( "/", "_" )

  end function clean_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_netcdf
