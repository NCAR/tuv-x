! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
module netcdf_util
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

  subroutine run( this, filespec, Hdr )
    ! Reads data from a NetCDF file

    use musica_assert,                    only : assert_msg
    use musica_io,                        only : io_t
    use musica_io_netcdf,                 only : io_netcdf_t
    use musica_string,                    only : string_t

    class(netcdf_t), intent(inout) :: this
    character(len=*), intent(in)   :: filespec
    character(len=*), intent(in)   :: Hdr

    integer, parameter :: noErr = 0
    character(len=*), parameter :: Iam = "read_netcdf_file: "

    class(io_t), pointer :: nc_file
    integer :: stat, nLambda
    integer, allocatable :: dims(:)
    type(string_t) :: var_name
    type(string_t) :: filespec_str

    !  open the netcdf file
    filespec_str = filespec
    nc_file => io_netcdf_t( filespec_str )

    !  parameter array must be in netcdf file, read it
    var_name = trim( Hdr ) // 'parameters'
    call assert_msg( 118377216, nc_file%exists( var_name, Iam ),              &
                     "File '"//filespec//"' must have 'parameters' field." )
    call nc_file%read( var_name, this%parameters, Iam )
    nLambda = size( this%parameters, 1 )

    ! if it exists, read wavelength array
    var_name = 'wavelength'
    if( nc_file%exists( var_name, Iam ) ) then
      call nc_file%read( var_name, this%wavelength, Iam )
      call assert_msg( 944197086, size( this%wavelength, 1 ) == nLambda,      &
                       "Wavelength and parameters array size mismatch in '"// &
                       filespec//"'" )
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

end module netcdf_util
