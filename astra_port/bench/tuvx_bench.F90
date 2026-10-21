! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! The column sweep of `turbotuv_bench --sweep`, run through the original
! TUV-x core: one model, column c at solar zenith angle 89 c / N degrees with
! unit radiation scale. Only core%run is timed, after one untimed column. The
! FNV-1a hash covers the six radiation fields and the photolysis, heating and
! dose rates, in the order and memory layout `turbotuv_bench --sweep` hashes,
! so equal hashes mean bit-identical results. Run it where TUV-x's relative
! data paths resolve, such as the reference test runtime directory.
!
!   tuvx_bench CONFIG [--columns N] [--repeat R]
program tuvx_bench
  use iso_fortran_env, only : int8, int64, real64
  use musica_config,   only : config_t
  use tuvx_core,       only : core_t
  use tuvx_grid,       only : grid_t
  use tuvx_solver,     only : radiation_field_t
  implicit none
  type(config_t) :: config
  class(core_t), pointer :: model
  class(grid_t), pointer :: grid
  type(radiation_field_t) :: field
  character(len=1024) :: path, arg
  character(len=16) :: text
  integer :: columns, repeat, nz, nw, np, nh, nd, i, c, r
  integer(int64) :: start, finish, rate, ticks
  ! FNV-1a 64 state as four 16-bit limbs, least significant first, so that
  ! the multiplication modulo 2**64 never overflows a signed integer.
  integer(int64) :: hash(0:3)
  real(real64), allocatable :: photo(:,:), heat(:,:), dose(:,:)
  real(real64) :: zenith, elapsed

  path = ''
  columns = 1000
  repeat = 1
  i = 1
  do while( i <= command_argument_count( ) )
    call get_command_argument( i, arg )
    if( arg == '--columns' .or. arg == '--repeat' ) then
      call get_command_argument( i + 1, text )
      if( arg == '--columns' ) then
        read( text, * ) columns
      else
        read( text, * ) repeat
      end if
      i = i + 2
    else
      path = arg
      i = i + 1
    end if
  end do
  if( path == '' .or. columns < 1 .or. repeat < 1 ) then
    write( *, '(a)' ) 'Usage: tuvx_bench CONFIG [--columns N] [--repeat R]'
    stop 2
  end if

  call config%from_file( trim( path ) )
  model => core_t( config )
  grid => model%get_grid( 'height', 'km' )
  nz = grid%size( ) + 1
  deallocate( grid )
  grid => model%get_grid( 'wavelength', 'nm' )
  nw = grid%size( )
  deallocate( grid )
  np = model%number_of_photolysis_reactions( )
  nh = model%number_of_heating_rates( )
  nd = model%number_of_dose_rates( )
  allocate( photo( nz, np ), heat( nz, nh ), dose( nz, nd ) )

  hash = [ int( z'2325', int64 ), int( z'8422', int64 ), int( z'9CE4', int64 ),     &
           int( z'CBF2', int64 ) ]
  ticks = 0
  call model%run( 0.0_real64, 1.0_real64, photo, dose, heat )
  do r = 1, repeat
    do c = 0, columns - 1
      zenith = 89.0_real64 * real( c, real64 ) / real( columns, real64 )
      call system_clock( start, rate )
      call model%run( zenith, 1.0_real64, photo, dose, heat )
      call system_clock( finish )
      ticks = ticks + ( finish - start )
      if( r > 1 ) cycle
      field = model%get_radiation_field( )
      call digest( reshape( field%fdr_, [ size( field%fdr_ ) ] ) )
      call digest( reshape( field%fdn_, [ size( field%fdn_ ) ] ) )
      call digest( reshape( field%fup_, [ size( field%fup_ ) ] ) )
      call digest( reshape( field%edr_, [ size( field%edr_ ) ] ) )
      call digest( reshape( field%edn_, [ size( field%edn_ ) ] ) )
      call digest( reshape( field%eup_, [ size( field%eup_ ) ] ) )
      call digest( reshape( photo, [ size( photo ) ] ) )
      call digest( reshape( heat, [ size( heat ) ] ) )
      call digest( reshape( dose, [ size( dose ) ] ) )
    end do
  end do
  elapsed = real( ticks, real64 ) / real( rate, real64 )

  write( text, '(4z4.4)' ) hash(3), hash(2), hash(1), hash(0)
  do i = 1, len( text )
    if( text(i:i) >= 'A' .and. text(i:i) <= 'F' )                               &
      text(i:i) = achar( iachar( text(i:i) ) + 32 )
  end do
  write( *, '(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0,a)' ) 'config=' // trim( path ) //   &
      ' columns=', columns * repeat, ' levels=', nz, ' wavelengths=', nw,       &
      ' photolysis=', np, ' heating=', nh, ' dose=', nd, ' execution=TUV-x'
  write( *, '(a,f0.3,a,f0.2,a)' ) 'elapsed_ms=', 1.0e3_real64 * elapsed,        &
      ' per_column_us=', 1.0e6_real64 * elapsed / real( columns * repeat, real64 ), &
      ' hash=' // text

contains

  ! Add the bytes of values, in memory order, to the FNV-1a hash.
  subroutine digest( values )
    real(real64), intent(in) :: values(:)
    integer(int8), allocatable :: bytes(:)
    integer :: b
    bytes = transfer( values, [ 0_int8 ] )
    do b = 1, size( bytes )
      hash(0) = ieor( hash(0), iand( int( bytes(b), int64 ), 255_int64 ) )
      call multiply( )
    end do
  end subroutine digest

  ! hash = hash * 1099511628211 modulo 2**64, the prime in 16-bit limbs.
  subroutine multiply( )
    integer(int64), parameter :: prime(0:3) = [ int( z'01B3', int64 ), 0_int64,    &
                                                int( z'0100', int64 ), 0_int64 ]
    integer(int64) :: product(0:3), carry
    integer :: k, j
    carry = 0
    do k = 0, 3
      product(k) = carry
      do j = 0, k
        product(k) = product(k) + hash(j) * prime(k - j)
      end do
      carry = ishft( product(k), -16 )
      product(k) = iand( product(k), 65535_int64 )
    end do
    hash = product
  end subroutine multiply

end program tuvx_bench
