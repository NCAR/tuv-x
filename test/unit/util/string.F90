! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_string module

!> Test module for the musica_string module
program test_util_string

  use musica_constants,                only : musica_ik, musica_rk, musica_dk
  use musica_assert
  use musica_mpi
  use musica_string
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

  character(len=256) :: failure_test_type

  call musica_mpi_init( )

  if( command_argument_count( ) .eq. 0 ) then
    call test_string_t( )
    if( musica_mpi_rank( MPI_COMM_WORLD ) == 0 ) then
      call replace_example( )
      call substring_example( )
      call split_example( )
      call table_test( )
    end if
  else if( command_argument_count( ) .eq. 1 ) then
    call get_command_argument( 1, failure_test_type )
    call failure_test( failure_test_type )
  else
    call die( 253391339 )
  end if

  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test string_t functionality
  subroutine test_string_t( )

    type(string_t) :: a, b, c, unalloced
    type(string_t), allocatable :: split_string(:)
    integer(kind=musica_ik) :: i
    real(kind=musica_rk) :: r
    logical :: l
    real(kind=musica_dk) :: d
    character(len=10) :: ca
    character(len=:), allocatable :: aca
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD

    if( musica_mpi_rank( comm ) == 0 ) then
      a = "an MPI test string"
      pack_size = a%pack_size( comm ) + b%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call a%mpi_pack( buffer, pos, comm )
      call b%mpi_pack( buffer, pos, comm )
    end if

    call musica_mpi_bcast( pack_size, comm )

    if( musica_mpi_rank( comm ) /= 0 ) allocate( buffer( pack_size ) )

    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) /= 0 ) then
      pos = 0
      call a%mpi_unpack( buffer, pos, comm )
      call b%mpi_unpack( buffer, pos, comm )
      call assert( 165108300, pos <= pack_size )
      call assert( 616287214, allocated( a%val_ ) )
      call assert( 619645987, a .eq. "an MPI test string" )
      call assert( 226757927, .not. allocated( b%val_ ) )
    end if

#ifdef MUSICA_USE_OPENMP
    write(*,*) "Testing string module on ", omp_get_max_threads( ), " threads"
#else
    write(*,*) "Testing string module without OpenMP support"
#endif

    !$omp parallel &
    !$omp   private( a, b, c, unalloced, split_string, i, r, l, d, ca, aca )

    ! string assignment

    a = "test string  "
    call assert( 814138261, a .eq. "test string" )

    b = a
    call assert( 240083225, a .eq. b )

    deallocate( a%val_ )
    b = a
    call assert( 325157354, .not. allocated( b%val_ ) )

    a = 1469
    call assert( 124406107, a .eq. "1469" )

    a = 13.4563
    call assert( 915898829, a%substring(1,6) .eq. "13.456" )

    a = 14.563d0
    call assert( 381828325, a%substring(1,5) .eq. "14.56" )

    a = .true.
    call assert( 827742932, a .eq. "true" )

    a = .false.
    call assert( 317272220, a .eq. "false" )

    ! string join to's

    a = "foo"
    b = "bar"
    c = a//b
    call assert( 938608038, c .eq. "foobar" )

    c = b//"foo"
    call assert( 817666613, c .eq. "barfoo" )

    c = a//123
    call assert( 984464744, c .eq. "foo123" )

    c = a//52.33
    call assert( 810949067, c%substring(1,7) .eq. "foo52.3" )

    c = a//53.43d0
    call assert( 419966541, c%substring(1,7) .eq. "foo53.4" )

    c = a//.true.
    call assert( 581500365, c .eq. "footrue" )

    c = a//.false.
    call assert( 971029652, c .eq. "foofalse" )

    c = "bar "//a
    call assert( 670923431, c .eq. "bar foo" )

    ! equality

    a = "foo"
    b = "foo"
    c = "bar"
    call assert( 160576005, a .eq. b )
    call assert( 667687944, .not. a .eq. c )

    call assert( 322109829, a .eq. "foo" )
    call assert( 264271270, .not. a .eq. "bar" )

    a = 134
    call assert( 325920897, a .eq. 134 )
    call assert( 315392283, .not. a .eq. 432 )

    a = 52.3
    call assert( 420993082, a .eq. 52.3 )
    call assert( 428162923, .not. a .eq. 762.4 )

    a = 87.45d0
    call assert( 307221498, a .eq. 87.45d0 )
    call assert( 696750785, .not. a .eq. 43.5d9 )

    a = .true.
    b = .false.
    call assert( 240759859, a .eq. .true. )
    call assert( 919934236, .not. a .eq. .false. )
    call assert( 179562527, b .eq. .false. )
    call assert( 969149715, .not. b .eq. .true. )

    ! not-equals

    a = "foo"
    b = "foo"
    c = "bar"
    call assert( 678503681, .not. a .ne. b )
    call assert( 173297276, a .ne. c )

    call assert( 903140371, .not. a .ne. "foo" )
    call assert( 732983467, a .ne. "bar" )

    a = 134
    call assert( 845301812, .not. a .ne. 134 )
    call assert( 957620157, a .ne. 432 )

    a = 52.3
    call assert( 787463253, .not. a .ne. 52.3 )
    call assert( 334831100, a .ne. 762.4 )

    a = 87.45d0
    call assert( 447149445, .not. a .ne. 87.45d0 )
    call assert( 894517291, a .ne. 43.5d9 )

    a = .true.
    b = .false.
    call assert( 389310886, .not. a .ne. .true. )
    call assert( 501629231, a .ne. .false. )
    call assert( 948997077, .not. b .ne. .false. )
    call assert( 778840173, b .ne. .true. )

    ! case convert

    a = "FoObAr 12 %"
    call assert( 500463115, a%to_lower( ) .eq. "foobar 12 %" )
    call assert( 614686994, a%to_upper( ) .eq. "FOOBAR 12 %" )

    ! substring

    call assert( 328852972, a%substring(1,6) .eq. "FoObAr" )
    call assert( 272919947, a%substring(4,5) .eq. "bAr 1" )
    call assert( 604610675, a%substring(7,20) .eq. " 12 %" )

    ! split

    a = "foobar1foofoobar2foofoo"
    b = "foo"
    split_string = a%split( b )
    call assert( 106051866, size( split_string ) .eq. 6 )
    call assert( 815260865, split_string(1) .eq. ""     )
    call assert( 432478287, split_string(2) .eq. "bar1" )
    call assert( 805184546, split_string(3) .eq. ""     )
    call assert( 381809569, split_string(4) .eq. "bar2" )
    call assert( 417108498, split_string(5) .eq. ""     )
    call assert( 742081680, split_string(6) .eq. ""     )

    split_string = a%split( b, compress = .true. )
    call assert( 413749725, size( split_string ) .eq. 2 )
    call assert( 238328514, split_string(1) .eq. "bar1" )
    call assert( 456247658, split_string(2) .eq. "bar2" )

    split_string = a%split( "bar" )
    call assert( 883657201, size( split_string ) .eq. 3 )
    call assert( 655661738, split_string(1) .eq. "foo" )
    call assert( 480240527, split_string(2) .eq. "1foofoo" )

    split_string = a%split( "bar", compress = .true. )
    call assert( 983657201, size( split_string ) .eq. 3 )
    call assert( 455661738, split_string(1) .eq. "foo" )
    call assert( 680240527, split_string(2) .eq. "1foofoo" )
    call assert( 104877217, split_string(3) .eq. "2foofoo" )

    split_string = a%split( "not in there" )
    call assert( 366468943, size( split_string ) .eq. 1 )
    call assert( 473522981, split_string(1) .eq. a )

    split_string = a%split( "" )
    call assert( 357845863, size( split_string ) .eq. 1 )
    call assert( 300007304, split_string(1) .eq. a )

    a = "foo bar"
    split_string = a%split( " " )
    call assert( 484519904, size( split_string ) .eq. 2 )
    call assert( 853182732, split_string(1) .eq. "foo" )
    call assert( 737505614, split_string(2) .eq. "bar" )

    deallocate( a%val_ )
    split_string = a%split( " " )
    call assert( 341895943, allocated( split_string ) )
    call assert( 966590534, size( split_string ) .eq. 0 )


    ! replace

    a = "foobar1foobar2foo1"
    b = a%replace( "foo", "bar" )
    call assert( 282451682, b .eq. "barbar1barbar2bar1" )
    b = a%replace( "bar", "foo" )
    call assert( 331667161, b .eq. "foofoo1foofoo2foo1" )

    ! convert to character array
    a = "string to convert"
    aca = a%to_char( )
    call assert( 476488677, aca .eq. "string to convert" )

    ! assignment from string

    ca = "XXXXXXXXXX"
    a = "foo"
    ca = a
    call assert( 189690040, trim( ca ) .eq. "foo" )

    ca = "XXXXXXXXXX"
    deallocate( a%val_ )
    ca = a
    call assert( 137411891, trim( ca ) .eq. "" )

    ca = "XXXXXXXXXX"
    a = "12345678901234567890"
    ca = a
    call assert( 811321961, trim( ca ) .eq. "1234567890" )

    a = "-12.02"
    r = a
    call assert( 179687753,                                                   &
                 almost_equal( real( r, kind=musica_dk ),                     &
                               real( -12.02, kind=musica_dk ) ) )

    a = "32.54"
    d = a
    call assert( 321521234, almost_equal( d, 32.54d0 ) )

    a = "-14"
    i = a
    call assert( 464068536, i .eq. -14 )

    a = "true"
    l = a
    call assert( 853597823, l )

    a = "false"
    l = a
    call assert( 237978607, .not. l )

    ! joins from strings

    ca = "foo"
    a = "bar"
    call assert( 511304449, trim( ca )//a .eq. "foobar" )

    i = 122
    call assert( 678841998, i//a .eq. "122bar" )

    r = 34.63
    b = r//a
    call assert( 165012513, b%substring(1,4) .eq. "34.6" )
    call assert( 610927120, b%substring( b%length( ) - 2, 3 ) .eq. "bar" )

    d = 43.63d0
    b = d//a
    call assert( 625841048, b%substring(1,4) .eq. "43.6" )
    call assert( 848572204, b%substring( b%length( ) - 2, 3 ) .eq. "bar" )

    call assert( 345271333, .true.//a .eq. "truebar" )
    call assert( 164585815, .false.//a .eq. "falsebar" )

    ! equality

    a = "foo"
    b = "foo"
    call assert( 719459994, "foo" .eq. a )
    call assert( 549303090, .not. "bar" .eq. b )

    a = 134
    call assert( 944096684, 134 .eq. a )
    call assert( 773939780, .not. 432 .eq. a )

    a = 52.3
    call assert( 603782876, 52.3 .eq. a )
    call assert( 433625972, .not. 762.4 .eq. a )

    a = 87.45d0
    call assert( 828419566, 87.45d0 .eq. a )
    call assert( 375787413, .not. 43.5d9 .eq. a )

    a = .true.
    b = .false.
    call assert( 153056257, .true. .eq. a )
    call assert( 882899352, .not. .false. .eq. a )
    call assert( 995217697, .false. .eq. b )
    call assert( 542585544, .not. .true. .eq. b )

    ! not-equals

    a = "foo"
    b = "foo"
    call assert( 597065330, .not. "foo" .ne. a )
    call assert( 426908426, "bar" .ne. a )

    a = 134
    call assert( 539226771, .not. 134 .ne. a )
    call assert( 369069867, 432 .ne. a )

    a = 52.3
    call assert( 146338711, .not. 52.3 .ne. a )
    call assert( 876181806, 762.4 .ne. a )

    a = 87.45d0
    call assert( 706024902, .not. 87.45d0 .ne. a )
    call assert( 535867998, 43.5d9 .ne. a )

    a = .true.
    b = .false.
    call assert( 648186343, .not. .true. .ne. a )
    call assert( 760504688, .false. .ne. a )
    call assert( 872823033, .not. .false. .ne. b )
    call assert( 702666129, .true. .ne. b )

    ca = "XXXXXXXXXX"
    ca = to_char( 345 )
    call assert( 278095873, trim( ca ) .eq. "345" )

    ca = "XXXXXXXXXX"
    ca = to_char( 482.53 )
    call assert( 876921224, ca(1:5) .eq. "482.5" )

    ca = "XXXXXXXXXX"
    ca = to_char( 873.453d0 )
    call assert( 989239569, ca(1:6) .eq. "873.45" )

    ca = "XXXXXXXXXX"
    ca = to_char( .true. )
    call assert( 201557915, trim( ca ) .eq. "true" )

    ca = "XXXXXXXXXX"
    ca = to_char( .false. )
    call assert( 931401010, trim( ca ) .eq. "false" )

    !$omp end parallel

  end subroutine test_string_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Replace example from documentation
  subroutine replace_example( )

type(string_t) :: my_string
my_string = "foo bar foobar"
my_string = my_string%replace( 'foo', 'bar' )
write(*,*) my_string%val_

  end subroutine replace_example

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Split example from documentation
  subroutine split_example( )

type(string_t) :: my_string
type(string_t), allocatable :: sub_strings(:)
integer :: i
my_string = "my original    string"
sub_strings = my_string%split( ' ' )
do i = 1, size( sub_strings )
  write(*,*) i, sub_strings( i )%val_
end do
sub_strings = my_string%split( ' ', .true. )
do i = 1, size( sub_strings )
  write(*,*) i, sub_strings( i )%val_
end do

  end subroutine split_example

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Substring example from documentation
  subroutine substring_example( )

type(string_t) :: my_string, sub_string
my_string = "Hi there!"
sub_string = my_string%substring( 4, 5 )
write(*,*) sub_string%val_
sub_string = my_string%substring( 9, 50 )
write(*,*) sub_string%val_

  end subroutine substring_example

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test table output
  subroutine table_test( )

    type(string_t) :: header(3)
    type(string_t) :: table(3,2)
    character(len=256) :: line

    header(1) = "foo"
    header(2) = "bar"
    header(3) = "baz"
    table(1,1) = "f1"
    table(1,2) = "f2"
    table(2,1) = "b1b1b1"
    table(2,2) = "b2"
    table(3,1) = "z1"
    table(3,2) = "z2"
    open( 12, file = "output_table.txt", status = "replace" )
    call output_table( header, table, 12 )
    close( 12 )
    open( 12, file = "output_table.txt", status = "old" )
    read( 12, '(A)' ) line
    call assert( 635926347, trim( line ) .eq. "----------------------" )
    read( 12, '(A)' ) line
    call assert( 804630012, trim( line ) .eq. "| foo | bar    | baz |" )
    read( 12, '(A)' ) line
    call assert( 351997859, trim( line ) .eq. "----------------------" )
    read( 12, '(A)' ) line
    call assert( 799365705, trim( line ) .eq. "| f1  | b1b1b1 | z1  |" )
    read( 12, '(A)' ) line
    call assert( 911684050, trim( line ) .eq. "| f2  | b2     | z2  |" )
    read( 12, '(A)' ) line
    call assert( 124002396, trim( line ) .eq. "----------------------" )
    close( 12 )

    table(2,1) = "1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9"//&
                " 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8"//&
                " 9 0 1 2 3 4 5 6 7 8 9 0"
    open( 12, file = "output_table_2.txt", status = "replace" )
    call output_table( header, table, 12 )
    close( 12 )
    open( 12, file = "output_table_2.txt", status = "old" )
    read( 12, '(A)' ) line
    call assert( 284539358, trim( line ) .eq. "----------------------------"//&
      "--------------------------------------------------------------------"//&
      "--------------------------------" )
    read( 12, '(A)' ) line
    call assert( 114382454, trim( line ) .eq. "| fo | bar                  "//&
      "                                                                    "//&
      "                           | ba |" )
    read( 12, '(A)' ) line
    call assert( 844225549, trim( line ) .eq. "----------------------------"//&
      "--------------------------------------------------------------------"//&
      "--------------------------------" )
    read( 12, '(A)' ) line
    call assert( 674068645, trim( line ) .eq. "| f1 | 1 2 3 4 5 6 7 8 9 0 1"//&
      " 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5"//&
      " 6 7 8 9 0 1 2 3 4 5 6 7 8 | z1 |" )
    read( 12, '(A)' ) line
    call assert( 503911741, trim( line ) .eq. "| f2 | b2                   "//&
      "                                                                    "//&
      "                           | z2 |" )
    read( 12, '(A)' ) line
    call assert( 333754837, trim( line ) .eq. "----------------------------"//&
      "--------------------------------------------------------------------"//&
      "--------------------------------" )
    close( 12 )

  end subroutine table_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Failure tests for string_t class
  subroutine failure_test( test_type )

    character(len=*), intent(in) :: test_type

    if( test_type .eq. "359920976" ) then
      call failure_test_359920976( )
    else
      call die( 592539031 )
    end if

  end subroutine failure_test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test invalid logical assignment from string
  subroutine failure_test_359920976( )

    type(string_t) :: string
    logical :: bar

    string = "foo"
    bar = string

  end subroutine failure_test_359920976

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_string
