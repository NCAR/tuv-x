! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_config module

!> Test module for the musica_config module
program test_config

  use musica_assert
  use musica_config
  use musica_mpi

  implicit none

  call musica_mpi_init( )
  call test_config_t_mpi( )
  call test_config_t( )
  call config_example( )
  call musica_mpi_finalize( )

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test config_t MPI functions
  subroutine test_config_t_mpi( )

    use musica_string,                 only : string_t

    type(config_t) :: a, b
    type(string_t) :: sa
    character, allocatable :: buffer(:)
    integer :: pos, pack_size
    integer, parameter :: comm = MPI_COMM_WORLD
    character(len=*), parameter :: my_name = "config tests"

    if( musica_mpi_rank( comm ) == 0 ) then
      a = '{ "foo": "bar" }'
      pack_size = a%pack_size( comm )
      allocate( buffer( pack_size ) )
      pos = 0
      call a%mpi_pack( buffer, pos, comm )
    end if

    call musica_mpi_bcast( pack_size, comm )

    if( musica_mpi_rank( comm ) /= 0 ) allocate( buffer( pack_size ) )

    call musica_mpi_bcast( buffer, comm )

    if( musica_mpi_rank( comm ) /= 0 ) then
      pos = 0
      call b%mpi_unpack( buffer, pos, comm )
      call b%get( "foo", sa, my_name )
      call assert( 529948470, sa == "bar" )
    end if

  end subroutine test_config_t_mpi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test config_t functionality
  subroutine test_config_t( )

    use musica_constants,              only : musica_rk, musica_dk, musica_ik
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t

    type(config_t) :: a, a_file, b, c, array
    type(config_t), allocatable :: orig_array(:), dest_array(:)
    real(kind=musica_rk) :: ra
    real(kind=musica_dk) :: da
    real(kind=musica_dk), allocatable :: daa(:), dab(:)
    integer(kind=musica_ik) :: ia
    logical :: la, found
    type(string_t) :: sa, sb
    type(string_t), allocatable :: saa(:), sab(:)
    character(len=*), parameter :: my_name = "config tests"
    class(iterator_t), pointer :: iterator

    ! constructors
    a = '{ "foo": "bar" }'
    call a%empty( )
    call a_file%from_file( "test/data/test_config.json" )
    if( musica_mpi_rank( MPI_COMM_WORLD ) .eq. 0 ) then
      call a_file%to_file( "temp_file.json" )
      call a_file%empty( )
      call a_file%from_file( "temp_file.json" )
    end if

    ! size
    a = '{ "foo": "bar", "baz": "qux" }'
    call assert( 917322918, a%number_of_children() .eq. 2 )

    ! get config
    call a_file%get( "my sub object", b, my_name, found = found )
    call assert( 169832207, found )

    call b%get( "sub real", da, my_name )
    call assert( 630635145, almost_equal( da, 87.3d0 ) )

    call b%get( "sub int", ia, my_name )
    call assert( 892957756, ia .eq. 42 )

    call b%get( "really?", la, my_name )
    call assert( 389656885, la )

    call b%get( "a bunch of strings", saa, my_name )
    call assert( 603764961, size( saa ) .eq. 3 )
    call assert( 210876901, saa(1) .eq. "bar" )
    call assert( 325100780, saa(2) .eq. "foo" )
    call assert( 202253821, saa(3) .eq. "barfoo" )

    call a_file%get( "not there", b, my_name, found = found )
    call assert( 430701579, .not. found )

    c = '{ "an int" : 13, "foo" : "bar" }'
    call a_file%get( "not there", b, my_name, default = c, found = found )
    call assert( 250468356, .not. found )
    call b%get( "foo", sa, my_name )
    call assert( 464576432, sa .eq. "bar" )
    call b%get( "an int", ia, my_name )
    call assert( 457145065, ia .eq. 13 )

    ! get string

    call a_file%get( "a string", sa, my_name )
    call assert( 651552798, sa .eq. "foo" )
    call a_file%get( "another string", sa, my_name, found = found )
    call assert( 411575482, found )
    call assert( 927310501, sa .eq. "bar" )
    call a_file%get( "a string", sa, my_name, default = "default value" )
    call assert( 292539591, sa .eq. "foo" )
    call a_file%get( "not there", sa, my_name, default = "default value", found = found )
    call assert( 968355195, .not. found )
    call assert( 345566138, sa .eq. "default value" )
    call a_file%get( "also not there", sa, my_name, found = found )
    call assert( 564491555, .not. found )

    ! get integer

    call a_file%get( "another int", ia, my_name )
    call assert( 851875875, ia .eq. 31 )
    call a_file%get( "my integer", ia, my_name, found = found )
    call assert( 338046390, found )
    call assert( 397790483, ia .eq. 12 )
    call a_file%get( "another int", ia, my_name, default = 42 )
    call assert( 271584751, ia .eq. 31 )
    call a_file%get( "not there", ia, my_name, default = 96, found = found )
    call assert( 440288416, .not. found )
    call assert( 382449857, ia .eq. 96 )
    call a_file%get( "also not there", ia, my_name, found = found )
    call assert( 395787890, .not. found )

    ! get real

    call a_file%get( "this real", ra, my_name )
    call assert( 821646918, almost_equal( ra, 23.4 ) )
    call a_file%get( "that real", ra, my_name, found = found )
    call assert( 425400085, found )
    call assert( 702611027, almost_equal( ra, 52.3e-4 ) )
    call a_file%get( "this real", ra, my_name, default = 432.5 )
    call assert( 901830772, almost_equal( ra, 23.4e0 ) )
    call a_file%get( "not there", ra, my_name, default = 643.78, found = found )
    call assert( 505583939, .not. found )
    call assert( 165270131, ra .eq. 643.78 )
    call a_file%get( "also not there", ra, my_name, found = found )
    call assert( 736101698, .not. found )

    ! get double

    call a_file%get( "this real", da, my_name )
    call assert( 155933230, almost_equal( da, 23.4d0 ) )
    call a_file%get( "that real", da, my_name, found = found )
    call assert( 550726824, found )
    call assert( 663045169, almost_equal( da, 52.3d-4 ) )
    call a_file%get( "this real", da, my_name, default = 432.5d0 )
    call assert( 775363514, almost_equal( da, 23.4d0 ) )
    call a_file%get( "not there", da, my_name, default = 643.78d0, found = found )
    call assert( 887681859, .not. found )
    call assert( 435049706, da .eq. 643.78d0 )
    call a_file%get( "also not there", da, my_name, found = found )
    call assert( 228989759, .not. found )

    ! get boolean

    call a_file%get( "is it?", la, my_name )
    call assert( 807245669, .not. la )
    call a_file%get( "is it really?", la, my_name, found = found )
    call assert( 405734529, found )
    call assert( 630371219, la )
    call a_file%get( "is it?", la, my_name, default = .false. )
    call assert( 511335328, .not. la )
    call a_file%get( "not there", la, my_name, default = .true., found = found )
    call assert( 672869152, .not. found )
    call assert( 227406840, la )
    call a_file%get( "also not there", la, my_name, found = found )
    call assert( 344666877, .not. found )

    ! get double array

    call a_file%get( "a bunch of doubles", daa, my_name )
    call assert( 302144795, size( daa ) .eq. 4 )
    call assert( 421632981, daa(1) .eq. 12.5_musica_dk )
    call assert( 976054865, daa(2) .eq. 13.2_musica_dk )
    call assert( 465584153, daa(3) .eq. 72.5_musica_dk )
    call assert( 972696092, daa(4) .eq. -142.64_musica_dk )
    call a_file%get( "another bunch of doubles", daa, my_name, found = found )
    call assert( 707754126, found )
    call assert( 460772141, size( daa ) .eq. 2 )
    call assert( 511893154, daa(1) .eq. 52.3_musica_dk )
    call assert( 401480343, daa(2) .eq. 0.0_musica_dk )
    allocate( dab( 2 ) )
    dab(1) = 83.32_musica_dk
    dab(2) = -64.23_musica_dk
    call a_file%get( "a bunch of doubles", daa, my_name, default = dab )
    call assert( 607417634, size( daa ) .eq. 4 )
    call assert( 826790017, daa(1) .eq. 12.5_musica_dk )
    call assert( 656633113, daa(2) .eq. 13.2_musica_dk )
    call assert( 204000960, daa(3) .eq. 72.5_musica_dk )
    call assert( 998852455, daa(4) .eq. -142.64_musica_dk )
    call a_file%get( "not there", daa, my_name, default = dab, found = found )
    call assert( 369345852, .not. found )
    call assert( 611515825, size( daa ) .eq. 2 )
    call assert( 441358921, daa(1) .eq.  83.32_musica_dk )
    call assert( 836152515, daa(2) .eq. -64.23_musica_dk )
    call a_file%get( "also not there", daa, my_name, found = found )
    call assert( 242877146, .not. found )

    ! get string array

    call a_file%get( "a bunch of strings", saa, my_name )
    call assert( 215424987, size( saa ) .eq. 3 )
    call assert( 834855271, saa(1) .eq. "foo" )
    call assert( 376958811, saa(2) .eq. "bar" )
    call assert( 884070750, saa(3) .eq. "foobar" )
    call a_file%get( "another bunch of strings", saa, my_name, found = found )
    call assert( 821420179, found )
    call assert( 533680623, size( saa ) .eq. 2 )
    call assert( 875899965, saa(1) .eq. "boo" )
    call assert( 135528256, saa(2) .eq. "far" )
    allocate( sab(2) )
    sab(1) = "default 1"
    sab(2) = "default 2"
    call a_file%get( "a bunch of strings", saa, my_name, default = sab )
    call assert( 802720780, size( saa ) .eq. 3 )
    call assert( 632563876, saa(1) .eq. "foo" )
    call assert( 127357471, saa(2) .eq. "bar" )
    call assert( 857200566, saa(3) .eq. "foobar" )
    call a_file%get( "not there", saa, my_name, default = sab, found = found )
    call assert( 801267541, .not. found )
    call assert( 513527985, size( saa ) .eq. 2 )
    call assert( 120639925, saa(1) .eq. "default 1" )
    call assert( 792644461, saa(2) .eq. "default 2" )
    call a_file%get( "also not there", saa, my_name, found = found )
    call assert( 354743196, .not. found ) 

    ! add config

    a = '{ "some int" : 1263 }'
    b = '{ "some real" : 14.3, "some string" : "foo" }'
    call a%add( "sub props", b, my_name )
    call b%add( "some string", "bar", my_name )
    call b%get( "some string", sa, my_name )
    call assert( 384683830, sa .eq. "bar" )
    call a%get( "some int", ia, my_name )
    call assert( 762415504, ia .eq. 1263 )
    call a%get( "sub props", c, my_name )
    call c%get( "some string", sa, my_name )
    call assert( 643379613, sa .eq. "foo" )
    call c%get( "some real", da, my_name )
    call assert( 252397087, almost_equal( da, 14.3d0 ) )

    ! add char array

    call a%add( "new char array", "new char array value", my_name )
    call a%get( "some int", ia, my_name )
    call assert( 575490332, ia .eq. 1263 )
    call a%get( "new char array", sa, my_name )
    call assert( 110876326, sa .eq. "new char array value" )

    ! add string

    sa = "new string value"
    call a%add( "new string", sa, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 428870436, ia .eq. 1263 )
    call a%get( "new string", sb, my_name )
    call assert( 258713532, sb .eq. "new string value" )

    ! add int

    call a%add( "new int", 432, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 601194400, ia .eq. 1263 )
    call a%get( "new int", ia, my_name )
    call assert( 827736624, ia .eq. 432 )

    ! add float

    call a%add( "new float", 12.75, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 313907139, ia .eq. 1263 )
    call a%get( "new float", ra, my_name )
    call assert( 875498864, almost_equal( ra, 12.75 ) )

    ! add double

    call a%add( "new double", 53.6d0, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 470628951, ia .eq. 1263 )
    call a%get( "new double", da, my_name )
    call assert( 468723417, almost_equal( da, 53.60d0 ) )

    ! add logical

    call a%add( "new logical", .true., my_name )
    call a%get( "some int", ia, my_name )
    call assert( 570965443, ia .eq. 1263 )
    call a%get( "new logical", la, my_name )
    call assert( 128861904, la )

    ! add double array

    if( allocated( daa ) ) deallocate( daa )
    if( allocated( dab ) ) deallocate( dab )
    allocate( daa(2) )
    daa(1) = -32.51_musica_dk
    daa(2) = 10.324_musica_dk
    call a%add( "new double array", daa, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 971982271, ia .eq. 1263 )
    call a%get( "new double array", dab, my_name )
    call assert( 456247252, size( dab ) .eq. 2 )
    call assert( 115933444, dab(1) .eq. -32.51_musica_dk )
    call assert( 570471131, dab(2) .eq. 10.324_musica_dk )

    ! add string array

    if( allocated( saa ) ) deallocate( saa )
    if( allocated( sab ) ) deallocate( sab )
    allocate( saa(2) )
    saa(1) = "foo"
    saa(2) = "bar"
    call a%add( "new string array", saa, my_name )
    call a%get( "some int", ia, my_name )
    call assert( 729592789, ia .eq. 1263 )
    call a%get( "new string array", sab, my_name )
    call assert( 225839623, size( sab ) .eq. 2 )
    call assert( 115426812, sab(1) .eq. "foo" )
    call assert( 275055102, sab(2) .eq. "bar" )

    ! assignment

    a = '{ "my favorite int" : 42 }'
    b = a
    call a%add( "my favorite int", 43, my_name )
    call a%get( "my favorite int", ia, my_name )
    call assert( 277177497, ia .eq. 43 )
    call b%get( "my favorite int", ia, my_name )
    call assert( 679211194, ia .eq. 42 )
    sa = '{ "another int" : 532 }'
    c = sa
    call c%get( "another int", ia, my_name )
    call assert( 842650552, ia .eq. 532 )

    ! iterator
    a = '{ "my int" : 2,'//&
        '  "my real" : 4.2,'//&
        '  "my double" : 5.2,'//&
        '  "my logical" : true,'//&
        '  "my string" : "foo bar",'//&
        '  "my sub config" : { "an int" : 3, "a double" : 6.7 },'//&
        '  "my string array" : [ "foo", "bar", "foobar" ] }'
    call assert( 494127713, a%number_of_children( ) .eq. 7 )
    iterator => a%get_iterator( )
    call assert( 909667855, iterator%next( ) )
    call assert( 432671110, a%key( iterator ) .eq. "my int" )
    call a%get( iterator, ia, my_name )
    call assert( 227587000, ia .eq. 2 )
    call assert( 217058386, iterator%next( ) )
    call a%get( iterator, ra, my_name )
    call assert( 391026358, almost_equal( ra, 4.2 ) )
    call assert( 270084933, iterator%next( ) )
    call a%get( iterator, da, my_name )
    call assert( 384308812, almost_equal( da, 5.2d0 ) )
    call assert( 826412351, iterator%next( ) )
    call a%get( iterator, la, my_name )
    call assert( 258103080, la )
    call assert( 147690269, iterator%next( ) )
    call a%get( iterator, sa, my_name )
    call assert( 361110121, sa .eq. "foo bar" )
    call assert( 468164159, iterator%next( ) )
    call a%get( iterator, b, my_name )
    call b%get( "a double", da, my_name )
    call assert( 749186169, almost_equal( da, 6.7d0 ) )
    call b%get( "an int", ia, my_name )
    call assert( 915984300, ia .eq. 3 )
    call assert( 182782432, iterator%next( ) )
    call a%get( iterator, saa, my_name )
    call assert( 902549208, saa(1) .eq. "foo" )
    call assert( 334239937, saa(2) .eq. "bar" )
    call assert( 164083033, saa(3) .eq. "foobar" )
    call assert( 441293975, .not. iterator%next( ) )
    call iterator%reset( )
    call assert( 102885701, iterator%next( ) )
    call a%get( iterator, ia, my_name )
    call assert( 162629794, ia .eq. 2 )
    deallocate( iterator )

    ! sequence iterator
    a = '[ 2, 3, "foo", { "bar": 4 } ]'
    call assert( 443487346, a%number_of_children( ) .eq. 4 )
    iterator => a%get_iterator( )
    call assert( 447298414, iterator%next( ) )
    call a%get( iterator, ia, my_name )
    call assert( 612191011, ia .eq. 2 )
    call assert( 442034107, iterator%next( ) )
    call a%get( iterator, ia, my_name )
    call assert( 889401953, ia .eq. 3 )
    call assert( 101720299, iterator%next( ) )
    call a%get( iterator, sa, my_name )
    call assert( 214038644, sa .eq. "foo" )
    call assert( 661406490, iterator%next( ) )
    call a%get( iterator, b, my_name )
    call b%get( "bar", ia, my_name )
    call assert( 208774337, ia .eq. 4 )
    call assert( 103625833, .not. iterator%next( ) )
    call iterator%reset( )
    call assert( 685807284, iterator%next( ) )
    call a%get( iterator, ia, my_name )
    call assert( 233175131, ia .eq. 2 )
    call assert( 410501876, iterator%next( ) )
    call a%get( iterator, ia, my_name )
    call assert( 857869722, ia .eq. 3 )
    call assert( 122762320, iterator%next( ) )
    call a%get( iterator, sa, my_name )
    call assert( 852605415, sa .eq. "foo" )
    call assert( 682448511, iterator%next( ) )
    call a%get( iterator, b, my_name )
    call b%get( "bar", ia, my_name )
    call assert( 229816358, ia .eq. 4 )
    call assert( 124667854, .not. iterator%next( ) )
    deallocate( iterator )

    ! empty object iterator
    a = ""
    call assert( 753171096, a%number_of_children( ) .eq. 0 )
    iterator => a%get_iterator( )
    call assert( 358377502, .not. iterator%next( ) )
    deallocate( iterator )

    ! merging
    a = '{ "a key" : 12,'//&
        '  "another key" : 14.2,'//&
        '  "sub stuff" : {'//&
        '    "orig key" : 72'//&
        '  },'//&
        '  "yet another key" : "hi" }'
    b = '{ "a new key" : true, '//&
        '  "sub stuff" : {'//&
        '    "new key" : "foo"'//&
        '  },'//&
        '  "another new key" : 51 }'
    call a%merge_in( b, my_name )
    call a%get( "a key", ia, my_name )
    call assert( 111746421, ia .eq. 12 )
    call a%get( "another key", da, my_name )
    call assert( 838230743, almost_equal( da, 14.2d0 ) )
    call a%get( "yet another key", sa, my_name )
    call assert( 259845153, sa .eq. "hi" )
    call a%get( "a new key", la, my_name )
    call assert( 879275437, la )
    call a%get( "another new key", ia, my_name )
    call assert( 756880773, ia .eq. 51 )
    call a%get( "sub stuff", c, my_name )
    call c%get( "orig key", ia, my_name )
    call assert( 172568249, ia .eq. 72 )
    call c%get( "new key", sa, my_name )
    call b%get( "a new key", la, my_name )
    call assert( 816624866, la )
    call b%get( "another new key", ia, my_name )
    call assert( 877822198, ia .eq. 51 )
    call b%get( "sub stuff", c, my_name )
    call c%get( "new key", sa, my_name )
    call assert( 597877923, sa .eq. "foo" )
    call c%get( "orig key", ia, my_name, found = found )
    call assert( 933379719, .not. found )
    call b%get( "a key", ia, my_name, found = found )
    call assert( 597164102, .not. found )
    call b%get( "another key", da, my_name, found = found )
    call assert( 293082976, .not. found )
    call b%get( "yet another key", sa, my_name, found = found )
    call assert( 907248953, .not. found )

    ! get and set config array
    allocate( orig_array( 3 ) )
    call orig_array( 1 )%empty( )
    call orig_array( 1 )%add( "a key", "a", my_name )
    call orig_array( 1 )%add( "same key", "same value", my_name )
    call orig_array( 2 )%empty( )
    call orig_array( 2 )%add( "b key", "b", my_name )
    call orig_array( 2 )%add( "same key", "same value", my_name )
    call orig_array( 3 )%empty( )
    call orig_array( 3 )%add( "c key", "c", my_name )
    call orig_array( 3 )%add( "same key", "same value", my_name )
    call array%empty( )
    call array%add( "my array", orig_array, my_name )
    deallocate( orig_array )
    call array%get( "my array", dest_array, my_name )
    call assert( 706554286, allocated( dest_array ) )
    call assert( 874805656, size( dest_array ) .eq. 3 )
    call dest_array( 1 )%get( "a key", sa, my_name )
    call assert( 308401919, sa .eq. "a" )
    call dest_array( 2 )%get( "b key", sa, my_name )
    call assert( 475200050, sa .eq. "b" )
    call dest_array( 3 )%get( "c key", sa, my_name )
    call assert( 640092647, sa .eq. "c" )
    deallocate( dest_array )
    call array%get( "my array", b, my_name )
    iterator => b%get_iterator( )
    call assert( 259072462, b%number_of_children( ) .eq. 3 )
    do while( iterator%next( ) )
      call b%get( iterator, a, my_name )
      call a%get( "same key", sa, my_name )
      call assert( 322175328, sa .eq. "same value" )
    end do
    deallocate( iterator )

    ! string assignment

    a = '{ "foo": 12, "bar": false }'
    sa = a
    b = sa
    call assert( 618824101, b%number_of_children( ) .eq. 2 )
    call b%get( "foo", ia, my_name )
    call assert( 733047980, ia .eq. 12 )
    call b%get( "bar", la, my_name )
    call assert( 787527766, .not. la )

    ! JSON/YAML validation
    a = '{ "a reqd key": 12.3,'// &
        '  "an optional key": "abcd",'// &
        '  "another reqd key": false,'// &
        '  "__a user key": { "foo": "bar" } }'
    if( allocated( saa ) ) deallocate( saa )
    if( allocated( sab ) ) deallocate( sab )
    allocate( saa( 2 ) )
    allocate( sab( 2 ) )
    saa(1) = "a reqd key"
    saa(2) = "another reqd key"
    sab(1) = "an optional key"
    sab(2) = "another optional key"
    call assert( 645591305, a%validate( saa, sab ) )
    deallocate( saa )
    allocate( saa( 1 ) )
    saa(1) = "a reqd key"
    call assert( 264571120, .not. a%validate( saa, sab ) )

  end subroutine test_config_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test the \c config_t example code
  subroutine config_example( )

use musica_config,                   only : config_t
use musica_constants,                only : musica_dk, musica_ik
use musica_iterator,                 only : iterator_t
use musica_string,                   only : string_t
 
character(len=*), parameter :: my_name = "config file example"
type(config_t) :: main_config, sub_config, sub_real_config
real(musica_dk) :: my_real
integer(musica_ik) :: my_int
type(string_t) :: my_string
class(iterator_t), pointer :: iter
logical :: found

call main_config%from_file( 'test/data/config_example.json' )

! this would fail with an error if 'a string' is not found
call main_config%get( "a string", my_string, my_name )
write(*,*) "a string value: ", my_string%val_
 
! add the found argument to avoid failure if the pair is not found
call main_config%get( "my int", my_int, my_name, found = found )
if( found ) then
  write(*,*) "my int value: ", my_int
else
  write(*,*) "'my int' was not found"
end if
 
! when you get a subset of the properties, a new config_t object is
! created containing the subset data. The two config_t objects are
! independent of one another after this point.
call main_config%get( "other props", sub_config, my_name )
call sub_config%get( "an int", my_int, my_name )
write(*,*) "other props->an int value: ", my_int
 
! you can iterate over a set of key-value pairs. but remember that
! the order is always arbitrary. you also must provide the right type
! of variable for the values.
call main_config%get( "real props", sub_real_config, my_name )
iter => sub_real_config%get_iterator( )
do while( iter%next( ) )
  my_string = sub_real_config%key( iter )
  call sub_real_config%get( iter, my_real, my_name )
  write(*,*) my_string%val_, " value: ", my_real
end do
 
! you can also get the number of child objects before iterating over
! them, if you want to allocate an array or something first
write(*,*) "number of children: ", sub_real_config%number_of_children( )
 
! you can add key-value pairs with the add function
call main_config%add( "my new int", 43, my_name )
call main_config%get( "my new int", my_int, my_name )
write(*,*) "my new int value: ", my_int
 
! clean up memory
deallocate( iter )

  end subroutine config_example

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_config
