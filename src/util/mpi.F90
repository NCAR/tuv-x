! Copyright (C) 2007-2021 Barcelona Supercomputing Center and University of
! Illinois at Urbana-Champaign
! SPDX-License-Identifier: MIT
module musica_mpi
  ! Wrapper functions for MPI.
  !
  ! This module was adapted from CAMP (https://github.com/open-atmos/camp).
  !
  ! All of these functions can be called irrespective of whether MPI
  ! support was compiled in or not. If MPI support is not enabled then
  ! they do the obvious trivial thing (normally nothing).

#ifdef MUSICA_USE_MPI
  use mpi
#endif

  use musica_constants,                only : dp => musica_dk

  implicit none

  private
  public :: musica_mpi_support, musica_mpi_init, musica_mpi_abort,            &
            musica_mpi_finalize, musica_mpi_barrier, musica_mpi_rank,         &
            musica_mpi_size, musica_mpi_bcast, musica_mpi_pack_size,          &
            musica_mpi_pack, musica_mpi_unpack, MPI_COMM_WORLD

#ifndef MUSICA_USE_MPI
  ! Parameter to make a communicator available when MPI support is not
  ! compiled in (to avoid a lot of preprocessor flags in tests)
  integer, parameter :: MPI_COMM_WORLD = 0
#endif

  integer, parameter :: dc = dp ! kind for double-precision complex numbers

  ! Broadcasts a variable from the primary process to all other processes
  interface musica_mpi_bcast
    procedure :: musica_mpi_bcast_integer
    procedure :: musica_mpi_bcast_string
    procedure :: musica_mpi_bcast_packed
  end interface musica_mpi_bcast

  ! Returns the size of a character buffer needed to pack a given variable
  interface musica_mpi_pack_size
    procedure :: musica_mpi_pack_size_integer
    procedure :: musica_mpi_pack_size_string
    procedure :: musica_mpi_pack_size_real
    procedure :: musica_mpi_pack_size_logical
    procedure :: musica_mpi_pack_size_complex
    procedure :: musica_mpi_pack_size_integer_array
    procedure :: musica_mpi_pack_size_string_array
    procedure :: musica_mpi_pack_size_real_array
    procedure :: musica_mpi_pack_size_real_array_2d
    procedure :: musica_mpi_pack_size_real_array_3d
  end interface musica_mpi_pack_size

  ! Packs the given variable onto a character buffer
  interface musica_mpi_pack
    procedure :: musica_mpi_pack_integer
    procedure :: musica_mpi_pack_string
    procedure :: musica_mpi_pack_real
    procedure :: musica_mpi_pack_logical
    procedure :: musica_mpi_pack_complex
    procedure :: musica_mpi_pack_integer_array
    procedure :: musica_mpi_pack_string_array
    procedure :: musica_mpi_pack_real_array
    procedure :: musica_mpi_pack_real_array_2d
    procedure :: musica_mpi_pack_real_array_3d
  end interface musica_mpi_pack

  ! Unpacks a variable from a character buffer
  interface musica_mpi_unpack
    procedure :: musica_mpi_unpack_integer
    procedure :: musica_mpi_unpack_string
    procedure :: musica_mpi_unpack_real
    procedure :: musica_mpi_unpack_logical
    procedure :: musica_mpi_unpack_complex
    procedure :: musica_mpi_unpack_integer_array
    procedure :: musica_mpi_unpack_string_array
    procedure :: musica_mpi_unpack_real_array
    procedure :: musica_mpi_unpack_real_array_2d
    procedure :: musica_mpi_unpack_real_array_3d
  end interface musica_mpi_unpack

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function musica_mpi_support( )
    ! Whether MPI support is compiled in.

#ifdef MUSICA_USE_MPI
    musica_mpi_support = .true.
#else
    musica_mpi_support = .false.
#endif

  end function musica_mpi_support

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_check_ierr( ierr )
    ! Dies if ``ierr`` is not ok.

    integer, intent(in) :: ierr ! MPI status code

#ifdef MUSICA_USE_MPI
    if( ierr /= MPI_SUCCESS )then
       call musica_mpi_abort(1)
    end if
#endif

  end subroutine musica_mpi_check_ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_init( )
    ! Initialize MPI.

#ifdef MUSICA_USE_MPI
    integer :: ierr

    call mpi_init( ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_abort( status )
    ! Abort the program.

    integer, intent(in) :: status ! Status flag to abort with

#ifdef MUSICA_USE_MPI
    integer :: ierr

    call mpi_abort( MPI_COMM_WORLD, status, ierr )
#else
    call assert( status, .false. )

#endif

  end subroutine musica_mpi_abort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine musica_mpi_finalize( )

    ! Shut down MPI.

#ifdef MUSICA_USE_MPI
    integer :: ierr

    call mpi_finalize( ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Synchronize all processes.
  subroutine musica_mpi_barrier( comm )

    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: ierr

    call mpi_barrier( comm, ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_barrier

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_rank( comm )
    ! Returns the rank of the current process.

    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: rank, ierr

    call mpi_comm_rank( comm, rank, ierr )
    call musica_mpi_check_ierr( ierr )
    musica_mpi_rank = rank
#else
    musica_mpi_rank = 0
#endif

  end function musica_mpi_rank

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_size( comm )
    ! Returns the total number of processes.

    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: size, ierr

    call mpi_comm_size( comm, size, ierr )
    call musica_mpi_check_ierr( ierr )
    musica_mpi_size = size
#else
    musica_mpi_size = 1
#endif

  end function musica_mpi_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_bcast_integer( val, comm )
    ! Broadcast the given value from process 0 to all other processes.

    integer, intent(inout) :: val ! value to broadcast
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast( val, 1, MPI_INTEGER, root, comm, ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_bcast_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_bcast_string( val, comm )
    ! Broadcast the given value from process 0 to all other processes.

    character(len=*), intent(inout) :: val ! value to broadcast
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast( val, len( val ), MPI_CHARACTER, root, comm, ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_bcast_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_bcast_packed( val, comm )
    ! Broadcast the given value from process 0 to all other processes.

    character, intent(inout) :: val(:) ! value to be broadcast
    integer,   intent(in)    :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: root, ierr

    root = 0 ! source of data to broadcast
    call mpi_bcast( val, size( val ), MPI_CHARACTER, root, comm, ierr )
    call musica_mpi_check_ierr( ierr )
#endif

  end subroutine musica_mpi_bcast_packed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_integer( val, comm )
    ! Determines the number of bytes required to pack the given value.

    integer, intent(in) :: val ! value to be packed
    integer, intent(in) :: comm ! MPI communicator

    integer :: ierr

#ifdef MUSICA_USE_MPI

    call mpi_pack_size( 1, MPI_INTEGER, comm,                               &
                        musica_mpi_pack_size_integer, ierr )
    call musica_mpi_check_ierr( ierr )
#else
    musica_mpi_pack_size_integer = 0
#endif

  end function musica_mpi_pack_size_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_real( val, comm )
    ! Determines the number of bytes required to pack the given value.

    real(kind=dp), intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

    integer :: ierr

#ifdef MUSICA_USE_MPI

    call mpi_pack_size( 1, MPI_DOUBLE_PRECISION, comm,                      &
                        musica_mpi_pack_size_real, ierr )
    call musica_mpi_check_ierr( ierr )
#else
    musica_mpi_pack_size_real = 0
#endif

  end function musica_mpi_pack_size_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_string( val, comm )
    ! Determines the number of bytes required to pack the given value.

    character(len=*), intent(in) :: val ! value to be packed
    integer, intent(in) :: comm ! MPI communicator

    integer :: ierr

#ifdef MUSICA_USE_MPI

    call mpi_pack_size( len_trim( val ), MPI_CHARACTER, comm,               &
                        musica_mpi_pack_size_string, ierr )
    call musica_mpi_check_ierr( ierr )
    musica_mpi_pack_size_string = musica_mpi_pack_size_string                 &
         + musica_mpi_pack_size_integer( len_trim( val ), comm )
#else
    musica_mpi_pack_size_string = 0
#endif

  end function musica_mpi_pack_size_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_logical( val, comm )
    ! Determines the number of bytes required to pack the given value.

    logical, intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

    integer :: ierr

#ifdef MUSICA_USE_MPI

    call mpi_pack_size( 1, MPI_LOGICAL, comm,                               &
                        musica_mpi_pack_size_logical, ierr )
    call musica_mpi_check_ierr( ierr )
#else
    musica_mpi_pack_size_logical = 0
#endif

  end function musica_mpi_pack_size_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_complex( val, comm )
    ! Determines the number of bytes required to pack the given value.

    complex(kind=dc), intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

    integer :: ierr

#ifdef MUSICA_USE_MPI

    call mpi_pack_size( 1, MPI_DOUBLE_COMPLEX, comm,                        &
                        musica_mpi_pack_size_complex, ierr )
    call musica_mpi_check_ierr( ierr )
#else
    musica_mpi_pack_size_complex = 0
#endif

  end function musica_mpi_pack_size_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_integer_array( val, comm )
    ! Determines the number of bytes required to pack the given value.

    integer, allocatable, intent(in) :: val(:) ! value to be packed
    integer, intent(in) :: comm ! MPI communicator

    integer :: total_size, ierr

#ifdef MUSICA_USE_MPI
    logical :: is_allocated


    total_size = 0
    is_allocated = allocated( val )
    if( is_allocated ) then
       call mpi_pack_size( size( val ), MPI_INTEGER, comm, total_size, ierr )
       call musica_mpi_check_ierr( ierr )
       total_size = total_size +                                              &
                    musica_mpi_pack_size_integer( size( val ), comm )
    end if
    total_size = total_size +                                                 &
                 musica_mpi_pack_size_logical( is_allocated, comm )
#else
    total_size = 0
#endif

    musica_mpi_pack_size_integer_array = total_size

  end function musica_mpi_pack_size_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_real_array( val, comm )
    ! Determines the number of bytes required to pack the given value.

    real(kind=dp), allocatable, intent(in) :: val(:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

    integer :: total_size, ierr

#ifdef MUSICA_USE_MPI
    logical :: is_allocated


    total_size = 0
    is_allocated = allocated( val )
    if( is_allocated ) then
       call mpi_pack_size( size( val ), MPI_DOUBLE_PRECISION, comm,         &
                           total_size, ierr )
       call musica_mpi_check_ierr( ierr )
       total_size = total_size +                                              &
                    musica_mpi_pack_size_integer( size( val ), comm )
    end if
    total_size = total_size +                                                 &
                 musica_mpi_pack_size_logical( is_allocated, comm )
#else
    total_size = 0
#endif

    musica_mpi_pack_size_real_array = total_size

  end function musica_mpi_pack_size_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_string_array( val, comm )
    ! Determines the number of bytes required to pack the given value.

    character(len=*), allocatable, intent(in) :: val(:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

    integer :: i, total_size
#ifdef MUSICA_USE_MPI
    logical :: is_allocated


    is_allocated = allocated( val )
    if( is_allocated ) then
       total_size = musica_mpi_pack_size_integer( size( val ), comm )
       do i = 1, size( val )
          total_size = total_size +                                           &
            musica_mpi_pack_size_string( val( i ), comm )
       end do
    end if
    total_size = total_size +                                                 &
                 musica_mpi_pack_size_logical( is_allocated, comm )
    musica_mpi_pack_size_string_array = total_size
#else
    total_size = 0
#endif

  end function musica_mpi_pack_size_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_real_array_2d( val, comm )
    ! Determines the number of bytes required to pack the given value.

    real(kind=dp), allocatable, intent(in) :: val(:,:) ! value to pack
    integer, intent(in) :: comm ! MPI Communicator

    integer :: total_size, ierr

#ifdef MUSICA_USE_MPI
    logical :: is_allocated


    total_size = 0
    is_allocated = allocated( val )
    if( is_allocated ) then
       call mpi_pack_size( size( val ), MPI_DOUBLE_PRECISION, comm,         &
                           total_size, ierr )
       call musica_mpi_check_ierr( ierr )
       total_size = total_size                                                &
            + musica_mpi_pack_size_integer( size( val, 1 ), comm )          &
            + musica_mpi_pack_size_integer( size( val, 2 ), comm )
    end if
    total_size = total_size +                                                 &
                 musica_mpi_pack_size_logical( is_allocated, comm )
#else
    total_size = 0
#endif

    musica_mpi_pack_size_real_array_2d = total_size

  end function musica_mpi_pack_size_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function musica_mpi_pack_size_real_array_3d( val, comm )
    ! Determines the number of bytes required to pack the given value.

    real(kind=dp), allocatable, intent(in) :: val(:,:,:) ! value to pack
    integer, intent(in) :: comm ! MPI Communicator

    integer :: total_size, ierr

#ifdef MUSICA_USE_MPI
    logical :: is_allocated


    total_size = 0
    is_allocated = allocated( val )
    if( is_allocated ) then
       call mpi_pack_size( size( val ), MPI_DOUBLE_PRECISION, comm,           &
                           total_size, ierr )
       call musica_mpi_check_ierr( ierr )
       total_size = total_size                                                &
            + musica_mpi_pack_size_integer( size( val, 1 ), comm )            &
            + musica_mpi_pack_size_integer( size( val, 2 ), comm )            &
            + musica_mpi_pack_size_integer( size( val, 3 ), comm )
    end if
    total_size = total_size +                                                 &
                 musica_mpi_pack_size_logical( is_allocated, comm )
#else
    total_size = 0
#endif

    musica_mpi_pack_size_real_array_3d = total_size

  end function musica_mpi_pack_size_real_array_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_integer( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! current buffer position
    integer, intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_pack( val, 1, MPI_INTEGER, buffer, size( buffer ),               &
                   position, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 913495993,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_integer( val, comm ) )
#endif

  end subroutine musica_mpi_pack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_real( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! current buffer position
    real(kind=dp), intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_pack( val, 1, MPI_DOUBLE_PRECISION, buffer, size( buffer ),      &
                   position, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 395354132,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real( val, comm ) )
#endif

  end subroutine musica_mpi_pack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_string( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    character(len=*), intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, length, ierr


    prev_position = position
    length = len_trim( val )
    call musica_mpi_pack_integer( buffer, position, length, comm )
    call mpi_pack( val, length, MPI_CHARACTER, buffer, size( buffer ),        &
                   position, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 607212018,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_string( val, comm ) )
#endif

  end subroutine musica_mpi_pack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_logical( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    logical, intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_pack( val, 1, MPI_LOGICAL, buffer, size( buffer ),               &
                   position, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 104535200,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_logical( val, comm ) )
#endif

  end subroutine musica_mpi_pack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_complex( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    complex(kind=dc), intent(in) :: val ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_pack( val, 1, MPI_DOUBLE_COMPLEX, buffer, size( buffer ),        &
                   position, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 640416372,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_complex( val, comm ) )
#endif

  end subroutine musica_mpi_pack_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_integer_array( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    integer, allocatable, intent(in) :: val(:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n, ierr
    logical :: is_allocated


    prev_position = position
    is_allocated = allocated( val )
    call musica_mpi_pack_logical( buffer, position, is_allocated, comm )
    if( is_allocated ) then
       n = size( val )
       call musica_mpi_pack_integer( buffer, position, n, comm )
       call mpi_pack( val, n, MPI_INTEGER, buffer, size( buffer ),            &
                      position, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 698601296,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_integer_array( val, comm ) )
#endif

  end subroutine musica_mpi_pack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_real_array( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(in) :: val(:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n, ierr
    logical :: is_allocated


    prev_position = position
    is_allocated = allocated( val )
    call musica_mpi_pack_logical( buffer, position, is_allocated, comm )
    if( is_allocated ) then
       n = size( val )
       call musica_mpi_pack_integer( buffer, position, n, comm )
       call mpi_pack( val, n, MPI_DOUBLE_PRECISION, buffer, size( buffer ),   &
                      position, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 825718791,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real_array( val, comm ) )
#endif

  end subroutine musica_mpi_pack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_string_array( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    character(len=*), allocatable, intent(in) :: val(:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, i, n
    logical :: is_allocated


    prev_position = position
    is_allocated = allocated( val )
    call musica_mpi_pack_logical( buffer, position, is_allocated, comm )
    if( is_allocated) then
       n = size( val )
       call musica_mpi_pack_integer( buffer, position, n, comm )
       do i = 1, n
          call musica_mpi_pack_string( buffer, position, val( i ), comm )
       end do
    end if
    call assert( 630900704,                                                   &
                 position - prev_position <= &
                 musica_mpi_pack_size_string_array( val, comm ) )
#endif

  end subroutine musica_mpi_pack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_real_array_2d( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(in) :: val(:,:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n1, n2, ierr
    logical :: is_allocated


    prev_position = position
    is_allocated = allocated( val )
    call musica_mpi_pack_logical( buffer, position, is_allocated, comm )
    if( is_allocated ) then
       n1 = size( val, 1 )
       n2 = size( val, 2 )
       call musica_mpi_pack_integer( buffer, position, n1, comm )
       call musica_mpi_pack_integer( buffer, position, n2, comm )
       call mpi_pack( val, n1 * n2, MPI_DOUBLE_PRECISION, buffer,             &
                      size( buffer ), position, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 567349745,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real_array_2d( val, comm ) )
#endif

  end subroutine musica_mpi_pack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_pack_real_array_3d( buffer, position, val, comm )
    ! Packs the given value into the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(in) :: val(:,:,:) ! value to pack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n1, n2, n3, ierr
    logical :: is_allocated


    prev_position = position
    is_allocated = allocated( val )
    call musica_mpi_pack_logical( buffer, position, is_allocated, comm )
    if( is_allocated ) then
       n1 = size( val, 1 )
       n2 = size( val, 2 )
       n3 = size( val, 3 )
       call musica_mpi_pack_integer( buffer, position, n1, comm )
       call musica_mpi_pack_integer( buffer, position, n2, comm )
       call musica_mpi_pack_integer( buffer, position, n3, comm )
       call mpi_pack( val, n1 * n2 * n3, MPI_DOUBLE_PRECISION, buffer,        &
                      size( buffer ), position, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 851684870,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real_array_3d( val, comm ) )
#endif

  end subroutine musica_mpi_pack_real_array_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_integer( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    integer, intent(out) :: val ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_unpack( buffer, size( buffer ), position, val, 1, MPI_INTEGER,   &
                     comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 890243339,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_integer( val, comm ) )
#else
    val = 0
#endif

  end subroutine musica_mpi_unpack_integer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_real( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), intent(out) :: val ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_unpack( buffer, size( buffer ), position, val, 1,                &
                     MPI_DOUBLE_PRECISION, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 570771632,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real( val, comm ) )
#else
    val = real( 0.0, kind = dp )
#endif

  end subroutine musica_mpi_unpack_real

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_string( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    character(len=*), intent(out) :: val ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, length, ierr


    prev_position = position
    call musica_mpi_unpack_integer( buffer, position, length, comm )
    call assert(946399479, length <= len( val ) )
    val = ''
    call mpi_unpack( buffer, size( buffer ), position, val, length,           &
                     MPI_CHARACTER, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 503378058, &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_string( val, comm ) )
#else
    val = ''
#endif

  end subroutine musica_mpi_unpack_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_logical( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    logical, intent(out) :: val ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_unpack( buffer, size( buffer ), position, val, 1, MPI_LOGICAL,   &
                     comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 694750528,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_logical( val, comm ) )
#else
    val = .false.
#endif

  end subroutine musica_mpi_unpack_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_complex( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    complex(kind=dc), intent(out) :: val ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, ierr


    prev_position = position
    call mpi_unpack( buffer, size( buffer ), position, val, 1,                &
                     MPI_DOUBLE_COMPLEX, comm, ierr )
    call musica_mpi_check_ierr( ierr )
    call assert( 969672634,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_complex( val, comm ) )
#else
    val = cmplx( 0 )
#endif

  end subroutine musica_mpi_unpack_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_integer_array( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    integer, allocatable, intent(inout) :: val(:) ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n, ierr
    logical :: is_allocated


    prev_position = position
    call musica_mpi_unpack_logical( buffer, position, is_allocated, comm )
    if( allocated( val ) ) deallocate( val )
    if( is_allocated ) then
       call musica_mpi_unpack_integer( buffer, position, n, comm )
       allocate( val( n ) )
       call mpi_unpack( buffer, size( buffer ), position, val, n, MPI_INTEGER,&
                        comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 565840919,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_integer_array( val, comm ) )
#endif

  end subroutine musica_mpi_unpack_integer_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_real_array( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(inout) :: val(:) ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n, ierr
    logical :: is_allocated


    prev_position = position
    call musica_mpi_unpack_logical( buffer, position, is_allocated, comm )
    if( allocated( val ) ) deallocate( val )
    if( is_allocated ) then
       call musica_mpi_unpack_integer( buffer, position, n, comm )
       allocate( val( n ) )
       call mpi_unpack( buffer, size( buffer ), position, val, n,             &
                        MPI_DOUBLE_PRECISION, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 782875761,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_real_array( val, comm ) )
#endif

  end subroutine musica_mpi_unpack_real_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_string_array( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    character(len=*), allocatable, intent(inout) :: val(:) ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, i, n
    logical :: is_allocated


    prev_position = position
    call musica_mpi_unpack_logical( buffer, position, is_allocated, comm )
    if( allocated( val ) ) deallocate( val )
    if( is_allocated ) then
       call musica_mpi_unpack_integer( buffer, position, n, comm )
       allocate( val( n ) )
       do i = 1, n
          call musica_mpi_unpack_string( buffer, position, val( i ), comm )
       end do
    end if
    call assert( 320065648,                                                   &
                 position - prev_position <=                                  &
                 musica_mpi_pack_size_string_array( val, comm ) )
#endif

  end subroutine musica_mpi_unpack_string_array

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_real_array_2d( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(inout) :: val(:,:) ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n1, n2, ierr
    logical :: is_allocated


    prev_position = position
    call musica_mpi_unpack_logical( buffer, position, is_allocated, comm )
    if( allocated( val ) ) deallocate( val )
    if( is_allocated ) then
       call musica_mpi_unpack_integer( buffer, position, n1, comm )
       call musica_mpi_unpack_integer( buffer, position, n2, comm )
       allocate( val( n1, n2 ) )
       call mpi_unpack( buffer, size( buffer ), position, val, n1 * n2,       &
                                MPI_DOUBLE_PRECISION, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 781681739, position - prev_position                          &
                 <= musica_mpi_pack_size_real_array_2d( val, comm ) )
#endif

  end subroutine musica_mpi_unpack_real_array_2d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine musica_mpi_unpack_real_array_3d( buffer, position, val, comm )
    ! Unpacks the given value from the buffer, advancing position.

    character, intent(inout) :: buffer(:) ! memory buffer
    integer, intent(inout) :: position ! curent buffer position
    real(kind=dp), allocatable, intent(inout) :: val(:,:,:) ! value to unpack
    integer, intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_position, n1, n2, n3, ierr
    logical :: is_allocated


    prev_position = position
    call musica_mpi_unpack_logical( buffer, position, is_allocated, comm )
    if( allocated( val ) ) deallocate( val )
    if( is_allocated ) then
       call musica_mpi_unpack_integer( buffer, position, n1, comm )
       call musica_mpi_unpack_integer( buffer, position, n2, comm )
       call musica_mpi_unpack_integer( buffer, position, n3, comm )
       allocate( val( n1, n2, n3 ) )
       call mpi_unpack( buffer, size( buffer ), position, val, n1 * n2 * n3,  &
                                MPI_DOUBLE_PRECISION, comm, ierr )
       call musica_mpi_check_ierr( ierr )
    end if
    call assert( 162434174, position - prev_position                          &
                 <= musica_mpi_pack_size_real_array_3d( val, comm ) )
#endif

  end subroutine musica_mpi_unpack_real_array_3d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !> Local assert
  subroutine assert( code, condition )

    !> Unique code for the assertion
    integer, intent(in) :: code
    !> Condition to evaluate
    logical, intent(in) :: condition

    character(len=50) :: str_code
    integer, parameter :: kErrorId = 0
    integer, parameter :: kErrorFileId = 10

    if( .not. condition ) then
      write(str_code,'(i30)') code
      write(kErrorId,*) "ERROR (Musica-"//trim( adjustl( str_code ) )//"): "  &
                        //"assertion failed"
      open( unit = kErrorFileId, file = "error.json", action = "WRITE" )
      write(kErrorFileId,'(A)') '{'
      write(kErrorFileId,'(A)') '  "code" : "'//trim( adjustl( str_code ) )//'",'
      write(kErrorFileId,'(A)') '  "message" : "assertion failed"'
      write(kErrorFileId,'(A)') '}'
      close(kErrorFileId)
      stop 3
    end if

  end subroutine assert

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module musica_mpi
