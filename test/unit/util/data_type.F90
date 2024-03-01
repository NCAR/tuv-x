! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests for the musica_data_type module

!> Tests for the data_type_t type
program test_util_data_type

  use musica_assert
  use musica_constants,                only : musica_dk, musica_ik,           &
                                              musica_lk, musica_rk
  use musica_data_type
#ifdef MUSICA_USE_OPENMP
  use omp_lib
#endif

  implicit none

#ifdef MUSICA_USE_OPENMP
  write(*,*) "Testing with ", omp_get_max_threads( ), " threads"
#else
  write(*,*) "Testing without OpenMP support"
#endif
  !$omp parallel
  call test_data_type_t( )
  !$omp end parallel

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Test data_type_t functionality
  subroutine test_data_type_t( )

    use musica_string,                 only : string_t

    type(data_type_t) :: dt_b, dt_d, dt_f, dt_i
    type(string_t) :: temp_str
    class(*), pointer :: temp_val

    ! constructors
    dt_b = data_type_t( "BooLean" )
    dt_d = data_type_t( "DOUBLE"  )
    dt_f = data_type_t( "float"   )
    temp_str = "iNtEgEr"
    dt_i = data_type_t( temp_str  )

    ! equality
    call assert( 322458369, dt_b .eq. kBoolean )
    call assert( 824306001, dt_b .ne. kTypeUnknown )
    call assert( 596310538, dt_b .ne. kFloat )
    call assert( 433323475, dt_b .ne. kDouble )
    call assert( 375484916, dt_b .ne. kInteger )
    call assert( 765014203, dt_d .eq. kDouble )
    call assert( 307117743, dt_f .eq. kFloat )
    call assert( 926548027, dt_i .eq. kInteger )
    call assert( 868709468, dt_d .ne. dt_b )
    call assert( 423247156, dt_f .ne. dt_i )

    ! names
    call assert( 300400197, dt_b%name( ) .eq. "boolean" )
    call assert( 467198328, dt_d%name( ) .eq. "double"  )
    call assert( 297041424, dt_f%name( ) .eq. "float"   )
    call assert( 126884520, dt_i%name( ) .eq. "integer" )

  end subroutine test_data_type_t

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end program test_util_data_type
