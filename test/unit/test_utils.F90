module tuvx_test_utils

  implicit none
  public

  interface check_values
    procedure :: check_values_1D, check_values_2D, check_values_1D_no_code,   &
                 check_values_2D_no_code
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_1D( code, results, expected_results,                &
      relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    integer,       intent(in) :: code
    real(kind=dk), intent(in) :: results(:)
    real(kind=dk), intent(in) :: expected_results(:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_elem

    call assert( code, size( results ) == size( expected_results ) )
    do i_elem = 1, size( results )
      call assert( code, almost_equal(                                       &
        results( i_elem ),                                                   &
        expected_results( i_elem ),                                          &
        relative_tolerance))
    end do

  end subroutine check_values_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_2D( code, results, expected_results,                &
      relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    integer,       intent(in) :: code
    real(kind=dk), intent(in) :: results(:,:)
    real(kind=dk), intent(in) :: expected_results(:,:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_level, i_wavelength

    call assert( code, &
      size( results, dim = 1 ) == size( expected_results, dim = 1) )
    call assert( code, &
      size( results, dim = 2 ) == size( expected_results, dim = 2) )
    do i_wavelength = 1, size( results, dim = 2 )
      do i_level = 1, size( results, dim = 1 )
        call assert( code, almost_equal(                                      &
          results( i_level, i_wavelength ),                                   &
          expected_results( i_level, i_wavelength ),                          &
          relative_tolerance))
      end do
    end do

  end subroutine check_values_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_1D_no_code( results, expected_results,              &
      relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    real(kind=dk), intent(in) :: results(:)
    real(kind=dk), intent(in) :: expected_results(:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_elem

    call assert( 966913792, size( results ) == size( expected_results ) )
    do i_elem = 1, size( results )
      call assert( 796756888, almost_equal(                                  &
        results( i_elem ),                                                   &
        expected_results( i_elem ),                                          &
        relative_tolerance))
    end do

  end subroutine check_values_1D_no_code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine check_values_2D_no_code( results, expected_results,              &
      relative_tolerance )

    use musica_assert,                 only : assert, almost_equal
    use musica_constants,              only : dk => musica_dk

    real(kind=dk), intent(in) :: results(:,:)
    real(kind=dk), intent(in) :: expected_results(:,:)
    real(kind=dk), intent(in) :: relative_tolerance

    integer :: i_level, i_wavelength

    call assert( 516187173, &
      size( results, dim = 1 ) == size( expected_results, dim = 1) )
    call assert( 963555019, &
      size( results, dim = 2 ) == size( expected_results, dim = 2) )
    do i_wavelength = 1, size( results, dim = 2 )
      do i_level = 1, size( results, dim = 1 )
        call assert( 228447617, almost_equal(                                 &
          results( i_level, i_wavelength ),                                   &
          expected_results( i_level, i_wavelength ),                          &
          relative_tolerance))
      end do
    end do

  end subroutine check_values_2D_no_code

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_test_utils
