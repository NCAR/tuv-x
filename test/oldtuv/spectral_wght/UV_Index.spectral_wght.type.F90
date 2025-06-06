! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This uv_index_spectral_wght module

!> The spectral_wght type and related functions
module micm_uv_index_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: uv_index_spectral_wght_t

  !> Calculator for uv_index_spectral_wght
  type, extends(base_spectral_wght_t) :: uv_index_spectral_wght_t
    !> The spectral wght array
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type uv_index_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> uv_index type
    class(uv_index_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'uv_index calculate: '

    spectral_wght = 40._musica_dk*fery( this%mdl_lambda_center )

  end function run

  FUNCTION fery(w)

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the action spectrum value for erythema at a given wavelength   =*
!* Webb, A.R., H. Slaper, P. Koepke, and A. W. Schmalwieser, 
!* Know your standard: Clarifying the CIE erythema action spectrum,
!* Photochem. Photobiol. 87, 483-486, 2011.
!=  Value at 300 nm = 0.6486                                                 =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  W - REAL, wavelength (nm)                                             (I)=*
!-----------------------------------------------------------------------------*

      REAL(musica_dk), intent(in) :: w(:)

      REAL(musica_dk) :: fery(size(w))

      where( w <= 298._musica_dk )
        fery = 1._musica_dk
      elsewhere( w > 298._musica_dk .and. w <= 328._musica_dk )
        fery = 10._musica_dk**(.094_musica_dk*(298._musica_dk - w))
      elsewhere( w > 328._musica_dk .and. w <= 400._musica_dk )
        fery = 10._musica_dk**(.015_musica_dk*(140._musica_dk - w))
      elsewhere
        fery = 1.E-36_musica_dk
      endwhere

  END FUNCTION fery

end module micm_uv_index_spectral_wght_type
