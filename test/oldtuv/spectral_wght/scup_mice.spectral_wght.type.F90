! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This scup_mice_spectral_wght module

!> The spectral_wght type and related functions
module micm_scup_mice_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: scup_mice_spectral_wght_t

  !> Calculator for scup_mice_spectral_wght
  type, extends(base_spectral_wght_t) :: scup_mice_spectral_wght_t
    !> The spectral wght array
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type scup_mice_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> scup_mice type
    class(scup_mice_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    real(kind=musica_dk)        :: factor(1)
    character(len=*), parameter :: Iam = 'scup_mice calculate: '

    write(*,*) Iam,'entering'

    factor = 1._musica_dk/sw_futr( (/300._musica_dk/) )
    spectral_wght = sw_futr( this%mdl_lambda_center ) * factor(1)

    write(*,*) Iam,'exiting'

  end function run

  FUNCTION sw_futr(w)
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Calculate the action spectrum value for skin cancer of albino hairless   =*
!=  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
!=  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
!=  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
!=  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
!=  pp. 53-60, 1993                                                          =*
!=  (Action spectrum for carcinomas)                                         =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  W  - REAL, wavelength (nm)                                            (I)=*
!-----------------------------------------------------------------------------*

      REAL(musica_dk), intent(in) :: w(:)

      REAL(musica_dk) :: sw_futr(size(w))

      real(musica_dk), allocatable :: t1(:), t2(:), t3(:), t4(:), t5(:)
      real(musica_dk), allocatable :: p(:)

      real(musica_dk), parameter :: a1 = -10.91_musica_dk
      real(musica_dk), parameter :: a2 = - 0.86_musica_dk
      real(musica_dk), parameter :: a3 = - 8.60_musica_dk
      real(musica_dk), parameter :: a4 = - 9.36_musica_dk
      real(musica_dk), parameter :: a5 = -13.15_musica_dk

      real(musica_dk), parameter :: x1 = 270._musica_dk
      real(musica_dk), parameter :: x2 = 302._musica_dk
      real(musica_dk), parameter :: x3 = 334._musica_dk
      real(musica_dk), parameter :: x4 = 367._musica_dk
      real(musica_dk), parameter :: x5 = 400._musica_dk

      real(musica_dk), parameter :: b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
      real(musica_dk), parameter :: b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)
      real(musica_dk), parameter :: b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
      real(musica_dk), parameter :: b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
      real(musica_dk), parameter :: b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

      real(musica_dk), parameter :: w1 = a1/b1
      real(musica_dk), parameter :: w2 = a2/b2
      real(musica_dk), parameter :: w3 = a3/b3
      real(musica_dk), parameter :: w4 = a4/b4
      real(musica_dk), parameter :: w5 = a5/b5

      t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
      t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
      t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
      t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
      t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

      p = w1*t1 + w2*t2 + w3*t3 + w4*t4 + w5*t5

      sw_futr  = EXP(p)

      END FUNCTION sw_futr

end module micm_scup_mice_spectral_wght_type
