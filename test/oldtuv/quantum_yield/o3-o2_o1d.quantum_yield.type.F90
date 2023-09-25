! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This o3+hv->o2+o1d quantum yield module

!> The o3+hv->o2+o1d quantum yield type and related functions
module micm_o3_o2_o1d_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: o3_o2_o1d_quantum_yield_t

  integer(musica_ik), parameter :: iONE = 1_musica_ik
  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for o3+hv->o2+o1d quantum yield
  type, extends(base_quantum_yield_t) :: o3_o2_o1d_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type o3_o2_o1d_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(o3_o2_o1d_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: a(3)  = (/ 0.8036_musica_dk, 8.9061_musica_dk, 0.1192_musica_dk /)
    real(musica_dk), parameter :: x(3)  = (/ 304.225_musica_dk, 314.957_musica_dk, 310.737_musica_dk /)
    real(musica_dk), parameter :: om(3) = (/ 5.576_musica_dk, 6.601_musica_dk, 2.187_musica_dk /)

    character(len=*), parameter :: Iam = 'o3+hv->o2+o1d calculate: '

    integer(musica_ik) :: wNdx
    real(musica_dk) ::  kt, q1, q2, T300, lambda
    real(musica_dk) ::  qfac1, qfac2

    write(*,*) Iam,'entering'

!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
! function to calculate the quantum yield O3 + hv -> O(1D) + O2,             =*
! according to:                                                             
! Matsumi, Y., F. J. Comes, G. Hancock, A. Hofzumanhays, A. J. Hynes,
! M. Kawasaki, and A. R. Ravishankara, QUantum yields for production of O(1D)
! in the ultraviolet photolysis of ozone:  Recommendation based on evaluation
! of laboratory data, J. Geophys. Res., 107, 10.1029/2001JD000510, 2002.
!-----------------------------------------------------------------------------*

      associate( w => this%mdl_lambda_center, Temp => environment%temperature )

      kt = 0.695_musica_dk * Temp
      q1 = rONE
      q2 = exp( -825.518_musica_dk/kt )
      qfac1 = q1/(q1 + q2)
      qfac2 = q2/(q1 + q2)

      T300 = Temp/300._musica_dk

      quantum_yield = rZERO
      where( w(:) <= 305._musica_dk )
        quantum_yield = 0.90_musica_dk
      elsewhere( w(:) > 328._musica_dk .and. w(:) <= 340._musica_dk )
        quantum_yield = 0.08_musica_dk
      endwhere
      do wNdx = iONE,size(w)
        lambda = w(wNdx)
        if( lambda > 305._musica_dk .and. lambda <= 328._musica_dk ) then
          quantum_yield(wNdx) = 0.0765_musica_dk &
          + a(1)          *qfac1*EXP( -((x(1) - lambda)/om(1))**4 ) &
          + a(2)*T300*T300*qfac2*EXP( -((x(2) - lambda)/om(2))**2 ) &
          + a(3)*T300**1.5_musica_dk*EXP( -((x(3) - lambda)/om(3))**2 )
        endif
      enddo

      end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_o3_o2_o1d_quantum_yield_type
