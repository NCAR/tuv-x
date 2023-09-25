! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch2chcho+hv->products quantum yield module

!> The ch2chcho+hv->prodcuts quantum yield type and related functions
module micm_ch2chcho_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch2chcho_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch2chcho+hv->oh+h quantum yield
  type, extends(base_quantum_yield_t) :: ch2chcho_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ch2chcho_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ch2chcho_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: phiL = .004_musica_dk
    real(musica_dk), parameter :: phiU = .086_musica_dk

    character(len=*), parameter :: Iam = 'ch2chcho+hv->products calculate: '
    real(musica_dk) :: phi0

    write(*,*) Iam,'entering'

    associate( M => environment%number_density_air )
      if( M > 2.6e19_musica_dk ) then
        quantum_yield = phiL
      else
        if( M <= 8.e17_musica_dk ) then
          phi0 = phiU + 1.613e-17_musica_dk*8.e17_musica_dk
        else
          phi0 = phiU + 1.613e-17_musica_dk*M
        endif
        quantum_yield = phiL + rONE/phi0
      endif
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch2chcho_quantum_yield_type
