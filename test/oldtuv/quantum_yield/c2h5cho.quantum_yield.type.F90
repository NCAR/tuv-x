! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This c2h5cho+hv->c2h5+hco quantum yield module

!> The c2h5cho+hv->c2h5+hco quantum yield type and related functions
module micm_c2h5cho_c2h5_hco_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: c2h5cho_c2h5_hco_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for c2h5cho+hv->c2h5+hco quantum yield
  type, extends(base_quantum_yield_t) :: c2h5cho_c2h5_hco_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type c2h5cho_c2h5_hco_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(c2h5cho_c2h5_hco_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: largest=1.E+36_musica_dk
    real(musica_dk), parameter :: pzero = 10._musica_dk/largest
    character(len=*), parameter :: Iam = 'c2h5cho+hv->c2h5+hco calculate: '
    real(musica_dk) :: air_den_fac

    write(*,*) Iam,'entering'

    air_den_fac = environment%number_density_air/2.45e19_musica_dk

! quantum yields:
! use Stern-Volmer pressure dependence:
    where( this%quantum_yield(1)%array(:,1) < pzero )
      quantum_yield = rZERO
    elsewhere
      quantum_yield = rONE/(rONE + (rONE/this%quantum_yield(1)%array(:,1) - rONE)*air_den_fac)
      quantum_yield = min( rONE,quantum_yield )
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_c2h5cho_c2h5_hco_quantum_yield_type
