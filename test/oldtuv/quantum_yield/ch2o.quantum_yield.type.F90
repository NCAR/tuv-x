! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch2o+hv->h2+co quantum yield module

!> The ch2o+hv->h2+co quantum yield type and related functions
module micm_ch2o_h2_co_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch2o_h2_co_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch2o+hv->h2+co quantum yield
  type, extends(base_quantum_yield_t) :: ch2o_h2_co_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ch2o_h2_co_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ch2o_h2_co_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter  :: lambdaL = 330._musica_dk
    real(musica_dk), parameter  :: lambdaU = 360._musica_dk
    character(len=*), parameter :: Iam = 'ch2o+hv->h2_co quantum yield calculate: '
    integer(musica_ik) :: m
    real(musica_dk) :: air_den_factor, Tfactor
    real(musica_dk), allocatable :: quantum_yield_chnl1(:)
    real(musica_dk), allocatable :: quantum_yield_chnl2(:)
    real(musica_dk), allocatable :: quantum_yield_tmp(:)
    real(musica_dk), allocatable :: quantum_yield_wrk(:)

    write(*,*) Iam,'entering'

    quantum_yield_chnl1 = this%quantum_yield(1)%array(:,1)
    quantum_yield_chnl2 = this%quantum_yield(1)%array(:,2)
    quantum_yield_tmp   = rONE - quantum_yield_chnl1
    allocate( quantum_yield_wrk(size(this%mdl_lambda_center)) )
    Tfactor = (300._musica_dk - environment%temperature)/80._musica_dk
    where( this%mdl_lambda_center >= lambdaL .and. this%mdl_lambda_center < lambdaU &
                                             .and. quantum_yield_chnl2 > rZERO )
      quantum_yield_wrk = (rONE - (quantum_yield_chnl1 + quantum_yield_chnl2)) &
                          /(2.45e19_musica_dk*quantum_yield_chnl2*quantum_yield_tmp)
      quantum_yield_wrk = quantum_yield_wrk*(rONE &
                        + .05_musica_dk*(this%mdl_lambda_center - 329._musica_dk)*Tfactor)
      quantum_yield = rONE/(rONE/quantum_yield_tmp + quantum_yield_wrk*environment%number_density_air)
    elsewhere
      quantum_yield = quantum_yield_chnl2
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch2o_h2_co_quantum_yield_type
