! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3cho+hv->ch3+hco quantum yield module

!> The ch3cho+hv->ch3+hco quantum yield type and related functions
module micm_ch3cho_ch3_hco_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3cho_ch3_hco_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch3cho+hv->ch3+hco quantum yield
  type, extends(base_quantum_yield_t) :: ch3cho_ch3_hco_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ch3cho_ch3_hco_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ch3cho_ch3_hco_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'ch3cho+hv->ch3_hco quantum yield calculate: '
    integer(musica_ik) :: m
    real(musica_dk) :: air_den_factor
    real(musica_dk), allocatable :: quantum_yield_chnl1(:)
    real(musica_dk), allocatable :: quantum_yield_chnl2(:)
    real(musica_dk), allocatable :: quantum_yield_wrk(:)

    write(*,*) Iam,'entering'

    quantum_yield_chnl1 = this%quantum_yield(1)%array(:,2)
    quantum_yield_chnl2 = rONE - this%quantum_yield(1)%array(:,1)
    quantum_yield_wrk = (/ (rZERO,m=1,size(this%quantum_yield(1)%array,dim=1)) /)
    where( quantum_yield_chnl1(:) > rZERO )
      quantum_yield_wrk = quantum_yield_chnl2/quantum_yield_chnl1 - rONE
    elsewhere
      quantum_yield_wrk = rZERO
    endwhere
    air_den_factor = environment%number_density_air/2.465e19_musica_dk
    quantum_yield = quantum_yield_chnl1*(rONE + quantum_yield_wrk) &
                    /(rONE + quantum_yield_wrk*air_den_factor)
    quantum_yield = min( rONE,max(rZERO,quantum_yield) )

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3cho_ch3_hco_quantum_yield_type
