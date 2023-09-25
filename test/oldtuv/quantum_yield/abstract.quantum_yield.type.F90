! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The abstract photolysis quantum yield module

!> The abstract quantum yield type and related functions
module micm_abs_quantum_yield_type

  use musica_constants,                only : musica_dk, musica_ik
  use micm_environment,                only : environment_t

  implicit none
  private

  public :: abs_quantum_yield_t, abs_quantum_yield_ptr

  !> Photo rate quantum yield abstract type
  type, abstract :: abs_quantum_yield_t
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    real(musica_dk), allocatable :: mdl_lambda_center(:)
  contains
    procedure(initial),   deferred :: initialize
    !> Calculate the photo rate quantum yield
    procedure(calculate), deferred :: calculate
    procedure                      :: addpnts
  end type abs_quantum_yield_t

  !> Pointer type for building sets of photo rate constants
  type :: abs_quantum_yield_ptr
    class(abs_quantum_yield_t), pointer :: val_ => null( )
  end type abs_quantum_yield_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of conditions
  function calculate( this, environment ) result( quantum_yield )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk
!   use micm_photolysis_wavelength_grid,  only : wavelength_grid
    import abs_quantum_yield_t

    !> Quantum yield calculator
    class(abs_quantum_yield_t), intent(in) :: this
    !> quantum yield on model photo grid
    real(kind=musica_dk)              :: quantum_yield(size(this%mdl_lambda_center))
    !> Environmental conditions
    class(environment_t), intent(in) :: environment
  end function calculate

  !> Initialize the base quantum yield type
  subroutine initial( this, config, mdlLambdaEdge )
    use musica_config,    only : config_t
    use musica_constants, only : musica_dk

    import abs_quantum_yield_t

    !> Quantum yield calculator
    class(abs_quantum_yield_t), intent(inout) :: this
    !> Environmental conditions
    type(config_t), intent(inout) :: config
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)
 end subroutine initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

  subroutine addpnts( this, config, data_lambda, data_parameter )
    use musica_config, only : config_t
    use musica_string, only : string_t
    use photo_utils,   only : addpnt

    class(abs_quantum_yield_t)    :: this
    type(config_t), intent(inout) :: config
    real(musica_dk), allocatable, intent(inout) :: data_lambda(:)
    real(musica_dk), allocatable, intent(inout) :: data_parameter(:)

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    character(len=*), parameter :: Iam = 'cross_section; addpnts: '

    integer(musica_ik) :: nRows
    real(musica_dk) :: lowerLambda, upperLambda
    real(musica_dk) :: addpnt_val_
    type(string_t)  :: addpnt_type_
    logical         :: found
    character(len=:), allocatable :: number

    write(*,*) Iam,'entering'

    !> add endpoints to data arrays; first the lower bound
    nRows = size(data_lambda)
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda(nRows)
    call config%get( 'lower extrapolation', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      addpnt_val_ = rZERO
    elseif( addpnt_type_ == 'boundary' ) then
      addpnt_val_ = data_parameter(1)
    else
      number = addpnt_type_%to_char()
      read( number, '(g30.20)' ) addpnt_val_
    endif

    call addpnt(x=data_lambda,y=data_parameter,xnew=rZERO,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE-deltax)*lowerLambda,ynew=addpnt_val_) 
    !> add endpoints to data arrays; now the upper bound
    call config%get( 'upper extrapolation', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      addpnt_val_ = rZERO
    elseif( addpnt_type_ == 'boundary' ) then
      addpnt_val_ = data_parameter(nRows)
    else
      number = addpnt_type_%to_char()
      read( number, '(g30.20)' ) addpnt_val_
    endif

    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE+deltax)*upperLambda,ynew=addpnt_val_) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=1.e38_musica_dk,ynew=addpnt_val_) 

    write(*,*) Iam,'exiting'

  end subroutine addpnts

end module micm_abs_quantum_yield_type
