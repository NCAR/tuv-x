! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The spectral weight module

!> The abstract spectral wght type and related functions
module micm_abs_spectral_wght_type

  use musica_constants,                only : musica_dk, musica_ik
  use micm_environment,                only : environment_t

  implicit none
  private

  public :: abs_spectral_wght_t, abs_spectral_wght_ptr

  !> spectral wght abstract type
  type, abstract :: abs_spectral_wght_t
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    real(musica_dk), allocatable :: mdl_lambda_center(:)
  contains
    !> Calculate the spectral weight
    procedure(initial),   deferred :: initialize
    procedure(calculate), deferred :: calculate
    procedure                      :: addpnts
  end type abs_spectral_wght_t

  !> Pointer type for building sets of spectral wght objects
  type :: abs_spectral_wght_ptr
    class(abs_spectral_wght_t), pointer :: val_ => null( )
  end type abs_spectral_wght_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of conditions
  function calculate( this, environment ) result( spectral_wght )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk

    import abs_spectral_wght_t

    !> Spectral weight calculator
    class(abs_spectral_wght_t), intent(in) :: this
    !> spectral weight on model photolysis grid
    real(musica_dk)                   :: spectral_wght(size(this%mdl_lambda_center))
    !> Environmental conditions
    class(environment_t), intent(in)  :: environment
  end function calculate

  !> Initialize the base spectral weight type
  subroutine initial( this, config, mdlLambdaEdge )
    use musica_config,    only : config_t
    use musica_constants, only : musica_dk

    import abs_spectral_wght_t

    !> Spectral weight calculator
    class(abs_spectral_wght_t), intent(inout) :: this
    !> Environmental conditions
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)
    type(config_t), intent(inout) :: config
 end subroutine initial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end interface

contains

  subroutine addpnts( this, config, data_lambda, data_parameter )
    use musica_config, only : config_t
    use musica_string, only : string_t
    use photo_utils,   only : addpnt

    class(abs_spectral_wght_t)    :: this
    type(config_t), intent(inout) :: config
    real(musica_dk), allocatable, intent(inout) :: data_lambda(:)
    real(musica_dk), allocatable, intent(inout) :: data_parameter(:)

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    real(musica_dk), parameter :: large  = 1.e36_musica_dk
    character(len=*), parameter :: Iam = 'spectral_wght; addpnts: '

    integer(musica_ik) :: nRows
    real(musica_dk) :: lowerLambda, upperLambda
    real(musica_dk) :: addpnt_val_lower, addpnt_val_upper
    type(string_t)  :: addpnt_type_
    logical         :: found
    character(len=:), allocatable :: number

    write(*,*) Iam,'entering'

    !> add endpoints to data arrays; first the lower bound
    nRows = size(data_lambda)
    lowerLambda = data_lambda(1) ; upperLambda = data_lambda(nRows)
    call config%get( 'lower extrapolation', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      addpnt_val_lower = rZERO
    elseif( addpnt_type_ == 'boundary' ) then
      addpnt_val_lower = data_parameter(1)
    else
      number = addpnt_type_%to_char()
      read( number, '(g30.20)' ) addpnt_val_lower
    endif

    !> add endpoints to data arrays; now the upper bound
    call config%get( 'upper extrapolation', addpnt_type_, Iam, found=found )
    if( .not. found ) then
      addpnt_val_upper = rZERO
    elseif( addpnt_type_ == 'boundary' ) then
      addpnt_val_upper = data_parameter(nRows)
    else
      number = addpnt_type_%to_char()
      read( number, '(g30.20)' ) addpnt_val_upper
    endif

    call addpnt(x=data_lambda,y=data_parameter,xnew=-large,ynew=addpnt_val_lower) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE-deltax)*lowerLambda,ynew=addpnt_val_lower) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE+deltax)*upperLambda,ynew=addpnt_val_upper) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=large,ynew=addpnt_val_upper) 

    write(*,*) Iam,'exiting'

  end subroutine addpnts

end module micm_abs_spectral_wght_type
