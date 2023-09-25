! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The photolysis cross section module

!> The abstract cross section type and related functions
module micm_abs_cross_section_type

  use musica_constants,                only : musica_dk, musica_ik

  implicit none
  private

  public :: abs_cross_section_t, abs_cross_section_ptr

  !> Photo rate cross section abstract type
  type, abstract :: abs_cross_section_t
    real(musica_dk), allocatable :: mdl_lambda_edge(:)
    real(musica_dk), allocatable :: mdl_lambda_center(:)
  contains
    !> Calculate the photo rate cross section
    procedure(initial),   deferred :: initialize
    procedure(calculate), deferred :: calculate
    procedure                      :: addpnts
  end type abs_cross_section_t

  !> Pointer type for building sets of photo rate constants
  type :: abs_cross_section_ptr
    class(abs_cross_section_t), pointer :: val_ => null( )
  end type abs_cross_section_ptr

interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the rate constant for a given set of conditions
  function calculate( this, environment ) result( cross_section )
    use micm_environment,              only : environment_t
    use musica_constants,              only : musica_dk

    import abs_cross_section_t

    !> Cross section calculator
    class(abs_cross_section_t), intent(in) :: this
    !> cross section on model photo grid
    real(musica_dk)                   :: cross_section(size(this%mdl_lambda_center))
    !> Environmental conditions
    class(environment_t), intent(in)  :: environment
  end function calculate

  !> Initialize the base cross section type
  subroutine initial( this, config, mdlLambdaEdge )
    use musica_config,    only : config_t
    use musica_constants, only : musica_dk

    import abs_cross_section_t

    !> Cross section calculator
    class(abs_cross_section_t), intent(inout) :: this
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

    class(abs_cross_section_t)    :: this
    type(config_t), intent(inout) :: config
    real(musica_dk), allocatable, intent(inout) :: data_lambda(:)
    real(musica_dk), allocatable, intent(inout) :: data_parameter(:)

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    real(musica_dk), parameter :: deltax = 1.e-5_musica_dk
    character(len=*), parameter :: Iam = 'cross_section; addpnts: '

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

    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE-deltax)*lowerLambda,ynew=addpnt_val_lower) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=rZERO,ynew=addpnt_val_lower) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=(rONE+deltax)*upperLambda,ynew=addpnt_val_upper) 
    call addpnt(x=data_lambda,y=data_parameter,xnew=1.e38_musica_dk,ynew=addpnt_val_upper) 

    write(*,*) Iam,'exiting'

  end subroutine addpnts

end module micm_abs_cross_section_type
