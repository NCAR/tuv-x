! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_photo_kinetics module

!> The photo_kinetics_t type and related functions
!!
!! \todo Consider rewriting micm_kinetics to follow musica style conventions
module micm_photo_kinetics

  use micm_abs_cross_section_type,     only : abs_cross_section_ptr
  use micm_abs_quantum_yield_type,     only : abs_quantum_yield_ptr
  use musica_constants,                only : musica_dk, musica_ik
  use musica_string,                   only : string_t
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_constants,              only : musica_rk
    use musica_assert,                 only : die_msg
    use micm_cross_section_factory,    only : cross_section_builder
    use micm_quantum_yield_factory,    only : quantum_yield_builder

  implicit none

  private
  public :: photo_kinetics_t

  !> Kinetics
  type :: photo_kinetics_t
    !> Photo rate constant calculators
    type(abs_cross_section_ptr), allocatable :: cross_section_objs_(:)
    type(abs_quantum_yield_ptr), allocatable :: quantum_yield_objs_(:)
    !> Current photo rate constant values
    real(kind=musica_dk), allocatable :: cross_section_values_(:,:)
    real(kind=musica_dk), allocatable :: quantum_yield_values_(:,:)
    real(kind=musica_dk), allocatable :: rate_constant_alias_factor_(:)
    type(string_t), allocatable       :: reaction_key(:)
  contains
    !> Update the object for new environmental conditions
    procedure :: update_for_new_environmental_state
    !> Finalize the object
    final :: finalize
  end type photo_kinetics_t

  !> photo_kinetics_t constructor
  interface photo_kinetics_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of photo_kinetics_t objects
  function constructor( config,mdlLambdaEdge ) result( new_photo_kinetics_obj )


    real(musica_dk), intent(in)      :: mdlLambdaEdge(:)
    !> Kinetics configuration data
    type(config_t), intent(inout) :: config
    !> New kinetics
    class(photo_kinetics_t), pointer :: new_photo_kinetics_obj

    !> local variables
    integer(musica_ik), parameter :: noErr = 0

    real(musica_dk) :: rate_aliasing_factor
    integer :: nSize, ndx
    character(len=*), parameter :: Iam = "Photo kinetics constructor: "
    type(config_t) :: reaction_set, reaction_config
    type(config_t) :: cross_section_config, quantum_yield_config
    class(iterator_t), pointer :: iter
    type(string_t) :: netcdfFile, Object
    type(abs_cross_section_ptr) :: cross_section_ptr
    type(abs_quantum_yield_ptr) :: quantum_yield_ptr
    character(:), allocatable   :: jsonkey
    character(len=32)           :: keychar
    type(string_t)              :: areaction_key
    type(string_t), allocatable :: netcdfFiles(:)
    logical :: found

    allocate( new_photo_kinetics_obj )

    associate(new_obj=>new_photo_kinetics_obj)

    allocate( string_t :: new_obj%reaction_key(0) )

    allocate( new_obj%cross_section_objs_(0) )
    allocate( new_obj%quantum_yield_objs_(0) )
    allocate( new_obj%rate_constant_alias_factor_(0) )

    jsonkey = 'photo dissociation rate constants'
    call config%get( jsonkey, reaction_set, Iam )
    iter => reaction_set%get_iterator( )
!-----------------------------------------------------------------------------
!> iterate over photo reactions
!-----------------------------------------------------------------------------
    do while( iter%next( ) )
      keychar = reaction_set%key(iter)
      areaction_key = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      new_obj%reaction_key = [new_obj%reaction_key,areaction_key]
      call reaction_set%get( iter, reaction_config, Iam )
!-----------------------------------------------------------------------------
!> cross section first
!-----------------------------------------------------------------------------
      call reaction_config%get( "cross section", cross_section_config, Iam )
      cross_section_ptr%val_ => cross_section_builder( cross_section_config,mdlLambdaEdge )
      new_obj%cross_section_objs_ = [new_obj%cross_section_objs_,cross_section_ptr]
!-----------------------------------------------------------------------------
!> now quantum yield
!-----------------------------------------------------------------------------
      call reaction_config%get( "quantum yield", quantum_yield_config, Iam )
      quantum_yield_ptr%val_ => quantum_yield_builder( quantum_yield_config,mdlLambdaEdge )
      new_obj%quantum_yield_objs_ = [new_obj%quantum_yield_objs_,quantum_yield_ptr]
!-----------------------------------------------------------------------------
!> finally "aliasing" factor
!-----------------------------------------------------------------------------
      call reaction_config%get( "rate constant alias factor", rate_aliasing_factor, Iam, found=found )
      if( .not. found ) then
        rate_aliasing_factor = 1.0_musica_dk
      endif
      new_obj%rate_constant_alias_factor_ = [new_obj%rate_constant_alias_factor_,rate_aliasing_factor]
    end do

    deallocate( iter )

    nSize = size(new_obj%cross_section_objs_)
    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' cross sections'')') Iam,nSize
    nSize = size(new_obj%quantum_yield_objs_)
    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' quantum yields'')') Iam,nSize

!-----------------------------------------------------------------------------
!> setup cross section, quantum yield arrays
!-----------------------------------------------------------------------------
    allocate( new_obj%cross_section_values_(0,0) )
    allocate( new_obj%quantum_yield_values_(0,0) )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the object for new environmental conditions
  subroutine update_for_new_environmental_state( this, environment, nwave )

    use micm_environment, only : environment_t
    use musica_assert,    only : die_msg

    !> Kinetics
    class(photo_kinetics_t), intent(inout) :: this
    integer(musica_ik), intent(in)         :: nwave
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'update_for_new_environmental_state: '
    integer(kind=musica_ik) :: ndx
    real(musica_dk), allocatable :: a_cross_section(:)
    real(musica_dk), allocatable :: cross_section_tray(:)
    real(musica_dk), allocatable :: a_quantum_yield(:)
    real(musica_dk), allocatable :: quantum_yield_tray(:)

    write(*,*) ' '
    write(*,*) Iam,'entering'

    allocate(cross_section_tray(0))
    do ndx = 1, size(this%cross_section_objs_)
      associate( calc_ftn => this%cross_section_objs_(ndx)%val_ )
        a_cross_section = calc_ftn%calculate( environment )
      end associate
      cross_section_tray = [cross_section_tray,a_cross_section]
    end do

    this%cross_section_values_ = reshape( cross_section_tray, &
                                          (/nwave,size(this%cross_section_objs_) /) )

    write(*,*) Iam,'size of cross section values = ',&
        size(this%cross_section_values_,dim=1), size(this%cross_section_values_,dim=2)

    allocate(quantum_yield_tray(0))
    do ndx = 1, size(this%quantum_yield_objs_)
      associate( calc_ftn => this%quantum_yield_objs_(ndx) )
        if( .not. associated(calc_ftn%val_) ) then
          call die_msg( 200000009,Iam//'quantum yield is NOT associated!!!' )
        endif
        a_quantum_yield = calc_ftn%val_%calculate( environment )
      end associate
      quantum_yield_tray = [quantum_yield_tray,a_quantum_yield]
    end do

    this%quantum_yield_values_ = reshape( quantum_yield_tray, &
                                          (/nwave,size(this%quantum_yield_objs_) /) )

    write(*,*) Iam,'size of quantum_yield values = ',&
        size(this%quantum_yield_values_,dim=1), size(this%quantum_yield_values_,dim=2)

    write(*,*) Iam,'exiting'

  end subroutine update_for_new_environmental_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the kinetics
  subroutine finalize( this )

    !> Kinetics
    type(photo_kinetics_t), intent(inout) :: this

    integer(kind=musica_ik) :: ndx
    character(len=*), parameter :: Iam = 'photo_kinetics finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%cross_section_values_ ) ) then
      deallocate( this%cross_section_values_ )
    endif
    if( allocated( this%quantum_yield_values_ ) ) then
      deallocate( this%quantum_yield_values_ )
    endif
    if( allocated( this%cross_section_objs_ ) ) then
      do ndx = 1,size(this%cross_section_objs_)
        if( associated( this%cross_section_objs_(ndx)%val_ ) ) then
          deallocate( this%cross_section_objs_(ndx)%val_ )
        endif
      enddo
      deallocate( this%cross_section_objs_ )
    end if

    if( allocated( this%quantum_yield_objs_ ) ) then
      do ndx = 1,size(this%quantum_yield_objs_)
        if( associated( this%quantum_yield_objs_(ndx)%val_ ) ) then
          deallocate( this%quantum_yield_objs_(ndx)%val_ )
        endif
      enddo
      deallocate( this%quantum_yield_objs_ )
    end if
    if( allocated( this%rate_constant_alias_factor_ ) ) then
      deallocate( this%rate_constant_alias_factor_ )
    end if
    if( allocated( this%reaction_key ) ) then
      deallocate( this%reaction_key )
    end if

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_photo_kinetics
