! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_radiator warehouse module

!> The radiator_warehouse_t type and related functions
!!
module micm_radiator_warehouse

  use musica_constants,                only : musica_dk, musica_ik
  use micm_abs_radiator_type,          only : radiator_ptr
  use musica_string,                   only : string_t

  implicit none

  private
  public :: radiator_warehouse_t

  !> radiator warehouse type
  type :: radiator_warehouse_t
    !> radiators
    type(radiator_ptr), allocatable :: radiators_(:)
    !> radiator "handles"
    type(string_t), allocatable     :: handle_(:)
  contains
    !> Get a copy of a specific radiator
    procedure :: get_radiator
    !> Finalize
    final     :: finalize
  end type radiator_warehouse_t

  !> radiator_warehouse_t constructor
  interface radiator_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Construct radiator_warehouse
  function constructor( config ) result( radiator_warehouse )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_assert,                 only : die_msg
    use micm_radiator_factory,         only : radiator_builder
    use musica_string,                 only : string_t

    !> arguments
    !> radiator configuration
    type(config_t), intent(inout) :: config
    !> warehouse type
    class(radiator_warehouse_t), pointer :: radiator_warehouse

    !> local variables
    character(len=*), parameter :: Iam = "Radiator warehouse constructor: "
    type(config_t)              :: radiator_config_set, radiator_config
    class(iterator_t), pointer  :: iter
    type(radiator_ptr)          :: aRadiator
    character(len=32)           :: keychar
    type(string_t)              :: keyString

    write(*,*) Iam // 'entering'

    call config%get( 'Radiators', radiator_config_set, Iam )
    iter => radiator_config_set%get_iterator()

    allocate( radiator_warehouse )
    allocate( radiator_warehouse%radiators_(0) )
    allocate( radiator_warehouse%handle_(0) )

    do while( iter%next() )
      keychar = radiator_config_set%key(iter)
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      keyString = keychar
      call radiator_config_set%get( iter, radiator_config, Iam )
      call radiator_config%add( 'Handle', keyString, Iam )
!-----------------------------------------------------------------------------
!> build and store the radiator
!-----------------------------------------------------------------------------
      aRadiator%val_ => radiator_builder( radiator_config )
      radiator_warehouse%radiators_ = [radiator_warehouse%radiators_,aRadiator]
      radiator_warehouse%handle_ = [radiator_warehouse%handle_,keyString]
    end do

    deallocate( iter )

    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' radiators'')') Iam,size(radiator_warehouse%radiators_)

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a specific radiator object
  function get_radiator( this, radiator_handle ) result( radiator_ptr )

    use micm_abs_radiator_type, only : abs_radiator_t
    use musica_string,      only : string_t
    use musica_constants,   only : lk => musica_lk, ik => musica_ik
    use musica_assert,      only : die_msg

    !> Arguments
    class(radiator_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)             :: radiator_handle

    class(abs_radiator_t), pointer         :: radiator_ptr

    !> Local variables
    character(len=*), parameter :: Iam = 'radiator warehouse get_radiator: '
    integer(ik) :: ndx
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do ndx = 1,size(this%handle_)
      if( radiator_handle .eq. this%handle_(ndx) ) then
        found = .true._lk
        exit
      endif
    end do

    if( found ) then
      allocate( radiator_ptr, source = this%radiators_(ndx)%val_ )
    else
      call die_msg( 460768324, "Invalid radiator handle: '"// radiator_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_radiator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radiator warehouse
  subroutine finalize( this )

    !> radiator warehouse
    type(radiator_warehouse_t), intent(inout) :: this

    integer(kind=musica_ik) :: ndx
    character(len=*), parameter :: Iam = 'radiator_warehouse finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%radiators_ ) ) then
      deallocate( this%radiators_ )
    endif

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_radiator_warehouse
