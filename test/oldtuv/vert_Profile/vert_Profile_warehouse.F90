! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_Profile_warehouse module

!> The Profile warehouse type and related functions
!!
module micm_Profile_warehouse

  use micm_Profile, only : abs_Profile_ptr

  implicit none

  private
  public :: Profile_warehouse_t

  !> Vert Profile warehouse type
  type :: Profile_warehouse_t
    private
    !> Vert Profile objects
    type(abs_Profile_ptr), allocatable :: Profile_objs_(:)
  contains
    !> get a copy of a Profile object
    procedure :: get_Profile
    !> Finalize the object
    final :: finalize
  end type Profile_warehouse_t

  !> Grid warehouse_t constructor
  interface Profile_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Vert Profile warehouse constructor
  function constructor( config, gridwarehouse ) result( Profile_warehouse_obj )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_Profile_factory,     only : Profile_builder

    !> Arguments
    !> Vert Profile configuration data
    type(config_t), intent(inout) :: config
    type(grid_warehouse_t), intent(inout) :: gridwarehouse

    !> New Profile_warehouse_obj
    class(Profile_warehouse_t), pointer :: Profile_warehouse_obj

    !> local variables
    character(len=*), parameter :: Iam = "Vert Profile warehouse constructor: "
    type(config_t)              :: Profile_set, Profile_config
    class(iterator_t), pointer  :: iter
    class(Profile_warehouse_t), pointer :: Profile_warehouse_ptr
    type(abs_Profile_ptr)            :: Profile_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey

    write(*,*) Iam // 'entering'

    allocate( Profile_warehouse_obj )

    associate(new_obj=>Profile_warehouse_obj)

    allocate( new_obj%Profile_objs_(0) )

    call config%get( 'Vertical Profiles', Profile_set, Iam )
    iter => Profile_set%get_iterator()
!-----------------------------------------------------------------------------
!> iterate over Profiles
!-----------------------------------------------------------------------------
    do while( iter%next() )
      keychar = Profile_set%key(iter)
      aswkey  = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      call Profile_set%get( iter, Profile_config, Iam )
      call Profile_config%add( 'Handle', aswkey, Iam )
!-----------------------------------------------------------------------------
!> Build Profile objects
!-----------------------------------------------------------------------------
      Profile_obj%ptr_ => Profile_builder( Profile_config, gridwarehouse )
      new_obj%Profile_objs_ = [new_obj%Profile_objs_,Profile_obj]
    end do

    deallocate( iter )

    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' Profile objects'')') Iam,size(new_obj%Profile_objs_)

    end associate

    write(*,*) Iam // 'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a Profile object
  function get_Profile( this, Profile_handle ) result( Profile_ptr )

    use musica_string,     only : string_t
    use musica_constants,  only : lk => musica_lk, ik => musica_ik
    use musica_assert,     only : die_msg
    use micm_Profile, only : base_profile_t

    !> Arguments
    class(Profile_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)                     :: Profile_handle

    class(base_profile_t), pointer          :: Profile_ptr

    !> Local variables
    character(len=*), parameter :: Iam = 'Profile warehouse get_Profile: '
    integer(ik) :: ndx
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do ndx = 1,size(this%Profile_objs_)
      if( Profile_handle .eq. this%Profile_objs_(ndx)%ptr_%handle_ ) then
        found = .true._lk
        exit
      endif
    end do

    if( found ) then
      allocate( Profile_ptr, source = this%Profile_objs_(ndx)%ptr_ )
    else
      call die_msg( 460768214, "Invalid Profile handle: '"// Profile_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_Profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize Profile warehouse
  subroutine finalize( this )

    use musica_constants, only : ik => musica_ik

    !> Arguments
    type(Profile_warehouse_t), intent(inout) :: this

    !> Local variables
    integer(kind=ik) :: ndx
    character(len=*), parameter :: Iam = 'Profile warehouse finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%Profile_objs_ ) ) then
      deallocate( this%Profile_objs_ )
    endif

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_Profile_warehouse
