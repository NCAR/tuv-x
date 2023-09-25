! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_vert_Profile_warehouse module

!> The vert Profile warehouse type and related functions
!!
module micm_vert_Profile_warehouse

  use micm_vert_Profile, only : abs_vert_Profile_ptr

  implicit none

  private
  public :: vert_Profile_warehouse_t

  !> Vert Profile warehouse type
  type :: vert_Profile_warehouse_t
    !> Vert Profile objects
    type(abs_vert_Profile_ptr), allocatable :: vert_Profile_objs_(:)
  contains
    !> get a copy of a vert Profile object
    procedure :: get_vert_Profile
    !> Finalize the object
    final :: finalize
  end type vert_Profile_warehouse_t

  !> Grid warehouse_t constructor
  interface vert_Profile_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Vert Profile warehouse constructor
  function constructor( config, gridwarehouse ) result( vert_Profile_warehouse_obj )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_vert_Profile_factory,     only : vert_Profile_builder

    !> Arguments
    !> Vert Profile configuration data
    type(config_t), intent(inout) :: config
    type(grid_warehouse_t), intent(inout) :: gridwarehouse

    !> New vert_Profile_warehouse_obj
    class(vert_Profile_warehouse_t), pointer :: vert_Profile_warehouse_obj

    !> local variables
    character(len=*), parameter :: Iam = "Vert Profile warehouse constructor: "
    type(config_t)              :: vert_Profile_set, vert_Profile_config
    class(iterator_t), pointer  :: iter
    class(vert_Profile_warehouse_t), pointer :: vert_Profile_warehouse_ptr
    type(abs_vert_Profile_ptr)            :: vert_Profile_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey

    write(*,*) Iam // 'entering'

    allocate( vert_Profile_warehouse_obj )

    associate(new_obj=>vert_Profile_warehouse_obj)

    allocate( new_obj%vert_Profile_objs_(0) )

    call config%get( 'Vertical Profiles', vert_Profile_set, Iam )
    iter => vert_Profile_set%get_iterator()
!-----------------------------------------------------------------------------
!> iterate over vert Profiles
!-----------------------------------------------------------------------------
    do while( iter%next() )
      keychar = vert_Profile_set%key(iter)
      aswkey  = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      call vert_Profile_set%get( iter, vert_Profile_config, Iam )
      call vert_Profile_config%add( 'Handle', aswkey, Iam )
!-----------------------------------------------------------------------------
!> Build vert Profile objects
!-----------------------------------------------------------------------------
      vert_Profile_obj%ptr_ => vert_Profile_builder( vert_Profile_config, gridwarehouse )
      new_obj%vert_Profile_objs_ = [new_obj%vert_Profile_objs_,vert_Profile_obj]
    end do

    deallocate( iter )

    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' vert Profile objects'')') Iam,size(new_obj%vert_Profile_objs_)

    end associate

    write(*,*) Iam // 'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a vert Profile object
  function get_vert_Profile( this, vert_Profile_handle ) result( vert_Profile_ptr )

    use musica_string,     only : string_t
    use musica_constants,  only : lk => musica_lk, ik => musica_ik
    use musica_assert,     only : die_msg
    use micm_vert_Profile, only : abs_vert_Profile_t

    !> Arguments
    class(vert_Profile_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)                     :: vert_Profile_handle

    class(abs_vert_Profile_t), pointer          :: vert_Profile_ptr

    !> Local variables
    character(len=*), parameter :: Iam = 'vert Profile warehouse get_vert_Profile: '
    integer(ik) :: ndx
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do ndx = 1,size(this%vert_Profile_objs_)
      if( vert_Profile_handle .eq. this%vert_Profile_objs_(ndx)%ptr_%handle_ ) then
        found = .true._lk
        exit
      endif
    end do

    if( found ) then
      allocate( vert_Profile_ptr, source = this%vert_Profile_objs_(ndx)%ptr_ )
    else
      call die_msg( 460768214, "Invalid vert Profile handle: '"// vert_Profile_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_vert_Profile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize vert Profile warehouse
  subroutine finalize( this )

    use musica_constants, only : ik => musica_ik

    !> Arguments
    type(vert_Profile_warehouse_t), intent(inout) :: this

    !> Local variables
    integer(kind=ik) :: ndx
    character(len=*), parameter :: Iam = 'vert Profile warehouse finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%vert_Profile_objs_ ) ) then
      deallocate( this%vert_Profile_objs_ )
    endif

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_vert_Profile_warehouse
