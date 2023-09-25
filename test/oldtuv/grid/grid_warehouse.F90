! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_grid_warehouse module

!> The grid warehouse type and related functions
!!
module micm_grid_warehouse

  use micm_1d_grid, only : base_grid_ptr
! use micm_1d_grid, only : base_grid_t

  implicit none

  private
  public :: grid_warehouse_t

  !> Grid warehouse type
  type :: grid_warehouse_t
    private
    !> grid objects
!   class(base_grid_ptr), allocatable :: grid_objs_(:)
    type(base_grid_ptr), allocatable :: grid_objs_(:)
  contains
    !> get a copy of a grid object
    procedure :: get_grid
    !> Finalize the object
    final :: finalize
  end type grid_warehouse_t

  !> Grid warehouse_t constructor
  interface grid_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Grid warehouse constructor
  function constructor( config ) result( grid_warehouse_obj )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_string,                 only : string_t
    use micm_grid_factory,             only : grid_builder

    !> Arguments
    !> grid configuration data
    type(config_t), intent(inout) :: config

    !> New grid_warehouse_obj
    class(grid_warehouse_t), pointer :: grid_warehouse_obj

    !> local variables
    character(len=*), parameter :: Iam = "Grid warehouse constructor: "
    type(config_t)              :: grid_set, grid_config
    class(iterator_t), pointer  :: iter
    class(grid_warehouse_t), pointer :: grid_warehouse_ptr
    type(base_grid_ptr)            :: grid_obj
    character(len=32)           :: keychar
    type(string_t)              :: aswkey

    write(*,*) Iam // 'entering'

    allocate( grid_warehouse_obj )

    associate(new_obj=>grid_warehouse_obj)

    allocate( new_obj%grid_objs_(0) )

    call config%get( 'Grids', grid_set, Iam )
    iter => grid_set%get_iterator()
!-----------------------------------------------------------------------------
!> iterate over grids
!-----------------------------------------------------------------------------
    do while( iter%next() )
      keychar = grid_set%key(iter)
      aswkey  = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      call grid_set%get( iter, grid_config, Iam )
      call grid_config%add( 'Handle', aswkey, Iam )
!-----------------------------------------------------------------------------
!> Build grid objects
!-----------------------------------------------------------------------------
      grid_obj%ptr_ => grid_builder( grid_config )
      new_obj%grid_objs_ = [new_obj%grid_objs_,grid_obj]
    end do

    deallocate( iter )

    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' grid objects'')') Iam,size(new_obj%grid_objs_)

    end associate

    write(*,*) Iam // 'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a grid object
  function get_grid( this, grid_handle ) result( grid_ptr )

    use micm_1d_grid,      only : base_grid_t
    use musica_string,     only : string_t
    use musica_constants,  only : lk => musica_lk, ik => musica_ik
    use musica_assert,     only : die_msg

    !> Arguments
    class(grid_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)             :: grid_handle

    class(base_grid_t), pointer          :: grid_ptr

    !> Local variables
    character(len=*), parameter :: Iam = 'grid warehouse get_grid: '
    integer(ik) :: ndx
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do ndx = 1,size(this%grid_objs_)
      if( grid_handle .eq. this%grid_objs_(ndx)%ptr_%handle_ ) then
        found = .true._lk
        exit
      endif
    end do

    if( found ) then
      allocate( grid_ptr, source = this%grid_objs_(ndx)%ptr_ )
    else
      call die_msg( 460768214, "Invalid grid handle: '"// grid_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize grid warehouse
  subroutine finalize( this )

    use musica_constants, only : ik => musica_ik

    !> Arguments
    type(grid_warehouse_t), intent(inout) :: this

    !> Local variables
    integer(kind=ik) :: ndx
    character(len=*), parameter :: Iam = 'grid warehouse finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%grid_objs_ ) ) then
      deallocate( this%grid_objs_ )
    endif

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_grid_warehouse
