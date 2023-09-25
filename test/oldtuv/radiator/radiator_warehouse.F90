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
  use musica_iterator,                 only : iterator_t

  implicit none

  private
  public :: radiator_warehouse_t, warehouse_iterator_t

  !> radiator warehouse type
  type :: radiator_warehouse_t
    private
    !> radiators
    type(radiator_ptr), allocatable :: radiators_(:)
    !> radiator "handles"
    type(string_t), allocatable     :: handle_(:)
  contains
    procedure, private :: get_radiator_from_handle
    procedure, public  :: get_radiator_ndx_from_handle
    procedure, private :: get_radiator_from_iterator
    !> Get a copy of a specific radiator
    generic   :: get_radiator => get_radiator_from_handle, get_radiator_from_iterator
    !> is radiator in warehouse?
    procedure :: in_warehouse
    !> Get radiator iterator
    procedure :: get_iterator
    !> Finalize
    final     :: finalize
  end type radiator_warehouse_t

  !>  radiator warehouse iterator
  type :: warehouse_iterator_t
    !> Pointer to the radiator warehouse
    type(radiator_warehouse_t), pointer :: warehouse_ => null( )
    !> Current index in the data set
    integer(kind=musica_ik) :: id_ = 0_musica_ik
  contains
    !> Advances to the next key-value pair
    procedure :: next => iterator_next
    !> Resets the iterator
    procedure :: reset => iterator_reset
  end type warehouse_iterator_t

  !> radiator_warehouse_t constructor
  interface radiator_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Construct radiator_warehouse
  function constructor( config, gridWareHouse ) result( radiator_warehouse )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_assert,                 only : die_msg
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_radiator_factory,         only : radiator_builder
    use musica_string,                 only : string_t

    !> arguments
    !> radiator configuration
    type(config_t), intent(inout)         :: config
    type(grid_warehouse_t), intent(inout) :: gridWareHouse
    !> warehouse type
    class(radiator_warehouse_t), pointer :: radiator_warehouse

    !> local variables
    character(len=*), parameter :: Iam = "Radiator warehouse constructor: "
    integer(musica_ik)          :: ndx
    type(config_t)              :: radiator_config_set, radiator_config
    class(iterator_t), pointer  :: iter
    type(radiator_ptr)          :: aRadiator
    character(len=32)           :: keychar
    type(string_t)              :: keyString

    write(*,*) Iam // 'entering'

    call config%get( 'Radiators', radiator_config_set, Iam )

    allocate( radiator_warehouse )
    allocate( radiator_warehouse%radiators_(0) )
    allocate( radiator_warehouse%handle_(0) )

    iter => radiator_config_set%get_iterator()
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
      aRadiator%val_ => radiator_builder( radiator_config, gridWareHouse )
      radiator_warehouse%radiators_ = [radiator_warehouse%radiators_,aRadiator]
      radiator_warehouse%handle_ = [radiator_warehouse%handle_,keyString]
    end do

    deallocate( iter )

    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' radiators'')') Iam,size(radiator_warehouse%radiators_)
    write(*,*) 'radiator handles'
    do ndx = 1,size(radiator_warehouse%handle_)
      write(*,'(a)') radiator_warehouse%handle_(ndx)%to_char()
    enddo

    write(*,*) Iam // 'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a specific radiator object
  function get_radiator_from_handle( this, radiator_handle ) result( radiator )

    use micm_abs_radiator_type, only : abs_radiator_t
    use musica_string,      only : string_t
    use musica_constants,   only : lk => musica_lk, ik => musica_ik
    use musica_assert,      only : die_msg

    !> Arguments
    class(radiator_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)                 :: radiator_handle
    class(abs_radiator_t), pointer             :: radiator

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
      radiator => this%radiators_(ndx)%val_
    else
      call die_msg( 460768324, "Invalid radiator handle: '"// radiator_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_radiator_from_handle

  !> Get index of a specific radiator object
  function get_radiator_ndx_from_handle( this, radiator_handle ) result( Index )

    use micm_abs_radiator_type, only : abs_radiator_t
    use musica_string,      only : string_t
    use musica_constants,   only : lk => musica_lk, ik => musica_ik
    use musica_assert,      only : die_msg

    !> Arguments
    class(radiator_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)                 :: radiator_handle
    integer(ik)                                :: Index


    !> Local variables
    character(len=*), parameter :: Iam = 'radiator warehouse get_radiator_ndx: '
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do Index = 1,size(this%handle_)
      if( radiator_handle .eq. this%handle_(Index) ) then
        found = .true._lk
        exit
      endif
    end do

    if( .not. found ) then
!     call die_msg( 460768324, "radiator '"// radiator_handle%to_char()//"' not found" )
      Index = -1
    endif

    write(*,*) Iam,'exiting'

  end function get_radiator_ndx_from_handle

  !> Get copy of a radiator object using an iterator
  function get_radiator_from_iterator( this, iterator ) result( radiator )

    use micm_abs_radiator_type, only : abs_radiator_t
    use musica_constants,   only : lk => musica_lk, ik => musica_ik
    use musica_assert,      only : die_msg

    !> Arguments
    class(radiator_warehouse_t), intent(inout) :: this
    type(warehouse_iterator_t),  intent(in)    :: iterator
    class(abs_radiator_t), pointer             :: radiator

    !> Local variables
    character(len=*), parameter :: Iam = 'radiator warehouse get_radiator from iterator: '
    integer(ik) :: ndx

    write(*,*) ' '
    write(*,*) Iam,'entering'

    ndx = iterator%id_
    write(*,*) Iam,'radiator handle = ',this%radiators_(ndx)%val_%handle_%to_char()

    radiator => this%radiators_(ndx)%val_

    write(*,*) Iam,'radiator diagnostics'
    write(*,*) Iam,'radiator handle = ',radiator%handle_%to_char()
    write(*,*) Iam,'radiator state OD is allocated = ',allocated(radiator%state_%layer_OD_)
    write(*,*) Iam,'radiator state SSA is allocated = ',allocated(radiator%state_%layer_SSA_)
    write(*,*) Iam,'radiator state G is allocated = ',allocated(radiator%state_%layer_G_)

    write(*,*) ' '
    write(*,*) Iam,'exiting'

!   stop 'debugging'

  end function get_radiator_from_iterator

  !> Is a radiator in the warehouse?
  function in_warehouse( this, radiator_handle )

    use micm_abs_radiator_type, only : abs_radiator_t
    use musica_string,      only : string_t
    use musica_constants,   only : lk => musica_lk, ik => musica_ik
    use musica_assert,      only : die_msg

    !> Arguments
    class(radiator_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)             :: radiator_handle

    logical(lk)                            :: in_warehouse

    !> Local variables
    character(len=*), parameter :: Iam = 'radiator in warehouse: '
    integer(ik) :: ndx

    write(*,*) ' '
    write(*,*) Iam,'entering'

    in_warehouse = .false._lk
    do ndx = 1,size(this%handle_)
      if( radiator_handle == this%handle_(ndx) ) then
        in_warehouse = .true._lk
        exit
      endif
    end do

    write(*,*) Iam,'exiting'

  end function in_warehouse

  !> Gets an interator for the radiator warehouse
  function get_iterator( this )

    use musica_assert,                 only : assert

    !> Pointer to the iterator
    class(warehouse_iterator_t), pointer :: get_iterator
    !> Radiator warehouse
    class(radiator_warehouse_t), intent(in), target :: this

    call assert( 753334333, allocated( this%radiators_ ) )
    allocate( warehouse_iterator_t :: get_iterator )
    select type( iter => get_iterator )
      type is( warehouse_iterator_t )
        iter%warehouse_ => this
    end select

  end function get_iterator

  !> Advances the iterator
  !> Returns false if the end of the collection has been reached
  function iterator_next( this ) result( Continue )

    use musica_constants, only : lk => musica_lk, ik => musica_ik

    !> Iterator
    class(warehouse_iterator_t), intent(inout) :: this

    logical(lk) :: Continue

    this%id_ = this%id_ + 1_ik
    Continue =  this%id_ <= size(this%warehouse_%radiators_)

  end function iterator_next

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Resets the iterator
  subroutine iterator_reset( this )

    use musica_constants, only : ik => musica_ik

    !> Iterator
    class(warehouse_iterator_t), intent(inout) :: this

    this%id_ = 0_ik

  end subroutine iterator_reset

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
