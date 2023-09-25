! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_radXfer_xsect_warehouse module

!> The radXfer_xsect_warehouse_t type and related functions
!!
module micm_radXfer_xsect_warehouse

  use micm_radXfer_abs_cross_section_type,     only : abs_cross_section_ptr
  use micm_radXfer_cross_section_factory, only : cross_section_builder
  use musica_constants,                only : musica_dk, musica_ik
  use musica_string,                   only : string_t

  implicit none

  private
  public :: radXfer_xsect_warehouse_t

  !> Radiative xfer cross section type
  type :: radXfer_xsect_warehouse_t
    private
    !> cross section calculators
    type(abs_cross_section_ptr), allocatable :: cross_section_objs_(:)
    !> cross section "handle"
    type(string_t), allocatable              :: handles_(:)
  contains
    !> Get a copy of a specific radXfer cross section
    procedure :: get_radXfer_cross_section
    !> Finalize the object
    final :: finalize
  end type radXfer_xsect_warehouse_t

  !> radXfer_xsect_warehouse_t constructor
  interface radXfer_xsect_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of radXfer_xsect_warehouse_t objects
  function constructor( config, gridWareHouse, ProfileWareHouse ) result( radXfer_xsect_obj )

    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_constants,              only : musica_rk, musica_lk
    use musica_assert,                 only : die_msg
    use micm_grid_warehouse,           only : grid_warehouse_t
    use micm_Profile_warehouse,        only : Profile_warehouse_t

    !> Arguments
    !> radXfer configuration data
    type(config_t), intent(inout) :: config
    !> The warehouses
    type(grid_warehouse_t), intent(inout)         :: gridWareHouse
    type(Profile_warehouse_t), intent(inout) :: ProfileWareHouse

    !> New radiative xfer cross section obj
    class(radXfer_xsect_warehouse_t), pointer :: radXfer_xsect_obj

    !> local variables
    integer :: nSize, ndx
    character(len=*), parameter :: Iam = "Radiative xfer xsect constructor: "
    type(config_t) :: reaction_set, reaction_config
    type(config_t) :: cross_section_config
    class(iterator_t), pointer :: iter
    type(abs_cross_section_ptr) :: cross_section_ptr
    character(:), allocatable   :: jsonkey
    character(len=32)           :: keychar
    type(string_t)              :: areaction_key
    type(string_t), allocatable :: netcdfFiles(:)
    logical(musica_lk)          :: found

    write(*,*) ' '
    write(*,*) Iam // 'entering'

    allocate( radXfer_xsect_obj )

    associate(new_obj=>radXfer_xsect_obj)

    allocate( string_t :: new_obj%handles_(0) )

    allocate( new_obj%cross_section_objs_(0) )

    jsonkey = 'Radiative xfer cross sections'
    call config%get( jsonkey, reaction_set, Iam, found=found )
    !> radXfer cross section objects not required
has_radXfer_xsects: &
    if( found ) then
      iter => reaction_set%get_iterator( )
!-----------------------------------------------------------------------------
!> iterate over cross sections
!-----------------------------------------------------------------------------
      do while( iter%next( ) )
        keychar = reaction_set%key(iter)
        areaction_key = keychar 
        write(*,*) ' '
        write(*,*) Iam,'key = ',trim(keychar)
        new_obj%handles_ = [new_obj%handles_,areaction_key]
        call reaction_set%get( iter, reaction_config, Iam )
!-----------------------------------------------------------------------------
!> build and store cross section object
!-----------------------------------------------------------------------------
        call reaction_config%get( "cross section", cross_section_config, Iam )
        cross_section_ptr%val_ => cross_section_builder( cross_section_config, gridWareHouse, ProfileWareHouse )
        new_obj%cross_section_objs_ = [new_obj%cross_section_objs_,cross_section_ptr]
      end do
      deallocate( iter )
    endif has_radXfer_xsects

    nSize = size(new_obj%cross_section_objs_)
    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' cross sections'')') Iam,nSize

    end associate

    write(*,*) ' '
    write(*,*) Iam // 'exiting'

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Get copy of a specific radXfer cross section object
  function get_radXfer_cross_section( this, radXfer_cross_section_handle ) result( radXfer_cross_section_ptr )

    use micm_radXfer_abs_cross_section_type, only : abs_cross_section_t
    use musica_string,     only : string_t
    use musica_constants,  only : lk => musica_lk, ik => musica_ik
    use musica_assert,     only : die_msg

    !> Arguments
    class(radXfer_xsect_warehouse_t), intent(inout) :: this
    type(string_t), intent(in)             :: radXfer_cross_section_handle

    class(abs_cross_section_t), pointer    :: radXfer_cross_section_ptr

    !> Local variables
    character(len=*), parameter :: Iam = 'radXfer cross section warehouse get_radXfer_cross_section: '
    integer(ik) :: ndx
    logical(lk) :: found

    write(*,*) ' '
    write(*,*) Iam,'entering'

    found = .false._lk
    do ndx = 1,size(this%handles_)
      if( radXfer_cross_section_handle .eq. this%handles_(ndx) ) then
        found = .true._lk
        exit
      endif
    end do

    if( found ) then
      allocate( radXfer_cross_section_ptr, source = this%cross_section_objs_(ndx)%val_ )
    else
      call die_msg( 460768224, "Invalid radXfer_cross_section_handle: '"// radXfer_cross_section_handle%to_char()//"'" )
    endif

    write(*,*) Iam,'exiting'

  end function get_radXfer_cross_section

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the radXfer xsect warehouse
  subroutine finalize( this )

    !> radXfer xsect warehouse
    type(radXfer_xsect_warehouse_t), intent(inout) :: this

    integer(kind=musica_ik) :: ndx
    character(len=*), parameter :: Iam = 'radXfer_xsect finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%cross_section_objs_ ) ) then
      do ndx = 1,size(this%cross_section_objs_)
        if( associated( this%cross_section_objs_(ndx)%val_ ) ) then
          deallocate( this%cross_section_objs_(ndx)%val_ )
        endif
      enddo
      deallocate( this%cross_section_objs_ )
    end if

    if( allocated( this%handles_ ) ) then
      deallocate( this%handles_ )
    end if

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_radXfer_xsect_warehouse
