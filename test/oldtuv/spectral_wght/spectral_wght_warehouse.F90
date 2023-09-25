! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_spectral_wght_warehouse module

!> The spectral_wght_t type and related functions
!!
module micm_spectral_wght_warehouse

  use micm_abs_spectral_wght_type,     only : abs_spectral_wght_ptr
    use micm_environment, only : environment_t
  use micm_spectral_wght_factory,      only : spectral_wght_builder
  use musica_constants,                only : musica_dk, musica_ik
  use musica_string,                   only : string_t
    use musica_config,                 only : config_t
    use musica_iterator,               only : iterator_t
    use musica_constants,              only : musica_rk
    use musica_assert,                 only : die_msg

  implicit none

  private
  public :: spectral_wght_warehouse_t

  !> Spectral weight warehouse type
  type :: spectral_wght_warehouse_t
    !> spectral wght objects
    type(abs_spectral_wght_ptr), allocatable :: spectral_wght_objs_(:)
    !> Current spectral wght values
    real(kind=musica_dk), allocatable :: spectral_wght_values_(:,:)
    type(string_t), allocatable       :: spectral_wght_key(:)
  contains
    !> Update the object for new environmental conditions
    procedure :: update_for_new_environmental_state
    !> Finalize the object
    final :: finalize
  end type spectral_wght_warehouse_t

  !> radXfer_xsect_t constructor
  interface spectral_wght_warehouse_t
    module procedure :: constructor
  end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Constructor of spectral_wght_t objects
  function constructor( config,mdlLambdaEdge ) result( spectral_wght_warehouse_obj )


    real(musica_dk), intent(in)      :: mdlLambdaEdge(:)
    !> spectral wght configuration data
    type(config_t), intent(inout) :: config
    !> New spectral_wght_warehouse_obj
    class(spectral_wght_warehouse_t), pointer :: spectral_wght_warehouse_obj

    !> local variables
    integer :: nSize, ndx
    character(len=*), parameter :: Iam = "spectral wght warehouse constructor: "
    type(config_t) :: spectral_weight_set, spectrum_config
    type(config_t) :: spectral_wght_config
    class(iterator_t), pointer :: iter
    type(abs_spectral_wght_ptr) :: spectral_wght_ptr
    character(:), allocatable   :: jsonkey
    character(len=32)           :: keychar
    type(string_t)              :: aswkey
    type(string_t), allocatable :: netcdfFiles(:)

    allocate( spectral_wght_warehouse_obj )

    associate(new_obj=>spectral_wght_warehouse_obj)

    allocate( string_t :: new_obj%spectral_wght_key(0) )

    allocate( new_obj%spectral_wght_objs_(0) )

    jsonkey = 'spectral weights'
    call config%get( jsonkey, spectral_weight_set, Iam )
    iter => spectral_weight_set%get_iterator()
!-----------------------------------------------------------------------------
!> iterate over spectral weights
!-----------------------------------------------------------------------------
    do while( iter%next() )
      keychar = spectral_weight_set%key(iter)
      aswkey  = keychar 
      write(*,*) ' '
      write(*,*) Iam,'key = ',trim(keychar)
      new_obj%spectral_wght_key = [new_obj%spectral_wght_key,aswkey]
      call spectral_weight_set%get( iter, spectrum_config, Iam )
!-----------------------------------------------------------------------------
!> get the iterator spectrum
!-----------------------------------------------------------------------------
      call spectrum_config%get( "weights", spectral_wght_config, Iam )
      spectral_wght_ptr%val_ => spectral_wght_builder( spectral_wght_config,mdlLambdaEdge )
      new_obj%spectral_wght_objs_ = [new_obj%spectral_wght_objs_,spectral_wght_ptr]
    end do

    deallocate( iter )

    nSize = size(new_obj%spectral_wght_objs_)
    write(*,*) ' '
    write(*,'(a,''There are '',i3,'' spectral wghts'')') Iam,nSize

!-----------------------------------------------------------------------------
!> setup spectral weight arrays
!-----------------------------------------------------------------------------
    allocate( new_obj%spectral_wght_values_(0,0) )

    end associate

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Update the object for new environmental conditions
  subroutine update_for_new_environmental_state( this, environment, nwave )

    use musica_assert,    only : die_msg

    !> Kinetics
    class(spectral_wght_warehouse_t), intent(inout) :: this
    integer(musica_ik), intent(in)                  :: nwave
    !> Environmental conditions
    class(environment_t), intent(in) :: environment

    character(len=*), parameter :: Iam = 'update_for_new_environmental_state: '
    integer(kind=musica_ik) :: ndx
    real(musica_dk), allocatable :: a_spectral_wght(:)
    real(musica_dk), allocatable :: spectral_wght_tray(:)

    write(*,*) ' '
    write(*,*) Iam,'entering'

    allocate(spectral_wght_tray(0))
    do ndx = 1, size(this%spectral_wght_objs_)
      associate( calc_ftn => this%spectral_wght_objs_(ndx)%val_ )
        a_spectral_wght = calc_ftn%calculate( environment )
      end associate
      spectral_wght_tray = [spectral_wght_tray,a_spectral_wght]
    end do

    this%spectral_wght_values_ = reshape( spectral_wght_tray, &
                                          (/nwave,size(this%spectral_wght_objs_) /) )

    write(*,*) Iam,'size of spectral weight values = ',&
        size(this%spectral_wght_values_,dim=1), size(this%spectral_wght_values_,dim=2)

    write(*,*) Iam,'exiting'

  end subroutine update_for_new_environmental_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Finalize the kinetics
  subroutine finalize( this )

    !> Kinetics
    type(spectral_wght_warehouse_t), intent(inout) :: this

    integer(kind=musica_ik) :: ndx
    character(len=*), parameter :: Iam = 'spectral_wght_warehouse finalize: '

    write(*,*) Iam,'entering'

    if( allocated( this%spectral_wght_values_ ) ) then
      deallocate( this%spectral_wght_values_ )
    endif
    if( allocated( this%spectral_wght_objs_ ) ) then
      do ndx = 1,size(this%spectral_wght_objs_)
        if( associated( this%spectral_wght_objs_(ndx)%val_ ) ) then
          deallocate( this%spectral_wght_objs_(ndx)%val_ )
        endif
      enddo
      deallocate( this%spectral_wght_objs_ )
    end if

    if( allocated( this%spectral_wght_key ) ) then
      deallocate( this%spectral_wght_key )
    end if

    write(*,*) Iam,'exiting'

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_spectral_wght_warehouse
