! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!one dimension, equally spaced  grid type
module micm_1d_from_csv_file_grid

  use musica_constants, only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_1d_grid,     only : base_grid_t

  implicit none

  public :: fromCsvFile_t

  type, extends(base_grid_t) :: fromCsvFile_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type fromCsvFile_t

contains
  !> Initialize grid
  subroutine initialize( this, grid_config )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : die_msg

    !> arguments
    class(fromCsvFile_t), intent(inout) :: this
    type(config_t), intent(inout)       :: grid_config

    !> local variables
    integer(ik), parameter :: Ok = 0_ik
    integer(ik), parameter :: inUnit = 20_ik
    character(len=*), parameter :: Iam = 'From_csv_file grid initialize: '
 
    integer(ik) :: istat
    real(dk)    :: Value
    logical(lk) :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec
    

    !> Get the configuration settings
    call grid_config%get( 'Filespec', Filespec, Iam )
    call grid_config%get( 'Handle', this%handle_, Iam, default = 'None' )

    !> Does input grid file exist?
    inquire( file=Filespec%to_char(), exist=found )
    if( found ) then
      open(unit=inUnit,file=Filespec%to_char(),iostat=istat)
      if( istat == Ok ) then
        read(inUnit,*,iostat=istat) InputLine
        if( istat == Ok ) then
          allocate( this%edge_(0) )
          do
            read(inUnit,*,iostat=istat) Value
            if( istat /= Ok ) then
              exit
            endif
            this%edge_ = [this%edge_,Value]
          enddo
        else
          call die_msg( 560768225, "Error reading " // Filespec%to_char() )
        endif
      endif
    else
      call die_msg( 560768215, "File " // Filespec%to_char() // " not found" )
    endif

    this%ncells_ = size(this%edge_)-1_ik
    allocate( this%mid_(this%ncells_) )
    allocate( this%delta_(this%ncells_) )
    this%mid_(:) = .5_dk &
                   *(this%edge_(1_ik:this%ncells_) + this%edge_(2_ik:this%ncells_+1_ik))
    this%delta_(:) = (this%edge_(2_ik:this%ncells_+1_ik) - this%edge_(1_ik:this%ncells_))

    close(unit=inUnit)

  end subroutine initialize

end module micm_1d_from_csv_file_grid
