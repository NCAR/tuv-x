! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!from csv file type
module micm_from_csv_file_vert_Profile

  use musica_constants,  only : dk => musica_dk, ik => musica_ik, lk => musica_lk
  use micm_vert_Profile, only : abs_vert_Profile_t

  implicit none

  public :: fromCsvFile_t

  type, extends(abs_vert_Profile_t) :: fromCsvFile_t
  contains
    !> Initialize grid
    procedure :: initialize
  end type fromCsvFile_t

contains
  !> Initialize grid
  subroutine initialize( this, profile_config, zGrid )
      
    use musica_config, only : config_t
    use musica_string, only : string_t
    use musica_assert, only : die_msg
    use interpolation

    !> arguments
    class(fromCsvFile_t), intent(inout) :: this
    type(config_t), intent(inout)       :: profile_config
    real(dk), intent(in)                :: zGrid(:)

    !> local variables
    integer(ik), parameter :: Ok = 0_ik
    integer(ik), parameter :: inUnit = 20_ik
    character(len=*), parameter :: Iam = 'From_csv_file vert_Profile initialize: '
 
    integer(ik) :: istat
    real(dk)    :: zd, Value
    real(dk), allocatable :: zdata(:)
    real(dk), allocatable :: Profile(:)
    logical(lk) :: found
    character(len=132) :: InputLine
    type(string_t)     :: Filespec

    write(*,*) Iam // 'entering'

    !> Get the configuration settings
    call profile_config%get( 'Filespec', Filespec, Iam )
    call profile_config%get( 'Handle', this%handle_, Iam, default = 'None' )

    !> Does input grid file exist?
    inquire( file=Filespec%to_char(), exist=found )
    file_exists: if( found ) then
      open(unit=inUnit,file=Filespec%to_char(),iostat=istat)
      if( istat == Ok ) then
        !> Skip the header
        do
          read(inUnit,'(a)',iostat=istat) InputLine
          if( istat /= Ok ) then
            exit
          elseif( verify( InputLine(1:1),'#!$%*' ) /= 0 ) then
            exit
          endif
        enddo
        if( istat == Ok ) then
          allocate( Profile(0) )
          allocate( zdata(0) )
          !> Read the data
          do
            read(InputLine,*,iostat=istat) zd, Value
            if( istat /= Ok ) then
              call die_msg( 560768229, "Invalid data format in " // Filespec%to_char() )
            endif
            Profile = [Profile,Value]
            zdata = [zdata,zd]
            read(inUnit,'(a)',iostat=istat) InputLine
            if( istat /= Ok ) then
              exit
            endif
          enddo
        else
          call die_msg( 560768227, "Error reading " // Filespec%to_char() )
        endif
      else
        call die_msg( 560768231, "Error opening " // Filespec%to_char() )
      endif
    else file_exists
      call die_msg( 560768215, "File " // Filespec%to_char() // " not found" )
    endif file_exists

    this%ncells_ = size(zGrid) - 1_ik
    allocate( this%edge_val_(this%ncells_+1_ik) )
    this%edge_val_ = this%inter1( zGrid, zdata,Profile )

    write(*,*) ' '
    write(*,*) Iam // 'data z grid'
    write(*,'(1p10g15.7)') zdata
    write(*,*) ' '
    write(*,*) Iam // this%handle_%to_char() // ' on data z grid'
    write(*,'(1p10g15.7)') Profile
    write(*,*) ' '
    write(*,*) Iam // this%handle_%to_char() // ' @ mdl z grid edges'
    write(*,'(1p10g15.7)') this%edge_val_

    allocate( this%mid_val_(this%ncells_) )
    allocate( this%delta_val_(this%ncells_) )
    this%mid_val_(:) = .5_dk &
                   *(this%edge_val_(1_ik:this%ncells_) + this%edge_val_(2_ik:this%ncells_+1_ik))
    this%delta_val_(:) = (this%edge_val_(2_ik:this%ncells_+1_ik) - this%edge_val_(1_ik:this%ncells_))

    close(unit=inUnit)

    write(*,*) Iam // 'exiting'

  end subroutine initialize

end module micm_from_csv_file_vert_Profile
