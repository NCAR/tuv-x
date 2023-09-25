! Copyright (C) 2021 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
!> \file
!> Tests the photo_decomp radXfer module

!> Test module for the radXfer_core_t type
program radXfer_test

  use musica_string,    only : string_t
  use radXfer_core,     only : radXfer_core_t

  implicit none

  class(radXfer_core_t), pointer :: core

  !> Command-line options
  character(len=256) :: argument
  type(string_t)     :: configFileSpec

  !> Get the model configuration file and options from the command line
  argument = 'test/data/radiative_transfer.test.config.json'

  configFileSpec = argument

  !> instatiate and initialize radXfer core object
  core => radXfer_core_t( configFileSpec )

  !> set radXfer cross sections
  call core%test()

  deallocate( core )

  contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Fail run and print usage info
  subroutine fail_run( )

    write(*,*) "Usage: ./radiative_transfer_test configuration_file.json"
    stop 3

  end subroutine fail_run

end program radXfer_test
