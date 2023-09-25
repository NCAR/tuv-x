! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This oclo_cross_section module

!> The oclo_cross_section type and related functions
module micm_oclo_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: oclo_cross_section_t

  !> Calculator for oclo_cross_section
  type, extends(base_cross_section_t) :: oclo_cross_section_t
    !> The cross section array
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type oclo_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t

    !> oclo cross section
    class(oclo_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'oclo cross section calculate: '
    integer(musica_ik) :: ndx, nParms
    real(musica_dk)    :: Tfac

    write(*,*) Iam,'entering'

    associate( Temp => environment%temperature, Xsection => this%cross_section )
      nParms = size(Xsection)
      if( Temp <= Xsection(1)%temperature(1) ) then
        cross_section = Xsection(1)%array(:,1)
      elseif( Temp >= Xsection(nParms)%temperature(1) ) then
        cross_section = Xsection(nParms)%array(:,1)
      else
        do ndx = 2,nParms
          if( Xsection(ndx)%temperature(1) > Temp ) then
            exit
          endif
        enddo
        ndx = ndx - 1
        Tfac = (Temp - Xsection(ndx)%temperature(1)) &
               /(Xsection(ndx+1)%temperature(1) - Xsection(ndx)%temperature(1))
        cross_section = Xsection(ndx)%array(:,1) &
                      + Tfac*(Xsection(ndx+1)%array(:,1) - Xsection(ndx)%array(:,1))
      endif
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_oclo_cross_section_type
