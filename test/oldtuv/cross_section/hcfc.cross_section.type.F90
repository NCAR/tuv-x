! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This hcfc+hv->products cross_section module

!> The hcfc+hv->products_cross_section type and related functions
module micm_hcfc_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: hcfc_cross_section_t

  !> Calculator for acetone cross_section
  type, extends(base_cross_section_t) :: hcfc_cross_section_t
  contains
    !> Initialize the cross section
    procedure :: calculate => run
  end type hcfc_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )
!     qyacet - q.y. for acetone, based on Blitz et al. (2004)
! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!       (2004), Pressure and temperature-dependent quantum yields for the 
!       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.

    use micm_environment,  only : environment_t

    !> base cross section
    class(hcfc_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter  :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter  :: LBar  = 206.214_musica_dk
    character(len=*), parameter :: Iam = 'ch3coch3+hv->ch3co+ch3 cross section calculate: '
    integer(musica_ik) :: wNdx, polyNdx
    real(musica_dk)    :: Tadj, sigma, uLambda

    uLambda = this%cross_section(1)%temperature(2)
    Tadj = min( 295._musica_dk,max( 203._musica_dk,environment%temperature ) ) &
           - this%cross_section(1)%temperature(1)
    do wNdx = 1,size(this%mdl_lambda_center)
      if( this%mdl_lambda_center(wNdx) >= 190._musica_dk &
          .and. this%mdl_lambda_center(wNdx) <= uLambda ) then
        sigma = rZERO
        associate( coefficient => this%cross_section(1)%array )
        do polyNdx = 1,size(coefficient,dim=1)
          sigma = sigma &
                + (coefficient(polyNdx,1) + Tadj*(coefficient(polyNdx,2) + Tadj*coefficient(polyNdx,3))) &
                  * (this%mdl_lambda_center(wNdx) - LBar)**(polyNdx-1)
        enddo
        end associate
        sigma = exp( sigma )
      else
        sigma = rZERO
      endif
      cross_section(wNdx) = sigma
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_hcfc_cross_section_type
