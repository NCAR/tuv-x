! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_spectral_wght_factory module

!> Builder of cross section calculators
module micm_spectral_wght_factory

  use micm_abs_spectral_wght_type,          only : abs_spectral_wght_t
  use micm_base_spectral_wght_type,         only : base_spectral_wght_t
  use micm_uv_b_280_315_nm_spectral_wght_type, only : uv_b_280_315_nm_spectral_wght_t
  use micm_uv_b_280_320_nm_spectral_wght_type, only : uv_b_280_320_nm_spectral_wght_t
  use micm_uv_a_315_400_nm_spectral_wght_type, only : uv_a_315_400_nm_spectral_wght_t
  use micm_visplus_spectral_wght_type,         only : visplus_spectral_wght_t
  use micm_gaussian_305_nm_10_nm_FWHM_spectral_wght_type,  only : gaussian_305_nm_10_nm_FWHM_spectral_wght_t
  use micm_gaussian_320_nm_10_nm_FWHM_spectral_wght_type,  only : gaussian_320_nm_10_nm_FWHM_spectral_wght_t
  use micm_gaussian_340_nm_10_nm_FWHM_spectral_wght_type,  only : gaussian_340_nm_10_nm_FWHM_spectral_wght_t
  use micm_gaussian_380_nm_10_nm_FWHM_spectral_wght_type,  only : gaussian_380_nm_10_nm_FWHM_spectral_wght_t
  use micm_eppley_uv_photometer_spectral_wght_type,        only : eppley_uv_photometer_spectral_wght_t
  use micm_par_400_700nm_spectral_wght_type,               only : par_400_700nm_spectral_wght_t
  use micm_exponential_decay_spectral_wght_type,           only : exponential_decay_spectral_wght_t
  use micm_scup_mice_spectral_wght_type,                   only : scup_mice_spectral_wght_t
  use micm_standard_human_erythema_spectral_wght_type,     only : standard_human_erythema_spectral_wght_t
  use micm_uv_index_spectral_wght_type,                    only : uv_index_spectral_wght_t
  use micm_phytoplankton_boucher_spectral_wght_type,       only : phytoplankton_boucher_spectral_wght_t
  use micm_plant_damage_spectral_wght_type,                only : plant_damage_spectral_wght_t
  use micm_plant_damage_flint_caldwell_spectral_wght_type, only : plant_damage_flint_caldwell_spectral_wght_t
  use micm_plant_damage_flint_caldwell_ext_spectral_wght_type, only : plant_damage_flint_caldwell_ext_spectral_wght_t

  implicit none

  private
  public :: spectral_wght_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function spectral_wght_builder( config, mdlLambdaEdge ) result( new_spectral_wght_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use musica_constants,              only : musica_dk

    !> New rate constant calculator
    class(abs_spectral_wght_t), pointer :: new_spectral_wght_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)

    type(string_t) :: spectral_wght_type
    character(len=*), parameter :: Iam = 'spectral wght builder: '

    write(*,*) Iam,'entering'
    new_spectral_wght_t => null()
    call config%get( 'spectral wght type', spectral_wght_type, Iam )

    select case( spectral_wght_type%to_char() )
      case( 'base spectral wght' )
        allocate( base_spectral_wght_t :: new_spectral_wght_t )
      case( 'UV-B,280-315nm' )
        allocate( uv_b_280_315_nm_spectral_wght_t :: new_spectral_wght_t )
      case( 'UV-B*,280-320nm' )
        allocate( uv_b_280_320_nm_spectral_wght_t :: new_spectral_wght_t )
      case( 'UV-A,315-400nm' )
        allocate( uv_a_315_400_nm_spectral_wght_t :: new_spectral_wght_t )
      case( 'vis+,> 400nm' )
        allocate( visplus_spectral_wght_t :: new_spectral_wght_t )
      case( 'Gaussian, 305 nm, 10 nm FWHM' )
        allocate( gaussian_305_nm_10_nm_FWHM_spectral_wght_t :: new_spectral_wght_t )
      case( 'Gaussian, 320 nm, 10 nm FWHM' )
        allocate( gaussian_320_nm_10_nm_FWHM_spectral_wght_t :: new_spectral_wght_t )
      case( 'Gaussian, 340 nm, 10 nm FWHM' )
        allocate( gaussian_340_nm_10_nm_FWHM_spectral_wght_t :: new_spectral_wght_t )
      case( 'Gaussian, 380 nm, 10 nm FWHM' )
        allocate( gaussian_380_nm_10_nm_FWHM_spectral_wght_t :: new_spectral_wght_t )
      case( 'Eppley UV Photometer' )
        allocate( eppley_uv_photometer_spectral_wght_t :: new_spectral_wght_t )
      case( 'PAR, 400-700nm, umol m-2 s-1' )
        allocate( par_400_700nm_spectral_wght_t :: new_spectral_wght_t )
      case( 'Exponential decay' )
        allocate( exponential_decay_spectral_wght_t :: new_spectral_wght_t )
      case( 'SCUP-mice (de Gruijl et al., 1993)' )
        allocate( scup_mice_spectral_wght_t :: new_spectral_wght_t )
      case( 'Standard human erythema(Webb et al., 2011)' )
        allocate( standard_human_erythema_spectral_wght_t :: new_spectral_wght_t )
      case( 'UV Index(WMO, 1994; Webb et al., 2011)' )
        allocate( uv_index_spectral_wght_t :: new_spectral_wght_t )
      case( 'Phytoplankton (Boucher et al., 1994)' )
        allocate( phytoplankton_boucher_spectral_wght_t :: new_spectral_wght_t )
      case( 'Plant damage (Caldwell, 1971)' )
        allocate( plant_damage_spectral_wght_t :: new_spectral_wght_t )
      case( 'Plant damage,Flint&Caldwell,2003,orig' )
        allocate( plant_damage_flint_caldwell_spectral_wght_t :: new_spectral_wght_t )
      case( 'Plant damage,Flint&Caldwell,2003,ext390' )
        allocate( plant_damage_flint_caldwell_ext_spectral_wght_t :: new_spectral_wght_t )
      case default
        call die_msg( 450768215, "Invalid spectral wght type: '"//              &
                                 spectral_wght_type%to_char()//"'" )
    end select
    call new_spectral_wght_t%initialize( config, mdlLambdaEdge )
    write(*,*) Iam,'exiting'

  end function spectral_wght_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_spectral_wght_factory
