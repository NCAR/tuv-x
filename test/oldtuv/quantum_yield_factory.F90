! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> The micm_quantum_yield_factory module

!> Builder of quantum yield calculators
module micm_quantum_yield_factory

  use micm_abs_quantum_yield_type,              only : abs_quantum_yield_t
  use micm_base_quantum_yield_type,             only : base_quantum_yield_t
  use micm_tint_quantum_yield_type,             only : tint_quantum_yield_t
  use micm_no2_tint_quantum_yield_type,         only : no2_tint_quantum_yield_t
  use micm_o3_o2_o1d_quantum_yield_type,        only : o3_o2_o1d_quantum_yield_t
  use micm_o3_o2_o3p_quantum_yield_type,        only : o3_o2_o3p_quantum_yield_t
  use micm_ho2_oh_o_quantum_yield_type,         only : ho2_oh_o_quantum_yield_t
  use micm_no3m_aq_quantum_yield_type,          only : no3m_aq_quantum_yield_t
  use micm_ch2o_h2_co_quantum_yield_type,       only : ch2o_h2_co_quantum_yield_t
  use micm_ch3cho_ch3_hco_quantum_yield_type,   only : ch3cho_ch3_hco_quantum_yield_t
  use micm_c2h5cho_c2h5_hco_quantum_yield_type, only : c2h5cho_c2h5_hco_quantum_yield_t
  use micm_ch2chcho_quantum_yield_type,         only : ch2chcho_quantum_yield_t
  use micm_mvk_quantum_yield_type,              only : mvk_quantum_yield_t
  use micm_ch3coch3_ch3co_ch3_quantum_yield_type, only : ch3coch3_ch3co_ch3_quantum_yield_t
  use micm_ch3coch2ch3_quantum_yield_type,      only : ch3coch2ch3_quantum_yield_t
  use micm_ch3cocho_ch3co_hco_quantum_yield_type, only : ch3cocho_ch3co_hco_quantum_yield_t
  use micm_clo_cl_o1d_quantum_yield_type,       only : clo_cl_o1d_quantum_yield_t
  use micm_clo_cl_o3p_quantum_yield_type,       only : clo_cl_o3p_quantum_yield_t
  use micm_clono2_cl_no3_quantum_yield_type,    only : clono2_cl_no3_quantum_yield_t
  use micm_clono2_clo_no2_quantum_yield_type,   only : clono2_clo_no2_quantum_yield_t

  implicit none

  private
  public :: quantum_yield_builder

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Builder of quantum yield calculators
  function quantum_yield_builder( config, mdlLambdaEdge ) result( new_quantum_yield_t )

    use musica_assert,                 only : die_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use musica_constants,              only : musica_dk

    !> New rate constant calculator
    class(abs_quantum_yield_t), pointer :: new_quantum_yield_t
    !> cross section configuration data
    type(config_t), intent(inout) :: config
    real(musica_dk), intent(in)   :: mdlLambdaEdge(:)

    type(string_t) :: quantum_yield_type
    character(len=*), parameter :: Iam = 'quantum yield builder: '

    write(*,*) Iam,'entering'
    new_quantum_yield_t => null()
    call config%get( 'quantum yield type', quantum_yield_type, Iam )

    select case( quantum_yield_type%to_char() )
      case( 'base quantum yield' )
        allocate( base_quantum_yield_t :: new_quantum_yield_t )
      case( 'tint quantum yield' )
        allocate( tint_quantum_yield_t :: new_quantum_yield_t )
      case( 'no2 tint quantum yield' )
        allocate( no2_tint_quantum_yield_t :: new_quantum_yield_t )
      case( 'O3+hv->O2+O(1D) quantum yield' )
        allocate( o3_o2_o1d_quantum_yield_t :: new_quantum_yield_t )
      case( 'O3+hv->O2+O(3P) quantum yield' )
        allocate( o3_o2_o3p_quantum_yield_t :: new_quantum_yield_t )
      case( 'HO2 quantum yield' )
        allocate( ho2_oh_o_quantum_yield_t :: new_quantum_yield_t )
      case( 'NO3-_(aq)+hv->NO2(aq)+O- quantum yield' )
        allocate( no3m_aq_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH2O quantum yield' )
        allocate( ch2o_h2_co_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH3CHO+hv->CH3+HCO quantum yield' )
        allocate( ch3cho_ch3_hco_quantum_yield_t :: new_quantum_yield_t )
      case( 'C2H5CHO quantum yield' )
        allocate( c2h5cho_c2h5_hco_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH2CHCHO+hv->Products quantum yield' )
        allocate( ch2chcho_quantum_yield_t :: new_quantum_yield_t )
      case( 'MVK+hv->Products quantum yield' )
        allocate( mvk_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH3COCH3+hv->CH3CO+CH3 quantum yield' )
        allocate( ch3coch3_ch3co_ch3_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH3COCH2CH3+hv->CH3CO+CH2CH3 quantum yield' )
        allocate( ch3coch2ch3_quantum_yield_t :: new_quantum_yield_t )
      case( 'CH3COCHO+hv->CH3CO+HCO quantum yield' )
        allocate( ch3cocho_ch3co_hco_quantum_yield_t :: new_quantum_yield_t )
      case( 'ClO+hv->Cl+O(1D) quantum yield' )
        allocate( clo_cl_o1d_quantum_yield_t :: new_quantum_yield_t )
      case( 'ClO+hv->Cl+O(3P) quantum yield' )
        allocate( clo_cl_o3p_quantum_yield_t :: new_quantum_yield_t )
      case( 'ClONO2+hv->Cl+NO3 quantum yield' )
        allocate( clono2_cl_no3_quantum_yield_t :: new_quantum_yield_t )
      case( 'ClONO2+hv->ClO+NO2 quantum yield' )
        allocate( clono2_clo_no2_quantum_yield_t :: new_quantum_yield_t )
      case default
        call die_msg( 450768214, "Invalid quantum yield type: '"//              &
                                 quantum_yield_type%to_char( )//"'" )
    end select
    call new_quantum_yield_t%initialize( config, mdlLambdaEdge )
    write(*,*) Iam,'exiting'

  end function quantum_yield_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module micm_quantum_yield_factory
