// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Delta-Eddington solver for radiative transfer
#pragma once

#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/util/array2d.hpp"
#include "tuvx/util/array3d.hpp"

#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/radiator.hpp>

#include <cassert>
#include <functional>
#include <map>
#include <vector>
namespace tuvx
{

  /// @brief Data structure holding approximation variables for the two stream approximation. These variables are used for
  /// computing the source functions and assembling the tridiagonal system. The different types of approximations can be
  /// found in Table 1 of the paper. For now, only the Eddington approximation is implemented see
  template<typename ArrayPolicy>
  struct ApproximationVariables
  {
    // Two stream approximation coeffcients
    ArrayPolicy gamma1_;
    ArrayPolicy gamma2_;
    ArrayPolicy gamma3_;
    ArrayPolicy gamma4_;
    ArrayPolicy mu_;
    ArrayPolicy lambda_;
    ArrayPolicy gamma_;
  };

  /// @brief Solve Function that solves the radiative transfer equation using two stream approximation
  /// The function first initializes the solver variables (delta scaling, source function definitions), assembles and solves
  /// the tridiagonal linear system and computes the radiation field.
  ///
  /// @param solar_zenith_angles Zenith angles for each column
  /// @param grids Domain grids - wavelength and altitude
  template<
      typename T,
      typename ArrayPolicy,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy>
  void Solve(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::function<void(const RadiatorStatePolicy&, const ArrayPolicy&, const std::vector<T>)> ApproximationFunction,
      const RadiatorStatePolicy& accumulated_radiator_states,
      RadiationFieldPolicy& radiation_field);

  /// @brief Initialize the solver variables
  template<typename T, typename GridPolicy, typename RadiatorStatePolicy, typename SourceFunctionPolicy>
  void InitializeSolver(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      RadiatorStatePolicy& accumulated_radiator_states,
      SourceFunctionPolicy& source_functions);

  /// @brief Compute the Eddington parameters
  template<typename T, typename ArrayPolicy, typename RadiatorStatePolicy>
  void EddingtonApproximation(
      const RadiatorStatePolicy& accumulated_radiator_states,
      const std::vector<T>& solar_zenith_angles,
      ApproximationVariables<ArrayPolicy>& approximation_variables);

  template<typename T, typename GridPolicy, typename RadiatorStatePolicy>
  void ScaleVariables(
      const std::map<std::string, GridPolicy>& grids,
      const std::vector<T>& solar_zenith_angles,
      RadiatorStatePolicy& radiator_state);

  /// @brief BuildSourceFunctions Source functions $C_1$ and $C_2$ functions from the paper
  template<
      typename T,
      typename ArrayPolicy,
      typename GridPolicy,
      typename RadiatorStatePolicy,
      typename SourceFunctionPolicy>
  void BuildSourceFunctions(
      const std::size_t& number_of_columns,
      const std::size_t& number_of_wavelengths,
      const std::size_t& number_of_layers,
      const std::vector<T>& solar_zenith_angles,
      const RadiatorStatePolicy& accumulated_radiator_states,
      const ApproximationVariables<ArrayPolicy> approximation_variables,
      const ArrayPolicy& surface_reflectivity,
      std::map<std::string, SourceFunctionPolicy> source_functions);

};  // namespace tuvx

#include "delta_eddington.inl"
#include "initialize_solver.inl"
/*
include "initialize_solver.inl";
include "assemble_tridiagonal_system.inl";
include "compute_radiation_field.inl";
include "solve.inl"
*/
