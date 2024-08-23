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

  /// @brief Compute the Eddington parameters
  template<typename T, typename RadiatorStatePolicy, typename ArrayPolicy>
  void DeltaEddingtonApproximation(const RadiatorStatePolicy& radiator_state, const ArrayPolicy& parameters);

  template<typename T, typename GridPolicy, typename ArrayPolicy, typename RadiatorStatePolicy>
  void InitializeSolver(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, std::vector<T>>& solver_variables,
      const RadiatorStatePolicy& accumulated_radiator_states);

  template<typename T, typename GridPolicy, typename ArrayPolicy, typename SourceFunctionPolicy>
  void BuildSourceFunctions(
      GridPolicy grids,
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, ArrayPolicy>& solver_variables,
      std::map<std::string, SourceFunctionPolicy> source_functions);

  template<typename T, typename ArrayPolicy, typename GridPolicy>
  void
  ScaleVariables(GridPolicy grids, std::vector<T> solar_zenith_angles, std::map<std::string, ArrayPolicy>& solver_variables);

  /*
  template<typename T>
  void AssembleTridiagonalMatrix(
      std::size_t number_of_layers,
      const std::map<std::string, T> solver_variables,
      TridiagonalMatrix<T>& coeffcient_matrix);

  template<typename T>
  void AssembleCoeffcientVector(
      std::size_t number_of_layers,
      const std::map<std::string, std::vector<T>> solver_variables,
      std::map<std::string, std::function<T(T)>> source_functions,
      std::vector<T>& coeffcient_vector);

  template<typename T, typename RadiationFieldPolicy>
  void ComputeRadiationField(
      const std::vector<T>& solar_zenith_angles,
      std::map<std::string, std::vector<T>> solver_variables,
      const TridiagonalMatrix<T>& coeffcient_matrix,
      const std::vector<T>& coeffcient_vector,
      const RadiationFieldPolicy& radiation_field);
  */

};  // namespace tuvx
#include "delta_eddington.inl";
#include "initialize_solver.inl";
/*
include "delta_eddington.inl";
include "initialize_solver.inl";
include "assemble_tridiagonal_system.inl";
include "compute_radiation_field.inl";
include "solve.inl"
*/
