// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Delta-Eddington solver for radiative transfer
#pragma once

#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/util/array2d.hpp"

#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/radiator.hpp>

#include <cassert>
#include <cmath>
#include <functional>
#include <map>
#include <vector>
namespace tuvx
{

  /// @brief Radiative flux calculator that applies the delta-Eddington Approximation.
  ///
  /// [DEV NOTES] We can determine whether this should be a class or a set of functions

  /// @brief Solve the radiative transfer equation for a collection of columns
  /// @param solar_zenith_angles Solar zenith angles for each column [radians].
  /// @param grids Grids available for the radiative transfer calculation.
  /// @param profiles Profiles available for the radiative transfer calculation.
  /// @param radiation_field The calculated radiation field.
  ///
  /// Solves two-stream equations for multiple layers. These routines are based
  /// on equations from: Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.
  /// DOI: https://doi.org/10.1029/JD094iD13p16287
  /// It contains 9 two-stream methods to choose from. A pseudo-spherical
  /// correction has also been added.
  ///
  /// The original delta-Eddington paper is:
  /// Joseph and Wiscombe, J. Atmos. Sci., 33, 2453-2459, 1976
  /// DOI: https://doi.org/10.1175/1520-0469(1976)033%3C2452:TDEAFR%3E2.0.CO;2
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
      const std::map<std::string, ProfilePolicy>& profiles,
      const RadiatorStatePolicy& accumulated_radiator_states,
      const std::function<void(const RadiatorStatePolicy&, const ArrayPolicy&)> approximation_function,
      RadiationFieldPolicy& radiation_field);

  /// @brief Compute the Eddington parameters
  template<typename T, typename RadiatorStatePolicy, typename ArrayPolicy>
  void DeltaEddingtonApproximation(const RadiatorStatePolicy& radiator_state, const ArrayPolicy& parameters);

  template<typename T, typename GridPolicy, typename ArrayPolicy, typename ProfilePolicy, typename RadiatorStatePolicy>
  void InitializeVariables(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::map<std::string, std::vector<T>>& solver_variables,
      const RadiatorStatePolicy& accumulated_radiator_states);

  template<typename T>
  void BuildSourceFunctions(
      std::map<std::string, std::vector<T>> solver_variables,
      std::map<std::string, std::function<T(T)>> source_functions);

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

  template<
      typename T,
      typename ArrayPolicy,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy>
  inline void Solve(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::function<void(const RadiatorStatePolicy&, const ArrayPolicy&, const std::vector<T>)> ApproximationFunction,
      const RadiatorStatePolicy& accumulated_radiator_state,
      RadiationFieldPolicy& radiation_field)

};  // namespace tuvx
include "delta_eddington.inl";
include "initialize_solver.inl";
include "assemble_tridiagonal_system.inl";
include "compute_radiation_field.inl";
include "solve.inl"
