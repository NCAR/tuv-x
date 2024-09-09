#include "tuvx/grid.hpp"
#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"
#include "tuvx/util/array2d.hpp"
#include "tuvx/util/array3d.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <gtest/gtest.h>

#include <functional>
#include <vector>

TEST(TUVX_INITIALIZE_SOLVER, SCALING)
{
  const std::size_t n_columns = 10;
  const std::size_t n_wavelengths = 3;
  const std::size_t n_layers = 4;

  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid("m", n_wavelengths, n_columns);
  tuvx::Grid<tuvx::Array2D<double>> vertical_grid("m", n_layers, n_columns);
  std::vector<double> solar_zenith_angles(n_columns, 90);

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>> grids{ { "altitude [m]", vertical_grid },
                                                                  { "wavelength [m]", wavelength_grid } };

  // initialize random radiator state
  tuvx::RadiatorState<tuvx::Array3D<double>> accumulated_radiator_states(n_columns, vertical_grid, wavelength_grid);
  tuvx::FillRandom(accumulated_radiator_states.optical_depth_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.single_scattering_albedo_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.asymmetry_parameter_.AsVector());

  // Scale variables for the Delta-Eddington approximation
  tuvx::ScaleVariables(grids, solar_zenith_angles, accumulated_radiator_states);

  EXPECT_EQ(0, 0);
}

TEST(TUVX_INITIALIZE_SOLVER, DeltaEddingtonApproximationVariables)
{
  const std::size_t number_of_layers = 10;
  const std::size_t number_of_columns = 3;
  const std::size_t number_of_wavelengths = 4;

  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid("m", number_of_wavelengths, number_of_columns);
  tuvx::Grid<tuvx::Array2D<double>> vertical_grid("m", number_of_layers, number_of_columns);
  std::vector<double> solar_zenith_angles(number_of_columns, 90);

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>> grids{ { "altitude [m]", vertical_grid },
                                                                  { "wavelength [m]", wavelength_grid } };

  // source functions (C1, C2 from the paper)
  std::map<std::string, tuvx::Array3D<std::function<double(double)>>> source_functions;

  // initialize random radiator state
  tuvx::RadiatorState<tuvx::Array3D<double>> accumulated_radiator_states(number_of_columns, vertical_grid, wavelength_grid);
  tuvx::FillRandom(accumulated_radiator_states.optical_depth_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.single_scattering_albedo_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.asymmetry_parameter_.AsVector());

  // compute the delta eddington approximation variables
  tuvx::ApproximationVariables<tuvx::Array3D<double>> approximation_variables;
  tuvx::EddingtonApproximation<double, tuvx::Array3D<double>, tuvx::RadiatorState<tuvx::Array3D<double>>>(
      number_of_columns,
      number_of_wavelengths,
      number_of_layers,
      accumulated_radiator_states,
      solar_zenith_angles,
      approximation_variables);

  EXPECT_EQ(0, 0);
}

TEST(TUVX_INITIALIZE_SOLVER, SourceFunctions)
{
  const std::size_t number_of_layers = 10;
  const std::size_t number_of_columns = 3;
  const std::size_t number_of_wavelengths = 4;

  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid("m", number_of_wavelengths, number_of_columns);
  tuvx::Grid<tuvx::Array2D<double>> vertical_grid("m", number_of_layers, number_of_columns);
  std::vector<double> solar_zenith_angles(number_of_columns, 90);

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>> grids{ { "altitude [m]", vertical_grid },
                                                                  { "wavelength [m]", wavelength_grid } };

  // source functions (C1, C2 from the paper)
  std::map<std::string, tuvx::Array3D<std::function<double(double)>>> source_functions;
  source_functions["C_upwelling"] =
      tuvx::Array3D<std::function<double(double)>>(number_of_columns, number_of_wavelengths, number_of_layers);
  source_functions["C_downwelling"] =
      tuvx::Array3D<std::function<double(double)>>(number_of_columns, number_of_wavelengths, number_of_layers);

  // initialize random radiator state
  tuvx::RadiatorState<tuvx::Array3D<double>> accumulated_radiator_states(number_of_columns, vertical_grid, wavelength_grid);
  tuvx::FillRandom(accumulated_radiator_states.optical_depth_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.single_scattering_albedo_.AsVector());
  tuvx::FillRandom(accumulated_radiator_states.asymmetry_parameter_.AsVector());

  // surface reflectivity
  tuvx::Array3D<double> surface_reflectivity(number_of_columns, number_of_wavelengths, number_of_layers);
  tuvx::FillRandom(surface_reflectivity.AsVector());

  // compute the delta eddington approximation variables
  tuvx::ApproximationVariables<tuvx::Array3D<double>> approximation_variables;
  tuvx::EddingtonApproximation<double, tuvx::Array3D<double>, tuvx::RadiatorState<tuvx::Array3D<double>>>(
      number_of_columns,
      number_of_wavelengths,
      number_of_layers,
      accumulated_radiator_states,
      solar_zenith_angles,
      approximation_variables);

  // Compute source function
  tuvx::BuildSourceFunctions<
      double,
      tuvx::Array3D<double>,
      tuvx::RadiatorState<tuvx::Array3D<double>>,
      tuvx::Array3D<std::function<double(double)>>>(
      number_of_columns,
      number_of_wavelengths,
      number_of_layers,
      solar_zenith_angles,
      accumulated_radiator_states,
      approximation_variables,
      surface_reflectivity,
      source_functions);

  // print a runtime generated source function
  std::cout << source_functions.at("C_upwelling")(0, 0, 0)(0.5) << ", " << source_functions.at("C_downwelling")(0, 0, 0)(0.5)
            << std::endl;

  EXPECT_TRUE(true);
}
