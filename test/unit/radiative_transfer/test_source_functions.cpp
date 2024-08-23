#include "tuvx/grid.hpp"
#include "tuvx/radiative_transfer/radiation_field.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"
#include "tuvx/util/array2d.hpp"
#include "tuvx/util/array3d.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

TEST(SOURCE_FUNCTION, SCALING)
{
  const std::size_t n_columns = 10;
  const std::size_t n_wavelengths = 10;
  const std::size_t n_layers = 10;

  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid("m", n_wavelengths, n_columns);
  tuvx::Grid<tuvx::Array2D<double>> vertical_grid("m", n_layers, n_columns);
  std::vector<double> solar_zenith_angles(n_columns, 90);

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>> grids{ { "altitude [m]", vertical_grid },
                                                                  { "wavelength [m]", wavelength_grid } };

  tuvx::RadiatorState<tuvx::Array3D<double>> accumulated_radiator_states(n_columns, vertical_grid, wavelength_grid);

  std::map<std::string, tuvx::Array3D<double>> solver_variables;
  std::map<std::string, std::function<double(double)>> source_functions;

  tuvx::InitializeVariables<
      double,
      tuvx::Grid<tuvx::Array2D<double>>,
      tuvx::RadiatorState<tuvx::Array3D<double>>,
      tuvx::RadiationField<tuvx::Array3D<double>>,
      std::function<double(double)>,
      tuvx::Array3D<double>>(solar_zenith_angles, grids, accumulated_radiator_states, solver_variables, source_functions);

  EXPECT_EQ(0, 0);
}
