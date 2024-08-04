
// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace tuvx
{

  template<typename T>
  void AssembleTridiagonalMatrix(
      std::size_t number_of_layers,
      const std::map<std::string, std::vector<T>> solver_variables,
      TridiagonalMatrix<T>& coeffcient_matrix)
  {
    // get linear system size
    std::size_t matrix_size = 2 * number_of_layers;
    {
      // LEFT HAND SIDE coeffcient matrix diagonals
      std::vector<T>& upper_diagonal = coeffcient_matrix.upper_diagonal_;
      std::vector<T>& main_diagonal = coeffcient_matrix.main_diagonal_;
      std::vector<T>& lower_diagonal = coeffcient_matrix.lower_diagonal_;

      // extract internal variables to build the matrix
      const std::vector<T>& e1 = solver_variables.at("e1");
      const std::vector<T>& e2 = solver_variables.at("e2");
      const std::vector<T>& e3 = solver_variables.at("e3");
      const std::vector<T>& e4 = solver_variables.at("e4");

      //  surface reflectivity
      const T& R_sfc = solver_variables.at("Surface Reflectivity");

      // first row
      upper_diagonal.front() = 0;
      main_diagonal.front() = e1.front();
      lower_diagonal.front() = -e2.front();

      // odd rows
      for (std::size_t n = 1; n < matrix_size - 1; n += 2)
      {
        upper_diagonal[n] = e2[n + 1] * e1[n] - e3[n] * e4[n + 1];
        main_diagonal[n] = e2[n] * e2[n + 1] - e3[n] * e4[n + 1];
        lower_diagonal[n] = e3[n] * e4[n + 1] - e1[n + 1] * e2[n + 1];
      }

      // even rows
      for (std::size_t n = 2; n < matrix_size - 2; n += 2)
      {
        upper_diagonal[n] = e2[n] * e3[n] - e4[n] * e1[n];
        main_diagonal[n] = e1[n] * e1[n + 1] - e3[n] * e3[n + 1];
        lower_diagonal[n] = e3[n] * e4[n + 1] - e1[n + 1] * e2[n + 1];
      }

      // last row
      lower_diagonal.back() = e1.back() - R_sfc * e3.back();
      main_diagonal.back() = e2.back() - R_sfc * e4.back();
      upper_diagonal.back() = 0;
    }
  }

  template<typename T>
  void AssembleCoeffcientVector(
      std::size_t number_of_layers,
      const std::map<std::string, std::vector<T>> solver_variables,
      const std::map<std::string, std::vector<T>> source_terms,
      std::vector<T>& coeffcient_vector)
  {
    // set up tridiagonal matrix
    std::size_t matrix_size = 2 * number_of_layers;
    {
      //  internal variables to build the matrix
      const std::vector<T>& e1 = solver_variables.at("e1");
      const std::vector<T>& e2 = solver_variables.at("e2");
      const std::vector<T>& e3 = solver_variables.at("e3");
      const std::vector<T>& e4 = solver_variables.at("e4");

      //  surface reflectivity and flux source
      const std::vector<T>& R_sfc = solver_variables.at("Surface Reflectivity");
      const std::vector<T>& f_0 = solver_variables.at("Initial Source Flux");
      const std::vector<T>& tau = solver_variables.at("Optical Depth");
      const std::vector<T>& S_sfc = solver_variables.at("Source Reflectivity");

      //  source functions
      // this is a vector of std::functions
      const auto& C_upwelling = source_terms.at("C_upwelling");
      const auto& C_downwelling = source_terms.at("C_downwelling");

      // first row
      coeffcient_vector[0] = f_0 - C_downwelling[0](tau[0]);

      // odd rows
      for (std::size_t n = 1; n < matrix_size - 1; n += 2)
      {
        coeffcient_vector[n] = e3[n] * (C_upwelling[n + 1](0) - C_upwelling[n](tau[n])) +
                               e1[n] * (C_downwelling[n](tau[n]) - C_downwelling[n](0));
      }

      // even rows
      for (std::size_t n = 2; n < matrix_size - 2; n += 2)
      {
        coeffcient_vector[n] = e2[n + 1] * (C_upwelling[n + 1](0) - C_upwelling[n](tau[n])) +
                               e4[n + 1] * (C_downwelling[n + 1](0) - C_downwelling[n + 1](tau[n]));
      }

      // last row
      coeffcient_vector.back() = S_sfc - C_upwelling.back()(tau.back()) + R_sfc * C_downwelling.back()(tau.back());
    }
  }
}  // namespace tuvx
