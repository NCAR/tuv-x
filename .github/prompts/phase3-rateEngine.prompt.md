# Phase 3: Rate Calculation Engine (Steps 13–15)

See [master plan](plan-tuvXCppSolverRewrite.prompt.md) for overall architecture and key decisions.

## Context

Once the radiation field is solved (Phase 1) and transforms are defined (Phase 2), rates are computed by a three-step pipeline: (1) **calculate** each transform (cross-section, quantum yield, or spectral weight) to produce weight matrices `[λ × z × col]`, (2) **apply** the weights to the radiation field via element-wise multiplication, (3) **reduce** by summing over wavelengths to produce rates `[z × col]`. These are all weighted spectral integrals — highly parallelizable across columns and reactions.

## Step 13: Implement photolysis rate calculator

Port `src/photolysis_rates.F90`. The core operation:

$$J_i(z) = \sum_\lambda F(\lambda, z) \cdot \sigma_i(\lambda, z) \cdot \phi_i(\lambda, z) \cdot \Delta\lambda$$

where $F$ is actinic flux, $\sigma_i$ is cross-section, $\phi_i$ is quantum yield, and $\Delta\lambda$ is wavelength bin width.

Implementation:
- **Calculate**: evaluate each reaction's cross-section and quantum yield transforms at the current atmospheric state to produce weight matrices
- **Apply**: compute element-wise product of actinic flux × σ weights × φ weights × Δλ using `ForEachRow` and column views
- **Reduce**: sum over wavelength dimension (reduction per height layer per column)
- Output: 3D array `[reaction × height × column]`
- Use `ArrayPolicy::Function` to pre-compile the per-timestep rate calculation for optimal vectorization across columns

API:
```cpp
template<typename ArrayPolicy>
class PhotolysisCalculator {
    void add_reaction(std::string name, TransformFunc<ArrayPolicy> cross_section, TransformFunc<ArrayPolicy> quantum_yield);
    Array3D<typename ArrayPolicy::value_type> calculate(const RadiationField<ArrayPolicy>& field, const AtmosphericState<ArrayPolicy>& state);
    // calculate() internally: (1) calculates each transform's weights from `state`,
    // (2) applies weights to `field` (element-wise multiply), (3) reduces over wavelengths → returns rates [reaction × height × column]
};
```

Parallelization: reactions are independent — can be computed in parallel across reactions and columns.

### Fortran reference
- `src/photolysis_rates.F90` — `photolysis_rates_t` class, `calculate()` method uses `dot_product(actinicFlux, cross_section * quantum_yield)` per height layer

## Step 14: Implement dose rate calculator

Port `src/dose_rates.F90`. Similar structure:

$$D(z) = \sum_\lambda I_{irrad}(\lambda, z) \cdot W(\lambda)$$

where $I_{irrad}$ combines direct, upwelling, and downwelling spectral irradiance (converted to $W \cdot m^{-2}$), and $W(\lambda)$ is a spectral weighting function.

Implementation uses `matmul(irradiance_matrix, spectral_weight_vector)` per height layer.

### Fortran reference
- `src/dose_rates.F90` — `dose_rates_t` class

## Step 15: Implement heating rate calculator

Port `src/heating_rates.F90`.

### Fortran reference
- `src/heating_rates.F90` — `heating_rates_t` class

### Documentation

Document `PhotolysisCalculator`, `DoseRateCalculator`, and `HeatingRateCalculator` with Doxygen `///` comments. Explain the calculate/apply/reduce pipeline in the class-level doc. Reference the Fortran source.
