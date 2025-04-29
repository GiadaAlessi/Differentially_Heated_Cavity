# Differentially Heated Cavity (DHC)

This project implements a **C++ numerical solver** for the **Differentially Heated Cavity** problem (DHC), based on the **Fractional Step Method** and using a **staggered mesh** for improved stability and accuracy.

## Problem Description

The differentially heated cavity problem involves a two-dimensional square cavity, where vertical walls are kept at different constant temperatures, and horizontal walls are thermally insulated. Natural convection is driven by the temperature gradient, resulting in a complex coupled flow and temperature field.

## Methodology

- **Discretization**:  
  - Finite Volume Method (FVM) on a staggered grid  
  - Central Difference Scheme (CDS) for diffusion terms  
  - Central Difference Scheme (CDS) or Upwind Scheme (UDS) for convective terms (switchable)
  
- **Time Advancement**:  
  - Fractional Step Method (projection method) for pressure-velocity coupling  
  - Adaptive time stepping based on residual evolution

- **Numerical Stabilization**:  
  - Under-relaxation for velocity **u**, **v** and temperature **T** fields  
  - Temperature "freezing" when residual stagnation is detected

- **Convergence Monitoring**:  
  - Residual tracking for velocity and temperature fields  
  - Automatic stopping criteria based on residual thresholds

- **Post-Processing Outputs**:  
  - Field data for temperature, velocities, and pressure  
  - Key quantities for Grid Independence Study:
    - Average Nusselt number **Nu**
    - Maximum velocity components **u_max**, **v_max**

- **Post-Processing and Plotting**:  
  - All colormaps, residual plots, and grid independence graphs were generated using **Python** via the script `Diff_Heat_Cav_Plots.py`

## References
G. de Vahl Davis, "Natural Convection of Air in a Square Cavity: A Benchmark Numerical Solution", International Journal for Numerical Methods in Fluids, 1983.

## Author
**Giada Alessi**  
Master in Thermal Engineering  
Universitat Politècnica de Catalunya
