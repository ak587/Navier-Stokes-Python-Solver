# Lid-Driven Cavity

A Python implementation of the classic **lid-driven cavity flow** problem, a benchmark case in computational fluid dynamics (CFD). This solver simulates incompressible Navier-Stokes flow in a square cavity with the top wall moving at constant velocity.

## Overview

The lid-driven cavity is a fundamental test case for validating CFD algorithms. It features:
- A unit square domain (1 × 1)
- A moving top wall (lid) at velocity U = 1
- Three stationary bottom, left, and right walls
- Viscous incompressible flow governed by the Navier-Stokes equations

The code solves the problem across multiple Reynolds numbers (100, 400, 1000, 3200) and validates results against the classical **Ghia et al. benchmark data**.

## Implementation Details

### Mathematical Approach

The solver uses the **stream function-vorticity (ψ-ω) formulation** of the Navier-Stokes equations:

1. **Vorticity Transport Equation**: Solved explicitly using an upwind scheme for convective terms and central differences for diffusive terms
2. **Poisson Equation for Stream Function**: Solved iteratively using the Jacobi iteration method
3. **Velocity Recovery**: Velocities are computed from the stream function using finite differences

### Key Features

- **Adaptive Time Stepping**: CFL-based stability control to ensure convergence
- **Boundary Conditions**: Proper treatment of vorticity and stream function at walls
- **Validation**: Results compared with Ghia et al. (1982) benchmark data for U and V velocity profiles
- **Multiple Reynolds Numbers**: Supports Re = 100, 400, 1000, and 3200

### Numerical Methods

| Component | Method |
|-----------|--------|
| Convective Terms | Upwind scheme |
| Diffusive Terms | Central difference |
| Poisson Solver | Jacobi iteration |
| Time Integration | Explicit (forward Euler) |
| Time Stepping | Adaptive (CFL-based) |

## Files

- **`Lid-Driven_Cavity.py`** - Main solver implementation
- **`Lid_Driven_Cavity_Re_*.png`** - Output plots showing streamlines and velocity comparisons for different Reynolds numbers
  - Re = 100 (t = 8.40)
  - Re = 400 (t = 44.16)
  - Re = 1000 (t = 95.47)
  - Re = 3200 (t = 400.00)

## Requirements
