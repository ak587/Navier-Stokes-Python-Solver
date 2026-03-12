# 2D Lid-Driven Cavity Flow Solver (CFD)
[![Language](https://img.shields.io/badge/Language-Python-blue.svg)](https://www.python.org/)
[![Library](https://img.shields.io/badge/Library-NumPy-orange.svg)](https://numpy.org/)
[![Library](https://img.shields.io/badge/Library-Matplotlib-green.svg)](https://matplotlib.org/)

## 📌 Project Overview
This project implements a numerical solver for the **2D Lid-Driven Cavity Flow**, a classic benchmark problem in Computational Fluid Dynamics (CFD). The solver simulates the behavior of an incompressible fluid inside a square cavity where the top wall moves at a constant velocity, creating a primary vortex and smaller secondary eddies.

The simulation solves the **Navier-Stokes equations** using the **Vorticity-Stream Function ($\omega-\psi$) formulation**.

### Key Features
*   **Numerical Schemes:** 1st-order Upwind scheme for convective terms (stability) and 2nd-order Central Difference for diffusive terms.
*   **Iterative Solver:** Jacobi Iteration method for solving the Poisson equation for the stream function.
*   **Stability:** Adaptive time-stepping based on the **CFL (Courant–Friedrichs–Lewy) condition**.
*   **Validation:** Results are validated against the landmark study by **Ghia et al. (1982)**.
*   **Performance:** Capable of simulating high Reynolds numbers (up to $Re=3200$) on a $129 \times 129$ grid.

---

## 🧪 Mathematical Formulation
The governing equations for incompressible flow in the $\omega-\psi$ form are:

1.  **Vorticity Transport Equation:**
    $$\frac{\partial \omega}{\partial t} + u \frac{\partial \omega}{\partial x} + v \frac{\partial \omega}{\partial y} = \frac{1}{Re} \left( \frac{\partial^2 \omega}{\partial x^2} + \frac{\partial^2 \omega}{\partial y^2} \right)$$
2.  **Stream Function (Poisson Equation):**
    $$\frac{\partial^2 \psi}{\partial x^2} + \frac{\partial^2 \psi}{\partial y^2} = -\omega$$
3.  **Velocity Definitions:**
    $$u = \frac{\partial \psi}{\partial y}, \quad v = -\frac{\partial \psi}{\partial x}$$

---

## 📊 Results & Validation
The solver computes the steady-state velocity profiles and compares them directly with experimental benchmark data.

### Reynolds Number: 1000
*(Insert your generated image here. Example path: `assets/results_re1000.png`)*
![Lid Driven Cavity Re 1000](https://via.placeholder.com/800x300.png?text=Placeholder:+Upload+your+Lid_Driven_Cavity_Re_1000.png+here)

**Observations:**
*   **Streamlines:** The primary vortex shifts toward the center as $Re$ increases.
*   **U-Velocity Profile:** The horizontal velocity across the vertical centerline matches the Ghia et al. data points with high accuracy.
*   **V-Velocity Profile:** The vertical velocity across the horizontal centerline captures the peak magnitudes correctly.

---

## 💻 Implementation Details
The code is written in Python, prioritizing readability and scientific accuracy:
*   **`poisson_solver_for_psi`**: Solves the elliptic Poisson equation using Jacobi iteration with a convergence tolerance of $10^{-3}$.
*   **`explicit_solver_for_omega`**: Updates the vorticity field using an explicit time-marching scheme.
*   **Adaptive Time-Stepping**: 
    ```python
    convective_limit = 1.0 / ((np.max(np.abs(u))/dx) + (np.max(np.abs(v))/dy) + 1e-12)
    diffusive_limit  = Re * dx**2 / 4
    dt = 0.2 * min(convective_limit, diffusive_limit)
    ```

---

## 🚀 How to Run
1.  **Clone the repository:**
    ```bash
    git clone https://github.com/yourusername/lid-driven-cavity-cfd.git
    ```
2.  **Install dependencies:**
    ```bash
    pip install numpy matplotlib
    ```
3.  **Execute the solver:**
    ```bash
    python solver.py
    ```

---

## 🎓 Author
**[Akash Mishra]**
*   LinkedIn: [linkedin.com/in/ak587](https://linkedin.com/in/ak587)
