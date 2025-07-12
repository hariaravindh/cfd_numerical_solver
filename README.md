<h1 align="center"> Numerical Solver for 1-D Unsteady Heat Conduction </h1>


##  Project Overview

This repository contains the Python implementation for **Computational Fluid Dynamics** , which solves the **One-dimensional unsteady heat conduction equation** in a homogeneous wall.

---

##  Objective

Implement and compare **Three finite difference methods**:
- ✅ **Explicit FTCS**
- ✅ **Implicit Laasonen**
- ✅ **Crank-Nicolson**

All methods are benchmarked against the **analytical solution** for six test cases to analyze **accuracy**, **stability**, and **efficiency**.

---

## Features

-  **Three Numerical Methods:** Explicit FTCS, Implicit Laasonen, Crank-Nicolson.
-  **Analytical Benchmark:** Compares numerical results to an analytical series solution.
-  **Multiple Test Cases:** Six configurations with varying mesh sizes and Courant numbers.
-  **Automated Plots:** Generates solution plots for visual comparison.
-  **Simple, Modular Python Code:** Easy to extend or adapt.

---

##  Governing Equation

<p align="center">
∂T/∂t = α ∂²T/∂x²
</p>

- **T(x, t)**: Temperature distribution  
- **α = 0.02 m²/hr**: Thermal diffusivity

**Initial & Boundary Conditions:**
- Initial: Triangular profile — T = 0 K at x = 0 and x = 1 m; T = 100 K at x = 0.5 m
- Boundaries: T(0, t) = T(1, t) = 0 K

---

## ⚙️ Numerical Methods

### Explicit FTCS  
- Forward time, centered space.
###  Implicit Laasonen  
- Solves a tridiagonal system.
###  Crank-Nicolson  
- Averages time levels.

---

<p align="center">
  <img src="https://img.shields.io/badge/Language-Python-blue?style=flat-square" alt="Language">
  <img src="https://img.shields.io/badge/NumPy-013243?style=flat&logo=numpy" />
  <img src="https://img.shields.io/badge/Matplotlib-F7931E?style=flat&logo=matplotlib" />
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" alt="License">
</p>

---
