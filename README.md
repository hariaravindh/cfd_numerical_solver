# Numerical Investigation of the 1D Transient Heat Conduction Equation

Comparative finite-difference study of the transient 1D heat equation using Python, NumPy and Matplotlib.

The repository investigates numerical stability, transient diffusion behavior and solution accuracy across multiple finite-difference schemes through comparison with the analytical solution.

---

## Governing Equation

```text
∂T/∂t = α ∂²T/∂x²
```

where:
- `T(x,t)` → temperature distribution
- `α = 0.02 m²/hr` → thermal diffusivity

Initial condition:
- triangular temperature profile
- `T = 0 K` at `x = 0` and `x = 1 m`
- `T = 100 K` at `x = 0.5 m`

Boundary condition:
- `T(0,t) = T(1,t) = 0 K`

---

## Numerical Methods

The following finite-difference schemes are implemented and compared:

- **Explicit FTCS**
  - forward-time centered-space formulation
  - conditionally stable

- **Implicit Laasonen**
  - implicit tridiagonal formulation
  - unconditionally stable

- **Crank–Nicolson**
  - second-order implicit scheme
  - improved transient accuracy

---

## Validation & Investigation

All numerical schemes are benchmarked against the analytical solution for multiple test cases involving varying:
- mesh sizes
- timestep values
- Fourier number conditions

Current investigations include:
- timestep sensitivity
- numerical stability behavior
- transient diffusion response
- explicit scheme instability
- comparative scheme accuracy

---

## Repository Structure

```text
project/
│
├── core/
│   ├── analytical_solution.py
│   └── intial_comditions.py
│
├── plots/
│   ├── test_case_1.png
│   ├── test_case_2.png
│   ├── test_case_3.png
│   ├── test_case_4.png
│   ├── test_case_5.png
│   ├── test_case_6.1.png
│   ├── test_case_6.2.png
│   └── test_case_6.png
│
├── solvers/
│   ├── crank_nicolson.py
│   ├── explicit_ftcs.py
│   └── implicit_laasonen.py
│
├── studies/
│   ├── case_1.py
│   ├── case_2.py
│   ├── case_3.py
│   ├── case_4.py
│   └── case_5&6.py
│
└── docs/
│   
├── LICENSE
├── README.md

```

---

## Current Limitations

- restricted to one-dimensional conduction
- constant thermal diffusivity assumption
- no multidimensional extension
- no adaptive meshing
- convergence study not yet implemented

---

## Engineering Relevance

This study investigates foundational numerical behaviors that appear throughout computational heat transfer and broader PDE-based simulation methods.

---

## Tools Used

<p align="left">

<img src="https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white" />

<img src="https://img.shields.io/badge/NumPy-013243?style=flat&logo=numpy&logoColor=white" />

<img src="https://img.shields.io/badge/Matplotlib-11557C?style=flat" />

</p>

---

![License: MIT](https://img.shields.io/badge/license-MIT-green?style=flat-square)
