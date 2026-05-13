# Numerical Investigation of the 1D Transient Heat Conduction Equation

Finite-difference implementation and comparative numerical study of the transient 1D heat equation using Python, NumPy and Matplotlib.

The repository investigates the numerical behavior, stability characteristics and solution accuracy of three classical finite-difference schemes through comparison with the analytical solution.

---

## Governing Equation

Transient heat conduction is governed by:

```text
∂T/∂t = α ∂²T/∂x²
```

where:

- `T(x,t)` → temperature distribution
- `α = 0.02 m²/hr` → thermal diffusivity


### Initial Condition
Triangular temperature profile:
- \(T = 0 \; K\) at \(x = 0\)
- \(T = 100 \; K\) at \(x = 0.5 \; m\)
- \(T = 0 \; K\) at \(x = 1 \; m\)

### Boundary Conditions
\[
T(0,t) = T(1,t) = 0
\]

---

## Numerical Methods

- **Explicit FTCS**
  - forward-time centered-space discretization
  - conditionally stable

- **Implicit Laasonen**
  - tridiagonal implicit formulation
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
- stability behavior
- transient diffusion response
- comparative scheme accuracy
- numerical diffusion effects

---

## Repository Structure

```text
project/
├── explicit_ftcs/
├── implicit_laasonen/
├── crank_nicolson/
├── analytical_solution/
├── plots/
└── results/
```

---

## Current Limitations

- restricted to one-dimensional conduction
- constant thermal diffusivity assumption
- no convergence study implemented yet
- no multidimensional extension

---

## ## Engineering Relevance

Transient diffusion problems provide foundational understanding for:
- PDE discretization
- numerical stability analysis
- scientific computing workflows
- broader CFD and thermal simulation methodologies

---

## Tools Used

<p align="left">

<img src="https://img.shields.io/badge/Python-3776AB?style=flat&logo=python&logoColor=white" />

<img src="https://img.shields.io/badge/NumPy-013243?style=flat&logo=numpy&logoColor=white" />

<img src="https://img.shields.io/badge/Matplotlib-11557C?style=flat" />

</p>

---

## License

MIT License
