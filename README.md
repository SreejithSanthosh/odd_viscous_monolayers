# Odd viscous monolayers: flow and stress around nematic defects

MATLAB code accompanying

> **Achiral Odd Mechanics in Cell Monolayers**
> Sreejith Santhosh and Mattia Serra, Department of Physics, University of California San Diego.

The code computes the steady-state velocity, strain-rate and stress fields of a two-dimensional viscous monolayer with a prescribed nematic texture (an isolated ±1/2 topological defect), for a linear constitutive relation that contains an *achiral odd* viscous modulus $`\lambda^o`$. It reproduces the defect flow and stress profiles in Fig. 2, Fig. 3 and Fig. S2 of the paper.

## Model

The monolayer is a 2D viscous continuum on the unit disk with a prescribed nematic tensor $`\mathbf{Q}(\mathbf{x}) = S(\mathbf{x})\,(\mathbf{n}\mathbf{n}^\top - \mathbf{I}/2)`$. Writing the stress and strain-rate tensors in the basis $`\{\mathbf{D}_1, \mathbf{D}_2, \mathbf{D}_3\}`$ (isotropic part, and the two traceless-symmetric parts), the nondimensional constitutive relation is

```math
\begin{pmatrix}\sigma_1\\ \sigma_2\\ \sigma_3\end{pmatrix}
=
\begin{pmatrix} m_0\\ 0\\ 0\end{pmatrix}
+
\begin{pmatrix}
1 & (\lambda^s-\lambda^o)\,q_1 & (\lambda^s-\lambda^o)\,q_2\\
(\lambda^s+\lambda^o)\,q_1 & \nu + \delta\nu^1(\mathbf{Q}) & \delta\nu^2(\mathbf{Q})\\
(\lambda^s+\lambda^o)\,q_2 & \delta\nu^2(\mathbf{Q}) & \nu - \delta\nu^1(\mathbf{Q})
\end{pmatrix}
\begin{pmatrix} e_1\\ e_2\\ e_3\end{pmatrix},
\qquad
\nabla\cdot\boldsymbol{\sigma} = 0,
\qquad
\boldsymbol{\sigma}\,\hat{\mathbf{n}} = -\kappa\,\mathbf{u}\ \text{ on } r=1,
```

with $`(q_1, q_2) = (Q_{11}, Q_{12})`$, $`e_1 = \tfrac12\mathrm{Tr}(\nabla\mathbf{u})`$, $`e_2 = \tfrac12(\partial_x u_x - \partial_y u_y)`$, $`e_3 = \tfrac12(\partial_x u_y + \partial_y u_x)`$ (Eq. S17 of the paper). Stress is nondimensionalized by the bulk viscosity, so the (1,1) entry is 1. The parameters are

| symbol | meaning |
|---|---|
| $`m_0`$ | isotropic prestress (contractility) |
| $`\nu`$ | shear viscosity |
| $`\lambda^s`$ | symmetric (reciprocal) coupling between isotropic and deviatoric modes |
| $`\lambda^o`$ | **odd** (non-reciprocal) coupling between isotropic and deviatoric modes |
| $`\delta\nu`$ | anisotropy of the shear viscosity |
| $`\kappa`$ | boundary drag; $`\kappa \ll 1`$ makes the boundary approximately traction-free while removing the rigid-body null space |

The nematic texture of a defect of charge $`s`$ at the origin is $`\phi = s\theta + \phi_0`$ and $`S(r) = 1 - \exp\left(-r^2/2K^2\right)`$, where $`K`$ is the defect core size and $`(r,\theta)`$ are polar coordinates.

The problem is solved with the MATLAB Partial Differential Equation Toolbox in its standard form $`-\nabla\cdot(\mathbf{c}\otimes\nabla\mathbf{u}) = \mathbf{f}`$ with $`\mathbf{f}=0`$ and generalized Neumann condition $`\hat{\mathbf{n}}\cdot(\mathbf{c}\otimes\nabla\mathbf{u}) + \mathbf{q}\,\mathbf{u} = \mathbf{g}`$. The 16-entry coefficient tensor $`\mathbf{c}`$ is the constitutive matrix above rotated into the lab frame (Eq. S12), $`\mathbf{q} = \kappa\mathbf{I}`$ and $`\mathbf{g} = -m_0\hat{\mathbf{n}}`$, which together give $`\boldsymbol{\sigma}\hat{\mathbf{n}} = -\kappa\mathbf{u}`$ for the total stress including the prestress.

## Requirements

* MATLAB R2019b or newer (uses `tiledlayout`).
* Partial Differential Equation Toolbox (`createpde`, `circleg`, `solvepde`, `interpolateSolution`, `evaluateGradient`).

No other toolboxes or external dependencies are needed.

## Quick start

1. Clone the repository and open the folder in MATLAB.
2. Edit `get_parameters.m` to set the moduli, and `orient_order.m` to choose the defect charge.
3. Run

   ```matlab
   main
   ```

`main.m` builds the unit-disk geometry, assembles the coefficients, meshes the domain (`Hmax = 0.01`), solves for the steady-state velocity, interpolates the solution and its gradient onto a 100 × 100 Cartesian grid, and produces a six-panel figure. A run takes on the order of a minute on a laptop; increase `hmax` in `main.m` for a quicker, coarser solve.

To save the figure, add for example `exportgraphics(fig, 'defect.png', 'Resolution', 300)` at the end of `main.m`. The handle `ax(5)` is the mean-stress panel if you want to fix its color limits (`clim(ax(5), [...])`) when comparing runs.

## Files

| file | purpose |
|---|---|
| `main.m` | Driver script: geometry, coefficients, boundary conditions, mesh, solve, post-processing and plotting. |
| `get_parameters.m` | All physical parameters (moduli, boundary drag $`\kappa`$, prestress $`m_0`$). Edit this file to change the material. |
| `orient_order.m` | Prescribed nematic field: returns the order parameter $`S`$ and director angle $`\phi`$ at given points. Set the defect charge (`n = ±0.5`), the core size and the defect position/orientation here. |
| `coeff.m` | The 16-entry PDE-Toolbox coefficient tensor $`\mathbf{c}(\mathbf{x})`$, i.e. the lab-frame viscous-moduli tensor for the local $`\mathbf{Q}(\mathbf{x})`$. |
| `bc_drag.m` | Boundary coefficient $`\mathbf{q} = \kappa\mathbf{I}`$ (weak boundary drag). |
| `bc_sfbc.m` | Boundary coefficient $`\mathbf{g} = -m_0\hat{\mathbf{n}}`$ (stress-free condition for the total stress with isotropic prestress). |
| `compute_stress.m` | Evaluates the constitutive relation: strain-rate components $`(e_1,e_2,e_3)`$ → viscous stress components $`(\sigma_1,\sigma_2,\sigma_3)`$ in the $`\mathbf{D}_i`$ basis. |
| `plot_fields.m` | Plots the six panels described below. |

### Parameter names in `get_parameters.m`

| variable | symbol | default |
|---|---|---|
| `bulk_visc` | bulk viscosity (normalization) | 1 |
| `shear_visc` | $`\nu`$ | 1 |
| `lambda_0` | $`\lambda^o`$ (odd modulus) | 0 |
| `lambda_s` | $`\lambda^s`$ | 0.2 |
| `del_nu` | $`\delta\nu`$ | 0 |
| `kappa` | $`\kappa`$ | $`10^{-4}`$ |
| `m_iso` | $`m_0`$ | 1 |

## Output

`plot_fields.m` draws one row of six panels on the unit disk:

1. $`\mathbf{Q}(\mathbf{x})`$: order parameter $`S`$ (color) with the director $`\mathbf{n}`$ (magenta rods).
2. $`\mathbf{u}(\mathbf{x})`$: speed (color) and velocity vectors.
3. $`\mathrm{Tr}[\nabla\mathbf{u}]/2 = e_1`$: isotropic strain rate (positive for expansion, negative for compression).
4. $`e_d,\ \mathbf{p}`$: magnitude and principal axis of the deviatoric strain rate.
5. $`\mathrm{Tr}[\boldsymbol{\sigma}]/2 = \sigma_1 + m_0`$: mean (isotropic) stress including the prestress.
6. $`\sigma_d,\ \mathbf{p}_\sigma`$: magnitude and principal axis of the deviatoric stress.

The super-title lists the moduli used for the run.

## Reproducing the paper figures

All simulations in the paper use $`R = 1`$, $`K = 0.05`$, $`m_0 = 1`$, $`\nu = 1`$, $`\delta\nu = 0`$ and $`\kappa = 10^{-4}`$ unless stated otherwise. Set `lambda_0` ($`\lambda^o`$) and `lambda_s` ($`\lambda^s`$) in `get_parameters.m` and `n` in `orient_order.m`:

| figure | $`\lambda^o`$ | $`\lambda^s`$ | defect |
|---|---|---|---|
| Fig. 2B | 0.5 | 0.5 | +1/2 and −1/2 |
| Fig. 2C | −0.5 | −0.5 | +1/2 and −1/2 |
| Fig. 2D, Fig. 3A | 1 | 0 | +1/2 and −1/2 |
| Fig. 3B | 0 | −0.5 | +1/2 and −1/2 |
| Fig. 3C | 0 | 0.5 | +1/2 and −1/2 |
| Fig. S2 | see caption | see caption | +1/2 and −1/2 |

Robustness with respect to the mesh size and to $`\kappa`$ (varied over $`10^{-5}`$–$`10^{-3}`$) was checked by changing `hmax` in `main.m` and `kappa` in `get_parameters.m`.

## Citation

If you use this code, please cite the paper:

```bibtex
@article{santhosh2026odd,
  title   = {Achiral Odd Mechanics in Cell Monolayers},
  author  = {Santhosh, Sreejith and Serra, Mattia},
  year    = {2026},
  note    = {University of California San Diego}
}
```

## Contact

Sreejith Santhosh (ssanthos@ucsd.edu) and Mattia Serra (mserra@ucsd.edu), Department of Physics, UC San Diego.
