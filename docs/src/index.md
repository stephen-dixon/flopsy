# Flopsy.jl Documentation

## What is Flopsy?

Flopsy is a Julia package for experimenting with stiff ODE solvers applied to one-dimensional reaction-diffusion (RD) problems. It is specifically designed for hydrogen transport through materials, leveraging the Palioxis Fortran library (via the Hotgates wrapper) to compute trapping and detrapping reaction rates with high fidelity.

The package provides a composable architecture for building complex coupled systems: reaction operators, diffusion operators, and constraint operators combine flexibly to form different problem formulations, which can then be solved using robust implicit and semi-implicit time integrators from SciML.

## Target Problems

Flopsy targets stiff, spatially-distributed hydrogen transport problems where:

- **Stiffness** arises from disparate timescales between mobile hydrogen diffusion and multi-occupancy trapping dynamics.
- **Heterogeneity** stems from spatial variation in material properties, defect densities, and temperatures.
- **Multi-occupancy trapping** requires faithful representation of defect interactions, accessible through the Palioxis library.
- **1D geometry** is sufficient for many practical materials-science applications (thin films, implantation, permeation).

## Quick Links

- **[Architecture](architecture.md)** — Understand the code structure: variable layouts, operators, formulations, and the solve pipeline.
- **[Hotgates Interface](hotgates.md)** — Learn how Flopsy bridges to Palioxis for accurate reaction-rate calculations.
- **API Reference** — Full API documentation (generated from docstrings).
