# Silicon Photovoltaic (PV) Solar Cell Module

This directory contains the computational framework for simulating a high-fidelity single-diode equivalent circuit PV array coupled with an explicit dynamic lumped capacitance thermal model.

## 📋 Contents
* `PV.ipynb`: Interactive Jupyter Notebook featuring scholarly Markdown equations, inline code annotations, and characteristic diagnostic sweeps.
* `PV.py`: Standard standalone production runtime execution script.

## 🔬 Scholarly Reference & Validation Source
The physical constants, governing mathematical matrices, and specific transient parameters map directly onto the solar-harvest baseline engine derived in:

> **Boodaghi, H.**, et al. (2026). *"Achieving holistic sustainability in solar-hydrogen systems: A 6E-based multi-objective optimization of a PV–PEMFC–PEME–ORC integrated framework."* **Thermal Science and Engineering Progress**, 72, 104773. [DOI: 10.1016/j.tsep.2026.104773](https://doi.org/10.1016/j.tsep.2026.104773)

## 📊 Core Governing Physics
* **Implicit Electrical Domain:** Resolves p-n junction non-linear shock losses via an implicit transcendental current loop solved dynamically using Newton-Raphson approximations (`scipy.optimize.fsolve`).
* **Transient Thermal Domain:** Quantifies real-time array operating temperature fluctuations ($C_t \frac{dT}{dt}$) by tracking localized solar radiative gains against simultaneous electrical extraction and utility convective heat losses via adaptive Runge-Kutta numerical integration (`solve_ivp`).