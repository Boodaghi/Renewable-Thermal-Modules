# Thermal Energy Storage (TES) indirect Two-Tank Module

This directory houses the transient numerical simulation engine for an indirect two-tank sensible Thermal Energy Storage (TES) loop using Therminol VP-1 as the working storage fluid.

## 📋 Contents
* `TES.ipynb`: Interactive Jupyter Notebook featuring explicit transient tank heat balances, storage volume depletions, and state-of-charge (SoC) graphics.
* `TES.py`: Standalone production-ready Python execution script.

## 🔬 Scholarly Reference & Validation Source
The physical configurations, temperature thresholds, and storage fluid properties map directly onto the backup thermal storage system evaluated in:

> **Boodaghi, H.**, et al. (2023). *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering**, 48. [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)

## 📊 Core Operational Metrics
* **Storage Fluid Medium:** Therminol VP-1 Synthetic Liquid Heat Transfer Fluid.
* **Hot Tank Charging Upper Limit:** 390.0 °C 
* **Cold Tank Return Lower Limit:** 290.0 °C
* **Numerical Method:** Runge-Kutta Adaptive Step Initial Value Integration (`scipy.integrate.solve_ivp`).