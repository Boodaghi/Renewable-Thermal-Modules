# Recuperative Organic Rankine Cycle (RORC) Module

This directory contains the variable, fluid-flexible numerical simulation framework for a Recuperative Organic Rankine Cycle (RORC) equipped with an explicit condenser cooling utility loop.

## 📋 Contents
* `RORC.ipynb`: Interactive Jupyter Notebook featuring high-fidelity code comments, Markdown descriptions, and validation logs.
* `RORC.py`: Standalone production-ready Python execution script.

## 🔬 Scholarly Reference & Validation Source
The governing thermodynamic states, boundary constraints, and fluid equations are directly mapped onto the design parameters published in:

> **Boodaghi, H.**, et al. (2023). *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering**. [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)

## 📊 Core Simulation Parameters
* **Target Net Output:** $850 \text{ kW}$ ($W_{\text{net}}$)
* **Working Medium:** Fully fluid-agnostic lookups via CoolProp (Default: Toluene)
* **High-Side Pressure Limit:** Dynamically bound to $90\%$ of the fluid's critical pressure ($P_{\text{crit}}$) to maintain stable, subcritical cycle operation.