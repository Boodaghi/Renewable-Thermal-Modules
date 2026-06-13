# Supercritical CO2 Brayton Cycle (SCBC) Topping Module

This directory houses the fluid-property closed simulation framework for a Supercritical $CO_2$ ($sCO_2$) Brayton Topping Power Cycle.

## 📋 Contents
* `SCBC.ipynb`: Interactive Jupyter Notebook featuring high-fidelity thermodynamic property maps, isentropic turbine/compressor loss algorithms, and self-correcting T-s visualizations.
* `SCBC.py`: Standalone production-ready Python execution script.

## 🔬 Scholarly Reference & Validation Source
The cycle pressure boundaries, component state parameters, and performance results are verified against the topping thermodynamic layout published in:

> **Boodaghi, H.**, et al. (2023). *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering**, 48. [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)

## 📊 Core Operational Layout
* **Working Medium:** Supercritical Carbon Dioxide (`sCO2`) handled via high-fidelity Helmholz energy equations of state (CoolProp reference engine).
* **High-Pressure Stage:** 259.0 bar (Gas Turbine Inlet)
* **Low-Pressure Stage:** 74.0 bar (Compressor Inlet)
* **Net Target Output:** 1250.0 kW (1.25 MW)