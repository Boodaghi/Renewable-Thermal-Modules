# Parabolic Trough Concentrating (PTC) Solar Collector Module

This directory contains the computational model for a linear parabolic concentrating solar collector loop. The engine resolves dynamic optical and thermal losses across varying solar geometries.

## 📋 Contents
* `PTC.ipynb`: Interactive Jupyter Notebook featuring comprehensive inline Markdown documentation, flow loop corrections, and daily parametric performance curves.
* `PTC.py`: Standalone production-ready Python execution script.

## 🔬 Scholarly Reference & Validation Source
The physical correction metrics, structural array factors ($R_1$ and $R_2$), and performance baselines are mapped onto the model components deployed in:

> **Boodaghi, H.**, et al. (2023). *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering**, 48. [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)

## 📊 Core Governing Physics
* **Incidence Angle Modifier (IAM):** Dynamically penalizes nominal optical absorption properties based on the localized solar incidence angle ($\theta$).
* **Flow Divergence Corrections ($R_1$):** Adjusts the manufacturer-reported collector loss metrics ($F_R U_L$) when actual operational mass flows vary from test benchmarks ($g_{\text{test}}$).
* **Series Cascade Modifiers ($R_2$):** Corrects heat removal metrics when tracking multiple collectors connected sequentially in series, factoring in fluid temperature rises along the focal loop.