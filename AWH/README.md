# Adsorption-Based Atmospheric Water Harvesting (AWH) Module

This directory houses the transient, variable fluid numerical simulation framework for a solar adsorption-based solid desiccant Atmospheric Water Harvesting (AWH) array loop.

## 📋 Contents
* `AWH.ipynb`: Primary interactive Jupyter Notebook featuring academic Markdown equations, finite-step forward-Euler solvers, and diurnal trajectory graphs.
* `AWH.py`: Stan-alone production-ready Python execution script.

## 🔬 Scholarly Reference & Validation Source
The governing dynamic thermal energy matrices, desorption mass splits, and Antoine saturation variables are directly validated against the empirical states published in:

>*"Experimental study of a solar adsorption-based atmospheric water harvesting system for off-grid cogeneration."* **Applied Thermal Engineering**, 264, 127744. [DOI: 10.1016/j.applthermaleng.2025.127744](https://doi.org/10.1016/j.applthermaleng.2025.127744)

## 📊 Core Performance Metrics
* **Array Configuration:** 27 standalone modular adsorption beds wired in a localized collection pattern.
* **Kinetic Evaluation Engine:** Couples macro-hourly changing climate conditions with second-by-second transient mass-transfer evaluations ($\Delta t = 1\text{ s}$).
* **Performance Indicators:** Tracks absolute water mass flow rates ($\dot{m}_{\text{water}}$), silica temperature profiles ($T_{\text{silica}}$), and Specific Moisture Production (SMP).