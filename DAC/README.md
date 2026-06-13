# Direct Air Capture (DAC) Solid Sorbent Module

This directory contains the computational simulation framework for a solid-sorbent Direct Air Capture (DAC) module designed to extract carbon dioxide directly from atmospheric air using low-temperature thermal regeneration.

## 📋 Contents
* `DAC.ipynb`: Interactive Jupyter Notebook featuring mass transfer kinetics, fan power evaluation matrices, and daily capture yield profiles.
* `DAC.py`: Production-ready standalone Python execution script.

## 🔬 Core Modeling Framework
* **Adsorption/Desorption Thermal Swing:** Captures ambient $\text{CO}_2$ and utilizes low-temperature thermal energy ($85\text{--}100^\circ\text{C}$) for sorbent regeneration, matching the heat rejection profile of low-grade industrial/cycle utility streams.
* **Parasitic Load Metrics:** Quantifies specific thermal energy consumption ($\text{GJ/t}_{\text{CO2}}$) and mechanical fan power requirements based on contactor pressure drops.