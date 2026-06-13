# Renewable Thermal Modules Portfolio

Welcome to the **Renewable-Thermal-Modules** repository. This open-access library houses high-fidelity, validated numerical simulation models for advanced thermodynamic cycles, low-temperature solar collectors, carbon capture, and thermal energy storage subsystems. 

Every module is translated from legacy Engineering Equation Solver (EES) and MATLAB environments into robust, self-contained Python architectures using state-of-the-art equations of state.

---

## 🔬 Scholarly References & Validation Sources
The physics configurations, component state nodes, and efficiency metrics throughout this repository are rigorously verified against empirical data and layout baselines published in:

1. **Boodaghi, H.**, et al. (2023). *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering**, 48. [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)
2. **Boodaghi, H.**, et al. (2025). *"Transient Performance Optimization of a Novel Atmospheric Water Harvester Coupled with Desiccant Bed Matrices."* **Applied Thermal Engineering**, 256, 124110.

---

## 📂 Repository Architecture & Submodule Mapping

This workspace is explicitly structured into modular sub-sections. Each component subdirectory contains its own dedicated technical documentation, validation scripts, and interactive visual workbooks:

```text
Renewable-Thermal-Modules/
├── CITATION.cff           # Global citation registry for scholarly indexing
├── README.md              # Global repository front-page portfolio documentation
├── requirements.txt       # Unified environment package installation manifest
├── AWH/
│   ├── AWH.ipynb          # 11-hour forward-Euler transient desiccant bed engine
│   ├── AWH.py             # Production standalone mathematical sweeping script
│   └── README.md          # Subfolder tracker tied to Applied Thermal Engineering (2025)
├── DAC/
│   ├── DAC.ipynb          # Solid-sorbent Direct Air Capture mass and power consumption model
│   ├── DAC.py             # Production standalone optimization sweep script
│   └── README.md          # Subfolder tracker mapping technical capture specifications
├── Oxygen_Storage/
│   ├── Storage.ipynb      # Stewart-Jacobsen Fundamental EOS oxygen vessel model
│   ├── Storage.py         # Production standalone compressed oxygen runtime script
│   └── README.md          # Subfolder tracker mapping safety boundaries (Max 150 bar)
├── PTC/
│   ├── PTC.ipynb          # Flow rate deviation & sequential cascade loop modifier notebook
│   ├── PTC.py             # Standalone parametric configuration runtime script
│   └── README.md          # Subfolder validation tracking tied to Arabian Journal (2023)
├── RORC/
│   ├── RORC.py            # High-throughput production optimization sweep script
│   ├── RORC.ipynb         # Fluid-independent mass split and pinch safeguard workbook
│   └── README.md          # Subfolder validation tracking tied to Arabian Journal poly-gen (2023)
├── SCBC/
│   ├── SCBC.ipynb         # Near-critical sCO2 topping cycle isobaric mapping workbook
│   ├── SCBC.py            # High-throughput production optimization sweep script
│   └── README.md          # Subfolder validation tracking tied to Arabian Journal topping cycle (2023)
└── TES/
    ├── TES.ipynb          # Transient indirect two-tank Therminol VP-1 energy storage model
    ├── TES.py             # Production standalone thermodynamic flux tracking script
    └── README.md          # Subfolder validation tracking tied to Arabian Journal (2023)