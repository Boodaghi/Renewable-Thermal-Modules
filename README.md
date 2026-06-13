# Renewable Thermal Modules

A high-fidelity, open-source Python suite for simulating, validating, and optimizing advanced thermodynamic and power cycle conversion configurations. This repository provides verified, fluid-agnostic numerical engineering engines cross-referenced with real-fluid thermodynamic lookup properties via CoolProp.

---

## 📚 Associated Publications

These models serve as the verified computational baselines for the following peer-reviewed research papers:

1. **H. Boodaghi**, et al. *"Design and Performance Assessment of a Novel Poly-generation System with Stable Production of Electricity, Hydrogen, and Hot Water: Energy and Exergy Analyses."* **Arabian Journal for Science and Engineering** (2023). [DOI: 10.1007/s13369-023-08410-7](https://doi.org/10.1007/s13369-023-08410-7)

---

## 🔬 Component & Cycle Frameworks

### 1. Regenerative Organic Rankine Cycle (`RORC/`)
An optimized numerical simulation engine for a Recuperative Organic Rankine Cycle (RORC) equipped with an explicit condenser cooling loop. 

Key architectural elements implemented:
* **Fluid-Agnostic Flex Stability:** Directly integrated with CoolProp real-fluid datasets. High-side pressure targets adjust dynamically to $90\%$ of the fluid's true critical pressure ($P_{\text{crit}}$).
* **Second-Law Thermal Boundary Safeguard:** Computes cycle mass splits dynamically based on targeted net power goals ($W_{\text{net}} = 850 \text{ kW}$) while enforcing minimum pinch point limits ($\Delta T_{\text{rec,min}} = 20\text{ K}$) to prevent internal temperature crosses.

---

## 🛠️ Dependencies & Local Execution

To run these thermodynamic cycle computations locally, install the core scientific dependencies via your package manager:

```bash
pip install numpy coolprop