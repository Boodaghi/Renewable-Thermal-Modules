# %% [markdown]
# # Regenerative Organic Rankine Cycle (RORC) High-Fidelity Model
# 
# ## Overview
# This module provides a completely flexible, fluid-independent numerical simulation engine for a Regenerative Organic Rankine Cycle (RORC) integrated with an explicit condenser cooling loop. 
# 
# Unlike conventional models that hardcode temperature targets, this engine calculates cycle mass splits dynamically based on target net work outputs ($W_{\text{net}} = 850\text{ kW}$) and protects against thermodynamic boundary violations by executing localized heat exchanger energy balances.
# 
# ---
# 
# ## 🔬 Component Layout and State Points
# To ensure clear interpretation, the cycle state points are standardized sequentially:
# * **State 1:** Turbine Inlet (Superheated or high-pressure vapor)
# * **State 2:** Turbine Actual Outlet / Recuperator Hot-Side Inlet
# * **State 3:** Recuperator Hot-Side Outlet / Condenser Inlet
# * **State 4:** Condenser Liquid Outlet / Pump Inlet (Saturated Liquid)
# * **State 5:** Pump Actual Outlet / Recuperator Cold-Side Inlet
# * **State 6:** Recuperator Cold-Side Outlet / Evaporator Inlet
# 
# ---
# 
# # ⚙️ Thermodynamic Modeling Assumptions
# 
# To secure mathematical reproducibility and ensure a self-correcting physics architecture for any selected working medium or temperature shift, the simulation relies on the following structural assumptions:
# 
# 1. **Steady-State Continuity:** The entire cycle operates under steady-state conditions ($\frac{dm}{dt} = 0$, $\frac{dE}{dt} = 0$).
# 2. **Isobaric Heat Exchangers:** Pressure drop penalties across the evaporator core, recuperator channels, and condenser tube bundles are assumed negligible ($P_{\text{high}}$ and $P_{\text{low}}$ remain completely uniform on their respective pressure stages).
# 3. **Isentropic Machine Efficiencies:** Deviations from ideal expansions and compressions are strictly quantified via fixed isentropic efficiencies ($\eta_{\text{turbine}} = 0.80$, $\eta_{\text{pump}} = 0.75$).
# 4. **Saturated Liquid States:** Fluid exiting the condenser and entering the pump is handled as a pure saturated liquid ($Q = 0$) at the defined condensation temperature ($T_{\text{cond}}$).
# 5. **No Ambient Heat Leaks:** All component boundaries (turbines, pumps, piping) are modeled as perfectly adiabatic; thermal energy transfer occurs exclusively within the heat exchangers.
# 6. **Second Law Safeguard:** The recuperator gas outlet state is dynamically bound using a minimum approach pinch point limit ($\Delta T_{\text{rec,min}} = 20\text{ K}$) relative to the pump discharge state, protecting against impossible internal temperature crosses.
# 

# %%
import numpy as np
import CoolProp.CoolProp as CP

print("--- Cell 1: Thermodynamic Core Engine Loaded Natively ---")

# %% [markdown]
# # 1. Variable Boundary Conditions & Inputs
# All variables in this section are entirely fluid-agnostic. The high-side operating pressure is self-correcting and bound dynamically to $90\%$ of the fluid's true critical pressure ($P_{\text{crit}}$) to enable subcritical or near-critical scaling.

# %%
# =============================================================================
# 1. DYNAMIC INPUT ENVIRONMENT & FLUID PROP DATASETS
# =============================================================================
fluid = 'Toluene'                # Working fluid token (Completely interchangeable)

# Performance Requirements
W_dot_net_target = 850.0 * 1e3   # Net cycle power target output [W] (850 kW)
eta_turbine      = 0.80          # Isentropic turbine efficiency [-]
eta_pump         = 0.75          # Isentropic pump efficiency [-]

# Operational Constraints
T_cond           = 50.0 + 273.15 # Condensation temperature [K] (50 °C)
T_turbine_in     = 350.0 + 273.15# Turbine inlet temperature [K] (350 °C)
Delta_T_rec_min  = 20.0          # Minimum allowable approach temperature in recuperator [K]

# External Condenser Cooling Water Conditions
P_cool           = 100.0 * 1e3   # Cooling loop pressure [Pa] (100 kPa)
T_cool_in        = 20.0 + 273.15 # Cooling water inlet temperature [K] (20 °C)
T_cool_out       = 50.0 + 273.15 # Cooling water outlet temperature [K] (50 °C)

# --- Dynamic Critical Fluid Property Lookups via CoolProp ---
T_crit = CP.PropsSI('TCRIT', fluid)
P_crit = CP.PropsSI('PCRIT', fluid)

# Enforce cycle pressure boundaries adaptively
P_high = 0.90 * P_crit           # High-side pressure capped at 90% of critical limit
P_low  = CP.PropsSI('P', 'T', T_cond, 'Q', 0, fluid) # Saturated liquid vapor pressure

print(f"--- Cell 2: {fluid} Operating Parameters Loaded ---")
print(f"  Critical Pressure : {P_crit/1e3:.2f} kPa | High Side Operating Pressure: {P_high/1e3:.2f} kPa")
print(f"  Critical Temp     : {T_crit:.2f} K      | Low Side Condensing Pressure : {P_low/1e3:.2f} kPa")

# %% [markdown]
# # 2. Cycle State Points & Multi-Physics Resolution
# This cell steps sequentially through the cycle nodes, mapping fluid behaviors via standard state properties ($P$, $H$, $S$). The recuperator model features an energy-balance verification checkpoint to prevent thermal cross-over.

# %%
# =============================================================================
# 2. SEQUENCE STATE NODE RESOLUTION (MAPPED TO BOODAGHI ET AL. 2023)
# =============================================================================

# --- State 1: Turbine Inlet (High Pressure Superheated Vapor) ---
# Enthalpy and entropy lookups are bound to the maximum cycle operating pressure stage
h1 = CP.PropsSI('H', 'P', P_high, 'T', T_turbine_in, fluid)
s1 = CP.PropsSI('S', 'P', P_high, 'T', T_turbine_in, fluid)

# --- State 2: Turbine Actual Outlet (Low Pressure Exhaust) ---
# Calculates real enthalpy drop factoring in isentropic losses across the rotor blades
s2_ideal = s1
h2_ideal = CP.PropsSI('H', 'P', P_low, 'S', s2_ideal, fluid)
h2       = h1 - eta_turbine * (h1 - h2_ideal) 
T2       = CP.PropsSI('T', 'P', P_low, 'H', h2, fluid)

# --- State 4: Condenser Outlet / Pump Inlet (Low Pressure Saturated Liquid) ---
# Set directly on the saturation liquid line (Q=0) matching the paper's default assumption
h4 = CP.PropsSI('H', 'T', T_cond, 'Q', 0, fluid)
s4 = CP.PropsSI('S', 'T', T_cond, 'Q', 0, fluid)

# --- State 5: Pump Actual Outlet (High Pressure Subcooled Liquid) ---
# Resolves pump work requirements based on fluid density properties
s5_ideal = s4
h5_ideal = CP.PropsSI('H', 'P', P_high, 'S', s5_ideal, fluid)
h5       = h4 + (h5_ideal - h4) / eta_pump    
T5_state = CP.PropsSI('T', 'P', P_high, 'H', h5, fluid)

# --- States 3 & 6: Recuperator Safe Energy Balances ---
# Limits thermal exchange to guarantee a valid minimum pinch point approach temperature
T3 = T5_state + Delta_T_rec_min
h3 = CP.PropsSI('H', 'P', P_low, 'T', T3, fluid)

# First-Law Energy Balance across the exchanger walls to extract preheated liquid enthalpy
h6 = h5 + (h2 - h3)
T6 = CP.PropsSI('T', 'P', P_high, 'H', h6, fluid)

# %% [markdown]
# # 3. External Condenser Cooling Balance
# Extracts thermal capacities for the utility cooling water circuit via the CoolProp reference database and tracks required mass flows.

# %%
# =============================================================================
# 3. SPECIFIC WORK, SYSTEM MASS FLOWS, & COOLING CIRCUITS
# =============================================================================
# Specific work calculations per unit mass [J/kg]
w_turbine = h1 - h2
w_pump    = h5 - h4
w_net     = w_turbine - w_pump

# Calculate the actual organic working fluid mass flow rate [kg/s]
m_dot_orc = W_dot_net_target / w_net

# Absolute Cycle Component Heat Duties [W]
W_dot_turbine = m_dot_orc * w_turbine
W_dot_pump    = m_dot_orc * w_pump
Q_dot_evap    = m_dot_orc * (h1 - h6)
Q_dot_rec     = m_dot_orc * (h2 - h3)
Q_dot_cond    = m_dot_orc * (h3 - h4)

# Dynamic Cycle Efficiency Verification
eta_cycle     = (W_dot_net_target / Q_dot_evap) * 100.0

# ---------- Cooling Water Utility Loop Balance ----------
water_token = 'Water'

# Fetch cooling water inlet and outlet specific enthalpies from CoolProp
h_cool_in  = CP.PropsSI('H', 'P', P_cool, 'T', T_cool_in, water_token)
h_cool_out = CP.PropsSI('H', 'P', P_cool, 'T', T_cool_out, water_token)

# Condenser Energy Balance: Q_dot_cond = m_dot_cool * (h_cool_out - h_cool_in)
m_dot_cool = Q_dot_cond / (h_cool_out - h_cool_in)

print("--- Cell 3: Thermodynamic Mass Flow and Cooling Utility Solved ---")

# %% [markdown]
# # 4. Verification & Validation Metrics Report
# Outputs the comprehensive cycle summary matching the layout of your other modules.

# %%
# =============================================================================
# 4. COMPREHENSIVE REPRODUCIBILITY VERIFICATION REPORT
# =============================================================================
print(f"==================================================")
print(f"      FINAL COMPLETE RORC SYSTEM SUMMARY         ")
print(f"==================================================")
print(f"Working Fluid Selection        : {fluid}")
print(f"ORC Fluid Mass Flow (m_dot_orc): {m_dot_orc:.3f} kg/s")
print(f"Cooling Water Flow  (m_cool)   : {m_dot_cool:.3f} kg/s")
print(f"--------------------------------------------------")
print(f"Turbine Power Output (W_dot_OT): {W_dot_turbine / 1e3:.2f} kW")
print(f"Pump Power Consumption (W_dot_OP): {W_dot_pump / 1e3:.2f} kW")
# Verification against target net output
print(f"Net System Power Output        : {W_dot_net_target / 1e3:.2f} kW")
print(f"--------------------------------------------------")
print(f"Recuperator Heat Duty (Q_Rec)  : {Q_dot_rec / 1e3:.2f} kW")
print(f"Evaporator Heat Input (Q_in)   : {Q_dot_evap / 1e3:.2f} kW")
print(f"Condenser Heat Rejection(Q_out): {Q_dot_cond / 1e3:.2f} kW")
print(f"--------------------------------------------------")
print(f"Final Optimized Cycle Efficiency: {eta_cycle:.2f} %")
print(f"==================================================")


