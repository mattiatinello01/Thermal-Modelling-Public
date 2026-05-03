"""
Scroll Compressor - Semi-empirical model (Winandy-Lebrun)
=========================================================
Inputs:  T_suc, P_suc, P_dis, frequency N
Outputs: mass flow rate, electric power, discharge temperature, COP

Reference: Winandy E., Saavedra C., Lebrun J. (2002)
           Simplified modelling of an open-type reciprocating compressor
           International Journal of Thermal Sciences
"""

import numpy as np
import CoolProp.CoolProp as CP

# ── Fluid ─────────────────────────────────────────────────────────────────────
FLUID = 'R290'  # Propane

# ── Compressor parameters (calibrated on manufacturer datasheet) ──────────────
V_sw   = 13.2e-6   # Swept volume [m³/rev]
N_nom  = 50.0      # Nominal frequency [Hz]
eta_vol_nom = 0.85 # Nominal volumetric efficiency [-]
eta_is_nom  = 0.72 # Nominal isentropic efficiency [-]
W_mec  = 50.0      # Mechanical losses [W]
Q_amb  = 30.0      # Heat loss to ambient [W]
dT_sh  = 10.0      # Suction superheating [K]

# ── Core functions ─────────────────────────────────────────────────────────────

def get_state(fluid, T, P):
    """Return enthalpy and specific volume at given T [K] and P [Pa]."""
    h = CP.PropsSI('H', 'T', T, 'P', P, fluid)
    v = CP.PropsSI('D', 'T', T, 'P', P, fluid)  # density
    v = 1 / v                                     # specific volume [m³/kg]
    s = CP.PropsSI('S', 'T', T, 'P', P, fluid)
    return h, v, s


def compressor(T_suc, P_suc, P_dis, N=N_nom, fluid=FLUID):
    """
    Semi-empirical scroll compressor model.

    Parameters
    ----------
    T_suc : float  Suction temperature [K]
    P_suc : float  Suction pressure [Pa]
    P_dis : float  Discharge pressure [Pa]
    N     : float  Frequency [Hz]

    Returns
    -------
    dict with mass flow, electric power, discharge temperature, COP
    """

    # 1. Suction state (with superheating correction)
    T_suc_real = T_suc + dT_sh
    h_suc, v_suc, s_suc = get_state(fluid, T_suc_real, P_suc)

    # 2. Isentropic discharge state
    h_dis_is = CP.PropsSI('H', 'S', s_suc, 'P', P_dis, fluid)

    # 3. Mass flow rate (volumetric model)
    m_dot = (V_sw * N) / v_suc * eta_vol_nom

    # 4. Electric power
    W_is  = m_dot * (h_dis_is - h_suc)   # isentropic work [W]
    W_el  = W_is / eta_is_nom + W_mec    # electric power [W]

    # 5. Real discharge enthalpy (energy balance)
    h_dis_real = h_suc + (W_el - W_mec - Q_amb) / m_dot

    # 6. Discharge temperature
    T_dis = CP.PropsSI('T', 'H', h_dis_real, 'P', P_dis, fluid)

    # 7. COP (compressor only)
    Q_evap = m_dot * (h_suc - CP.PropsSI('H', 'Q', 0, 'P', P_suc, fluid))
    COP = Q_evap / W_el

    return {
        'm_dot'  : m_dot,
        'W_el'   : W_el,
        'T_dis'  : T_dis - 273.15,   # °C
        'COP'    : COP,
        'W_is'   : W_is,
        'eta_is' : W_is / (W_el - W_mec),
    }


# ── Quick test ─────────────────────────────────────────────────────────────────
if __name__ == '__main__':

    # Typical R290 chiller operating point
    T_evap = -5.0 + 273.15   # K
    T_cond = 40.0 + 273.15   # K

    P_suc = CP.PropsSI('P', 'T', T_evap, 'Q', 1, FLUID)
    P_dis = CP.PropsSI('P', 'T', T_cond, 'Q', 0, FLUID)

    res = compressor(T_suc=T_evap, P_suc=P_suc, P_dis=P_dis)

    print(f"\n{'='*40}")
    print(f"  Scroll Compressor - R290")
    print(f"{'='*40}")
    print(f"  Mass flow rate  : {res['m_dot']*3600:.2f} kg/h")
    print(f"  Electric power  : {res['W_el']:.1f} W")
    print(f"  Discharge temp  : {res['T_dis']:.1f} °C")
    print(f"  Isentropic eff  : {res['eta_is']:.3f}")
    print(f"  COP             : {res['COP']:.2f}")
    print(f"{'='*40}\n")