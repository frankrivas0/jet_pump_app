import streamlit as st
import numpy as np
import pandas as pd
from utils.jet import jet_pump
from utils.jet import ipr
import matplotlib.pyplot as plt


st.set_page_config(page_title="Jet Pump Calculator", layout="wide")

well_data = False
range_data = False
simulation_completed = False

with st.form("well_data_form"):
    col1, col2 = st.columns(2)
    with col1:
        st.header("Well Data")
        well_name = st.text_input("Well name:")
        h_perforations = st.number_input("Perforations depth (ft):", format="%.2f")
        total_depth = st.number_input("Pump depth (ft):", format="%.2f")
        dtbgID = st.number_input("Tubing ID (in):", format="%.2f")
        dtbgOD = st.number_input("Tubing OD (in):", format="%.2f")
        dcsg = st.number_input("Casing ID (in):", format="%.2f")
        p_inj = st.number_input("Injection pressure (psi):", format="%.2f")
        p_wellhead = st.number_input("Wellhead pressure (psi):", format="%.2f")
        t_wellhead = st.number_input("Wellhead temperature (°F):", format="%.2f")
        st.divider()

        st.header("Fluid Properties")
        api = st.number_input("Oil API gravity:", format="%.2f")
        pb = st.number_input("Oil Bubble point pressure (psi):", format="%.2f")
        yg = st.number_input("Gas specific gravity:", format="%.2f")
        yw = st.number_input("Water specific gravity:", format="%.2f")
        mu_inj = st.number_input("Injection fluid viscosity (cP):", format="%.2f")
    with col2:

        st.header("Production Data")
        q_prod = st.number_input("Production flow rate (STB/D):", format="%.2f")
        bsw = st.number_input("BSW (%):", format="%.2f")/100
        GOR = st.number_input("Gas-Oil Ratio (SCF/STB):", format="%.2f")
        st.divider()

        st.header("Reservoir Data")
        p_res = st.number_input("Reservoir pressure (psi):", format="%.2f")
        t_bottom = st.number_input("Bottomhole temperature (°F):", format="%.2f")
        st.divider()

        st.header("Test Data")
        q_test = st.number_input("Test flow rate (STB/D):", format="%.2f")
        pwf_test = st.number_input("Test bottomhole pressure (psi):", format="%.2f")
        st.divider()

        st.header("Jet Pump Properties")
        pump_model = st.text_input("Pump name:")
        aj = st.number_input("Nozzle area (in²):", format="%.5f")
        at = st.number_input("Throat area (in²):", format="%.5f")

    data = st.form_submit_button("Submit")

if data:
    if (
        well_name.strip() == "" or
        pump_model.strip() == "" or
        any(v == 0.0 for v in [
            h_perforations, total_depth, dtbgID, dtbgOD, dcsg,
            p_inj, p_wellhead, t_wellhead, api, pb, yg, yw, mu_inj,
            q_prod, bsw, GOR, p_res, t_bottom, q_test, pwf_test, aj, at
        ])
    ):
        st.error("Please fill in all fields with non-zero values.")
    else:
        st.success("All fields are filled!")
        well_data = True
        st.session_state['well_data'] = {'well_data': well_data}

with st.form("Bottomhole Pressure range"):
    st.header("Bottomhole Pressure Range for Simulation")
    pwfmin = st.number_input("Minimum Bottomhole Pressure (psi):", format="%.2f")
    pwfmax = st.number_input("Maximum Bottomhole Pressure (psi):", format="%.2f")
    range_data = st.form_submit_button("Set Range")

if range_data:
    if pwfmin <= 0 or pwfmax <= 0 or pwfmin >= pwfmax:
        st.error("Please enter a valid pressure range where min < max and both are greater than 0.")
    else:
        st.success("Pressure range set!")
        range_data = True
        st.session_state['range_data'] = {'range_data': range_data}

simulate = st.button("Run Simulation")
if simulate:
    if not (st.session_state["well_data"] and st.session_state["range_data"]):
        st.error("Please complete the forms before simulating.")
    else:
        pwf = np.linspace(pwfmin, pwfmax, 25)
        values = []
        for pwfi in pwf:
            q_prod, pwf, q_inj, p_discharge, p_nozzle, power, efficiency, acm, qcav = jet_pump(pwfi, q_prod, aj, at, p_inj, p_wellhead,
                                                                                            pb, t_bottom, t_wellhead, total_depth, dcsg,
                                                                                            dtbgID, dtbgOD, api, GOR, bsw, yg, yw, mu_inj)
            values.append([q_prod, pwf, q_inj, p_discharge, p_nozzle, power, efficiency, acm, qcav])

        values = pd.DataFrame(values, columns = ['Production bbl/d', 'Pressure@pwf', 'Injection bbl/d', 'Pressure@Discharge', 'Pressure@Nozzle',
                                            'HP', 'Eff %', 'Cavitation area sq-in', 'Cavitation flow bbl/d'])
        x = ipr(q_test, pwf_test, p_res, pb)
        pwf2 = np.round(np.linspace(0, p_res, 50))
        qo = []
        for i in pwf2:
            qo.append(x.voguel(i))

        pwf2 = pwf2 - 0.433*(141.5/(api+131.5))*(h_perforations - total_depth)
        simulation_completed = True
        st.success("Simulation completed!")


if simulation_completed:
    st.title("Simulation Results")
    st.dataframe(values)
    st.divider()
    st.title("IPR and VLP Curves")
    fig, ax = plt.subplots()
    ax.plot(qo, pwf2, color='red', label='IPR CURVE')
    ax.plot(values['Production bbl/d'], values['Pressure@pwf'],  color='green', label='VLP CURVE')
    ax.set_xlim(0)
    ax.set_ylim(0)
    ax.set_ylabel('Pwf, psi')
    ax.set_xlabel('Fluid Rate, bpd')
    ax.legend()
    from matplotlib.ticker import MultipleLocator
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(40))
    ax.grid(which='both')
    plt.title(f'Well: {well_name} - Jet Pump: {pump_model}')
    st.pyplot(fig)
   
    

    
            



    




