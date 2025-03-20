import streamlit as st
from sympy import symbols, Eq, solve

# Title of the app
st.title("Bernoulli's Equation Solver")
st.write("""
This app solves Bernoulli's equation for one unknown, including head loss, pump work, and turbine work.
""")

# Sidebar for user inputs
st.sidebar.header("Input Parameters")
st.sidebar.write("Enter known values and leave one field blank for the unknown.")

# Input fields
P1 = st.sidebar.number_input("Pressure at point 1 (P1, Pa)", value=None, format="%f")
P2 = st.sidebar.number_input("Pressure at point 2 (P2, Pa)", value=None, format="%f")
v1 = st.sidebar.number_input("Velocity at point 1 (v1, m/s)", value=None, format="%f")
v2 = st.sidebar.number_input("Velocity at point 2 (v2, m/s)", value=None, format="%f")
h1 = st.sidebar.number_input("Height at point 1 (h1, m)", value=None, format="%f")
h2 = st.sidebar.number_input("Height at point 2 (h2, m)", value=None, format="%f")
rho = st.sidebar.number_input("Fluid density (rho, kg/m³)", value=None, format="%f")
g = st.sidebar.number_input("Gravitational acceleration (g, m/s²)", value=9.81, format="%f")
hl = st.sidebar.number_input("Head loss (hl, J/kg)", value=0.0, format="%f")
Wp = st.sidebar.number_input("Pump work (Wp, J/kg)", value=0.0, format="%f")
Wt = st.sidebar.number_input("Turbine work (Wt, J/kg)", value=0.0, format="%f")

# Solve button
if st.sidebar.button("Solve"):
    # Validate inputs
    inputs = [P1, P2, v1, v2, h1, h2, rho, hl, Wp, Wt]
    if inputs.count(None) != 1:
        st.error("Exactly one variable must be left blank (unknown).")
    elif rho is None:
        st.error("Density (rho) must be provided.")
    else:
        # Define the unknown variable
        unknowns = [P1, P2, v1, v2, h1, h2, rho, hl, Wp, Wt]
        missing_vars = [i for i, val in enumerate(unknowns) if val is None]
        missing_index = missing_vars[0]
        missing_var = symbols('x')
        unknowns[missing_index] = missing_var

        # Assign variables
        P1, P2, v1, v2, h1, h2, rho, hl, Wp, Wt = unknowns

        # Define Bernoulli's equation
        equation = Eq(P1 + 0.5 * rho * v1**2 + rho * g * h1 + Wp, 
                      P2 + 0.5 * rho * v2**2 + rho * g * h2 + hl + Wt)

        # Solve for the unknown variable
        solution = solve(equation, missing_var)

        # Filter valid solutions
        valid_solution = [sol.evalf() for sol in solution if sol.is_real and sol >= 0]

        if not valid_solution:
            st.error("No physically meaningful solution found.")
        else:
            # Define unit based on missing variable
            units = ["Pa", "Pa", "m/s", "m/s", "m", "m", "kg/m³", "J/kg", "J/kg", "J/kg"]
            unit = units[missing_index]

            # Display the result
            st.success(f"The solved value is: **{round(float(valid_solution[0]), 7)} {unit}**")

# Instructions
st.write("### Instructions:")
st.write("""
1. Fill in all known values in the sidebar.
2. Leave one field blank for the unknown variable.
3. Click the 'Solve' button to calculate the unknown.
""")
