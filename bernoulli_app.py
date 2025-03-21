import streamlit as st
from sympy import symbols, Eq, solve, log

# Title of the app
st.title("Bernoulli's Equation Solver")
st.write("""
This app solves Bernoulli's equation for one unknown, including head loss, pump work, and turbine work.
It also calculates the Reynolds Number, Moody's Friction Factor, Power from heads, Major Friction Loss, and Iterative Friction Loss.
""")

# Sidebar for user inputs
st.sidebar.header("Input Parameters")
st.sidebar.write("Enter known values and leave one field blank for the unknown.")

# Input fields for Bernoulli's equation
P1 = st.sidebar.number_input("Pressure at point 1 (P1, Pa)", value=None, format="%f")
P2 = st.sidebar.number_input("Pressure at point 2 (P2, Pa)", value=None, format="%f")
v1 = st.sidebar.number_input("Velocity at point 1 (v1, m/s)", value=None, format="%f")
v2 = st.sidebar.number_input("Velocity at point 2 (v2, m/s)", value=None, format="%f")
h1 = st.sidebar.number_input("Height at point 1 (h1, m)", value=None, format="%f")
h2 = st.sidebar.number_input("Height at point 2 (h2, m)", value=None, format="%f")
rho = st.sidebar.number_input("Fluid density (rho, kg/m³)", value=None, format="%f")
g = st.sidebar.number_input("Gravitational acceleration (g, m/s²)", value=9.81, format="%f")
hl_meters = st.sidebar.number_input("Head loss (hl, m)", value=None, format="%f", help="Head loss in meters of fluid.")
Wp_meters = st.sidebar.number_input("Pump head (Wp, m)", value=None, format="%f", help="Pump head in meters of fluid.")
Wt_meters = st.sidebar.number_input("Turbine head (Wt, m)", value=None, format="%f", help="Turbine head in meters of fluid.")

# Solve button for Bernoulli's equation
if st.sidebar.button("Solve Bernoulli's Equation"):
    # Validate inputs
    inputs = [P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters]
    if inputs.count(None) != 1:
        st.error("Exactly one variable must be left blank (unknown).")
    elif rho is None:
        st.error("Density (rho) must be provided.")
    else:
        # Define the unknown variable
        unknowns = [P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters]
        missing_vars = [i for i, val in enumerate(unknowns) if val is None]
        missing_index = missing_vars[0]
        missing_var = symbols('x')
        unknowns[missing_index] = missing_var

        # Assign variables
        P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters = unknowns

        # Convert head loss, pump work, and turbine work from meters to pressure (Pa)
        hl = rho * g * hl_meters if hl_meters is not None else None
        Wp = rho * g * Wp_meters if Wp_meters is not None else None
        Wt = rho * g * Wt_meters if Wt_meters is not None else None

        # Define Bernoulli's equation
        equation = Eq(P1 + 0.5 * rho * v1**2 + rho * g * h1 + (Wp if Wp else 0), 
                      P2 + 0.5 * rho * v2**2 + rho * g * h2 + (hl if hl else 0) + (Wt if Wt else 0))

        # Solve for the unknown variable
        solution = solve(equation, missing_var)

        # Filter valid solutions
        valid_solution = [sol.evalf() for sol in solution if sol.is_real]

        if not valid_solution:
            st.error("No physically meaningful solution found.")
        else:
            # Define unit based on missing variable
            units = ["Pa", "Pa", "m/s", "m/s", "m", "m", "kg/m³", "m", "m", "m"]
            unit = units[missing_index]

            # Check if the solution is physically meaningful
            value = float(valid_solution[0])
            if value < 0:  # Non-negative check for most parameters
                st.warning(f"The calculated value is **{round(value, 7)} {unit}**, which is physically impossible.")
            else:
                st.success(f"The solved value is: **{round(value, 7)} {unit}**")

# New section for Power Calculation
st.header("Power Calculation from Heads")

# Input fields for Power Calculation
st.write("### Power from Pump or Turbine Head")
flow_rate = st.number_input("Flow rate (Q, m³/s)", value=None, format="%f", help="Volumetric flow rate of the fluid.")
head = st.number_input("Head (m)", value=None, format="%f", help="Pump head or turbine head in meters.")

# Calculate Power
if st.button("Calculate Power"):
    if rho is None or g is None or flow_rate is None or head is None:
        st.error("Density (rho), gravitational acceleration (g), flow rate (Q), and head must be provided.")
    else:
        power = rho * g * head * flow_rate
        st.success(f"The power is: **{round(float(power), 7)} W**")

# New section for Reynolds Number and Moody's Friction Factor
st.header("Reynolds Number and Moody's Friction Factor Calculator")

# Input fields for Reynolds Number and Moody's Friction Factor
st.write("### Reynolds Number (Re)")
velocity = st.number_input("Fluid velocity (v, m/s)", value=1.0, format="%f", help="Velocity of the fluid in the pipe.")
diameter = st.number_input("Pipe diameter (D, m)", value=0.1, format="%f", help="Inner diameter of the pipe.")
viscosity = st.number_input("Dynamic viscosity (μ, Pa·s)", value=0.001, format="%f", help="Measure of fluid resistance to flow.")

# Calculate Reynolds Number
if st.button("Calculate Reynolds Number"):
    if rho is None:
        st.error("Density (rho) must be provided.")
    elif velocity <= 0 or diameter <= 0 or viscosity <= 0:
        st.error("Velocity, diameter, and viscosity must be positive.")
    else:
        Re = (rho * velocity * diameter) / viscosity
        st.success(f"The Reynolds Number is: **{round(float(Re), 7)}**")

# Input fields for Moody's Friction Factor
st.write("### Moody's Friction Factor (f)")
roughness = st.number_input("Pipe roughness (ε, m)", value=0.0001, format="%f", help="Roughness of the pipe's inner surface.")

# Calculate Moody's Friction Factor
if st.button("Calculate Moody's Friction Factor"):
    if rho is None:
        st.error("Density (rho) must be provided.")
    elif velocity <= 0 or diameter <= 0 or viscosity <= 0 or roughness < 0:
        st.error("Velocity, diameter, viscosity, and roughness must be positive.")
    else:
        Re = (rho * velocity * diameter) / viscosity
        if Re < 2000:  # Laminar flow
            f = 64 / Re
            st.success(f"The flow is laminar. Friction factor (f) is: **{round(float(f), 7)}**")
        elif Re >= 4000:  # Turbulent flow (Colebrook-White equation)
            # Initial guess for f
            f = 0.02
            tolerance = 1e-6
            max_iterations = 1000
            for i in range(max_iterations):
                f_new = (-2 * log((roughness / diameter) / 3.7 + 2.51 / (Re * f**0.5)) / log(10))**-2
                if abs(f_new - f) < tolerance:
                    f = f_new
                    break
                f = f_new
            st.success(f"The flow is turbulent. Friction factor (f) is: **{round(float(f), 7)}**")
        else:  # Transitional flow
            st.warning("The flow is in the transitional region. Friction factor cannot be accurately calculated.")

# New section for Major Friction Loss Calculation
st.header("Major Friction Loss Calculator")

# Input fields for Major Friction Loss
st.write("### Major Friction Loss (Darcy-Weisbach Equation)")
pipe_length = st.number_input("Pipe length (L, m)", value=None, format="%f", help="Length of the pipe.")
pipe_diameter = st.number_input("Pipe diameter (D, m)", value=None, format="%f", help="Inner diameter of the pipe.")
pipe_roughness = st.number_input("Pipe roughness (ε, m)", value=None, format="%f", help="Roughness of the pipe's inner surface.")
flow_velocity = st.number_input("Flow velocity (v, m/s)", value=None, format="%f", help="Velocity of the fluid in the pipe.")

# Calculate Major Friction Loss
if st.button("Calculate Major Friction Loss"):
    if rho is None or pipe_length is None or pipe_diameter is None or pipe_roughness is None or flow_velocity is None:
        st.error("All inputs (rho, L, D, ε, v) must be provided.")
    else:
        # Calculate Reynolds Number
        Re = (rho * flow_velocity * pipe_diameter) / viscosity

        # Determine friction factor (f)
        if Re < 2000:  # Laminar flow
            f = 64 / Re
        elif Re >= 4000:  # Turbulent flow (Colebrook-White equation)
            # Initial guess for f
            f = 0.02
            tolerance = 1e-6
            max_iterations = 1000
            for i in range(max_iterations):
                f_new = (-2 * log((pipe_roughness / pipe_diameter) / 3.7 + 2.51 / (Re * f**0.5)) / log(10))**-2
                if abs(f_new - f) < tolerance:
                    f = f_new
                    break
                f = f_new
        else:  # Transitional flow
            st.warning("The flow is in the transitional region. Friction factor cannot be accurately calculated.")
            f = None

        if f is not None:
            # Calculate head loss due to friction
            h_f = f * (pipe_length / pipe_diameter) * (flow_velocity**2 / (2 * g))
            st.success(f"The major friction loss is: **{round(float(h_f), 7)} m**")
        else:
            st.error("Cannot calculate friction loss for transitional flow.")

# New section for Iterative Friction Loss Calculation
st.header("Iterative Friction Loss Calculator")

# Input fields for Iterative Friction Loss
st.write("### Iterative Friction Loss (Unknown Flow Velocity)")
pipe_length_iter = st.number_input("Pipe length (L, m)", value=None, format="%f", help="Length of the pipe.")
pipe_diameter_iter = st.number_input("Pipe diameter (D, m)", value=None, format="%f", help="Inner diameter of the pipe.")
pipe_roughness_iter = st.number_input("Pipe roughness (ε, m)", value=None, format="%f", help="Roughness of the pipe's inner surface.")
flow_rate_iter = st.number_input("Flow rate (Q, m³/s)", value=None, format="%f", help="Volumetric flow rate of the fluid.")

# Calculate Iterative Friction Loss
if st.button("Calculate Iterative Friction Loss"):
    if rho is None or pipe_length_iter is None or pipe_diameter_iter is None or pipe_roughness_iter is None or flow_rate_iter is None:
        st.error("All inputs (rho, L, D, ε, Q) must be provided.")
    else:
        # Initial guess for flow velocity
        v_guess = 1.0  # Initial guess for velocity (m/s)
        tolerance = 1e-6  # Convergence tolerance
        max_iterations = 1000  # Maximum number of iterations
        converged = False

        for i in range(max_iterations):
            # Calculate Reynolds Number
            Re = (rho * v_guess * pipe_diameter_iter) / viscosity

            # Determine friction factor (f)
            if Re < 2000:  # Laminar flow
                f = 64 / Re
            elif Re >= 4000:  # Turbulent flow (Colebrook-White equation)
                # Initial guess for f
                f = 0.02
                for j in range(max_iterations):
                    f_new = (-2 * log((pipe_roughness_iter / pipe_diameter_iter) / 3.7 + 2.51 / (Re * f**0.5)) / log(10))**-2
                    if abs(f_new - f) < tolerance:
                        f = f_new
                        break
                    f = f_new
            else:  # Transitional flow
                st.warning("The flow is in the transitional region. Friction factor cannot be accurately calculated.")
                f = None
                break

            if f is not None:
                # Calculate head loss due to friction
                h_f = f * (pipe_length_iter / pipe_diameter_iter) * (v_guess**2 / (2 * g))

                # Update flow velocity using the flow rate
                v_new = flow_rate_iter / (3.1416 * (pipe_diameter_iter / 2)**2)

                # Check for convergence
                if abs(v_new - v_guess) < tolerance:
                    converged = True
                    break

                # Update guess for velocity
                v_guess = v_new
            else:
                break

        if converged:
            st.success(f"The flow velocity is: **{round(float(v_guess), 7)} m/s**")
            st.success(f"The major friction loss is: **{round(float(h_f), 7)} m**")
        else:
            st.error("The solution did not converge. Check inputs or try a different initial guess.")

# Instructions
st.write("### Instructions:")
st.write("""
1. Fill in all known values in the sidebar.
2. Leave one field blank for the unknown variable.
3. Click the 'Solve' button to calculate the unknown.
4. Use the sections below to calculate Reynolds Number, Moody's Friction Factor, Power, Major Friction Loss, and Iterative Friction Loss.
""")
