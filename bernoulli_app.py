import streamlit as st
from sympy import symbols, Eq, solve, log

# Title of the app
st.title("Fluid Mechanics Problem Solver")
st.write("""
This app solves various fluid mechanics problems, including:
- Bernoulli's Equation
- Power Calculation
- Reynolds Number and Moody's Friction Factor
- Major and Iterative Friction Loss
- Flow in Porous Media (Darcy's Law, Forchheimer Equation, Ergun Equation, etc.)
""")

# Sidebar for user inputs
st.sidebar.header("Input Parameters")
st.sidebar.write("Enter known values and leave one field blank for the unknown.")

# Input fields for Bernoulli's equation
st.sidebar.write("### Bernoulli's Equation")
P1 = st.sidebar.number_input("Pressure at point 1 (P1, Pa)", value=None, format="%f", key="P1")
P2 = st.sidebar.number_input("Pressure at point 2 (P2, Pa)", value=None, format="%f", key="P2")
v1 = st.sidebar.number_input("Velocity at point 1 (v1, m/s)", value=None, format="%f", key="v1")
v2 = st.sidebar.number_input("Velocity at point 2 (v2, m/s)", value=None, format="%f", key="v2")
h1 = st.sidebar.number_input("Height at point 1 (h1, m)", value=None, format="%f", key="h1")
h2 = st.sidebar.number_input("Height at point 2 (h2, m)", value=None, format="%f", key="h2")
rho = st.sidebar.number_input("Fluid density (rho, kg/m³)", value=None, format="%f", key="rho")
g = st.sidebar.number_input("Gravitational acceleration (g, m/s²)", value=9.81, format="%f", key="g")
hl_meters = st.sidebar.number_input("Head loss (hl, m)", value=None, format="%f", help="Head loss in meters of fluid.", key="hl_meters")
Wp_meters = st.sidebar.number_input("Pump head (Wp, m)", value=None, format="%f", help="Pump head in meters of fluid.", key="Wp_meters")
Wt_meters = st.sidebar.number_input("Turbine head (Wt, m)", value=None, format="%f", help="Turbine head in meters of fluid.", key="Wt_meters")

# Function to solve Bernoulli's equation
def solve_bernoulli(P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters, g):
    unknowns = [P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters]
    missing_vars = [i for i, val in enumerate(unknowns) if val is None]
    if len(missing_vars) != 1:
        return None, None, "Exactly one variable must be left blank (unknown)."
    
    missing_index = missing_vars[0]
    missing_var = symbols('x')
    unknowns[missing_index] = missing_var

    P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters = unknowns

    hl = rho * g * hl_meters if hl_meters is not None else None
    Wp = rho * g * Wp_meters if Wp_meters is not None else None
    Wt = rho * g * Wt_meters if Wt_meters is not None else None

    equation = Eq(P1 + 0.5 * rho * v1**2 + rho * g * h1 + (Wp if Wp else 0), 
                  P2 + 0.5 * rho * v2**2 + rho * g * h2 + (hl if hl else 0) + (Wt if Wt else 0))

    solution = solve(equation, missing_var)
    valid_solution = [sol.evalf() for sol in solution if sol.is_real]

    if not valid_solution:
        return None, None, "No physically meaningful solution found."
    else:
        units = ["Pa", "Pa", "m/s", "m/s", "m", "m", "kg/m³", "m", "m", "m"]
        unit = units[missing_index]
        value = float(valid_solution[0])
        if value < 0:
            return value, unit, f"The calculated value is **{round(value, 7)} {unit}**, which is physically impossible."
        else:
            return value, unit, f"The solved value is: **{round(value, 7)} {unit}**"

# Solve button for Bernoulli's equation
if st.sidebar.button("Solve Bernoulli's Equation"):
    if rho is None:
        st.error("Density (rho) must be provided.")
    else:
        value, unit, message = solve_bernoulli(P1, P2, v1, v2, h1, h2, rho, hl_meters, Wp_meters, Wt_meters, g)
        if value is None:
            st.error(message)
        elif value < 0:
            st.warning(message)
        else:
            st.success(message)

# New section for Power Calculation
st.header("Power Calculation from Heads")
st.write("### Power from Pump or Turbine Head")
flow_rate = st.number_input("Flow rate (Q, m³/s)", value=None, format="%f", help="Volumetric flow rate of the fluid.", key="flow_rate")
head = st.number_input("Head (m)", value=None, format="%f", help="Pump head or turbine head in meters.", key="head")

# Calculate Power
if st.button("Calculate Power"):
    if rho is None or g is None or flow_rate is None or head is None:
        st.error("Density (rho), gravitational acceleration (g), flow rate (Q), and head must be provided.")
    else:
        power = rho * g * head * flow_rate
        st.success(f"The power is: **{round(float(power), 7)} W**")

# New section for Reynolds Number and Moody's Friction Factor
st.header("Reynolds Number and Moody's Friction Factor Calculator")
st.write("### Reynolds Number (Re)")
velocity = st.number_input("Fluid velocity (v, m/s)", value=1.0, format="%f", help="Velocity of the fluid in the pipe.", key="velocity")
diameter = st.number_input("Pipe diameter (D, m)", value=0.1, format="%f", help="Inner diameter of the pipe.", key="diameter")
viscosity = st.number_input("Dynamic viscosity (μ, Pa·s)", value=0.001, format="%f", help="Measure of fluid resistance to flow.", key="viscosity")

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
roughness = st.number_input("Pipe roughness (ε, m)", value=0.0001, format="%f", help="Roughness of the pipe's inner surface.", key="roughness")

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
st.write("### Major Friction Loss (Darcy-Weisbach Equation)")
pipe_length = st.number_input("Pipe length (L, m)", value=None, format="%f", help="Length of the pipe.", key="pipe_length")
pipe_diameter = st.number_input("Pipe diameter (D, m)", value=None, format="%f", help="Inner diameter of the pipe.", key="pipe_diameter")
pipe_roughness = st.number_input("Pipe roughness (ε, m)", value=None, format="%f", help="Roughness of the pipe's inner surface.", key="pipe_roughness")
flow_velocity = st.number_input("Flow velocity (v, m/s)", value=None, format="%f", help="Velocity of the fluid in the pipe.", key="flow_velocity")

# Calculate Major Friction Loss
if st.button("Calculate Major Friction Loss"):
    if rho is None or pipe_length is None or pipe_diameter is None or pipe_roughness is None or flow_velocity is None:
        st.error("All inputs (rho, L, D, ε, v) must be provided.")
    else:
        Re = (rho * flow_velocity * pipe_diameter) / viscosity
        if Re < 2000:  # Laminar flow
            f = 64 / Re
        elif Re >= 4000:  # Turbulent flow (Colebrook-White equation)
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
            h_f = f * (pipe_length / pipe_diameter) * (flow_velocity**2 / (2 * g))
            st.success(f"The major friction loss is: **{round(float(h_f), 7)} m**")
        else:
            st.error("Cannot calculate friction loss for transitional flow.")

# New section for Iterative Friction Loss Calculation
st.header("Iterative Friction Loss Calculator")
st.write("### Iterative Friction Loss (Unknown Flow Velocity)")
pipe_length_iter = st.number_input("Pipe length (L, m)", value=None, format="%f", help="Length of the pipe for iterative calculation.", key="pipe_length_iter")
pipe_diameter_iter = st.number_input("Pipe diameter (D, m)", value=None, format="%f", help="Inner diameter of the pipe for iterative calculation.", key="pipe_diameter_iter")
pipe_roughness_iter = st.number_input("Pipe roughness (ε, m)", value=None, format="%f", help="Roughness of the pipe's inner surface for iterative calculation.", key="pipe_roughness_iter")
flow_rate_iter = st.number_input("Flow rate (Q, m³/s)", value=None, format="%f", help="Volumetric flow rate of the fluid for iterative calculation.", key="flow_rate_iter")

# Calculate Iterative Friction Loss
if st.button("Calculate Iterative Friction Loss"):
    if rho is None or pipe_length_iter is None or pipe_diameter_iter is None or pipe_roughness_iter is None or flow_rate_iter is None:
        st.error("All inputs (rho, L, D, ε, Q) must be provided.")
    else:
        v_guess = 1.0  # Initial guess for velocity (m/s)
        tolerance = 1e-6  # Convergence tolerance
        max_iterations = 1000  # Maximum number of iterations
        converged = False

        for i in range(max_iterations):
            Re = (rho * v_guess * pipe_diameter_iter) / viscosity
            if Re < 2000:  # Laminar flow
                f = 64 / Re
            elif Re >= 4000:  # Turbulent flow (Colebrook-White equation)
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
                h_f = f * (pipe_length_iter / pipe_diameter_iter) * (v_guess**2 / (2 * g))
                v_new = flow_rate_iter / (3.1416 * (pipe_diameter_iter / 2)**2)
                if abs(v_new - v_guess) < tolerance:
                    converged = True
                    break
                v_guess = v_new
            else:
                break

        if converged:
            st.success(f"The flow velocity is: **{round(float(v_guess), 7)} m/s**")
            st.success(f"The major friction loss is: **{round(float(h_f), 7)} m**")
        else:
            st.error("The solution did not converge. Check inputs or try a different initial guess.")

# New section for Porous Media Flow
st.header("Flow in Porous Media")
st.write("This section solves problems related to flow in porous media using:")
st.write("- Darcy's Law")
st.write("- Forchheimer Equation")
st.write("- Blake–Kozeny Equation")
st.write("- Ergun Equation")
st.write("- Burke–Plummer Equation")

# Input fields for porous media flow
st.sidebar.header("Porous Media Flow Parameters")
st.sidebar.write("Enter known values and leave one field blank for the unknown.")

# Common inputs
epsilon = st.sidebar.number_input("Porosity (ε)", value=None, format="%f", help="Porosity of the porous medium.", key="epsilon")
dp = st.sidebar.number_input("Particle diameter (dp, m)", value=None, format="%f", help="Diameter of particles in the porous medium.", key="dp")
k = st.sidebar.number_input("Permeability (k, m²)", value=None, format="%f", help="Permeability of the porous medium.", key="k")
beta = st.sidebar.number_input("Inertial resistance coefficient (β, m⁻¹)", value=None, format="%f", help="Coefficient for inertial resistance in Forchheimer Equation.", key="beta")
L_porous = st.sidebar.number_input("Length of porous medium (L, m)", value=None, format="%f", help="Length of the porous medium.", key="L_porous")
A_porous = st.sidebar.number_input("Cross-sectional area (A, m²)", value=None, format="%f", help="Cross-sectional area of the porous medium.", key="A_porous")

# Function to calculate Reynolds number for porous media
def calculate_reynolds_porous(rho, v, dp, mu, epsilon):
    return (rho * v * dp) / (mu * (1 - epsilon))

# Function to solve Darcy's Law
def solve_darcy(Q, k, A, mu, deltaP, L):
    if Q is None:
        Q = - (k * A * deltaP) / (mu * L)
        return Q, "Volumetric flow rate (Q)", f"The volumetric flow rate is: **{round(float(Q), 7)} m³/s**"
    elif k is None:
        k = - (Q * mu * L) / (A * deltaP)
        return k, "Permeability (k)", f"The permeability is: **{round(float(k), 7)} m²**"
    elif deltaP is None:
        deltaP = - (Q * mu * L) / (k * A)
        return deltaP, "Pressure drop (ΔP)", f"The pressure drop is: **{round(float(deltaP), 7)} Pa**"
    else:
        return None, None, "Exactly one variable must be left blank (unknown)."

# Function to solve Forchheimer Equation
def solve_forchheimer(deltaP, L, mu, k, beta, rho, v):
    if deltaP is None:
        deltaP = L * ( (mu / k) * v + beta * rho * v**2 )
        return deltaP, "Pressure drop (ΔP)", f"The pressure drop is: **{round(float(deltaP), 7)} Pa**"
    elif v is None:
        a = beta * rho
        b = mu / k
        c = -deltaP / L
        discriminant = b**2 - 4 * a * c
        if discriminant < 0:
            return None, None, "No real solution for velocity."
        v1 = (-b + discriminant**0.5) / (2 * a)
        v2 = (-b - discriminant**0.5) / (2 * a)
        v = max(v1, v2)  # Take the positive root
        return v, "Velocity (v)", f"The velocity is: **{round(float(v), 7)} m/s**"
    else:
        return None, None, "Exactly one variable must be left blank (unknown)."

# Function to solve Blake–Kozeny Equation
def solve_blake_kozeny(deltaP, L, mu, epsilon, dp, v):
    if deltaP is None:
        deltaP = 150 * ((1 - epsilon)**2 / epsilon**3) * (mu * v / dp**2) * L
        return deltaP, "Pressure drop (ΔP)", f"The pressure drop is: **{round(float(deltaP), 7)} Pa**"
    elif v is None:
        v = (deltaP / L) * (epsilon**3 / (150 * (1 - epsilon)**2)) * (dp**2 / mu)
        return v, "Velocity (v)", f"The velocity is: **{round(float(v), 7)} m/s**"
    else:
        return None, None, "Exactly one variable must be left blank (unknown)."

# Function to solve Ergun Equation
def solve_ergun(deltaP, L, mu, epsilon, dp, rho, v):
    if deltaP is None:
        deltaP = L * (150 * ((1 - epsilon)**2 / epsilon**3) * (mu * v / dp**2) + 1.75 * ((1 - epsilon) / epsilon**3) * (rho * v**2 / dp))
        return deltaP, "Pressure drop (ΔP)", f"The pressure drop is: **{round(float(deltaP), 7)} Pa**"
    elif v is None:
        return None, None, "Velocity (v) cannot be solved directly. Use an iterative approach."
    else:
        return None, None, "Exactly one variable must be left blank (unknown)."

# Function to solve Burke–Plummer Equation
def solve_burke_plummer(deltaP, L, epsilon, dp, rho, v):
    if deltaP is None:
        deltaP = L * 1.75 * ((1 - epsilon) / epsilon**3) * (rho * v**2 / dp)
        return deltaP, "Pressure drop (ΔP)", f"The pressure drop is: **{round(float(deltaP), 7)} Pa**"
    elif v is None:
        v = ((deltaP / L) * (epsilon**3 / (1.75 * (1 - epsilon))) * (dp / rho))**0.5
        return v, "Velocity (v)", f"The velocity is: **{round(float(v), 7)} m/s**"
    else:
        return None, None, "Exactly one variable must be left blank (unknown)."

# Automated solver for porous media flow
def automated_porous_media_solver(rho, v, dp, mu, epsilon, L, A, k, beta, deltaP=None, Q=None):
    Re_p = calculate_reynolds_porous(rho, v, dp, mu, epsilon)
    
    if Re_p < 1:
        st.write("**Flow Regime: Laminar (Re < 1)**")
        st.write("Using **Darcy's Law**.")
        return solve_darcy(Q, k, A, mu, deltaP, L)
    elif 1 <= Re_p <= 10:
        st.write("**Flow Regime: Transitional (1 ≤ Re ≤ 10)**")
        st.write("Using **Forchheimer Equation**.")
        return solve_forchheimer(deltaP, L, mu, k, beta, rho, v)
    elif Re_p > 10:
        st.write("**Flow Regime: Turbulent (Re > 10)**")
        st.write("Using **Ergun Equation**.")
        return solve_ergun(deltaP, L, mu, epsilon, dp, rho, v)

# Solve button for porous media flow
if st.button("Solve Porous Media Flow"):
    if rho is None or mu is None or epsilon is None or dp is None:
        st.error("Fluid density (ρ), viscosity (μ), porosity (ε), and particle diameter (dp) must be provided.")
    else:
        result, unit, message = automated_porous_media_solver(
            rho, v, dp, mu, epsilon, L_porous, A_porous, k, beta, deltaP, Q
        )
        if result is None:
            st.error(message)
        else:
            st.success(message)

# Instructions
st.write("### Instructions:")
st.write("""
1. Fill in all known values.
2. Leave one field blank for the unknown variable.
3. Click the 'Solve' button to calculate the unknown.
""")
