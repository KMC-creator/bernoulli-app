# New section for Major Friction Loss Calculation
st.header("Major Friction Loss Calculator")

# Input fields for Major Friction Loss
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
pipe_length_iter = st.number_input("Pipe length (L, m)", value=None, format="%f", help="Length of the pipe for iterative calculation.", key="pipe_length_iter")
pipe_diameter_iter = st.number_input("Pipe diameter (D, m)", value=None, format="%f", help="Inner diameter of the pipe for iterative calculation.", key="pipe_diameter_iter")
pipe_roughness_iter = st.number_input("Pipe roughness (ε, m)", value=None, format="%f", help="Roughness of the pipe's inner surface for iterative calculation.", key="pipe_roughness_iter")
flow_rate_iter = st.number_input("Flow rate (Q, m³/s)", value=None, format="%f", help="Volumetric flow rate of the fluid for iterative calculation.", key="flow_rate_iter")

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
