import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Define the parameters
π = 0.01
δ = 0.005
β = 0.009
α = 0.003
γ = 0.002
ρ = 0.01

# Define the initial conditions
u0 = [1000, 10, 0, 0]

# Define the differential equations
def zombie_outbreak(t, u):
    H, I, Z, R = u
    du = np.zeros(4)
    
    du[0] = π - β*H*Z - δ*H
    du[1] = β*H*Z - (ρ + δ)*I
    du[2] = ρ*I + γ*R - α*H*Z
    du[3] = δ*(H + I) + α*H*Z - γ*R
    
    return du

# Define the time span
tspan = (0.0, 100.0)

# Solve the differential equations
sol = solve_ivp(zombie_outbreak, tspan, u0, method="RK45", args=(), dense_output=False)

# Plot the results
plt.plot(sol.t, sol.y[0,:], label="Humans")
plt.plot(sol.t, sol.y[2,:], label="Zombies")
plt.xlabel("Time")
plt

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np

# Define the model function
def zombie_outbreak(t, y, π, δ, β, α, γ, ρ):
    H, I, Z, Rh, Rz = y

    # Calculate the derivatives of each state variable
    dHdt = π - β * H * Z - δ * H
    dIdt = β * H * Z - (ρ + δ) * I
    dZdt = ρ * I + γ * Rh - α * H * Z
    dRhdt = δ * (H + I) - γ * Rh
    dRzdt = α * H * Z
    
    return [dHdt, dIdt, dZdt, dRhdt, dRzdt]

# Set the initial conditions and parameters
y0 = [500, 1, 1, 0, 0]
π, δ, β, α, γ, ρ = 0, 0, 0.005, 0.001, 0.002, 1
params = (π, δ, β, α, γ, ρ)

# Set the time span for the simulation
t_span = (0.0, 500.0)

# Solve the differential equations and store the results
sol = solve_ivp(zombie_outbreak, t_span, y0, args=params, method='Radau', rtol=1e-8, atol=1e-8)

# Plot the results
plt.plot(sol.t, sol.y[0], label="Healthy People")
plt.plot(sol.t, sol.y[1], label="Infected People")
plt.plot(sol.t, sol.y[2], label="Zombies")
plt.plot(sol.t, sol.y[3], label="Removed People")
plt.plot(sol.t, sol.y[4], label="Removed Zombies")
plt.xlabel('Time')
plt.ylabel('Population')
plt.legend()
plt.show()

# Solve the differential equations and store the results
sol = solve_ivp(zombie_outbreak, t_span, y0, args=params, method='Radau', rtol=1e-8, atol=1e-8)

# Print the values of H and Z at each time step
for i in range(len(sol.t)):
    print("t = ", sol.t[i], ", H = ", sol.y[0, i], ", Z = ", sol.y[2, i])
