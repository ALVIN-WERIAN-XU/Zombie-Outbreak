using DifferentialEquations, Plots

# Define the model parameters
params = (π=0.02, δ=0.01, β=0.05, α=0.03, γ=0.01, ρ=5.0)

# Define the differential equations
function zombie_eqs(du, u, p, t)
    H, I, Z, Rh, Rz = u
    π, δ, β, α, γ, ρ = p

    du[1] = π - β*H*Z - δ*H
    du[2] = β*H*Z - (ρ+δ)*I
    du[3] = ρ*I + γ*Rh - α*H*Z
    du[4] = δ*(H+I) - γ*Rh
    du[5] = α*H*Z
end

# Define the initial conditions
u0 = [500.0, 1.0, 1.0, 0.0, 0.0]

# Define the time span
tspan = (0.0, 100.0)

# Define the problem and solve it using the Tsit5 method
prob = ODEProblem(zombie_eqs, u0, tspan, params)
sol = solve(prob, Tsit5(), saveat=0.1)

# Plot the results
plot(sol.t, sol[1,:], label="Humans")
plot!(sol.t, sol[3,:], label="Zombies")
xlabel!("Time")
ylabel!("Population")
title!("Zombie Outbreak Model")

