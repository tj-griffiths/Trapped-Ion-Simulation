using IonSim
using QuantumOptics: timeevolution, expect
using Plots

# Same as Setup.jl
ca = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
laser = Laser()

chain = LinearChain(
    ions = [ca],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1])
)

chamber = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [laser])
wavelength_from_transition!(laser, ca, ("S", "D"), chamber)
polarization!(laser, (x̂ - ẑ)/√2)
wavevector!(laser, (x̂ + ẑ)/√2)

intensity_from_pitime!(laser, 2e-6, ca, ("S", "D"), chamber)
mode = modes(chamber)[1]
ν = frequency(mode) # Vibrational mode frequency set to 1 MHz from comfrequencies.ẑ

# The Lamb-Dicke Parameter: How strongly the ion/mode/laser geometry couples to internal and motional states
η = abs(lambdicke(mode, ca, laser)) 
println("Lamb-Dicke Parameter η = ", η)

# Bare Carrier Rabi Frequency:
Ω = π/2

T_sideband = 2π / (η * Ω) # Time for a full sideband Rabi oscillation
tspan = 0:0.05:T_sideband # Time span for simulation

# Blue and Red Sidebands - Starting in Motional Groundstate |n = 0⟩

# Blue sideband drives |S,0⟩ ↔ |D,1⟩ transition
detuning!(laser, ν) # Detune laser to blue sideband (Δ = +ν)
h_blue = hamiltonian(chamber, timescale = 1e-6)
ψ0_ground = ca["S"] ⊗ mode[0] # Initial state: ion in S state and mode in ground state
tout,sol_bluefrom0 = timeevolution.schroedinger_dynamic(tspan, ψ0_ground, h_blue)
ex_blue_from0 = real.(expect(ionprojector(chamber, "D"), sol_bluefrom0)) # Expectation value of the D state population (|D⟩⟨D|)

# Red Sideband drives |S,0⟩ ↔ |D,0⟩ transition
detuning!(laser, -ν) # Detune laser to red sideband (Δ = -ν)
h_red = hamiltonian(chamber, timescale = 1e-6)
tout, sol_redfrom0 = timeevolution.schroedinger_dynamic(tspan, ψ0_ground, h_red)
ex_red_from0 = real.(expect(ionprojector(chamber, "D"), sol_redfrom0)) # Expectation value of the D state population (|D⟩⟨D|)

# Red and Blue Sidebands - Starting from one phonon state |n = 1⟩

ψ0_excited = ca["S"] ⊗ mode[1] # Initial state: ion in S state and mode in first excited state

# Blue sideband drives |S,1⟩ ↔ |D,2⟩ transition
detuning!(laser, ν) # Detune laser to blue sideband (Δ = +ν)
h_blue = hamiltonian(chamber, timescale = 1e-6)
tout, sol_bluefrom1 = timeevolution.schroedinger_dynamic(tspan, ψ0_excited, h_blue)
ex_blue_from1 = real.(expect(ionprojector(chamber, "D"), sol_bluefrom1)) 

# Red sideband drives |S,1⟩ ↔ |D,0⟩ transition
detuning!(laser, -ν) # Detune laser to red sideband (Δ = -ν)
h_red = hamiltonian(chamber, timescale = 1e-6)
tout, sol_redfrom1 = timeevolution.schroedinger_dynamic(tspan, ψ0_excited, h_red)
ex_red_from1 = real.(expect(ionprojector(chamber, "D"), sol_redfrom1))

# Compare leading-order Lamb-Dicke approximation with full Hamiltonian simulation
Ω_eff_01 = Ω * η * exp(-η^2/2) * √1 # Effective Rabi frequency for |S,0⟩ ↔ |D,1⟩ transition
Ω_eff_12 = Ω * η * exp(-η^2/2) * √2 # Effective Rabi frequency for |S,1⟩ ↔ |D,2⟩ transition

analytic_01 = sin.(Ω_eff_01/2 .* tout).^2 # Analytical solution for |S,0⟩ ↔ |D,1⟩ transition
analytic_02 = sin.(Ω_eff_12/2 .* tout).^2 # Analytical solution for |S,1⟩ ↔ |D,2⟩ transition

# Plotting

p1 = plot(tout, ex_blue_from0, label = "blue, from n=0", color =:blue, lw = 2.5)
plot!(p1, tout, ex_red_from0, label = "red, from n=0", color =:red, lw = 2.5)
plot!(p1, tout, analytic_01, label = "analytic |0> ↔ |1>", color =:black, lw = 2.5, ls = :dash)
xlabel!(p1, "Time (μs)")
ylabel!(p1, "Population in |D⟩")
title!(p1, "Sideband Rabi Flopping from |n=0⟩")
ylims!(p1, 0, 1)

p2 = plot(tout, ex_blue_from1, label = "blue, from n=1", color =:red, lw = 2.5)
plot!(p2, tout, ex_red_from1, label = "red, from n=1", color =:blue, lw = 2.5)
plot!(p2, tout, analytic_02, label = "analytic |1> ↔ |2>", color =:black, lw = 2.5, ls = :dash)
xlabel!(p2, "Time (μs)")
ylabel!(p2, "Population in |D⟩")
title!(p2, "Sideband Rabi Flopping from |n=1⟩")
ylims!(p2, 0, 1)

plot(p1, p2, layout = (1, 2) , size = (1000, 420), left_margin = 8mm, bottom_margin = 8mm, top_margin = 5mm)
savefig("IonTrap/Results/Sideband_Rabi_Flopping.png")
println("Saved plot to Sideband_Rabi_Flopping.png")

# Looking at the plot we see that the blue sideband only reaches ~ 0.27 population and fails to match the analytic curve's amplitude.
# The difference in detuning and bare carrier Rabi frequency (1 MHz, 250 kHz respectively)  produces an AC Stark shift that is comparable to the sideband Rabi frequency (Ω_eff/2π  = 17 kHz). 
# This detunes the Rabi oscillation. To fix we drop the off-resonant carrier term from the Hamiltonian:

detuning!(laser, ν) 
h_blue = hamiltonian(chamber, timescale = 1e-6, rwa_cutoff = 1e5)
detuning!(laser, -ν)
h_red = hamiltonian(chamber, timescale = 1e-6, rwa_cutoff = 1e5)

tout, sol_bluefrom0 = timeevolution.schroedinger_dynamic(tspan, ψ0_ground, h_blue)
ex_blue_from0 = real.(expect(ionprojector(chamber, "D"), sol_bluefrom0)) 
tout, sol_redfrom0 = timeevolution.schroedinger_dynamic(tspan, ψ0_ground, h_red)
ex_red_from0 = real.(expect(ionprojector(chamber, "D"), sol_redfrom0))

tout, sol_bluefrom1 = timeevolution.schroedinger_dynamic(tspan, ψ0_excited, h_blue)
ex_blue_from1 = real.(expect(ionprojector(chamber, "D"), sol_bluefrom1))
tout, sol_redfrom1 = timeevolution.schroedinger_dynamic(tspan, ψ0_excited, h_red)
ex_red_from1 = real.(expect(ionprojector(chamber, "D"), sol_redfrom1))

p1 = plot(tout, ex_blue_from0, label = "blue, from n=0", color =:blue, lw = 2.5)
plot!(p1, tout, ex_red_from0, label = "red, from n=0", color =:red, lw = 2.5)
plot!(p1, tout, analytic_01, label = "analytic |0> ↔ |1>", color =:black, lw = 2.5, ls = :dash)
xlabel!(p1, "Time (μs)")
ylabel!(p1, "Population in |D⟩")
title!(p1, "Sideband Rabi Flopping from |n=0⟩ (RWA)")
ylims!(p1, 0, 1)

p2 = plot(tout, ex_blue_from1, label = "blue, from n=1", color =:blue, lw = 2.5)
plot!(p2, tout, ex_red_from1, label = "red, from n=1", color =:red, lw = 2.5)
plot!(p2, tout, analytic_02, label = "analytic |1> ↔ |2>", color =:black, lw = 2.5, ls = :dash)
xlabel!(p2, "Time (μs)")
ylabel!(p2, "Population in |D⟩")
title!(p2, "Sideband Rabi Flopping from |n=1⟩ (RWA)")
ylims!(p2, 0, 1)

plot(p1, p2, layout = (1, 2) , size = (1000, 420), left_margin = 8mm, bottom_margin = 8mm, top_margin = 5mm)
savefig("IonTrap/Results/Sideband_Rabi_Flopping_RWA.png")
