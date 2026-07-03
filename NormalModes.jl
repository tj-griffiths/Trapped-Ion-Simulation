using IonSim
using QuantumOptics: timeevolution, expect
using Plots
using StatsPlots


# Build a two-ion chain with two identical Ca40 ions
ca1 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
ca2 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])

laser = Laser()

# Selected modes: mode 1 = COM, mode 2 = Stretch
chain  = LinearChain(
    ions = [ca1, ca2],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1,2]) # Modes along z-axis
)

chamber = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [laser])
wavelength_from_transition!(laser, ca1, ("S", "D"), chamber)
polarization!(laser, (x̂ - ẑ)/√2)
wavevector!(laser, (x̂ + ẑ)/√2)
intensity_from_pitime!(laser, 2e-6, ca1, ("S", "D"), chamber)

# Extract and validate normal mode frequencies
mode_com, mode_stretch = modes(chamber)

ν_com = frequency(mode_com) # COM mode frequency
ν_stretch = frequency(mode_stretch) # Stretch mode frequency


# Analytic predictions:
# ω_COM = ω_z (Coulumb shifts cancel for in-phase motion)
# ω_stretch = √3 * ω_z (Coulumb adds restoring force for out-of-phase motion)

ν_z = 1e6
analytic_com = ν_z
analytic_stretch = √3 * ν_z

println("=== Normal Mode Frequencies ===")
print("COM     Simulated: ", round(ν_com/1e6, digits = 6), " MHz")
print("COM     Analytic: ", round(analytic_com/1e6, digits = 6), " MHz")
print("Stretch Simulated: ", round(ν_stretch/1e6, digits = 6), " MHz")
print("Stretch Analytic: ", round(analytic_stretch/1e6, digits = 6), " MHz")
print("Stretch/COM ratio: ", round(ν_stretch/ν_com, digits = 6), " (Simulated)")


# Lamb-Dicke Parameters: Each ion × each mode
# For two equal-mass ions in a harmonic trap: the analytic mode vectors are:
# COM mode: b = [1/√2, 1/√2] (in-phase)
# Stretch mode: b = [1/√2, -1/√2] (out-of-phase)
# The stretch mode being opposite phase for ion 2 is needed for MS gate

η_ion1_com = lambdicke(mode_com, ca1, laser)
η_ion2_com = lambdicke(mode_com, ca2, laser)
η_ion1_stretch = lambdicke(mode_stretch, ca1, laser)
η_ion2_stretch = lambdicke(mode_stretch, ca2, laser)

println("\n=== Lamb-Dicke Parameters ===")
println("η(ion1, COM)):   ", round(η_ion1_com, digits = 6))
println("η(ion2, COM)):   ", round(η_ion2_com, digits = 6))
println("η(ion1, Stretch): ", round(η_ion1_stretch, digits = 6))
println("η(ion2, Stretch): ", round(η_ion2_stretch, digits = 6))

# Analytic η_0 for the COM mode:
#η = k_laser * cos(θ) * √(ħ/(2 m ω_COM))

println("\nRatio check (expect ≈ 1):")
println("η(ion1,COM)/η(ion2,COM) = ", round(η_ion1_com/η_ion2_com, digits = 4))
println("η(ion1,Stretch)/-η(ion2,Stretch) = ", round(η_ion1_stretch/-η_ion2_stretch, digits = 4))
println("Stretch η / COM η scaling:   ", round(abs(η_ion1_stretch)/abs(η_ion1_com), digits=4),
        " (expect ν_z/ν_stretch scaling ~ ", round(√(ν_z/analytic_stretch), digits=4), ")")



# Visualizing Mode Structure with Bar Chart

# To see how much each ion participates in each mode, we extract the normalized eigenvector component b_{ik} by dividing out η_0
η0 = abs(η_ion1_com) #Use ion1 COM as our reference η_0

b_com = [η_ion1_com, η_ion2_com] ./ η0 
b_stretch = [η_ion1_stretch, η_ion2_stretch] ./ η0

p_mode = groupedbar(
    ["Ion 1", "Ion 2"],
    hcat(b_com, b_stretch),   # each column is one mode
    label=["COM (ω = ν_z)" "Stretch (ω = √3 ν_z)"],
    ylabel="Normalized eigenvector component bᵢₖ",
    title="Two-ion normal mode eigenvectors",
    color=[:steelblue :crimson],
    lw=0,
    ylims=(-1.3, 1.3),
    legend=:topright,
    size=(500, 380)
)

hline!(p_mode, [0], lw=1, ls=:dash, color=:black, label = "")

# Sideband Rabi Flopping for COM and Stretch Modes

# Initial state: both ions in |S⟩ and mode in ground state |n=0⟩
# In IonSim, the full state is |ψ⟩ = |ion1⟩ ⊗ |ion2⟩ ⊗ |modes⟩
ψ0 = ca1["S"] ⊗ ca2["S"] ⊗ mode_com[0] ⊗ mode_stretch[0] # Both modes in ground state

T_com = 2π / (abs(η_ion1_com) * π/2) # Time for a full Rabi oscillation on COM mode
tspan_com = 0:0.1:2*T_com # Time span for simulation

# COM Blue Sideband
detuning!(laser, ν_com)
h_com = hamiltonian(chamber, timescale = 1e-6, rwa_cutoff = 1e5)
tout_com, sol_com = timeevolution.schroedinger_dynamic(tspan_com, ψ0, h_com)

ex_ion1_com = real.(expect(ionprojector(chamber, "D", 1), sol_com)) # Expectation value of the D state population for ion 1
ex_ion2_com = real.(expect(ionprojector(chamber, "D", 2), sol_com)) # Expectation value of the D state population for ion 2 

# Stretch Blue Sideband
T_stretch = 2π / (abs(η_ion1_stretch) * π/2) # Time for a full Rabi oscillation on Stretch mode
tspan_stretch = 0:0.1:2*T_stretch # Time span for simulation

detuning!(laser, ν_stretch)
h_stretch = hamiltonian(chamber, timescale = 1e-6, rwa_cutoff = 1e5)
tout_stretch, sol_stretch = timeevolution.schroedinger_dynamic(tspan_stretch, ψ0, h_stretch)

ex_ion1_stretch = real.(expect(ionprojector(chamber, "D", 1), sol_stretch))
ex_ion2_stretch = real.(expect(ionprojector(chamber, "D", 2), sol_stretch))

p_com = plot(tout_com, ex_ion1_com, label = "Ion 1", color = :blue, lw = 2.5)
plot!(p_com, tout_com, ex_ion2_com, label = "Ion 2", color = :red, lw = 2.5)
xlabel!(p_com, "Time (μs)")
ylabel!(p_com, "Population in |D⟩")
title!(p_com, "COM Mode Rabi Flopping")
ylims!(p_com, 0, 1)

p_stretch = plot(tout_stretch, ex_ion1_stretch, label = "Ion 1", color = :blue, lw = 2.5)
plot!(p_stretch, tout_stretch, ex_ion2_stretch, label = "Ion 2", color = :red, lw = 2.5)
xlabel!(p_stretch, "Time (μs)")
ylabel!(p_stretch, "Population in |D⟩")
title!(p_stretch, "Stretch Mode Rabi Flopping")
ylims!(p_stretch, 0, 1)


p_dynamics = plot(p_com, p_stretch, layout = (1,2), size=(1000,400), left_margin = 8mm, bottom_margin = 8mm, top_margin = 5mm, legend=:topright)

savefig(p_mode, "IonTrap/Results/TwoIon_NormalMode_Eigenvectors.png")
savefig(p_dynamics, "IonTrap/Results/TwoIon_NormalMode_RabiFlopping.png")
println("\n Plots Saved")

