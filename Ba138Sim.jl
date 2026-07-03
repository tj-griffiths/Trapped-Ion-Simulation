# Ba138 Simulation - Carrier Rabi Flopping, Sidebands, and Ca40 Comparison

using IonSim
using QuantumOptics: timeevolution, expect
using Plots
using Plots.PlotMeasures

# Ba138 Carrier Rabi Flopping

ba = Ba138([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
laser_ba = Laser()

chain_ba = LinearChain(
    ions = [ba],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1])
)

chamber_ba = Chamber(iontrap=chain_ba, B=4e-4, Bhat=ẑ, δB=0, lasers=[laser_ba])

wavelength_from_transition!(laser_ba, ba, ("S", "D"), chamber_ba)
polarization!(laser_ba, (x̂ - ẑ)/√2)
wavevector!(laser_ba, (x̂ + ẑ)/√2)

# Use 5 μs pi-time for Ba138 as labs due to different matrix elements at 1762 nm
τ_π_ba = 5e-6
intensity_from_pitime!(laser_ba, τ_π_ba, ba, ("S", "D"), chamber_ba)

mode_ba = modes(chamber_ba)[1]
η_ba = abs(lambdicke(mode_ba, ba, laser_ba))
Ω_ba = π / (τ_π_ba * 1e6)   # rad/μs

println("=== Ba138 Parameters ===")
println("Qubit wavelength: ", round(laser_ba.λ * 1e9, digits=2), " nm")
println("Lamb-Dicke η:     ", round(η_ba, digits=5))
println("Carrier Ω:        ", round(Ω_ba, digits=4), " rad/μs  (π-time = ", τ_π_ba*1e6, " μs)")

h_ba_carrier = hamiltonian(chamber_ba, timescale=1e-6)
ψ0_ba = ba["S"] ⊗ mode_ba[0]

tspan_carrier = 0:0.1:15   # 3 full pi-pulses
tout_c, sol_c = timeevolution.schroedinger_dynamic(tspan_carrier, ψ0_ba, h_ba_carrier)

ex_ba_carrier = real.(expect(ionprojector(chamber_ba, "D"), sol_c))

# Analytic: P_D(t) = sin²(Ωt/2) 
analytic_ba = sin.(Ω_ba .* tout_c ./ 2).^2

# Ba138 Sideband Driving

ν_ba = frequency(mode_ba)
T_sb_ba = 2π / (η_ba * Ω_ba)   # expected sideband period (μs)
tspan_sb = 0:0.5:2*T_sb_ba

# Blue sideband from |n=0⟩ — expect slow oscillation due to small η
detuning!(laser_ba, ν_ba)
h_ba_blue = hamiltonian(chamber_ba, timescale=1e-6, rwa_cutoff=1e5)
tout_sb, sol_sb_blue = timeevolution.schroedinger_dynamic(tspan_sb, ψ0_ba, h_ba_blue)
ex_ba_blue = real.(expect(ionprojector(chamber_ba, "D"), sol_sb_blue))

# Red sideband from |η=0⟩ — expect near-zero (no phonon to remove)
detuning!(laser_ba, -ν_ba)
h_ba_red = hamiltonian(chamber_ba, timescale=1e-6, rwa_cutoff=1e5)
_, sol_sb_red = timeevolution.schroedinger_dynamic(tspan_sb, ψ0_ba, h_ba_red)
ex_ba_red = real.(expect(ionprojector(chamber_ba, "D"), sol_sb_red))

# Analytic Lab-Dicke sideband formula: Ω_eff = Ω × η × exp(-η²/2) × √1
Ω_eff_ba = Ω_ba * η_ba * exp(-η_ba^2/2)
analytic_ba_sb = sin.(Ω_eff_ba/2 .* tout_sb).^2

println("\n=== Ba138 Sideband Parameters ===")
println("Sideband Ω_eff: ", round(Ω_eff_ba, digits=4), " rad/μs")
println("Sideband period: ", round(T_sb_ba, digits=4), " μs")

# Ca40 Comparison

ca = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
laser_ca = Laser()
chain_ca = LinearChain(
    ions = [ca],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1])
)

chamber_ca = Chamber(iontrap=chain_ca, B=4e-4, Bhat=ẑ, δB=0, lasers=[laser_ca])
wavelength_from_transition!(laser_ca, ca, ("S", "D"), chamber_ca)
polarization!(laser_ca, (x̂ - ẑ)/√2)
wavevector!(laser_ca, (x̂ + ẑ)/√2)
intensity_from_pitime!(laser_ca, 5e-6, ca, ("S", "D"), chamber_ca)

mode_ca = modes(chamber_ca)[1]
η_ca = abs(lambdicke(mode_ca, ca, laser_ca))
Ω_ca = π / 5   # rad/μs

# Ca40 Blue Sideband
ν_ca = frequency(mode_ca)
T_sb_ca = 2π / (η_ca * Ω_ca)  
tspan_ca_sb = 0:0.5:2*T_sb_ca
ψ0_ca = ca["S"] ⊗ mode_ca[0]

detuning!(laser_ca, ν_ca)
h_ca_blue = hamiltonian(chamber_ca, timescale=1e-6, rwa_cutoff=1e5)
tout_ca_sb, sol_ca_sb_blue = timeevolution.schroedinger_dynamic(tspan_ca_sb, ψ0_ca, h_ca_blue)
ex_ca_blue = real.(expect(ionprojector(chamber_ca, "D"), sol_ca_sb_blue))
Ω_eff_ca = Ω_ca * η_ca * exp(-η_ca^2/2)
analytic_ca_sb = sin.(Ω_eff_ca/2 .* tout_ca_sb).^2

println("\n=== Ca40 vs Ba138 Comparison ===")
println("                  Ca40        Ba138")
println("Wavelength:       729 nm      ", round(laser_ba.λ*1e9, digits=1), " nm")
println("Mass:             40 amu      138 amu")
println("η:                ", round(η_ca, digits=4), "      ", round(η_ba, digits=4))
println("Sideband Ω_eff:   ", round(Ω_eff_ca, digits=4), "      ", round(Ω_eff_ba, digits=4), " rad/μs")
println("Sideband period:  ", round(T_sb_ca, digits=4), " μs  ", round(T_sb_ba, digits=4), " μs")

# Plots

# Plot 1: Carrier Rabi Flopping
p1 = plot(tout_c, ex_ba_carrier, label="Numeric (IonSim)", color =:steelblue, lw=2.5)
plot!(tout_c, analytic_ba, label="Analytic", linestyle=:dash, color =:black, lw = 1.5)
xlabel!(p1, "Time (μs)")
ylabel!(p1, "Population in |D⟩")
title!(p1, "¹³⁸Ba⁺ Carrier Rabi Flopping (λ = 1762 nm)")
ylims!(p1, 0, 1)

# Plot 2: Sideband Driving
p2 = plot(tout_sb, ex_ba_blue, label="Ba138 BSB from n=0", color=:steelblue, lw=2.5)
plot!(p2, tout_sb, ex_ba_red, label="Ba138 RSB from n=0", color=:indianred, lw=2.5)
plot!(p2, tout_sb, analytic_ba_sb, label="Analytic |0⟩↔|1⟩", linestyle=:dot, color=:black, lw=1.5)
xlabel!(p2, "Time (μs)")
ylabel!(p2, "Population in |D⟩")
title!(p2, "¹³⁸Ba⁺ COM Sideband Rabi Flopping (η = $(round(η_ba, digits=4)))")
ylims!(p2, 0, 1)

# Plot 3: Ca40 vs Ba138 Sideband Comparison
p3 = plot(tout_ca_sb ./T_sb_ca, ex_ca_blue,  label="⁴⁰Ca⁺  η=$(round(η_ca,digits=3))", color=:steelblue, lw=4.5)
plot!(p3, tout_sb ./T_sb_ba, ex_ba_blue, label="¹³⁸Ba⁺  η=$(round(η_ba,digits=3))", color=:indianred, lw=2.5)
plot!(p3, 0:0.01:2, t-> sin(π*t)^2, label = "Ideal sin²(πt/2)", color=:black, lw=1.5, linestyle=:dash)
xlabel!(p3, "Time/Sideband Period")
ylabel!(p3, "Population in |D⟩")
title!(p3, "Ca40 vs Ba138 Sideband Rabi Flopping")

plot(p1, p2, p3, layout=(1,3), size=(1600, 450), left_margin = 8mm, bottom_margin = 8mm, top_margin = 5mm)
savefig("IonTrap/Results/Ba138_Ca40_Comparison.png")
println("\nPlots saved to IonTrap/Results/Ba138_Ca40_Comparison.png")
println("Max carrier deviation from analytic: ", round(maximum(abs.(ex_ba_carrier - analytic_ba)), digits=5))