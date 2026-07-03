# Now will explore Molmer-Sorensen (MS) Entangling Gates
# Using bichromatic drive (laser_blue + laser_red)

using IonSim
using QuantumOptics: timeevolution, expect, dm, normalize, dagger, embed, create, destroy, ptrace, FockBasis, DenseOperator
using Plots
using Plots.PlotMeasures

# Setup
ca1 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
ca2 = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])

# IonSim adds the Hamiltonians automatically 
laser_blue = Laser()
laser_red = Laser()

chain = LinearChain(
    ions = [ca1, ca2],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1,2])
)

chamber = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [laser_blue, laser_red])

# Set both lasers to the same wavelength, geometry, and intensity
for laser in [laser_blue, laser_red]
    wavelength_from_transition!(laser, ca1, ("S", "D"), chamber)
    polarization!(laser, (x̂ - ẑ)/√2)
    wavevector!(laser, (x̂ + ẑ)/√2)
    intensity_from_pitime!(laser, 2e-6, ca1, ("S", "D"), chamber)
end

mode_com, mode_stretch = modes(chamber)
ν_com = frequency(mode_com)
η = abs(lambdicke(mode_com, ca1, laser_blue)) # Lamb-Dicke parameter for stretch mode
Ω = π/2 # Bare carrier Rabi freq, rad/μs (pi-time 2 μs -> Ω = π/τ_π)

# Gate Parameters

# Gate Condition (one phase-space loop, maximally entangling)

Ω_eff = η * Ω # Effective Rabi frequency for the stretch mode
δ_MS_rpu = 2 * Ω_eff # Detuning for the MS gate (red/blue sideband detuning), rpu = rad/μs
δ_MS_Hz = δ_MS_rpu * 1e6/ (2π) # Convert to Hz for detuning
t_gate = 2π/ δ_MS_rpu # Gate time for one phase-space loop

println("=== MS Gate Parameters ===")
println("η (COM)   = ", round(η, digits=6))
println("Ω_eff (rad/μs) = ", round(Ω_eff, digits=6), " kHz")
println("δ_MS (rad/μs) = ", round(δ_MS_rpu, digits=6), " kHz")
println("t_gate    = ", round(t_gate, digits=6), " μs")
println("Verify φ=π/4:    π(ηΩ/δ)² = ",
        round(π*(Ω_eff/δ_MS_rpu)^2, digits=4), "  (want ", round(π/4, digits=4), ")")

# Build Bichromatic Hamiltonian
detuning!(laser_blue, ν_com + δ_MS_Hz) # Blue sideband detuning
detuning!(laser_red, -(ν_com + δ_MS_Hz)) # Red sideband detuning
h_MS = hamiltonian(chamber, timescale = 1e-6, rwa_cutoff = 1e5)

# Ideal gate: Schrodinger Evolution
ψ0 = ca1["S"] ⊗ ca2["S"] ⊗ mode_com[0] ⊗ mode_stretch[0] 

# Run for 2 gate times so we can see the state return and repeat

tspan = 0:0.2:2 * t_gate # Time span for simulation
tout, sol_ideal = timeevolution.schroedinger_dynamic(tspan, ψ0, h_MS)

# Observables: Population and Bell-State Fidelity

# Target Bell State in Full Hilbert Space (Motional State must return to |0,0⟩)

ψ_bell = normalize(ca1["S"] ⊗ ca2["S"] ⊗ mode_com[0] ⊗ mode_stretch[0] + im * ca1["D"] ⊗ ca2["D"] ⊗ mode_com[0] ⊗ mode_stretch[0])

# Bell Fidelity at each time: F(t) = |<ϕ'|ψ(t)>|²
fidelity_ideal = [abs2(dagger(ψ_bell) * ψ_t) for ψ_t in sol_ideal]

# Two-qubit populations 

function two_qubit_pops(sol)
    ρ2 = [ptrace(dm(ψt), [3,4]) for ψt in sol] # Trace out motional modes (subsystems 3 and 4)
    SS = real.([ρ.data[1,1] for ρ in ρ2]) # |S>_1|S>_2
    SD = real.([ρ.data[2,2] for ρ in ρ2]) # |S>_1|D>_2
    DS = real.([ρ.data[3,3] for ρ in ρ2]) # |D>_1|S>_2
    DD = real.([ρ.data[4,4] for ρ in ρ2]) # |D>_1|D>_2
    return SS, SD, DS, DD
end

ex_SS, ex_SD, ex_DS, ex_DD = two_qubit_pops(sol_ideal)

# Print Gate Result at t_gate
idx_gate = argmin(abs.(tout .- t_gate)) # Find index closest to t_gate

println("\n=== MS Gate Result at t = ", round(tout[idx_gate], digits=1), " μs ===")
println("P(SS) = ", round(ex_SS[idx_gate], digits=4), "(expect ≈0.5)")
println("P(DD) = ", round(ex_DD[idx_gate], digits=4), "(expect ≈0.5)")
println("P(SD) = ", round(ex_SD[idx_gate], digits=4), "(expect ≈0)")
println("P(DS) = ", round(ex_DS[idx_gate], digits=4), "(expect ≈0)")

println("Bell Fidelity = ", round(fidelity_ideal[idx_gate], digits=4), "(expect ≈1)")

# Noisy Gate: Lindblad Master Equation with Motional Heating

# Motional Heating adds phonons via L = √Γ â†_COM
# Where Γ = n_dot * 1e-6 converts phonons/s -> phonons/μs
# Build â†_COM and â_COM in the full 4-subsystem Hilbert space
# Subsystem order: [ca1, ca2, mode_com, mode_stretch] → COM is index 3


b_full = ψ0.basis
vm_basis = mode_com[0].basis # VibrationalMode basis
fb = FockBasis(vm_basis.N) # Fock basis for COM mode with same dimensions

adag_com = embed(b_full, 3, DenseOperator(vm_basis, vm_basis, create(fb).data)) # â†_COM in full Hilbert space
a_com = embed(b_full, 3, DenseOperator(vm_basis, vm_basis, destroy(fb).data)) # â_COM in full Hilbert space
ρ0 = dm(ψ0) # Initial density matrix
ρ_bell = dm(ψ_bell) # Target Bell state density matrix

# Compute Noisy Fidelity for Two Heating Rates

n_dot_values = [1e3, 1e4] # phonons/s
fidelities_noisy = []

for n_dot in n_dot_values
    fid = let n_dot = n_dot # Ensure the closure captures the right n_dot_values
        Γ = n_dot * 1e-6
        J = [√Γ * adag_com] # Lindblad jump operators
        Jdag = [√Γ * a_com] # Its adjoint

        # Master_dynamic expects f(t, ρ) -> (H, J_list, Jdag_list)
        fmaster(t, ρ) = (h_MS(t, ρ), J, Jdag)
        _, sol = timeevolution.master_dynamic(tspan, ρ0, fmaster)
        real.(expect(ρ_bell, sol))
    end
    push!(fidelities_noisy, fid) # Store the fidelity for this heating rate
end

# Plots

# Plot 1: Two-Qubit Populations (Ideal)
p_pop = plot(tout, ex_SS, label = "|SS⟩", color = :steelblue, lw = 2.5)
plot!(p_pop, tout, ex_DD, label = "|DD⟩", color = :indianred, lw = 2.5)
plot!(p_pop, tout, ex_SD, label = "|SD⟩", color = :seagreen, lw = 2.5)
plot!(p_pop, tout, ex_DS, label = "|DS⟩", color = :orange, lw = 2.5)
vline!(p_pop, [t_gate], label = "Gate Time", color = :black, lw = 2, ls = :dash)
xlabel!(p_pop, "Time (μs)")
ylabel!(p_pop, "Population")
title!(p_pop, "Molmer-Sorensen Gate: Two-Qubit Populations (Ideal)")
ylims!(p_pop, -0.05, 1.05)

# Plot 2: Bell State Fidelity (Ideal and Noisy)
p_fid = plot(tout, fidelity_ideal, label = "Ideal", color = :black, lw = 2.5)
plot!(p_fid, tout, fidelities_noisy[1], label = "Noisy (ṅ = 1e3 phonons/s)", color = :steelblue, lw = 2.5)
plot!(p_fid, tout, fidelities_noisy[2], label = "Noisy (ṅ = 1e4 phonons/s)", color = :indianred, lw = 2.5)
vline!(p_fid, [t_gate], label = "Gate Time", color = :black, lw = 2, ls = :dash)
hline!(p_fid, [π/4 / (π/2)], color =:gray, lw = 1, ls =:dash, label = "") # Visual Guide
xlabel!(p_fid, "Time (μs)")
ylabel!(p_fid, "Bell State Fidelity F = |⟨ϕ'|ψ(t)⟩|²")
title!(p_fid, "Molmer-Sorensen Gate: Bell State Fidelity")
ylims!(p_fid, -0.05, 1.05)

plot(p_pop, p_fid, layout = (1, 2), size = (1100, 430), left_margin = 8mm, bottom_margin = 8mm, top_margin = 5mm)
savefig("IonTrap/Results/MolmerSorensenGate.png")
print("\n Saved MS Gate Plots")
println("Peak Ideal Bell Fidelity: ", round(maximum(fidelity_ideal), digits=4))
