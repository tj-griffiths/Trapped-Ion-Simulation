# Simulating Trapped-Ion Quantum Systems with Laser-Ion Interactions, Molmer-Sorensen Gates, and Motional Mode Coupling

using Pkg
# Pkg.add("IonSim", "QuantumOptics", "Plots"]
using IonSim
using QuantumOptics: timeevolution, expect
using Plots

# Step 1: Single Ion, Carrier Rabi Flopping

# Define Ion:
ca = Ca40([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])
# Energy levels and magnetic sublevel 'm'

# Define Laser
laser = Laser()
# No wavelength set, will set to transition S -> D

# Define the Ion Trap Configuration

# Computes normal modes for a Coulomb crystal with center-of-mass trap freq: x = 3Mhz, y = 3Mhz, and z = 1Mhz
chain = LinearChain(
    ions = [ca],
    comfrequencies = (x=3e6, y=3e6, z=1e6),
    selectedmodes = (;z=[1]) # Keep first selected mode along z-field
)

# Assembling Chamber

# Creates  a chamber with the ties between the ion trap and the laser, using a magnetic field of 4e-4 Tesla along the z-axis
chamber = Chamber(iontrap = chain, B = 4e-4, Bhat = ẑ, δB = 0, lasers = [laser])

# Laser Tuning
# Tunes the laser wavelength to that of the transition between the S and D states of the Ca40 ion
wavelength_from_transition!(laser, ca, ("S", "D"), chamber)

# Need specific  laser polarization and wavelength to account for selection rules and B-field direction
polarization!(laser, (x̂ - ẑ)/√2) # Polarization along x-z diagonal
wavevector!(laser, (x̂ + ẑ)/√2) # Wavevector along x-z diagonal

# Set laser intensity for population swap (pi pulse) - takes 2ms as Rabi Frequency
intensity_from_pitime!(laser, 2e-6, ca, ("S", "D"), chamber)

# Build Hamiltonian
h = hamiltonian(chamber, timescale = 1e-6) 

# Initial state
mode = modes(chamber)[1] # Get the first mode
ψ0 = ca["S"] ⊗ mode[0]  # Initial state: ion in S state and mode in ground state (mode[0] is Fock State |0⟩)

# Time-Evolve
tspan = 0:0.05:6 # 0 to 6 ms in steps of 0.05 ms (3 full pi-pulses)
tout, sol = timeevolution.schroedinger_dynamic(tspan, ψ0, h)

ex = real.(expect(ionprojector(chamber, "D"), sol)) # Expectation value of the D state population (|D⟩⟨D|)

# Compare with analytical solution for Rabi flopping
# With pi-time τ_π = 2 ms, the Rabi frequency Ω = π / τ_π = π / 2e-3 s
τ_π = 2
Ω = π / τ_π
analytic = sin.(Ω .* tout ./2).^2 # Analytical solution for Rabi flopping

# Plot

plot(tout, ex, label = "Numerical (IonSim)", lw = 2.5)
plot!(tout, analytic, label = "Analytical sin²(Ωt/2)", lw = 2.5, ls = :dash)
xlabel!("Time (μs)")
ylabel!("Population in |D⟩")
title!("Single ⁴⁰Ca⁺ Carrier Rabi Flopping, τ_π = 2μs")
ylims!(0, 1)

savefig("IonTrap/Results/Carrier_Rabi_Flopping.png")
println("Saved plot to Carrier_Rabi_Flopping.png")
println("Max Deviation from Analysic Curve: ", maximum(abs.(ex - analytic))) 
