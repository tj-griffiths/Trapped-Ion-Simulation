# Test if Ba138 was properly added to IonSim

using IonSim
using QuantumOptics: timeevolution, expect
using Plots

ba = Ba138([("S1/2", -1/2, "S"), ("D5/2", -1/2, "D")])

laser = Laser()
chain = LinearChain(ions=[ba], comfrequencies = (x=3e6, y=3e6, z =1e6), selectedmodes=(;z=[1]))
chamber = Chamber(iontrap=chain, B=4e-4, Bhat = ẑ, δB = 0, lasers = [laser])
wavelength_from_transition!(laser, ba, ("S", "D"), chamber)
println("Qubit wavelength: ", round(laser.λ * 1e9, digits = 2), " nm (expect ≈ 1762 nm")

mode = modes(chamber)[1]
η_ba = abs(lambdicke(mode, ba, laser))
print("Lamb-Dicke η(Ba138): " , round(η_ba, digits = 5), " (expect < 0.0685 from Ca40)")

intensity_from_pitime!(laser, 5e-6,  ba, ("S", "D"), chamber)
h = hamiltonian(chamber, timescale = 1e-6)
ψ0 = ba["S"] ⊗ mode[0]
tout, sol = timeevolution.schroedinger_dynamic(0:0.1:10, ψ0, h)
println("\nTime evolution: ", length(sol), " steps completed  ✓")
