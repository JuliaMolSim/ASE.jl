using JuLIP
using JuLIP.Potentials: StillingerWeber
using JuLIP.Solve: minimise!
using JuLIP.Constraints: VariableCell
using ASE: ASECalculator, ASEAtoms
using PyCall

@pyimport ase.units as units
@pyimport quippy.potential as quippy_potential

# Define Stillinger-Weber potential from quippy
sw_pot = quippy_potential.Potential("IP SW")
# Wrap it into a calculator
sw_calc_Q = ASECalculator(sw_pot)

# Define Stillinger-Weber potential from JuLIP
sw_calc_J = StillingerWeber()

at_J = bulk(:Si) * (4, 4, 4)
for k=1:100
    # Define perturbation of the positions
    at = rattle!(at_J, 0.1)

    # Convert into ASE format
    at_ASE = ASEAtoms(at)

    # Compare the two implementations
    @assert energy(sw_calc_J, at_ASE) - energy(sw_calc_Q, at_ASE) < 1e-6
end
