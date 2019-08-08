using JuLIP, ASE, JuLIP.Testing, Test, LinearAlgebra

emt = pyimport("ase.calculators.emt")        # import the EMT model
println("Test Direct use of ASECalculator")
at = bulk(:Cu, cubic=true) * 2    # generate periodic Cu supercell
deleteat!(at, 1)                 # vacancy defect
at = ASEAtoms(at)
try
   calc = ASECalculator(emt.EMT())        # wrap it into a Julia Object
   @show energy(calc, at)                 # compute the energy
   @show maximum(norm.(forces(calc, at)))
   println(@test true)

   atj = Atoms(at)
   println(@test energy(calc, at) ≈ energy(calc, atj))
catch
   println("failed ASECalculator Test")
   println(@test false)
end


println("Test ASE's EMT calculator implementation (ASE.EMTCalculator)")
emt = ASE.Models.EMT()
at = rattle!( set_pbc!( bulk(:Cu, cubic=true) * 2, (true,false,false) ), 0.1 )
println(@test JuLIP.Testing.fdtest(emt, at, verbose=true))
