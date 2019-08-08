h1("JuLIP vs ASE Calculator Tests")
using Test, JuLIP.Testing, ASE, PyCall, JuLIP, LinearAlgebra

@info("These tests to compare JuLIP vs ASE implementations of some potentials")

h3("Compare JuLIP vs ASE: EMT")
# JuLIP's EMT implementation
at = set_pbc!( bulk(:Cu, cubic=true) * (2,2,2), (true,false,false) )
rattle!(at, 0.1)
emt = EMT(at)

@info("Test JuLIP vs ASE EMT implementation")
pyemt = ASE.Models.EMT()
print("   energy: ")
println(@test abs(energy(emt, at) - energy(pyemt, at)) < 1e-10)
print("   forces: ")
println(@test norm(forces(pyemt, at) - forces(emt, at), Inf) < 1e-10)
# ------------------------------------------------------------------------


h3("Compare JuLIP vs ASE: EMT - Multi-species")
# JuLIP's EMT implementation
at = set_pbc!( bulk(:Cu, cubic=true) * (2,2,2), (true,false,false) )
at.Z[5:10] .= atomic_number(:Al)
@show unique(chemical_symbols(at))
rattle!(at, 0.1)
emt = EMT(at)

@info("Test JuLIP vs ASE EMT implementation")
pyemt = ASE.Models.EMT()
print("   energy: ")
println(@test abs(energy(emt, at) - energy(pyemt, at)) < 1e-10)
print("   forces: ")
println(@test norm(forces(pyemt, at) - forces(emt, at), Inf) < 1e-10)


# ------------------------------------------------------------------------


h3("Compare JuLIP vs ASE: EAM")


# pot_file = joinpath(dirname(pathof(JuLIP)), "..", "data", "w_eam4.fs")
pot_file = @__DIR__() * "/w_eam4.fs"
url = "https://www.ctcms.nist.gov/potentials/Download/2013--Marinica-M-C-Ventelon-L-Gilbert-M-R-et-al--W-4/1/w_eam4.fs"
if !isfile(pot_file)
   @info("Download an Example EAM potential")
   run(`curl -o $pot_file $url`)
end

@info("Generate the ASE potential")
eam4_ase = ASE.Models.EAM(pot_file)

@info("Generate low-, med-, high-accuracy JuLIP potential")
eam4_jl1 = EAM(pot_file)
eam4_jl2 = EAM(pot_file; s = 1e-4)
eam4_jl3 = EAM(pot_file; s = 1e-6)

at1 = rattle!(bulk(:W, cubic=true) * 3, 0.1)
at2 = deleteat!(bulk(:W, cubic=true) * 3, 1)
at1_ase = ASEAtoms(at1)
at2_ase = ASEAtoms(at2)

for (i, (at, at_ase)) in enumerate(zip([at1, at2], [at1_ase, at2_ase]))
   @info("Test $i")
   err_low = (energy(eam4_ase, at_ase) - energy(eam4_jl1, at)) / length(at)
   err_med = (energy(eam4_ase, at_ase) - energy(eam4_jl2, at)) / length(at)
   err_hi = (energy(eam4_ase, at_ase) - energy(eam4_jl3, at)) / length(at)
   println("      Low Accuracy energy error:", err_low, "; ",
         (@test abs(err_low) < 0.03))
   println("   Medium Accuracy energy error:", err_med, "; ",
         (@test abs(err_med) < 0.006))
   println("     High Accuracy energy error:", err_hi, "; ",
         (@test abs(err_hi) < 0.0005))
end

# ------------------------------------------------------------------------

h3("Compare JuLIP vs Quippy: SW")
if ASE.Models.has_quippy()
   swq = ASE.Models.QuippySW()
   sw = StillingerWeber()
   at0 = bulk(:Si) * (4, 4, 4)
   for k = 1:30
      at = rattle!(deepcopy(at0), 0.1)
      at_ase = ASEAtoms(at)
      err = (energy(swq, at_ase) - energy(sw, at))/length(at)
      print_tf(@test abs(err) < 1e-6)
   end
   println()
else
   @info("quippy was not found")
end
