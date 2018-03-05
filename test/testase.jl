

using ASE, NeighbourLists


# ======================================================================

println("Check that the cubic unit cell for Al is reproduced correctly")
a0 = 2.025
p = [ [0.0;0.0;0.0] [0.0;a0;a0] [a0;0.0;a0] [a0;a0;0.0] ]
at = bulk("Al", cubic=true)
@test (positions(at) |> mat) == p

# ======================================================================

println("Check neighbourlist without periodicity")
# TODO: implement a test with periodicity as well???
println("   ... assemble neighbour list ...")
at = bulk("Al", cubic=true) * 3
set_pbc!(at, (false,false,false))
rcut = 1.7 * a0
nlist = neighbourlist(at, rcut)
# create a neighbourlist via a naive double-loop
simple = zeros(length(at), length(at))
X = positions(at)
for n = 2:length(at), m = 1:n-1
   if norm(X[m]-X[n]) <= rcut
      simple[n,m] = simple[m,n] = 1
   end
end
println("   ... check the bond-iterator ... ")
pass_bonds_test = true
for (i,j,r,R) in pairs(nlist)
   if !( (simple[i,j] == 1) && (abs(norm(X[i]-X[j]) - r) < 1e-12) &&
         (norm(X[j]-X[i] - R) < 1e-12)  )
      pass_bonds_test = false
      break
   end
   # switch the flag
   simple[i,j] = -1
end
@test pass_bonds_test
# check that all pairs have been found
@test maximum(simple) == 0
# revert to original
simple *= -1
println("   ... check the site iterator ... ")
pass_site_test = true
for (i,j,r,R) in sites(nlist)
   for n = 1:length(j)
      if simple[i,j[n]] != 1
         pass_site_test = false
         break
      end
      simple[i,j[n]] = -1
   end
   if maximum(simple[i,:]) != 0
      pass_site_test = false
      break
   end
end
@test pass_site_test
@test maximum(simple) == 0

# ======================================================================

println("Checking `***_data`, `***_array`, `***_info`, `***_transient`")

at = bulk("Cu") * 3
N = length(at)
# set an array and test that it is read back correctly
z = rand(N)
set_array!(at, "z", z)
@test get_array(at, "z") == z
# set some info and test reading
i = "some info"
set_info!(at, "i", i)
@test get_info(at, "i") == i
# ***_transient should be tested automatically via the calculators.
# test the has_***
@test has_array(at, "z")
@test has_info(at, "i")
@test !has_array(at, "i")
@test !has_info(at, "z")


println("Checking momenta and velocites")
p = rand(3,N) |> vecs
set_momenta!(at, p)
@test isapprox(momenta(at) |> mat, p |> mat)
v = p ./ masses(at)
@test isapprox(velocities(at) |> mat, v |> mat)
set_velocities!(at, v/2.0)
@test isapprox(velocities(at) |> mat, v/2.0 |> mat)
@test isapprox(momenta(at) |> mat, p/2.0 |> mat)

println("Checking chemical_symbols")
@test chemical_symbols(at) == ["Cu" for i in 1:length(at)]

# ======================================================================

println("testing static_neighbourlist")
nlist = static_neighbourlist(at, rnn("Cu") * 1.1)
rattle!(at, 0.05)
nlist1 = static_neighbourlist(at, rnn("Cu") * 1.1)
@test nlist === nlist1

# ======================================================================

println("testing read and write")
fname = "test.xyz"
at = bulk("Cu") * 2
write(fname, at)
at2 = read(fname)
@assert positions(at) == positions(at2)
rm(fname)



# ---------
println("Conversion between JuLIP.Atoms and ASEAtoms")
at1 = JuLIP.bulk(:Si)    # JuLIP.Atoms
at2 = ASEAtoms(at1)      # ASE.ASEAtoms
at3 = Atoms(at2)         # JuLIP.Atoms
@assert at1 == at3



# IMPORTED FROM TESTPOTENTIALS:

# # [2] ASE's EMT calculator
# emt = JuLIP.ASE.EMTCalculator()
# at = rattle!( set_pbc!( bulk("Cu", cubic=true) * 2, (true,false,false) ), 0.1 )
# set_calculator!(at, emt)
# if notCI
#    push!(calculators, (emt, at))
# end

# println("--------------------------------------------------")
# println(" EMT Consistency test: ")
# println("--------------------------------------------------")
# println(" E_ase - E_jl = ", energy(at) - energy(at2))
# println(" |Frc_ase - Frc_jl = ", maxnorm(forces(at) - forces(at2)))
# println("--------------------------------------------------")
# @test abs(energy(at) - energy(at2)) < 1e-10


println("attempt to work with an ASECalculator")
at = bulk("Cu", cubic=true) * 2    # generate periodic Cu supercell
deleteat!(at, 1)                       # vacancy defect
try
   @pyimport ase.calculators.emt as emt   # import the EMT model
   calc = ASECalculator(emt.EMT())        # wrap it into a Julia Object
   @show energy(calc, at)                 # compute the energy
   @show maximum(norm.(forces(calc, at)))
   @test true
catch
   println("failed ASECalculator Test")
   @test false
end
