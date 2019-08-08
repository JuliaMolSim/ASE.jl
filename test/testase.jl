

using Test, ASE, NeighbourLists, JuLIP, JuLIP.Testing
using LinearAlgebra: norm

h1("Misc ASE Tests")

# ---------
h3("Conversion between JuLIP.Atoms and ASEAtoms")
at1 = bulk(:Si)          # JuLIP.Atoms
at2 = ASEAtoms(at1)      # ASE.ASEAtoms
@show cell(at2)
at3 = Atoms(at2)         # JuLIP.Atoms
println(@test at1 == at3)

# ---------
h3("Conversion between JuLIP.Atoms and ASEAtoms, Multi-species")
at1 = bulk(:Al) * 3         # JuLIP.Atoms
at1.Z[1:(length(at1)÷2)] .= atomic_number(:Ti)
at1.M[:] .= 1.0 .+ rand(length(at1))
at2 = ASEAtoms(at1)      # ASE.ASEAtoms
@show cell(at2)
at3 = Atoms(at2)         # JuLIP.Atoms
println(@test at1 == at3)


# ======================================================================

h3("Check that the cubic unit cell for Al is reproduced correctly")
a0 = 2.025
p = [ [0.0;0.0;0.0] [0.0;a0;a0] [a0;0.0;a0] [a0;a0;0.0] ]
at = bulk("Al", cubic=true)
println(@test (positions(at) |> mat) == p)

# ======================================================================

h2("Check neighbourlist without periodicity")
# TODO: implement a test with periodicity as well???
at = ASEAtoms(bulk(:Al, cubic=true) * 3)
set_pbc!(at, (false,false,false))
rcut = 1.7 * a0
nlist = neighbourlist(at, rcut)

@info "   ... assemble neighbour list ..."
# create a neighbourlist via a naive double-loop
simple = zeros(length(at), length(at))
X = positions(at)
for n = 2:length(at), m = 1:n-1
   if norm(X[m]-X[n]) <= rcut
      simple[n,m] = simple[m,n] = 1
   end
end

@info "   ... check the bond-iterator ... "
pass_bonds_test = true
for (i,j,r,R) in pairs(nlist)
   if !( (simple[i,j] == 1) && (abs(norm(X[i]-X[j]) - r) < 1e-12) &&
         (norm(X[j]-X[i] - R) < 1e-12)  )
      global pass_bonds_test = false
      break
   end
   # switch the flag
   simple[i,j] = -1
end

h3("check that no false bonds have been found")
println(@test pass_bonds_test)
h3("check that all bonds have been found")
println(@test maximum(simple) == 0)

# revert to original
simple *= -1
h3("   ... check the site iterator ... ")
pass_site_test = true
for (i,j,r,R) in sites(nlist)
   for n = 1:length(j)
      if simple[i,j[n]] != 1
         global pass_site_test = false
         break
      end
      simple[i,j[n]] = -1
   end
   if maximum(simple[i,:]) != 0
      global pass_site_test = false
      break
   end
end
println(@test pass_site_test)
println(@test maximum(simple) == 0)

# ======================================================================

h3("Checking momenta and velocites")
N = length(at)
p = rand(3,N) |> vecs
set_momenta!(at, p)
println(@test isapprox(momenta(at) |> mat, p |> mat))
v = p ./ masses(at)
println(@test isapprox(velocities(at) |> mat, v |> mat))
set_velocities!(at, v/2.0)
println(@test isapprox(velocities(at) |> mat, v/2.0 |> mat))
println(@test isapprox(momenta(at) |> mat, p/2.0 |> mat))

println("Checking chemical_symbols")
println(@test chemical_symbols(at) == ["Al" for i in 1:length(at)])

# ======================================================================

h3("testing read and write")
fname = "test.xyz"
at = ASEAtoms(bulk(:Cu) * 2)
write_xyz(fname, at)
at2 = read_xyz(fname)
println(@test positions(at) == positions(at2))
rm(fname)
