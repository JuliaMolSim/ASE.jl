
Base.eltype(::ASEAtoms) = Float64

pyobject(at::ASEAtoms) = at.po


"""
`update!(calc::ASECalculator, at::Atoms)`

updates species, positions, cell, pbc. (only those features that
will affect energy, forces, virial...)
"""
function update!(calc::ASECalculator, at::Atoms)
   # if the composition or length has changed, then we just start
   # from scratch
   if calc.at.Z != at.Z
      refresh!(calc, at)
   end
   if at.X != calc.at.X
      set_positions!(calc.at, at.X)
      set_positions!(calc.aseat, at.X)
   end
   if at.cell != calc.at.cell
      set_cell!(calc.at, at.cell)
      set_cell!(calc.aseat, at.cell)
   end
   if at.pbc != calc.at.pbc
      set_pbc!(calc.at, at.pbc)
      set_cell!(calc.aseat, at.pbc)
   end
end

function refresh!(calc::ASECalculator, at::Atoms)
   calc.at = deepcopy(at)
   calc.aseat = ASEAtoms(at)
end


ASEAtoms(syms::AbstractVector{<:String}, X::AbstractVector{<: JVec}) =
      ASEAtoms(ase_atoms.Atoms(collect(syms), collect(mat(X)')))

positions(at::ASEAtoms) = at.po.get_positions()' |> vecs |> collect
set_positions!(at::ASEAtoms, X::AbstractMatrix) =
      (at.po.set_positions(collect(X')); at)
set_positions!(at::ASEAtoms, X::AbstractVector{JVecF}) =
      set_positions!(at, mat(X))

length(at::ASEAtoms) = length(positions(at))

set_pbc!(at::ASEAtoms, val::Bool) = set_pbc!(at, (val,val,val))
set_pbc!(at::ASEAtoms, val::NTuple{3,Bool}) = (at.po.pbc = val; at)
set_pbc!(at::ASEAtoms, val::AbstractVector{Bool}) = set_pbc!(at, tuple(val...))
pbc(a::ASEAtoms) = a.po.pbc


set_cell!(at::ASEAtoms, p::Matrix) = (at.po.set_cell(p); at)
cell(at::ASEAtoms) = _convert_cell(at.po.get_cell())
# the _convert_cell is needed since on later ASE versions the cell is
# no longer just a numpy array but its own cell class.
_convert_cell(C::Array) = C
_convert_cell(C::PyObject) = C.__array__()


# special arrays: momenta, velocities, masses, chemical_symbols
"get the momenta array"
momenta(at::ASEAtoms) = at.po.get_momenta()' |> vecs
"set the momenta array"
set_momenta!(at::ASEAtoms, p::AbstractVector{<:JVec}) = at.po.set_momenta(p |> mat |> PyReverseDims)
"get the velocities array (convert from momenta)"
velocities(at::ASEAtoms) = at.po.get_velocities()' |> vecs
"convert to momenta, then set the momenta array"
set_velocities!(at::ASEAtoms, v::AbstractVector{<:JVec}) = at.po.set_velocities(v |> mat |> PyReverseDims)
"get Vector of atom masses"
masses(at::ASEAtoms) = at.po.get_masses()
"set atom mass array as Vector{Float64}"
set_masses!(at::ASEAtoms, m::Vector{Float64}) = at.po.set_masses(m)
"return vector of chemical symbols as strings"
chemical_symbols(at::ASEAtoms) = pyobject(at).get_chemical_symbols()
"set the chemical symbols"
set_chemical_symbols!(at::ASEAtoms, s::Vector{T}) where {T <: AbstractString} =
   pyobject(at).set_chemical_symbols(s)
"return vector of atomic numbers"
atomic_numbers(at::ASEAtoms) = pyobject(at).get_atomic_numbers()



# ===================== Translations / Conversion ======================

# TODO: write tests for consistency of these conversions
#       for the new test_ase code  \\  X, P, M, Z, cell, pbc
#       also the |> collect should not be needed, probably too strict
#       typing in JuLIP
Atoms(at_ase::ASE.ASEAtoms) =
   Atoms( positions(at_ase) |> collect,
          momenta(at_ase) |> collect,
          masses(at_ase) |> collect,
          atomic_numbers(at_ase) |> collect,
          JMat{Float64}(cell(at_ase)),
          pbc(at_ase),
          nothing )

function ASEAtoms(at::Atoms)
   syms = string.(chemical_symbols(at))
   at_ase = ASEAtoms(syms, positions(at))
   set_momenta!(at_ase, momenta(at))
   set_masses!(at_ase, masses(at))
   set_cell!(at_ase, Matrix(cell(at)))
   set_pbc!(at_ase, tuple(pbc(at)...))
   set_calculator!(at_ase, calculator(at))
end


# ===================== BUILD ATOMS =================================

bulk(s::AbstractString, args...; pbc=true, kwargs...) =
   set_pbc!(ASEAtoms(ase_build.bulk(s, args...; kwargs...)), pbc)
