
__precompile__(false)

module ASE

using Reexport

@reexport using JuLIP
@reexport using PyCall

# the functions to be implemented
import JuLIP:
      positions, set_positions!,
      cell, set_cell!,             # ✓
      pbc, set_pbc!,               # ✓
      calculator, set_calculator!, # ✓
      neighbourlist,                # ✓
      energy, forces, virial, stress,
      momenta, set_momenta!,
      masses, set_masses!,
      set_transient!,
      atomic_numbers,
      Atoms, chemical_symbols,
      get_data, has_data, set_data!,
      bulk

import Base.length, Base.deleteat!, Base.deepcopy

# from arrayconversions:
using JuLIP: mat, vecs, JVecF, JMatF,
      AbstractAtoms,
      AbstractCalculator, maxdist, SVec,
      Dofs, set_dofs!

using LinearAlgebra: det

# extra ASE functionality:
import Base: repeat         # ✓

export ASEAtoms,      # ✓
      AbstractASECalculator, ASECalculator,
      extend!, get_info, set_info!, get_array, set_array!, has_array, has_info,
      get_transient, has_transient,
      velocities, set_velocities!,
      static_neighbourlist,
      read_xyz, write_xyz

ase_build = pyimport("ase.build")
ase_atoms = pyimport("ase.atoms")
ase_io    = pyimport("ase.io")


"""
`type ASEAtoms <: AbstractAtoms`

Julia wrapper for the ASE `Atoms` class.

### Constructors:
```julia
ASEAtoms(po::PyObject)         # from a given ASE object
ASEAtoms(s::AbstractString)    # e.g., "Si2" for a Si cluster containing 2 atoms
ASEAtoms(at::Atoms)            # from a JuLIP Atoms object
```

"""
mutable struct ASEAtoms <: AbstractAtoms{Float64}
   po::PyObject       # ase.Atoms instance
   calc::Union{AbstractCalculator, Nothing}
end

Base.eltype(::ASEAtoms) = Float64

ASEAtoms(po::PyObject) = ASEAtoms(po, nothing)
pyobject(a::ASEAtoms) = a.po

# this one is needed e.g. for JuLIP conversions
ASEAtoms(s::AbstractString) = ASEAtoms(ase_atoms.Atoms(s))

ASEAtoms(syms::AbstractVector{<:String}, X::AbstractVector{<: JVec}) =
      ASEAtoms(ase_atoms.Atoms(collect(syms), collect(mat(X)')))

set_calculator!(at::ASEAtoms, calc::Union{AbstractCalculator, Nothing}) = (at.calc = calc; at)
calculator(at::ASEAtoms) = at.calc

positions(at::ASEAtoms) = at.po.get_positions()' |> vecs |> collect
set_positions!(at::ASEAtoms, X::Matrix) = (at.po.set_positions(convert(Matrix, X')); at)
set_positions!(at::ASEAtoms, X::AbstractVector{JVecF}) = set_positions!(at, Matrix(mat(X)))
# TODO: revert to "clever" set_positions!
# TODO: how to get rid of the `convert`? There is an awful lot of copying going on here!
# function set_positions!(a::ASEAtoms, p::JVecsF)
#    pold = positions(a)
#    r = maxdist(pold, p)
#    p_py = PyReverseDims(mat(p))
#    a.po.set_positions(p_py)
#    update_transient_data!(a, r)
#    return a
# end

length(at::ASEAtoms) = length(positions(at))
# TODO: could be made efficient by staying in Python
#       e.g. size(at.po.get_positions(), 1)

set_pbc!(at::ASEAtoms, val::Bool) = set_pbc!(at, (val,val,val))
set_pbc!(at::ASEAtoms, val::NTuple{3,Bool}) = (at.po.pbc = val; at)
set_pbc!(at::ASEAtoms, val::AbstractVector{Bool}) = set_pbc!(at, tuple(val...))
pbc(a::ASEAtoms) = a.po.pbc


cell(at::ASEAtoms) = _convert_cell(at.po.get_cell())
_convert_cell(C::Array) = C
_convert_cell(C::PyObject) = C.__array__()
set_cell!(at::ASEAtoms, p::Matrix) = (at.po.set_cell(p); at)
# TODO: if we want to return to transient data
# A = pinv(cell(a)) * p
# update_transient_data!(a, vecnorm(A))


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

JuLIP.Chemistry.rnn(s::AbstractString) = rnn(Symbol(s))

#
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
          calculator(at_ase) )

function ASEAtoms(at::Atoms)
   # simplify by assuming there is only one species
   syms = string.(chemical_symbols(at))
   at_ase = ASEAtoms(syms, positions(at))
   set_momenta!(at_ase, momenta(at))
   set_masses!(at_ase, masses(at))
   set_cell!(at_ase, Matrix(cell(at)))
   set_pbc!(at_ase, tuple(pbc(at)...))
   set_calculator!(at_ase, calculator(at))
end


# ===================== BUILD ATOMS =================================

# @doc ase_build.bulk.__doc__ ->
bulk(s::AbstractString, args...; pbc=true, kwargs...) =
   set_pbc!(ASEAtoms(ase_build.bulk(s, args...; kwargs...)), pbc)

build(sym::Symbol, args...; kwargs...) = ase_build[sym](args...; kwargs...)

repeat(a::ASEAtoms, n::NTuple{3, Int64}) = ASEAtoms(a.po.repeat(n))
repeat(a::ASEAtoms, n::AbstractArray) = repeat(a, tuple(n...))

import Base.*
*(at::ASEAtoms, n) = repeat(at, n)
*(n::NTuple{3, Int64}, at::ASEAtoms) = repeat(at, n)
*(at::ASEAtoms, n::Integer) = repeat(at, (n,n,n))
*(n::Integer, at::ASEAtoms) = repeat(at, (n,n,n))

function deleteat!(at::ASEAtoms, n::Integer)
   at.po.__delitem__(n-1) # delete in the actual array
   # update_transient_data!(at, Inf)
   return at
end

# function deleteat!(at::ASEAtoms, nn)
#    deleteat!(at, collect(nn) - 1)
#    return at
# end

# ======================== CALCULATORS ===============================

######################################################
#    Attaching an ASE-style calculator
#    ASE-style aliases
######################################################

"""
abstract type for all calculators that interface to ASE
"""
abstract type AbstractASECalculator <: AbstractCalculator end

"""
Concrete subtype of ASECalculator for classical potentials
"""
mutable struct ASECalculator <: AbstractASECalculator
   po::PyObject
end

function set_calculator!(at::ASEAtoms, calc::AbstractASECalculator)
   at.po.set_calculator(calc.po)
   at.calc = calc
   return at
end

set_calculator!(at::ASEAtoms, po::PyObject) =
      set_calculator!(at, AbstractASECalculator(po))

forces(calc::AbstractASECalculator, at::ASEAtoms) = calc.po.get_forces(at.po)' |> vecs
energy(calc::AbstractASECalculator, at::ASEAtoms) = calc.po.get_potential_energy(at.po)

function virial(calc::AbstractASECalculator, at::ASEAtoms)
    s = calc.po.get_stress(at.po)
    vol = det(cell(at))
    if size(s) == (6,)
      # unpack stress from compressed Voigt vector form
      s11, s22, s33, s23, s13, s12 = s
      return -JMatF([s11 s12 s13;
                     s12 s22 s23;
                     s13 s23 s33]) * vol
    elseif size(s) == (3,3)
      return -JMatF(s) * vol
    else
      error("got unxpected size(stress) $(size(stress)) from ASE")
    end
end



################### extra ASE functionality  ###################
# TODO: we probably want more of this
#       and a more structured way to translate
#

"""
* `extend!(at::ASEAtoms, atadd::ASEAtoms)::ASEAtoms`
* `extend!(at::ASEAtoms, S::AbstractString, x::JVecF)::ASEAtoms`

add `atadd` atoms to `at` and returns `at`; only `at` is modified.

A short variant is
```julia
#  extend!(at, (s, x))  <<<< deprecated
extend!(at, s, x)
```
where `s` is a string, `x::JVecF` a position
"""
function extend!(at::ASEAtoms, atadd::ASEAtoms)::ASEAtoms
   at.po.extend(atadd.po)
   return at
end

function extend!(at::ASEAtoms, atnew::Tuple{S,JVecF}) where S <: AbstractString
   @warn("`extend!(at, (s,x))` is deprecated; use `extend!(at, s, x)` instead")
   extend!(at, atnew[1], atnew[2])
end

extend!(at::ASEAtoms, S::AbstractString, x::JVecF) = extend!(at, ASEAtoms(S, [x]))

"""
* `write_xyz(filename, at, mode=:write)` : write atoms object to `filename`
* `write_xyz(filehandle, at)` : write atoms object as xyz file
* `write_xyz(filename, ats::Vector{ASEAtoms}, mode=:write)` : write a time series to a file
* `write_xyz(filename, at, x::Vector{Dofs}, mode=:write)` : write a time series to a file

to append to an existing file, use `:append` or `"a"` instead of `:write`.
"""
write_xyz(filename::AbstractString, at::ASEAtoms, mode=:write) =
   mode == :write ? ase_io.write(filename, at.po) : write_xyz(filename, [at], mode)

write_xyz(filehandle::PyObject, at::ASEAtoms) = ase_io.write(filehandle, at.po, format="xyz")

# open and close files from Python (to get a python filehandle)
pyopenf(filename::AbstractString, mode::AbstractString) = py"open($(filename), $(mode))"
pyclosef(filehandle) = filehandle.close()

function write_xyz(filename::AbstractString, at::ASEAtoms, xs::AbstractVector{Dofs}, mode=:write)
   x0 = dofs(at) # save the dofs
   filehandle = pyopenf(filename, string(mode)[1:1])
   for x in xs
     write_xyz(filehandle, set_dofs!(at, x))
   end
   pyclosef(filehandle)
   set_dofs!(at, x0)   # revert to original configuration
end

function write_xyz(filename::AbstractString, ats::AbstractVector{ASEAtoms}, mode=:write)
   filehandle = pyopenf(filename, string(mode)[1:1])
   for at in ats
      write_xyz(filehandle, at)
   end
   pyclosef(filehandle)
end


read_xyz(filename::AbstractString) = ASEAtoms(ase_io.read(filename))





# ====================== TODO REWRITE =================================
#
# * deepcopy
# *

# ASEAtoms(s::AbstractString, X::JVecsF) = ASEAtoms(s, mat(X))
# ASEAtoms(s::AbstractString, X::Matrix) = ASEAtoms(ase_atoms.Atoms(s, X'))
# ASEAtoms(s::AbstractString) = ASEAtoms(ase_atoms.Atoms(s))

# function deleteat!(at::ASEAtoms, n::Integer)
#    at.po.__delitem__(n-1) # delete in the actual array
#    update_transient_data!(at, Inf)
#    return at
# end
#
# function deleteat!(at::ASEAtoms, nn)
#    deleteat!(at, collect(nn) - 1)
#    return at
# end


# "nanotube(n, m, length=1, bond=1.42, symbol=\"C\", verbose=False)"
# function nanotube(args...; kwargs...)
#    at = ASEAtoms(ase_build.nanotube(args...; kwargs...))
#    C = Matrix(cell(at))
#    if det(C) ≈ 0.0
#       C[1,1] = 1.0
#       C[2,2] = 1.0
#       @assert !(det(C) ≈ 0.0)
#    end
#    set_cell!(at, C)
#    return at
# end


include("nlist.jl")

include("models.jl")


# -------------- Calling a JuLIP Calculator on ASEAtoms --------------

# energy(V::AbstractCalculator, at::ASEAtoms) = energy(V, Atoms(at))
# forces(V::AbstractCalculator, at::ASEAtoms) = forces(V, Atoms(at))
# virial(V::AbstractCalculator, at::ASEAtoms) = virial(V, Atoms(at))
# stress(V::AbstractCalculator, at::ASEAtoms) = stress(V, Atoms(at))

# -------------- and ASE Calculator on JuLIP Atoms --------------

energy(V::AbstractASECalculator, at::Atoms) = energy(V, ASEAtoms(at))
forces(V::AbstractASECalculator, at::Atoms) = forces(V, ASEAtoms(at))
virial(V::AbstractASECalculator, at::Atoms) = virial(V, ASEAtoms(at))
# stress(V::AbstractASECalculator, at::Atoms) = stress(V, ASEAtoms(at))


end
