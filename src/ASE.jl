
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
       ASECalculator,
       read_xyz, write_xyz


build = ase_build = pyimport("ase.build")
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
end


mutable struct ASECalculator <: AbstractCalculator
   pycalc::PyObject
   aseat::ASEAtoms
   at::Atoms{Float64}
end

# main functionality for ASEAtoms
include("atoms.jl")

# main functionality for ASECalculator
include("calc.jl")

# xyz file io
include("fio.jl")

# wrapper for the matscipy neighbourlist
include("nlist.jl")

# some simple wrappers for some ase models
include("models.jl")


end
