
# ASE.jl

Julia Bindings for the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)

[![Build Status](https://travis-ci.org/libAtoms/ASE.jl.svg?branch=master)](https://travis-ci.org/libAtoms/ASE.jl)

### Summary

Provides Julia wrappers for a limited subset of ASE's functionality, to be used within
[JuLIP.jl](https://github.com/libAtoms/JuLIP.jl). On top of `JuLIP.jl`, which is
a pure Julia library, `ASE.jl` also provides an interface to ASE, via [PyCall.jl](https://github.com/JuliaPy/PyCall.jl).

### Getting Started

To install
```julia
Pkg.add("ASE")
```

Quick test
```Julia
using ASE
at = bulk("Cu", cubic=true) * 2        # generate periodic Cu supercell
deleteat!(at, 1)                       # vacancy defect
emt = pyimport("ase.calculators.emt")  # import the EMT model
calc = ASECalculator(emt.EMT())        # wrap it into a Julia Object
@show energy(calc, at)                 # compute the energy
# -------------------------------------------
#  or to use more of the JuLIP framework:
# -------------------------------------------
set_calculator!(at, calc)
set_constraint!(at, FixedCell(at))
minimise!(at)
@show energy(at)
```


Note that in `ASE.jl` a bulk cell is generated using `bulk("Cu")` while
in `JuLIP.jl` it is generated using `bulk(:Cu)`. Conversion between
`ASE.ASEAtoms` and `JuLIP.Atoms` is via
```
at1 = bulk(:Si)          # JuLIP.Atoms
at2 = ASEAtoms(at1)      # ASE.ASEAtoms
at3 = Atoms(at2)         # JuLIP.Atoms
@assert at1 == at3
```

### TODO

* better integration with JuLIP, i.e. work with `JuLIP.Atoms` instead of
`ASEAtoms`, but this will require rewriting the `ASECalculators` a bit
* provide more convenience functions to call ASE functionality
