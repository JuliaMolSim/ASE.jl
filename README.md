
# ASE.jl

Julia Bindings for the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/)

[![Build Status](https://travis-ci.org/libAtoms/ASE.jl.svg?branch=master)](https://travis-ci.org/libAtoms/ASE.jl)

### Summary

Provides minimal Julia wrappers for a limited subset of ASE's functionality, to
be used within [JuLIP.jl](https://github.com/libAtoms/JuLIP.jl). On top of
`JuLIP.jl`, which is a pure Julia library, `ASE.jl` also provides an interface
to ASE, via [PyCall.jl](https://github.com/JuliaPy/PyCall.jl).

### Getting Started

To install, simply calling
```julia
Pkg.add("ASE")
```
from the REPL should suffice. 

Quick test
```Julia
using ASE
at = bulk(:Cu, cubic=true) * 2         # generate periodic Cu supercell
deleteat!(at, 1)                       # vacancy defect
emt = pyimport("ase.calculators.emt")  # import the EMT model
calc = ASECalculator(emt.EMT())        # wrap it into a Julia Object
@show energy(calc, at)                 # compute the energy
set_calculator!(at, calc)
minimise!(at)
@show energy(at)
```

### Useful functions

(See their inline documentation for more details.)

* `ASECalculator`: wrap a python object into a `JuLIP.jl` calculator to
compute energies, forces, etc...
* `read_xyz, write_xyz`: functionality for reading and writing .xyz files

`ASE.Models` implements some convenience functions to create ASE calculators:
* `ASE.Models.EMT()`
* `ASE.Models.EAM(filename)`
* `ASE.Models.QuippySW()`
Please file an issue or PR to add more models to this sub-module.

### Internals

The most important type defined in `ASE.jl` is the `ASECalculator`. It
stores the `PyObject` defining the calculator, an `ASEAtoms` object wrapping
an `ase.Atoms` `PyObject` and a `JuLIP.Atoms` object. Each call to `energy`,
`forces`, etc, first updates the two atoms objects and then calls the
the ase calculator.

In `ASE.jl` a bulk cell is generated using `bulk("Cu")` while
in `JuLIP.jl` it is generated using `bulk(:Cu)`. Conversion between
`ASE.ASEAtoms` and `JuLIP.Atoms` is via
```
at1 = bulk(:Si)          # JuLIP.Atoms
at2 = ASEAtoms(at1)      # ASE.ASEAtoms
at3 = Atoms(at2)         # JuLIP.Atoms
@assert at1 == at3
```
A "user" should normally not need to use `ASEAtoms` directly.
