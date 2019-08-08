
ASECalculator(po::PyObject) = ASECalculator(po, nothing, nothing)

"""
`update!(calc::ASECalculator, at::Atoms)`

updates species, positions, cell, pbc. (only those features that
will affect energy, forces, virial...)
"""
function update!(calc::ASECalculator, at::Atoms)
   # if the `at, aseat` fields haven't been set yet, then
   # start from scratch
   if calc.at == nothing || calc.aseat == nothing
      refresh!(calc, at)
   end
   # if the composition or length has changed, then we also start
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


forces(calc::ASECalculator, at::ASEAtoms) = calc.po.get_forces(at.po)' |> vecs

energy(calc::ASECalculator, at::ASEAtoms) = calc.po.get_potential_energy(at.po)

function virial(calc::ASECalculator, at::ASEAtoms)
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

function forces(calc::ASECalculator, at::Atoms)
   update!(calc, at)
   return forces(calc, calc.aseat)
end

function energy(calc::ASECalculator, at::Atoms)
   update!(calc, at)
   return energy(calc, calc.aseat)
end

function virial(calc::ASECalculator, at::Atoms)
   update!(calc, at)
   return virial(calc, calc.aseat)
end



# -------------- Calling a JuLIP Calculator on ASEAtoms --------------
# This one is tricky!

# energy(V::AbstractCalculator, at::ASEAtoms) = energy(V, Atoms(at))
# forces(V::AbstractCalculator, at::ASEAtoms) = forces(V, Atoms(at))
# virial(V::AbstractCalculator, at::ASEAtoms) = virial(V, Atoms(at))
# stress(V::AbstractCalculator, at::ASEAtoms) = stress(V, Atoms(at))
