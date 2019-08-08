
# ======================== CALCULATORS ===============================


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
