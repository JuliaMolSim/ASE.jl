
module Models

using PyCall
using ASE: ASECalculator

"""
`EMT()`

Creates an `ASECalculator` that uses `ase.calculators.emt` to compute
energy and forces. This is very slow and is only included for
demonstration purposes.
"""
function EMT()
   emt = pyimport("ase.calculators.emt")
   return ASECalculator(emt.EMT())
end


"""
`EAM(potfile)`

Creates an `ase.calculators.eam` calculator to compute
energy and forces. This is very slow and is only included for
demonstration purposes.
"""
function EAM(potfile)
   eam = pyimport("ase.calculators.eam")
   return ASECalculator( eam.EAM(potential=potfile) )
end


function QuippyPot(id)
   quippy_potential = pyimport("quippy.potential")
   return ASECalculator(quippy_potential.Potential(id))
end

QuippySW() = QuippyPot("IP SW")

function has_quippy()
   try
      pyimport("quippy.potential")
      return true
   catch
      return false
   end
end


end
