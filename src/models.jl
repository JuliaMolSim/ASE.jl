
module Models

using PyCall

"""
Creates an `ASECalculator` that uses `ase.calculators.emt` to compute
energy and forces. This is very slow and is only included for
demonstration purposes.
"""
function EMTCalculator()
   @pyimport ase.calculators.emt as emt
   return ASECalculator(emt.EMT())
end



end
