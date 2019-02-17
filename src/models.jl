
module Models

using PyCall

# TODO: figure out how to load this dynamically 
@pyimport ase.calculators.emt as emt

"""
Creates an `ASECalculator` that uses `ase.calculators.emt` to compute
energy and forces. This is very slow and is only included for
demonstration purposes.
"""
function EMTCalculator()
   return ASECalculator(emt.EMT())
end



end
