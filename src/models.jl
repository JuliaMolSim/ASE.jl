
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


function QuippyPot(id; kwargs...)
   QuippyPotential = pyimport("quippy.potential").Potential
   return ASECalculator(QuippyPotential(id; kwargs...))
end

function QuippySW()
   QuippyPotential = pyimport("quippy.potential").Potential
   qpot = QuippyPotential("IP SW", param_str = """<SW_params n_types="3" label="PRB_31_plus_H_Ge">
   <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements.  Ge and Si-Ge from Ethier and Lewis '92 </comment>
   <per_type_data type="1" atomic_num="1" />
   <per_type_data type="2" atomic_num="14" />
   <per_type_data type="3" atomic_num="32" />
   <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
         p="0" q="0" a="1.0" sigma="1.0" eps="0.0" />
   <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
         p="4" q="0" a="1.25" sigma="2.537884" eps="2.1672" />
   <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

   <!-- pair data for Ge and Si-Ge -->
   <per_pair_data atnum_i="32" atnum_j="32" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.181" eps="1.92590365783410138248" />
   <per_pair_data atnum_i="14" atnum_j="32" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.138" eps="2.04313391101890694121" />

   <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
   <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

   <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

   <!-- triplet data for Ge and Si-Ge -->
   <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="32"
      lambda="23.1" gamma="1.20" eps="2.10444772465437788017" />
   <per_triplet_data atnum_c="14" atnum_j="32" atnum_k="32"
      lambda="25.5" gamma="1.20" eps="2.04326828917050691242" />
   <per_triplet_data atnum_c="32" atnum_j="14" atnum_k="14"
      lambda="25.5" gamma="1.20" eps="2.04326828917050691242" />
   <per_triplet_data atnum_c="32" atnum_j="14" atnum_k="32"
      lambda="28.1" gamma="1.20" eps="1.98396169354838709676" />
   <per_triplet_data atnum_c="32" atnum_j="32" atnum_k="32"
      lambda="31.0" gamma="1.20" eps="1.92590365783410138247" />

   </SW_params>""")
   return ASECalculator(qpot)
end

function has_quippy()
   try
      pyimport("quippy")
      return true
   catch
      return false
   end
end


end
