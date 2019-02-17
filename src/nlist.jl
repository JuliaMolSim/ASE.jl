

############################################################
# matscipy neighbourlist functionality
############################################################

include("MatSciPy.jl")

function neighbourlist(at::ASEAtoms, cutoff::Float64)::MatSciPy.NeighbourList
   # TODO: also recompute if rcut is different !!!!!
   # if no previous neighbourlist is available, compute a new one
   # if !has_transient(at, (:nlist, cutoff)) || recompute
   #    # this nlist will be destroyed as soon as positions change
   #    set_transient!(at, (:nlist, cutoff), MatSciPy.NeighbourList(at, cutoff))
   # end
   # return get_transient(at, (:nlist, cutoff))
   return MatSciPy.NeighbourList(at, cutoff)
end

# TODO: no longer needed in the new ASE design
# """
# `static_neighbourlist(at::ASEAtoms, rcut::Float64)`
#
# This function first checks whether a static neighbourlist already exists
# with cutoff `rcut` and if it does then it returns the existing list.
# If it does not, then it computes a new neighbour list with the current
# configuration, stores it for later use and returns it.
# """
# function static_neighbourlist(at::ASEAtoms, rcut::Float64)
#    # if !has_transient(at, (:snlist, rcut))
#    #    set_transient!( at, (:snlist, rcut),
#    #                        MatSciPy.NeighbourList(at, rcut),
#    #                        Inf )    # Inf means this is never deleted!
#    # end
#    # return get_transient(at, (:snlist, rcut))
#    return MatSciPy.NeighbourList(at, rcut)
# end
