
# ------------------------------------------------------
# Expose ASE FIO functionality

"""
* `write_xyz(filename, at, mode=:write)` : write atoms object to `filename`
* `write_xyz(filehandle, at)` : write atoms object as xyz file
* `write_xyz(filename, ats::Vector{ASEAtoms}, mode=:write)` : write a time series to a file
* `write_xyz(filename, at, x::Vector{Dofs}, mode=:write)` : write a time series to a file

to append to an existing file, use `:append` or `"a"` instead of `:write`.
"""
write_xyz(filename::AbstractString, at::ASEAtoms, mode=:write) =
   mode == :write ? ase_io.write(filename, at.po) : write_xyz(filename, [at], mode)

write_xyz(filehandle::PyObject, at::ASEAtoms) = ase_io.write(filehandle, at.po, format="xyz")

# open and close files from Python (to get a python filehandle)
pyopenf(filename::AbstractString, mode::AbstractString) = py"open($(filename), $(mode))"
pyclosef(filehandle) = filehandle.close()

function write_xyz(filename::AbstractString, at::ASEAtoms, xs::AbstractVector{Dofs}, mode=:write)
   x0 = dofs(at) # save the dofs
   filehandle = pyopenf(filename, string(mode)[1:1])
   for x in xs
     write_xyz(filehandle, set_dofs!(at, x))
   end
   pyclosef(filehandle)
   set_dofs!(at, x0)   # revert to original configuration
end

function write_xyz(filename::AbstractString, ats::AbstractVector{ASEAtoms}, mode=:write)
   filehandle = pyopenf(filename, string(mode)[1:1])
   for at in ats
      write_xyz(filehandle, at)
   end
   pyclosef(filehandle)
end


read_xyz(filename::AbstractString) = ASEAtoms(ase_io.read(filename))
