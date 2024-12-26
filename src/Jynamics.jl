module Jynamics

include("Cylind.jl")
include("File.jl")
include("Tracing.jl")
include("PeriodicSearch.jl")
include("NBodyParticle.jl")
include("NBodyODE.jl")
using .Cylind
using .File
using .Tracing
using .PeriodicSearch
using .NBodyParticle
using .NBodyODE

# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

for name in names(Cylind, all=true) ∪ names(File, all=true) ∪ names(Tracing, all=true) ∪ names(NBodyParticle, all=true) ∪ names(NBodyODE, all=true)  
    # only export the functions that are not built-in
    if !(name in builtin_functions)
        @eval export $name
    end
end
export lookfor_cycle

end
