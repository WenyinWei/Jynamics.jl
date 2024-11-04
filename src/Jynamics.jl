module Jynamics

include("DiscreteUtils.jl")
include("Cylind.jl")
include("File.jl")
include("Tracing.jl")
include("PeriodicSearch.jl")
using .DiscreteUtils
using .Cylind
using .File
using .Tracing
using .PeriodicSearch

# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

for name in names(DiscreteUtils, all=true) ∪ names(Cylind, all=true) ∪ names(File, all=true) ∪ names(Tracing, all=true) 
    # only export the functions that are not built-in
    if !(name in builtin_functions)
        @eval export $name
    end
end
export lookfor_cycle

end
