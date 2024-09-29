module Jynamics

include("Cylind.jl")
include("File.jl")
include("Tracing.jl")
using .Cylind
using .File
using .Tracing

# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

for name in names(Cylind, all=true) ∪ names(File, all=true) ∪ names(Tracing, all=true)
    # only export the functions that are not built-in
    if !(name in builtin_functions)
        @eval export $name
    end
end

end
