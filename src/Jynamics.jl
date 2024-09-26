module Jynamics

include("Cylind.jl")

using .Cylind

# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

for name in names(Cylind, all=true)
    # only export the functions that are not built-in
    if !(name in builtin_functions)
        @eval export $name
    end
end

end
