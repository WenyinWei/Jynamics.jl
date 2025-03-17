module Jynamics

include("DiscreteUtils.jl")
include("Cylind.jl")
include("File.jl")
include("Tracing.jl")
include("PeriodicSearch.jl")
<<<<<<< HEAD
include("FourierUtils.jl")
using .DiscreteUtils
=======
include("NBodyParticle.jl")
include("NBodyODE.jl")
>>>>>>> b9f324b2c3147f7d0968404dce2862b958d39b5b
using .Cylind
using .File
using .Tracing
using .PeriodicSearch
<<<<<<< HEAD
using .FourierUtils
=======
using .NBodyParticle
using .NBodyODE
>>>>>>> b9f324b2c3147f7d0968404dce2862b958d39b5b

# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

<<<<<<< HEAD
for name in names(DiscreteUtils, all=true) ∪ names(Cylind, all=true) ∪ names(File, all=true) ∪ names(Tracing, all=true) ∪ names(FourierUtils, all=true) 
=======
for name in names(Cylind, all=true) ∪ names(File, all=true) ∪ names(Tracing, all=true) ∪ names(NBodyParticle, all=true) ∪ names(NBodyODE, all=true)  
>>>>>>> b9f324b2c3147f7d0968404dce2862b958d39b5b
    # only export the functions that are not built-in
    if !(name in builtin_functions)
        @eval export $name
    end
end
export lookfor_cycle

end
