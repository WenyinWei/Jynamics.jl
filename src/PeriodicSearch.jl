module PeriodicSearch

using ..Cylind: CylindricalVectorField, VR_VZ_VPhi_interp, RVpoloVPhi_pRpZ_interp

# import DifferentialEquations as DE
# using DifferentialEquations: ODEProblem, solve
# using OrdinaryDiffEq
import DifferentialEquations as DE
import LinearAlgebra as LA


function lookfor_cycle(v::CylindricalVectorField, init_rzphi, tor_turn=1, h=0.98, delta_threshold=1e-6; cb=nothing)
    if isnothing(cb) 
        R, Z = v.R, v.Z
        condition(x,t,integrator) = x[1]>R[end-3] || x[1] < R[3] || x[2] >Z[end-3] || x[2] < Z[3]
        affect!(integrator) = DE.terminate!(integrator)
        cb = DE.DiscreteCallback(condition,affect!)
    end

    phi = init_rzphi[3]
    tspan = (phi, phi + tor_turn*2pi)
    x0rz = init_rzphi[1:2]
    BR_interp, BZ_interp, BPhi_interp = VR_VZ_VPhi_interp(v)
    
    function dXpol_dphi!(dX,X,p,phi)
        phimod = mod( phi, 2pi/v.nSym )
        dX[1] = X[1] * BR_interp(X[1], X[2], phimod) / BPhi_interp(X[1], X[2], phimod)
        dX[2] = X[1] * BZ_interp(X[1], X[2], phimod) / BPhi_interp(X[1], X[2], phimod)
    end

    FLT_A = RVpoloVPhi_pRpZ_interp(v)

    while true
        _prob_Xpol = DE.ODEProblem(dXpol_dphi!, x0rz, tspan)
        _sol_Xpol = DE.solve(_prob_Xpol, DE.RK4(), 
            abstol=1e-9, reltol=1e-6, maxiters=1e9, dt=pi/128, callback=cb) 
        
        function dDXpol_dphi!(dDXpol,DXpol,p,phi)
            phimod = mod( phi, 2pi/v.nSym )
            r,z = _sol_Xpol(phi)
            dDXpol[:,:] = FLT_A(r,z,phimod) * DXpol
        end
        _prob_DXpol = DE.ODEProblem(dDXpol_dphi!, [1 0; 0 1], tspan, )
        
        _sol_DXpol = DE.solve(_prob_DXpol, DE.RK4(), 
            abstol=1e-9, reltol=1e-6, maxiters=1e9, dt=pi/128,)
        
        delta =  - h * LA.inv(_sol_DXpol(phi + tor_turn*2pi) - LA.I(2) ) * (_sol_Xpol(phi + tor_turn*2pi) - x0rz) 
        x0rz += delta
        # println(x0rz)
        if LA.norm(delta) < delta_threshold
            break
        end
    end
    return push!(x0rz, phi)
end
export lookfor_cycle

# for name in names(@__MODULE__; all=true)
#     if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
#         @eval export $name
#     end
# end


end