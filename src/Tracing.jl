module Tracing

using ..Cylind: CylindricalVectorField, VR_VZ_VPhi_interp, RVpoloVPhi_pRpZ_interp, RVpoloVPhi_pRpZ_delta_interp, RVpoloVPhi_pRpZ_dirderivative_interp, delta_RVpoloVPhi_interp


using Memoization
@memoize function get_tracing_sols(v::CylindricalVectorField, r0z0, phi0::Number)
    dict = Dict{String, Any}()
    dict["v"] = v  # 存储 CylindricalVectorField 实例的引用
    dict["r0z0"] = r0z0  # 存储 r0z0 的引用
    dict["phi0"] = phi0  # 存储 phi0 的引用
    return dict
end


function dXpol_dphi!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    BR_interp, BZ_interp, BPhi_interp = VR_VZ_VPhi_interp(v)

    function dXpol_dphi!(dX,X,p,phi)
        phimod = mod( phi, 2pi/v.nSym )
        dX[1] = X[1] * BR_interp(X[1], X[2], phimod) / BPhi_interp(X[1], X[2], phimod)
        dX[2] = X[1] * BZ_interp(X[1], X[2], phimod) / BPhi_interp(X[1], X[2], phimod)
    end
    return dXpol_dphi!
end

function dDXpol_dphi!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    sol_Xpol = sols_dict["sol_Xpol"]
    FLT_A =  RVpoloVPhi_pRpZ_interp(sols_dict["v"])

    function dDXpol_dphi!(dDXpol,DXpol,p,phi)
        phimod = mod( phi, 2pi/v.nSym )
        r,z = sol_Xpol(phi)
        dDXpol[:,:] = FLT_A(r,z,phimod) * DXpol
    end
    return dDXpol_dphi!
end

function dDPm_dphi!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    sol_Xpol = sols_dict["sol_Xpol"]
    FLT_A =  RVpoloVPhi_pRpZ_interp(v)

    function dDPm_dphi!(dDPm,DPm,p,phi)
        phimod = mod( phi, 2pi/v.nSym )
        r,z = sol_Xpol(phi)
        dDPm[:,:] = FLT_A(r,z,phimod) * DPm - DPm * FLT_A(r,z,phimod) 
    end
    return dDPm_dphi!
end

function d_delta_Xpol_dphi!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    sol_Xpol = sols_dict["sol_Xpol"]
    FLT_A =  RVpoloVPhi_pRpZ_interp(v)
    RBpoloBPhi_delta_DeltaB = delta_RVpoloVPhi_interp(v, sols_dict["delta_v"])

    function d_delta_Xpol_dphi!(dXpol_delta, Xpol_delta, p, phi)
        phimod = mod( phi, 2pi/v.nSym )
        r,z = sol_Xpol(phi)
        dXpol_delta[:] = FLT_A(r,z,phimod) * Xpol_delta + RBpoloBPhi_delta_DeltaB(r,z,phimod) 
    end
    return d_delta_Xpol_dphi!    
end
function d_delta_DXpol_dphi_init_Xcyc!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    v_pert = sols_dict["delta_v"]
    sol_Xpol = sols_dict["sol_Xpol"]
    sol_DXpol = sols_dict["sol_DXpol"]
    sol_delta_Xcyc = sols_dict["sol_delta_Xcyc"]
    FLT_A =  RVpoloVPhi_pRpZ_interp(v)
    A_delta = RVpoloVPhi_pRpZ_delta_interp(v, v_pert)
    A_dirderivative = RVpoloVPhi_pRpZ_dirderivative_interp(v)

    function d_delta_Xpol_dphi_init_Xcyc!(dXpol_delta, Xpol_delta, p, phi)

        phimod = mod( phi, 2pi/v.nSym )
        r,z = sol_Xpol(phi)
        dr, dz = sol_delta_Xcyc(phi)
        # println( A_delta(r,z,phimod)  )
        # println( A_dirderivative(r,z,phimod,dr,dz) ) 
        dXpol_delta[:,:] = (A_delta(r,z,phimod) + A_dirderivative(r,z,phimod,dr,dz) ) * sol_DXpol(phi) + FLT_A(r,z,phimod) * Xpol_delta
    end
    return d_delta_Xpol_dphi_init_Xcyc!
end

function d_delta_DPm_dphi_init_Xcyc!_generator(sols_dict::Dict{String, Any})
    v = sols_dict["v"]
    v_pert = sols_dict["delta_v"]
    sol_Xpol = sols_dict["sol_Xpol"]
    sol_DPm = sols_dict["sol_DPm"]
    sol_delta_Xcyc = sols_dict["sol_delta_Xcyc"]
    FLT_A =  RVpoloVPhi_pRpZ_interp(v)
    A_delta = RVpoloVPhi_pRpZ_delta_interp(v, v_pert)
    A_dirderivative = RVpoloVPhi_pRpZ_dirderivative_interp(v)

    function d_delta_DPm_dphi_init_Xcyc!(dDPm_delta, DPm_delta, p, phi)
        phimod = mod( phi, 2pi/v.nSym )
        r,z = sol_Xpol(phi)
        dr, dz = sol_delta_Xcyc(phi)
        dDPm_delta[:,:] = (
            (A_delta(r,z,phimod) + A_dirderivative(r,z,phimod,dr,dz) ) * sol_DPm(phi) 
            - sol_DPm(phi) * (A_delta(r,z,phimod) + A_dirderivative(r,z,phimod,dr,dz) ) 
            + FLT_A(r,z,phimod) * DPm_delta 
            - DPm_delta * FLT_A(r,z,phimod)
        )
    end
    return d_delta_DPm_dphi_init_Xcyc!
end


for name in names(@__MODULE__; all=true)
    if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $name
    end
end



end
