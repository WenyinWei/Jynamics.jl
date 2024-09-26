module Cylind

struct CylindricalScalarField{T}
    R::Vector{T}
    Z::Vector{T}
    Phi::Vector{T}
    value::Array{T,3}
    nSym::Int
end

function CylindricalScalarField{T}(R::Vector{T}, Z::Vector{T}, Phi::Vector{T}, value::Array{T,3}; nSym::Int = 1)
    return CylindricalScalarField(R, Z, Phi, value, nSym)
end

function Base.:+(s1::CylindricalScalarField, s2::CylindricalScalarField)
    return CylindricalScalarField(s.R, s.Z, s.Phi, s1.value + s2.value)
end
function Base.:*(a::Number, s::CylindricalScalarField)
    return CylindricalScalarField(s.R, s.Z, s.Phi, a * s.value)
end
function Base.:/(s::CylindricalScalarField, a::Number)
    return CylindricalScalarField(s.R, s.Z, s.Phi, s.value / a)
end

using Memoization
@memoize function pRpZ(s::Array, Rord::Int, Zord::Int, R::Vector, Z::Vector)
    dR = R[2]-R[1]
    dZ = Z[2]-Z[1]

    if Rord==0 && Zord==0
        return s
    elseif Rord > 0
        lastord_field = pRpZ(s, Rord-1, Zord, R, Z)
        thisord_field = similar(s)
        thisord_field[:,:,:] .= NaN
        thisord_field[2:end-1,:,:] = (lastord_field[3:end,:,:]- lastord_field[1:end-2,:,:]) / (2dR)
        return thisord_field
    elseif Zord > 0 
        lastord_field = pRpZ(s, Rord, Zord-1, R, Z)
        thisord_field = similar(s)
        thisord_field[:,:,:] .= NaN
        thisord_field[:,2:end-1,:] = (lastord_field[:,3:end,:]- lastord_field[:,1:end-2,:]) / (2dZ)
        return thisord_field
    end
end
@memoize function pRpZ(s::CylindricalScalarField, Rord::Int, Zord::Int)
    return CylindricalScalarField( s.R, s.Z, s.Phi, pRpZ(s.value, Rord, Zord, R, Z) )
end

using Interpolations
@memoize function pRpZ_interp(s::CylindricalScalarField, Rord::Int, Zord::Int)
    return linear_interpolation( 
            (s.R, s.Z, s.Phi), pRpZ(s, Rord, Zord).value )
end








struct CylindricalVectorField{T}
    R::Vector{T}
    Z::Vector{T}
    Phi::Vector{T}
    VR::Array{T,3}
    VZ::Array{T,3}
    VPhi::Array{T,3}
    nSym::Int
end

function CylindricalScalarField{T}(R::Vector{T}, Z::Vector{T}, Phi::Vector{T}, value::Array{T,3}; nSym::Int=1)
    return CylindricalScalarField{T}(R, Z, Phi, value, nSym)
end

function Base.:+(v1::CylindricalVectorField, v2::CylindricalVectorField)
    return CylindricalVectorField(v1.R, v1.Z, v1.Phi, v1.VR + v2.VR, v1.VZ + v2.VZ, v1.VPhi + v2.VPhi)
end
function Base.:*(a::Number, v::CylindricalVectorField)
    return CylindricalVectorField(v.R, v.Z, v.Phi, a * v.VR, a * v.VZ, a * v.VPhi)
end
function Base.:/(v::CylindricalVectorField, a::Number)
    return CylindricalVectorField(v.R, v.Z, v.Phi, v.VR / a, v.VZ / a, v.VPhi / a)
end

function get_VR(v::CylindricalVectorField)
    return CylindricalScalarField(v.R, v.Z, v.Phi, v.VR)
end
function get_VZ(v::CylindricalVectorField)
    return CylindricalScalarField(v.R, v.Z, v.Phi, v.VZ)
end
function get_VPhi(v::CylindricalVectorField)
    return CylindricalScalarField(v.R, v.Z, v.Phi, v.VPhi)
end


using Memoization
@memoize function RVpoloVPhi_pRpZ(v::CylindricalVectorField)
    R, Z = v.R, v.Z
    BR, BZ, BPhi = v.VR, v.VZ, v.VPhi

    A11 = BR./BPhi + R.*pRpZ( BR./BPhi,1,0,R,Z) ;
    A12 =            R.*pRpZ( BR./BPhi,0,1,R,Z) ;
    A21 = BZ./BPhi + R.*pRpZ( BZ./BPhi,1,0,R,Z) ;
    A22 =            R.*pRpZ( BZ./BPhi,0,1,R,Z) ;

    return A11, A12, A21, A22
end
@memoize function RVpoloVPhi_pRpZ_interp(v::CylindricalVectorField)
    R, Z, Phi = v.R, v.Z, v.Phi
    A11, A12, A21, A22 = RVpoloVPhi_pRpZ(v)
    A11_intp = linear_interpolation( (R,Z,Phi), A11 );
    A12_intp = linear_interpolation( (R,Z,Phi), A12 );
    A21_intp = linear_interpolation( (R,Z,Phi), A21 );
    A22_intp = linear_interpolation( (R,Z,Phi), A22 );

    return (r,z,phi) -> [A11_intp(r,z,phi)  A12_intp(r,z,phi); A21_intp(r,z,phi)  A22_intp(r,z,phi)]
end


using TensorCast
@memoize function RVpoloVPhi_pRpZ_delta_v_pert(v::CylindricalVectorField, v_pert::CylindricalVectorField)
    R, Z = v.R, v.Z
    BR, BZ, BPhi = v.VR, v.VZ, v.VPhi
    BR_pR = pRpZ( v.VR, 1, 0, R, Z);
    BR_pZ = pRpZ( v.VR, 0, 1, R, Z);
    BZ_pR = pRpZ( v.VZ, 1, 0, R, Z); 
    BZ_pZ = pRpZ( v.VZ, 0, 1, R, Z);
    BPhi_pR = pRpZ( v.VPhi, 1, 0, R, Z);
    BPhi_pZ = pRpZ( v.VPhi, 0, 1, R, Z);
    BR_pert = v_pert.VR
    BZ_pert = v_pert.VZ
    BPhi_pert = v_pert.VPhi
    BR_pert_pR = pRpZ( v_pert.VR, 1, 0, R, Z);
    BR_pert_pZ = pRpZ( v_pert.VR, 0, 1, R, Z);
    BZ_pert_pR = pRpZ( v_pert.VZ, 1, 0, R, Z);
    BZ_pert_pZ = pRpZ( v_pert.VZ, 0, 1, R, Z);
    BPhi_pert_pR = pRpZ( v_pert.VPhi, 1, 0, R, Z);
    BPhi_pert_pZ = pRpZ( v_pert.VPhi, 0, 1, R, Z);

    temp1 = BR_pert ./ BPhi - BR .* BPhi_pert ./ (BPhi.^2);
    temp2 = BR_pert_pR ./ BPhi - BR_pR ./ (BPhi.^2) .* BPhi_pert - BPhi_pR ./ (BPhi.^2) .* BR_pert - BR .* (  BPhi_pert_pR .* (BPhi.^2) - 2 * BPhi_pR .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A11_delta[iR,iZ,iPhi] := temp1[iR,iZ,iPhi] + R[iR] * temp2[iR,iZ,iPhi];

    temp1 = BZ_pert ./ BPhi - BZ .* BPhi_pert ./ (BPhi.^2);
    temp2 = BZ_pert_pR ./ BPhi - BZ_pR ./ (BPhi.^2) .* BPhi_pert - BPhi_pR ./ (BPhi.^2) .* BZ_pert - BZ .* (  BPhi_pert_pR .* (BPhi.^2) - 2 * BPhi_pR .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A21_delta[iR,iZ,iPhi] := temp1[iR,iZ,iPhi] + R[iR] * temp2[iR,iZ,iPhi];

    temp2 = BR_pert_pZ ./ BPhi - BR_pZ ./ (BPhi.^2) .* BPhi_pert - BPhi_pZ ./ (BPhi.^2) .* BR_pert - BR .* (  BPhi_pert_pZ .* (BPhi.^2) - 2 * BPhi_pZ .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A12_delta[iR,iZ,iPhi] :=                   + R[iR] * temp2[iR,iZ,iPhi];

    temp2 = BZ_pert_pZ ./ BPhi - BZ_pZ ./ (BPhi.^2) .* BPhi_pert - BPhi_pZ ./ (BPhi.^2) .* BZ_pert - BZ .* (  BPhi_pert_pZ .* (BPhi.^2) - 2 * BPhi_pZ .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A22_delta[iR,iZ,iPhi] :=                   + R[iR] * temp2[iR,iZ,iPhi];
    return A11_delta, A12_delta, A21_delta, A22_delta
end
@memoize function RVpoloVPhi_pRpZ_delta_v_pert_interp(v::CylindricalVectorField, v_pert::CylindricalVectorField)
    R, Z, Phi = v.R, v.Z, v.Phi

    A11_delta, A12_delta, A21_delta, A22_delta = RVpoloVPhi_pRpZ_delta_v_pert(v::CylindricalVectorField, v_pert::CylindricalVectorField)
    A11_delta_interp = linear_interpolation((R,Z,Phi), A11_delta);
    A12_delta_interp = linear_interpolation((R,Z,Phi), A12_delta);
    A21_delta_interp = linear_interpolation((R,Z,Phi), A21_delta);
    A22_delta_interp = linear_interpolation((R,Z,Phi), A22_delta);
    return (r,z,phi) -> [A11_delta_interp(r,z,phi) A12_delta_interp(r,z,phi); A21_delta_interp(r,z,phi) A22_delta_interp(r,z,phi);]
end


@memoize function divergence(v::CylindricalVectorField)
    dR = v.R[2]-v.R[1];
    dZ = v.Z[2]-v.Z[1];
    dPhi = v.Phi[2]-v.Phi[1];
    ans = Array{eltype(v.VR),3}(undef, length(v.R), length(v.Z), length(v.Phi) );
    ans[2:end-1,2:end-1,:] .= (v.VR[3:end,2:end-1,:].-v.VR[1:end-2,2:end-1,:])./(2*dR) .+ (v.VZ[2:end-1,3:end,:].-v.VZ[2:end-1,1:end-2,:])./(2*dZ);
    for iR in 2:length(v.R)-1
        ans[iR,:,2:end-1] .+= (v.VR[iR,:,2:end-1] .+ (v.VPhi[iR,:,3:end].-v.VPhi[iR,:,1:end-2])./(2*dPhi) ) ./ R[iR];
        ans[iR,:,1] .+= (v.VR[iR,:,1] .+ (v.VPhi[iR,:,2].-v.VPhi[iR,:,end-1])./(2*dPhi) ) ./ R[iR];
        ans[iR,:,end] .+= (v.VR[iR,:,1] .+ (v.VPhi[iR,:,2].-v.VPhi[iR,:,end-1])./(2*dPhi) ) ./ R[iR];
    end
    return CylindricalScalarField(v.R, v.Z, v.Phi, ans )
end

@memoize function magnitude(v::CylindricalVectorField)
    return CylindricalScalarField(v.R, v.Z, v.Phi, sqrt.(v.VR.^2 .+ v.VZ.^2 .+ v.VPhi.^2) )
end

@memoize function cross(v1::CylindricalVectorField, v2::CylindricalVectorField)
    # TODO: check R, Z, Phi grid are identical
    return CylindricalVectorField(
        v1.R, v1.Z, v1.Phi,
        v1.VPhi*v2.VZ - v1.VZ*v2.VPhi, 
        v1.VR*v2.VPhi - v1.VPhi*v2.VR,
        v1.VZ*v2.VR   - v1.VR*v2.VZ)
end

function directional_derivative_along_v_of_s(v::CylindricalVectorField, s::CylindricalScalarField)
    R, Z, Phi = v1.R, v1.Z, v1.Phi
    
    dR = R[2]-R[1]
    dZ = Z[2]-Z[1]
    dPhi = Phi[2]-Phi[1]
    scal = s.value
    ans = Array{eltype(v.VR),3}(undef, length(v.R), length(v.Z), length(v.Phi) );
    ans[2:end-1,2:end-1,:] .= v.VR[2:end-1,2:end-1,:] .* ( (scal[3:end,2:end-1,:].-scal[1:end-2,2:end-1,:])./(2*dR) );
    ans[2:end-1,2:end-1,:].+= v.VZ[2:end-1,2:end-1,:] .* ( (scal[2:end-1,3:end,:].-scal[2:end-1,1:end-2,:])./(2*dZ) );
    for iR in 2:length(v.R)-1
        ans[iR,:,2:end-1] .+=  v.VPhi[iR,:,2:end-1]./R[iR] .* (scal[iR,:,3:end].-scal[iR,:,1:end-2])./(2*dPhi)  ;
        ans[iR,:,1] .+= v.VPhi[iR,:,1]./R[iR] .* (scal[iR,:,2].-scal[iR,:,end-1])./(2*dPhi) ;
        ans[iR,:,end] .+= v.VPhi[iR,:,1]./R[iR]  .* (scal[iR,:,2].-scal[iR,:,end-1])./(2*dPhi) ;
    end
    return CylindricalScalarField(R, Z, Phi, ans)
end

function directional_derivative_along_v1_of_v2(v1::CylindricalVectorField, v2::CylindricalVectorField)
    """v1 is the direction, v2 is the field
    
    ```math
    \\vec{v}\\cdot\\nabla (u) = v_{R} \\dfrac{\\partial (u)}{\\partial R} + v_{Z} \\dfrac{\\partial (u)}{\\partial Z} +  \\dfrac{v_{\\phi}}{R} \\dfrac{\\partial (u)}{\\partial \\phi}

    \\vec{v}_{1}\\cdot\\nabla (\\vec{v}_{2}) = v_{1,R} \\dfrac{\\partial (\\vec{v}_{2})}{\\partial R} + v_{1,Z} \\dfrac{\\partial (\\vec{v}_{2})}{\\partial Z} +  \\dfrac{v_{1,\\phi}}{R} \\dfrac{\\partial (\\vec{v}_{2})}{\\partial \\phi}
    =v_{1,R}(\\dfrac{\\partial v_{2,R} }{\\partial R} \\hat{\\vec{e}}_{R} + \\dfrac{\\partial v_{2,Z} }{\\partial R} \\hat{\\vec{e}}_{Z}   +\\dfrac{\\partial v_{2,\\phi} }{\\partial R} \\hat{\\vec{e}}_{\\phi}    ) 
    +v_{1,Z}(\\dfrac{\\partial v_{2,R} }{\\partial Z} \\hat{\\vec{e}}_{R} + \\dfrac{\\partial v_{2,Z} }{\\partial Z} \\hat{\\vec{e}}_{Z}   +\\dfrac{\\partial v_{2,\\phi} }{\\partial Z} \\hat{\\vec{e}}_{\\phi}    ) 
    +\\dfrac{v_{1,\\phi}}{R}(\\dfrac{\\partial v_{2,R} }{\\partial \\phi} \\hat{\\vec{e}}_{R} + \\dfrac{\\partial v_{2,Z} }{\\partial \\phi} \\hat{\\vec{e}}_{Z}   +\\dfrac{\\partial v_{2,\\phi} }{\\partial \\phi} \\hat{\\vec{e}}_{\\phi}    ) 
    +\\dfrac{v_{1,\\phi}}{R}(v_{2,R} \\hat{\\vec{e}}_{\\phi} + \\qquad \\cdots\\qquad   -  v_{2,\\phi} \\hat{\\vec{e}}_{R}    ) 

    \\hat{\\vec{e}}_{R}= \\cos\\phi \\hat{\\vec{e}}_{x} + \\sin\\phi \\hat{\\vec{e}}_{y}

    \\hat{\\vec{e}}_{\\phi}= -\\sin\\phi \\hat{\\vec{e}}_{x} + \\cos\\phi \\hat{\\vec{e}}_{y}

    \\dfrac{\\partial \\hat{\\vec{e}}_{R} }{\\partial \\phi} = \\hat{\\vec{e}}_{\\phi} 

    \\dfrac{\\partial \\hat{\\vec{e}}_{\\phi} }{\\partial \\phi} = -\\hat{\\vec{e}}_{R} 
    ```
    """
    R, Z, Phi = v1.R, v1.Z, v1.Phi
    
    dR = R[2]-R[1]
    dZ = Z[2]-Z[1]
    dPhi = Phi[2]-Phi[1]
    ans_VR = Array{eltype(v1.VR),3}(undef, length(R), length(Z), length(Phi) );
    ans_VZ = Array{eltype(v1.VR),3}(undef, length(R), length(Z), length(Phi) );
    ans_VPhi = Array{eltype(v1.VR),3}(undef, length(R), length(Z), length(Phi) );
    ans_VR[2:end-1,2:end-1,:] .= v1.VR[2:end-1,2:end-1,:] .* ( (v2.VR[3:end,2:end-1,:].-v2.VR[1:end-2,2:end-1,:])./(2*dR) );
    ans_VZ[2:end-1,2:end-1,:] .= v1.VR[2:end-1,2:end-1,:] .* ( (v2.VZ[3:end,2:end-1,:].-v2.VZ[1:end-2,2:end-1,:])./(2*dR) );
    ans_VPhi[2:end-1,2:end-1,:] .= v1.VR[2:end-1,2:end-1,:] .* ( (v2.VPhi[3:end,2:end-1,:].-v2.VPhi[1:end-2,2:end-1,:])./(2*dR) );
    
    ans_VR[2:end-1,2:end-1,:] .+= v1.VZ[2:end-1,2:end-1,:] .* ( (v2.VR[2:end-1,3:end,:].-v2.VR[2:end-1,1:end-2,:])./(2*dZ) );
    ans_VZ[2:end-1,2:end-1,:] .+= v1.VZ[2:end-1,2:end-1,:] .* ( (v2.VZ[2:end-1,3:end,:].-v2.VZ[2:end-1,1:end-2,:])./(2*dZ) );
    ans_VPhi[2:end-1,2:end-1,:] .+= v1.VZ[2:end-1,2:end-1,:] .* ( (v2.VPhi[2:end-1,3:end,:].-v2.VPhi[2:end-1,1:end-2,:])./(2*dZ) );
    
    for iR in 2:length(R)-1
        # For the 0 < \phi < 2\pi/n sections
        ans_VR[iR,:,2:end-1] .+= v1.VPhi[iR,:,2:end-1]./R[iR] .* ( (v2.VR[iR,:,3:end].-v2.VR[iR,:,1:end-2])./(2*dPhi) ) ;
        ans_VZ[iR,:,2:end-1] .+= v1.VPhi[iR,:,2:end-1]./R[iR] .* ( (v2.VZ[iR,:,3:end].-v2.VZ[iR,:,1:end-2])./(2*dPhi) );
        ans_VPhi[iR,:,2:end-1] .+= v1.VPhi[iR,:,2:end-1]./R[iR] .* ( (v2.VPhi[iR,:,3:end].-v2.VPhi[iR,:,1:end-2])./(2*dPhi) );
        
        ans_VPhi[iR,:,2:end-1] .+= v1.VPhi[iR,:,2:end-1]./R[iR] .* v2.VR[iR,:,2:end-1];
        ans_VR[iR,:,2:end-1] .-= v1.VPhi[iR,:,2:end-1]./R[iR] .* v2.VPhi[iR,:,2:end-1];
        
        # For the \phi=0 section
        ans_VR[iR,:,1] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VR[iR,:,2].-v2.VR[iR,:,end-1])./(2*dPhi) );
        ans_VZ[iR,:,1] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VZ[iR,:,2].-v2.VZ[iR,:,end-1])./(2*dPhi) );
        ans_VPhi[iR,:,1] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VPhi[iR,:,2].-v2.VPhi[iR,:,end-1])./(2*dPhi) );
        
        ans_VPhi[iR,:,1] .+= v1.VPhi[iR,:,1]./R[iR] .* v2.VR[iR,:,1];
        ans_VR[iR,:,1] .-= v1.VPhi[iR,:,1]./R[iR] .* v2.VPhi[iR,:,1];
        
        # For the \phi=2\pi/n section
        ans_VR[iR,:,end] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VR[iR,:,2].-v2.VR[iR,:,end-1])./(2*dPhi) ) ;
        ans_VZ[iR,:,end] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VZ[iR,:,2].-v2.VZ[iR,:,end-1])./(2*dPhi) );
        ans_VPhi[iR,:,end] .+= v1.VPhi[iR,:,1]./R[iR] .* ( (v2.VPhi[iR,:,2].-v2.VPhi[iR,:,end-1])./(2*dPhi) );
        
        ans_VPhi[iR,:,end] .+= v1.VPhi[iR,:,1]./R[iR] .* v2.VR[iR,:,1];
        ans_VR[iR,:,end] .-= v1.VPhi[iR,:,1]./R[iR] .* v2.VPhi[iR,:,1];
    end
    
    return CylindricalVectorField(R,Z,Phi, ans_VR, ans_VZ, ans_VPhi)
end



# pick out the functions that are built-in
builtin_functions = Set([:eval, :include])

for name in names(@__MODULE__, all=true)
    if !(name in builtin_functions)
        @eval export $name
    end
end


end
