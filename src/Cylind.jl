module Cylind

export CylindricalScalarField, CylindricalVectorField, divergence, magnitude, cross, directional_derivative_along_v_of_s, directional_derivative_along_v1_of_v2

struct CylindricalScalarField{T}
    R::Vector{T}
    Z::Vector{T}
    Phi::Vector{T}
    value::Array{T,3}
end
function Base.:+(s1::CylindricalScalarField, s2::CylindricalScalarField)
    return CylindricalVectorField(s.R, s.Z, s.Phi, s1.value + s2.value)
end
function Base.:*(a::Number, s::CylindricalScalarField)
    return CylindricalVectorField(s.R, s.Z, s.Phi, a * s.value)
end

struct CylindricalVectorField{T}
    R::Vector{T}
    Z::Vector{T}
    Phi::Vector{T}
    VR::Array{T,3}
    VZ::Array{T,3}
    VPhi::Array{T,3}
end
function Base.:+(v1::CylindricalVectorField, v2::CylindricalVectorField)
    return CylindricalVectorField(v1.R, v1.Z, v1.Phi, v1.VR + v2.VR, v1.VZ + v2.VZ, v1.VPhi + v2.VPhi)
end
function Base.:*(a::Number, v::CylindricalVectorField)
    return CylindricalVectorField(v.R, v.Z, v.Phi, a * v.VR, a * v.VZ, a * v.VPhi)
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

@memoize function pRpZ(s::CylindricalScalarField, Rord, Zord)
    dR = s.R[2]-s.R[1]
    dZ = s.Z[2]-s.Z[1]

    if Rord==0 && Zord==0
        return s
    elseif Rord > 0
        lastord_field = pRpZ(s, Rord-1, Zord).value
        thisord_field = similar(s.value)
        thisord_field[:,:,:] .= NaN
        thisord_field[2:end-1,:,:] = (lastord_field[3:end,:,:]- lastord_field[1:end-2,:,:]) / (2dR)
        return CylindricalScalarField( s.R, s.Z, s.Phi, thisord_field )
    elseif Zord > 0 
        lastord_field = pRpZ(s, Rord, Zord-1).value
        thisord_field = similar(s.value)
        thisord_field[:,:,:] .= NaN
        thisord_field[:,2:end-1,:] = (lastord_field[:,3:end,:]- lastord_field[:,1:end-2,:]) / (2dZ)
        return CylindricalScalarField( s.R, s.Z, s.Phi, thisord_field )
    end
end
using Interpolations
@memoize function pRpZ_interp(s::CylindricalScalarField, Rord, Zord)
    return linear_interpolation( 
            (s.R, s.Z, s.Phi), pRpZ(s, Rord, Zord).value )
end

using TensorCast
function A()
    
    BR_pR = pRpZ(BR,1,0);
    BR_pZ = pRpZ(BR,0,1);
    BZ_pR = pRpZ(BZ,1,0);
    BZ_pZ = pRpZ(BZ,0,1);
    BPhi_pR = pRpZ(BPhi,1,0);
    BPhi_pZ = pRpZ(BPhi,0,1);
    BR_pert_pR = pRpZ(BR_pert,1,0);
    BR_pert_pZ = pRpZ(BR_pert,0,1);
    BZ_pert_pR = pRpZ(BZ_pert,1,0);
    BZ_pert_pZ = pRpZ(BZ_pert,0,1);
    BPhi_pert_pR = pRpZ(BPhi_pert,1,0);
    BPhi_pert_pZ = pRpZ(BPhi_pert,0,1);

    temp1 = BR_pert ./ BPhi - BR .* BPhi_pert ./ (BPhi.^2);
    temp2 = BR_pert_pR ./ BPhi - BR_pR ./ (BPhi.^2) .* BPhi_pert - BPhi_pR ./ (BPhi.^2) .* BR_pert - BR .* (  BPhi_pert_pR .* (BPhi.^2) - 2 * BPhi_pR .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A11_delta_DeltaB[iR,iZ,iPhi] := temp1[iR,iZ,iPhi] + R[iR] * temp2[iR,iZ,iPhi];

    temp1 = BZ_pert ./ BPhi - BZ .* BPhi_pert ./ (BPhi.^2);
    temp2 = BZ_pert_pR ./ BPhi - BZ_pR ./ (BPhi.^2) .* BPhi_pert - BPhi_pR ./ (BPhi.^2) .* BZ_pert - BZ .* (  BPhi_pert_pR .* (BPhi.^2) - 2 * BPhi_pR .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A21_delta_DeltaB[iR,iZ,iPhi] := temp1[iR,iZ,iPhi] + R[iR] * temp2[iR,iZ,iPhi];

    temp2 = BR_pert_pZ ./ BPhi - BR_pZ ./ (BPhi.^2) .* BPhi_pert - BPhi_pZ ./ (BPhi.^2) .* BR_pert - BR .* (  BPhi_pert_pZ .* (BPhi.^2) - 2 * BPhi_pZ .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A12_delta_DeltaB[iR,iZ,iPhi] :=                   + R[iR] * temp2[iR,iZ,iPhi];

    temp2 = BZ_pert_pZ ./ BPhi - BZ_pZ ./ (BPhi.^2) .* BPhi_pert - BPhi_pZ ./ (BPhi.^2) .* BZ_pert - BZ .* (  BPhi_pert_pZ .* (BPhi.^2) - 2 * BPhi_pZ .* BPhi .* BPhi_pert ) ./ (BPhi.^4);
    @cast A22_delta_DeltaB[iR,iZ,iPhi] :=                   + R[iR] * temp2[iR,iZ,iPhi];


    A11_delta_DeltaB_interp = linear_interpolation((R,Z,Phi), A11_delta_DeltaB);
    A12_delta_DeltaB_interp = linear_interpolation((R,Z,Phi), A12_delta_DeltaB);
    A21_delta_DeltaB_interp = linear_interpolation((R,Z,Phi), A21_delta_DeltaB);
    A22_delta_DeltaB_interp = linear_interpolation((R,Z,Phi), A22_delta_DeltaB);
    @inline A_delta_DeltaB(r,z,phi) = [A11_delta_DeltaB_interp(r,z,phi) A12_delta_DeltaB_interp(r,z,phi); A21_delta_DeltaB_interp(r,z,phi) A22_delta_DeltaB_interp(r,z,phi);]
end


function divergence(v::CylindricalVectorField)
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

function magnitude(v::CylindricalVectorField)
    return CylindricalScalarField(v.R, v.Z, v.Phi, sqrt.(v.VR.^2 .+ v.VZ.^2 .+ v.VPhi.^2) )
end

function cross(v1::CylindricalVectorField, v2::CylindricalVectorField)
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

end
