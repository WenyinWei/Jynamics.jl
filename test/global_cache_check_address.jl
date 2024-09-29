using Jynamics
using Test

@testset "Global cache aided by Memoization" begin
    rminb = 4.0;
    rmaxb = 7.0;
    zminb = -1.5;
    zmaxb = 1.5;
    nSym = 5;

    R = collect( range(rminb, stop=rmaxb, length=128) );
    Z = collect( range(zminb, stop=zmaxb, length=128) );
    Phi = collect( range(0.0, stop=2pi/nSym, length=65) );

    BR = Array{Float64,3}(undef, 128, 128, 65);
    BZ = Array{Float64,3}(undef, 128, 128, 65);
    BPhi = Array{Float64,3}(undef, 128, 128, 65);

    fieldB = CylindricalVectorField(R, Z, Phi, BR, BZ, BPhi,  nSym);
    interp1 = RVpoloVPhi_pRpZ_interp(fieldB);
    interp2 = RVpoloVPhi_pRpZ_interp(fieldB);
    
    A11, A12, A21, A22 = RVpoloVPhi_pRpZ(fieldB)
    A11_, A12_, A21_, A22_ = RVpoloVPhi_pRpZ(fieldB)
    @test pointer_from_objref(A11) == pointer_from_objref(A11_)
    @test pointer_from_objref(A12) == pointer_from_objref(A12_)
    @test pointer_from_objref(A21) == pointer_from_objref(A21_)
    @test pointer_from_objref(A22) == pointer_from_objref(A22_)
    @test interp1 === interp2 # immutable objects can not be fetched addressed by pointer_from_objref  
end
