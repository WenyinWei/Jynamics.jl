module NBodyParticle

# 定义粒子种类和属性
struct ParticleSpecies
    name::String
    mass::Float64    # 静止质量，单位 kg
    charge::Float64  # 电荷量，单位 C
end

# 定义粒子
NaS = ParticleSpecies("NaS", NaN, NaN)
proton = ParticleSpecies("Proton", 1.6726219e-27, 1.602176634e-19)
deuteron = ParticleSpecies("Deuteron", 3.3435837724e-27, 1.602176634e-19)
triton = ParticleSpecies("Triton", 5.0073567446e-27, 1.602176634e-19)
electron = ParticleSpecies("Electron", 9.10938356e-31, -1.602176634e-19)

species_list = [NaS, proton, deuteron, triton, electron]






for name in names(@__MODULE__; all=true)
    if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $name
    end
end



end
