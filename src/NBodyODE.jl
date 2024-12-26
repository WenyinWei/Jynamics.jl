module NBodyODE

import LinearAlgebra: norm, cross, dot, I
import ..NBodyParticle: NaS

ϵ₀ = 8.854187817e-12   # 真空介电常数
μ₀ = 4π * 1e-7         # 真空磁导率 
c_real = 1 / sqrt(ϵ₀ * μ₀)  # 真空光速
c_initial = 3e8        # 初始光速，用于无延迟近似
c_transition_time = 1e-15  # 光速平滑过渡时间

function nbody_relativistic!(du::Matrix{Float64}, u::Matrix{Float64}, h, p, t)
    EBfield, species, u0 = p
    alive_mask = [s != NaS for s in species]
    N = length(alive_mask)
    N_alive = count(alive_mask)

    # 提取状态矩阵和延迟时间矩阵
    r_γv = view(u, 1:N, 1:6)
    delay_times = view(u, 1:N, 6+1:6+N)

    du[:] .= 0.0
    dr_γv_dt = view(du, 1:N, 1:6)
    delay_times_dt = view(du, 1:N, 6+1:6+N)

    # 当前光速随时间平滑过渡
    # c = c_real + (c_initial - c_real) * exp(-t / c_transition_time)
    # dc_dt = - (c_initial - c_real) / c_transition_time * exp(-t / c_transition_time)

    alive_indices = findall(alive_mask)
    for i in alive_indices
        species_i = species[i]
        m0_i = species_i.mass
        q_i = species_i.charge

        r_i = r_γv[i, 1:3]
        γv_i = r_γv[i, 4:6]

        γ_i = sqrt(1 + norm(γv_i / c_real)^2)
        v_i = γv_i / γ_i

        E_total, B_total = EBfield(r_i, t)
        for j in alive_indices
            if i == j
                continue  # 跳过对角元
            end

            species_j = species[j]
            m0_j = species_j.mass
            q_j = species_j.charge

            t_ji = delay_times[j, i] # the time particle j imposes its field on i, that is t_{j->i}
            r_j = h(p, t_ji, idxs=[j, j+N, j+2N])
            γv_j = h(p, t_ji, idxs=[j+3N, j+4N, j+5N])
            dγv_j_dt = h(p, t_ji, Val{1}, idxs=[j+3N, j+4N, j+5N])
            println("r_j = ", r_j)
            println("γv_j = ", γv_j)
            println("dγv_j_dt = ", dγv_j_dt)
            println("t_ji = ", t_ji, " i = ", i, " j = ", j)

            γ_j = sqrt(1 + norm(γv_j / c_real)^2)
            v_j = γv_j / γ_j

            r_ji = r_i - r_j
            distance_ji = norm(r_ji)
            hatr_ji = r_ji / distance_ji

            beta_j = v_j / c_real
            beta_j_dt = 1/m0_j * inv(γ_j^3 * v_j * beta_j' + γ_j * c_real * I(3)) * dγv_j_dt 

            hatr_dot_beta = dot(hatr_ji, beta_j)
            one_minus_hatr_dot_beta = (1 - hatr_dot_beta)

            E_ji = q_j / (4π * ϵ₀) * (
                (hatr_ji - beta_j) / (γ_j^2 * distance_ji^2)
                + cross(hatr_ji, cross(hatr_ji - beta_j, beta_j_dt))  / (c_real * distance_ji)
                ) / (one_minus_hatr_dot_beta^3)
            B_ji = cross(hatr_ji, E_ji) / c_real

            E_total += E_ji
            B_total += B_ji

            numerator = 1 - dot(hatr_ji, v_i) / c_real
            denominator = 1 - dot(hatr_ji, v_j) / c_real

            # delay_times_dt[j, i] = (numerator + dc_dt * distance_ji / c^2) / denominator
            delay_times_dt[j, i] = (numerator) / denominator
        end

        dr_γv_dt[i, 1:3] .= v_i
        dr_γv_dt[i, 4:6] .= (q_i / m0_i) * (E_total + cross(v_i, B_total))
    end
    
end




function nbody_h(p, t, deriv::Union{Nothing, Type{Val{1}}} = nothing; idxs::Union{Nothing, Vector{Int}} = nothing)
    EBfield, species, u0 = p
    N = length(species)
    if idxs === nothing
        if deriv === Val{1}
            u0_new = copy(u0)
            u0_new[1:N, 4:6] .= 0.0
            u0_new[1:N, 6+1:6+N] .= 1.0
            return u0_new
        else
            return u0
        end
    else
        if deriv === Val{1}
            return [idx <= 6N ? (idx <= 3N ? u0[idx+3N] : 0.0) : 1.0 for idx in idxs]
        else
            return [u0[idx] for idx in idxs] 
        end
    end   
    
end



for name in names(@__MODULE__; all=true)
    if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $name
    end
end



end
