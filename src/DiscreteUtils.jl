module DiscreteUtils


function find_Delta_theta(scatxy, axisxy)
    Nmap = length(scatxy[1,:])-1

    init_rel_xy = scatxy[:,1] - axisxy
    init_theta = atan(init_rel_xy[2], init_rel_xy[1])
    last_theta = init_theta
    last_turn = 0 
    Delta_theta_upper_bound = +Inf
    Delta_theta_lower_bound = -Inf
    for imap in 1:Nmap
        rel_xy = scatxy[:,imap+1] - axisxy
        this_theta = atan(rel_xy[2], rel_xy[1])
        this_turn = this_theta<last_theta ? last_turn+1 : last_turn # how many times it passed the semi-axis from axisxy to axisxy-(infty,0).
        if this_theta < init_theta
            Delta_theta_upper_bound = min(Delta_theta_upper_bound, 2pi*(this_turn-1+1)/imap )
            Delta_theta_lower_bound = max(Delta_theta_lower_bound, 2pi*(this_turn-1)/imap )
        else
            Delta_theta_upper_bound = min(Delta_theta_upper_bound, 2pi*(this_turn+1)/imap )
            Delta_theta_lower_bound = max(Delta_theta_lower_bound, 2pi*(this_turn)/imap )
        end
        last_theta = this_theta
        last_turn = this_turn
    end
    
    if Delta_theta_upper_bound < Delta_theta_lower_bound
        return NaN
    else
        return (Delta_theta_upper_bound + Delta_theta_lower_bound) / 2
    end
end

for name in names(@__MODULE__; all=true)
    if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $name
    end
end



end
