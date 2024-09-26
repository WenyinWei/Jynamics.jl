module File 

function VectorVector2Array(vec1_of_vec2)
    return permutedims( reduce(hcat, vec1_of_vec2), (2,1) );
end
function VectorMatrix2Array(vec_of_mat)
    vec_size = size(vec_of_mat)[1]
    mat_size = size(vec_of_mat[1])
    arr3D = Array{Float64,3}(undef, vec_size, mat_size[1], mat_size[2])
    for i in 1:mat_size[1]
        for j in 1:mat_size[2]
            arr3D[:,i,j] = [mat[i,j]  for mat in vec_of_mat ]
        end
    end
    return arr3D
end

# export all
for name in names(@__MODULE__; all=true)
    if Base.isidentifier(name) && name ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $name
    end
end

end