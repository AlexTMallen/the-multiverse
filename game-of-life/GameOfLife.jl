cd(@__DIR__)  # this sets the directory pointer to this file location
using Pkg  # Pkg.<command> is equivalent of `]` in command line julia
Pkg.activate("..")  # this specifies the environment being used and activates it


using Plots
using Random
using BenchmarkTools

seed = 633
Random.seed!(seed)

xdim, ydim = 320, 400
x = 1:xdim
y = 1:ydim
density = 0.1
cells = Random.rand(xdim, ydim) .< density

function update!(cells::BitMatrix, tmp::BitMatrix)
    xdim, ydim = size(cells)
    for j in 2:ydim - 1, i in 2:xdim - 1  # ignore edges
        neighbors = @view cells[i - 1:i + 1, j - 1:j + 1]
        num_neigh = sum(neighbors) - cells[i, j]
        if cells[i, j] 
            if num_neigh < 2 || num_neigh > 3
                tmp[i, j] = false
            else
                tmp[i, j] = true
            end
        elseif num_neigh == 3
            tmp[i, j] = true
        else
            tmp[i, j] = false
        end
    end
    tmp, cells
end

iterations = 1000
name = string(xdim, "x", ydim, " seed=", seed, " p=", density, " iters=", iterations)
tmp = similar(cells)

function frame!(cells::BitMatrix, tmp::BitMatrix, name::String)
    cells, tmp = update!(cells, tmp)
    # NOTE: x is vertical, y is horizontal, so swap the arguments for reasonable results
    heatmap(5:ydim - 5, 5:xdim - 5, view(cells, 5:xdim - 5, 5:ydim - 5),
        colorbar=false, grid=false, ticks=false,
        showaxis=false, title=name, aspect_ratio=1, dpi=150)
    cells, tmp
end

function animate!(cells::BitMatrix, tmp::BitMatrix, name::String)
    anim = @animate for i in 1:iterations
        cells, tmp = frame!(cells, tmp, name)
    end
    anim
end

anim = animate!(cells, tmp, name)
@btime frame!(cells, tmp, name)

gif(anim, string("gifs\\", name, ".gif"), fps=15)
