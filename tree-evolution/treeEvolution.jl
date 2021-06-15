using Base:num_bit_chunks
""" 
Tree evolution
Author: Alex Mallen (atmallen@uw.edu)

Can trees grow down?
    
Trees are tall because they compete for light.
    In a forest where all trees are at least 10 meters tall, with 8 meters of trunk, each tree
    could be 8 meters shorter and recieve the same amount of light, but
    spend much less energy on growth.

    However, the dominant strategy, paradoxically, is to grow tall because regardless
    of how tall the neighboring trees are, it is more advantageous to grow tall. This
    is because for every millimeter of height gained by the tree, y units of energy
    are used to grow, while z units of energy are gained due to increased light,
    where y < z.
    
    If trees were rational human beings, they could all be better off by
    coming to an agreement that changes their dominant strategy: A tree gives its
    neighbor (in the direction of the sun) x units of energy for every millimeter shorter it becomes, where
    z - y < x < y < z. The reason x < y is because the x delivered to one tree
    must come from the y gained by its neighbor which is now able to be shorter

    z - y < y ⟹ y > z/2 is one important constraint on the system.

    More precisely, "advantage" for a tree (or its genes) is defined as likelihood 
    of persisting long-term, and we assume that energy is a good proxy for this.

    (Is this green-beard altruism?) I don't think so because it doesn't distinguish like from unlike
"""

cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")

using Plots
using Random
using BenchmarkTools

seed = 633
Random.seed!(seed)

abstract type Tree end

mutable struct NormalTree <: Tree
    id::Int32
    height::Int32
    growth_cost::Float32  # y
    growth_gain::Float32  # z
    energy::Float32
end

mutable struct AltruisticTree <: Tree
    id::Int32
    height::Int32
    donation_rate::Float32  # x
    growth_cost::Float32  # y
    growth_gain::Float32  # z
    energy::Float32
end

modp1(x::Integer, m::Integer) = x < 1 ? m + x % m : 1 + (x - 1) % m

function both_animate!(trees::Vector{Tree}, grid::Matrix{Int32}, trees_cache::Vector{Tree}, generations::Int64, name::String)
    ydim, xdim = size(grid)
    avg_energies = Vector{Float32}()
    avg_heights = Vector{Float32}()
    anim = @animate for gen in 1:generations
        # first update each tree's energy
        avg_energy = 0
        for i in 1:length(trees)
            neigh1, neigh2 = trees[modp1(i - 1, length(trees))], trees[modp1(i + 1, length(trees))]  # neighbor in the direction of the sun
            tree = trees[i]
            tree.energy = tree.growth_gain * (tree.height - neigh1.height + 45) - tree.growth_cost * tree.height
            if typeof(neigh2) === AltruisticTree
                tree.energy += neigh2.donation_rate * (neigh2.height - tree.height)  # TODO this 45 ÷ 2 might need to be changed to avg_height or neigh2.height
            end
            avg_energy += tree.energy
            if tree.energy <= 0
                println(tree.energy)
            end
        end
        avg_energy = round(avg_energy / length(trees), digits=1)
        push!(avg_energies, avg_energy)

        avg_height = 0
        for i in 1:length(trees)
            # draw trees[i] onto current grid
            # grid[2:trees[i].height + 1, 4 * (i - 1) + 2(gen % 2 + 1)] = trees[i].id
            grid[2:ydim - 4, 4 * (i - 1) + 2] .= 0
            grid[2:trees[i].height + 1, 4 * (i - 1) + 2] .= trees[i].id
            avg_height += trees[i].height

            # now update trees
            j1, j2 = modp1(i - 1, length(trees)), modp1(i + 1, length(trees))
            trees_cache[i] = rand() < trees[j1].energy / (trees[j1].energy + trees[j2].energy) ? trees[j1] : trees[j2]
            trees_cache[i].height = min(max(trees_cache[i].height + rand(-1:1), 1), 45)  # random mutation: height must be between 1 and 45
        end
        avg_height = round(avg_height / length(trees), digits=1)
        push!(avg_heights, avg_height)
        trees, trees_cache = trees_cache, trees
        # plot
        heatmap(1:xdim, 1:ydim, grid,
                colorbar=false, grid=false, ticks=false,
                showaxis=false, title=string(gen, " ", avg_height, " both; ", name), dpi=150)
    end
    anim, avg_energies, avg_heights
end

# 1 normal trees
# 2 altruistic trees  (this by itself won't be interesting)
# 3 normal vs altruistic trees
x, y, z = 1.5, 3, 4
num_trees = 125
xdim = 4 * num_trees
ydim = 50
grid = zeros(Int32, ydim, xdim)
grid[1, :] .= -2 * num_trees  # ground
grid[ydim - 3:ydim - 1, 2:2 + num_trees ÷ 7] .= 2 * num_trees  # sunshine!
generations = 200
ratio = 0
name = string("gens=", generations, " trees=", num_trees, " x,y,z=", x, ",", y, ",", z, " ratio=", ratio)
    
function both_trees_init()
    trees = Vector{Tree}()  # from l to r
    for i in 1:num_trees
        height = rand(10:20)
        if rand() < ratio
            push!(trees, NormalTree(i + 10, height, y, z, 0))
        else
            push!(trees, AltruisticTree(-i - 10, height, x, y, z, 0))
        end
    end

    results = both_animate!(trees, grid, similar(trees), generations, name)
results
end

anim, avg_energies, avg_heights = both_trees_init()    
gif(anim, string("gifs\\", name, ".gif"), fps=3)
p = plot(avg_energies / maximum(avg_energies), lab="energy")
plot!(p, avg_heights / maximum(avg_heights), lab="height", xaxis="time", legend=:top)
savefig(p, string("gifs\\", name, ".png"))
p