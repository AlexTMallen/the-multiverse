using Base:num_bit_chunks
""" 
Tree evolution
Author: Alex Mallen (atmallen@uw.edu)

Can trees grow down?
    
Trees are tall because they compete for light.
    In a forest where all trees are at least 10 meters tall, each tree
    could be 10 meters shorter and recieve the same amount of light, but
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

    (Is this green-beard altruism?)
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

struct NormalTree
    id::Int8
    height::UInt16
    growthCost::Float32  # y
    growthGain::Float32  # z
    energy::Float32
end

struct AltruisticTree
    id::Int8
    height::UInt16
    growthCost::Float32  # y
    growthGain::Float32  # z
    donationRate::Float32  # x
    energy::Float32
end

modp1(x::Integer, m::Integer) = x < 1 ? m + x % m : 1 + (x - 1) % m

function normal_animate!(trees::Vector{NormalTree}, grid::Matrix{Int8}, trees_cache::Vector{NormalTree}, generations::Int64, name::String)
    ydim, xdim = size(grid)
    anim = @animate for gen in 1:generations
        for i in 1:length(trees)
            # draw trees[i] onto current grid
            # grid[2:trees[i].height + 1, 4 * (i - 1) + 2(gen % 2 + 1)] = i
            grid[2:ydim - 4, 4 * (i - 1) + 2] .= 0
            grid[2:trees[i].height + 1, 4 * (i - 1) + 2] .= trees[i].id

            # now update trees
            j1, j2 = modp1(i - 1, length(trees)), modp1(i + 1, length(trees))
            trees_cache[i] = rand() < trees[j1].energy / (trees[j1].energy + trees[j2].energy) ? trees[j1] : trees[j2]
        end
        trees, trees_cache = trees_cache, trees
        # plot
        heatmap(1:xdim, 1:ydim, grid,
                colorbar=false, grid=false, ticks=false,
                showaxis=false, title=string("normal ", name), aspect_ratio=1, dpi=150)
    end
    anim
end

# 1 normal trees
# 2 altruistic trees
# 3 normal vs altruistic trees
x, y, z = 2, 3, 4
num_trees = 15
xdim = 4 * num_trees
ydim = 50
grid = zeros(Int8, ydim, xdim)
grid[1, :] .= -num_trees  # ground
grid[ydim - 3:ydim - 1, 2:4] .= 2 * num_trees  # sunshine!
generations = 20
name = string("gens=", generations, " x,y,z=", x, ",", y, ",", z)


function normal_trees_init()
    trees = Vector{NormalTree}()  # from l to r
    for i in 1:num_trees
        height = rand(1:40)
        push!(trees, NormalTree(i, height, y, z, height * (y - z)))
    end

    anim = normal_animate!(trees, grid, similar(trees), generations, name)
    anim
end

anim = normal_trees_init()    
gif(anim, string("gifs\\", name, ".gif"), fps=2)