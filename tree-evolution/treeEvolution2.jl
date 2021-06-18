"""
As opposed to treeEvolution.jl, this simulation gives trees a trunk height 
seperate from the remaining height, and the trees start with 0 trunk height
"""


cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")

using Plots
using Random
using Statistics

seed = rand(1:1000)
Random.seed!(seed)
println(string("seed: ", seed))

abstract type Tree end

mutable struct AltruisticTree <: Tree
    id::Int32
    height::Int32 # = trunk_height + top_heights
    top_height::Int32 
    donation_rate::Float32  # x
    growth_cost::Float32  # y
    growth_gain::Float32  # z
    energy::Float32
    prev_energy::Float32  # memory of how much energy I had 1 unit of time ago
    age::Int64
end

modp1(x::Integer, m::Integer) = x < 1 ? m + x % m : 1 + (x - 1) % m


function continuous_time_evolve_no_animate!(trees::Vector{AltruisticTree}, occ_idxs::Set{Int64}, free_idxs::Set{Int64}, generations::Int64, age_cap::Int64=100, spawnrate::Float64=0.01)
    for gen in 1:generations
        
        # update energies and ages
        for idx in occ_idxs
            tree = trees[idx]
            idx_left = get_left_neigh_idx(trees, occ_idxs, idx)
            neigh_left = trees[idx_left]
            idx_right = get_right_neigh_idx(trees, occ_idxs, idx)
            neigh_right = trees[idx_right]

            if tree.energy > tree.prev_energy
                tree.height = min(45, tree.height + 1)
                tree.prev_energy = tree.energy
            end
            tree.energy = tree.growth_gain * (tree.height - neigh_left.height / modp1(idx - idx_left, length(trees)) + 10) - tree.growth_cost * tree.height
            offset = 0  # diff between tree heights at which altruistic trees will start giving #TODO make this evolve
            if neigh_right.height > tree.height + offset
                tree.energy += neigh_right.donation_rate * (neigh_right.height - tree.height + offset)
                # TODO decrease this tree so that its trunk height is its neighbor's height
            end
            if tree.height > neigh_left.height + offset
                tree.energy -= tree.donation_rate * (tree.height - neigh_left.height + offset)
                # TODO decrease height of tree --I think this needs a trunk_height parameter
            end
            tree.age += 1
        end

        # birth and death
        for idx in occ_idxs  # concurrent modification seems to be ok in this case
            tree = trees[idx]
            if tree.energy <= 0 || tree.age > age_cap  # death
                pop!(occ_idxs, idx)
                push!(free_idxs, idx)
            elseif length(free_idxs) > 0 && rand() < spawnrate * tree.energy  # birth
                new_x = max(0, tree.donation_rate + randn() / 10)
                new_tree = AltruisticTree(tree.id, tree.top_height, tree.top_height, new_x, tree.growth_cost, tree.growth_gain, 0, -1, 0)
                new_idx = rand(free_idxs)
                pop!(free_idxs, new_idx)
                push!(occ_idxs, new_idx)
                trees[new_idx] = new_tree
            end
        end
    end

    heights = Vector{Float32}()
    energies = Vector{Float32}()
    xs = Vector{Float32}()
    count = 0
    for idx in occ_idxs
        tree = trees[idx]
        push!(heights, tree.height)
        push!(energies, tree.energy)
        push!(xs, tree.donation_rate)
    count += 1
    end
    heights, energies, xs, count
end

function get_left_neigh_idx(trees::Vector{AltruisticTree}, occ_idxs::Set{Int64}, idx::Int64)
    curr = modp1(idx - 1, length(trees))
    while !(curr in occ_idxs)
        curr =  modp1(curr - 1, length(trees))
    end
    curr
end

function get_right_neigh_idx(trees::Vector{AltruisticTree}, occ_idxs::Set{Int64}, idx::Int64)
    curr = modp1(idx + 1, length(trees))
    while !(curr in occ_idxs)
        curr =  modp1(curr + 1, length(trees))
    end
    curr
end

x, y, z = 1, 3, 4
top_height = 5  # TODO: to make it more realistic, this should be the maximum height difference that matters
num_trees = 200
xdim = 4 * num_trees
ydim = 50  # ⟹ 45 is max height
grid = zeros(Int32, ydim, xdim)
grid[1, :] .= -2 * num_trees  # ground
grid[ydim - 3:ydim - 1, 2:2 + num_trees ÷ 7] .= 2 * num_trees  # sunshine!
generations = 50_000
age_cap = 100
spawnrate = 0.1
name = string("lifecycle--gens=", generations, " trees=", num_trees, " x,y,z=", x, ",", y, ",", z,
              " sr=", spawnrate, " ac=", age_cap, " seed=", seed)
    
function init_test_x()
    n_trials = 10
    n_quantiles = 5
    std = 0.1
    d = 1 / (n_quantiles + 1)
    quantiles = d:d:1 - d
    x_inits = [0:0.25:2; 5; 10; 100]
    s = (length(x_inits), n_quantiles)
    q_heights = Matrix{Float32}(undef, s)
    q_energies = Matrix{Float32}(undef, s)
    q_end_xs = Matrix{Float32}(undef, s)
    q_counts = Matrix{Float32}(undef, s)

    for pair in enumerate(x_inits)
        i, x = pair
        println(x)
        heights = Vector{Float32}()
        energies = Vector{Float32}()
        end_xs = Vector{Float32}()
        counts = Vector{Float32}()
        for trial in 1:n_trials
            trees = Vector{AltruisticTree}(undef, 2 * num_trees)
            perm = randperm(length(trees))
            idxs = Set(perm[1:num_trees])
            free_idxs = Set(perm[num_trees + 1:end])
            for i in idxs
                age = rand(0:age_cap - 1)
                height = min(45, age + top_height)
                trees[i] = AltruisticTree(-i - 10, height, top_height, max(x + randn() * std, 0), y, z, 0, - 1, age)
            end
            hs, es, xs, count = continuous_time_evolve_no_animate!(trees, idxs, free_idxs, generations, age_cap, spawnrate)
            push!(heights, mean(hs))
            push!(energies, mean(es))
            push!(end_xs, mean(xs))
            push!(counts, count)
        end
        sort!(heights)
        sort!(energies)
        sort!(end_xs)
        sort!(counts)
        for pair in enumerate(quantiles)
            j, q = pair
            q_heights[i, j] = get_quantile(heights, q)
            q_energies[i, j] = get_quantile(energies, q)
            q_end_xs[i, j] = get_quantile(end_xs, q)
            q_counts[i, j] = get_quantile(counts, q)

        end
    end
    labels = reshape([string("x--", round(q, digits=2)) for q in quantiles], 1, n_quantiles)
    p = plot(x_inits, q_end_xs, label=labels, xaxis="initial x", legend=:top)
    savefig(p, string("gifs\\", name, ".png"))
    p, q_heights, q_energies, q_end_xs, q_counts
end

get_quantile(list::Vector{Float32}, q::Float64) = list[Int(floor(length(list) * q + 1))]

# anim, avg_energies, avg_heights, counts, avg_xs = continuous_time_init()    
# g = gif(anim, string("gifs\\", name, ".gif"), fps=15)
# p = plot(avg_xs, lab="x")
# plot!(p, counts ./ maximum(counts), lab="# trees")
# plot!(p, avg_energies ./ maximum(avg_energies), lab="energy/tree")
# plot!(p, avg_heights ./ maximum(avg_heights), lab="height", xaxis="time", legend=:top, ylims=(-0.1, maximum(avg_xs) + 0.1))
# savefig(p, string("gifs\\", name, ".png"))
p, q_heights, q_energies, q_end_xs, q_counts = init_test_x()
p