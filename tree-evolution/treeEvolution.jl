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
    The tree actually needs to shrink when it donates
    The difference between economics and evolution is that capitalists are profit maximizing while genes are *relative* profit maximizing
"""

cd(@__DIR__)
cd("..")
using Pkg
Pkg.activate(".")

using Statistics
using Plots
using Random

seed = rand(1:1000)
Random.seed!(seed)
println(string("seed", seed))

abstract type Tree end

mutable struct AltruisticTree <: Tree
    id::Int32
    height::Int32
    donation_rate::Float32  # x
    growth_cost::Float32  # y
    growth_gain::Float32  # z
    energy::Float32
    age::Int64
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
            if typeof(neigh2) === AltruisticTree && neigh2.height > tree.height
                tree.energy += neigh2.donation_rate * (neigh2.height - tree.height)  # TODO this 45 ÷ 2 might need to be changed to avg_height or neigh2.height
            end
            if typeof(tree) == AltruisticTree && tree.height > neigh1.height
                tree.energy -= tree.donation_rate * (tree.height - neigh1.height)  # TODO this 45 ÷ 2 might need to be changed to avg_height or neigh2.height
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

function continuous_time_evolve!(trees::Vector{AltruisticTree}, occ_idxs::Set{Int64}, free_idxs::Set{Int64}, generations::Int64, age_cap::Int64=100, spawnrate::Float64=0.01)
    avg_energies = Vector{Float32}()
    avg_heights = Vector{Float32}()
    avg_xs = Vector{Float32}()
    counts = Vector{Int64}()

        anim = @animate for gen in 1:generations
        heights = Vector{Float32}()
        energies = Vector{Float32}()
        xs = Vector{Float32}()
        avg_height = 0
        avg_energy = 0
        count = 0
        avg_x = 0
        
        # update energies and ages
        for idx in occ_idxs
            tree = trees[idx]
            idx_left = get_left_neigh_idx(trees, occ_idxs, idx)
            neigh_left = trees[idx_left]
            idx_right = get_right_neigh_idx(trees, occ_idxs, idx)
            neigh_right = trees[idx_right]
            tree.energy = tree.growth_gain * (tree.height - neigh_left.height / modp1(idx - idx_left, length(trees)) + 10) - tree.growth_cost * tree.height
            offset = 0  # diff between tree heights at which altruistic trees will start giving
            if neigh_right.height > tree.height + offset
                tree.energy += neigh_right.donation_rate * (neigh_right.height - tree.height + offset)  # TODO this 45 ÷ 2 might need to be changed to avg_height or neigh2.height
            end
            if tree.height > neigh_left.height + offset
                tree.energy -= tree.donation_rate * (tree.height - neigh_left.height + offset)  # TODO this 45 ÷ 2 might need to be changed to avg_height or neigh2.height
            end
            count += 1
            tree.age += 1
            push!(heights, tree.height)
            push!(xs, tree.donation_rate)
            push!(energies, tree.energy)
            avg_height += tree.height
            avg_energy += tree.energy
            avg_x += tree.donation_rate
        end
        
        avg_height = round(avg_height / length(occ_idxs), digits=1)
        push!(avg_heights, avg_height)
        avg_energy = round(avg_energy / length(occ_idxs), digits=1)
        push!(avg_energies, avg_energy)
        avg_x = round(avg_x / length(occ_idxs), digits=2)
        push!(avg_xs, avg_x)
        push!(counts, count)
        
        print("correlation between height and energy: ")
        println(round(cor(heights, energies), digits=3))
        print("between x and energy: ")
        println(round(cor(xs, energies), digits=3))
        p1 = histogram(heights, bins=0:2:46, title=gen, xaxis="height", ylims=(-0.1, length(trees) / 4), leg=false)
        p2 = histogram(xs, bins=0:0.2:5, title=gen, xaxis="donation rate (x)", ylims=(-0.1, length(trees) / 4), leg=false)
        p3 = histogram(energies, bins=-50:4:60, title=gen, xaxis="energy", ylims=(-0.1, length(trees) / 4), leg=false)
        plot(p1, p2, p3)

        # birth and death
        for idx in occ_idxs  # concurrent modification seems to be ok in this case
            tree = trees[idx]
            if tree.energy <= 0 || tree.age > age_cap  # death
                pop!(occ_idxs, idx)
                push!(free_idxs, idx)
            elseif length(free_idxs) > 0 && rand() < spawnrate * tree.energy  # birth
                new_height = min(45, max(1, tree.height + rand(-1:1)))
                new_x = max(0, tree.donation_rate + randn() / 10)
                new_tree = AltruisticTree(tree.id, new_height, new_x, tree.growth_cost, tree.growth_gain, 0, 0)
                new_idx = rand(free_idxs)
                pop!(free_idxs, new_idx)
                push!(occ_idxs, new_idx)
                trees[new_idx] = new_tree
            end
        end
    end
    anim, avg_energies, avg_heights, counts, avg_xs
end

function continuous_time_evolve_no_animate!(trees::Vector{AltruisticTree}, occ_idxs::Set{Int64}, free_idxs::Set{Int64}, generations::Int64, age_cap::Int64=100, spawnrate::Float64=0.01)
    for gen in 1:generations
        
        # update energies and ages
        for idx in occ_idxs
            tree = trees[idx]
            idx_left = get_left_neigh_idx(trees, occ_idxs, idx)
            neigh_left = trees[idx_left]
            idx_right = get_right_neigh_idx(trees, occ_idxs, idx)
            neigh_right = trees[idx_right]
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
                new_height = min(45, max(1, tree.height + rand(-1:1)))
                new_x = max(0, tree.donation_rate + randn() / 10)
                new_tree = AltruisticTree(tree.id, new_height, new_x, tree.growth_cost, tree.growth_gain, 0, 0)
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

# 1 normal trees
# 2 altruistic trees
# 3 normal vs altruistic trees
x, y, z = 1, 3, 4
num_trees = 200
xdim = 4 * num_trees
ydim = 50  # ⟹ 45 is max height
grid = zeros(Int32, ydim, xdim)
grid[1, :] .= -2 * num_trees  # ground
grid[ydim - 3:ydim - 1, 2:2 + num_trees ÷ 7] .= 2 * num_trees  # sunshine!
generations = 10_000
age_cap = 100
spawnrate = 0.1
name = string("gens=", generations, " trees=", num_trees, " x,y,z=", x, ",", y, ",", z,
              " sr=", spawnrate, " ac=", age_cap, " seed=", seed)
    
function both_trees_init()
    trees = Vector{Tree}()  # from l to r
    for i in 1:num_trees
        height = rand(10:20)
        if rand() < ratio
            push!(trees, NormalTree(i + 10, height, y, z, 0, 0))
        else
            push!(trees, AltruisticTree(-i - 10, height, x, y, z, 0, 0))
        end
    end
    results = both_animate!(trees, grid, similar(trees), generations, name)
    results
end

function continuous_time_init()
    trees = Vector{AltruisticTree}(undef, 2 * num_trees)
    perm = randperm(length(trees))
    idxs = Set(perm[1:num_trees])
    free_idxs = Set(perm[num_trees + 1:end])
    for i in idxs
        height = rand(10:20)
        age = rand(0:age_cap - 1)
        trees[i] = AltruisticTree(-i - 10, height, max(x + randn(), 0), y, z, 0, age)
    end
    results = continuous_time_evolve!(trees, idxs, free_idxs, generations, age_cap, spawnrate)
    results
end
    
function init_test_x()
    n_trials = 10
    n_quantiles = 5
    std = 0.1
    d = 1 / (n_quantiles + 1)
    quantiles = d:d:1 - d
    x_inits = 0:0.25:3
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
                height = rand(10:20)
                age = rand(0:age_cap - 1)
                trees[i] = AltruisticTree(-i - 10, height, max(x + randn() * std, 0), y, z, 0, age)
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