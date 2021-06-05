using Base:Integer
cd(@__DIR__)  # this sets the directory pointer to this file location
using Pkg  # Pkg.<command> is equivalent of `]` in command line julia
Pkg.activate("..")  # this specifies the environment being used and activates it

using Plots
# using TrackingHeaps
using Random
import Base.getindex, Base.size, Base.setindex!
using BenchmarkTools

"""
fundamentally in evolution every pattern (organism) needs a source of information (matter and energy) to succeed

self, discrete and deterministic, reasonable way of gaining information to allow reproduction
    - as soon as they have enough bits to replicate, they replicate (into a new grid)

    - or an incomplete copy of them idly travels around the next grid and is built up as the pattern finds more bits (in the form of other incomplete patterns)
        - motion: one cell is moved to a edge location which rotates counterclockwise each frame (so motion is quasi-circular)
        - or, the most recently moved cell is the starting point for counting (+x or cckmost margin), count n cckwise, where n is size of self

fractal evolution (replication is at smaller scales)? Then info is exponential
"""

seed = 633
Random.seed!(seed)

const AGE_CAP = 16
xdim, ydim = 60, 80
x = 1:xdim
y = 1:ydim

struct Pattern
    coords::Set{CartesianIndex}  # push!(coords, CartesianIndex(1,1))
end

"""
array whose indexing wraps
must be used with the below function definitions
"""
struct Grid
    M::Matrix
end

# wraps numbers to stay between 1 <= y <= m
modp1(x::Integer, m::Integer) = x < 1 ? m + x % m : 1 + (x - 1) % m

size(G::Grid) = size(G.M)

function getindex(G::Grid, c::CartesianIndex)
    mx, my = size(G.M)
    G.M[modp1(c.I[1], mx), modp1(c.I[2], my)]
end

function setindex!(G::Grid, val::Number, c::CartesianIndex)
    mx, my = size(G.M)
    G.M[modp1(c.I[1], mx), modp1(c.I[2], my)] = val
end

function updategrid!(patterns::Set{Pattern}, cells::Grid, tmp::Grid)
    xdim, ydim = size(cells)
    tmp.M .= (cells.M .> 0) .+ cells.M
    for pattern in patterns
        if length(pattern.coords) == 0
            println("empty pattern")
        end
        newcoords = Set{CartesianIndex}()
        rightmost = argmax(x -> x.I[2], pattern.coords)
        start = CartesianIndex(rightmost.I[1], rightmost.I[2] + 1)
        for coord in pattern.coords
            if cells[coord] == AGE_CAP
                current = start
                step = 1
                # the second condition actually does nothing as is because the edge-finding algo does not know friend from foe
                while step <= length(pattern.coords) - 1 || cells[current] > 0  # TODO stepping needs to be a bit smarter to handle weird shapes
                    step += 1
                    neighs = [  CartesianIndex(current.I[1], current.I[2] + 1)
                                CartesianIndex(current.I[1] + 1, current.I[2] + 1)
                                CartesianIndex(current.I[1] + 1, current.I[2])
                                CartesianIndex(current.I[1] + 1, current.I[2] - 1)
                                CartesianIndex(current.I[1], current.I[2] - 1)
                                CartesianIndex(current.I[1] - 1, current.I[2] - 1)
                                CartesianIndex(current.I[1] - 1, current.I[2])
                                CartesianIndex(current.I[1] - 1, current.I[2] + 1)]
                    done = false
                    started = false
                    # println(current)
                    # println("neighs ", neighs, "\n")
                    i = 1
                    while i <= 16
                        i += 1
                        idx = i % 8 + 1  # increment idx cyclically
                        if !started
                            if cells[neighs[idx]] == 0  # is the cell empty
                                started = true
                            end
                        else
                            if cells[neighs[idx]] > 0  
                                break  # we update current to the previous neighbor once we encounter a live cell
                            end
                        end
                        current = neighs[idx]
                    end
                end
                new = current
                tmp[coord] = 0
                tmp[new] = 1  # birth new cell TODO check if this is vacant
                push!(newcoords, new)
                start = new
            else
                push!(newcoords, coord)
            end
        end
        updatecoords!(pattern, newcoords)
    end
    patterns, tmp, cells
end

function updatecoords!(pattern::Pattern, coords::Set{CartesianIndex})
    for i in 1:length(pattern.coords)
        pop!(pattern.coords)
    end
    for i in 1:length(coords)
        push!(pattern.coords, pop!(coords))
    end
end

# pattern_region1 = [ 1 5 6 7 2 3 0;
#                     9 8 2 6 9 2 11; 
#                     0 4 5 3 4 2 9; 
#                     0 0 0 1 2 1 1]
# pattern_region1 = [ 0 0 0 0 0 0 0;
#                     1 3 5 7 9 11 13; 
#                     0 0 1 0 0 0 0; 
#                     0 0 0 0 0 0 0]
# pattern_region2 = [15 13 10 9 9 5 0;
#                    14 11 0  8 7 2 1; 
#                     0 8  5  3 4 1 1; 
#                     0 0  0  0 1 2 0]
# cells[(1:4) .+ 20, (1:7) .+ 25] = pattern_region1
# cells[(1:4) .+ 35, (1:7) .+ 25] = pattern_region2
p1 = Pattern(Set(CartesianIndex[CartesianIndex(37, 26), CartesianIndex(38, 27), CartesianIndex(36, 31), CartesianIndex(37, 30), CartesianIndex(36, 29), CartesianIndex(38, 30), CartesianIndex(38, 28), CartesianIndex(39, 30), CartesianIndex(37, 32), CartesianIndex(38, 32), CartesianIndex(37, 31), CartesianIndex(38, 31), CartesianIndex(36, 
27), CartesianIndex(39, 31), CartesianIndex(36, 26), CartesianIndex(37, 29), CartesianIndex(38, 29), CartesianIndex(36, 30), CartesianIndex(36, 28), CartesianIndex(37, 27)]))
p2 = Pattern(Set(CartesianIndex[CartesianIndex(23, 28), CartesianIndex(22, 32), CartesianIndex(22, 29), CartesianIndex(22, 27), CartesianIndex(22, 26), CartesianIndex(22, 31), CartesianIndex(22, 30), CartesianIndex(22, 28)]))
patterns = Set{Pattern}([p1, p2])

cells = zeros(UInt8, xdim, ydim)
for pattern in patterns
    for c in pattern.coords
        cells[c] = rand(1:AGE_CAP)
    end
end

cells = Grid(cells)

iterations = 500
name = string("main ", xdim, "x", ydim, " seed=", seed, " iters=", iterations)
tmp = Grid(similar(cells.M))
bsize = 0

anim = @animate for i in 1:iterations
    # NOTE: x is vertical, y is horizontal, so swap the arguments for reasonable results
    heatmap(bsize + 1:ydim - bsize, bsize + 1:xdim - bsize, view(cells.M, bsize + 1:xdim - bsize, bsize + 1:ydim - bsize),
    colorbar=false, grid=false, ticks=false,
    showaxis=false, title=name, aspect_ratio=1, dpi=150)
    patterns, cells, tmp = updategrid!(patterns, cells, tmp)
end

gif(anim, string("gifs\\", name, ".gif"), fps=8)
