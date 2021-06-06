using Base:Integer
cd(@__DIR__)  # this sets the directory pointer to this file location
using Pkg  # Pkg.<command> is equivalent of `]` in command line julia
Pkg.activate("..")  # this specifies the environment being used and activates it

using Plots
# using TrackingHeaps
using Random
import Base.getindex, Base.size, Base.setindex!, Base.length
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
length(p::Pattern) = length(p.coords)

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
    newcells = Set{CartesianIndex}()
    println()
    println("sum ", sum(length.(patterns)))
    for pattern in patterns
        if length(pattern) > 0
            newcoords = Set{CartesianIndex}()
            rightmost = argmax(x -> x.I[2], pattern.coords)
            start = CartesianIndex(rightmost.I[1], rightmost.I[2] + 1)
            for coord in pattern.coords
                if cells[coord] == AGE_CAP
                    current = start
                    step = 1
                    # the second condition actually does nothing as is because the edge-finding algo does not know friend from foe
                    while (step <= length(pattern.coords) - 1 || cells[current] > 0 || current in newcells) && step < 4 * length(pattern.coords)  # TODO stepping needs to be a bit smarter to handle weird shapes
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
                                if cells[neighs[idx]] < 1  # is the cell empty
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
                    push!(newcells, new)
                    tmp[new] = AGE_MIN
                    push!(newcoords, new)
                    start = new
                else
                    push!(newcoords, coord)
                end
            end
            updatecoords!(pattern, newcoords)
        else
            prinln("DEAD PATTERN")
        end
    end
    patterns, tmp, cells
end

function updatecoords!(pattern::Pattern, coords::Set{CartesianIndex})
    for _ in 1:length(pattern.coords)
        pop!(pattern.coords)
    end
    for _ in 1:length(coords)
        push!(pattern.coords, pop!(coords))
    end
end

const AGE_MIN = 16
const AGE_CAP = 32
xdim, ydim = 100, 140
x = 1:xdim
y = 1:ydim

# initialize patterns
pattern_locs = [CartesianIndex(i, j) for i in 0:20:xdim - 1 for j in 0:20:ydim - 1]
# pattern_locs = CartesianIndex.([(5, 10), (35, 10), (40, 25), (20, 30), (10, 55), (5, 50), (35, 40), (40, 60)])
# pattern_locs = CartesianIndex.([(5, 10)])
patterns = Set{Pattern}()
for ploc in pattern_locs
    xd = rand(1:20)
    yd = rand(1:20)
    n = rand(xd * yd ÷ 2 + 1:xd * yd)
    possible_coords = Set{CartesianIndex}([CartesianIndex(i + ploc.I[1], j + ploc.I[2]) for i in 1:xd for j in 1:yd])
    coords = Set{CartesianIndex}(rand(possible_coords, n))
    push!(patterns, Pattern(coords))
end

cells = zeros(Int8, xdim, ydim)
for pattern in patterns
    for c in pattern.coords
        cells[c] = rand(AGE_MIN:AGE_CAP)
    end
end

cells = Grid(cells)

iterations = 1000
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

gif(anim, string("gifs\\", name, ".gif"), fps=15)
