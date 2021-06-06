using Base:Integer, UInt
cd(@__DIR__)  # this sets the directory pointer to this file location
using Pkg  # Pkg.<command> is equivalent of `]` in command line julia
Pkg.activate("..")  # this specifies the environment being used and activates it

using Plots
# using TrackingHeaps
using Random
import Base.getindex, Base.size, Base.setindex!, Base.length, Base.+, Base.-, Base.*, Base.÷, Base.show, Base.print, Base.string, Base.println, Base.<, Base.>, Base.<=, Base.>=
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

# wraps numbers to stay between 1 <= y <= m
modp1(x::Integer, m::Integer) = x < 1 ? m + x % m : 1 + (x - 1) % m

"""
modular integer
where 1 <= val <= mod
"""
struct MInt <: Integer
    val::UInt8
    mod::UInt8
    MInt(val, mod) = 1 <= val <= mod ? new(val, mod) : new(modp1(val, mod), mod)
end

getindex(A::Matrix, i1::MInt, i2::MInt) = A[i1.val, i2.val]
function setindex!(A::Matrix, v::Number, i1::MInt, i2::MInt)
    A[i1.val, i2.val] = v
end
+(a::MInt, b::Integer) = MInt(a.val + b, a.mod)
+(b::Integer, a::MInt) = MInt(a.val + b, a.mod)
-(b::Integer, a::MInt) = MInt(b - a.val, a.mod)
-(a::MInt, b::Integer) = MInt(a.val - b, a.mod)
*(a::MInt, b::Integer) = MInt(a.val * b, a.mod)
*(b::Integer, a::MInt) = MInt(a.val * b, a.mod)
÷(b::Integer, a::MInt) = MInt(b ÷ a.val, a.mod)
÷(a::MInt, b::Integer) = MInt(a.val ÷ b, a.mod)
<(a::MInt, b::Integer) = a.val < b
<(b::Integer, a::MInt) = a.val > b
<(a::MInt, b::MInt) = a.val > b.val
>(a::MInt, b::Integer) = a.val > b
>(b::Integer, a::MInt) = a.val < b
>(a::MInt, b::MInt) = a.val > b.val
<=(a::MInt, b::Integer) = a.val <= b
<=(b::Integer, a::MInt) = a.val >= b
<=(a::MInt, b::MInt) = a.val <= b.val
>=(a::MInt, b::Integer) = a.val >= b
>=(b::Integer, a::MInt) = a.val <= b
>=(a::MInt, b::MInt) = a.val >= b.val

show(io::IO, a::MInt) = print(io, a.val, "%", a.mod)

"""
matrix indexer whose indexing wraps
"""
struct GridIndex
    I::Tuple{MInt, MInt}
    GridIndex(x::MInt, y::MInt) = new(Tuple((x, y)))
    GridIndex(tup::Tuple{MInt, MInt}) = new(tup)
end

getindex(A::Matrix, I::GridIndex) = A[I.I[1], I.I[2]]
function setindex!(A::Matrix, v::Number, I::GridIndex)
    A[I.I[1], I.I[2]] = v
end

struct Pattern
    coords::Set{GridIndex}
end

length(p::Pattern) = length(p.coords)


### frame to frame updates
function updategrid!(patterns::Set{Pattern}, cells::Matrix{Int8}, tmp::Matrix{Int8})
    xdim, ydim = size(cells)
    tmp .= (cells .> 0) .+ cells
    newcells = Set{GridIndex}()
    println()
    println("sum ", sum(length.(patterns)))
    for pattern in patterns
        pat_size = length(pattern.coords)
        if pat_size > 0
            newcoords = Set{GridIndex}()
            if pat_size % 4 < 2
                extremum = argmax(x -> x.I[pat_size % 2 + 1], pattern.coords)
            else
                extremum = argmin(x -> x.I[pat_size % 2 + 1], pattern.coords)
            end
            start = GridIndex(extremum.I[1], extremum.I[2] + 1)
            for coord in pattern.coords
                if cells[coord] == AGE_CAP
                    prev = extremum
                    current = start
                    step = 1
                    cutoff = 4 * pat_size  # max possible surface area for squares
                    while (step <= pat_size - 1 || current in newcells) && step <= cutoff
                        step += 1
                        prev, current = makestep(prev, current, pattern)
                    end
                    if step == cutoff + 1
                        current = coord  # if you can't find anywhere to grow, fill the space of the cell the just died
                    end
                    new = current
                    tmp[coord] = 0
                    push!(newcells, new)
                    tmp[new] = AGE_MIN
                    push!(newcoords, new)
                    if AGE_MIN <= cells[current] < AGE_CAP  # overwriting foreign cell
                        popped = false  # TODO remove
                        for p in patterns
                            if current in p.coords
                                pop!(p.coords, current)
                                if popped
                                    println("popped twice???")
                                end
                                popped = true
                            end
                        end
                        tmp[coord] = 1
                        push!(newcoords, coord)
                    end
                    start = new
                else
                    push!(newcoords, coord)
                end
            end
            updatecoords!(pattern, newcoords)
        else
            println("dead pattern")
        end
    end
    patterns, tmp, cells
end

function indexof(vec::Vector{T}, val::T) where T
    for i in 1:length(vec)
        if vec[i] == val
            return i
        end
    end
    -1
end

function makestep(prev::GridIndex, current::GridIndex, pattern::Pattern)
    neighs = GridIndex.([   (current.I[1], current.I[2] + 1)
                            (current.I[1] + 1, current.I[2] + 1)
                            (current.I[1] + 1, current.I[2])
                            (current.I[1] + 1, current.I[2] - 1)
                            (current.I[1], current.I[2] - 1)
                            (current.I[1] - 1, current.I[2] - 1)
                            (current.I[1] - 1, current.I[2])
                            (current.I[1] - 1, current.I[2] + 1)])
    done = false
    started = false
    # println(current)
    # println("neighs ", neighs, "\n")
    i = indexof(neighs, prev)
    prev = current
    while i <= 16  # 2 len neighs
        i += 1
        idx = i % 8 + 1  # increment idx cyclically
        if !started
            if !(neighs[idx] in pattern.coords)  # is the location part of my body
                started = true
            end
        else
            if neighs[idx] in pattern.coords 
                break  # we let current be the previous neighbor once we encounter a live cell
            end
        end
        current = neighs[idx]
    end
    prev, current
end

function updatecoords!(pattern::Pattern, coords::Set{GridIndex})
    for _ in 1:length(pattern.coords)
        pop!(pattern.coords)
    end
    for _ in 1:length(coords)
        push!(pattern.coords, pop!(coords))
    end
end

const AGE_MIN = 16
const AGE_CAP = 32
xdim, ydim = 40, 60
x = 1:xdim
y = 1:ydim

# initialize patterns
pattern_locs = [GridIndex(MInt(i, xdim), MInt(j, ydim)) for i in 0:20:xdim - 1 for j in 0:20:ydim - 1]
# pattern_locs = CartesianIndex.([(5, 10), (35, 10), (40, 25), (20, 30), (10, 55), (5, 50), (35, 40), (40, 60)])
# pattern_locs = CartesianIndex.([(5, 10)])
patterns = Set{Pattern}()
for ploc in pattern_locs
    xd = rand(1:20)
    yd = rand(1:20)
    n = rand(xd * yd ÷ 2 + 1:xd * yd)
    possible_coords = Set{GridIndex}([GridIndex(i + ploc.I[1], j + ploc.I[2]) for i in 1:xd for j in 1:yd])
    coords = Set{GridIndex}(rand(possible_coords, n))
    push!(patterns, Pattern(coords))
end

cells = zeros(Int8, xdim, ydim)
for pattern in patterns
    for c in pattern.coords
        cells[c] = rand(AGE_MIN:AGE_CAP)
    end
end

iterations = 1000
name = string("main ", xdim, "x", ydim, " seed=", seed, " iters=", iterations)
tmp = similar(cells)
bsize = 0

anim = @animate for i in 1:iterations
    # NOTE: x is vertical, y is horizontal, so swap the arguments for reasonable results
    heatmap(bsize + 1:ydim - bsize, bsize + 1:xdim - bsize, view(cells, bsize + 1:xdim - bsize, bsize + 1:ydim - bsize),
    colorbar=false, grid=false, ticks=false,
    showaxis=false, title=name, aspect_ratio=1, dpi=150)
    patterns, cells, tmp = updategrid!(patterns, cells, tmp)
end

gif(anim, string("gifs\\", name, ".gif"), fps=60)
