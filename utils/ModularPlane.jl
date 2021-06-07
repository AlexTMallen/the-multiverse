module ModularPlane

using ExportAll
import Base.getindex, Base.size, Base.setindex!, Base.length, Base.+, Base.-, Base.*, Base.÷, Base.show, Base.print, Base.string, Base.println, Base.<, Base.>, Base.<=, Base.>=

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
@exportAll()

end