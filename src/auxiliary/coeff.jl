# I support time-dependent quantum operator by using this class
abstract type AbstractCoefficient end

"""
	struct Coefficient{T} <: AbstractCoefficient

T is either a scalar or an unary function
"""
struct Coefficient{T} <: AbstractCoefficient
	value::T

function Coefficient{T}(f) where {T<:Function}
	m = f(1.)
	isa(m, Number) || throw(ArgumentError("f should be either a scalar or a unary function"))
	new{T}(f)
end
Coefficient{T}(f) where {T<:Number} = new{T}(f)
end

Coefficient(f::Union{Number, Function}) = Coefficient{typeof(f)}(f)
Coefficient(x::Coefficient) = Coefficient(value(x))
const AllowedCoefficient = Union{Number, Function, Coefficient}
const ScalarCoefficient = Coefficient{<:Number}
const FunctionCoefficient = Coefficient{<:Function}

# value(x::Number) = x
value(x::Coefficient) = x.value
TK.scalar(x::ScalarCoefficient) = value(x)
TK.scalar(x::FunctionCoefficient) = throw(ArgumentError("cannot convert a function into a scalar"))
Base.copy(x::Coefficient) = Coefficient(value(x))
DMRG.coeff(x::AllowedCoefficient) = Coefficient(x)

Base.:+(x::ScalarCoefficient, y::ScalarCoefficient) = Coefficient(value(x) + value(y))
Base.:+(x::FunctionCoefficient, y::ScalarCoefficient) = x + value(y)
Base.:+(y::ScalarCoefficient, x::FunctionCoefficient) = x + y
Base.:+(x::FunctionCoefficient, y::FunctionCoefficient) = Coefficient(z->value(x)(z) + value(y)(z))

Base.:-(x::ScalarCoefficient, y::ScalarCoefficient) = Coefficient(value(x) - value(y))
Base.:-(x::FunctionCoefficient, y::ScalarCoefficient) = x - value(y)
Base.:-(y::ScalarCoefficient, x::FunctionCoefficient) = -x + y
Base.:-(x::FunctionCoefficient, y::FunctionCoefficient) = Coefficient(z->value(x)(z) - value(y)(z))

Base.:*(x::ScalarCoefficient, y::ScalarCoefficient) = Coefficient(value(x) * value(y))
Base.:*(x::FunctionCoefficient, y::ScalarCoefficient) = x * value(y)
Base.:*(y::ScalarCoefficient, x::FunctionCoefficient) = x * y
Base.:*(x::FunctionCoefficient, y::FunctionCoefficient) = Coefficient(z->value(x)(z) * value(y)(z))

Base.:/(x::ScalarCoefficient, y::ScalarCoefficient) = Coefficient(value(x) / value(y))
Base.:/(x::FunctionCoefficient, y::ScalarCoefficient) = x / value(y)
Base.:/(y::ScalarCoefficient, x::FunctionCoefficient) = Coefficient(z -> y(value(x)(z)))
Base.:/(x::FunctionCoefficient, y::FunctionCoefficient) = Coefficient(z->value(x)(z) / value(y)(z))


Base.:+(x::ScalarCoefficient, y::Number) = Coefficient(value(x) + y)
Base.:+(x::FunctionCoefficient, y::Number) = Coefficient(z->value(x)(z+y))
Base.:+(y::Number, x::Coefficient) = x + y

Base.:-(x::ScalarCoefficient, y::Number) = Coefficient(value(x) - y)
Base.:-(x::FunctionCoefficient, y::Number) = Coefficient(z->value(x)(z-y))
Base.:-(y::Number, x::Coefficient) = -x + y

Base.:*(x::ScalarCoefficient, y::Number) = Coefficient(value(x) * y)
Base.:*(x::FunctionCoefficient, y::Number) = Coefficient(z->value(x)(z)*y)
Base.:*(y::Number, x::Coefficient) = x * y

Base.:/(x::ScalarCoefficient, y::Number) = Coefficient(value(x) / y)
Base.:/(x::FunctionCoefficient, y::Number) = Coefficient(z->value(x)(z)/y)
Base.:/(x::Number, y::ScalarCoefficient) = Coefficient(x / value(y))
Base.:/(x::Number, y::FunctionCoefficient) = Coefficient(z->x/(value(y)(z)))

Base.:-(x::ScalarCoefficient) = Coefficient(-value(x))
Base.:-(x::FunctionCoefficient) = Coefficient(z->-value(x)(z))

# only implemented those I need, not aimed for general unary operations
Base.conj(x::ScalarCoefficient) = Coefficient(conj(value(x)))
Base.conj(x::FunctionCoefficient) = Coefficient(z->conj(value(x)(z)))
Base.sqrt(x::ScalarCoefficient) = Coefficient(sqrt(value(x)))
Base.sqrt(x::FunctionCoefficient) = Coefficient(z->sqrt(value(x)(z)))

(op::ScalarCoefficient)(t::Number) = value(op)
(op::FunctionCoefficient)(t::Number) = value(op)(t)

"""
	isconstant(x)
	check if the object x is a pure constant or contrans function 
"""
isconstant(x::ScalarCoefficient) = true
isconstant(x::FunctionCoefficient) = false


Base.:(==)(x::Coefficient, y::Coefficient) = false
Base.:(==)(x::ScalarCoefficient, y::ScalarCoefficient) = value(x) == value(y)
Base.:(==)(x::FunctionCoefficient, y::FunctionCoefficient) = value(x) === value(y)

Base.isapprox(x::Coefficient, y::Coefficient; kwargs...) = false
Base.isapprox(x::ScalarCoefficient, y::ScalarCoefficient; kwargs...) = isapprox(x, y; kwargs...)
Base.isapprox(x::FunctionCoefficient, y::FunctionCoefficient; kwargs...) = x == y

"""
	Base.eltype(x)

Return the scale type of x, which must be a subtype of Number.
"""
TK.scalartype(::Type{Coefficient{T}}) where {T<:Number} = T
TK.scalartype(x::ScalarCoefficient) = scalartype(typeof(x))
# Base.eltype(x::FunctionCoefficient) = begin
#     m = value(x)(1.)
#     return typeof(m)
# end

"""
	iszero(x) 
	check if x is zero.
"""
Base.iszero(x::ScalarCoefficient) = iszero(value(x))
Base.iszero(x::FunctionCoefficient) = false

Base.iszero(x::AbstractTensorMap) = isapprox(x, zero(x))
