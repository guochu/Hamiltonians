abstract type AbstractInfiniteQuantumOperator{S} end

# storage should be implemented, which is a Dict type

TK.spacetype(x::AbstractInfiniteQuantumOperator) = spacetype(typeof(x))
TK.spacetype(::Type{<:AbstractInfiniteQuantumOperator{S}}) where {S} = S
Base.isempty(s::AbstractInfiniteQuantumOperator) = isempty(storage(s))
Base.length(x::AbstractInfiniteQuantumOperator) = length(physical_spaces(x))
unitcell_size(x::AbstractInfiniteQuantumOperator) = length(x)

isconstant(x::AbstractInfiniteQuantumOperator) = all(isconstant, storage(x))

# interactionrange(x::AbstractInfiniteQuantumOperator) = maximum([interactionrange(k) for k in storage(x)])

Base.:*(y::AllowedCoefficient, x::AbstractInfiniteQuantumOperator) = x * y
Base.:/(x::AbstractInfiniteQuantumOperator, y::AllowedCoefficient) = x * (1 / Coefficient(y))
Base.:+(x::AbstractInfiniteQuantumOperator) = x
Base.:-(x::AbstractInfiniteQuantumOperator) = (-1) * x 

function Base.keys(x::AbstractInfiniteQuantumOperator)
	r = Set{Vector{Int}}()
	for item in qterms(x)
		push!(r, positions(item))
	end
	return r
end