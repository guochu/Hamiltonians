abstract type AbstractQuantumOperator{S} end

# storage should be implemented, which is a Dict type

TK.spacetype(x::AbstractQuantumOperator) = spacetype(typeof(x))
TK.spacetype(::Type{<:AbstractQuantumOperator{S}}) where {S} = S
Base.isempty(s::AbstractQuantumOperator) = isempty(storage(s))
Base.length(x::AbstractQuantumOperator) = length(physical_spaces(x))

isconstant(x::AbstractQuantumOperator) = all(isconstant, storage(x))

interactionrange(x::AbstractQuantumOperator) = maximum([interactionrange(k) for k in storage(x)])

Base.:*(y::AllowedCoefficient, x::AbstractQuantumOperator) = x * y
Base.:/(x::AbstractQuantumOperator, y::AllowedCoefficient) = x * (1 / Coefficient(y))
Base.:+(x::AbstractQuantumOperator) = x
Base.:-(x::AbstractQuantumOperator) = (-1) * x 

function Base.keys(x::AbstractQuantumOperator)
	r = Set{Vector{Int}}()
	for item in qterms(x)
		push!(r, positions(item))
	end
	return r
end