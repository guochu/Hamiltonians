abstract type AbstractQuantumTerm end


positions(x::AbstractQuantumTerm) = x.positions
op(x::AbstractQuantumTerm) = x.op
DMRG.coeff(x::AbstractQuantumTerm) = x.coeff

DMRG.space_l(x::AbstractQuantumTerm) = space(op(x)[1], 1)
DMRG.space_r(x::AbstractQuantumTerm) = space(op(x)[end], 3)

TK.spacetype(x::AbstractQuantumTerm) = spacetype(typeof(x))


Base.isempty(x::AbstractQuantumTerm) = isempty(op(x))

Base.:*(m::AllowedCoefficient, s::AbstractQuantumTerm) = s * m
Base.:/(s::AbstractQuantumTerm, m::AllowedCoefficient) = s * (1 / Coefficient(m))
Base.:+(s::AbstractQuantumTerm) = s
Base.:-(s::AbstractQuantumTerm) = (-1) * s 



nterms(s::AbstractQuantumTerm) = length(op(s))
isconstant(s::AbstractQuantumTerm) = isconstant(coeff(s))
TK.scalartype(x::AbstractQuantumTerm) = scalartype(typeof(x))

function _interaction_range(x::Union{Vector{Int}, Tuple})::Int
	(length(x) == 0) && return 0
	(length(x)==1) && return 1
	return x[end] - x[1] + 1
end

interactionrange(x::AbstractQuantumTerm) = _interaction_range(positions(x))


function Base.iszero(x::AbstractQuantumTerm) 
	iszero(coeff(x)) && return true
	isempty(x) && return true
	for item in op(x)
	    iszero(item) && return true
	end
	return false
end

function DMRG.isstrict(t::AbstractQuantumTerm)
	isempty(t) && return true
	iszero(t) && return true
	return isoneunit(space_r(t)) 
end

DMRG.bond_dimension(h::AbstractQuantumTerm, bond::Int) = begin
	((bond >= 1) && (bond <= nterms(h))) || throw(BoundsError())
	dim(space(op(h)[bond], 3))
end 
DMRG.bond_dimensions(h::AbstractQuantumTerm) = [bond_dimension(h, i) for i in 1:nterms(h)]
DMRG.bond_dimension(h::AbstractQuantumTerm) = maximum(bond_dimensions(h))

DMRG.ophysical_spaces(psi::AbstractQuantumTerm) = [space(item, 2) for item in op(psi)]
DMRG.iphysical_spaces(psi::AbstractQuantumTerm) = [space(item, 4) for item in op(psi)]

function DMRG.physical_spaces(psi::AbstractQuantumTerm)
	xs = ophysical_spaces(psi)
	(xs == adjoint.(iphysical_spaces(psi))) || error("i and o physical dimension mismatch.")
	return xs
end

