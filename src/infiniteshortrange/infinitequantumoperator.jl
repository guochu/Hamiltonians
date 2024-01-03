"""
	struct InfiniteQuantumOperator{S <: ElementarySpace}

Stores all the terms with its first position in the first unitcell
"""
struct InfiniteQuantumOperator{S <: ElementarySpace} <: AbstractInfiniteQuantumOperator{S}
	physpaces::PeriodicArray{S, 1}
	data::Vector{Any}

function InfiniteQuantumOperator(physpaces::PeriodicArray{S, 1}, ms::Vector{Any}) where {S<:ElementarySpace}
	@assert all(x->isa(x, QTerm), ms)
	return new{S}(physpaces, ms)
end

end
InfiniteQuantumOperator(physpaces::AbstractVector{S}, ms::Vector{Any}) where {S<:ElementarySpace} = InfiniteQuantumOperator(PeriodicArray(physpaces), ms)
InfiniteQuantumOperator(physpaces::AbstractVector{S}) where {S<:ElementarySpace} = InfiniteQuantumOperator(physpaces, Vector{Any}())
InfiniteQuantumOperator(physpaces::AbstractVector{<:ElementarySpace}, ms::Vector{<:QTerm}) = InfiniteQuantumOperator(physpaces, convert(Vector{Any}, ms))


storage(x::InfiniteQuantumOperator) = x.data
DMRG.physical_spaces(x::InfiniteQuantumOperator) = x.physpaces

function Base.push!(x::InfiniteQuantumOperator, m::QTerm) 
	(spacetype(x) == spacetype(m)) || throw(SpaceMismatch())
	iszero(m) && return x
	(0 < positions(m)[1] <= unitcell_size(x)) || throw(DimensionMismatch("first position should be within cell"))
	isstrict(m) || throw(SpaceMismatch("only strict QTerm is allowed"))
	push!(x.data, m)
	return x
end  
Base.push!(x::InfiniteQuantumOperator, m::AdjointQTerm) = push!(x, convert(QTerm, m))

TK.scalartype(x::InfiniteQuantumOperator) = compute_scalartype(storage(x))
Base.copy(x::InfiniteQuantumOperator) = InfiniteQuantumOperator(copy(physical_spaces(x)), copy(storage(x)))

(x::InfiniteQuantumOperator)(t::Number) = InfiniteQuantumOperator(physical_spaces(x), [item(t) for item in storage(x)])
Base.:*(x::InfiniteQuantumOperator, y::AllowedCoefficient) = InfiniteQuantumOperator(physical_spaces(x), [item * y for item in storage(x)])
function Base.:+(x::InfiniteQuantumOperator, y::InfiniteQuantumOperator)
	(physical_spaces(x) == physical_spaces(y)) || throw(SpaceMismatch("only infinite operator with same unit cell can be added"))
	return InfiniteQuantumOperator(physical_spaces(x), vcat(storage(x), storage(y)) )
end 
Base.:-(x::InfiniteQuantumOperator, y::InfiniteQuantumOperator) = x + (-y)


qterms(x::InfiniteQuantumOperator) = storage(x)

DMRG.isstrict(x::InfiniteQuantumOperator) = all(isstrict, storage(x))

function changeunitcell(m::InfiniteQuantumOperator; unitcellsize::Int)
	L = unitcell_size(m)
	(unitcellsize == L) && return m
	(unitcellsize % L == 0) || error("new unitcellsize $unitcellsize mismatch with curren unitcell size $L")
	nperoid = div(unitcellsize, L)
	physpaces = repeat(physical_spaces(m), nperoid)
	r = InfiniteQuantumOperator(physpaces)
	for i in 1:nperoid
		start = (i - 1) * L
		for t in qterms(m)
			push!(r, shift(t, start))
		end
	end
	return r
end

function absorb_one_bodies(x::InfiniteQuantumOperator) 
	L = length(x)
	physpaces = physical_spaces(x)
	y = InfiniteQuantumOperator(physpaces)
	if L == 1
		iden = id(oneunit(physpaces[1]) ⊗ physpaces[1])
		for item in storage(x)
			if length(positions(item)) == 1
				i = positions(item)[1]
				@assert i == 1
				t = QTerm([i, i+1], [op(item)[1], iden], coeff=0.5*coeff(item))
				push!(y, t)
				t = QTerm([i, i+1], [iden, op(item)[1]], coeff=0.5*coeff(item))
				push!(y, t)
			else
				push!(y, item)
			end
		end
	else
		for item in storage(x)
			if length(positions(item)) == 1
				i = positions(item)[1]
				if i < L
					iden = id(oneunit(physpaces[i+1]) ⊗ physpaces[i+1])
					t = QTerm([i, i+1], [op(item)[1], iden], coeff=coeff(item))
				else
					iden = id(oneunit(physpaces[i-1]) ⊗ physpaces[i-1])
					t = QTerm([i-1, i], [iden, op(item)[1]], coeff=coeff(item))
				end
				push!(y, t)
			else
				push!(y, item)
			end
		end
	end
	return y
end

"""
	oneperiod(m::InfiniteQuantumOperator)
Get one period of Hamiltonian as a finite QuantumOperator
this function is used for TEBD
"""
function oneperiod(m::InfiniteQuantumOperator)
	L = maximum(x->maximum(positions(x)), qterms(m))
	physpaces = physical_spaces(m)[1:L]
	return QuantumOperator(physpaces, storage(m))
end
