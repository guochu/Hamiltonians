struct QuantumOperator{S <: ElementarySpace} <: AbstractQuantumOperator{S}
	physpaces::Vector{S}
	data::Vector{Any}

function QuantumOperator(physpaces::Vector{S}, ms::Vector{Any}) where {S<:ElementarySpace}
	@assert all(x->isa(x, QTerm), ms)
	return new{S}(physpaces, ms)
end

end

storage(x::QuantumOperator) = x.data
DMRG.physical_spaces(x::QuantumOperator) = x.physpaces

# I assert the physical spaces to be set here
QuantumOperator(physpaces::Vector{S}) where {S<:ElementarySpace} = QuantumOperator(physpaces, Vector{Any}())
QuantumOperator(physpaces::Vector{<:ElementarySpace}, ms::Vector{<:QTerm}) = QuantumOperator(physpaces, convert(Vector{Any}, ms))

function QuantumOperator(ms::Vector{Any})
	@assert !isempty(ms)
	@assert all(x->isa(x, AbstractQuantumTerm), ms)
	ms = [isa(item, QTerm) ? item : convert(QTerm, item) for item in ms]
	# compute spaces
	S = spacetype(ms[1])
	physpaces = Vector{Union{Missing, S}}()
	u = oneunit(S)
	for item in ms
		pos = positions(item)
		L = length(physpaces)
		(space_r(item)' == u) || throw(SpaceMismatch("only strict QTerm is allowed"))
		if pos[end] > L
			resize!(physpaces, pos[end])
			for i in L+1:pos[end]
				physpaces[i] = missing
			end
		end
		for (posj, mj) in zip(pos, op(item))
			if !ismissing(physpaces[posj])
				(space(mj, 2) == physpaces[posj]) || throw(SpaceMismatch("physical space mismatch"))
			else
				physpaces[posj] = space(mj, 2)
			end
		end
	end
	any(ismissing, physpaces) && error("physical space is missing")
	physpaces = convert(Vector{S}, physpaces)
	return QuantumOperator(physpaces, ms)
end


"""
	push!(x::QuantumOperator{M}, m::QTerm{M}) where {M <: MPOTensor}
	adding a new term into the quantum operator
"""
function Base.push!(x::QuantumOperator, m::QTerm) 
	(spacetype(x) == spacetype(m)) || throw(SpaceMismatch())
	(positions(m)[end] <= length(x)) || throw(DimensionMismatch("QTerm out of Hamiltonian range"))
	isstrict(m) || throw(SpaceMismatch("only strict QTerm is allowed"))
	iszero(m) && return x
	push!(x.data, m)
	return x
end  
# m has to be strict
Base.push!(x::QuantumOperator, m::AdjointQTerm) = push!(x, convert(QTerm, m))


TK.scalartype(x::QuantumOperator) = compute_scalartype(storage(x))
Base.copy(x::QuantumOperator) = QuantumOperator(copy(physical_spaces(x)), copy(storage(x)))

(x::QuantumOperator)(t::Number) = QuantumOperator(physical_spaces(x), [item(t) for item in storage(x)])
Base.:*(x::QuantumOperator, y::AllowedCoefficient) = QuantumOperator(physical_spaces(x), [item * y for item in storage(x)])
Base.:+(x::QuantumOperator, y::QuantumOperator) = QuantumOperator(_merge_spaces(physical_spaces(x), physical_spaces(y)), vcat(storage(x), storage(y)))
Base.:-(x::QuantumOperator, y::QuantumOperator) = x + (-y)

function _merge_spaces(x, y)
	# assume length(x) >= length(y)
	if length(x) < length(y)
		return _merge_spaces(y, x)
	end
	L = length(x)
	for i in 1:length(y)
		(x[i] == y[i]) || throw(SpaceMismatch())
	end
	return copy(x)
end

DMRG.isstrict(x::QuantumOperator) = all(isstrict, storage(x))

qterms(x::QuantumOperator) = storage(x)
qterms(x::QuantumOperator, k::Vector{Int}) = [item for item in storage(x) if positions(item) == k]

function todict(x::QuantumOperator{S}) where {S}
	T = scalartype(x)
	M = mpotensortype(S, T)
	data = Dict{Tuple{Int, Vararg{Int, N} where N}, Vector{Tuple{Vector{M}, AbstractCoefficient}}}()
	for item in storage(x)
		pos = Tuple(positions(item))
		v = get!(data, pos, Vector{Tuple{Vector{M}, AbstractCoefficient}}())
		push!(v, (op(item), coeff(item)))
	end
	return data
end
function absorb_one_bodies(x::QuantumOperator) 
	L = length(x)
	physpaces = physical_spaces(x)
	y = QuantumOperator(physpaces)
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
	return y
end


# function _join_ops(m::QTerm)
# 	isconstant(m) || error(ArgumentError("functional term not allowed"))
# 	iszero(m) && return nothing
# 	isstrict(m) || error(ArgumentError("only strict QTerm allowed"))
# 	return _join(op(m)) * value(coeff(m))
# end
get_trivial_leg(m::AbstractTensorMap) = TensorMap(ones,scalartype(m),oneunit(space(m,1)), one(space(m,1)))