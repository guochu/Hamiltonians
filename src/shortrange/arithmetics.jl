function apply!(t::QTerm, mps::AbstractMPS)
	_start, _end = positions(t)[1], positions(t)[end]
	S = spacetype(mps)
	T = scalartype(mps)
	physpaces = [space(mps[i], 2) for i in _start:_end]
	mpo = prodmpo(T, physpaces, positions(t) .- (_start-1), op(t)) * scalar(coeff(t))
	mpo_scale = coeff(mpo)

	M = tensormaptype(S, 2, 3, T)
	r = Vector{M}(undef, _end - _start + 1)
	for (i, pos) in enumerate(_start:_end)
		r[i] = @tensor tmp[-1 -2; -3 -4 -5] := mpo_scale * mpo[i][-1, -3, -4, 1] * mps[pos][-2, 1, -5]
	end
	fusion_ts = [isomorphism(space(item, 4)' ⊗ space(item, 5)', fuse(space(item, 4)', space(item, 5)')) for item in r]
	left = isomorphism(fuse(space(r[1], 1), space(r[1], 2)), space(r[1], 1) ⊗ space(r[1], 2))
	mps[_start] = @tensor tmp[1,4;7] := left[1,2,3] * r[1][2,3,4,5,6] * fusion_ts[1][5,6,7]
	for (i, pos) in enumerate(_start+1:_end)
		mps[pos] = @tensor tmp[3,4;7] := conj(fusion_ts[i][1,2,3]) * r[i+1][1,2,4,5,6] * fusion_ts[i+1][5,6,7]
	end
	return mps
end
Base.:*(t::QTerm, mps::AbstractMPS) = apply!(t, copy(mps))


function mult(x::QTerm, y::AdjointQTerm)
	pos, opx, opy = _coerce_qterms(x, y.parent)
	v = DMRG._mult_n_a(opx, opy)
	return QTerm(pos, v, coeff=coeff(x)*coeff(y))
end
function mult(x::AdjointQTerm, y::QTerm)
	pos, opx, opy = _coerce_qterms(x.parent, y)
	v = DMRG._mult_a_n(opx, opy)
	return QTerm(pos, v, coeff=coeff(x)*coeff(y))
end
function mult(x::QTerm, y::QTerm)
	pos, opx, opy = _coerce_qterms(x, y)
	v = DMRG._mult_n_n(opx, opy)
	return QTerm(pos, v, coeff=coeff(x)*coeff(y))
end
mult(x::AdjointQTerm, y::AdjointQTerm)  = adjoint(y.parent * x.parent)

"""
	Base.:*(x::QTerm, y::AdjointQTerm)
multiplication between two QTerms
"""
Base.:*(x::AbstractQuantumTerm, y::AbstractQuantumTerm) = mult(x, y)

function Base.:+(x::QTerm, y::QTerm)
	(isconstant(x) && isconstant(y)) || throw(ArgumentError("only constant QTerms can be added"))
	pos, hA, hB = _coerce_qterms(x, y)
	c1, c2 = scalar(coeff(x)), scalar(coeff(y))
	hA[1] *= c1
	hB[1] *= c2
	hC = MPO(hA) + MPO(hB)
	return QTerm(pos, hC.data, coeff=1)
end
Base.:-(x::QTerm, y::QTerm) = x + (-y)

function Base.:*(x::QuantumOperator, y::QuantumOperator)
	a = qterms(x)
	b = qterms(y)
	return QuantumOperator(_merge_spaces(physical_spaces(x), physical_spaces(y)), [aj * bj for aj in a for bj in b])
end

"""
	simplify(m::QTerm; atol::Real=1.0e-12) 

Simplifying QTerm by removing the identities
"""
simplify(m::QTerm; atol::Real=1.0e-12) = remove_identities(m; atol=atol)
function remove_identities(m::QTerm; kwargs...)
	_pos = positions(m)
	_ops = op(m)

	new_pos = Int[]
	new_ops = []
	scale = 1.
	for (a, b) in zip(_pos, _ops)
		_is_id, scal = DMRG.isid(b; kwargs...)
		if _is_id
			scale *= scal
		else
			push!(new_pos, a)
			push!(new_ops, b)
		end
	end
	if isempty(new_pos)
		@warn "QTerm is identity"
		return scale * coeff(m)
	else
		return QTerm(new_pos, [new_ops...], coeff=scale * coeff(m))
	end
end

simplify(m::AdjointQTerm; kwargs...) = AdjointQTerm(simplify(m.parent; kwargs...))

function _coerce_qterms(x::QTerm, y::QTerm)
	T = promote_type(scalartype(x), scalartype(y))
	S = spacetype(x)
	M = mpotensortype(S, T)
    opx = convert(Vector{M}, copy(op(x))) 
    opy = convert(Vector{M}, copy(op(y))) 
    pos = positions(x)
    if !(positions(x) == positions(y))
    	x_left = oneunit(spacetype(x))
    	y_left = oneunit(spacetype(y))

    	# new_pos = sort([Set(vcat(positions(x), positions(y)))...])
    	new_pos = sort(union(positions(x), positions(y)))
    	new_opx = M[]
    	new_opy = M[]
    	for pos in new_pos
    		pos_x = findfirst(a->a==pos, positions(x))
    		pos_y = findfirst(a->a==pos, positions(y))
    		if isnothing(pos_x) && !(isnothing(pos_y))
    			push!(new_opx, id(Matrix{T}, x_left ⊗ space(opy[pos_y], 2) ))
    			push!(new_opy, opy[pos_y])
    			y_left = space(opy[pos_y], 3)'
    		elseif !(isnothing(pos_x)) && isnothing(pos_y)
    			push!(new_opx, opx[pos_x])
    			push!(new_opy, id(Matrix{T}, y_left ⊗ space(opx[pos_x], 2)))
    			x_left = space(opx[pos_x], 3)'
    		elseif !(isnothing(pos_x)) && !(isnothing(pos_y))
    			push!(new_opx, opx[pos_x])
    			push!(new_opy, opy[pos_y])
    			x_left = space(opx[pos_x], 3)'
    			y_left = space(opy[pos_y], 3)'
    		else
    			throw(ArgumentError("why here?"))
    		end
    	end
    	opx = new_opx
    	opy = new_opy
    	pos = new_pos
    end
    return pos, opx, opy
end