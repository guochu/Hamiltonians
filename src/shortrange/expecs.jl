
"""
	assume the underlying state is canonical
"""
function _expectation_canonical(m::QTerm, psi)
	isstrict(m) || throw(ArgumentError("QTerm should conserve quantum number"))
	isconstant(m) || throw(ArgumentError("only constant QTerm allowed"))
	iszero(m) && return 0.
	pos = positions(m)
	ops = op(m)
	pos_end = pos[end]
	util = get_trivial_leg(psi[1])
	@tensor hold[-3 -2; -1] := conj(psi[pos_end][-1, 1, 2]) * psi[pos_end][-3, 4, 2] * ops[end][-2, 1, 3, 4] * util[3] 
	for j in pos_end-1:-1:pos[1]
		pj = findfirst(x->x==j, pos)
		if isnothing(pj)
			hold = updateright(hold, psi[j], pj, psi[j])
		else
			hold = updateright(hold, psi[j], ops[pj], psi[j])
		end
	end	 
	s = convert(TensorMap, psi.s[pos[1]]) 
	@tensor hnew[-1; -2] := conj(s[-1, 1]) * hold[3, 2, 1] * conj(util[2]) * s[-2, 3]
	return tr(hnew) * value(coeff(m))
end

function expectation_canonical(m::QTerm, psi::MPS)
	(positions(m)[end] <= length(psi)) || throw(BoundsError())
	return _expectation_canonical(m, psi)
end

function DMRG.expectation(psiA::AbstractMPS, m::QTerm, psiB::AbstractMPS, envs::OverlapCache=environments(psiA, psiB))
	cstorage = envs.cstorage
	(length(psiA) == length(psiB) == length(cstorage)-1) || throw(DimensionMismatch())
	isstrict(m) || throw(ArgumentError("onlt strict QTerm is allowed"))
	isconstant(m) || throw(ArgumentError("only constant QTerm allowed"))
	iszero(m) && return 0.
	L = length(psiA)
	pos = positions(m)
	ops = op(m)
	pos_end = pos[end]
	(pos_end <= L) || throw(BoundsError())
	util = get_trivial_leg(psiA[1])
	@tensor hold[-3 -2; -1] := conj(psiA[pos_end][-1, 1, 2]) * cstorage[pos_end+1][3, 2] * psiB[pos_end][-3, 5, 3] * ops[end][-2, 1, 4, 5] * util[4]  
	for j in pos_end-1:-1:pos[1]
		pj = findfirst(x->x==j, pos)
		if isnothing(pj)
			hold = updateright(hold, psiA[j], pj, psiB[j])
		else
			hold = updateright(hold, psiA[j], ops[pj], psiB[j])
		end
	end
	@tensor hnew[-2; -1] := conj(util[1]) * hold[-2, 1, -1]
	for j in pos[1]-1:-1:1
		hnew = updateright(hnew, psiA[j], psiB[j])
	end
	return scalar(hnew) * value(coeff(m))
end

DMRG.expectation(m::QTerm, psi::MPS; iscanonical::Bool=false) = iscanonical ? expectation_canonical(m, psi) : expectation(psi, m, psi)
DMRG.expectation(m::QTerm, psi::ExactMPS) = expectation(psi, m, psi)

function DMRG.expectation(psiA::AbstractMPS, h::QuantumOperator, psiB::AbstractMPS)
	(length(h) <= length(psiA)) || throw(DimensionMismatch())
	envs = environments(psiA, psiB)
	r = 0.
	for m in qterms(h)
		r += expectation(psiA, m, psiB, envs)
	end
	return r
end
function DMRG.expectation(h::QuantumOperator, psi::MPS; iscanonical::Bool=false)
	if iscanonical
		r = 0.
		for m in qterms(h)
			r += expectation_canonical(m, psi)
		end
		return r
	else
		return expectation(psi, h, psi)
	end
end
DMRG.expectation(h::QuantumOperator, psi::ExactMPS) = expectation(psi, h, psi)

