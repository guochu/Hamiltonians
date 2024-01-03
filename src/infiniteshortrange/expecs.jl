"""
	expectation(m::QTerm, psi::InfiniteMPS)
The infinite MPS is assumed to be in the canonical form
"""
function DMRG.expectation(m::QTerm, psi::InfiniteMPS)
	DMRG.svectors_uninitialized(psi) && throw(ArgumentError("A canonical form infinite MPS is assumed"))
	return _expectation_canonical(m, psi)
end 


function DMRG.expectation(m::InfiniteQuantumOperator, psi::InfiniteMPS)
	m2 = changeunitcell(m, length(psi))
	r = 0.
	for t in qterms(m)
		r += expectation(t, psi)
	end
	return r
end