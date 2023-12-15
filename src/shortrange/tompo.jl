

function DMRG.prodmpo(physpaces::Vector{S}, m::ScalarQTerm) where {S <: ElementarySpace}
	(S == spacetype(m)) || throw(SpaceMismatch())
	isconstant(m) || throw(ArgumentError("only constant term allowed"))
	return prodmpo(scalartype(m), physpaces, positions(m), op(m)) * scalar(coeff(m))
end

DMRG.prodmpo(physpaces::Vector{<: ElementarySpace}, m::AdjointQTerm) = adjoint(prodmpo(physpaces, m.parent))


function DMRG.MPO(h::QuantumOperator, alg::MPOCompression) 
	physpaces = physical_spaces(h)
	local mpo
	compress_threshold = 20
	for m in qterms(h)
		if (@isdefined mpo) && (!isnothing(mpo))
			mpo += prodmpo(physpaces, m)
		else
			mpo = prodmpo(physpaces, m)
		end
		if bond_dimension(mpo) >= compress_threshold
			mpo = compress!(mpo, alg)
			compress_threshold += 5
		end
	end
	(!(@isdefined mpo)) && throw(ArgumentError("Hamiltonian is empty after compression."))
	mpo = compress!(mpo, alg)
	return mpo
end
DMRG.MPO(h::QuantumOperator; alg::MPOCompression = Deparallelise()) = MPO(h, alg)
