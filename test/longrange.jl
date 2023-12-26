println("------------------------------------")
println("-----| Long-range Hamiltonian |-----")
println("------------------------------------")


function longrange_xxz(J, Jzz, hz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	C = [sp, sm, z]
	B = [2*J * sp', 2*J * sm', Jzz * z]
	terms = []
	for (a1, a2) in zip(C, B)
		push!(terms, ExponentialDecayTerm(a1, a2, α=exp(-α)))
	end
	return SchurMPOTensor(hz * z, [terms...])
end

function longrange_xxz_ham(L, hz, J, Jzz, α, p)
	sp, sm, z = p["+"], p["-"], p["z"]
	mpo = prodmpo(L, [1], [hz * z])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * z])
	end
	compress!(mpo)
	for i in 1:L
	    for j in i+1:L
	    	coeff = exp(-α*(j-i))
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sp, sp'])
	    	mpo += prodmpo(L, [i, j], [2*J*coeff*sm, sm'])
	    	mpo += prodmpo(L, [i, j], [Jzz*coeff*z, z])
	    	compress!(mpo)
	    end
	end
	return mpo
end

function longrange_xxz_mpoham(L, hz, J, Jzz, α, p)
	# the last term of J and Jzz not used
	mpo = MPOHamiltonian([longrange_xxz(J, Jzz, hz, α, p) for i in 1:L])
end

function longrange_fermion_mpoham(L, hz, J, alpha, p)
	adag, n, JW = p["+"], p["n"], p["JW"]
	@tensor adagJW[1,2;3,5] := adag[1,2,3,4] * JW[4,5]
	m = ExponentialDecayTerm(adagJW, adag', middle=JW, α=alpha, coeff=-J)
	t = SchurMPOTensor(hz * n, [m, m'])
	return MPOHamiltonian([t for i in 1:L])
end

function longrange_fermion_ham(L, hz, J, alpha, p)
	adag, n, JW = p["+"], p["n"], p["JW"]
	@tensor adagJW[1,2;3,5] := adag[1,2,3,4] * JW[4,5]
	a = adag'
	mpo = prodmpo(L, [1], [hz * n])
	for i in 2:L
		mpo += prodmpo(L, [i], [hz * n])
	end
	compress!(mpo)

	for i in 1:L
	    for j in i+1:L
	    	coeff = alpha^(j-i) 
	    	pos = collect(i:j)
	    	op_v = vcat(vcat([-J * coeff * adagJW], [JW for k in (i+1):(j-1)]), [a])
	    	tmp = prodmpo(L, pos, [op_v...])
	    	mpo += tmp
	    	mpo += tmp'
	    	compress!(mpo)
	    end
	end
	return mpo
end

function initial_state_u1_su2(::Type{T}, L) where {T<:Number}
	physpace = Rep[U₁×SU₂]((-0.5, 0)=>1, (0.5, 0)=>1, (0, 0.5)=>1)

	init_state = [(-0.5, 0) for i in 1:L]
	for i in 2:2:L
		init_state[i] = (0.5, 0)
	end
	n = sum([item[1] for item in init_state])
	n2 = 0
	right = Rep[U₁×SU₂]((n, 0)=>1)
	state = prodmps(T, physpace, init_state, right=right )

	return state
end


function initial_state_u1_u1(::Type{T}, L) where {T<:Number}
	physpace = Rep[U₁×U₁]((0, 0)=>1, (0, 1)=>1, (1, 0)=>1, (1, 1)=>1)

	init_state = [(0, 0) for i in 1:L]
	for i in 2:2:L
		init_state[i] = (1, 1)
	end
	n1 = sum([item[1] for item in init_state])
	n2 = sum([item[2] for item in init_state])

	right = Rep[U₁×U₁]((n1, n2)=>1)
	state = prodmps(T, physpace, init_state, right=right )
	return state
end

function do_dmrg(dmrg, alg)
	dmrg_sweeps = 10
	# Evals, delta = compute!(dmrg, alg)
	Evals = Float64[]
	for i in 1:dmrg_sweeps
		Evals, delta = sweep!(dmrg, alg)
	end
	return Evals[end]
end


@testset "MPOHamiltonian: long-range XXZ        " begin
	p = spin_site_ops_u1()
	for L in (2, 3, 4)
		hz = 0.8
		J = 1
		Jzz = 1.2
		α = 0.9
		h1 = longrange_xxz_mpoham(L, hz, J, Jzz, α, p)
		@test space_l(h1) == oneunit(space_l(h1))
		@test space_r(h1)' == oneunit(space_r(h1))
		@test length(h1) == L
		h2 = longrange_xxz_ham(L, hz, J, Jzz, α, p)
		@test length(h2) == L
		@test physical_spaces(h1) == physical_spaces(h2)
		right = iseven(L) ? Rep[U₁](0=>1) : Rep[U₁](1=>1)
		state = randommps(Float64, physical_spaces(h1), right=right, D=4)
		state = canonicalize!(state)
		@test expectation(h1, state) ≈ expectation(h2, state) atol = 1.0e-12
		@test distance(MPO(h1), h2) ≈ 0. atol = 1.0e-5		
	end
end

@testset "MPOHamiltonian: long-range fermions    " begin
	for p in (spinal_fermion_site_ops_u1_u1(), spinal_fermion_site_ops_u1_su2())
		for L  in (2, 3, 4, 5)
			hz = 0.7
			J = 1.
			alpha = 0.8
			h1 = longrange_fermion_mpoham(L, hz, J, alpha, p)
			@test space_l(h1) == oneunit(space_l(h1))
			@test space_r(h1)' == oneunit(space_r(h1))
			@test length(h1) == L
			h2 = longrange_fermion_ham(L, hz, J, alpha, p)
			@test length(h2) == L
			@test physical_spaces(h1) == physical_spaces(h2)
			if spacetype(h1) == Rep[U₁×U₁]
				right = Rep[U₁×U₁]((div(L, 2), div(L, 2)+1)=>1) 
			else 
				right = iseven(L) ? Rep[U₁×SU₂]((0,0)=>1) : Rep[U₁×SU₂]((0, 0.5)=>1)
			end
			state = randommps(ComplexF64, physical_spaces(h1), right=right, D=4)
			state = canonicalize!(state)
			@test expectation(h1, state) ≈ expectation(h2, state) atol = 1.0e-12
			@test distance(MPO(h1), h2) ≈ 0. atol = 1.0e-5		
		end
	end
end

@testset "DMRG with MPOHamiltonian: comparison with ED" begin
	J = 1.
	U = 1.5
	alpha = exp(-0.45)
	for L in 4:5
		for p in (spinal_fermion_site_ops_u1_u1(), spinal_fermion_site_ops_u1_su2())
			mpo = longrange_fermion_mpoham(L, J, U, alpha, p)
			observers = [prodmpo(L, [i], [p["n"]]) for i in 1:L]

			if spacetype(mpo) == Rep[U₁×U₁]
				state = initial_state_u1_u1(Float64, L)
			else
				state = initial_state_u1_su2(Float64, L)
			end	
			state = randommps(scalartype(state), physical_spaces(state), right=space_r(state)', D=10)
			state = canonicalize!(state, alg=Orthogonalize(normalize=true))

			# ED energy
			E, _st = exact_diagonalization(mpo, right=space_r(state)', num=1, ishermitian=true)
			E = E[1]

			E1 = do_dmrg(environments(mpo, copy(state)), DMRG2())
			@test E ≈ E1 atol = 1.0e-6

			E2 = do_dmrg(environments(mpo, copy(state)), DMRG1S())
			@test E ≈ E2 atol = 1.0e-6

			# check excited state
			dmrg = environments(mpo, copy(state))
			do_dmrg(dmrg, DMRG2())
			gs_state = dmrg.state

			E3, _st = exact_diagonalization(mpo, right=space_r(state)', num=2, ishermitian=true)
			@test E3[1] ≈ E atol = 1.0e-12
			E3 = E3[2]

			E4 = do_dmrg(environments(mpo, copy(state), [gs_state]), DMRG2(trunc=truncdimcutoff(D=20, ϵ=1.0e-10)))
			@test E3 ≈ E4 atol = 1.0e-6

			E5 = do_dmrg(environments(mpo, copy(state), [gs_state]), DMRG1S(trunc=truncdimcutoff(D=20, ϵ=1.0e-10)))
			@test E3 ≈ E5 atol = 1.0e-6

		end
	end
end
