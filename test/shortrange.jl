println("------------------------------------")
println("|       Quantum operators          |")
println("------------------------------------")
@testset "Quantum operators: QTerm" begin
	# u1 symmetry
	p = spin_site_ops_u1()
	sp, sm, z = p["+"], p["-"], p["z"]
	ph = space(sp, 2)
	physpaces = [ph for i in 1:5]
	t2 = QTerm(2=>z, 3=>z, coeff=0.5)
	@test bond_dimension(t2) == 1
	t = QTerm(2=>sp, 4=>sp', coeff=0.6)
	@test positions(t) == [2, 4]
	@test scalar(coeff(t)) == 0.6
	@test scalartype(t) == Float64
	@test spacetype(t) == spacetype(sp)
	@test isconstant(t)
	@test !iszero(t)
	@test space_l(t) == oneunit(space_l(t))
	@test space_r(t)' == oneunit(space_r(t))
	@test bond_dimension(t) == 1
	h1 = prodmpo(physpaces, t)
	@test distance(h1, 0.6*prodmpo(physpaces, [2, 4], [sp, sp'])) ≈ 0. atol=1.0e-6
	t3 = t2 * t
	@test positions(t3) == [2,3,4]
	@test scalar(coeff(t3)) == 0.3
	@test isconstant(t3)
	# t3 = t ⊗ t2
	# @test positions(t3) == [2,3,4]
	# @test scalar(coeff(t3)) == 0.3
	# t3 = t ⊠ t2
	# @test positions(t3) == [2,3,4]
	# @test scalar(coeff(t3)) == 0.3

end

function twobody(pos1::Int, pos2::Int, c::Number, p)
	adag = a_dagger(pos1, p)
	a = a_dagger(pos2, p)'
	r = adag * a
	return simplify(r) * Coefficient(c)
end

function fourbody(pos1::Int, pos2::Int, pos3::Int, pos4::Int, v::Number, p)
	adag1 = a_dagger(pos1, p)
	adag2 = a_dagger(pos2, p)
	a3 = a_dagger(pos3, p)'
	a4 = a_dagger(pos4, p)'

	r = (adag1 * adag2) * (a3 * a4)
	return simplify(r) * (0.5*Coefficient(v))
end


function a_dagger(pos::Int, p)
	adag, JW = p["+"], p["JW"]
	_pos = collect(1:pos)
	_ops = Vector{Any}(undef, length(_pos)) 
	# _ops[1:end-1] .= JW
	for i in 1:length(_pos)-1
		_ops[i] = JW
	end
	_ops[end] = adag
	return QTerm(_pos, [_ops...])
end

@testset "QTerm: multiplication" begin
	for p in (spinal_fermion_site_ops_u1_u1(), spinal_fermion_site_ops_u1_su2())
		JW, adag = p["JW"], p["+"]
		a = QTerm([1,2], [JW, adag])
		b = QTerm([1,2,3,4,5,6], [JW, JW, JW, JW, JW, adag])
		c = QTerm([1,2,3,4], [JW, JW, JW, adag])'
		d = QTerm([1,2,3,4,5], [JW, JW, JW, JW, adag])'
		ph = space(adag, 2)
		for L in (6, 7)
			spaces = [ph for i in 1:L]
			ma = prodmpo(spaces, a)
			mb = prodmpo(spaces, b)
			mc = prodmpo(spaces, c)
			md = prodmpo(spaces, d)
			r1 = (ma * mb) * (mc * md)
			r2 = (ma * md) * (mb * mc)
			@test distance(r1, r2) < 1.0e-4
		end
	end

	for p in (spinal_fermion_site_ops_u1_u1(), spinal_fermion_site_ops_u1_su2())
		# twobody terms
		m = twobody(2,2, 1, p)
		@test positions(m) == [2]
		m = twobody(2,3, 0.5, p)
		@test positions(m) == [2,3]
		m = twobody(2,4, 0.5, p)
		@test positions(m) == [2,3,4]
		# four body terms
		m = fourbody(2,2,2,2, 0.4, p)
		@test positions(m) == [2]
		m = fourbody(2,4,4,2, 0.4, p)
		@test positions(m) == [2,4]
		m = fourbody(2,5,4,2, 0.4, p)
		@test positions(m) == [2,4,5]
		m = fourbody(2,5,5,4, 0.4, p)
		@test positions(m) == [2,3,4,5]
		m = fourbody(4,5,7,2, 0.4, p)
		@test positions(m) == [2,3,4,5,6,7]
	end
end

@testset "QTerm: mpo conversion" begin
	for p in (spinal_fermion_site_ops_u1_u1(), spinal_fermion_site_ops_u1_su2())
		_ph = space(p["+"], 2)
		physpaces = [_ph for i in 1:8]
		if spacetype(_ph) == Rep[U₁×U₁]
			mps = randommps(physpaces, D=8, right = Rep[U₁×U₁]((4,4)=>1))
		else
			mps = randommps(physpaces, D=8)
		end
		mps = canonicalize!(mps, alg=Orthogonalize(SVD(), normalize=true))
		ms = [twobody(2,2, 1, p), twobody(2,3, 0.5, p), twobody(2,4, 0.5, p), fourbody(2,2,2,2, 0.4, p), fourbody(2,4,4,2, 0.4, p), 
		fourbody(2,5,4,2, 0.4, p), fourbody(2,3,5,4, 0.4, p), fourbody(4,5,7,2, 0.4, p)]
		for m in ms
			mpo = prodmpo(physpaces, m)
			@test mpo * mps ≈ m * mps atol = 1.0e-6
		end

		@test prodmpo(physpaces, ms[1] + ms[end]) ≈ prodmpo(physpaces, ms[1]) + prodmpo(physpaces, ms[end]) atol = 1.0e-6
		for m in ms
			mi = id(m)
			@test prodmpo(physpaces, m + mi) ≈ prodmpo(physpaces, m) + prodmpo(physpaces, mi) atol = 1.0e-6
 		end
	end
end

# function bosonhubbard_ham(L::Int, J::Real, U::Real, mu::Real)
# 	p = boson_site_ops_u1(4)
# 	adag, a, n = p["+"], p["-"], p["n"]
# 	n2 = n * n
# 	terms = []
# 	for i in 1:L
# 		push!(terms, QTerm(i=>(n*(mu-U/2.)+n2*(U/2.))))
# 	end
# 	ham = QuantumOperator(terms)
# 	for i in 1:L-1
# 		t = QTerm(i=>adag, i+1=>adag', coeff=-J)
# 		add!(ham, t)
# 		add!(ham, t')
# 	end
# 	return ham
# end



