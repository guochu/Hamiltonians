println("------------------------------------")
println("|    Infinite Quantum operators    |")
println("------------------------------------")

function spin_site_ops_u1x()
    ph = Rep[U₁](-0.5=>1, 0.5=>1)
    vacuum = oneunit(ph)
    σ₊ = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
    blocks(σ₊)[Irrep[U₁](0.5)] = ones(1, 1)
    σ₋ = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph)
    blocks(σ₋)[Irrep[U₁](-0.5)] = ones(1, 1)
    σz = TensorMap(ones, ph ← ph)
    blocks(σz)[Irrep[U₁](-0.5)] = -ones(1, 1)
    return Dict("+"=>σ₊, "-"=>σ₋, "z"=>σz)
end

@testset "Infinite Expectation value: QTerm" begin

	ph = Rep[U₁](-0.5=>1, 0.5=>1)
	L = 4
	mps = randommps(Float64, [ph for i in 1:L], D=10)

	p = spin_site_ops_u1x()
	sp, sm, sz = p["+"], p["-"], p["z"]
	observers = [QTerm(i=>sp, (i+1)=>sp') for i in 1:L-1]

	tol = 1.0e-8

	obs1 = [expectation(ob, mps) for ob in observers]

	imps = convert(InfiniteMPS, mps)

	obs2 = [expectation(ob, imps) for ob in observers]

	@test max_error(obs1, obs2) < tol

	observers2 = [shift(item, L) for item in observers]

	obs3 = [expectation(ob, imps) for ob in observers2]

	@test max_error(obs1, obs3) < tol

end