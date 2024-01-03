
function boson_site_ops_u1(d::Int)
	@assert d > 1
	ph = Rep[U₁](i-1=>1 for i in 1:d)
	vacuum = oneunit(ph)
	adag = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
	for i in 1:d-1
		blocks(adag)[Irrep[U₁](i)] = sqrt(i) * ones(1,1)
	end
	a = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph )
	for i in 1:d-1
		blocks(a)[Irrep[U₁](i-1)] = sqrt(i) * ones(1, 1)
	end
	n = TensorMap(zeros, ph ← ph)
	for i in 1:d-1
		blocks(n)[Irrep[U₁](i)] = i * ones(1, 1)
	end
	return Dict("+"=>adag, "-"=>a, "n"=>n)
end

function boson_site_ops_u1_2(d::Int)
	@assert d > 1
	ph = Rep[U₁](i-1=>1 for i in 1:d)
	vacuum = oneunit(ph)
	adag = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
	for i in 1:d-1
		blocks(adag)[Irrep[U₁](i)] = sqrt(i) * ones(1,1)
	end
	a = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph )
	for i in 1:d-1
		blocks(a)[Irrep[U₁](i-1)] = sqrt(i) * ones(1, 1)
	end
	n = TensorMap(zeros, ph ← ph)
	for i in 1:d-1
		blocks(n)[Irrep[U₁](i)] = i * ones(1, 1)
	end
	return Dict("+"=>adag, "-"=>a, "n"=>n)
end

function spin_site_ops_u1()
    ph = Rep[U₁](0=>1, 1=>1)
    vacuum = oneunit(ph)
    σ₊ = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](1=>1) ⊗ ph)
    blocks(σ₊)[Irrep[U₁](1)] = ones(1, 1)
    σ₋ = TensorMap(zeros, vacuum ⊗ ph ← Rep[U₁](-1=>1) ⊗ ph)
    blocks(σ₋)[Irrep[U₁](0)] = ones(1, 1)
    σz = TensorMap(ones, ph ← ph)
    blocks(σz)[Irrep[U₁](0)] = -ones(1, 1)
    return Dict("+"=>σ₊, "-"=>σ₋, "z"=>σz)
end

"""
	The convention is that the creation operator on the left of the annihilation operator

By convention space_l of all the operators are vacuum
"""
function spinal_fermion_site_ops_u1_su2()
	ph = Rep[U₁×SU₂]((-0.5, 0)=>1, (0.5, 0)=>1, (0, 0.5)=>1)
	bh = Rep[U₁×SU₂]((0.5, 0.5)=>1)
	vh = oneunit(ph)
	adag = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
	blocks(adag)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
	blocks(adag)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = sqrt(2) * ones(1,1) 

	bh = Rep[U₁×SU₂]((-0.5, 0.5)=>1)
	a = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
	blocks(a)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
	blocks(a)[Irrep[U₁](-0.5) ⊠ Irrep[SU₂](0)] = -sqrt(2) * ones(1,1) 


	onsite_interact = TensorMap(zeros, Float64, ph ← ph)
	blocks(onsite_interact)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = ones(1, 1)

	JW = TensorMap(ones, Float64, ph ← ph)
	blocks(JW)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = -ones(1, 1)

	# adagJW = TensorMap(zeros, Float64, vh ⊗ ph ← bh ⊗ ph)
	# blocks(adagJW)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1,1)
	# blocks(adagJW)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = -sqrt(2) * ones(1,1) 

	# hund operators
	# c↑† ⊗ c↓†
	bhr = Rep[U₁×SU₂]((1, 0)=>1)
	adagadag = TensorMap(ones, Float64, vh ⊗ ph ← bhr ⊗ ph)

	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
	bhr = Rep[U₁×SU₂]((0, 1)=>1)
	adaga = TensorMap(zeros, Float64, vh ⊗ ph ← bhr ⊗ ph)
	blocks(adaga)[Irrep[U₁](0) ⊠ Irrep[SU₂](0.5)] = ones(1, 1) * (-sqrt(3) / 2)

	n = TensorMap(ones, Float64, ph ← ph)
	blocks(n)[Irrep[U₁](-0.5) ⊠ Irrep[SU₂](0)] = zeros(1, 1)
	blocks(n)[Irrep[U₁](0.5) ⊠ Irrep[SU₂](0)] = 2 * ones(1, 1)

	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, "JW"=>JW, "n"=>n)
end

function spinal_fermion_site_ops_u1_u1()
	ph = Rep[U₁×U₁]((1, 1)=>1, (1,0)=>1, (0,1)=>1, (0,0)=>1)
	vacuum = oneunit(ph)

	# adag
	adagup = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,0)=>1) ⊗ ph )
	blocks(adagup)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = ones(1,1)
	blocks(adagup)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = ones(1,1)

	adagdown = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,1)=>1) ⊗ ph)
	blocks(adagdown)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1)
	blocks(adagdown)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = -ones(1,1)

	adag = cat(adagup, adagdown, dims=3)

	# a
	aup = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,0)=>1) ⊗ ph)
	blocks(aup)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = ones(1,1)
	blocks(aup)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1)

	adown = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((0,-1)=>1) ⊗ ph)
	blocks(adown)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = ones(1,1)
	blocks(adown)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = -ones(1,1)

	a = cat(aup, - adown, dims=3)

	# hund operators
	adagadag = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,1)=>1) ⊗ ph)
	blocks(adagadag)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = ones(1, 1)

	# c↑† ⊗ c↓, this is a spin 1 sector operator!!!
	up = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((1,-1)=>1) ⊗ ph)
	blocks(up)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = ones(1,1) / (-sqrt(2))
	middle = TensorMap(zeros, Float64, vacuum ⊗ ph ← vacuum ⊗ ph )
	blocks(middle)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = 0.5 * ones(1,1)
	blocks(middle)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = -0.5 * ones(1,1)
	down = TensorMap(zeros, Float64, vacuum ⊗ ph ← Rep[U₁×U₁]((-1,1)=>1) ⊗ ph)
	blocks(down)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = ones(1,1) / sqrt(2)
	adaga = cat(cat(up, middle, dims=3), down, dims=3)

	onsite_interact = TensorMap(zeros, Float64, ph ← ph)
	blocks(onsite_interact)[Irrep[U₁](1) ⊠ Irrep[U₁](1)]= ones(1,1)

	JW = TensorMap(ones, Float64, ph ← ph)
	blocks(JW)[Irrep[U₁](1) ⊠ Irrep[U₁](0)] = -ones(1, 1)
	blocks(JW)[Irrep[U₁](0) ⊠ Irrep[U₁](1)] = -ones(1, 1)

	occupy = TensorMap(ones, Float64, ph ← ph)
	blocks(occupy)[Irrep[U₁](0) ⊠ Irrep[U₁](0)] = zeros(1, 1)
	blocks(occupy)[Irrep[U₁](1) ⊠ Irrep[U₁](1)] = 2 * ones(1, 1)
	return Dict("+"=>adag, "-"=>a, "++"=>adagadag, "+-"=>adaga, "n↑n↓"=>onsite_interact, 
		"JW"=>JW, "n"=>occupy)
end

max_error(a::Vector, b::Vector) = maximum(abs.(a - b))
