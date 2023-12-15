
# function compute_eltype(ms)
# 	T = Float64
# 	for v in ms
# 		T = promote_type(T, eltype(v))
# 	end
# 	return T	
# end

# function _is_abelian_sector(::Type{C}) where {C <: Sector}
# 	(C <: AbelianIrrep) && return true
# 	if (C <: ProductSector)
# 		symms = C.parameters[1].parameters
# 		for item in symms
# 			(item <: AbelianIrrep) || return false
# 		end
# 		return true
# 	end
# 	return false
# end

"""
	struct AbelianMatrix{S <: ElementarySpace, T <: Number, Ic <: Sector} 
Struct to hold matricies which are labeled with U(1) indices, note that the quantum numbers does not have to be conserved, compared to TensorMap.
This is for convenience to construct U(1) MPO tensors
"""
struct AbelianMatrix{S <: ElementarySpace, T <: Number, Ic <: Sector} 
	data::Dict{Tuple{Ic, Ic}, Matrix{T}}
	physpace::S
end

function AbelianMatrix{S, T, C}(physpace::S, data::AbstractDict) where {S <: GradedSpace, T <: Number, C <: Sector}
	isdual(physpace) && throw(ArgumentError("dual space not allowed"))
	(FusionStyle(C) isa UniqueFusion) || throw(ArgumentError("Abelian sector expected"))
	(sectortype(S) == C) || throw(SectorMismatch())
	new_data = Dict{Tuple{C, C}, Matrix{T}}()
	for (k, vo) in data
		k1 = convert(C, k[1])
		k2 = convert(C, k[2])
		v = convert(Matrix{T}, vo)
		(hassector(physpace, k1) && hassector(physpace, k2)) || throw(ArgumentError("sector does not exist."))
		((dim(physpace, k1) == size(v, 1)) && (dim(physpace, k2) == size(v, 1))) || throw(SpaceMismatch())
		if !iszero(v)
			new_data[(k1, k2)] = v
		end
		# new_data[(k1, k2)] = v
	end
	return AbelianMatrix{S, T, C}(new_data, physpace)
end 
AbelianMatrix(::Type{T}, physpace::S, data::AbstractDict) where {S <: GradedSpace, T <: Number} = AbelianMatrix{S, T, sectortype(S)}(physpace, data)
AbelianMatrix(physpace::S, data::AbstractDict) where {S <: GradedSpace} = AbelianMatrix(compute_scalartype(values(data)), physpace, data)


DMRG.physical_space(m::AbelianMatrix) = m.physpace
TK.spacetype(::Type{<:AbelianMatrix{S}}) where {S<:ElementarySpace} = S
TK.spacetype(x::AbelianMatrix) = spacetype(typeof(x))
TK.scalartype(::Type{<:AbelianMatrix{S, T}}) where {S<:ElementarySpace, T<:Number} = T
raw_data(m::AbelianMatrix) = m.data


Base.copy(m::AbelianMatrix) = AbelianMatrix(physical_space(m), deepcopy(raw_data(m)))
Base.similar(m::AbelianMatrix) = AbelianMatrix(physical_space(m), typeof(raw_data(m))() )
Base.eltype(::Type{AbelianMatrix{S, T, I}}) where {S, T, I} = T
Base.eltype(x::AbelianMatrix) = eltype(typeof(x))

Base.:*(m::AbelianMatrix, y::Number) = AbelianMatrix(Dict(k=>v*y for (k, v) in raw_data(m)), physical_space(m))
Base.:*(y::Number, m::AbelianMatrix) = m * y
Base.:/(m::AbelianMatrix, y::Number) = m * (1 / y)
Base.:+(m::AbelianMatrix) = m
Base.:-(m::AbelianMatrix) = m * (-1)
function Base.:+(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace}
	(physical_space(x) == physical_space(y)) || throw(SpaceMismatch())
	T = promote_type(eltype(x), eltype(y))
	return AbelianMatrix(T, physical_space(x), merge(+, raw_data(x), raw_data(y)))
end
Base.:-(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace} = x + (-1) * y

function Base.:*(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace}
	T = promote_type(eltype(x), eltype(y))
	Ic = sectortype(S)
	r = Dict{Tuple{Ic, Ic}, Matrix{T}}()
	for (kx, vx) in raw_data(x)
		for (ky, vy) in raw_data(y)
			if kx[2] == ky[1]
				kr = (kx[1], ky[2])
				m = get!(r, kr, zeros(T, size(vx, 1), size(vy, 2)))
				m .+= vx * vy
			end
		end
	end
	return AbelianMatrix(T, physical_space(x), r)
end
Base.adjoint(x::AbelianMatrix) = AbelianMatrix(physical_space(x), Dict((k[2], k[1])=>v' for (k, v) in raw_data(x)))

function Base.:(==)(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace} 
	(physical_space(x) == physical_space(y)) || throw(SpaceMismatch())
	return raw_data(x) == raw_data(y)
end

function TK.tr(x::AbelianMatrix)
	r = zero(eltype(x))
	for (k, v) in x
		if k[1] == k[2]
			r += tr(v)
		end
	end
	return r
end

TK.dot(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace} = tr(x' * y)
TK.norm(x::AbelianMatrix) = sqrt(real(dot(x, x)))


"""
	abelian_matrix_from_dense(m::AbstractMatrix{T}) where {T <: Number}

Convert a dense matrix into an Abelian matrix with U(1) symmetry
"""
function Base.convert(::Type{AbelianMatrix}, m::AbstractMatrix{T}) where {T <: Number}
	(size(m, 1) == size(m, 2)) || throw(ArgumentError("square matrix ecpected."))
	physpace = U1Space(i-1=>1 for i in 1:size(m, 1))
	data = Dict{Tuple{Int, Int}, Matrix{T}}()
	for i in 1:size(m, 1)
		for j in 1:size(m, 2)
			if m[i, j] != zero(T)
				data[(i-1, j-1)] = ones(1, 1) * m[i, j]
			end
		end
	end
	return AbelianMatrix(physpace, data)
end
function Base.convert(::Type{Array}, m::AbelianMatrix{S, T}) where {S, T}
	physpace = physical_space(m)
	L = dim(physpace)
	m2 = zeros(T, L, L)
	for (k, v) in raw_data(m)
		m2[axes(physpace, k[1]), axes(physpace, k[2])] = v
	end
	return m2
end

function Base.kron(x::AbelianMatrix{S}, y::AbelianMatrix{S}) where {S <: GradedSpace}
	physpace = fuse(physical_space(x) ⊠ physical_space(y))
	C = sectortype(physpace)
	T = promote_type(eltype(x), eltype(y))
	data = Dict{Tuple{C, C}, Matrix{T}}()
	for (kx, vx) in raw_data(x)
		for (ky, vy) in raw_data(y)
			k = (kx[1] ⊠ ky[1], kx[2] ⊠ ky[2])
			data[k] = kron(vx, vy)
		end
	end
	return AbelianMatrix(physpace, data)
end

Base.zero(x::AbelianMatrix) = similar(x)
function Base.iszero(x::AbelianMatrix)
	isempty(raw_data(x)) && return true
	for (k, v) in raw_data(x)
		iszero(v) || return false
	end
	return true
end 
function Base.one(x::AbelianMatrix)
	physpace = physical_space(x)
	data = typeof(raw_data(x))()
	T = eltype(x)
	for s in sectors(physpace)
		d = dim(physpace, s)
		data[(s, s)] =  one(zeros(T, d, d))
	end
	return AbelianMatrix(physpace, data)
end

function tompotensor(m::AbelianMatrix{S}; left::S=oneunit(S)) where {S <: GradedSpace}
	phy = physical_space(m)
	right_sectors = sectortype(S)[]
	for (k, v) in raw_data(m)
		ko, ki = k
		for _l in sectors(left)
			(dim(left, _l) == 1) || throw(SpaceMismatch("wrong leftind"))
			_r, = first(ko ⊗ _l) ⊗ conj(ki)
			if !(_r in right_sectors)
				push!(right_sectors, _r)
			end
		end
	end
	right = S(item=>1 for item in right_sectors)
	r = TensorMap(zeros, left ⊗ phy ← right ⊗ phy )
	for (k, v) in raw_data(m)
		ko, ki = k
		for _l in sectors(left)
			_r, = first(ko ⊗ _l) ⊗ conj(ki)
			copyto!(r[(_l, ko, conj(_r), conj(ki))], v)
		end
	end
	return r
end
