"""
	struct QTerm{M <: MPOTensor} <: AbstractQuantumTerm 
"""
struct QTerm{M <: MPOTensor, C <: Coefficient} <: AbstractQuantumTerm 
	positions::Vector{Int}
	op::Vector{M}
	coeff::C


function QTerm(pos::Vector{Int}, m::Vector{M}, v::C) where {M <: MPOTensor, C <: Coefficient}
	# checks
	@assert !isempty(m)
	(length(pos) == length(m)) || throw(DimensionMismatch())
	check_qterm_positions(pos)
	DMRG.check_mpo_spaces(m)
	return new{M, C}(pos, m, v)
end
end

const ScalarQTerm{M} = QTerm{M, C} where {M <: MPOTensor, C <: ScalarCoefficient}

"""
	QTerm(pos::Vector{Int}, m::Vector{<:AbstractTensorMap}; coeff::Union{Number, Function, Coefficient}=1.) 
 	entrance to construct a single quantum term
"""
# entrance point
function QTerm(pos::Vector{Int}, m::Vector, v::Coefficient)
	@assert !isempty(m)
	L = length(m)
	S = spacetype(m[1])
	T = compute_scalartype(m)
	M = mpotensortype(S, T)
	mpotensors = Vector{M}(undef, L)
	left = oneunit(S)
	for i in 1:L
		mj = m[i]
		if isa(mj, MPSBondTensor)
			tmp = id(Matrix{T}, left)
			@tensor mj[-1 -2; -3 -4] := tmp[-1, -3] * mj[-2, -4]
		elseif isa(mj, AbelianMatrix)
			mj = tompotensor(mj, left=left)
		else
			isa(mj, MPOTensor) || throw(ArgumentError("unsupported operator type $(typeof(mj)) for QTerm"))
		end
		# @assert isa(mj, MPOTensor)
		mpotensors[i] = mj
		left = space(mj, 3)'
	end	
	return QTerm(pos, mpotensors, v)
end 

QTerm(pos::Vector{Int}, m::Vector; coeff::AllowedCoefficient=1.) = QTerm(pos, m, Coefficient(coeff))
QTerm(pos::Tuple, m::Vector; coeff::AllowedCoefficient=1.) = QTerm([pos...], m; coeff=coeff)
(x::QTerm)(t::Number) = QTerm(positions(x), op(x), coeff=coeff(x)(t))

QTerm(pos::Int, m::SiteOperator; coeff::AllowedCoefficient=1.) = QTerm([pos], [m], coeff=coeff)

TK.scalartype(::Type{QTerm{M, C}}) where {M<:MPOTensor, C <: ScalarCoefficient} = promote_type(scalartype(M), scalartype(C))
TK.spacetype(::Type{QTerm{M, C}}) where {M <: MPOTensor, C} = spacetype(M)

function QTerm(x::Pair{Int}...; coeff::AllowedCoefficient=1.) 
	pos, ms = _parse_pairs(x...)
	return QTerm(pos, ms; coeff=coeff)
end 


Base.copy(x::QTerm) = QTerm(copy(positions(x)), copy(op(x)), copy(coeff(x)))
Base.:(==)(x::QTerm, y::QTerm)  = (positions(x) == positions(y)) && (op(x) == op(y)) && (coeff(x) == coeff(y))

Base.:*(s::QTerm, m::AllowedCoefficient) = QTerm(positions(s), op(s), coeff(s) * Coefficient(m))

shift(m::QTerm, i::Int) = QTerm(positions(m) .+ i, op(m), coeff(m))
TK.id(x::QTerm) = QTerm(positions(x), [id(Matrix{scalartype(x)}, oneunit(spacetype(x)) ⊗ space(m, 2)) for m in op(x)], coeff=1.)


# Qterm adjoint
struct AdjointQTerm{M <: MPOTensor, C <: Coefficient} <: AbstractQuantumTerm 
	parent::QTerm{M, C}
end

Base.adjoint(m::QTerm) = AdjointQTerm(m)
Base.adjoint(m::AdjointQTerm) = m.parent
positions(m::AdjointQTerm) = positions(m.parent)
DMRG.coeff(m::AdjointQTerm) = conj(coeff(m.parent))
op(m::AdjointQTerm) = error("op for AdjointQTerm type is not supported")
TK.scalartype(::Type{AdjointQTerm{M, C}}) where {M<:MPOTensor, C<:ScalarCoefficient} = promote_type(scalartype(M), scalartype(C))


function Base.convert(::Type{<:QTerm}, m::AdjointQTerm) 
	isstrict(m.parent) || throw(ArgumentError("can not convert non-strict QTerm adjoint into a QTerm."))
	QTerm(positions(m.parent), DMRG.unsafe_mpotensor_adjoint.(op(m.parent)), conj(coeff(m.parent)) )
end 

function check_qterm_positions(pos::Vector{Int})
	(length(Set(pos)) == length(pos)) || throw(ArgumentError("duplicate positions not allowed"))
	(sort(pos) == pos) || throw(ArgumentError("QTerm positions should be strictly ordered."))
end

function _parse_pairs(x::Pair{Int}...) 
	pos = Int[]
	ms = []
	for (a, b) in x
		push!(pos, a)
		push!(ms, b)
	end
	return pos, ms
end

