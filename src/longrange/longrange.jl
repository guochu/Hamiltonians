function generate_Fmat(fvec::Vector{<:Number}, n::Int)
    L = length(fvec)
    (L >= n) || error("number of sites must be larger than number of terms in expansion")
    F = zeros(eltype(fvec), L-n+1, n)
    for j in 1:n
        for i in 1:L-n+1
            F[i, j] = fvec[i + j - 1]
        end
    end
    return F
end

generate_Fmat(f, L::Int, n::Int) = generate_Fmat([f(k) for k in 1:L], n)

function exponential_expansion(fmat::AbstractMatrix)
    s1, n = size(fmat)
    (s1 >= n) || error("wrong input, try increase L, or decrease the tolerance")
    L = s1 - 1 + n
    _u, _v = qr(fmat)
    U = Matrix(_u)
    V = Matrix(_v)
    U1 = U[1:L-n, :]
    U2 = U[(s1-L+n+1):s1, :]
    m = pinv(U1) * U2
    lambdas = eigvals(m)
    (length(lambdas) == n) || error("something wrong")
    m = zeros(eltype(lambdas), L, n)
    for j in 1:n
        for i in 1:L
            m[i, j] = lambdas[j]^i
        end
    end
    fvec = zeros(eltype(fmat), L)
    for i in 1:n
        fvec[i] = fmat[1, i]
    end
    for i in n+1:L
        fvec[i] = fmat[i-n+1, n]
    end
    xs = m \ fvec
    # err = norm(m * xs - fvec)
    err = maximum(abs.(m * xs - fvec))
    return  xs, lambdas, err
end

exponential_expansion_n(f, L::Int, n::Int) = exponential_expansion(generate_Fmat(f, L, n))
function exponential_expansion(f, L::Int; atol::Real=1.0e-5)
    for n in 1:L
        xs, lambdas, err = exponential_expansion_n(f, L, n)
        if err <= atol
            # println("converged $n iterations, error is $err.")
            return xs, lambdas
        end
        if n >= L-n+1
            @warn "can not converge to $atol with size $L, try increase L, or decrease the tolerance"
            return xs, lambdas
        end
    end
    error("can not find a good approximation")
end


abstract type AbstractLongRangeTerm end

DMRG.space_l(x::AbstractLongRangeTerm) = isa(x.a, MPOTensor) ? space_l(x.a) : oneunit(spacetype(x))
DMRG.space_r(x::AbstractLongRangeTerm) = isa(x.b, MPOTensor) ? space_r(x.b) : oneunit(spacetype(x))'
TK.spacetype(x::AbstractLongRangeTerm) = spacetype(typeof(x))
DMRG.coeff(x::AbstractLongRangeTerm) = x.coeff


# coeff * α^n, α must be in [0, 1]
struct ExponentialDecayTerm{M1<:SiteOperator, M<:SiteOperator, M2, T <:Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    α::T
    coeff::T
end

function ExponentialDecayTerm(a::SiteOperator, b::SiteOperator; middle::MPSBondTensor=id(physical_space(a)), α::Number=1., coeff::Number=1.) 
    T = promote_type(typeof(α), typeof(coeff))
    check_a_b(a, b)
    return ExponentialDecayTerm(a, middle, b, convert(T, α), convert(T, coeff))
end

TK.scalartype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T)
TK.scalartype(x::ExponentialDecayTerm) = scalartype(typeof(x))
TK.spacetype(::Type{ExponentialDecayTerm{M1, M, M2, T}}) where {M1, M, M2, T} = spacetype(M1)

Base.adjoint(x::ExponentialDecayTerm) = ExponentialDecayTerm(_op_adjoint(x.a, x.m, x.b)..., conj(x.α), conj(coeff(x)))
_op_adjoint(a::MPSBondTensor, m::MPSBondTensor, b::MPSBondTensor) = (a', m', b')
_op_adjoint(a::MPOTensor, m::MPSBondTensor, b::MPOTensor) = (DMRG.unsafe_mpotensor_adjoint(a), m', DMRG.unsafe_mpotensor_adjoint(b))


struct GenericDecayTerm{M1, M<:MPSBondTensor, M2, F, T <: Number} <: AbstractLongRangeTerm
    a::M1
    m::M
    b::M2
    f::F
    coeff::T
end

function GenericDecayTerm(a::SiteOperator, b::SiteOperator; 
                            middle::MPSBondTensor = id(physical_space(a)), f, coeff::Number=1.) 
    check_a_b(a, b)
    GenericDecayTerm(a, middle, b, f, coeff)
end 

TK.scalartype(x::GenericDecayTerm{M1, M, M2, F, T}) where {M1, M, M2, F, T} = promote_type(scalartype(M1), scalartype(M), scalartype(M2), T, typeof(x.f(0.)))
Base.adjoint(x::GenericDecayTerm) = GenericDecayTerm(_op_adjoint(x.a, x.m, x.b)..., y->conj(x.f(y)), conj(coeff(x)))
TK.spacetype(::Type{GenericDecayTerm{M1, M, M2, F, T}}) where {M1, M, M2, F, T} = spacetype(M1)

"""
    coeff * n^α, α must be negative (diverging otherwise)
"""
PowerlawDecayTerm(a::M, m::MPSBondTensor, b::M; α::Number=1., coeff::Number=1.) where {M<:SiteOperator} = GenericDecayTerm(a, m, b, f=x->x^α, coeff=coeff)
PowerlawDecayTerm(a::M, b::M; α::Number=1., coeff::Number=1.) where {M<:SiteOperator} = GenericDecayTerm(a, b, f=x->x^α, coeff=coeff)



# L is the number of sites

"""
    exponential_expansion(x::GenericDecayTerm{M, T, F}; len::Int, atol::Real=1.0e-5)

Convert a general decaying term into a list of exponential decaying term
"""
function exponential_expansion(x::GenericDecayTerm{M, T, F}; len::Int, atol::Real=1.0e-5) where {M, T, F}
    xs, lambdas = exponential_expansion(x.f, len, atol=atol)
    r = []
    for (c, alpha) in zip(xs, lambdas)
        push!(r, ExponentialDecayTerm(x.a, x.m, x.b; α=alpha, coeff=c * coeff(x)))
    end
    return [r...]
end

function _longrange_schurmpo_util(h1, h2s::Vector{<:ExponentialDecayTerm})
    isempty(h2s) && throw(ArgumentError("empty interactions."))
	pspace = physical_space(h2s[1].a)
	N = length(h2s)
	T = Float64
	for item in h2s
		T = promote_type(T, scalartype(item))
	end
	cell = Matrix{Any}(undef, N+2, N+2)
	for i in 1:length(cell)
		cell[i] = zero(T)
	end
	# diagonals
	cell[1, 1] = 1
	cell[end, end] = 1
	cell[1, end] = h1
	for i in 1:N
        if isa(h2s[i].a, MPSBondTensor)
            b_iden = id(Matrix{T}, oneunit(spacetype(h2s[i].a)))
        else
            b_iden = id(Matrix{T}, space_r(h2s[i].a)')
        end
        m = h2s[i].m
        @tensor iden[-1 -2; -3 -4] := b_iden[-1, -3] * m[-2, -4]
		cell[i+1, i+1] = h2s[i].α * iden
		cell[1, i+1] = h2s[i].coeff * h2s[i].a
		cell[i+1, end] = h2s[i].α * h2s[i].b
	end
	return SchurMPOTensor(cell)
end

function check_a_b(a::SiteOperator, b::SiteOperator)
    if isa(a, MPOTensor)
        isa(b, MPOTensor) || throw(ArgumentError("a and b must both be MPOTensor or MPSBondTensor"))
        s_l = space_l(a)
        (s_l == space_r(b)' == oneunit(s_l)) || throw(ArgumentError("only strict MPOTensor is allowed"))
    else
        isa(b, MPSBondTensor) || throw(ArgumentError("a and b must both be MPOTensor or MPSBondTensor"))
    end    
end

"""
    SchurMPOTensor(h1::ScalarSiteOp, h2s::Vector{<:ExponentialDecayTerm})
    SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm})

Return an SchurMPOTensor, with outer matrix size (N+2)×(N+2) (N=length(h2s))

Algorithm reference: "Time-evolving a matrix product state with long-ranged interactions"
"""
DMRG.SchurMPOTensor(h1::MPSBondTensor, h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(h1, h2s)
DMRG.SchurMPOTensor(h2s::Vector{<:ExponentialDecayTerm}) = _longrange_schurmpo_util(0., h2s)
