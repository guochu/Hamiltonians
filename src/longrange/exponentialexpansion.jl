abstract type ExponentialExpansionAlgorithm end

# hankel expansion
struct HankelExpansion <: ExponentialExpansionAlgorithm 
    atol::Float64
    verbosity::Int
end
HankelExpansion(; atol::Real = 1.0e-8, verbosity::Int=1) = HankelExpansion(convert(Float64, atol), verbosity)
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

function hankel_exponential_expansion(fmat::AbstractMatrix)
    s1, n = size(fmat)
    (s1 >= n) || error("wrong input, try increase L, or decrease tol")
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
    err = norm(m * xs - fvec)
    # err = maximum(abs.(m * xs - fvec))
    # println("norm error is ", norm(m * xs - fvec), " ", err)
    return  xs, lambdas, err
end

hankel_exponential_expansion_n(f::Vector{<:Number}, n::Int) = hankel_exponential_expansion(generate_Fmat(f, n))
function exponential_expansion(f::Vector{<:Number}, alg::HankelExpansion)
    L = length(f)
    results = []
    errs = Float64[]
    atol = alg.atol
    verbosity = alg.verbosity
    for n in 1:L
        xs, lambdas, err0 = hankel_exponential_expansion_n(f, n)
        err = expansion_error(f, xs, lambdas)
        if err <= atol
            (verbosity > 1) && println("converged in $n iterations, error is $err.")
            return xs, lambdas
        else
            if (n > 1) && (err >= errs[end])
                (verbosity > 0) && println("stop at $n-th iteration due to error increase from $(errs[end]) to $err")
                return results[end]
            else
                push!(results, (xs, lambdas))
                push!(errs, err)
            end
        end
        if n >= L-n+1
            (verbosity > 0) && @warn "can not converge to $atol with size $L, try increase L, or decrease tol"
            return xs, lambdas
        end
    end
    error("can not find a good approximation")
end

# least square expansion
struct LsqExpansion <: ExponentialExpansionAlgorithm 
    atol::Float64
    verbosity::Int
end
LsqExpansion(; atol::Real = 1.0e-8, verbosity::Int=1) = LsqExpansion(convert(Float64, atol), verbosity)

function _predict(x, p)
    @assert length(p) % 2 == 0
    n = div(length(p), 2)
    L = length(x)
    T = eltype(p)
    r = zeros(T, L)
    for i in 1:L
        xi = x[i]
        @assert xi == i
        tmp = zero(T)
        for j in 1:n
            tmp += p[j] * p[n+j]^xi
        end
        r[i] = tmp
    end
    return r
end

function expansion_error(f::Vector{<:Number}, p::Vector{<:Number})
    T = eltype(f)
    xdata = [convert(T, i) for i in 1:length(f)]
    f_pred = _predict(xdata, p)
    return norm(f_pred - f)
end
expansion_error(f::Vector{<:Number}, coeffs::Vector{<:Number}, alphas::Vector{<:Number}) = expansion_error(f, vcat(coeffs, alphas))

function lsq_expansion_n(f::Vector{<:Real}, n::Int, coeffs::Vector{<:Real}, alphas::Vector{<:Real})
    @assert n == length(coeffs) == length(alphas)
    T = eltype(f)
    xdata = [convert(T, i) for i in 1:length(f)]
    p0 = vcat(coeffs, alphas)
    fit = curve_fit(_predict, xdata, f, p0, autodiff=:forwarddiff)
    # println("converged? ", fit.converged)
    p = fit.param
    err = norm(_predict(xdata, p) - f)
    return p[1:n], p[n+1:end], err
end
function exponential_expansion(f::Vector{<:Number}, alg::LsqExpansion)
    L = length(f)
    results = []
    errs = Float64[]
    atol = alg.atol
    verbosity = alg.verbosity
    for n in 1:L
        if n == 1
            xs, lambdas, err = lsq_expansion_n(f, n, [0.5], [0.5])
        else
            _xs, _lambdas = results[end]
            xs, lambdas, err = lsq_expansion_n(f, n, vcat(_xs, 0.5), vcat(_lambdas, 0.5))
        end
        if err <= atol
            (verbosity > 1) && println("converged in $n iterations, error is $err.")
            return xs, lambdas
        else
            push!(results, (xs, lambdas))
            push!(errs, err)
        end
        if n == L
            (verbosity > 0) && @warn "can not converge to $atol with size $L, try increase L, or decrease tol"
            return results[argmin(errs)]
        end
    end
    error("can not find a good approximation")
end

exponential_expansion(f::Vector{<:Number}; alg::ExponentialExpansionAlgorithm=HankelExpansion()) = exponential_expansion(f, alg)
exponential_expansion(f, L::Int, alg::ExponentialExpansionAlgorithm) = exponential_expansion([f(k) for k in 1:L], alg)
exponential_expansion(f, L::Int; alg::ExponentialExpansionAlgorithm=HankelExpansion()) = exponential_expansion(f, L, alg)
