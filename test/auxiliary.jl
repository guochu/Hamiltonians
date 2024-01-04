println("------------------------------------")
println("auxiliary")
println("------------------------------------")


@testset "Coefficients" begin
	a = Coefficient(1.2)
	b = a * 3.5
	@test scalartype(a) <: Number
	@test scalartype(b) <: Number
	@test typeof(a + 1) == typeof(a)
	@test scalartype(1.3 + a) <: Real
	@test scalartype(a * 3.11) <: Real
	@test scalartype(3.11 * a) <: Real
	@test scalartype(a / 3.1) <: Real
	@test scalartype(3.1 / a) <: Real
	@test value(b) == value(a) * 3.5
	@test value(-a) == -value(a)
	@test a(0.1) == 1.2
	@test value(a / 2) == value(a) / 2
	@test value(3.7 / a) == 3.7 / value(a)
	@test conj(Coefficient(3.1+5*im)) == Coefficient(3.1-5*im)
	@test iszero(Coefficient(0))
	@test isconstant(a)
	f = Coefficient(sin)
	@test Coefficient(sin) == f
	@test !(f == a)
	@test !isconstant(f)
	@test f(1.7) == sin(1.7)
	@test (1 / f)(2.3) == 1 / f(2.3)
	@test (f / 2.4)(1.5) == f(1.5) / 2.4
	@test sqrt(f)(1.6) ≈ sqrt(f(1.6))
	@test conj(f(3+4*im)) ≈ conj(f)(3+4*im)
	@test (2.4 / f)(1.7) == 2.4 / f(1.7)
	@test (-f)(1.4) == -(f(1.4))
	g = Coefficient(cos)
	@test (f + g)(1.3) == f(1.3) + g(1.3)
	@test (f / g)(2.7) ≈ f(2.7) / g(2.7)
	@test (f - g)(1.4) ≈ f(1.4) - g(1.4)
end


@testset "Abelian matrix" begin
	p = boson_matrices(d=5)
	adag, a, n = p["+"], p["-"], p["n"]
	adag′ = convert(AbelianMatrix, adag)
	a′ = convert(AbelianMatrix, a)
	n′ = convert(AbelianMatrix, n)
	@test convert(Array, adag′) == adag
	@test convert(Array, a) == a
	@test convert(Array, n′) == n

end