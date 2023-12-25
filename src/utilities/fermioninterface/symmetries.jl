

abstract type FermionicSymmetry end

abstract type SpinfulFermionicSymmetry <: FermionicSymmetry end
abstract type AbelianFermionicSymmetry <: SpinfulFermionicSymmetry end
abstract type NonAbelianFermionicSymmetry <: SpinfulFermionicSymmetry end


"""
	struct ChargeCharge

Two spins per site with U₁U₁ symmetry, namely the chagre of up and down spin are both conserved.
"""
struct ChargeCharge <: AbelianFermionicSymmetry end

"""
	struct SpinCharge

Two spins per site with U₁SU₂ symmetry, namely charge (U₁) and spin (SU₂) conservation.
"""
struct SpinCharge <: NonAbelianFermionicSymmetry end


abstract type FermionicSymmetrySector end

struct ChargeChargeSector <: FermionicSymmetrySector
	up::Int
	down::Int
end

ChargeChargeSector(;up::Int, down::Int) = ChargeChargeSector(up, down)

function Base.getproperty(x::ChargeChargeSector, s::Symbol)
	if s == :charge
		return x.up + x.down
	else
		getfield(x, s)
	end
end


struct SpinChargeSector <: FermionicSymmetrySector
	spin::Float64
	charge::Int

function SpinChargeSector(spin::Real, charge::Int)
	@assert _is_half_integer(spin)
	return new(convert(Float64, spin), charge)
end

end

SpinChargeSector(; spin::Real, charge::Int) = SpinChargeSector(spin, charge)

DMRG.space(sector::ChargeChargeSector) = Rep[U₁×U₁]((sector.up, sector.down)=>1)
DMRG.space(sector::SpinChargeSector) = Rep[U₁×SU₂]((sector.charge, sector.spin)=>1)


abstract type SpinlessFermionicSymmetry <: FermionicSymmetry end

"""
	struct Charge

Spinless fermion with charge conservation symmetry
"""
struct Charge <: SpinlessFermionicSymmetry end

abstract type SpinlessFermionicSymmetrySector end

struct ChargeSector <: SpinlessFermionicSymmetrySector
	charge::Int
end
ChargeSector(; charge::Int) = ChargeSector(charge)
DMRG.space(sector::ChargeSector) = Rep[U₁](sector.charge=>1)

function _is_half_integer(m::Real)
	mi = round(Int, 2*m)
	return abs(2*m - mi) < 1.0e-14
end

