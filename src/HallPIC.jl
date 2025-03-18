module HallPIC

using DataInterpolations: DataInterpolations, LinearInterpolation

const q_e = 1.608e-19
const N_A = 6.022e26
const m_0 = 1 / N_A
const q0_m0 = q_e / m_0
const phi_0 = 300.0
const x_0 = 0.05
const E_0 = phi_0 / x_0
const u_0 = sqrt(q_e * phi_0 / m_0)
const t_0 = x_0 / u_0

struct Particle{T<:AbstractFloat}
	weight::T
	pos::T
	vel::T
	acc::T
end

struct ParticleContainer{T<:AbstractFloat, I<:Integer}
    """
    T: Float32 or Float64
    I: UInt8 (if less than 256 rxns) or UInt16
    """
    weight::Vector{T}
    pos::Vector{T}
    vel::Vector{T}
    acc::Vector{T}
    reactions::Vector{I}
    mass::T
    charge::UInt8
end

function ParticleContainer(::Type{T}, N, mass, charge) where T
	weight = zeros(T, N)
	pos = zeros(T, N)
	vel = zeros(T, N)
	acc = zeros(T, N)
	rxn = UInt8[]
	return ParticleContainer{T, eltype(rxn)}(
		weight, pos, vel, acc, rxn, T(mass), eltype(rxn)(charge)
	)
end

function push_pos!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.pos)
        pc.pos[i] = muladd(pc.vel[i], dt, pc.pos[i])
    end
end

function push_vel!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.vel)
        pc.vel[i] = muladd(pc.acc[i], dt, pc.vel[i])
    end
end

function push!(pc::ParticleContainer, dt::AbstractFloat)
    """
    Push particle container to next timestep using Leapfrog scheme
    """
    push_vel!(pc, dt)
    push_pos!(pc, dt)
end

function gather!(pc::ParticleContainer, E_itp::LinearInterpolation)
	charge_to_mass = pc.charge / pc.mass
	@inbounds for i in eachindex(pc.pos)
		pc.acc[i] = charge_to_mass * E_itp(pc.pos[i])
	end
end

end
