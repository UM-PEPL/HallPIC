module HallPIC

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


function push_pos(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.pos)
        pc.pos[i] = muladd(pc.vel[i], dt, pc.pos[i])
    end
end

function push_vel(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.vel)
        pc.vel[i] = muladd(pc.acc[i], dt, pc.vel[i])
    end
end

function push(pc::ParticleContainer, dt::AbstractFloat)
    """
    Push particle container to next timestep using Leapfrog scheme
    """
    push_vel(pc, dt)
    push_pos(pc, dt)
end

end
