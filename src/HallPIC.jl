module HallPIC

using DataInterpolations: DataInterpolations, LinearInterpolation
using DocStringExtensions
using Random: Random


#===================================================
Nondimensionalization notes

# Base quantities
m_0 = 1 amu 		[mass]
q_0 = q_e 			[charge]
n_0 = 1e14			[volume -> length]
phi_0 = q T_0 = q E_0 = 1 V or 1 eV  [energy -> time]

x_0 = 1 / (n_0^1/3)
u_0 = sqrt(q_e * phi_0 / m_0)
t_0 = x_0 / u_0
E_0 = phi_0 / x_0
===================================================#

# Physical constants and base units
const q_e = 1.608e-19
const N_A = 6.022e26
const m_0 = 1 / N_A
const n_0 = 1e12
const phi_0 = 1.0

# Derived quantities
const x_0 = 1 / cbrt(n_0)
const u_0 = sqrt(q_e * phi_0 / m_0)
const t_0 = x_0 / u_0
const E_0 = phi_0 / x_0

#======================================================
Gas and Species definitions
======================================================#

@kwdef struct Gas
	mass::Float64
	name::Symbol
end

"""
$(TYPEDEF)
A single chemical species (charge + mass)
$(TYPEDFIELDS)
"""
@kwdef struct Species
	gas::Gas
    reactions::Vector{UInt8}
	charge::UInt8
end

@inline (g::Gas)(Z::Integer) = Species(gas=g, charge=UInt8(Z), reactions=UInt8[])

# Helper methods
@inline mass(s::Species) = s.gas.mass
@inline charge(s::Species) = s.charge
@inline charge_to_mass(s::Species) = s.charge / s.gas.mass

#======================================================
ParticleContainer definitions
======================================================#

"""
$(TYPEDEF)
Per-particle information for a single species.
`T` is expected to be either Float32 or Float64
$(TYPEDFIELDS)
"""
struct ParticleContainer{T<:AbstractFloat}
    weight::Vector{T}
    pos::Vector{T}
    vel::Vector{T}
    acc::Vector{T}
	species::Species
end

function ParticleContainer{T}(N, species::Species) where T<:AbstractFloat
	weight = zeros(T, N)
	pos = zeros(T, N)
	vel = zeros(T, N)
	acc = zeros(T, N)
	return ParticleContainer{T}(
		weight, pos, vel, acc, species 
	)
end

"""
$(TYPEDSIGNATURES)
Push particle container to next timestep using Leapfrog scheme
"""
function push!(pc::ParticleContainer, dt::AbstractFloat)
    push_vel!(pc, dt)
    push_pos!(pc, dt)
end

"""
$(TYPEDSIGNATURES)
Add particles to a `ParticleContainer`.
"""
function add_particles!(pc::ParticleContainer{T}, x::Vector{T}, v::Vector{T}, w::Vector{T}) where T
	# check that new arrays have same length
	M = length(x)
	N = length(pc.pos)
	@assert M == length(v) && M == length(w)
	# append position, velocity, weight to pc arrays
	append!(pc.pos, x)
	append!(pc.vel, v)
	append!(pc.weight, w)
	# extend acceleration array to correct size and fill with zeros
	resize!(pc.acc, N+M)
	for i in 1:M
		pc.acc[N+i] = zero(T)
	end
	return pc
end

function push_pos!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.pos)
        pc.pos[i] = muladd(pc.vel[i], dt, pc.pos[i])
    end
	return pc
end

function push_vel!(pc::ParticleContainer, dt::AbstractFloat)
    @inbounds for i in eachindex(pc.vel)
        pc.vel[i] = muladd(pc.acc[i], dt, pc.vel[i])
    end
	return pc
end


function gather!(pc::ParticleContainer{T}, E_itp::LinearInterpolation) where T
	q_m = charge_to_mass(pc.species)
	@inbounds for i in eachindex(pc.pos)
		pc.acc[i] = T(q_m * E_itp(pc.pos[i]))
	end
	return pc
end


#======================================================
SpeciesProperties and CellProperties definitions
======================================================#

"""
$(TYPEDEF)
T: Float32 or Float64
Holds fluid properties (density, velocity, temperature) for a single species.
$(TYPEDFIELDS)
"""
struct SpeciesProperties{T<:AbstractFloat}
    dens::Vector{T}
    vel::Vector{T}
    temp::Vector{T}
    avg_weight::Vector{T}
	species::Species
end

"""
$(TYPEDEF)
Contains fluid information for a single cell.
$(TYPEDFIELDS)
"""
struct CellProperties{T<:AbstractFloat}
	dens::T
	vel::T
	temp::T
	avg_weight::T
end

function SpeciesProperties{T}(num_cells, species::Species) where T<:AbstractFloat
	dens = zeros(T, num_cells)
	vel = zeros(T, num_cells)
	temp = zeros(T, num_cells)
	avg_weight = zeros(T, num_cells)
	return SpeciesProperties{T}(dens, vel, temp, avg_weight, species)
end


#======================================================
  Iterator interface definition for SpeciesProperties
======================================================#

Base.length(sp::SpeciesProperties) = length(sp.dens)

function Base.getindex(sp::SpeciesProperties{T}, i) where T
	return CellProperties{T}(
		sp.dens[i], sp.vel[i], sp.temp[i], sp.avg_weight[i]
	)
end

Base.isdone(sp::SpeciesProperties) = i > length(sp)
Base.eltype(::SpeciesProperties{T}) where T = CellProperties{T}

Base.iterate(sp::SpeciesProperties) = length(sp) > 0 ? (sp[1], 2) : nothing 
function Base.iterate(sp::SpeciesProperties, i)
	i > length(sp) && return nothing
	return sp[i], i+1
end

#======================================================
Functions using SpeciesProperties
======================================================#

"""
$(TYPEDSIGNATURES)
Create a particle container and fill it with particles matching the macroscopic properties
from a provided SpeciesProperties object.

Accounting note: the edges of cell i are edges[i] and edges[i+1]
"""
function initialize_particles(sp::SpeciesProperties{T}, edges, volumes, particles_per_cell) where T
	pos_buf = zeros(T, N_ppc)
	vel_buf = zeros(T, N_ppc)
	weight_buf = zeros(T, N_ppc)

	num_cells = length(cell_centers)

	@assert num_cells == length(sp)

	pc = ParticleContainer{T}(0, species)

	for (i, (V, cell)) in enumerate(zip(volumes, sp))
		z_L = edges[i]
		z_R = edges[i+1]
		dz = z_R - z_L

		Random.rand!(pos_buf)
		@. pos_buf = dz * pos_buf + z_L 

		Random.randn!(vel_buf)
		v_th = sqrt(cell.temp / pc.gas.m)
		@. vel_buf + vel_buf * v_th + cell.vel

		weight = cell.dens * V / particles_per_cell
		@. weight_buf = weight

		add_particles!(pc, pos_buf, vel_buf, weight_buf)
	end

	return pc
end

#=

function Deposit(N_cell, x, dx, w_bar, Particles::ParticleContainer)
    """
    Deposit bulk particle quantities to the grid 
    Inputs: 
        N_cell (int)
            number of cells on the grid
        x (N_cell+1 array of floats)
            positions of cell centers 
        dx (N_cell array of floats)
            cell widths
        w_bar (N_cell array of floats)
            time-averaged average particle weight
        Particles (Particle Object)
            object containing particle information 
    Outputs:
        density (N_cell array of floats)
            number density for each cell
        velocity (N_cell array of floats)
            velocity for each cell
        temperature (N_cell array of floats)
            temperature for each cell
        w_bar (N_cell array of floats)
            time-averaged average particle weight
        N (N_cell array of floats)
            number of particles in each cell 
    """

    #initialize output
    w_sum = zeros(N_cell)
    v_sum = zeros(N_cell)
    T_sum = zeros(N_cell)    
    N = zeros(N_cell)

    #pull quantities
    w_p = Particles.weight
    v_p = Particles.vel
    x_p = Particles.pos

    #loop over particles
    @inbounds for i in eachindex(Particles.pos)
        #calculate relative position 
        x_rel = abs.(x .- x_p[i]) ./ dx
    
        #contribute to cells that particles touch 
        w_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) * w_p[i]
        v_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) * w_p[i] * v_p[i]
        
        #count number of particles in the cell 
        N[x_rel .< 1] .+= 1 
    end

    #normalization
    density = w_sum ./ dx 
    velocity = v_sum ./ (dx .* density)

    #do the temperature calculation 
    @inbounds for i in eachindex(Particles.pos)
        #calculate relative position 
        x_rel = abs.(x .- x_p[i]) ./ dx

        T_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) .* w_p[i] .* (v_p[i] .-velocity[x_rel .< 1]).^2
    end
    
    temperature = sqrt.(1 * T_sum ./ ((N.-1)./N.* density .*dx))#hard coded Xe mass

    N_k = 50 
    w_bar = (w_bar * (N_k - 1 ) + (w_sum ./ N)) / N_k  

    return density, velocity, temperature, w_bar, N 
end
=#

end
