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

Base.isdone(sp::SpeciesProperties, i) = i > length(sp)
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
	pos_buf = zeros(T, particles_per_cell)
	vel_buf = zeros(T, particles_per_cell)
	weight_buf = zeros(T, particles_per_cell)

	num_cells = length(edges) - 1

	@assert num_cells == length(sp)

	pc = ParticleContainer{T}(0, sp.species)

	for (i, (V, cell)) in enumerate(zip(volumes, sp))
		z_L = edges[i]
		z_R = edges[i+1]
		dz = z_R - z_L

		Random.rand!(pos_buf)
		@. pos_buf = dz * pos_buf + z_L 

		Random.randn!(vel_buf)
		v_th = sqrt(cell.temp / pc.species.gas.mass)
		@. vel_buf = vel_buf * v_th + cell.vel

		weight = cell.dens * V / particles_per_cell
		@. weight_buf = weight

		add_particles!(pc, pos_buf, vel_buf, weight_buf)
	end

	return pc
end


function deposit(cell_centers, volumes, fluid_properties::SpeciesProperties{T}, particles::ParticleContainer) where T
    """
    Initialize a set of particles using a grid and properties 
    Inputs: 
        cell_centers (N_cell+1 array of floats)
            positions of cell centers 
        volumes (N_cell array of floats)
            cell widths/volume
        Fluid_Properties (Species properties object)
            object containing fluid properties for the species on the grid, to be updated
        Particles (Particle Object)
            object containing particle information 
    Outputs:
        Fluid_Properties (Species properties object)
            updated fluid properties for the species on the grid
    """

    #initialize sums 
    n_cell = length(cell_centers)
    w_sum = zeros(T, n_cell)
    v_sum = zeros(T, n_cell)
    temp_sum = zeros(T, n_cell)    
    n = zeros(T, n_cell)

    #pull quantities
    n_p = size(particles.pos)[1]
    #loop over particles
    #technically there's a loop over cells in here too, but that tis handled with logical indexing
    for i=1:n_p
        #calculate relative position 
        x_rel = abs.(cell_centers .- particles.pos[i]) ./ volumes
    
        #contribute to cells that particles touch 
        w_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) * particles.weight[i]
        v_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) * particles.weight[i] * particles.vel[i]
        
        #count number of particles in the cell 
        n[x_rel .< 1] .+= 1 
    end


    #normalization
    #assign the properties 
    fluid_properties.dens .= w_sum ./ volumes 
    fluid_properties.vel .= v_sum ./ (volumes .* fluid_properties.dens)

    #do the temperature calculation 
    for i=1:n_p 
        #calculate relative position 
        x_rel = abs.(cell_centers .- particles.pos[i]) ./ volumes
        #particle contribution to temperature 
        temp_sum[x_rel .< 1] += (1 .- x_rel[x_rel .< 1]) .* particles.weight[i] .* (particles.vel[i] .-fluid_properties.vel[x_rel .< 1]).^2
    end
    #normalize
    fluid_properties.temp .= sqrt.(1 * temp_sum ./ ((n.-1)./n.* fluid_properties.dens .*volumes))
    #calculation for average weight
    #hold time average interval to 50 for now (see Dominguez Vazquez thesis) can have as an input parameter later 
    n_k = 50 
    fluid_properties.avg_weight .= (fluid_properties.avg_weight * (n_k - 1 ) + (w_sum ./ n)) / n_k  

    
    return fluid_properties
end


end
