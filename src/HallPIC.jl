module HallPIC

using DataInterpolations: DataInterpolations, LinearInterpolation
using DocStringExtensions
using Random: Random
using DelimitedFiles 


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
    inds::Vector{Int}
	species::Species
    n_d::Int32
end

function ParticleContainer{T}(N, species::Species) where T<:AbstractFloat
	weight = zeros(T, N)
	pos = zeros(T, N)
	vel = zeros(T, N)
	acc = zeros(T, N)
    inds = zeros(Int, N)
    n_d = 500#hard code for now, should be a simulation hyperparamter
	return ParticleContainer{T}(
		weight, pos, vel, acc, inds, species, n_d 
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
function add_particles!(pc::ParticleContainer{T}, x::Vector{T}, v::Vector{T}, w::Vector{T}, i::Vector{Int}) where T
	# check that new arrays have same length
	M = length(x)
	N = length(pc.pos)
	@assert M == length(v) && M == length(w)
	# append position, velocity, weight, index to pc arrays
	append!(pc.pos, x)
	append!(pc.vel, v)
	append!(pc.weight, w)
    append!(pc.inds, i)
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
Geometry
======================================================#

"""
$(TYPEDSIGNATURES)
Create a particle container and fill it with particles matching the macroscopic properties
from a provided SpeciesProperties object.

Accounting note: the edges of cell i are edges[i] and edges[i+1]
"""
function initialize_particles(sp::SpeciesProperties{T}, grid, particles_per_cell) where T
	pos_buf = zeros(T, particles_per_cell)
	vel_buf = zeros(T, particles_per_cell)
	weight_buf = zeros(T, particles_per_cell)
    ind_buf = ones(Int, particles_per_cell)
	@assert length(grid.cell_centers) == length(sp)

	pc = ParticleContainer{T}(0, sp.species)

    n_itp = LinearInterpolation(sp.dens, grid.cell_centers)
    u_itp = LinearInterpolation(sp.vel, grid.cell_centers)
    T_itp = LinearInterpolation(sp.temp, grid.cell_centers)

    for i in 2:length(grid.cell_centers) - 1
        V = grid.cell_volumes[i]
		z_L = grid.face_centers[i]
		z_R = grid.face_centers[i+1]
		dz = z_R - z_L
 
		Random.rand!(pos_buf)
		@. pos_buf = dz * pos_buf + z_L 

		Random.randn!(vel_buf)
		@. vel_buf *= sqrt(T_itp(pos_buf) / pc.species.gas.mass)
		@. vel_buf = vel_buf + u_itp(pos_buf)

		@. weight_buf = n_itp(pos_buf) / particles_per_cell * V

		add_particles!(pc, pos_buf, vel_buf, weight_buf, ind_buf * i)
	end

	return pc
end

struct OpenBoundary end

struct WallBoundary
    temperature::Float32
end

const HeavySpeciesBoundary = Union{OpenBoundary, WallBoundary}

struct Grid
    cell_centers::Vector{Float64}
    cell_volumes::Vector{Float64}
    face_centers::Vector{Float64}
    face_areas::Vector{Float64}
    left_boundary::HeavySpeciesBoundary
    right_boundary::HeavySpeciesBoundary
end

function Grid(num_cells::Integer, left, right, area, left_boundary = OpenBoundary(), right_boundary = OpenBoundary())
    dz = (right - left) / num_cells

    face_centers = collect(range(left - dz, right + dz, step = dz))
    face_areas = fill(area, length(face_centers))

    cell_centers = zeros(num_cells+2)
    cell_volumes = zeros(num_cells+2)

    for i in eachindex(cell_centers)
        z_L, z_R = face_centers[i], face_centers[i+1]
        cell_centers[i] = 0.5 * (z_L + z_R)
        cell_volumes[i] = area * (z_R - z_L)

    end

    return Grid(cell_centers, cell_volumes, face_centers, face_areas, left_boundary, right_boundary)
end

"""
$(TYPEDSIGNATURES)
"""
function locate_particles!(pc::ParticleContainer{T}, grid::Grid) where T
    for (i, x) in enumerate(pc.pos)
        cell_index = searchsortedfirst(grid.face_centers, x) - 1
        center_pos = grid.cell_centers[cell_index]
        pc.inds[i] = copysign(cell_index, x - center_pos)
    end

    return pc.inds
end

"""
$(TYPEDSIGNATURES)

Deposit particle properties into a gridded SpeciesProperties object
"""
function deposit!(fluid_properties::SpeciesProperties{T}, particles::ParticleContainer, grid::Grid, avg_interval=50) where T

    # zero density, velocity, temperature
    fluid_properties.dens .= zero(T)
    fluid_properties.vel .= zero(T)
    fluid_properties.temp .= zero(T)

    # loop over particles, deposit density and momentum density
    for (ip, z_part) in enumerate(particles.pos)
        # The cell index is positive if the particle is right of the cell center,
        # or negative if it is left of the cell center
        cell_index = particles.inds[ip]
        s, ic = sign(cell_index), abs(cell_index)

        # don't deposit into ghost cells
        if cell_index == -2 || cell_index == length(grid.cell_centers)-1
            s, ic = 0, abs(cell_index)
        else
            s, ic = sign(cell_index), abs(cell_index)
        end

        # Grid information and cell weighting
        inv_vol_c = 1 / grid.cell_volumes[ic]
        inv_vol_s = 1 / grid.cell_volumes[ic+s]
        z_cell = grid.cell_centers[ic]
        dz = grid.face_centers[ic+1] - grid.face_centers[ic]
        t = abs((z_part - z_cell) / dz)

        # Particle velocity and weight
        up = particles.vel[ip]
        wp = particles.weight[ip]

        # Density contribution to each cell
        dens_c = (1 - t) * wp * inv_vol_c
        dens_s = t * wp * inv_vol_s
        fluid_properties.dens[ic]   += dens_c 
        fluid_properties.dens[ic+s] += dens_s

        # Momentum contribution to each cell
        fluid_properties.vel[ic]    += dens_c * up
        fluid_properties.vel[ic+s]  += dens_s * up

        # Use temperature as a buffer to count the number of particles in this cell
        fluid_properties.temp[ic]   += (1 - t)
        fluid_properties.temp[ic+s] += t
    end

    # calculation for average weight

    # Calculate cell-centered fluid properties
    for ic in 2:length(grid.cell_centers)-1 
        density = fluid_properties.dens[ic]
        num_particles = fluid_properties.temp[ic]

        # Get avg. velocity from momentum density
        fluid_properties.vel[ic] /= density

        # Compute smoothed average weight in cell
        new_avg_wt = density / num_particles * grid.cell_volumes[ic]
        old_avg_wt = fluid_properties.avg_weight[ic]
        fluid_properties.avg_weight[ic] = ((avg_interval - 1) * old_avg_wt + new_avg_wt) / avg_interval
    end

    # Zero temperature again
    fluid_properties.temp .= 0.0

    # Calculate energy density
    for (ip, z_part) in enumerate(particles.pos)
        # The cell index is positive if the particle is right of the cell center,
        # or negative if it is left of the cell center
        cell_index = particles.inds[ip]

        # don't deposit into ghost cells
        if cell_index == -2 || cell_index == length(grid.cell_centers)-1
            s, ic = 0, abs(cell_index)
        else
            s, ic = sign(cell_index), abs(cell_index)
        end

        # Grid information and cell weighting
        inv_vol_c = 1.0 / grid.cell_volumes[ic]
        inv_vol_s = 1.0 / grid.cell_volumes[ic+s]
        z_cell = grid.cell_centers[ic]
        dz = grid.face_centers[ic+1] - grid.face_centers[ic]
        t = abs((z_part - z_cell) / dz)

        # Particle velocity and weight
        up = particles.vel[ip]
        wp = particles.weight[ip]

        # Cell average velocity 
        vel_c = fluid_properties.vel[ic]
        vel_s = fluid_properties.vel[ic+s]

        # Accumulate energy density
        dens_c = (1 - t) * wp * inv_vol_c
        dens_s = t * wp * inv_vol_s
        fluid_properties.temp[ic] += dens_c * (up - vel_c)^2
        fluid_properties.temp[ic+s] += dens_s * (up - vel_s)^2
    end

    m = fluid_properties.species.gas.mass

    # Convert energy density to nondimensional temperature
    for ic in 2:length(grid.cell_centers)-1 
        fluid_properties.temp[ic] *= m / fluid_properties.dens[ic]
    end
    
    return fluid_properties
end

#======================================================
Reactions definitions
======================================================#

struct ReactingSpecies
    name::Symbol
    coefficient::Int8 
end
struct Reaction{T<:AbstractFloat}
    reactant::ReactingSpecies
    products::Vector{ReactingSpecies}
    threshold_energy::T
    energies::Vector{T}
    rate::Vector{T}
    delta_n::Vector{T}

end



function reaction_reduction(grid::Grid, reaction::Reaction, electron_properties::SpeciesProperties, fluid_properties::SpeciesProperties, reactant::ParticleContainer,  dt)

    #for each cell, calculate the delta n
    n_cell = length(grid.cells)
    for i=1:n_cell
        #calculate the number of particles produced 
        rate = reaction.rate[electron_properties.temp[i]]#need to fix this for a linear interpolation 
        reaction.delta_n[i] = fluid_properties.dens[i] * electron_properties.dens[i] * rate * dt 
        """
        #add to the daughter products 
        for p=1:n_products
            products[p].delta_n[i] += delta_n * reaction.product_multiplier[p]
        end
        """
    end

    #now that we have the delta_n, adjust particle weights 
    n_particles = length(reactant.pos)
    for i=1:n_particles
        #pull the index 
        ic = reactant.inds[i]
        s = sign(ic)
        ic = abs(ic)
        #count the first cell 
        reactant.weight[i] -= reactant.weight[i] * (reaction.reactant.coefficient * reaction.delta_n[ic]/fluid_properties.dens[ic]) * abs(reactant.pos[i] - grid.cells.center[ic]) / grid.cells.width[ic]
        #count the second cell 
        reactant.weight[i] -= reactant.weight[i] * (reaction.reactant.coefficient * reaction.delta_n[ic+s]/fluid_properties.dens[ic+s]) * abs(reactant.pos[i] - grid.cells.center[ic+s]) / grid.cells.width[ic+s]
    end

    return reaction, reactant 
end

"""
Note that the products vector should be sorted the same way as in the reaction structure 
"""

function daughter_particle_generation(grid::Grid, reaction::Reaction, reactant_properties::SpeciesProperties, product_properties::SpeciesProperties, products::Vector{ParticleContainer})



    n_products = length(products)
    for p = 1:n_products 
        n_d = products[p].n_d#number of desired particles touching cell, hyperparamter from the simualtion
        #for each cell, add particles
        n_cell = length(grid.cells)
        for i=2:n_cell-1
            #determine the number of partices to generate
            w_gen = product_properties.avg_weight[i] * ((sum(products[p].inds==i)+sum(products[p].inds==i-1)+sum(products[p].inds == -(i+1)))/n_d)#check for particles in the cell 
            n_gen = Int32(reaction.products[p].coefficient * reaction.delta_n[i] / w_gen)
            #update the generation weight slightly to consume the full delta_n 
            w_gen = reaction.delta_n[i] / n_gen 

            #generate the particles 
            #position 
            z_L = grid.cells[i].left_face.pos                
            pos_buf = grid.cells[i].width * Random.rand(n_gen) + z_L 

            #velocity 
            v_th = sqrt(reactant_properties.temp[i] / product.species.gas.mass)
            vel_buf = Random.randn(vel_buf) * v_th + reactant_properties.vel[i]

            #weight 
            weight_buf = w_gen.*ones(n_gen)
                
            #actual generation 
            add_particles!(product, pos_buf, vel_buf, weight_buf, ones(n_gen).*i)
        end
    end 

    return products
end 

function read_reaction_rates(filepath)

    energy, values = open(filepath) do file 
        energy = parse(Float64, strip(split(readline(file), ':')[2]))
        values = readdlm(file, skipstart = 1)
        energy, values 
    end

    return energy, values[:,1], values[:,2]
end

end
