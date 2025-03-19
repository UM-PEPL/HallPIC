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

function add_particle(pc::ParticleContainer, x, v, a, w)
    """
    Add a particle to a container 

    Inputs: 
    pc (Particle Container Object)
        the particles object to add the particle too
    x (float)
        the position of the particle 
    v (float)
        the velocity of the particle 
    a (float)
        the acceleration of the particle  
    w (float)
        the weight of the particle
    Outputs:
    pc (Particle Container Object)
        modified particles object 
    """

    #add to the particle 
    append!(pc.pos, x)
    append!(pc.vel, v)
    append!(pc.acc, a)
    append!(pc.weight, w)


    return pc 
end

struct SpeciesGridProperties{T<:AbstractFloat, I<:Integer}
    """
    T: Float32 or Float64
    I: UInt 
    """
    n::Matrix{T}
    v::Matrix{T}
    Temp::Matrix{T}
    w_bar::Matrix{T}
    w_tot::Matrix{T}
    N::Matrix{I}
end

function SpeciesGridProperties(::Type{T}, N_cell, N_species) where T
	n = zeros(T, N_cell, N_species)
	v = zeros(T, N_cell, N_species)
	Temp = zeros(T, N_cell, N_species)
	w_bar = zeros(T, N_cell, N_species)
	w_tot = zeros(T, N_cell, N_species)
    N = zeros(Int32, N_cell, N_species)
	return SpeciesGridProperties{T, Int32}(
		n, v, Temp, w_bar, w_tot, N)
end


function initialize_particles(::Type{T}, n_cell, x, dx, n_ppc, density, velocity, temperature, mass, charge) where T
    """
    Initialize a set of particles using a grid and properties 
    Inputs: 
        N_cell (int)
            number of cells on the grid
        x (N_cell+1 array of floats)
            positions of cell edges 
        dx (N_cell array of floats)
            cell widths
        N_per (int)
            number of particles per cell 
        density (N_cell array of floats)
            initial number density for each cell
        velocity (N_cell array of floats)
            initial velocity for each cell
        temperature (N_cell array of floats)
            initial temperature for each cell
        mass (float)
            mass for an individual particle for the set 
        Charge_Index(Int)
            Index to the state of the charge 
    Outputs:
        Particles (HallPIC particle Object)
            object containing particle information
        w_bar (N_cell array of floats)
            time-averaged average particle weight
    """

    #initialize particle object 
    particles = ParticleContainer(T, 0, mass, charge)
    #create storage arrays 
    n_p = n_ppc*n_cell
    cell_idx = 1
    for i = 1:n_p
        #sample uniformly in the cell 
        append!(particles.pos, rand()*dx[cell_idx] + x[cell_idx])
        #sample from Maxwellian and transform
        append!(particles.vel, randn() * temperature[cell_idx] + velocity[cell_idx])
        #evenly distribute the weights
        append!(particles.weight, density[cell_idx] .* dx[cell_idx] ./ n_ppc)
        #update the current cell 
        cell_idx = Int(floor(i/n_ppc) + 1)
    end

    #set average weight
    w_bar = density / n_ppc 

    return particles, w_bar 
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

end
