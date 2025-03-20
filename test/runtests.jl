using HallPIC: HallPIC as hp
using Documenter
using Test

Documenter.doctest(hp)

const Hydrogen = hp.Gas(name=:H, mass=1)
const Xenon = hp.Gas(name=:Xe, mass=131.293)

function test_leapfrog(::Type{T}) where T
	function leapfrog_gather!(pc::hp.ParticleContainer)
		@inbounds for i in eachindex(pc.pos)
			pc.acc[i] = -pc.pos[i]
		end
	end

	dt = 0.2
	x0 = 1.0
	v0 = 0
	tmax = 5*pi
	num_steps = ceil(Int, tmax / dt)
	t = range(0.0, tmax, length=num_steps)

	pos = zeros(num_steps)
	vel = zeros(num_steps)

	pc = hp.ParticleContainer{T}(1, Hydrogen(1))

	# Test 1: Initially-stationary particle remains at rest
	pc.pos[1] = 0.0
	pc.vel[1] = 0.0

	# No acceleration at x = 0
	leapfrog_gather!(pc)
	@test pc.acc[1] == 0

	hp.push_vel!(pc, -dt/2)

	# No movement
	hp.push!(pc, dt)
	@test pc.pos[1] == 0
	@test pc.vel[1] == 0

	# Test 2: simple harmonic oscillator starting at (x,v) = (1,0)
	# Check that leapfrog is centered as expected
	pc.pos[1] = x0
	pc.vel[1] = v0
	leapfrog_gather!(pc)

	# half step backward in velocity
	hp.push_vel!(pc, -dt/2)
	pos[1] = pc.pos[1]
	vel[1] = pc.vel[1]

	for i in 2:num_steps
		leapfrog_gather!(pc)
		hp.push!(pc, dt)

		pos[i] = pc.pos[1]
		vel[i] = pc.vel[1]
	end

	# Check that leapfrog is centered
	min, max = extrema(pos)
	@test min >= -1
	@test max <= 1

	min, max = extrema(vel)
	@test min >= -1
	@test max <= 1
end

@testset "Harmonic oscillator" begin
	for T in [Float32, Float64]
		test_leapfrog(T)
	end
end

function test_linear_drop(::Type{T}, V) where T
	pc = hp.ParticleContainer{T}(1, Hydrogen(1))
	M = 100
	L = 0.05
	x = LinRange(0, L, M)
	E = V / L * ones(M)

	x_norm = x ./ hp.x_0
	E_norm = E ./ hp.E_0

	extrap = hp.DataInterpolations.ExtrapolationType.Constant
	E_itp = hp.LinearInterpolation(E_norm, x_norm, extrapolation = extrap)

	q = hp.charge(pc.species) * hp.q_e
	m = hp.mass(pc.species) * hp.m_0

	u_expected = sqrt(2 * q * V / m)
	dt = 0.5 * L / u_expected / M
	dt_norm = dt / hp.t_0

	# pushback step
	hp.gather!(pc, E_itp)
	hp.push_vel!(pc, -dt_norm/2)

	t = 0.0
	while pc.pos[1] < x_norm[end] 
		t += dt_norm
		hp.gather!(pc, E_itp)
		hp.push!(pc, dt_norm)
	end

	hp.push_vel!(pc, dt_norm/2)

	a = hp.charge(pc.species) * E_norm[1] / hp.mass(pc.species)
	u_exact = T(a * t)
	x_exact = T(0.5 * a * t^2)

	@test x_exact ≈ pc.pos[1]
	@test u_exact ≈ pc.vel[1]
end

@testset "Linear potential drop" begin
	for T in [Float32, Float64]
		test_linear_drop(T, 600.0)
	end
end


function test_add_particles(::Type{T}) where T 
	# try initalizing the particles object 
	particles = hp.ParticleContainer{T}(0, Xenon(1))

	#check the initialization 
	@test hp.charge(particles.species) == 1
	@test hp.mass(particles.species) == Xenon.mass
	@test isempty(particles.pos)
	@test isempty(particles.vel)
	@test isempty(particles.acc)
	@test isempty(particles.weight)

	# try adding two particles 
	pos = [T(1.57), T(0.51)]
	vel = [T(19823.1041), T(981.471)]
	weight = [T(1.589e16), T(1.589e16)]

	particles = hp.add_particles!(particles, pos, vel, weight)

	#check them 
	@test length(particles.pos) == 2
	@test particles.pos[1] == pos[1]
	@test particles.pos[2] == pos[2]

	@test length(particles.vel) == 2
	@test particles.vel[1] == vel[1]
	@test particles.vel[2] == vel[2]

	@test length(particles.weight) == 2
	@test particles.weight[1] == weight[1]
	@test particles.weight[2] == weight[2]

	@test length(particles.acc) == 2
	@test particles.acc[1] == 0
	@test particles.acc[2] == 0
end

@testset "Add particles" begin
	for T in [Float32, Float64]
		test_add_particles(T)
	end
end

@testset "SpeciesProperties" begin
	for T in [Float32, Float64]
		N = 5
		sp = hp.SpeciesProperties{T}(N, Hydrogen(1))
		@test length(sp) == N

		base = collect(1.0:N)
		@. sp.dens = 1.0 * base
		@. sp.vel = 2.0 * base
		@. sp.temp = 3.0 * base
		@. sp.avg_weight = 4.0*base

		@test eltype(sp) == hp.CellProperties{T}

		for (i, cell) in enumerate(sp)
			@test cell.dens == sp.dens[i] && cell.dens == sp[i].dens
			@test cell.vel == sp.vel[i] && cell.vel == sp[i].vel
			@test cell.temp == sp.temp[i] && cell.temp == sp[i].temp
			@test cell.avg_weight == sp.avg_weight[i] == sp[i].avg_weight
			@test cell == sp[i]
		end
		c = collect(sp)
		@test eltype(c) == hp.CellProperties{T}
		@test length(c) == N

		sp_empty = hp.SpeciesProperties{T}(0, Hydrogen(1))
		for _ in sp_empty
		end
		c = collect(sp_empty)
		@test eltype(c) == hp.CellProperties{T}
		@test length(c) == 0
		@test length(sp_empty) == 0
	end
end


function test_initial_deposition(_::Type{T}) where T

	#define the cell properties
	n_cell = 4
	n_n = 50000#per cell 


	#seed properties/initialize cell arrays 
	neutral_properties = hp.SpeciesProperties{Float64}(n_cell, Xenon(0))
	base = ones(n_cell)
	neutral_properties.dens .= base .* 1e18#1/m^3
	neutral_properties.vel .= base .* 300
	neutral_properties.temp .= base .* (500 / 11604) #eV
	x = 0:(1/n_cell):1
	x_c = (x[2:end] + x[1:end-1]) / 2
	dx = x[2:end] - x[1:end-1]

	neutrals = hp.initialize_particles(neutral_properties, x, dx, n_n)

	#check that particles are in bounds and limits are correct
	@test size(neutrals.pos)[1] == n_cell * n_n
	min, max = extrema(neutrals.pos)
	@test max <= 1
	@test min >= 0

	min, max = extrema(neutrals.weight)
	@test max ≈ 5e12
	@test min ≈ 5e12

	min, max = extrema(neutrals.vel)
	thermal_speed = sqrt(2 * hp.q_e * (500 / 11604) / (neutrals.species.gas.mass * hp.m_0))
	@test max <= 300 + 10 * thermal_speed
	@test min >= 300 - 10 * thermal_speed

	
	#deposit to the cells 
	new_neutral_properties = hp.deposit(x_c, dx, neutral_properties, neutrals)

	#check that interior cells reproduce
	@test isapprox(new_neutral_properties.dens[2], neutral_properties.dens[2]; rtol = 0.01)
	@test isapprox(new_neutral_properties.vel[2], neutral_properties.vel[2]; rtol = 0.01)
	@test isapprox(new_neutral_properties.temp[2], neutral_properties.temp[2]; rtol = 0.01)

	@test isapprox(new_neutral_properties.dens[3], neutral_properties.dens[3]; rtol = 0.01)
	@test isapprox(new_neutral_properties.vel[3], neutral_properties.vel[3]; rtol = 0.01)
	@test isapprox(new_neutral_properties.temp[3], neutral_properties.temp[3]; rtol = 0.01)

end

#=
@testset "Initial Deposition" begin
	for T in [Float32, Float64]
		test_initial_deposition(T)
	end
end
=#

@testset "Grid construction" begin
	N = 100
	left_boundary = hp.OpenBoundary()
	right_boundary = hp.WallBoundary(0.5)
	area = 1.0
	x0 = 0
	x1 = 1
	grid = hp.Grid(N, x0, x1, area, left_boundary, right_boundary)

	@test length(grid.cell_centers) == N+2
	@test length(grid.face_centers) == N+3
	@show grid.face_centers
	@test grid.left_boundary == left_boundary
	@test grid.right_boundary == right_boundary
	dz = (x1 - x0) / N
	@test all(volume ≈ dz * area for volume in grid.cell_volumes)
	@test grid.face_centers[end-1] - grid.face_centers[2] == (x1 - x0)
end

@testset "Particle location" begin
	N = 100
	x0 = 0.
	x1 = 1.0*N
	grid = hp.Grid(N, x0, x1, 1.0)
	dx = (x1 - x0) / N
	pc = hp.ParticleContainer{Float32}(7, Hydrogen(1))

	pc.pos[1] = 50.5
	pc.pos[2] = 50.1
	pc.pos[3] = 50.7
	pc.pos[4] = x0
	pc.pos[5] = -0.6
	pc.pos[6] = dx
	pc.pos[7] = x1

	hp.locate_particles!(pc, grid)

	@test pc.inds[1] == 52
	@test pc.inds[2] == -52
	@test pc.inds[3] == 52
	@test pc.inds[4] == 1
	@test pc.inds[5] == -1
	@test pc.inds[6] == 2
	@test pc.inds[7] == N+1
end
function test_read_reaction_table()
	#load the reaction 
	filepath = "../reactions/ionization_Xe_Xe+.dat"
	threshold_energy, energy, rates = hp.read_reaction_rates(filepath)

	#ensure that the values are correct 
	@test threshold_energy ≈ 12.1298437

	@test energy[1] ≈ 0.3878E-01
	@test energy[100] ≈ 499.5

	@test rates[1] ≈ 0.000
	@test rates[100] ≈ 0.3526E-12

end


function test_initialize_reaction(::Type{T}) where T

	#load the rate table 
	filepath = "../reactions/ionization_Xe_Xe+.dat"
	threshold_energy, energy, rates = hp.read_reaction_rates(filepath)

	#define the species  
	reactant = hp.ReactingSpecies(Xenon(0).gas.name, 1)
	product = [hp.ReactingSpecies(Xenon(1).gas.name, 1)]

	#initialize the struct 
	Xe_ionization = hp.Reaction{T}(reactant, product, threshold_energy, energy, rates, [0.0, 0.0, 0.0, 0.0])

	#actually test 
	@test Xe_ionization.reactant == reactant
	@test Xe_ionization.products == product 
	@test Xe_ionization.threshold_energy ≈ threshold_energy
	@test Xe_ionization.energies ≈ energy 
	@test Xe_ionization.rate ≈ rates 
	@test Xe_ionization.delta_n ≈ zeros(4) 


end


function test_reaction_step(::Type{T}) where T

	#load the rate table 
	filepath = "../reactions/ionization_Xe_Xe+.dat"
	threshold_energy, energy, rates = hp.read_reaction_rates(filepath)

	#define the species  
	reactant = hp.ReactingSpecies(Xenon(0).gas.name, 1)
	product = [hp.ReactingSpecies(Xenon(1).gas.name, 1)]

	#initialize the reaction struct 
	Xe_ionization = hp.Reaction{T}(reactant, product, threshold_energy, energy, rates, [0.0, 0.0, 0.0, 0.0])

	#seed properties/initialize cell arrays 
	n_cell = 2 
	neutral_properties = hp.SpeciesProperties{Float64}(n_cell+2, Xenon(0))
	ion_properties = hp.SpeciesProperties{Float64}(n_cell+2, Xenon(1))
	base = ones(n_cell+2)
	neutral_properties.dens .= base .* 1e18#1/m^3
	neutral_properties.vel .= base .* 300
	neutral_properties.temp .= base .* (500 / 11604) #eV
	ion_properties.dens .= base .* 1e16#1/m^3
	ion_properties.vel .= base .* 5000
	ion_properties.temp .= base .* 0.1 #eV
	grid = hp.Grid(2, 0, 1, 3.14, hp.OpenBoundary(), hp.OpenBoundary())
	x = 0:(1/n_cell):1
	x_c = (x[2:end] + x[1:end-1]) / 2
	dx = x[2:end] - x[1:end-1]

	
	neutrals = hp.initialize_particles(neutral_properties, grid.cells[:].left_face, grid.cells.width, n_n)
	ions = hp.initialize_particles(ion_properties, grid.cells.left_face, grid.cells.width, n_n)

	#deposit to grid 
	neutral_properties = hp.deposit(grid.cells.center, grid.cells.width, neutral_properties, neutrals)
	ion_properties = hp.deposit(grid.cells.center, grid.cells.width, ion_properties, ions)

	#initialization for reaction properties 


	#reduce weights 

	#add particles 
end



@testset "Reactions" begin
	test_read_reaction_table()
	for T in [Float32, Float64]
		test_initialize_reaction(T)
		test_reaction_step(T)
	end
end