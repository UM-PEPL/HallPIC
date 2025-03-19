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


function test_Initial_Deposition(::Type{T}) where T
	#initialize species 
	Xenon = hp.Gas(name=:Xe, mass=131.293)

	#define the cell properties
	N_cell = 4
	N_n = 50000#per cell 


	#seed properties/initialize cell arrays 
	Neutral_Properties = hp.SpeciesProperties{Float64}(N_cell, Xenon(0))
	base = ones(N_cell)
	Neutral_Properties.dens .= base .* 1e18#1/m^3
	Neutral_Properties.vel .= base .* 300
	Neutral_Properties.temp .= base .* (500 / 11604) #eV
	x = 0:(1/N_cell):1
	x_c = (x[2:end] + x[1:end-1]) / 2
	dx = x[2:end] - x[1:end-1]

	Neutrals = hp.initialize_particles(Neutral_Properties, x, dx, N_n)

	#check that particles are in bounds and limits are correct
	@test size(Neutrals.pos)[1] == N_cell * N_n
	min, max = extrema(Neutrals.pos)
	@test max <= 1
	@test min >= 0

	min, max = extrema(Neutrals.weight)
	@test max ≈ 5e12
	@test min ≈ 5e12

	min, max = extrema(Neutrals.vel)
	thermal_speed = sqrt(2 * hp.q_e * (500 / 11604) / (Neutrals.species.gas.mass * hp.m_0))
	@test max <= 300 + 10 * thermal_speed
	@test min >= 300 - 10 * thermal_speed

	
	#deposit to the cells 
	New_Neutral_Properties = hp.deposit(x_c, dx, Neutral_Properties, Neutrals)

	#check that interior cells reproduce
	@test isapprox(New_Neutral_Properties.dens[2], Neutral_Properties.dens[2]; rtol = 0.01)
	@test isapprox(New_Neutral_Properties.vel[2], Neutral_Properties.vel[2]; rtol = 0.01)
	@test isapprox(New_Neutral_Properties.temp[2], Neutral_Properties.temp[2]; rtol = 0.01)

	@test isapprox(New_Neutral_Properties.dens[3], Neutral_Properties.dens[3]; rtol = 0.01)
	@test isapprox(New_Neutral_Properties.vel[3], Neutral_Properties.vel[3]; rtol = 0.01)
	@test isapprox(New_Neutral_Properties.temp[3], Neutral_Properties.temp[3]; rtol = 0.01)

end

@testset "Initial Deposition" begin
	for T in [Float32, Float64]
		test_Initial_Deposition(T)
	end
end
