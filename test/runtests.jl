using HallPIC: HallPIC as hp
using Documenter
using Test

#Documenter.doctest(hp)

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

	pc = hp.ParticleContainer(T, 1, 1, 1)

	# Test 1: Initially-stationary particle remains at rest
	pc.pos[1] = 0.0
	pc.vel[1] = 0.0

	# No acceleration at x = 0
	leapfrog_gather(pc)
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
	leapfrog_gather(pc)

	# half step backward in velocity
	hp.push_vel!(pc, -dt/2)
	pos[1] = pc.pos[1]
	vel[1] = pc.vel[1]

	for i in 2:num_steps
		leapfrog_gather(pc)
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
	pc = hp.ParticleContainer(T, 1, 1, 1)
	M = 100
	L = 0.05
	x = LinRange(0, L, M)
	E = V / L * ones(M)

	x_norm = x ./ hp.x_0
	E_norm = E ./ hp.E_0

	extrap = hp.DataInterpolations.ExtrapolationType.Constant
	E_itp = hp.LinearInterpolation(E_norm, x_norm, extrapolation = extrap)

	q = pc.charge * hp.q_e
	m = pc.mass * hp.m_0

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

	a = pc.charge * E_norm[1] / pc.mass
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

function test_initialization(::Type{T}) where T 
	#try initalizing the particles object 

	particles = hp.ParticleContainer(T, 0, 131.29*1.66e-27, 1)

	#check the initialization 
	@test particles.charge == 1
	@test particles.mass ≈ 131.29*1.66e-27
	@test isempty(particles.pos)
	@test isempty(particles.vel)
	@test isempty(particles.acc)
	@test isempty(particles.weight)

	#try adding two particles 
	particles = hp.add_particle(particles, 1.57, 19823.1041, 0.0, 1.589e16)
	particles = hp.add_particle(particles, 0.51, 981.471, 5.402e9, 1.589e16)

	#check them 
	@test particles.pos[1] ≈ 1.57
	@test particles.pos[2] ≈ 0.51 
	@test particles.vel[1] ≈ 19823.1041
	@test particles.vel[2] ≈ 981.471
	@test particles.acc[1] ≈ 0.0
	@test particles.acc[2] ≈ 5.402e9
	@test particles.weight[1] ≈ 1.589e16
	@test particles.weight[2] ≈ 1.589e16
end

@testset "Initialization" begin
	for T in [Float32, Float64]
		test_initialization(T)
	end
end

function test_Initial_Deposition(::Type{T}) where T
	#define the cell properties
	N_cell = 4
	N_n = 50000#per cell 


	#seed properties/initialize cell arrays 
	nn = ones(N_cell) * 1e18#1/m^3
	vn = ones(N_cell) * 300
	Tn = ones(N_cell) * (500 / 11604) #eV
	x = 0:(1/N_cell):1
	x_c = (x[2:end] + x[1:end-1]) / 2
	dx = x[2:end] - x[1:end-1]

	Neutrals, w_bar_n = hp.initialize_particles(T, N_cell, x, dx, N_n, nn, vn, Tn, 131.29*1.66e-27, 0)

	#check that particles are in bounds and limits are correct
	@test size(Neutrals.pos)[1] == N_cell * N_n
	min, max = extrema(Neutrals.pos)
	@test max <= 1
	@test min >= 0

	min, max = extrema(Neutrals.weight)
	@test max ≈ 5e12
	@test min ≈ 5e12

	min, max = extrema(Neutrals.vel)
	thermal_speed = sqrt(2 * hp.q_e * (500 / 11604) / Neutrals.mass)
	@test max <= 300 + 10 * thermal_speed
	@test min >= 300 - 10 * thermal_speed

	

	nn, vn, Tn, w_bar_n, N_n_cell = hp.Deposit(N_cell, x_c, dx, w_bar_n, Neutrals)

	#check that interior cells reproduce
	@test isapprox(nn[2], 1e18; rtol = 0.01)
	@test isapprox(vn[2], 300; rtol = 0.01)
	@test isapprox(Tn[2], (500 / 11604); rtol = 0.01)

	@test isapprox(nn[3], 1e18; rtol = 0.01)
	@test isapprox(vn[3], 300; rtol = 0.01)
	@test isapprox(Tn[3], (500 / 11604); rtol = 0.01)

end

@testset "Initial Deposition" begin
	for T in [Float32, Float64]
		test_Initial_Deposition(T)
	end
end