using HallPIC: HallPIC as hp
using Documenter
using Test
using CairoMakie: Makie as mk

#Documenter.doctest(hp)

@testset "Harmonic oscillator" begin

	function test_gather(pc::hp.ParticleContainer)
		@inbounds for i in eachindex(pc.pos)
			pc.acc[i] = -pc.pos[i]
		end
	end

	pc = hp.ParticleContainer(Float64, 1, 1, 1)

	# Test 1: Initially-stationary particle remains at rest
	x0 = 0.0
	v0 = 0.0
	dt = 0.2

	pc.pos[1] = 0.0
	pc.vel[1] = 0.0

	# No acceleration at x = 0
	test_gather(pc)
	@test pc.acc[1] == 0

	hp.push_vel(pc, -dt/2)

	# No movement
	hp.push(pc, dt)
	@test pc.pos[1] == 0
	@test pc.vel[1] == 0

	# Test 2: simple harmonic oscillator starting at (x,v) = (1,0)
	# Check that leapfrog is centered as expected
	x0 = 1.0
	v0 = 0
	tmax = 5*pi
	num_steps = ceil(Int, tmax / dt)
	t = range(0.0, tmax, length=num_steps)

	pos = zeros(num_steps)
	vel = zeros(num_steps)

	pc.pos[1] = x0
	pc.vel[1] = v0
	test_gather(pc)

	# half step backward in velocity
	hp.push_vel(pc, -dt/2)
	pos[1] = pc.pos[1]
	vel[1] = pc.vel[1]

	for i in 2:num_steps
		test_gather(pc)
		hp.push(pc, dt)

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
