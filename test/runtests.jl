using HallPIC: HallPIC as hp
using Documenter
using Test

Documenter.doctest(HallPIC)

@testset "HallPIC.jl" begin
	@test hp.add1(1) == 2
end
