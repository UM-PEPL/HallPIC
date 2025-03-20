var documenterSearchIndex = {"docs":
[{"location":"#Design-Document","page":"Design Document","title":"Design Document","text":"","category":"section"},{"location":"","page":"Design Document","title":"Design Document","text":"Here are some unordered thoughts for the design of this PIC module for HallThruster.jl","category":"page"},{"location":"#Modular-design","page":"Design Document","title":"Modular design","text":"","category":"section"},{"location":"","page":"Design Document","title":"Design Document","text":"We should strive to make the code as modular as possible, with only the minimum possible information shared between modules. Preliminarily, the modules that jump out to me are ","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"Heavy species\nElectron energy\nOhm's law","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"The electron energy and Ohm's law modules should know nothing about the implementation of the heavy species module (i.e. whether it is fluid or PIC). The same should go in reverse. The heavy species module should provide the following information to the other parts of the code:","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"Species densities/velocities/temperatures\nIonization/excitation energy losses","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"If these are sufficiently distinct, the heavy species module could be split further into one which advances particles/handles BCs, and one which performs reaction calculations.","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"We should take care to define the public APIs of each modules in such a way that implementation details are hidden as much as possible. This was a major mistake in HallThruster.jl, which exposed the internal details of the code to the user in a way that made refactoring without changing the user interface very challenging.","category":"page"},{"location":"#Temporal-order-of-accuracy-and-operator-splitting","page":"Design Document","title":"Temporal order of accuracy and operator splitting","text":"","category":"section"},{"location":"","page":"Design Document","title":"Design Document","text":"HallThruster.jl does not at present use a principled operator splitting scheme for updating its various subcomponents. This caps us to 1st order temporal accuracy and may be less than ideal from a stability standpoint. Using the highly modular design approach outlined above, we can instead adopt a more rigorous splitting strategy (likely Strang splitting) to perform the update at each timestep. This lets each module take care of its own timestepping in the way that makes then most sense for the numerics of that module, and gives us the possibility of second-order temporal accuracy. This will prove valuable for cases in which the transient behavior of the thruster is important.","category":"page"},{"location":"#Strong-emphasis-on-verification-and-testing","page":"Design Document","title":"Strong emphasis on verification and testing","text":"","category":"section"},{"location":"","page":"Design Document","title":"Design Document","text":"Types of test","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"Unit tests: These are tests of individual functions or \"units\" in the code.   They are the easiest kind of test to write, and are essential for making sure the nuts and bolts of the code work as anticipated.   As they grow to test larger chunks of the code, however, they can become brittle and flaky as more and more test setup is needed to \"mock\" parts of the code not directly under testing.\nProperty tests: These check to make sure specific properties hold.   When done well, these can automatically search an input space to try and find cases that violate the properties we have set.   Whenenver possible, we should try and identify properties that either 1) always hold or 2) hold under some known set of conditions and test for these automatically.   Some easy ones include mass conservation and one of momentum/energy conservation, depending on the PIC scheme we pic.   Should also check current conservation.   Smaller invariants are also worth checking, like verifying that the output of some function is always of a specific form or in a specific range\nOrder verification tests: These use artificial source terms and the method of manufactured solutions to check that a given discretized PDE is being solved to the desired order of accuracy.   We use these in HallThruster.jl for the ion fluid model and the electron energy equation, but we could also do this for other parts of the code, like gradient computations.   We should try to implement these in a way that requires as little mocking as possible, as the setup code for these tests needs to maintain and can become burdensome.\nRegression tests: These check that the observable behavior of the code has not changed unexpectedly.   There are three main types of regression test, each with their own challenges and uses.\nBehavior regression: This checks that the code has not changed its overall behavior.  For us, this means the predicted thruster performance and plasma properties.  These tests are easy to implement, but rely on you knowing exactly what observable behavior you care about, and what exactly constitutes a change.\nPerformance regression: This type of test checks that we haven't made the code slower.  This is quite challenging, as the performance of the code varies from machine to machine so we cannot easily specify a standard of performance that we can check against.  I have not yet figured out a good way to do these sorts of tests.\nBug regression: This sort of test verifies that a bug stays fixed.  When we fix a bug, a test that would have triggered the bug should be added to the code, so that we can check for all time afterward that the bug stays squashed","category":"page"},{"location":"#Axisymmetry-and-quasi-1D-flow","page":"Design Document","title":"Axisymmetry and quasi-1D flow","text":"","category":"section"},{"location":"","page":"Design Document","title":"Design Document","text":"It would be nice if particles could move in 2D so that we can more self-consistently model wall losses. This could move the code in a quasi-2D direction, given some assumptions. Axisymmetry complicates the PIC update significantly, however, so this needs to be done carefully.","category":"page"},{"location":"","page":"Design Document","title":"Design Document","text":"Modules = [HallPIC]","category":"page"},{"location":"#HallPIC.CellProperties","page":"Design Document","title":"HallPIC.CellProperties","text":"struct CellProperties{T<:AbstractFloat}\n\nContains fluid information for a single cell.\n\ndens::AbstractFloat\nvel::AbstractFloat\ntemp::AbstractFloat\navg_weight::AbstractFloat\n\n\n\n\n\n","category":"type"},{"location":"#HallPIC.ParticleContainer","page":"Design Document","title":"HallPIC.ParticleContainer","text":"struct ParticleContainer{T<:AbstractFloat}\n\nPer-particle information for a single species. T is expected to be either Float32 or Float64\n\nweight::Vector{T} where T<:AbstractFloat\npos::Vector{T} where T<:AbstractFloat\nvel::Vector{T} where T<:AbstractFloat\nacc::Vector{T} where T<:AbstractFloat\ninds::Vector{Int64}\nspecies::HallPIC.Species\n\n\n\n\n\n","category":"type"},{"location":"#HallPIC.Species","page":"Design Document","title":"HallPIC.Species","text":"struct Species\n\nA single chemical species (charge + mass)\n\ngas::HallPIC.Gas\nreactions::Vector{UInt8}\ncharge::Int8\n\n\n\n\n\n","category":"type"},{"location":"#HallPIC.SpeciesProperties","page":"Design Document","title":"HallPIC.SpeciesProperties","text":"struct SpeciesProperties{T<:AbstractFloat}\n\nT: Float32 or Float64 Holds fluid properties (density, velocity, temperature) for a single species.\n\ndens::Vector{T} where T<:AbstractFloat\nvel::Vector{T} where T<:AbstractFloat\ntemp::Vector{T} where T<:AbstractFloat\navg_weight::Vector{T} where T<:AbstractFloat\nspecies::HallPIC.Species\n\n\n\n\n\n","category":"type"},{"location":"#HallPIC.add_particles!-Union{Tuple{T}, Tuple{HallPIC.ParticleContainer{T}, Vector{T}, Vector{T}, Vector{T}}} where T","page":"Design Document","title":"HallPIC.add_particles!","text":"add_particles!(\n    pc::HallPIC.ParticleContainer{T},\n    x::Array{T, 1},\n    v::Array{T, 1},\n    w::Array{T, 1}\n) -> HallPIC.ParticleContainer{T} where T\n\n\nAdd particles to a ParticleContainer.\n\n\n\n\n\n","category":"method"},{"location":"#HallPIC.deposit!-Union{Tuple{T}, Tuple{HallPIC.SpeciesProperties{T}, HallPIC.ParticleContainer, HallPIC.Grid}, Tuple{HallPIC.SpeciesProperties{T}, HallPIC.ParticleContainer, HallPIC.Grid, Any}} where T","page":"Design Document","title":"HallPIC.deposit!","text":"deposit!(\n    fluid_properties::HallPIC.SpeciesProperties{T},\n    particles::HallPIC.ParticleContainer,\n    grid::HallPIC.Grid\n) -> HallPIC.SpeciesProperties{T} where T\ndeposit!(\n    fluid_properties::HallPIC.SpeciesProperties{T},\n    particles::HallPIC.ParticleContainer,\n    grid::HallPIC.Grid,\n    avg_interval\n) -> HallPIC.SpeciesProperties{T} where T\n\n\nDeposit particle properties into a gridded SpeciesProperties object\n\n\n\n\n\n","category":"method"},{"location":"#HallPIC.initialize_particles-Union{Tuple{T}, Tuple{HallPIC.SpeciesProperties{T}, Any, Any}} where T","page":"Design Document","title":"HallPIC.initialize_particles","text":"initialize_particles(\n    sp::HallPIC.SpeciesProperties{T},\n    grid,\n    particles_per_cell\n) -> HallPIC.ParticleContainer\n\n\nCreate a particle container and fill it with particles matching the macroscopic properties from a provided SpeciesProperties object.\n\nAccounting note: the edges of cell i are edges[i] and edges[i+1]\n\n\n\n\n\n","category":"method"},{"location":"#HallPIC.locate_particles!-Union{Tuple{T}, Tuple{HallPIC.ParticleContainer{T}, HallPIC.Grid}} where T","page":"Design Document","title":"HallPIC.locate_particles!","text":"locate_particles!(\n    pc::HallPIC.ParticleContainer{T},\n    grid::HallPIC.Grid\n) -> Vector{Int64}\n\n\n\n\n\n\n","category":"method"},{"location":"#HallPIC.push!-Tuple{HallPIC.ParticleContainer, AbstractFloat}","page":"Design Document","title":"HallPIC.push!","text":"push!(\n    pc::HallPIC.ParticleContainer,\n    dt::AbstractFloat\n) -> HallPIC.ParticleContainer\n\n\nPush particle container to next timestep using Leapfrog scheme\n\n\n\n\n\n","category":"method"}]
}
