# Design Document

Here are some unordered thoughts for the design of this PIC module for HallThruster.jl

## Modular design

We should strive to make the code as modular as possible, with only the minimum possible information shared between modules.
Preliminarily, the modules that jump out to me are 
- Heavy species
- Electron energy
- Ohm's law

The electron energy and Ohm's law modules should know nothing about the implementation of the heavy species module (i.e. whether it is fluid or PIC).
The same should go in reverse.
The heavy species module should provide the following information to the other parts of the code:
- Species densities/velocities/temperatures
- Ionization/excitation energy losses

If these are sufficiently distinct, the heavy species module could be split further into one which advances particles/handles BCs, and one which performs reaction calculations.

We should take care to define the public APIs of each modules in such a way that implementation details are hidden as much as possible.
This was a major mistake in HallThruster.jl, which exposed the internal details of the code to the user in a way that made refactoring without changing the user interface very challenging.

## Temporal order of accuracy and operator splitting

HallThruster.jl does not at present use a principled operator splitting scheme for updating its various subcomponents.
This caps us to 1st order temporal accuracy and may be less than ideal from a stability standpoint.
Using the highly modular design approach outlined above, we can instead adopt a more rigorous splitting strategy (likely [*Strang splitting*](https://en.wikipedia.org/wiki/Strang_splitting)) to perform the update at each timestep.
This lets each module take care of its own timestepping in the way that makes then most sense for the numerics of that module, and gives us the possibility of second-order temporal accuracy.
This will prove valuable for cases in which the transient behavior of the thruster is important.

## Strong emphasis on verification and testing
Types of test
- **Unit tests**: These are tests of individual functions or "units" in the code.
    They are the easiest kind of test to write, and are essential for making sure the nuts and bolts of the code work as anticipated.
    As they grow to test larger chunks of the code, however, they can become brittle and flaky as more and more test setup is needed to "mock" parts of the code not directly under testing.
- **Property tests**: These check to make sure specific properties hold.
    When done well, these can automatically search an input space to try and find cases that violate the properties we have set.
    Whenenver possible, we should try and identify properties that either 1) always hold or 2) hold under some known set of conditions and test for these automatically.
    Some easy ones include mass conservation and one of momentum/energy conservation, depending on the PIC scheme we pic.
    Should also check current conservation.
    Smaller invariants are also worth checking, like verifying that the output of some function is always of a specific form or in a specific range
- **Order verification tests**: These use artificial source terms and the method of manufactured solutions to check that a given discretized PDE is being solved to the desired order of accuracy.
    We use these in HallThruster.jl for the ion fluid model and the electron energy equation, but we could also do this for other parts of the code, like gradient computations.
    We should try to implement these in a way that requires as little mocking as possible, as the setup code for these tests needs to maintain and can become burdensome.
- **Regression tests**: These check that the observable behavior of the code has not changed unexpectedly.
    There are three main types of regression test, each with their own challenges and uses.
    1. *Behavior regression*: This checks that the code has not changed its overall behavior.
        For us, this means the predicted thruster performance and plasma properties.
        These tests are easy to implement, but rely on you knowing exactly what observable behavior you care about, and what exactly constitutes a change.
    2. *Performance regression*: This type of test checks that we haven't made the code slower.
        This is quite challenging, as the performance of the code varies from machine to machine so we cannot easily specify a standard of performance that we can check against.
        I have not yet figured out a good way to do these sorts of tests.
    3. *Bug regression*: This sort of test verifies that a bug stays fixed.
        When we fix a bug, a test that would have triggered the bug should be added to the code, so that we can check for all time afterward that the bug stays squashed

### Axisymmetry and quasi-1D flow

It would be nice if particles could move in 2D so that we can more self-consistently model wall losses.
This could move the code in a quasi-2D direction, given some assumptions.
Axisymmetry complicates the PIC update significantly, however, so this needs to be done carefully.

### Reactions 

#### Electron Impact Reactions

Electron impact reactions, where a reactant species collides with an electron (ionization, dissociation, etc.), are treated in the following way. For each reaction, we first compute a $\Delta n$ in each cell
```math
    \Delta n = n_e n_r k_r(T_e) \Delta t,
```
where $n_e$ is the electron number density, $n_r$ is the reactant density, $k_r$ is the reaction rate which is a function of the electron temperature $T_e$, and $\Delta t$ is the timestep. This quantity represents the change in number density in the timestep due to a given reaction. We then loop over all the products to determine the proportion of this created density that is actually consumed by the particles to be generated. The consumed density is calculated as 
```math 
    \Delta n_{consumed} = min_{products}(N_{gen} * w_{gen}),
```
where the minimum is over the product species of the reaction where the number of generated particles $N_{gen}$ for that species is given by
```math 
    N_{gen} = floor(\frac{\Delta n}{w_{gen}}),
```
and the weight the new particles are generated $w_{gen}$ at is given by
```math
    w_{gen} = \Bar{w} \frac{N_{cell}}{N_{desired}}
```
where $\Bar{w}$ is the time averaged average weight of the particles of the product in the cell, $N_{cell}$ is the number of product particles in the cell, and $N_{desired}$ is a hyperparameter that specifies the number of product particles that we want in the cell. Essentailly, this generation weight means that generated particles will not be far from the average particle weight, accounting for the fact that we do desire a reasonable number of particles per cell. The reason we do the loop over the products is that the fact that $\Bar{w}$ varies between species and we use a floor function can cause total weight added to each product species to differ betwen products, which does not conserve mass. By using the minimum to find the consumed $\Delta n$ we can then ensure that the $\Delta n$ is consistent across all species. We also calculate the remainder $\Delta n$ as
```math 
    \Delta n_{remainder} = \Delta n - \Delta n_{consumed},
```
which is tracked per reaction per cell and carried over to the next timestep. In this manner we can handle cases where the reaction does not generate enough weight to create a product particle with good statistics, but ensure mass conservation on average. 

For the actual particle generatation, we then use $\Delta n_{consumed}$ to calculate the number of particles to be generated using the floor method. We then modify the generation weight as 
```math 
    w_{gen} = \frac{\Delta n_{consumed}}{N_gen},
```
which essentially purturbes the generation weight to ensure that the $\Delta n_{consumed}$ is, in fact, consumed. 

Finally, we reduce the weight of each reactant particle  by 
```math 
    \frac{w_{p,cell}}{w_{total, cell}} \Delta n_{consumed},
```
where $w_{p,cell}$ is the weight the particle contributes to the cell and $w_{total, cell}$ is the total weight of the reactant speices in the cell. In this manner, all reactant particles have their weights reduced, in proportion to their size. This insures that larger particles see a larger reduction, which is consistent with the linear dependency of $\Delta n$ on the reactant density.  

```@autodocs

Modules = [HallPIC]
```