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
    1. *Behavior regression*
        This checks that the code has not changed its overall behavior.
        For us, this means the predicted thruster performance and plasma properties.
        These tests are easy to implement, but rely on you knowing exactly what observable behavior you care about, and what exactly constitutes a change.
    2. *Performance regression*
        This type of test checks that we haven't made the code slower.
        This is quite challenging, as the performance of the code varies from machine to machine so we cannot easily specify a standard of performance that we can check against.
        I have not yet figured out a good way to do these sorts of tests.
    3. *Bug regression*
        This sort of test verifies that a bug stays fixed.
        When we fix a bug, a test that would have triggered the bug should be added to the code, so that we can check for all time afterward that the bug stays squashed

### Axisymmetry and quasi-1D flow

It would be nice if particles could move in 2D so that we can more self-consistently model wall losses.
This could move the code in a quasi-2D direction, given some assumptions.
Axisymmetry complicates the PIC update significantly, however, so this needs to be done carefully.
