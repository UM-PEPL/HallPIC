# Reactions 

## Electron Impact Reactions

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