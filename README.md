# radVP : *Simulating the spherically symmetric gravitational Vlasov-Poisson system*
by Christopher Straub
> This is the code used for the numerical analysis conducted in Chapter 8 of the doctoral thesis 'Pulsating Galaxies' as well as for the paper 'Numerical experiments on stationary, oscillating, and damped galaxy models'. 

## Citation
If you would use, adapt, or extend this code, its author would be most pleased. To let him know that you have done so, please cite the doctoral thesis in which the code was originally conceived.
```
@article{Straub24
    title = {Pulsating Galaxies},
    author = {Christopher Straub},
    school = {University of Bayreuth},
    year = 2024,
    doi = {https://doi.org/10.15495/EPub_UBT_00007639},
}
```

## Uses
The code can be used to numerically compute steady states of the spherically symmetric gravitational Vlasov-Poisson system and their radial periods.
In addition, it includes Particle-In-Cell simulations of the linearised and non-linearised radial Vlasov-Poisson system.

## Style of the Code
The code is written in C++ and is kept as simple/non-sophisticated as possible. In particular, no specialised programming libraries are used, everything is written from scratch. This is done so that the code is also understandable for less programming experienced readers.

## Technical Instructions
The code was written in Code::Blocks 20.03 and compiled with C++17.
A makefile we used for compilation (under Ubuntu 18.04, g++ version 7.5.0) is included in a separate folder.

In the main folder, you can find five parameter files, all having the prefix "in_".
These files determine the parameters used for the particle-in-cell simulations of the Non-Linearised, Linearised, or Pure Transport system.
The different parameters are described in detail in these files; see also the explanations given in the thesis.

The five parameter files in the main folder lead to a simulation of the linearised VP for the steady state k=1=Rmax and initial condition w*phi'(E).
In the folder "Parameters", we have included further parameter files.
The first ones result in a similar simulation, but with better numerical parameters; these paramters were actually used to create the plots at the start of Section 8.3 of the thesis.
The second ones initiate a simulation of the non-linearised VP with initial condition close to the isotropic polytrope k=1=Rmax; these parameters were used to create the plots at the start of Section 8.4.

The outputs of the particle-in-cell simulations get printed in various output files, all having the prefix "out_". 

In the parallelised parts (in particular, in the PIC simulations), the program uses 200 threads in parallel.
This number can be changed in the ..."_barrier.h" Header Files.


The numerical simulations of a steady state itself and its radial periods have to be called by writing suitable function calls in main.cpp.
Examples for how to call these functions -- again in the case of the istropic polytrope k=1=Rmax -- are provided in main.cpp ; they are currently commented out.



