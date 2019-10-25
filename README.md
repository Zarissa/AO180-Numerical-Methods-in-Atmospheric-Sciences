# A&O-180-Numerical-Methods-in-Atmospheric-Sciences

Adrian Ling-Ho Lam 

This repository includes a series of 7 numerical simulation projects as part of the 10-week A&O 180 Course (aka Number-crunching in Hell) that I took in my first year at UCLA. 

What all this biscuit is about in English: 
Many complex phenomena in nature can only be modeled by (partial) differential equations, which often do not have analytical solutions. Numerical methods are essentially a set of computational schemes that allow a computer to obtain a "good enough" (sometimes, the best possible) solution of a set of equations. From these projects, I learned several basic approaches (e.g. Successive over Relaxation) in solving classical PDEs (e.g. Poisson equation). Numerical artifacts are covered in depth in this course as well. 

The final project: 
The previous projects culminate in a final project of a 2D Navier-Stokes Atmosphere Model, written in the vorticity-streamfunction formation. 

What you will find: 
Specifications, codes, results, project reports and everything surviving the hurricane that is A&O180 are all included in this repo in remembrance of my deceased neurones.  

****************************************************************************************************************************
*** If you are looking for a technical summary of what I did, please check out the individual project reports. If you are here for the cool pics, see the last few projects (Proj 5- Fin)***
****************************************************************************************************************************

Note: The first few projects are hosted in Jupyter notebooks, where the results are plotted inline. For the last few projects, only the codes and some images/animations are attached due to the prohibitive computing time (>5 hrs). 

------------------------------------------------------------------
Project 1:
Project 2: Modeling the Solution of a simple one-dimension linera wave equation via Forward in Time Backward in Space scheme
Project 3: Modeling the solution of a simple one-dimensional linear wave equation via Leap-Frog and 3rd order Runge-Kutta (RK3) scheme
Project 4: Advection Diffusion Problem in 2D
Project 5: 2D Poisson's Equation
Project 6: Quasi-Geotrophic Barotropic Vorticity Equation
Project 7: Rising Thermal 
Final Project: A Numerical Simulation of the 2D Kelvin-Helmholtz Instability Problem 

The final project investigates the Kelvin-Helmholtz Instability Problem, a fluid phenomena that develops at the interface of two flows with a non-zero shear velocity. (English: Strange patterns that develop between 2 parallel flows travelling at different velocities) I also experimented with different number crunching setups (e.g. grid sizes, etc) 

Best Setup: 2D Navier-Stokes Equation in Vorticity Streamfunction formulation. Implemented RK4, RK5 and AB3 integration schemes. 

------------------------------------------------------------------
