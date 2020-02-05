# A&O-180-Numerical-Methods-in-Atmospheric-Sciences

This repository includes a series of 7 numerical simulation projects as part of the A&O 180 Course that I took in my first year at UCLA. 

What you will find: 
Specifications, codes, results, project reports and everything surviving the hurricane that is A&O180 are all included in this repo. 

*************************************************************************************************************************

**Important Notes:**
                        
**For a technical summary of what I did, please check out the individual project reports.**

**For cool simulation pictures, see the last few projects (Project 5 till Final Project) and scroll down the Jupyter notebooks**

**Use https://nbviewer.jupyter.org/ to view the jupyter notebooks if you have issues displaying them**
*************************************************************************************************************************


What this is all about in non-technical speak: 

Many complex phenomena in nature can only be modeled by (partial) differential equations, which often do not have analytical solutions. Numerical methods are essentially a set of computational schemes that allow a computer to obtain a "good enough" (sometimes, the best possible) solution of a set of DEs. From these projects, several numerical algorithms (e.g. Successive-Over-Relaxation for solving the Poisson Equation) and numerical schemes (e.g.RK,AB families) for solving PDEs. 


The projects also touch on several mathematical constraints (e.g.stability) and numerical artifacts (e.g.virtual advection) that are intrinsic in some specific numerical methods used. 


The final project: 
The previous projects culminate in a final project of a 2D Navier-Stokes Atmosphere Model, written in the vorticity-streamfunction formation. It is based on the Project 6 QG code. 

************************************************************************************************************************


Note: The first few projects are hosted in Jupyter notebooks, where the results are plotted inline. For the last few projects, only the codes and some images/animations are attached due to the prohibitive computing time (>5 hrs). 


------------------------------------------------------------------

Project 1: Simulation of an air packet under the Buoyancy Equation (Integration Scheme: Euler Forward & Leap Frog :'))

Project 2: Modeling the Solution of a simple one-dimension linera wave equation via Forward in Time Backward in Space scheme

Project 3: Modeling the solution of a simple one-dimensional linear wave equation via Leap-Frog and 3rd order Runge-Kutta (RK3) scheme

Project 4: Advection Diffusion Problem in 2D

Project 5: 2D Poisson's Equation

Project 6: Quasi-Geotrophic Barotropic Vorticity Equation

Project 7: Rising Thermal 

Final Project: A Numerical Simulation of the 2D Kelvin-Helmholtz Instability Problem 

The final project investigates the Kelvin-Helmholtz Instability Problem, a flow instability that develops at the interface of two flows with a shear velocity. (Strange patterns that develop between 2 parallel flows travelling at different speeds) I also experimented with different number crunching setups (e.g. grid sizes, etc) 

Best Setup: 2D Navier-Stokes Equation in Vorticity Streamfunction formulation. Implemented RK4, RK5 and AB3 integration schemes. 


------------------------------------------------------------------
