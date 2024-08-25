# TraceGR
Submission for MIT-IEEE URTC '24 | Modular Spacetime Imaging using Null Geodesics in Raycasting Applications.

## Abstract

This study investigates the use of null geodesic photon trajectories as rays in a raycasting system to visualize spacetime warping (gravitational lensing) caused by celestial bodies. These geometries are defined via an isotropic line element or metric tensor, serving as a modular input for the raycasting program. Our approach enables the simulation and visualization of various theoretical spacetimes. We evaluate its performance on imaging a Schwarzschild black hole (stationary black hole) and the Alcubierre warp drive (geometry theorized to enable faster-than-light travel). Our results demonstrate the accuracy and versatility of this technique in rendering complex spacetime geometries and provide a valuable tool for theoretical and observational astrophysics. We expand on future work in the use of gravitational wave data to generate visualizations of real binary black hole systems, as well as the improvement of computational efficiency of the raycasting system.

## Results

Below are a few results,


### 3D Black Hole (Schwarzschild Isotropic Metric)

The line element for this render is

$ds^2 = \left( \frac{1 - \frac{r_s}{4R}}{1 + \frac{r_s}{4R}} \right)^2 dt^2 - \left( 1 + \frac{r_s}{4R} \right)^4 (dx^2 + dy^2 + dz^2)$




### 3D Alcubierre Warp Drive (Alcubierre Metric)

The line element for Alcubierre renders is

$ds^2 = -c^2 dt^2 + (dx - v_s f(r) dt)^2 + dy^2 + dz^2$



### 3D Moving Alcubierre Warp Drive (Alcubierre Metric)

This is a raytraced video for an Alcubierre Warp Drive

<div align="center">
  <img src="https://github.com/zanebeeai/TraceGR/blob/main/Results/AWD/bg4/bg4WARPDRIVE_final.gif" height="300">
</div>


## Algorithms

Key algorithms are given in the ```~/Functions/``` folder. Some examples include solving the metric tensor from a line element, solving the Christoffel symbols and solving geodesics.

## How to run

To run, create a new python virtual env and run

```pip install -r requirements.txt```

Then follow the notebooks in ```~/Notebooks/```
