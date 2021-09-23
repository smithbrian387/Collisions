# Collisions
Simulating elastic collisions between round objects. Collisions take into account particle mass.

Installation notes
Required libraries: pygame, decimal, math, sys, random

General notes:
Velocity and radii sliders affect new particles being added, not existing particles.

Full screens can cause glitches due to a series of particles overlapping as they get repositioned through collisions. Screen area fillable by particles is capped to reduce this.

Collisions are detected within sections with size based on the largest particle currently on the screen. Ex: 500 small particles will run much more efficiently than 500 small particles with 1 large particle.


