# Mass Spring Modal Analysis

Sample project simulating a 3d mass spring system using Backwards Euler time integration.

Includes code for linear modal analysis to find the vibration modes for the object.

The math for this project is here: https://paper.dropbox.com/doc/Mass-Spring-System-LmksbhzlJnmrIFlKN7nH6

Using:
 - Eigen for linear algebra
 - Spectra for the generalized eigenvalue solve
 - LibIGL for graphics

Make Instructions:

..../libigl

..../spectra-0.5.0

..../MassSpringModalAnalysis

Go into

 MassSpringModalAnalysis

$ mkdir build

$ cd build

$ cmake ..

$ make

