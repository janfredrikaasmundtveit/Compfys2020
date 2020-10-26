This repository contains all programs relevant to project 3 Fys-3150.

The file planet.h contains several classes needed for the program. theese clases are:
  coord: difines a coordinate with 1 x-component and 1 y component.'
  
  planet: defines a planet which contains 3 coord and 3 doubles. The doubles are the mass, the energy and the angular momentum.(as it is only used for 2D system angular momentum is a scalar).
  The 3 coords are: the posistion, the velocity and the total force.
  
  solarsystem: contains an integer size which is the number of planets, including the sun, and one array of 10 planets. 
  The solarsystem is initalized without any planets. Planets are included using the member function add planet which takes an interger j and a planet and adds the replaces the planet number j in the array with the planet included. This interger should be smaller than the size and planets should be added in order starting form 0.
  The number 10 only serves as an upper bound to the number of planets. 
  solarsystem all so contains member functions for updating the force of a member planet.
  
  solver: contains an interger method and a solarsystem. int method is set to 0 if the solver should use the Verlet method  or 1 if Eulers method is desired. 
  The member function step evolves the solarsystem by a single timestep using the spesified method.

The c++ file setup.cpp contains all the relevant setups of different subsystems, aswell as the full solarsystem with realistic inital condition.

The full program takes 5 inputs: first a filename, second the number of timsteps, third second a the number of years the system should be evolved over (double), fourth   the interger option, fifth the interger spesifiing the method th solver should usesee solver class.
the integer option  takes 6 different values: choose system 0 for full solarsystem, system 1 for sun-earth sytem, system 2 for sun-earth-jupiter system
		//, system 3 for sun mercury with GR correction, system 4 to vary beta and output E,l; option 5 to study escape velocities.

