/*******************************************************************************
 * Constants.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Parameters.h"

/******************************************************************************/

// If defined, verbose-level logging of system events to the console will
// be enabled
//#define ENABLE_LOGGING 1

// If defined, a simple test scene will be used for rendering
//#define SIMPLE_SCENE 1

// If defined, mesh spheres will be drawn for the particles, otherwise
// faster OpenGL points will be used:

//#define DRAW_PARTICLES_AS_SPHERES 1

/******************************************************************************/

namespace Constants {

/**
 * Number of iterations the Jacobi solver will run for
 */
const int SOLVER_ITERATIONS = 4;
   
/**
 * Default number of particles in the simulation
 */
const int DEFAULT_NUM_PARTICLES = 5000;

/**
 * Default time step = 1/30 of a second
 */
//const float DEFAULT_DT = 0.04166; // 1/24 s
//const float DEFAULT_DT = 0.033f; // 1/30 s
const float DEFAULT_DT = 0.025f; // 1/40 s
//const float DEFAULT_DT = 0.0166f; // 1/60 s

/**
 * Desired number of particles per cell in X. This value is used to 
 * compute the optimal number of cells in the spatial grid per axis
 */
const int PARTICLES_PER_CELL_X = 2;

/**
 * Desired number of particles per cell in Y This value is used to
 * compute the optimal number of cells in the spatial grid per axis
 */
const int PARTICLES_PER_CELL_Y = 2;

/**
 * Desired number of particles per cell in Z This value is used to
 * compute the optimal number of cells in the spatial grid per axis
 */
const int PARTICLES_PER_CELL_Z = 2;

/******************************************************************************/

/**
 * Parameters for a simulation with particles of radius 0.5 units
 */
const Parameters FOR_RADIUS_0_5 = Parameters(0.5f    // Particle radius
                                            ,1.15f   // Smoothing radius
                                            ,0.0033f // Relaxation
                                            ,0.1f    // Artifical Pressure K
                                            ,4.0f    // Artifical Pressure N
                                            ,0.1f    // Vorticity epsilon
                                            ,0.01f); // Viscosity coefficient

/**
 * Parameters for a simulation with particles of radius 0.25 units
 */
const Parameters FOR_RADIUS_0_25 = Parameters(0.25f   // Particle radius
                                             ,0.75f    // Smoothing radius
                                             ,0.1f    // Relaxation
                                             ,10.0f  // Artifical Pressure K
                                             ,4.0f    // Artifical Pressure N
                                             ,0.1f    // Vorticity epsilon
                                             ,0.01f); // Viscosity coefficient
    
}
    
#endif
