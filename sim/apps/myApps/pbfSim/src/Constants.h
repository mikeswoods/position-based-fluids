/*******************************************************************************
 * Constants.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef PBF_SIM_CONSTANTS_H
#define PBF_SIM_CONSTANTS_H

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

#if defined(SIMPLE_SCENE) && !defined(DRAW_PARTICLES_AS_SPHERES)
    #define DRAW_PARTICLES_AS_SPHERES 1
#endif

// If defined, the simulation will produce better results at the expense of
// speed

#define MORE_ACCURATE 1

/******************************************************************************/

namespace Constants {

/**
 * Number of iterations the Jacobi solver will run for
 */
const int SOLVER_ITERATIONS = 3;
   
/**
 * Default number of particles in the simulation
 */
const int DEFAULT_NUM_PARTICLES = 10000;

/**
 * Default time step
 */

#ifdef MORE_ACCURATE
const float DEFAULT_DT = 0.025f; // 1/40 s
#else
const float DEFAULT_DT = 0.033f; // 1/30 s
#endif

/**
 * Desired number of particles per cell in X. This value is used to 
 * compute the optimal number of cells in the spatial grid per axis
 */
#ifdef MORE_ACCURATE
const int PARTICLES_PER_CELL_X = 3;
#else
const int PARTICLES_PER_CELL_X = 2;
#endif

/**
 * Desired number of particles per cell in Y This value is used to
 * compute the optimal number of cells in the spatial grid per axis
 */
#ifdef MORE_ACCURATE
const int PARTICLES_PER_CELL_Y = 3;
#else
const int PARTICLES_PER_CELL_Y = 2;
#endif

/**
 * Desired number of particles per cell in Z This value is used to
 * compute the optimal number of cells in the spatial grid per axis
 */
#ifdef MORE_ACCURATE
const int PARTICLES_PER_CELL_Z = 3;
#else
const int PARTICLES_PER_CELL_Z = 2;
#endif

/******************************************************************************/

/**
 * Parameters for a simulation with particles of radius 0.5 units
 */
const Parameters DEFAULT_PARAMS = Parameters(0.5f    // Particle radius
                                            ,1.15f   // Smoothing radius
                                            ,0.0033f // Relaxation
                                            ,0.1f    // Artifical Pressure K
                                            ,5.0f    // Artifical Pressure N
                                            ,0.1f    // Vorticity epsilon
                                            ,0.01f); // Viscosity coefficient

}
    
#endif
