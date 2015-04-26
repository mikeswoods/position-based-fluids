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

/******************************************************************************/

namespace Constants {

/**
 * Number of iterations the Jacobi solver will run for
 */
const int SOLVER_ITERATIONS = 3;
   
/**
 * Default number of particles in the simulation
 */
const int DEFAULT_NUM_PARTICLES = 5000;

/**
 * Default time step = 1/30 of a second
 */
//const float DEFAULT_DT = 0.04166; // 1/24 s
const float DEFAULT_DT = 0.033f; // 1/30 s
//const float DEFAULT_DT = 0.025f; // 1/40 s
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
const Parameters DEFAULT_PARAMS = Parameters(0.5f    // Particle radius
                                            ,1.2f   // Smoothing radius
                                            ,0.0033f // Relaxation
                                            ,0.01f    // Artifical Pressure K
                                            ,4.0f    // Artifical Pressure N
                                            ,0.6f    // Vorticity epsilon
                                            ,0.01f); // Viscosity coefficient

}
    
#endif
