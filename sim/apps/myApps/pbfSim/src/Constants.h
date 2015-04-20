/*******************************************************************************
 * Constants.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

/******************************************************************************/

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
const int SOLVER_ITERATIONS = 3;
   
/**
 * Default number of particles in the simulation
 */
const int DEFAULT_NUM_PARTICLES = 5000;

/**
 * Default time step = 1/30 of a second
 */
const float DEFAULT_DT = 0.033f;

/**
 * Max particle radius
 */
const float MAX_PARTICLE_RADIUS = 3.0f;
    
/**
 * Default particle radius
 */
const float DEFAULT_PARTICLE_RADIUS = 0.1f;
    
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

}
    
#endif
