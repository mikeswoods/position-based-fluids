/*******************************************************************************
 * Constants.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace Constants {

/**
 * Default time step = 1/30 of a second
 */
const float DEFAULT_DT = 0.033f;

/**
 * Default particle mass
 */
const float DEFAULT_PARTICLE_MASS = 1.0f;

/**
 * Max particle mass
 */
const float MAX_PARTICLE_MASS = 3.0f;
    
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
