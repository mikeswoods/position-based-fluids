/*******************************************************************************
 * Simulation.h
 * CIS563: Physcially Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include "Definitions.h"
#include "AABB.h"
#include "MSAOpenCL.h"

/******************************************************************************/

typedef struct {
    float3 x;
    float3 v;
    float mass;
} Particle;

/**
 * This class encompasses the current statue of the 
 * Position-Based Fluids/Dynamics system at a given point in time. Much of the
 * code defining the implementation of this class was originally derived
 * from the second assignment in the class, which in turn, was based off of
 * Matthias Muller's "Position Based Dynamics" paper
 */
class Simulation
{
    protected:
        // Bounding volume
        AABB bounds;

        // Total number of particles in the system
        unsigned int numParticles;
    
        // Just as the name describes
        float massPerParticle;
    
        // A sparse (3k) x (3k) matrix encoding the mass of the i-th particle
        // at the (i,i)-th entry in the matrix
        SparseMatrix particleMass;

        // A sparse (3k) x (3k) matrix encoding the inverse mass of the i-th
        // particle at the (i,i)-th entry in the matrix
        SparseMatrix invParticleMass;

        // A k x 3 matrix, where the i-th block encodes the (x,y,z) components
        // of the position of the i-th particle
        VectorX particleX;

        // A k x 3 matrix, where the i-th block encodes the (x,y,z) components
        // of the velocity of the i-th particle
        VectorX particleV;

        void initialize();
        void integratePBD(VectorX& x, VectorX& v, float dt, unsigned int ns);
        void computeExternalForce();

    public:
        Simulation(AABB bounds, unsigned int numParticles, float massPerParticle);
        virtual ~Simulation();
    
        const AABB& getBounds() const { return this->bounds; }
    
        void reset();
        void update();
        void draw();
    
        friend std::ostream& operator<<(std::ostream& os, EigenVector3 v);
};

#endif
