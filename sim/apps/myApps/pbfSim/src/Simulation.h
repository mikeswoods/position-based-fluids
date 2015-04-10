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
    private:
        void loadKernels();
    
    protected:
        // OpenCL manager
        msa::OpenCL& openCL;
    
        // Bounding volume
        AABB bounds;
    
        // Timestep size:
        float dt;

        // Total number of particles in the system
        unsigned int numParticles;
    
        // Just as the name describes
        float massPerParticle;

        // Host managed buffer of particles; adapted from the of ofxMSAOpenCL
        // particle example
        msa::OpenCLBufferManagedT<Particle>	particles;
    
        void initialize();

        void computeExternalForce();
    
        void drawBounds() const;
        void drawParticles();

    public:
        Simulation(msa::OpenCL& openCL
                  ,AABB bounds
                  ,float dt
                  ,unsigned int numParticles
                  ,float massPerParticle);
        virtual ~Simulation();
    
        const AABB& getBounds() const { return this->bounds; }
    
        void reset();
        void update();
        void draw();
    
        friend std::ostream& operator<<(std::ostream& os, EigenVector3 v);
};

#endif
