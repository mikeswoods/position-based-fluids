/*******************************************************************************
 * Simulation.h
 * CIS563: Physically Based Animation final project
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
    float4 pos;
    float4 vel;
    float  mass;
    float  radius;
    /**
     * VERY IMPORTANT: This is needed so that the struct's size is aligned 
     * for x86 memory access along 16 byte intervals.
     *
     * If the size is not aligned, results WILL be screwed up!!! 
     * Don't be like me and waste hours trying to debug this issue. The
     * OpenCL compiler WILL NOT pad your struct to so that boundary aligned
     * like g++/clang will in the C++ world.
     *
     * See http://en.wikipedia.org/wiki/Data_structure_alignment
     */
    float  __dummy[2];
} Particle;

typedef struct {
    int particleIndex; // Index of particle in particle buffer
    int cellI;         // Corresponding grid index in the x-axis
    int cellJ;         // Corresponding grid index in the y-axis
    int cellK;         // Corresponding grid index in the z-axis
} ParticlePosition;

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
        // Count of the current frame number
        unsigned int frameNumber;
    
        // Load kernels and bind parameters
        void initializeKernels();
    
        // Moves data from GPU buffers back to the host
        void readFromGPU();
    
        // Writes data from the host to buffers on the GPU (device)
        void writeToGPU();
    
    protected:
        // OpenCL manager
        msa::OpenCL& openCL;
    
        // Bounding volume
        AABB bounds;
    
        // Timestep size:
        float dt;
    
        // Cells per axis for spatial subdivision:
        EigenVector3 cellsPerAxis;
    
        // Total number of particles in the system
        unsigned int numParticles;
    
        // Total number of cells in the system
        unsigned int numCells;
    
        // Just as the name describes
        float massPerParticle;

        // Host managed buffer of particles; adapted from the of ofxMSAOpenCL
        // particle example
        msa::OpenCLBufferManagedT<Particle>	        particles;
        msa::OpenCLBufferManagedT<ParticlePosition>	particleToCell;
        //msa::OpenCLBufferManagedT<int>	sortedParticleToCell;
        msa::OpenCLBufferManagedT<int>	            cellHistogram;

        // Initialization-related functions:
        void initialize();

        // Simulation state-related functions:
        void applyExternalForces();
        void predictPositions();
        void discretizeParticlePositions();
        void zeroCellHistogram();
        void computeCellHistogram();
        void sortParticlesByCell();
    
        // Drawing-related functions:
        void drawBounds() const;
        void drawParticles();
        void countingSort();

    public:
        Simulation(msa::OpenCL& openCL
                  ,AABB bounds
                  ,float dt = 0.025f
                  ,unsigned int numParticles = 1000
                  ,float massPerParticle = 1.0f);
        virtual ~Simulation();

        const unsigned int getFrameNumber() const { return this->frameNumber; }
        const AABB& getBounds() const { return this->bounds; }
        const unsigned int getNumberOfParticles() const { return this->numParticles; }
        const unsigned int getNumberOfCells() const { return this->numCells; }
    
        void reset();
        void step();
        void draw();
    
        friend std::ostream& operator<<(std::ostream& os, EigenVector3 v);
};

#endif
