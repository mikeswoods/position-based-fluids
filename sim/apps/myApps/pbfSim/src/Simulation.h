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

// A particle type:

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
    float  __dummy[2]; // Padding

} Particle;

// A type to represent the position of a given particle in the spatial
// grid the simulated world is divided into

typedef struct {

    int particleIndex; // Index of particle in particle buffer

    int cellI;         // Corresponding grid index in the x-axis

    int cellJ;         // Corresponding grid index in the y-axis

    int cellK;         // Corresponding grid index in the z-axis

} ParticlePosition;

// A type that encodes the start and length of a grid cell in sortedParticleToCell

typedef struct {
    
    int  start; // Start of the grid cell in sortedParticleToCell
    
    int length;
    
    int __dummy[2]; // Padding
    
} GridCellOffset;

/******************************************************************************/

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
    
        // Flag to draw the spatial grid
        bool flagDrawGrid;
    
        // Flag for visual debugging
        bool flagVisualDebugging;
    
        // Load kernels and bind parameters
        void initializeKernels();
    
        // Moves data from GPU buffers back to the host
        void readFromGPU();
    
        // Writes data from the host to buffers on the GPU (device)
        void writeToGPU();
    
        // Individual functions needed for grouping particles into bins/cells
        void initializeParticleSort();
        void computeCellHistogram();
        void discretizeParticlePositions();
        void sortParticlesByCell();
    
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
        int numParticles;
    
        // Total number of cells in the system
        int numCells;
    
        // Just as the name describes
        float massPerParticle;

        // Host managed buffer of particles; adapted from the of ofxMSAOpenCL
    
        // All particles in the simulation
        msa::OpenCLBufferManagedT<Particle>	particles;
    
        // An array of particle-to-cell mappings
        msa::OpenCLBufferManagedT<ParticlePosition>	particleToCell;
    
        // A cell count histogram used for particle neighbor finding
        msa::OpenCLBufferManagedT<int> cellHistogram;
    
        // A sorted version of particleToCell
        msa::OpenCLBufferManagedT<ParticlePosition>	sortedParticleToCell;
    
        // An array of cell start locations and spans in sortedParticleToCell
        msa::OpenCLBufferManagedT<GridCellOffset> gridCellOffsets;
    
        // Particle densities computed by SPH estimation
        msa::OpenCLBufferManagedT<float> density;
    
        // Initialization-related functions:
        void initialize();
        void resetParticleQuantities();

        // Simulation state-related functions:
        void applyExternalForces();
        void predictPositions();
        void groupParticlesByCell();
        void calculateDensity();
        void calculatePositionDelta();
        void handleCollisions();
    
        // Drawing-related functions:
        void drawBounds() const;
        void drawGrid() const;
        void drawParticles();

    public:
        Simulation(msa::OpenCL& openCL
                  ,AABB bounds
                  ,float dt = 0.025f
                  ,EigenVector3 _cellsPerAxis = EigenVector3(2, 2, 2)
                  ,int numParticles = 3
                  ,float massPerParticle = 1.0f);
        virtual ~Simulation();

        const unsigned int getFrameNumber() const { return this->frameNumber; }
        const AABB& getBounds() const { return this->bounds; }
        const unsigned int getNumberOfParticles() const { return this->numParticles; }
        const unsigned int getNumberOfCells() const { return this->numCells; }
    
        const bool drawGridEnabled() const { return this->flagDrawGrid; }
        void toggleDrawGrid() { this->flagDrawGrid = !this->flagDrawGrid; }
    
        const bool isVisualDebuggingEnabled() const { return this->flagVisualDebugging; }
        void toggleVisualDebugging() { this->flagVisualDebugging = !this->flagVisualDebugging; }
    
        void reset();
        void step();
        void draw();
    
        friend std::ostream& operator<<(std::ostream& os, EigenVector3 v);
};

#endif
