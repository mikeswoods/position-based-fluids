/*******************************************************************************
 * Simulation.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef SIMULATION_H
#define SIMULATION_H

#include <iostream>
#include "Parameters.h"
#include "Constants.h"
#include "AABB.h"
#include "MSAOpenCL.h"

/******************************************************************************/

// A particle type:

typedef struct {

    float4 pos;      // Current particle position (x)

    float4 posStar;  // Predicted particle position (x*)

    float4 vel;      // Current particle velocity (v)
    
    float4 velDelta; // Velocity delta from XSPH viscosity (v*), 4 words

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
    //float  __padding[2]; // Padding

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
    
    int __padding[2]; // Padding
    
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

        // Total number of cells in the system
        int numCells;
    
        // Flag to draw the spatial grid
        bool flagDrawGrid;
    
        // Flag for visual debugging
        bool flagVisualDebugging;
    
        // Given a particle count, particle radius and world bounds,
        // find the "ideal" cell count per axis
        ofVec3f findIdealParticleCount();

        // Moves data from GPU buffers back to the host
        void readFromGPU();
    
        // Writes data from the host to buffers on the GPU (device)
        void writeToGPU();
    
    protected:
        // Particle mesh sphere
        ofMesh particleMesh;
    
        // Particle vertices
        ofVbo particleVertices;
    
        // OpenCL manager
        msa::OpenCL& openCL;
    
        // Basic shader
        ofShader shader;
    
        // Bounding volume
        AABB bounds;
    
        // Timestep size:
        float dt;
    
        // Cells per axis for spatial subdivision:
        ofVec3f cellsPerAxis;

        // Total number of particles in the system
        int numParticles;

        // Simulation parameters to pass to the kernels
        Parameters parameterData;
        msa::OpenCLBufferManagedT<Parameters> parameters;
    
        // All particles in the simulation
        msa::OpenCLBufferManagedT<Particle>	particles;
    
        // An array of particle-to-cell mappings
        msa::OpenCLBuffer particleToCell;
    
        // A cell count histogram used for particle neighbor finding
        msa::OpenCLBuffer cellHistogram;
    
        // A sorted version of particleToCell, used to search for a given
        // particle's neighbors
        msa::OpenCLBuffer sortedParticleToCell;

        // An array of cell start locations and spans in sortedParticleToCell
        msa::OpenCLBuffer gridCellOffsets;
    
        // Particle densities computed by SPH estimation
        msa::OpenCLBuffer density;

        // Particle density lambda value from the section "Enforcing
        // Incompressibility" of "Position Based Fluids"
        msa::OpenCLBuffer lambda;
    
        // Vorticity curl force applied to each particle
        msa::OpenCLBuffer curl;
    
        // Position deltas
        msa::OpenCLBuffer posDelta;

        // Final render position for OpenCL <-> OpenGL instanced rendering
        msa::OpenCLBufferManagedT<float4> renderPos;
    
        // Initialization-related functions:
        void initialize();
        void initializeKernels();
        void initializeOpenGL();

        // Particle sorting functions:
        void computeCellHistogram();
        void discretizeParticlePositions();
        void sortParticlesByCell();
    
        // Simulation state-related functions:
        void resetQuantities();
        void applyExternalForces();
        void predictPositions();
        void groupParticlesByCell();
        void calculateDensity();
        void calculatePositionDelta();
        void applyPositionDelta();
        void handleCollisions();
        void updateVelocity();
        void applyXSPHViscosity();
        void applyVorticityConfinement();
        void updatePosition();
    
        // Drawing-related functions:
        void drawBounds(const ofCamera& camera) const;
        void drawGrid(const ofCamera& camera) const;
        void drawParticles(const ofCamera& camera);

    public:
        Simulation(msa::OpenCL& openCL
                  ,AABB bounds
                  ,int numParticles
                  ,Parameters parameters);
    
        Simulation(msa::OpenCL& openCL
                  ,AABB bounds
                  ,int numParticles
                  ,float dt
                  ,ofVec3f cellsPerAxis
                  ,Parameters parameters);

        virtual ~Simulation();

        const unsigned int getFrameNumber() const { return this->frameNumber; }
        const AABB& getBounds() const { return this->bounds; }
        const ofVec3f& getCellsPerAxis() const { return this->cellsPerAxis; }
        const unsigned int getNumberOfParticles() const { return this->numParticles; }
        const unsigned int getNumberOfCells() const { return this->numCells; }
        Parameters& getParameters() { return this->parameters[0]; }

        const bool drawGridEnabled() const { return this->flagDrawGrid; }
        void toggleDrawGrid() { this->flagDrawGrid = !this->flagDrawGrid; }
    
        const bool isVisualDebuggingEnabled() const { return this->flagVisualDebugging; }
        void toggleVisualDebugging() { this->flagVisualDebugging = !this->flagVisualDebugging; }
    
        void step();
        void draw(const ofCamera& camera);
};

#endif
