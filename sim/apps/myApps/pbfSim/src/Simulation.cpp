/*******************************************************************************
 * Simulation.cpp
 * - The heart of the position-based fluids simulator. This class encapsulates
 *   the current state of the simulation
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include "ofMain.h"
#include "Simulation.h"

/******************************************************************************/

const int SOLVER_ITERATIONS = 5;

/******************************************************************************/

using namespace std;

/******************************************************************************/

ostream& operator<<(ostream& os, EigenVector3 v)
{
    return os << "<" << v[0] << "," << v[1] << "," << v[2] << ">";
}

ostream& operator<<(ostream& os, Particle p)
{
    return os << "Particle {" << endl
              << "  pos: <" << p.pos.x << "," << p.pos.y << "," << p.pos.z << ">" << endl
              << "  vel: <" << p.vel.x << "," << p.vel.y << "," << p.vel.z << ">" << endl
              << "  mass: " << p.mass << endl
              << " radiusL " << p.radius << endl
              << "}";
}

/******************************************************************************/

/**
 * Constructs a new simulation instance
 *
 * @param [in] _openCL OpenCL manager instance
 * @param [in] _bounds Defines the boundaries of the simulation in world space
 * @param [in] _dt The time step (usually 0.1)
 * @param [in] _cellsPerAxis Cell spatial grid subdivisions per axis
 * @param [in] _numParticles The number of particles in the simulation
 * @param [in] _radiusPerParticle The mass per particle
 * @param [in] _particleMass The mass per particle
 */
Simulation::Simulation(msa::OpenCL& _openCL
                      ,AABB _bounds
                      ,float _dt
                      ,EigenVector3 _cellsPerAxis
                      ,int _numParticles
                      ,float _particleRadius
                      ,float _particleMass) :
    openCL(_openCL),
    bounds(_bounds),
    dt(_dt),
    cellsPerAxis(_cellsPerAxis),
    numParticles(_numParticles),
    particleRadius(_particleRadius),
    particleMass(_particleMass),
    frameNumber(0),
    flagDrawGrid(false),
    flagVisualDebugging(false)
{
    this->numCells =   static_cast<int>(this->cellsPerAxis[0])
                     * static_cast<int>(this->cellsPerAxis[1])
                     * static_cast<int>(this->cellsPerAxis[2]);

    this->initialize();
}

Simulation::~Simulation()
{
    
}

/**
 * Loads all of the OpenCL kernels that will be used for during the simulation
 */
void Simulation::initializeKernels()
{
    // Read the source files for the kernels:
    this->openCL.loadProgramFromFile("kernels/Simulation.cl");
    
    // KERNEL :: applyExternalForces
    this->openCL.loadKernel("applyExternalForces");
    this->openCL.kernel("applyExternalForces")->setArg(0, this->particles);
    this->openCL.kernel("applyExternalForces")->setArg(1, this->dt);
    
    // KERNEL :: predictPosition
    this->openCL.loadKernel("predictPosition");
    this->openCL.kernel("predictPosition")->setArg(0, this->particles);
    this->openCL.kernel("predictPosition")->setArg(1, this->dt);
    
    // KERNEL :: discretizeParticlePositions
    this->openCL.loadKernel("discretizeParticlePositions");
    this->openCL.kernel("discretizeParticlePositions")->setArg(0, this->particles);
    this->openCL.kernel("discretizeParticlePositions")->setArg(1, this->particleToCell);
    this->openCL.kernel("discretizeParticlePositions")->setArg(2, this->cellHistogram);
    this->openCL.kernel("discretizeParticlePositions")->setArg(3, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(4, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(5, static_cast<int>(this->cellsPerAxis[2]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(6, this->bounds.getMinExtent());
    this->openCL.kernel("discretizeParticlePositions")->setArg(7, this->bounds.getMaxExtent());
    
    // KERNEL :: sortParticlesByCell
    this->openCL.loadKernel("sortParticlesByCell");
    this->openCL.kernel("sortParticlesByCell")->setArg(0, this->particleToCell);
    this->openCL.kernel("sortParticlesByCell")->setArg(1, this->cellHistogram);
    this->openCL.kernel("sortParticlesByCell")->setArg(2, this->sortedParticleToCell);
    this->openCL.kernel("sortParticlesByCell")->setArg(3, this->gridCellOffsets);
    this->openCL.kernel("sortParticlesByCell")->setArg(4, this->numParticles);
    this->openCL.kernel("sortParticlesByCell")->setArg(5, this->numCells);
    this->openCL.kernel("sortParticlesByCell")->setArg(6, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("sortParticlesByCell")->setArg(7, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("sortParticlesByCell")->setArg(8, static_cast<int>(this->cellsPerAxis[2]));
    
    // KERNEL :: estimateDensity
    this->openCL.loadKernel("estimateDensity");
    this->openCL.kernel("estimateDensity")->setArg(0, this->particles);
    this->openCL.kernel("estimateDensity")->setArg(1, this->sortedParticleToCell);
    this->openCL.kernel("estimateDensity")->setArg(2, this->gridCellOffsets);
    this->openCL.kernel("estimateDensity")->setArg(3, this->numParticles);
    this->openCL.kernel("estimateDensity")->setArg(4, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("estimateDensity")->setArg(5, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("estimateDensity")->setArg(6, static_cast<int>(this->cellsPerAxis[2]));
    this->openCL.kernel("estimateDensity")->setArg(7, this->density);
    
    // KERNEL :: computeLambda
    this->openCL.loadKernel("computeLambda");
    this->openCL.kernel("computeLambda")->setArg(0, this->particles);
    this->openCL.kernel("computeLambda")->setArg(1, this->sortedParticleToCell);
    this->openCL.kernel("computeLambda")->setArg(2, this->gridCellOffsets);
    this->openCL.kernel("computeLambda")->setArg(3, this->density);
    this->openCL.kernel("computeLambda")->setArg(4, this->numParticles);
    this->openCL.kernel("computeLambda")->setArg(5, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("computeLambda")->setArg(6, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("computeLambda")->setArg(7, static_cast<int>(this->cellsPerAxis[2]));
    this->openCL.kernel("computeLambda")->setArg(8, this->lambda);
    
    // KERNEL :: computePositionDelta
    this->openCL.loadKernel("computePositionDelta");
    this->openCL.kernel("computePositionDelta")->setArg(0, this->particles);
    this->openCL.kernel("computePositionDelta")->setArg(1, this->sortedParticleToCell);
    this->openCL.kernel("computePositionDelta")->setArg(2, this->gridCellOffsets);
    this->openCL.kernel("computePositionDelta")->setArg(3, this->numParticles);
    this->openCL.kernel("computePositionDelta")->setArg(4, this->lambda);
    this->openCL.kernel("computePositionDelta")->setArg(5, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("computePositionDelta")->setArg(6, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("computePositionDelta")->setArg(7, static_cast<int>(this->cellsPerAxis[2]));
    this->openCL.kernel("computePositionDelta")->setArg(8, this->posDeltaX);
    this->openCL.kernel("computePositionDelta")->setArg(9, this->posDeltaY);
    this->openCL.kernel("computePositionDelta")->setArg(10, this->posDeltaZ);
    
    // KERNEL :: applyPositionDelta
    this->openCL.loadKernel("applyPositionDelta");
    this->openCL.kernel("applyPositionDelta")->setArg(0, this->posDeltaX);
    this->openCL.kernel("applyPositionDelta")->setArg(1, this->posDeltaY);
    this->openCL.kernel("applyPositionDelta")->setArg(2, this->posDeltaZ);
    this->openCL.kernel("applyPositionDelta")->setArg(3, this->particles);
    
    // KERNEL ::  updatePositionAndVelocity
    this->openCL.loadKernel("updatePositionAndVelocity");
    this->openCL.kernel("updatePositionAndVelocity")->setArg(0, this->particles);
    this->openCL.kernel("updatePositionAndVelocity")->setArg(1, this->dt);
    
    // KERNEL :: applyViscosity
    /*
    this->openCL.loadKernel("applyViscosity");
    this->openCL.kernel("applyViscosity")->setArg(0, this->particles);
    this->openCL.kernel("applyViscosity")->setArg(1, this->sortedParticleToCell);
    this->openCL.kernel("applyViscosity")->setArg(2, this->gridCellOffsets);
    this->openCL.kernel("applyViscosity")->setArg(3, this->density);
    this->openCL.kernel("applyViscosity")->setArg(4, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("applyViscosity")->setArg(5, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("applyViscosity")->setArg(6, static_cast<int>(this->cellsPerAxis[2]));
    */
    
    // KERNEL :: computeCurl
    /*
    this->openCL.loadKernel("computeCurl");
    this->openCL.kernel("computeCurl")->setArg(0, this->particles);
    this->openCL.kernel("computeCurl")->setArg(1, this->sortedParticleToCell);
    this->openCL.kernel("computeCurl")->setArg(2, this->gridCellOffsets);
    this->openCL.kernel("computeCurl")->setArg(3, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("computeCurl")->setArg(4, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("computeCurl")->setArg(5, static_cast<int>(this->cellsPerAxis[2]));
    */
    
    // KERNEL :: applyVorticity
    /*
    this->openCL.loadKernel("applyVorticity");
    this->openCL.kernel("applyVorticity")->setArg(0, this->particles);
    this->openCL.kernel("applyVorticity")->setArg(1, this->dt);
    this->openCL.kernel("applyVorticity")->setArg(2, this->sortedParticleToCell);
    this->openCL.kernel("applyVorticity")->setArg(3, this->gridCellOffsets);
    this->openCL.kernel("applyVorticity")->setArg(4, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("applyVorticity")->setArg(5, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("applyVorticity")->setArg(6, static_cast<int>(this->cellsPerAxis[2]));
    */
    
    // KERNEL :: resolveCollisions
    this->openCL.loadKernel("resolveCollisions");
    this->openCL.kernel("resolveCollisions")->setArg(0, this->particles);
    this->openCL.kernel("resolveCollisions")->setArg(1, this->bounds.getMinExtent());
    this->openCL.kernel("resolveCollisions")->setArg(2, this->bounds.getMaxExtent());
    
}

/**
 * Resets various particle quantities, like density, etc.
 */
void Simulation::resetParticleQuantities()
{
    for (int i = 0; i < this->numParticles; i++) {
        
        Particle &p = this->particles[i];
        
        // Clear out the particle's predicted position:
        p.posStar.x = p.posStar.y = p.posStar.z = 0.0f;
        
        this->density[i]   = 0.0f;
        this->lambda[i]    = 0.0f;
        this->posDeltaX[i] = 0.0f;
        this->posDeltaY[i] = 0.0f;
        this->posDeltaZ[i] = 0.0f;
    }
}

/**
 * Initializes the simulation
 */
void Simulation::initialize()
{
    auto p1 = this->bounds.getMinExtent();
    auto p2 = this->bounds.getMaxExtent();
    
    // Dimension the OpenCL buffer to hold the given number of particles:

    this->particles.initBuffer(this->numParticles);
    
    // particleToCell contains [0 .. this->numParticles - 1] entries, where
    // each ParticlePosition instance (index is not important) maps
    // a particle's index (ParticlePosition#particleIndex) to a spatial grid
    // cell (ParticlePosition#cellI, ParticlePosition#cellJ, ParticlePosition#cellK),
    // where 0 <= ParticlePosition#cellI < cellsPerAxis.x,
    // 0 <= ParticlePosition#cellJ < cellsPerAxis.y, and
    // 0 <= ParticlePosition#cellK < cellsPerAxis.z

    this->particleToCell.initBuffer(this->numParticles);

    // Where the sorted version of the above will be sorted per simulation
    // step. The ParticlePosition indices will be sorted in ascending order
    // according to the linearized index computed from
    // (ParticlePosition#cellI, ParticlePosition#cellJ, ParticlePosition#cellK
    //
    // See the kernel helper function sub2ind in kernels/Simulation.cl for
    // details
    
    this->sortedParticleToCell.initBuffer(this->numParticles);

    // An array containing [0 .. this->numCells - 1] entries, where the
    // i-th entry contains the offset information about the start of a
    // particular grid cell in sortedParticleToCell. gridCellOffsets entries
    // are structs of type GridCellOffset, and are considered valid if
    // GridCellOffset#start != -1. This is used to speed up the lookup for
    // particles that happen to be in the same cell, so for instance, given
    // a grid cell offset at index i, g_i, all of the particles in cell
    // i are in the range sortedParticleToCell[g_i.start .. (g_i.start + g_i.length)]
    
    this->gridCellOffsets.initBuffer(this->numCells);

    // A histogram (count table), where the i-th entry contains the number of
    // particles occupying that linearized grid cell. For a linear grid cell
    // z, z can be computed from subscripts (i, j, k) by way of
    // z = i + (j * GRIDWIDTH) + (k * GRIDWIDTH * GRIDHEIGHT)

    this->cellHistogram.initBuffer(this->numCells);

    // The density values associated with each particle. The i-th density
    // corresponds to the i-th particle in this->particles:
    
    this->density.initBuffer(this->numParticles);
    this->lambda.initBuffer(this->numParticles);
    
    // For particle position correction in the solver:
    this->posDeltaX.initBuffer(this->numParticles);
    this->posDeltaY.initBuffer(this->numParticles);
    this->posDeltaZ.initBuffer(this->numParticles);
    
    // Set up initial positions and velocities for the particles:
    
    for (int i = 0; i < this->numParticles; i++) {

        Particle &p = this->particles[i];

        // Random position in the bounding box:
        p.pos.x = ofRandom(p1[0], p2[0]);
        p.pos.y = ofRandom(p1[1], p2[1]);
        p.pos.z = ofRandom(p1[2], p2[2]);
        
        // No predicted position:
        p.posStar.x = p.posStar.y = p.posStar.z = 0.0f;
        
        // All particles have uniform mass (for now):
        p.mass = this->particleMass;
        
        // and a uniform radius (for now):
        p.radius = this->particleRadius;
        
        // and no initial velocity:
        p.vel.x = p.vel.y = p.vel.z = 0.0f;
    }
    
    // Set the initial state:

    this->resetParticleQuantities();
    
    // Needed for particle grouping/"binning" by cell later:

    this->initializeParticleSort();

    // Dump the initial quantities assigned to the particles to the GPU, so we
    // can use them in GPU-land/OpenCL

    this->writeToGPU();

    // Load the kernels:

    this->initializeKernels();
}

/**
 * This implementation is based off of the method described in
 * "￼FAST FIXED-RADIUS NEAREST NEIGHBORS: INTERACTIVE MILLION-PARTICLE FLUID" by
 * Hoetzlein, 2014 that uses counting sort as an alternative to radix sort.
 *
 * See http://on-demand.gputechconf.com/gtc/2014/presentations/S4117-fast-fixed-radius-nearest-neighbor-gpu.pdf
 */
void Simulation::groupParticlesByCell()
{
    this->discretizeParticlePositions();
    this->sortParticlesByCell();
}

/**
 * Initializes the datastructures needed to group particles by cell. This
 * mostly involves zeroing out the structures and/or setting marker values
 * like -1 to indicate unset entries in the requisite arrays, etc.
 */
void Simulation::initializeParticleSort()
{
    // Zero out the histogram counts and grid cell offsets:
    for (int i = 0; i < this->numCells; i++) {

        this->cellHistogram[i] = 0;

        this->gridCellOffsets[i].start  = -1;
        this->gridCellOffsets[i].length = -1;
    }

    // as well as the particle to cell mappings:
    for (int i = 0; i < this->numParticles; i++) {

        // Particle index; -1 indicates unset
        this->particleToCell[i].particleIndex       = -1;
        this->particleToCell[i].cellI               = -1;
        this->particleToCell[i].cellJ               = -1;
        this->particleToCell[i].cellK               = -1;
        
        this->sortedParticleToCell[i].particleIndex = -1;
        this->sortedParticleToCell[i].cellI         = -1;
        this->sortedParticleToCell[i].cellJ         = -1;
        this->sortedParticleToCell[i].cellK         = -1;
    }
}

/**
 * Resets the state of the simulation
 */
void Simulation::reset()
{
    this->frameNumber = 0;
    
    this->initialize();
}

/**
 * Moves data from GPU buffers back to the host
 */
void Simulation::readFromGPU()
{
    this->particles.readFromDevice();
    this->particleToCell.readFromDevice();
    this->cellHistogram.readFromDevice();
    this->sortedParticleToCell.readFromDevice();
    this->gridCellOffsets.readFromDevice();
    this->density.readFromDevice();
    this->lambda.readFromDevice();
    this->posDeltaX.readFromDevice();
    this->posDeltaY.readFromDevice();
    this->posDeltaZ.readFromDevice();
}

/**
 * Writes data from the host to buffers on the GPU (i.e. the "device" in 
 * OpenCL parlance)
 */
void Simulation::writeToGPU()
{
    this->particles.writeToDevice();
    this->particleToCell.writeToDevice();
    this->cellHistogram.writeToDevice();
    this->sortedParticleToCell.writeToDevice();
    this->gridCellOffsets.writeToDevice();
    this->density.writeToDevice();
    this->lambda.writeToDevice();
    this->posDeltaX.writeToDevice();
    this->posDeltaY.writeToDevice();
    this->posDeltaZ.writeToDevice();
}

/**
 * Moves the state of the simulation forward one time step according to the
 * time step value, dt, passed to the constrcutor
 *
 * In this method, the motion of the particles, as well as the various
 * quatities assigned to them are updated, as described in the paper
 * "Position Based Fluids" by Miles Macklin & Matthias Muller.
 */
void Simulation::step()
{
    // Solver iterations (this will be adjustable later)

    int N = SOLVER_ITERATIONS;
    
    // We need to perform this step (zeroing out the histogram array and some
    // other data structures needed) before we compute the particle cell groups
    // because the zeroed structures will be written back to the GPU on the
    // next line, i.e. writeToGPU()

    this->initializeParticleSort();

    this->resetParticleQuantities();
    
    // Dump whatever changes occurred in host-land to the GPU:

    this->writeToGPU();

    // Where the actual work is done: the sequence of substeps follows
    // more-or-less from the listing "Algorithm 1 Simulation Loop" in the
    // paper "Position Based Fluids". The main difference is that we are using
    // a different method than Macklin and Muller to compute the nearest
    // neighbors of a given particle. Whereas they use the method by
    // [Green 2008], we use the method described by Hoetzlein, 2014
    // in the slides
    // "￼FAST FIXED-RADIUS NEAREST NEIGHBORS: INTERACTIVE MILLION-PARTICLE FLUID"
    // that uses counting sort as an alternative to radix sort
    
    this->applyExternalForces();

    this->predictPositions();

    this->groupParticlesByCell();

    // Solver runs for N iterations:
    for (int i = 0; i < N; i++) {

        this->calculateDensity();
        
        this->calculatePositionDelta();
        
        this->applyPositionDelta();

        this->handleCollisions();
    }
    
    //this->applyXSPHViscosity();
    
    //this->applyVorticityConfinement();

    this->updatePositionAndVelocity();
    
    // Make sure the OpenCL work queue is empty before proceeding. This will
    // block until all the stuff in GPU-land is done before moving forward
    // and reading the results of the work we did on the GPU back into
    //host-land:
    
    this->openCL.finish();
    
    // Read the changes back from the GPU so we can manipulate the values
    // in our C++ program:

    this->readFromGPU();

    // Finally, bump up the frame counter:

    this->frameNumber++;
}

/**
 * Draws the cell grid
 */
void Simulation::drawGrid() const
{
    auto p1 = this->bounds.getMinExtent();
    auto p2 = this->bounds.getMaxExtent();
    
    float xCellWidth = (p2[0] - p1[0]) / static_cast<float>(this->cellsPerAxis[0]);
    float halfXWidth = xCellWidth * 0.5f;
    float yCellWidth = (p2[1] - p1[1]) / static_cast<float>(this->cellsPerAxis[1]);
    float halfYWidth = yCellWidth * 0.5f;
    float zCellWidth = (p2[2] - p1[2]) / static_cast<float>(this->cellsPerAxis[2]);
    float halfZWidth = zCellWidth * 0.5f;
    
    ofNoFill();
    ofSetColor(0, 255, 0);

    for (int i = 1; i < (2  * this->cellsPerAxis[0]); i += 2) {

        float xCorner = p1[0] + (static_cast<float>(i) * halfXWidth);

        for (int j = 1; j < (2  * this->cellsPerAxis[1]); j += 2) {
        
            float yCorner = p1[1] + (static_cast<float>(j) * halfYWidth);
        
            for (int k = 1; k < (2  * this->cellsPerAxis[2]); k += 2) {
                
                float zCorner = p1[2] + (static_cast<float>(k) * halfZWidth);

                ofDrawBox(xCorner, yCorner, zCorner, xCellWidth, yCellWidth, zCellWidth);
            }
        }
    }
}

/**
 * Draws the bounds of the simulated environment as a transparent with solid
 * lines indicating the edges of the bounding box.
 */
void Simulation::drawBounds() const
{
    // Draw the bounding box that will hold the particles:
    auto p1 = this->bounds.getMinExtent();
    auto p2 = this->bounds.getMaxExtent();
    auto x  = (p1[0] + p2[0]) * 0.5f;
    auto y  = (p1[1] + p2[1]) * 0.5f;
    auto z  = (p1[2] + p2[2]) * 0.5f;
    auto w  = p2[0] - p1[0];
    auto h  = p2[1] - p1[1];
    auto d  = p2[2] - p1[2];
    
    ofNoFill();
    ofSetColor(255, 255, 255);
    ofDrawBox(x, y, z, w, h, d);
}

/**
 * Currently, draws the positions of the particles using a fixed color. 
 * Later, this may be changed so that the color of the particle reflects
 * some quantity like velocity, mass, viscosity, etc.
 */
void Simulation::drawParticles()
{
    for (int i = 0; i < this->numParticles; i++) {
        
        Particle &p = this->particles[i];

        // The fill:
        ofSetColor(51, 153, 255);
        ofFill();
        ofDrawSphere(p.pos.x, p.pos.y, p.pos.z, p.radius);
        
        // Draw an outline of the particle:
        ofSetColor(0, 0, 255);
        ofNoFill();
        ofDrawSphere(p.pos.x, p.pos.y, p.pos.z, p.radius);
        
        if (this->isVisualDebuggingEnabled()) {

            // Label the particle with its number:
            ofSetColor(255, 255, 0);
            ofFill();
            ofPushMatrix();
                ofTranslate(0,0,p.pos.z);
                ofDrawBitmapString(ofToString(i), p.pos.x, p.pos.y);
            ofPopMatrix();
        }
    }
}

/**
 * This method is once per step of the simulation to render all graphical
 * output, including rendering the bounding box of the simulated environment,
 * all particles in the simulation, as well as any additional objects (meshs,
 * walls, etc.) that may exist.
 */
void Simulation::draw()
{
    //ofClear(0, 0, 0);
    this->drawBounds();

    if (this->drawGridEnabled()) {
        this->drawGrid();
    }

    this->drawParticles();
    
}

/******************************************************************************/

/**
 * Applies external forces to the particles in the simulation.
 *
 * @see kernels/Simulation.cl (applyExternalForces) for details
 */
void Simulation::applyExternalForces()
{
    this->openCL.kernel("applyExternalForces")->run1D(this->numParticles);
}

/**
 * Updates the predicted positions of the particles via an explicit Euler step
 *
 * @see kernels/Simulation.cl (predictPosition) for details
 */
void Simulation::predictPositions()
{
    this->openCL.kernel("predictPosition")->run1D(this->numParticles);
}

/**
 * Discretizes all of the particles to a grid cell, where the number of
 * grid cells along each axis in the simulated space is specified by 
 * cellsPerAxis, e.g. (4,5,6) specifies 4 cells in the x-axs, 5 in the y-axis, 
 * and 6 in the z-axis
 *
 * @see kernels/Simulation.cl (discretizeParticlePositions) for details
 */
void Simulation::discretizeParticlePositions()
{
    this->openCL.kernel("discretizeParticlePositions")->run1D(this->numParticles);
}

/**
 * Sorts the particles by the assigned grid cell. Following the run of this
 * function, sortedParticleToCell (after read back fro the GPU) will contain a 
 * listing of ParticlePosition, that are sorted by linearized cell indices, e.g.
 * particles that are in the same cell will be consecutive in 
 * sortedParticleToCell, making neighbor search quick.
 *
 * @see kernels/Simulation.cl (sortParticlesByCell) for details
 */
void Simulation::sortParticlesByCell()
{
    // Only 1 thread is needed to run this. The sorting operation is
    // sequential in nature, hence the invocation with 1 thread, e.g.
    // "kernel("sortParticlesByCell")->run1D(1)"!

    this->openCL.kernel("sortParticlesByCell")->run1D(1);
}

/**
 * Computes the density for each particle using the SPH density estimator
 * 
 * (*) Specifically, this function is part of the constraint solver loop
 *
 * @see kernels/Simulation.cl (estimateDensity) for details
 */
void Simulation::calculateDensity()
{
    this->openCL.kernel("estimateDensity")->run1D(this->numParticles);
}

/**
 * Computes the position delta
 *
 * (*) Specifically, this function is part of the constraint solver loop
 *
 * @see kernels/Simulation.cl (computeLambda, computePositionDelta) for details
 */
void Simulation::calculatePositionDelta()
{
    this->openCL.kernel("computeLambda")->run1D(this->numParticles);
    this->openCL.kernel("computePositionDelta")->run1D(this->numParticles);
}

/**
 * Apply the position delta
 *
 * @see kernels/Simulation.cl (applyPositionDelta) for details
 */
void Simulation::applyPositionDelta()
{
    this->openCL.kernel("applyPositionDelta")->run1D(this->numParticles);
}

/**
 * TODO: All this does now is clamp the particle positions to the simulation
 * bounding box
 *
 * @see kernels/Simulation.cl (resolveCollisions) for details
 */
void Simulation::handleCollisions()
{
    this->openCL.kernel("resolveCollisions")->run1D(this->numParticles);
}

/**
 * [TODO]
 *
 * Apply XSPH viscosity
 *
 *  f_i_viscosity = mu SUM_j(mass_j * ((v_j - v_i)/density_j) * Nabla^2 W(|x_i - x_j|)
 *
 *  where mu = 0.05f and Nabla^2 = (Laplace operator)
 *  particle i looks at neighbors j
 *
 */
void Simulation::applyXSPHViscosity()
{
    //this->openCL.loadKernel("applyViscosity")->run1D(this->numParticles);
}

/**
 * [TODO]
 *
 * Apply vorticity confinement
 */
void Simulation::applyVorticityConfinement()
{
    //this->openCL.loadKernel("computeCurl")->run1D(this->numParticles);
    //this->openCL.loadKernel("applyVorticity")->run1D(this->numParticles);
}

/**
 * Updates the actualm, final position and velocity of the particles
 * in the current simulation step
 *
 * @see kernels/Simulation.cl (updatePositionAndVelocity) for details
 */
void Simulation::updatePositionAndVelocity()
{
    this->openCL.kernel("updatePositionAndVelocity")->run1D(this->numParticles);
}

/******************************************************************************/
