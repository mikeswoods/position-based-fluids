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
 * @param [in] _numParticles The number of particles in the simulation
 * @param [in] _massPerParticle The mass per particle
 */
Simulation::Simulation(msa::OpenCL& _openCL
                      ,AABB _bounds
                      ,float _dt
                      ,int _numParticles
                      ,float _massPerParticle) :
    openCL(_openCL),
    bounds(_bounds),
    dt(_dt),
    cellsPerAxis(EigenVector3(4, 4, 4 )),
    numParticles(_numParticles),
    massPerParticle(_massPerParticle),
    frameNumber(0)
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
    
    //
    this->openCL.loadKernel("applyExternalForces");
    this->openCL.kernel("applyExternalForces")->setArg(0, this->particles);
    this->openCL.kernel("applyExternalForces")->setArg(1, this->dt);
    
    //
    this->openCL.loadKernel("predictPosition");
    this->openCL.kernel("predictPosition")->setArg(0, this->particles);
    this->openCL.kernel("predictPosition")->setArg(1, this->dt);
    
    //
    this->openCL.loadKernel("discretizeParticlePositions");
    this->openCL.kernel("discretizeParticlePositions")->setArg(0, this->particles);
    this->openCL.kernel("discretizeParticlePositions")->setArg(1, this->particleToCell);
    this->openCL.kernel("discretizeParticlePositions")->setArg(2, this->cellHistogram);
    this->openCL.kernel("discretizeParticlePositions")->setArg(3, static_cast<int>(this->cellsPerAxis[0]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(4, static_cast<int>(this->cellsPerAxis[1]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(5, static_cast<int>(this->cellsPerAxis[2]));
    this->openCL.kernel("discretizeParticlePositions")->setArg(6, this->bounds.getMinExtent());
    this->openCL.kernel("discretizeParticlePositions")->setArg(7, this->bounds.getMaxExtent());
    
    //
    this->openCL.loadKernel("sortParticlesByCell");
    this->openCL.kernel("sortParticlesByCell")->setArg(0, this->particleToCell);
    this->openCL.kernel("sortParticlesByCell")->setArg(1, this->cellHistogram);
    this->openCL.kernel("sortParticlesByCell")->setArg(2, this->sortedParticleToCell);
    this->openCL.kernel("sortParticlesByCell")->setArg(3, this->numParticles);
    this->openCL.kernel("sortParticlesByCell")->setArg(4, this->numCells);
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
    
    // A mapping of particle indices to linearized cell indices. Each entry
    // is a pair (i,j), where i refers to a particle in this->particles,
    // and j refers to a linear index in the range [0, N], where
    // N = (cellsPerAxis[0] * _cellsPerAxis[1] * _cellsPerAxis[2])

    this->particleToCell.initBuffer(this->numParticles);
    this->sortedParticleToCell.initBuffer(this->numParticles);

    // A histogram (count table), where the i-th entry contains the number of
    // particles occupying that linearized grid cell. For a linear grid cell
    // z, z can be computed from subscripts (i, j, k) by way of
    // z = i + (j * GRIDWIDTH) + (k * GRIDWIDTH * GRIDHEIGHT)

    this->cellHistogram.initBuffer(this->numCells);

    for (int i = 0; i < this->numParticles; i++) {

        Particle &p = this->particles[i];

        // Random position in the bounding box:
        p.pos.x = ofRandom(p1[0], p2[0]);
        p.pos.y = ofRandom(p1[1], p2[1]);
        p.pos.z = ofRandom(p1[2], p2[2]);
        
        // All particles have uniform mass (for now):
        p.mass = this->massPerParticle;
        
        // and a uniform radius (for now):
        p.radius = 0.1f;
        
        // and no initial velocity:
        p.vel.x = p.vel.y = p.vel.z = 0.0f;
    }
    
    // Needed for particle to cell binning later:
    this->initializeParticleSort();
    
    // Dump the initial quantities assigned to the particles to the GPU, so we
    // can use them in GPU-land/OpenCL

    this->writeToGPU();
    
    // Load the kernels:

    this->initializeKernels();
}

/**
 *
 */
void Simulation::initializeParticleSort()
{
    // Zero out the histogram counts:
    for (int i = 0; i < this->numCells; i++) {
        this->cellHistogram[i] = 0;
    }
    
    // as well as the particle to cell mappings:
    for (int i = 0; i < this->numParticles; i++) {
        this->particleToCell[i].particleIndex       = -1; // Particle index; -1 indicates unset
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
}

/**
 * Writes data from the host to buffers on the GPU (device
 */
void Simulation::writeToGPU()
{
    this->particles.writeToDevice();
    this->particleToCell.writeToDevice();
    this->cellHistogram.writeToDevice();
    this->sortedParticleToCell.writeToDevice();
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
    this->initializeParticleSort();
    
    // Dump whatever changes occurred in host-land to the GPU
    this->writeToGPU();

    this->applyExternalForces();
    this->predictPositions();
    this->discretizeParticlePositions();
    this->sortParticlesByCell();

    this->openCL.finish();
    
    // Read the changes back from the GPU so we can manipulate the values
    // in our C++ program:
    this->readFromGPU();

    // Finally, bump up the frame counter;
    this->frameNumber++;
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
    ofSetColor(0, 0, 255);
    ofFill();

    for (int i = 0; i < this->numParticles; i++) {
        Particle &p = this->particles[i];
        ofDrawSphere(p.pos.x, p.pos.y, p.pos.z, p.radius);
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
    this->drawParticles();
}

/******************************************************************************/

/**
 * Applies external forces to the particles in the simulation.
 *
 * @see kernels/UpdatePositions.cl for details
 */
void Simulation::applyExternalForces()
{
    this->openCL.kernel("applyExternalForces")->run1D(this->numParticles);
}

/**
 * Updates the predicted positions of the particles via an explicit Euler step
 *
 * @see kernels/UpdatePositions.cl for details
 */
void Simulation::predictPositions()
{
    this->openCL.kernel("predictPosition")->run1D(this->numParticles);
}

/**
 *
 */
void Simulation::discretizeParticlePositions()
{
    this->openCL.kernel("discretizeParticlePositions")->run1D(this->numParticles);
}

/**
 *
 */
void Simulation::sortParticlesByCell()
{
    this->openCL.kernel("sortParticlesByCell")->run1D(1);
}

/******************************************************************************/
