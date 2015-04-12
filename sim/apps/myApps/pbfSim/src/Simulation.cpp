/*******************************************************************************
 * Simulation.cpp
 * - The heart of the position-based fluids simulator. This class encapsulates
 *   the current state of the simulation
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include "ofMain.h"
#include "ofx3DModelLoader.h"
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
                      ,unsigned int _numParticles
                      ,float _massPerParticle) :
    openCL(_openCL),
    bounds(_bounds),
    dt(_dt),
    numParticles(_numParticles),
    massPerParticle(_massPerParticle),
    frameNumber(0)
{
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
    
    // Read the source:
    this->openCL.loadProgramFromFile("kernels/UpdatePositions.cl");
    
    // Initialize the specified kernels and bind the parameters that are
    // the same across invocations of the kernel:

    // 0) applyExternalForces
    this->openCL.loadKernel("applyExternalForces");
    this->openCL.kernel("applyExternalForces")->setArg(0, this->particles);
    this->openCL.kernel("applyExternalForces")->setArg(1, this->dt);
    
    // 1) predictPosition
    this->openCL.loadKernel("predictPosition");
    this->openCL.kernel("predictPosition")->setArg(0, this->particles);
    this->openCL.kernel("predictPosition")->setArg(1, this->dt);

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
    
    // Dump the initial quantities assigned to the particles to the GPU, so we
    // can use them in GPU-land/OpenCL
    this->particles.writeToDevice();
    
    // Load the kernels:
    this->initializeKernels();
}

/**
 * Resets the state of the simulation
 */
void Simulation::reset()
{
    
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
    // Dump whatever changes occurred in host-land to the GPU
    this->particles.writeToDevice();
 
    //cout << "STEP " << this->frameNumber << " " << this->particles[0] << endl;
    
    // Apply external forces to the particles, like gravity:
    this->applyExternalForces();
 
    // Predict the next position of the particles:
    this->predictPositions();

    // Make sure no more work remains in the OpenCL work queue. This will
    // block until all OpenCL related stuff in the step() function has
    // finished running:
    this->openCL.finish();
    
    // Read the changes back from the GPU so we can manipulate the values
    // in our C++ program:
    this->particles.readFromDevice();

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

/**
 * This method sorts the buckets by 
 *
 */
void countingSort(int arr[], int sz)
{
    int i, j, k, min, max, idx = 0;
    
    min = max = arr[0];
    
    for(i = 1; i < sz; i++)
    {
        min = (arr[i] < min) ? arr[i] : min;
        max = (arr[i] > max) ? arr[i] : max;
    }
    k = max - min + 1; /* creates k buckets */
    int *B = new int [k];
    
    for(i = 0; i < k; i++)
        B[i] = 0;
    
    for(i = 0; i < sz; i++)
        B[arr[i] - min]++;
    
    for(i = min; i <= max; i++)
        for(j = 0; j < B[i - min]; j++)
            arr[idx++] = i;
    
    delete [] B;

}

/*

#define H               1.5f  // smoothing radius
#define H_9             (H*H*H*H*H*H*H*H*H) // h^9
#define H_6             (H*H*H*H*H*H) // h^6

float poly6Kernel(Vector3DS p_i, Vector3DS p_j){
    Vector3DS diff = p_i - p_j;
    
    float r = diff.length();
    if (H > r && r > 0) {
        float h_minus_r = (H * H - r * r);
        float div = 64.0 * M_PI * H_9 * h_minus_r * h_minus_r * h_minus_r;
        return 315.0f / div;
    }
    return 0;
}

float spikyKernel(Vector3DS p_i, Vector3DS p_j){
    Vector3DS diff = p_i - p_j;
    float r = diff.length();
    if (H > r && r > 0) {
        float h_minus_r = H - r;
        float div = M_PI * H_6 * h_minus_r * h_minus_r * h_minus_r;
        return 15.0f / div;
    }
    return 0;
}
*/

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

/******************************************************************************/
