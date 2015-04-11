/*******************************************************************************
 * Simulation.cpp
 * - The heart of the position-based fluids simulator. This class encapsulates
 *   the current state of the simulation
 *
 * CIS563: Physcially Based Animation final project
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
    massPerParticle(_massPerParticle)
{
    this->initialize();
}

Simulation::~Simulation()
{
    
}

/**
 * Loads all of the OpenCL kernels that will be used for during the simulation
 */
void Simulation::loadKernels()
{
    
    this->openCL.loadProgramFromFile("kernels/UpdatePositions.cl");
    this->openCL.loadKernel("updatePositions");
    
    this->openCL.kernel("updateParticle")->setArg(0, particles);

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
        p.x.x = ofRandom(p1[0], p2[0]);
        p.x.y = ofRandom(p1[1], p2[1]);
        p.x.z = ofRandom(p1[2], p2[2]);
        
        // All particles have uniform mass:
        p.mass = 1.0f;
        
        // and no initial velocity:
        p.v.x = p.v.y = p.v.z = 0.0f;
    }
    
    // Dump the particles to the GPU:
    this->particles.writeToDevice();
    
    // Load the kernels:
    this->loadKernels();
}

/**
 * Resets the state of the simulation
 */
void Simulation::reset()
{
    
}

/**
 * Runs once per step of the simulation. The update() method is run before the
 * accompanying draw method to change the state of the simulation.
 *
 * In this method, the motion of the particles, as well as the various
 * quatities assigned to them are updated, as described in the paper
 * "Position Based Fluids" by Miles Macklin & Matthias Muller.
 */
void Simulation::update()
{
    // Apply external forces to the particles, like gravity:
    this->computeExternalForce();
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

        ofSetColor(0, 0, 255);
        ofFill();
        ofDrawSphere(p.x.x, p.x.y, p.x.z, 0.1f);
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
    this->drawBounds();
    this->drawParticles();
    
    glColor3f(1, 1, 1);
    string info = "fps: " + ofToString(ofGetFrameRate()) + "\nnumber of particles: " + ofToString(this->numParticles);
    ofDrawBitmapString(info, 20, 20);

}

/**
 * This method sorts the buckets
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

/******************************************************************************/

/**
 * Applies external forces to the
 */
void Simulation::computeExternalForce()
{
    
}

/******************************************************************************/
