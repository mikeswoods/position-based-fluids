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

Simulation::Simulation(AABB _bounds, unsigned int _numParticles, float _massPerParticle) :
    bounds(_bounds),
    numParticles(_numParticles),
    massPerParticle(_massPerParticle)
{
    this->initialize();
}

Simulation::~Simulation()
{
    
}

void Simulation::initialize()
{
    // We multiply by 3, since we need to store the (x,y,z) components for each
    unsigned int K = 3 * this->numParticles;

    // Dimension all of the vectors / matrices based on the number of
    // particles to be simulated:
    this->particleX.resize(K);
    this->particleV.resize(K);
    this->particleMass.resize(K, K);
    this->invParticleMass.resize(K, K);
    
    // Set the initial positions of the particles:
    // particleX ...
    
    
    // All particle velocities are initially zero:
    this->particleV.setZero();
}

void Simulation::reset()
{
    
}

void Simulation::update()
{
    float dt = 0.01;
    
    // Apply external forces to the particles, like gravity:
    this->computeExternalForce();

    // Perform PDB integration over all particles:
    this->integratePBD(this->particleX, this->particleV, dt, 10);

    //vector<CollisionInstance> collisions;
    //detectCollision(x, collisions);
    //resolveCollision(x, v, collisions);
}

void Simulation::draw()
{
    // Draw the bounding box that will hold the particles:
    auto p1 = this->bounds.getMinExtent();
    auto p2 = this->bounds.getMaxExtent();
    auto w  = p2[0] - p1[0];
    auto h  = p2[1] - p1[1];
    auto d  = p2[2] - p1[2];

    ofNoFill();
    ofSetColor(255, 255, 255);
    ofDrawBox(p1[0], p1[1], p1[2], w, h, d);
}

/******************************************************************************/

void Simulation::computeExternalForce()
{
    
}

void Simulation::integratePBD(VectorX& x, VectorX& v, float dt, unsigned int ns)
{
    /*
    // Explicit Euler step, with only gravity acting as an external force:
    v = v + (dt * (m_mesh->m_inv_mass_matrix * this->m_external_force));
    VectorX p = x + (dt * v);
    
    for (unsigned int i=0; i<this->m_iterations_per_frame; i++) {
        if ((i % 2) == 0) {
            // Forward
            for (auto k = this->m_constraints.begin(); k != this->m_constraints.end(); k++) {
                (*k)->PBDProject(p, this->m_mesh->m_inv_mass_matrix, ns);
            }
        } else {
            // Backward
            for (auto k = this->m_constraints.rbegin(); k != this->m_constraints.rend(); k++) {
                (*k)->PBDProject(p, this->m_mesh->m_inv_mass_matrix, ns);
            }
        }
    }
    
    v = (p - x) / dt;
    x = p;
    */
}

/******************************************************************************/
