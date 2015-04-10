/*******************************************************************************
 * UpdatePositions.cl
 * - The OpenCL kernel responsible for apply external forces, like gravity
 *   for instance to each particle in the simulation, and subsequenly updating
 *   the predicted position of each particle using a simple explicit Euler
 *   step
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

/**
 * This needs to be re-defined here, so we have an inline definition of a
 * Particle. Since the original definition sits in Simulation.h, it's in
 * host-land, and we're in GPU-world. We can pass data from host <-> GPU,
 * but not type definitions (as far as I know).
 */
typedef struct {
    float3 pos;
    float3 vel;
    float  mass;
    float  radius;
} Particle;

/**
 * Acceleration due to gravity: 9.8 m/s
 */
#define G 9.8f

/**
 * Currently, only applies gravity to the y component of the velocity.
 * Additional forces may be added later like wind and other forms of
 * turbulence, etc.
 *
 *   v_i = v_i + dt + f_external(x_i)
 */
__kernel void applyExternalForces(__global Particle* particles, float dt)
{
    int id = get_global_id(0);
    
    // Apply the force of gravity along the y-axis:
    particles[id].vel.y += (dt * -G);
}

/**
 * Updates the predicted position of the particle using an explicit Euler step
 *
 *   x_i = x_i + (dt * v_i)
 */
__kernel void predictPosition(__global Particle* particles, float dt)
{
    int id = get_global_id(0);

    // Explicit Euler step:
    //particles[id].pos += (dt * particles[id].vel);
    
    particles[id].pos.y += -0.5f;
}
