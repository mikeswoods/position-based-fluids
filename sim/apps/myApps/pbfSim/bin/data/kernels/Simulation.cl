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

/*******************************************************************************
 * Types
 ******************************************************************************/

typedef struct {
    
    float4 pos;    // 4 words
    
    float4 vel;    // 4 words
    
    float  mass;   // 1 word

    float  radius; // 1 word

    /**
     * VERY IMPORTANT: This is needed so that the struct's size is aligned
     * for x86 memory access along 4/word 16 byte intervals.
     *
     * If the size is not aligned, results WILL be screwed up!!!
     * Don't be like me and waste hours trying to debug this issue. The
     * OpenCL compiler WILL NOT pad your struct to so that boundary aligned
     * like g++/clang will in the C++ world.
     *
     * See http://en.wikipedia.org/wiki/Data_structure_alignment
     */
    float  __dummy[2]; // 2 words

} Particle; // total = 12 words = 64 bytes

typedef struct {
    int particleIndex; // Index of particle in particle buffer
    int cellI;         // Corresponding grid index in the x-axis
    int cellJ;         // Corresponding grid index in the y-axis
    int cellK;         // Corresponding grid index in the z-axis
} ParticlePosition;

/*******************************************************************************
 * Constants
 ******************************************************************************/

/**
 * Acceleration due to gravity: 9.8 m/s
 */
#define G 9.8f

/*******************************************************************************
 * Helper functions
 ******************************************************************************/

/**
 * A helper function that scales a value x in the range [a0,a1] to a new
 * range [b0,b1]
 */
float rescale(float x, float a0, float a1, float b0, float b1)
{
    return ((x - a0) / (a1 - a0)) * (b1 - b0) + b0;
}

/**
 * A function that converts a 3D subscript (i,j,k) into a linear index
 *
 * @param [in] int i x component of subscript
 * @param [in] int j y component of subscript
 * @param [in] int k z component of subscript
 * @param [in] int w grid width
 * @param [in] int h grid height
 */
int sub2ind(int i, int j, int k, int w, int h)
{
    return i + (j * w) + k * (w * h);
}

/**
 * [HELPER]
 *
 * A function that converts a linear index x into a 3D subscript (i,j,k)
 *
 * @param [in] int x The linear index x
 * @param [in] int w grid width
 * @param [in] int h grid height
 */
int3 ind2sub(int x, int w, int h)
{
    return (int3)(x % w, (x / w) % h, x / (w * h));
}

/*******************************************************************************
 * Kernels
 ******************************************************************************/

/**
 * Currently, only applies gravity to the y component of the velocity.
 * Additional forces may be added later like wind and other forms of
 * turbulence, etc.
 *
 *   v_i = v_i + dt + f_external(x_i)
 */
kernel void applyExternalForces(global Particle* particles, float dt)
{
    int i = get_global_id(0);
    
    // Apply the force of gravity along the y-axis:
    particles[i].vel.y += (dt * -G);
}

/**
 * Updates the predicted position of the particle using an explicit Euler step
 *
 *   x_i = x_i + (dt * v_i)
 */
kernel void predictPosition(global Particle* particles, float dt)
{
    int i = get_global_id(0);

    // Explicit Euler step:
    particles[i].pos += (dt * particles[i].vel);
    
}

/**
 * Discretizes particles based on their positions to a cell in
 * a grid consisting of (cellsPerAxis[0], cellsPerAxis[1], cellsPerAxis[2])
 * along the (x, y, z) axes.
 *
 * @param [in] Particle* particles The particles to assign to cells
 * @param [out] int2* particleToCell Each entry contains a int2 pair
 * (i,j), where i is the particle in the i-th entry of particles, and j is
 * the linear index of the corresponding linear bin (j_x, j_y, j_z), where
 * 0 <= j_x < cellsPerAxis.x, 0 <= j_y < cellsPerAxis.y,
 * and 0 <= j_z < cellsPerAxis.z
 * @param [in] int3 cellsPerAxis The number of cells per axis in the grid
 * @param [in] float3 minExtent The minimum extent of the simulation's
 *             bounding box in world space
 * @param [in] float3 maxExtent The maximum extent of the simulation's
 *             bounding box in world space
 */
kernel void discretizeParticlePositions(global Particle* particles
                                       ,global ParticlePosition* particleToCell
                                       ,int3 cellsPerAxis
                                       ,float3 minExtent
                                       ,float3 maxExtent)
{
    int i = get_global_id(0);
    global Particle *p = &particles[i];
    
    // Now we have the discretized cell at (i, j, k):
    particleToCell[i].cellI = (int)(rescale(p->pos.x, minExtent.x, maxExtent.x, 0.0, (float)(cellsPerAxis.x - 1)));
    particleToCell[i].cellJ = (int)(rescale(p->pos.y, minExtent.y, maxExtent.y, 0.0, (float)(cellsPerAxis.y - 1)));
    particleToCell[i].cellK = (int)(rescale(p->pos.z, minExtent.z, maxExtent.z, 0.0, (float)(cellsPerAxis.z - 1)));
}

/**
 *
 *
 */
kernel void computeCellHistogram(global ParticlePosition* particleToCell
                                ,global int* cellHistogram
                                ,int3 cellsPerAxis)
{
    int i = get_global_id(0);
    int cellI = particleToCell[i].cellI;
    int cellJ = particleToCell[i].cellJ;
    int cellK = particleToCell[i].cellK;
    
    cellHistogram[sub2ind(cellI, cellJ, cellK, cellsPerAxis.x, cellsPerAxis.y)]++;
}

/**
 *
 *
 */
/*
kernel void sortParticlesByCell(global ParticlePosition* particleToCell
                               ,global int* cellHistogram
                               ,global ParticlePosition* sortedParticleToCell
                               ,int n)
{
    
}
 */
