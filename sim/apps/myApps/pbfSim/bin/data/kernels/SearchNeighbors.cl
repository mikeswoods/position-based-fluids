/*******************************************************************************
 * SearchNeighbors.cl
 * - TODO
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

/**
 * A helper function that scales a value x in the range [a0,a1] to a new
 * range [b0,b1]
 */
float rescale(float x, float a0, float a1, float b0, float b1)
{
    return ((x - a0) / (a1 - a0)) * (b1 - b0) + b0;
}

/**
 * [HELPER]
 *
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

/**
 * [KERNEL] Discretizes particles based on their positions to a cell in
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
kernel void discretize(global Particle* particles
                      ,global int2* particleToCell
                      ,int3 cellsPerAxis
                      ,float3 minExtent
                      ,float3 maxExtent)
{
    int i = get_global_id(0);
    global Particle *p = &particles[i];
    
    // Now we have the discretized cell at (i, j, k):
    int iCell = (int)(rescale(p->pos.x, minExtent.x, maxExtent.x, 0.0, (float)(cellsPerAxis.x - 1)));
    int jCell = (int)(rescale(p->pos.y, minExtent.y, maxExtent.y, 0.0, (float)(cellsPerAxis.y - 1)));
    int kCell = (int)(rescale(p->pos.z, minExtent.z, maxExtent.z, 0.0, (float)(cellsPerAxis.z - 1)));
    
    // Convert the cell subscript (i,j,k) to a linear index z:
    int z = sub2ind(iCell, jCell, kCell, cellsPerAxis.x, cellsPerAxis.y);
    
    // Store the particle's index i, and the linearized and
    // discretized cell index z
    particleToCell[i] = (int2)(i, z);
}
