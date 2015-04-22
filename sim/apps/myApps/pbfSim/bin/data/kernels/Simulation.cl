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
 * Preprocessor directives
 ******************************************************************************/

//#define USE_BRUTE_FORCE_SEARCH 1

/*******************************************************************************
 * Constants
 ******************************************************************************/

/**
 * A small epislon value
 */
const constant float EPSILON = 1.0e-4f;

/**
 * The maximum number of neighbors to examine for a given particle:
 */
const constant int CHECK_MAX_NEIGHBORS = 8;

/**
 * Acceleration force due to gravity: 9.8 m/s
 */
const constant float G = 9.8f;

/**
 * Particle rest density: 1000kg/m^3 = 10,000g/m^3
 */
const constant float REST_DENSITY     = 10000.0f;
const constant float INV_REST_DENSITY = 1.0f / REST_DENSITY;

/*******************************************************************************
 * Types
 ******************************************************************************/

// Tuneable parameters for the simulation:

typedef struct {
    
    float particleRadius;      // 1. Particle radius

    float smoothingRadius;     // 2. Kernel smoothing radius
    
    float relaxation;          // 3. Pressure relaxation coefficient (epsilon)
    
    float artificialPressureK; // 4. Artificial pressure coefficient K
    
    float artificialPressureN; // 5. Artificial pressure coefficient N
    
    float epsilonVorticity;    // 6. Vorticity coefficient
    
    float viscosityCoeff;      // 7. Viscosity coefficient
    
    float __padding[1];
    
} Parameters;

// A particle type:

typedef struct {
    
    float4 pos;     // Current particle position (x), 4 words
    
    float4 posStar; // Predicted particle position (x*), 4 words
    
    float4 vel;     // Current particle velocity (v), 4 words
    
    float4 velDelta; // Velocity delta from XSPH viscosity (v*), 4 words

    /**
     * VERY IMPORTANT: This is needed so that the struct's size is aligned
     * for x86 memory access along 4/word 16 byte intervals.
     *
     * If the size is not aligned, results WILL be screwed up!!!
     * Don't be like me and waste hours trying to debug this issue. The
     * OpenCL compiler WILL NOT pad your struct to so that boundary aligned
     * like g++/clang will in host (C++) land!!!.
     *
     * See http://en.wikipedia.org/wiki/Data_structure_alignment
     */
    //float  __padding[2]; // 2 words

} Particle;

// A type to represent the position of a given particle in the spatial
// grid the simulated world is divided into

typedef struct {

    int particleIndex; // Index of particle in particle buffer (1 word)

    int cellI;         // Corresponding grid index in the x-axis (1 word)
    
    int cellJ;         // Corresponding grid index in the y-axis (1 word)
    
    int cellK;         // Corresponding grid index in the z-axis (1 word)

    int key;           // Linearized index key computed from the subscript
                       // (cellI, cellJ, cellK)
    int __padding[3];
    
} ParticlePosition;

// A type that encodes the start and length of a grid cell in sortedParticleToCell

typedef struct {
    
    int  start; // Start of the grid cell in sortedParticleToCell
    
    int length;
    
    int __padding[2]; // Padding
    
} GridCellOffset;

/*******************************************************************************
 * Forward declarations
 ******************************************************************************/

global float rescale(float x, float a0, float a1, float b0, float b1);

global int sub2ind(int i, int j, int k, int w, int h);

global int3 ind2sub(int x, int w, int h);

global int getKey(const global Particle* p
                 ,int cellsX
                 ,int cellsY
                 ,int cellsZ
                 ,float3 minExtent
                 ,float3 maxExtent);

global int3 getSubscript(const global Particle* p
                        ,int cellsX
                        ,int cellsY
                        ,int cellsZ
                        ,float3 minExtent
                        ,float3 maxExtent);

float poly6(float4 r, float h);

float poly6_scalar(float q, float h);

float4 spiky(float4 r, float h);

void callback_SPHDensityEstimator_i(const global Parameters* parameters
                                   ,int i
                                   ,const global Particle* p_i
                                   ,int j
                                   ,const global Particle* p_j
                                   ,void* dataArray
                                   ,void* accum);

 void callback_SPHGradient_i(const global Parameters* parameters
                            ,int i
                            ,const global Particle* p_i
                            ,int j
                            ,const global Particle* p_j
                            ,void* dataArray
                            ,void* accum);

 void callback_SquaredSPHGradientLength_j(const global Parameters* parameters
                                         ,int i
                                         ,const global Particle* p_i
                                         ,int j
                                         ,const global Particle* p_j
                                         ,void* dataArray
                                         ,void* accum);

 void callback_PositionDelta_i(const global Parameters* parameters
                              ,int i
                              ,const global Particle* p_i
                              ,int j
                              ,const global Particle* p_j
                              ,void* dataArray
                              ,void* accum);

void callback_Curl_i(const global Parameters* parameters
                    ,int i
                    ,const global Particle* p_i
                    ,int j
                    ,const global Particle* p_j
                    ,void* dataArray
                    ,void* accum);

void callback_Vorticity_i(const global Parameters* parameters
                         ,int i
                         ,const global Particle* p_i
                         ,int j
                         ,const global Particle* p_j
                         ,void* dataArray
                         ,void* accum);

void callback_XPSHViscosity_i(const global Parameters* parameters
                             ,int i
                             ,const global Particle* p_i
                             ,int j
                             ,const global Particle* p_j
                             ,void* dataArray
                             ,void* accum);

global int getNeighboringCells(const global ParticlePosition* sortedParticleToCell
                              ,const global GridCellOffset* gridCellOffsets
                              ,int cellsX
                              ,int cellsY
                              ,int cellsZ
                              ,int3 cellSubscript
                              ,int* neighborCells);

global void forAllNeighbors(const global Parameters* parameters
                           ,const global Particle* particles
                           ,const global ParticlePosition* sortedParticleToCell
                           ,const global GridCellOffset* gridCellOffsets
                           ,int numParticles
                           ,int cellsX
                           ,int cellsY
                           ,int cellsZ
                           ,float3 minExtent
                           ,float3 maxExtent
                           ,int particleId
                           ,void* dataArray
                           ,void (*callback)(const global Parameters*
                                            ,int
                                            ,const global Particle*
                                            ,int
                                            ,const global Particle*
                                            ,void* dataArray
                                            ,void* currentAccum)
                           ,void* initialAccum);

/*******************************************************************************
 * Utility functions
 ******************************************************************************/

/**
 * A helper function that scales a value x in the range [a0,a1] to a new
 * range [b0,b1]
 */
global float rescale(float x, float a0, float a1, float b0, float b1)
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
global int sub2ind(int i, int j, int k, int w, int h)
{
    return i + (j * w) + k * (w * h);
}

/**
 * A function that converts a linear index x into a 3D subscript (i,j,k)
 *
 * @param [in] int x The linear index x
 * @param [in] int w grid width
 * @param [in] int h grid height
 */
global int3 ind2sub(int x, int w, int h)
{
    return (int3)(x % w, (x / w) % h, x / (w * h));
}

/**
 * Given a Particle, this function returns the 3D cell subscript of the cell the
 * particle is contained in
 *
 * @param Particle* p The particle to find the cell key of
 * @param [in] int cellsX The number of cells in the x axis of the spatial
 *             grid
 * @param [in] int cellsY The number of cells in the y axis of the spatial
 *             grid
 * @param [in] int cellsZ The number of cells in the z axis of the spatial
 *             grid
 * @param [in] float3 minExtent The minimum extent of the simulation's
 *             bounding box in world space
 * @param [in] float3 maxExtent The maximum extent of the simulation's
 *             bounding box in world space
 * @returns int3 The 3D subscript (i,j,k) of the cell the particle is
 *               contained in
 */
global int3 getSubscript(const global Particle* p
                       ,int cellsX
                       ,int cellsY
                       ,int cellsZ
                       ,float3 minExtent
                       ,float3 maxExtent)
{
    // Find the discretized cell the particle will be in according to its
    // predicted position:
    int i = (int)round((rescale(p->posStar.x, minExtent.x, maxExtent.x, 0.0f, (float)(cellsX - 1))));
    int j = (int)round((rescale(p->posStar.y, minExtent.y, maxExtent.y, 0.0f, (float)(cellsY - 1))));
    int k = (int)round((rescale(p->posStar.z, minExtent.z, maxExtent.z, 0.0f, (float)(cellsZ - 1))));
    
    return (int3)(i, j, k);
}

/**
 * Given a Particle, this function returns 1D cell index of the cell the 
 * particle is contained in
 *
 * @param Particle* p The particle to find the cell key of
 * @param [in] int cellsX The number of cells in the x axis of the spatial
 *             grid
 * @param [in] int cellsY The number of cells in the y axis of the spatial
 *             grid
 * @param [in] int cellsZ The number of cells in the z axis of the spatial
 *             grid
 * @param [in] float3 minExtent The minimum extent of the simulation's
 *             bounding box in world space
 * @param [in] float3 maxExtent The maximum extent of the simulation's
 *             bounding box in world space
 * @returns int The 1D key of the cell the particle is contained in
 */
global int getKey(const global Particle* p
                 ,int cellsX
                 ,int cellsY
                 ,int cellsZ
                 ,float3 minExtent
                 ,float3 maxExtent)
{
    int3 subscript = getSubscript(p, cellsX, cellsY, cellsZ, minExtent, maxExtent);

    // Compute the linear index as the key:
    return sub2ind(subscript.x, subscript.y, subscript.z, cellsX, cellsY);
}

/**
 * Given the subscript (i,j,k) as an int3 of a cell to search the vicinity of,
 * this function will return a count of valid neighboring cells (including
 * itself) in the range [1,27], e.g. between 1 and 27 neighboring cells are
 * valid and need to be searched for neighbors. The indices from 
 * [0 .. neighborCount-1] will be populated with the indices of neighboring 
 * cells in gridCellOffsets, such that for each neighboring grid cell
 * (i', j', k'), 0 <= i' < cellX, 0 <= j' < cellY, 0 <= k' < cellZ, and the
 * corresponding entry for cell (i',j',k') in gridCellOffsets has a cell 
 * start index != -1.
 *
 * @param [in]  ParticlePosition* sortedParticleToCell
 * @param [in]  GridCellOffset* gridCellOffsets
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  int3 cellSubscript
 * @param [out] int* neighborCells
 */
global int getNeighboringCells(const global ParticlePosition* sortedParticleToCell
                              ,const global GridCellOffset* gridCellOffsets
                              ,int cellsX
                              ,int cellsY
                              ,int cellsZ
                              ,int3 cellSubscript
                              ,int* neighborCells)
{
    int i = cellSubscript.x;
    int j = cellSubscript.y;
    int k = cellSubscript.z;
    
    // Count of valid neighbors:

    int neighborCellCount = 0;

    // We need to search the following potential 27 cells about (i,j,k):
    // (i + [-1,0,1], j + [-1,0,1], k + [-1,0,1]):

    int offsets[3] = { -1, 0, 1 };
    int I, J, K;
    
    // -1 indicates an invalid/non-existent neighbor:

    for (int i = 0; i < 27; i++) {
        neighborCells[i] = -1;
    }

    I = J = K = -1;
    
    for (int u = 0; u < 3; u++) {

        I = i + offsets[u]; // I = i-1, i, i+1

        for (int v = 0; v < 3; v++) {
        
            J = j + offsets[v]; // J = j-1, j, j+1

            for (int w = 0; w < 3; w++) {
            
                K = k + offsets[w]; // K = k-1, k, k+1
                
                if (   (I >= 0 && I < cellsX)
                    && (J >= 0 && J < cellsY)
                    && (K >= 0 && K < cellsZ))
                {
                    int key = sub2ind(I, J, K, cellsX, cellsY);

                    // The specified grid cell offset has a valid starting
                    // index, so we can return it as a valid neighbor:

                    if (gridCellOffsets[key].start != -1) {
                        neighborCells[neighborCellCount++] = key;
                    }
                }
            }
        }
    }
    
    return neighborCellCount;
}

/**
 * For all neighbors p_j of a particle p_i, this function will apply the given
 * function to all particle pairs (p_i, p_j), accumulating the result and
 * returning it
 *
 * @param [in]  Parameters* parameters
 * @param [in]  Particle* particles
 * @param [in]  ParticlePosition* sortedParticleToCell
 * @param [in]  GridCellOffset* gridCellOffsets
 * @param [in]  int numParticles The total number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [in]  int particleId The ID (index) of the particle to find the n
 *              neighbors of. This value corresponds to the position of the 
 *              particle in the array particles, and must be in the range
 *              [0 .. numParticles - 1]
 * @param [in]  (*callback)(int, const global Particle*, int, const global Particle*, void* accum)
 * @param [out] void* accum The accumulated result, passed to and update by apply
 *              for every neighbor pair of particles
 */
global void forAllNeighbors(const global Parameters* parameters
                           ,const global Particle* particles
                           ,const global ParticlePosition* sortedParticleToCell
                           ,const global GridCellOffset* gridCellOffsets
                           ,int numParticles
                           ,int cellsX
                           ,int cellsY
                           ,int cellsZ
                           ,float3 minExtent
                           ,float3 maxExtent
                           ,int particleId
                           ,void* dataArray
                           ,void (*callback)(const global Parameters*
                                            ,int
                                            ,const global Particle*
                                            ,int
                                            ,const global Particle*
                                            ,void* dataArray
                                            ,void* currentAccum)
                           ,void* initialAccum)
{
    // Sanity check:
    if (particleId < 0 || particleId >= numParticles) {
        return;
    }

    const global Particle *p_i = &particles[particleId];
    
    // Neighbor radius threshold:
    
    float R2 = 2.0f * parameters->particleRadius;
    
    // Keep track of the number of neighbors seen so far. If/when this
    // number exceeds CHECK_MAX_NEIGHBORS, we can bail out of the neighbor
    // search loop

    int neighborsSeen = 0;
    
#ifdef USE_BRUTE_FORCE_SEARCH

    // Exhaustively search the space for neighbors by computing the distance
    // between p_i and for all j in [0 .. numParticles - 1], p_j:

    for (int j = 0; j < numParticles; j++) {

        // Skip instances in which we'd be comparing a particle to itself:
        if (j == particleId) {
            continue;
        }

        const global Particle* p_j = &particles[j];
        float d = distance(p_i->posStar, p_j->posStar);

        if (R2 >= d) {
            
            if (neighborsSeen >= CHECK_MAX_NEIGHBORS) {
                return;
            }
            
            callback(parameters, particleId, p_i, j, p_j, dataArray, initialAccum);

            neighborsSeen++;
        }
    }

#else // Use fixed-radius neighbor search:

    // Given a particle, find the cell it's in based on its position:

    int3 cellSubscript = getSubscript(p_i, cellsX, cellsY, cellsZ, minExtent, maxExtent);

    // 27 (3x3x3) possible neighbors to search:

    int neighborCells[27];
    int neighborCellCount = getNeighboringCells(sortedParticleToCell
                                               ,gridCellOffsets
                                               ,cellsX
                                               ,cellsY
                                               ,cellsZ
                                               ,cellSubscript
                                               ,neighborCells);
    
    // For all neighbors found for the given cell at grid subscript (i,j, k):

    for (int j = 0; j < neighborCellCount; j++) {
        
        // We fetch the all indices returned in neighbors and check that
        // the corresponding entries in gridCellOffsets (if neighbors[j]
        // is valid):
        
        if (neighborCells[j] != -1) {
            
            const global GridCellOffset* g = &gridCellOffsets[neighborCells[j]];
            
            // If the start index of the grid-cell is valid, we iterate over
            // every particle we find in the cell:
            
            if (g->start != -1) {
                
                int start = g->start;
                int end   = start + g->length;
                
                for (int k = start; k < end; k++) {
                    
                    int J = sortedParticleToCell[k].particleIndex;
                    
                    // Skip instances in which we'd be comparing a particle to itself:
                    
                    if (particleId == J) {
                        continue;
                    }
                    
                    // The current potentially neighboring particle to check
                    // the distance of:
                    
                    const global Particle* p_j = &particles[J];
                    
                    // To determine if p_j is actually a neighbor of p_i, we
                    // test if the position delta is less then the sum of the
                    // radii of both particles. If p_j is a neighbor of p_i,
                    // we invoke the specified callback and accumulate the
                    // result:

                    float d = distance(p_i->posStar, p_j->posStar);

                    if (R2 >= d) {
                        
                        if (neighborsSeen >= CHECK_MAX_NEIGHBORS) {
                            return;
                        }

                        // Invoke the callback function to the particle pair
                        // (p_i, p_j), along with their respective indices,
                        // and accumulate the result into accum:

                        callback(parameters, particleId, p_i, J, p_j, dataArray, initialAccum);
                    
                        neighborsSeen++;
                    }
                }
            }
        }
    }
#endif
}

/*******************************************************************************
 * Density estimation functions
 ******************************************************************************/

/**
 * Computed the poly6 scalar smoothing kernel
 *
 * From the PBF slides SIGGRAPH 2013, pg. 13
 *
 * @param [in] float4 r distance
 * @param [in] float h Smoothing kernel radius
 * @returns float The computed scalar value
 */
float poly6(float4 r, float h)
{
    float rBar = length(r);

    if (rBar > h) {
        return 0.0f;
    }
    
    // (315 / (64 * PI * h^9)) * (h^2 - |r|^2)^3
    float h9 = (h * h * h * h * h * h * h * h * h);
    if (h9 == 0.0) {
        return 0.0f;
    }
    float A  = 1.566681471061f / h9;
    float B  = (h * h) - (rBar * rBar);

    return A * (B * B * B);
}

/**
 * Computes poly6 using a scalar value in place of r
 *
 * From the PBF slides SIGGRAPH 2013, pg. 13
 *
 * @param [in] float r distance
 * @param [in] float h Smoothing kernel radius
 * @returns float The computed scalar value
 */
float poly6_scalar(float q, float h)
{
    if (q > h) {
        return 0.0f;
    }
    
    // (315 / (64 * PI * h^9)) * (h^2 - |r|^2)^3
    float h9 = (h * h * h * h * h * h * h * h * h);
    if (h9 == 0.0) {
        return 0.0f;
    }
    float A  = 1.566681471061f / h9;
    float B  = (h * h) - (q * q);
    
    return A * (B * B * B);
}

/**
 * Computes the spiky smoothing kernel gradient
 *
 * From the PBF slides SIGGRAPH 2013, pg. 13
 *
 * @param [in] float4 r distance
 * @param [in] float h Smoothing kernel radius
 * @returns float4 The computed gradient in (x,y,z,w)
 */
float4 spiky(float4 r, float h)
{
    float rBar = length(r);
    
    if (rBar > h) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }

    // (45 / (PI * h^6)) * (h - |r|)^2 * (r / |r|)
    float h6   = (h * h * h * h * h * h);
    if (h6 == 0.0f) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }
    float A    = 14.323944878271f / h6;
    float B    = (h - rBar);
    float4 out = A * (B * B) * (r / (rBar + EPSILON));
    out[3] = 0.0f;
    return out;
}

/**
 * SPH density estimator for a pair of particles p_i and p_j for use as a 
 * callback function with forAllNeighbors()
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_SPHDensityEstimator_i(const global Parameters* parameters
                                   ,int i
                                   ,const global Particle* p_i
                                   ,int j
                                   ,const global Particle* p_j
                                   ,void* dataArray
                                   ,void* accum)
{
    // Cast the void pointer to the type we expect, so we can update the
    // variable accordingly:
    
    float* accumDensity = (float*)accum;

    (*accumDensity) += poly6(p_i->posStar - p_j->posStar, parameters->smoothingRadius);
}

/**
 * A callback function that computes the SPH gradient of a constraint 
 * function C_i, w.r.t a particle p_j for the case when i = j
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_SPHGradient_i(const global Parameters* parameters
                           ,int i
                           ,const global Particle* p_i
                           ,int j
                           ,const global Particle* p_j
                           ,void* dataArray
                           ,void* accum)
{
    // Cast the void pointer to the type we expect, so we can update the
    // variable accordingly:

    float4* gradVector = (float4*)accum;

    (*gradVector) += spiky(p_i->posStar - p_j->posStar, parameters->smoothingRadius);
}

/**
 * A callback function that computes the squared length of the SPH gradient 
 * of a constraint function C_i, w.r.t a particle p_j for the case when i != j
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_SquaredSPHGradientLength_j(const global Parameters* parameters
                                        ,int i
                                        ,const global Particle* p_i
                                        ,int j
                                        ,const global Particle* p_j
                                        ,void* dataArray
                                        ,void* accum)
{
    // Cast the void pointer to the type we expect, so we can update the
    // variable accordingly:
    
    float* totalGradLength = (float*)accum;

    float4 gradVector      = (INV_REST_DENSITY * -spiky(p_i->posStar - p_j->posStar, parameters->smoothingRadius));
    float gradVectorLength = length(gradVector);
    
    (*totalGradLength) += (gradVectorLength * gradVectorLength);
}

/**
 * A callback function that computes the position delta of a particle p_i 
 * given a neighbor particle p_j
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
 void callback_PositionDelta_i(const global Parameters* parameters
                              ,int i
                              ,const global Particle* p_i
                              ,int j
                              ,const global Particle* p_j
                              ,void* dataArray
                              ,void* accum)
{
    // Cast the void pointer to the type we expect, so we can update the
    // variable accordingly:
    
    float* lambda    = (float*)dataArray;
    float4* posDelta = (float4*)accum;

    float lambda_i = lambda[i];
    float lambda_j = lambda[j];

    // Introduce the artificial pressure corrector:
    
    float h         = parameters->smoothingRadius;
    float deltaQ    = 0.3f * h;
    float4 r        = p_i->posStar - p_j->posStar;
    float4 gradient = spiky(r, h);
    float n         = poly6(r, h);
    //float d         = poly6_scalar(deltaQ, h);
    float d         = poly6_scalar(0.0f, h);
    float nd        = d <= EPSILON ? 0.0f : n / d;
    float sCorr     = -parameters->artificialPressureK * pow(nd, parameters->artificialPressureN);

    (*posDelta) += ((lambda_i + lambda_j + sCorr) * gradient);
}

/**
 * A callback function that computes the curl force acting on a given
 * particle, p_i
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_Curl_i(const global Parameters* parameters
                    ,int i
                    ,const global Particle* p_i
                    ,int j
                    ,const global Particle* p_j
                    ,void* dataArray
                    ,void* accum)
{
    float4* omega_i = (float4*)accum;

    float4 v_ij        = p_i->vel - p_j->vel;
    float4 gradient_ij = spiky(p_i->posStar - p_j->posStar, parameters->smoothingRadius);

    (*omega_i) += cross(v_ij, gradient_ij);
}

/**
 * A callback function that computes the vorticity force acting on a given
 * particle, p_i
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_Vorticity_i(const global Parameters* parameters
                         ,int i
                         ,const global Particle* p_i
                         ,int j
                         ,const global Particle* p_j
                         ,void* dataArray
                         ,void* accum)
{
    float4* curl          = (float4*)dataArray;
    float4* omegaGradient = (float4*)accum;
    
    float4 r = p_i->posStar - p_j->posStar;
    float4 omegaBar = length(curl[i] - curl[j]);
    
    (*omegaGradient) += (omegaBar / r);
}

/**
 * A callback function that computes the XSPH viscosity acting on a given
 * particle, p_i
 *
 * @param [in] Parameters* parameters Simulation parameters
 * @param [in] int i The fixed index of particle i
 * @param [in] Particle* p_i The i-th (fixed) particle, particle p_i
 * @param [in] int j The varying index of particle j
 * @param [in] Particle* p_j The j-th (varying) particle, particle p_j
 * @param [in] void* dataArray An auxiliary readonly source data array to access
 * @param [in] void* accum An accumulator value to update
 */
void callback_XPSHViscosity_i(const global Parameters* parameters
                             ,int i
                             ,const global Particle* p_i
                             ,int j
                             ,const global Particle* p_j
                             ,void* dataArray
                             ,void* accum)
{
    float4* v_ij_sum = (float4*)accum;
    
    float4 v_ij = p_i->vel - p_j->vel;
    float W_ij  = poly6(p_i->posStar - p_j->posStar, parameters->smoothingRadius);
    
    (*v_ij_sum) += (W_ij * v_ij);
}

/*******************************************************************************
 * Kernels
 ******************************************************************************/

/**
 * For all particles p_i in particles, this kernel resets all associated
 * quantities, like density, etc.
 */
kernel void resetParticleQuantities(global Particle* particles
                                   ,global ParticlePosition* particleToCell
                                   ,global ParticlePosition* sortedParticleToCell
                                   ,global float* density
                                   ,global float* lambda
                                   ,global float4* posDelta)
{
    int id = get_global_id(0);
    global Particle *p           = &particles[id];
    global ParticlePosition *pp  = &particleToCell[id];
    global ParticlePosition *spp = &sortedParticleToCell[id];

    // Particle index; -1 indicates unset
    pp->particleIndex = pp->cellI = pp->cellJ = pp->cellK = pp->key = -1;
    spp->particleIndex = spp->cellI = spp->cellJ = spp->cellK = spp->key = -1;

    p->posStar   = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    p->velDelta  = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    density[id]  = 0.0f;
    lambda[id]   = 0.0f;
    posDelta[id] = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
}

/**
 * For all cells in the spatial grid, this kernel resets all associated 
 * quantities
 */
kernel void resetCellQuantities(global int* cellHistogram
                               ,global GridCellOffset* gridCellOffsets)
{
    int id = get_global_id(0);

    cellHistogram[id] = 0;

    gridCellOffsets[id].start  = -1;
    gridCellOffsets[id].length = -1;
}

/**
 * For all particles p_i in particles, this kernel applies external forces to the
 * velocity of p_i
 *
 * Currently, only applies gravity to the y component of the velocity.
 * Additional forces may be added later like wind and other forms of
 * turbulence, etc.
 *
 *   v_i = v_i + dt + f_external(x_i)
 *
 * @param [in/out] Particle* particles The particles to update
 * @param [in]     float dt The timestep
 */
kernel void applyExternalForces(global Particle* particles
                               ,float dt)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    // Apply the force of gravity along the y-axis:
    p->vel.y += (dt * -G);
}

/**
 * For all particles p_i in particles, this kernel updates the predicted 
 * position of p_i using an explicit Euler step like so:
 *
 * x_i = x_i + (dt * v_i), where x_i is the position of p_i and v_i is
 * the velocity of p_i
 *
 * @param [in/out] Particle* particles The particles to update
 * @param [in]     float dt The timestep
 */
kernel void predictPosition(global Particle* particles
                           ,float dt)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];

    // Explicit Euler step on the predicted particle position posStar (x*)
    p->posStar = p->pos + (dt * p->vel);
}

/**
 * For all particles p_i in particles, this kernel discretizes each p_i's
 * position into a grid of cells with dimensions specified by cellsPerAxis.
 *
 * @param [in] Particle* particles The particles to assign to cells
 * @param [out] int2* particleToCell Each entry contains a int2 pair
 * (i,j), where i is the particle in the i-th entry of particles, and j is
 * the linear index of the corresponding linear bin (j_x, j_y, j_z), where
 * 0 <= j_x < cellsPerAxis.x, 0 <= j_y < cellsPerAxis.y,
 * and 0 <= j_z < cellsPerAxis.z
 * @param [out] int* cellHistogram A histogram of counts of particles per cell
 * @param [in] int cellsX The number of cells in the x axis of the spatial
 *             grid
 * @param [in] int cellsY The number of cells in the y axis of the spatial
 *             grid
 * @param [in] int cellsZ The number of cells in the z axis of the spatial
 *             grid
 * @param [in] float3 minExtent The minimum extent of the simulation's
 *             bounding box in world space
 * @param [in] float3 maxExtent The maximum extent of the simulation's
 *             bounding box in world space
 */
kernel void discretizeParticlePositions(const global Particle* particles
                                       ,global ParticlePosition* particleToCell
                                       ,global int* cellHistogram
                                       ,int cellsX
                                       ,int cellsY
                                       ,int cellsZ
                                       ,float3 minExtent
                                       ,float3 maxExtent)
{
    int id                       = get_global_id(0);
    const global Particle *p     = &particles[id];
    global ParticlePosition *p2c = &particleToCell[id];
    
    // Convert the particle's position (x,y,z) to a grid cell subscript (i,j,k):
    int3 subscript = getSubscript(p, cellsX, cellsY, cellsZ, minExtent, maxExtent);
    
    // Convert the particle's position (x,y,z) to a linear index key:
    int key = getKey(p, cellsX, cellsY, cellsZ, minExtent, maxExtent);

    p2c->particleIndex = id;
    
    // Set the (i,j,k) index of the cell:

    p2c->cellI = subscript.x;
    p2c->cellJ = subscript.y;
    p2c->cellK = subscript.z;
    p2c->key   = key;

    // This is needed; "cellHistogram[z] += 1" won't work here as multiple
    // threads are modifying cellHistogram simultaneously:

    atomic_add(&cellHistogram[key], 1);
}

/**
 * NOTE: This kernel is meant to be run with 1 thread. This is necessary
 * since we have to perform a sort and perform some other actions which are
 * inherently sequential in nature
 *
 * This kernel basically performs a counting sort 
 * (http://en.wikipedia.org/wiki/Counting_sort) on the particles, sorting
 * them by the grid cell they were each assigned to. Rather than sorting by
 * a 3 dimensional subscript (i,j,k), we linearize the subscript, and sort by
 * that
 *
 * @see discretizeParticlePositions
 *
 * @param [in]     ParticlePosition* particleToCell
 * @param [in/out] int* cellHistogram
 * @param [out]    ParticlePosition* sortedParticleToCell
 * @param [out]    GridCellOffset* gridCellOffsets An array of size 
 *                 [0 .. numCells-1], where each index i contains the start and 
 *                 length of the i-th cell in the grid as it occurs in 
 *                 sortedParticleToCell
 * @param [in] int numParticles The total number of particles in the simulation
 * @param [in] int numCells The total number of cells in the spatial grid
 * @param [in] int cellsX The number of cells in the x axis of the spatial
 *             grid
 * @param [in] int cellsY The number of cells in the y axis of the spatial
 *             grid
 * @param [in] int cellsZ The number of cells in the z axis of the spatial
 *             grid
 */
kernel void countSortParticles(const global ParticlePosition* particleToCell
                              ,global int* cellHistogram
                              ,global ParticlePosition* sortedParticleToCell
                              ,global GridCellOffset* gridCellOffsets
                              ,int numParticles
                              ,int numCells
                              ,int cellsX
                              ,int cellsY
                              ,int cellsZ)
{
    // ==== Sorting code =======================================================
    
    // First step of counting sort is done already, since we calculated
    //the histogram (cellHistogram) in the discretizeParticlePositions kernel:
    
    int prefixSum = 0;
    int totalSum  = 0;

    // Second step of counting sort:
    for (int i = 0; i < numCells; i++) {
        prefixSum        = cellHistogram[i];
        cellHistogram[i] = totalSum;
        totalSum        += prefixSum;
    }

    // Final step of counting sort:
    for (int i = 0; i < numParticles; i++) {

        const global ParticlePosition* pp = &particleToCell[i];

        //int key = sub2ind(pp->cellI, pp->cellJ, pp->cellK, cellsX, cellsY);
        int key = pp->key;
        int j   = cellHistogram[key];

        sortedParticleToCell[j] = *pp;
        
        cellHistogram[key] += 1;
    }
    
    // ==== Binning code =======================================================
    
    // Now, the ParticlePosition entries of sortedParticleToCell are sorted in
    // ascending order by the value sub2ind(pp[i].cellI, pp[i].cellJ, pp[i].cellK, cellsX, cellsY),
    // where pp is an instance of ParticlePosition  at index i, such that
    // 0 <= i < numParticles.

    // Record the offsets per grid cell:
    // The i-th entry of the gridCellOffsets contains the start and length
    // of the i-th linearized grid cell in sortedParticleToCell

    int lengthCount = 1;
    int cellStart   = 0;
    int currentKey  = -1;
    int nextKey     = -1;
    
    // We traverse the list to find sequences of consecutive particles that
    // are assigned the same cell. We record the start and length to these
    // sequences and store the results in gridCellOffsets, so we can
    // quickly find all of the particles in a given cell quickly.

    for (int i = 0; i < (numParticles - 1); i++) {

        // Compare the particle position at index i and i+1:
        
        const global ParticlePosition* currentP = &sortedParticleToCell[i];
        const global ParticlePosition* nextP    = &sortedParticleToCell[i+1];

        // If two particles p and q have cell subscripts (p_x, p_y, p_z) and
        // (q_x, q_y, q_z), then the keys are the linearized indices for p and
        // q, p_key and q_key.

        //currentKey = sub2ind(currentP->cellI, currentP->cellJ, currentP->cellK, cellsX, cellsY);
        //nextKey    = sub2ind(nextP->cellI, nextP->cellJ, nextP->cellK, cellsX, cellsY);
        
        currentKey = currentP->key;
        nextKey    = nextP->key;
        
        // If p_key and q_key are equal, increase the length of the span:

        if (currentKey == nextKey) {

            lengthCount++;

        } else {
            
            // We hit a new key. Record this grid cell offset and continue;

            gridCellOffsets[currentKey].start  = cellStart;
            gridCellOffsets[currentKey].length = lengthCount;
            
            cellStart   = i + 1;
            lengthCount = 1;
        }
    }
    
    // For the last particle, since we iterate up to, but not including
    // the particle at index (numParticles - 1):

    if (nextKey != -1) {
        gridCellOffsets[nextKey].start  = cellStart;
        gridCellOffsets[nextKey].length = lengthCount;
    }
}

/**
 * From the Macklin & Muller paper: SPH density estimation
 * 
 * The SPH density estimator calculates \rho_i = \sum_j * m_j * W(p_i - p_j, h),
 * where \rho_i is the density of the i-th particle, m_j is the mass of the 
 * j-th particle, p_i - p_j is the position delta between the particles p_i and
 * p_j and h is the smoothing radius
 *
 * @param [in]  Parameters* parameters
 * @param [in]  Particle* particles
 * @param [in]  ParticlePosition* sortedParticleToCell
 * @param [in]  GridCellOffset* gridCellOffsets
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [out] float* density
 */
void kernel estimateDensity(const global Parameters* parameters
                           ,const global Particle* particles
                           ,const global ParticlePosition* sortedParticleToCell
                           ,const global GridCellOffset* gridCellOffsets
                           ,int numParticles
                           ,int cellsX
                           ,int cellsY
                           ,int cellsZ
                           ,float3 minExtent
                           ,float3 maxExtent
                           ,global float* density)
{
    int id = get_global_id(0);

    // For all neighboring particles p_j of the current particle (specified
    // by particles[id], aka p_i), apply the function estimateDensity for
    // all (p_i, p_j), accumulating the result into the density variable:

    float estDensity = 0.0f;
    
    forAllNeighbors(parameters
                   ,particles
                   ,sortedParticleToCell
                   ,gridCellOffsets
                   ,numParticles
                   ,cellsX
                   ,cellsY
                   ,cellsZ
                   ,minExtent
                   ,maxExtent
                   ,id
                   ,(void*)particles
                   ,callback_SPHDensityEstimator_i
                   ,(void*)&estDensity);

    density[id] = estDensity;
}

/**
 * For all particles p_i in particles, this kernel computes the density
 * constraint lambda value, defined as
 *
 *   \lambda_i = -C_i(p_1, ..., p_n) / \sum_k |\nabla(p_k) C_i|^2
 *
 * where,
 * 
 *   1) C_i(p_1, ..., p_n) = (\rho_i / \rho_0) - 1 = 0,
 *
 *   2) \rho_0 is the rest density, and
 *
 *   3) \rho_i is the density for particle p_i
 *
 * NOTE:
 * This corresponds to Figure (1) in the section "Enforcing Incompressibility"
 *
 * @param [in]  Parameters* parameters Simulation parameters
 * @param [in]  const Particle* particles The particles in the simulation
 * @param [in]  const ParticlePosition* sortedParticleToCell
 * @param [in]  const GridCellOffset* gridCellOffsets
 * @param [in]  const float* density The density per particle. The i-th entry
 *              contains the density for the i-th particle
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [out] float* lambda The constraint lambda value
 */
kernel void computeLambda(const global Parameters* parameters
                         ,const global Particle* particles
                         ,const global ParticlePosition* sortedParticleToCell
                         ,const global GridCellOffset* gridCellOffsets
                         ,const global float* density
                         ,int numParticles
                         ,int cellsX
                         ,int cellsY
                         ,int cellsZ
                         ,float3 minExtent
                         ,float3 maxExtent
                         ,global float* lambda)
{
    int id = get_global_id(0);

    // Compute the constraint value C_i(p_1, ... p_n) for all neighbors [1..n]
    // of particle i:

    float C_i = (density[id] * INV_REST_DENSITY) - 1.0f;
    
    float gradientSum = 0.0f;
    
    // ==== Case (2) k = i =====================================================

    float4 gv_i = (float4)(0.0f, 0.0f, 0.0f, 0.0f);

    forAllNeighbors(parameters
                   ,particles
                   ,sortedParticleToCell
                   ,gridCellOffsets
                   ,numParticles
                   ,cellsX
                   ,cellsY
                   ,cellsZ
                   ,minExtent
                   ,maxExtent
                   ,id
                   ,(void*)particles
                   ,callback_SPHGradient_i
                   ,(void*)&gv_i);
    
    float gv_iLength = length(INV_REST_DENSITY * gv_i);
    
    gradientSum += (gv_iLength * gv_iLength);
    
    // ==== Case (2) k = j =====================================================
    
    float gv_sLengths = 0.0f;
    
    forAllNeighbors(parameters
                   ,particles
                   ,sortedParticleToCell
                   ,gridCellOffsets
                   ,numParticles
                   ,cellsX
                   ,cellsY
                   ,cellsZ
                   ,minExtent
                   ,maxExtent
                   ,id
                   ,(void*)particles
                   ,callback_SquaredSPHGradientLength_j
                   ,(void*)&gv_sLengths);
    
    gradientSum += gv_sLengths;
    
    // ==== lambda_i ===========================================================

    if (gradientSum == 0.0f) {
        gradientSum = EPSILON;
    }
    
    lambda[id] = -(C_i / ((gradientSum + parameters->relaxation)));
}

/**
 * For all particles p_i in particles, this kernel computes the position
 * delta of p_i, p_i*
 *
 * @param [in]  Parameters* parameters Simulation parameters
 * @param [in]  const Particle* particles The particles in the simulation
 * @param [in]  const ParticlePosition* sortedParticleToCell
 * @param [in]  const GridCellOffset* gridCellOffsets
 * @param [in]  const float* density The density per particle. The i-th entry
 *              contains the density for the i-th particle
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [out] float4* posDelta position changes
 */
kernel void computePositionDelta(const global Parameters* parameters
                                ,const global Particle* particles
                                ,const global ParticlePosition* sortedParticleToCell
                                ,const global GridCellOffset* gridCellOffsets
                                ,int numParticles
                                ,const global float* lambda
                                ,int cellsX
                                ,int cellsY
                                ,int cellsZ
                                ,float3 minExtent
                                ,float3 maxExtent
                                ,global float4* posDelta)
{
    int id = get_global_id(0);

    float4 posDelta_i = (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    
    forAllNeighbors(parameters
                   ,particles
                   ,sortedParticleToCell
                   ,gridCellOffsets
                   ,numParticles
                   ,cellsX
                   ,cellsY
                   ,cellsZ
                   ,minExtent
                   ,maxExtent
                   ,id
                   ,(void*)lambda
                   ,callback_PositionDelta_i
                   ,(void*)&posDelta_i);
    
    posDelta[id] = INV_REST_DENSITY * posDelta_i;
}

/**
 * For all particles p_i in particles, this kernel applies the computed 
 * position delta to the predicted positions p_i.posStar.(x|y|z), e.g. "x_i*"
 * in the Position Based Fluids paper
 *
 * @param [in]  float4* posDelta The predicted delta position
 * @param [out] Particle* particles The particles in the simulation to be updated
 */
kernel void applyPositionDelta(const global float4* posDelta
                              ,global Particle* particles)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    p->posStar += posDelta[id];
}

/**
 * Tests for collisions between particles and objects/bounds and projects
 * the positions of the particles accordingly
 *
 * TODO: For now, this just clamps the particle to the world bounds
 *
 * @param [in]     Parameters* parameters Simulation parameters
 * @param [in/out] const Particle* particles The particles in the simulation
 * @param [in]     float3 minExtent The minimum extent of the simulation's
 *                 bounding box in world space
 * @param [in]     float3 maxExtent The maximum extent of the simulation's
 *                 bounding box in world space
 */
kernel void resolveCollisions(const global Parameters* parameters
                             ,global Particle* particles
                             ,float3 minExtent
                             ,float3 maxExtent)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    // Clamp predicted and actual positions:
    
    float R = parameters->particleRadius;
    
    // Actual positions:
    p->pos.x = clamp(p->pos.x, minExtent.x + R, maxExtent.x - R);
    p->pos.y = clamp(p->pos.y, minExtent.y + R, maxExtent.y - R);
    p->pos.z = clamp(p->pos.z, minExtent.z + R, maxExtent.z - R);
    
    // Predicted positions:
    p->posStar.x = clamp(p->posStar.x, minExtent.x + R, maxExtent.x - R);
    p->posStar.y = clamp(p->posStar.y, minExtent.y + R, maxExtent.y - R);
    p->posStar.z = clamp(p->posStar.z, minExtent.z + R, maxExtent.z - R);
}

/**
 * For all particles p_i in particles, this kernel computes the final, true
 * velocity of the particle in the current simulation step
 *
 * @param [in/out] Particle* particles The particles in the simulation to be updated
 * @param [in]     float dt The timestep
 */
kernel void updateVelocity(global Particle* particles
                          ,float dt)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    // Update the particle's final velocity based on the actual (x) and
    // predicted (x*) positions:
    
    p->vel = (1.0f / dt) * (p->posStar - p->pos);
}

/**
 * Computes the curl associated with each particle
 *
 * @param [in]  Parameters* parameters Simulation parameters
 * @param [in]  const Particle* particles The particles in the simulation
 * @param [in]  const ParticlePosition* sortedParticleToCell
 * @param [in]  const GridCellOffset* gridCellOffsets
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [out] float4 curl The curl associated with each particle
 */
kernel void computeCurl(const global Parameters* parameters
                       ,const global Particle* particles
                       ,const global ParticlePosition* sortedParticleToCell
                       ,const global GridCellOffset* gridCellOffsets
                       ,int numParticles
                       ,int cellsX
                       ,int cellsY
                       ,int cellsZ
                       ,float3 minExtent
                       ,float3 maxExtent
                       ,global float4* curl)
{
    int id = get_global_id(0);
    
    // Curl for particle i:
    float4 omega_i = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    
    forAllNeighbors(parameters
                    ,particles
                    ,sortedParticleToCell
                    ,gridCellOffsets
                    ,numParticles
                    ,cellsX
                    ,cellsY
                    ,cellsZ
                    ,minExtent
                    ,maxExtent
                    ,id
                    ,(void*)particles
                    ,callback_Curl_i
                    ,(void*)&omega_i);
    
    curl[id] = omega_i;
}

/**
 * Computes and applies the vorticity confinement force
 *
 * @param [in]  Parameters* parameters Simulation parameters
 * @param [in]  const Particle* particles The particles in the simulation
 * @param [in]  const ParticlePosition* sortedParticleToCell
 * @param [in]  const GridCellOffset* gridCellOffsets
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  float4* curl Computed curl for the particle p_i
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 * @param [out] float4* extForces Accumulated external forces acting on 
 *              particle p_i
 */
kernel void applyVorticity(const global Parameters* parameters
                          ,const global Particle* particles
                          ,const global ParticlePosition* sortedParticleToCell
                          ,const global GridCellOffset* gridCellOffsets
                          ,int numParticles
                          ,const global float4* curl
                          ,int cellsX
                          ,int cellsY
                          ,int cellsZ
                          ,float3 minExtent
                          ,float3 maxExtent
                          ,global float4* extForces)
{
    int id = get_global_id(0);
    
    // Curl force for particle i:
    float4 omegaGradient = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    
    forAllNeighbors(parameters
                    ,particles
                    ,sortedParticleToCell
                    ,gridCellOffsets
                    ,numParticles
                    ,cellsX
                    ,cellsY
                    ,cellsZ
                    ,minExtent
                    ,maxExtent
                    ,id
                    ,(void*)curl
                    ,callback_Vorticity_i
                    ,(void*)&omegaGradient);
    
    float n = length(omegaGradient);
    float4 N = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    
    if (n > EPSILON) {
        N = normalize(omegaGradient);
    }
    
    float4 f_curl = parameters->relaxation * cross(N, curl[id]);
}

/**
 * Computes and applies the XSPH viscosity to each particle as described
 * by Schechter and Bridson, 2012.
 *
 * @param [in]  Parameters* parameters Simulation parameters
 * @param [in]  const Particle* particles The particles in the simulation
 * @param [in]  const ParticlePosition* sortedParticleToCell
 * @param [in]  const GridCellOffset* gridCellOffsets
 * @param [in]  int numParticles The number of particles in the simulation
 * @param [in]  int cellsX The number of cells in the x axis of the spatial
 *              grid
 * @param [in]  int cellsY The number of cells in the y axis of the spatial
 *              grid
 * @param [in]  int cellsZ The number of cells in the z axis of the spatial
 *              grid
 * @param [in]  float3 minExtent The minimum extent of the simulation's
 *              bounding box in world space
 * @param [in]  float3 maxExtent The maximum extent of the simulation's
 *              bounding box in world space
 */
kernel void applyXSPHViscosity(const global Parameters* parameters
                              ,const global Particle* particles
                              ,const global ParticlePosition* sortedParticleToCell
                              ,const global GridCellOffset* gridCellOffsets
                              ,int numParticles
                              ,int cellsX
                              ,int cellsY
                              ,int cellsZ
                              ,float3 minExtent
                              ,float3 maxExtent)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    float4 v_i_sum = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    
    forAllNeighbors(parameters
                   ,particles
                   ,sortedParticleToCell
                   ,gridCellOffsets
                   ,numParticles
                   ,cellsX
                   ,cellsY
                   ,cellsZ
                   ,minExtent
                   ,maxExtent
                   ,id
                   ,(void*)particles
                   ,callback_XPSHViscosity_i
                   ,(void*)&v_i_sum);

    p->velDelta = (parameters->viscosityCoeff * v_i_sum);
}

/**
 * For all particles p_i in particles, this kernel computes the final, true
 * position of the particle in the current simulation step
 *
 * @param [in/out] Particle* particles The particles in the simulation to be updated
 * @param [out]    float4* The final render position of the particle to use in
 *                 OpenGL for rendering the instanced position of the particle
 * @param [in]     float dt The timestep
 */
kernel void updatePosition(global Particle* particles
                          ,global float4* renderPos
                          ,float dt)
{
    int id = get_global_id(0);
    global Particle *p = &particles[id];
    
    // Update the particle's final velocity based on the actual (x) and
    // predicted (x*) positions:

    p->vel += p->velDelta;
    
    // And finally the position:
    
    p->pos = p->posStar;
    
    renderPos[id] = p->pos;
    
    // We need this since we're using a 4 dimensional homogenous representation
    // of a 3D point, e.g. we need to show a point in (x,y,z), but our
    // representation is in (x,y,z,w), which OpenGL converts to
    // (x/w,y/z,z/w) so that it can be displayed in 3D. If w is zero, then
    // we'll never see any output, so if we set this explicitly to 1.0, then
    // everything will work correctly:
    
    renderPos[id].w = 1.0f;
}

