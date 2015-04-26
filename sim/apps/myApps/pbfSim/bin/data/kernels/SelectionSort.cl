
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

/**
 * Parallel selection sort
 * Eric Bainville, June 2011
 * http://www.bealto.com/gpu-sorting_parallel-selection.html
 */
kernel void selectionSort(global ParticlePosition* in, global ParticlePosition* out)
{
    int i = get_global_id(0); // current thread
    int n = get_global_size(0); // input size
    global ParticlePosition* iData = &in[i];
    
    int iKey = iData->key;
    
    // Compute position of in[i] in output
    
    int pos = 0;
    
    for (int j=0; j<n; j++) {
        int jKey     = in[j].key; // broadcasted
        bool smaller = (jKey < iKey) || (jKey == iKey && j < i);  // in[j] < in[i] ?
        pos += (smaller) ? 1 : 0;
    }
    
    out[pos] = *iData;
}
