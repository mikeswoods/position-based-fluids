
#define DAMP			0.95f
#define CENTER_FORCE	0.007f
#define MOUSE_FORCE		300.0f
#define MIN_SPEED		0.1f



typedef struct {
    float3 x;
    float3 v;
    float mass;
} Particle;


__kernel void updateParticle(__global Particle* particles)  //, __global float2* posBuffer, const float2 mousePos, const float2 dimensions)
{
    int id = get_global_id(0);
    __global Particle *p = &particles[id];
    
    p->v += float3(0,9.8f,0);
    
    /*
    float2 diff = mousePos - posBuffer[id];
    float invDistSQ = 1.0f / dot(diff, diff);
    diff *= MOUSE_FORCE * invDistSQ;
    
    p->v += (dimensions*0.5f - posBuffer[id]) * CENTER_FORCE - diff* p->mass;
    
    float speed2 = dot(p->v, p->v);
    if(speed2<MIN_SPEED) posBuffer[id] = mousePos + diff * (1.0f + p->mass);
    
    posBuffer[id] += p->vel;
    p->v *= DAMP;
     */
}

#define H               1.5f  // smoothing radius
#define H_9             (H*H*H*H*H*H*H*H*H) // h^9
#define H_6             (H*H*H*H*H*H) // h^6

__kernel float poly6Kernel(float3 p_i, float3 p_j){
    float3 diff = p_i - p_j;
    float r = sqrt(diff * diff);
    if (H > r && r > 0) {
        float h_minus_r = (H * H - r * r);
        float div = 64.0 * PI * H_9 * h_minus_r * h_minus_r * h_minus_r;
        return 315.0f / div;
    }
    return 0;
}

__kernel float spikyKernel(float3 p_i, float3 p_j){
    float3 diff = p_i - p_j;
    float r = sqrt(diff * diff);
    if (H > r && r > 0) {
        float h_minus_r = H - r;
        float div = PI * H_6 * h_minus_r * h_minus_r * h_minus_r;
        return 15.0f / div;
    }
    return 0;
}
