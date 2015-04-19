#define DAMP			0.95f
#define CENTER_FORCE	0.007f
#define MOUSE_FORCE		300.0f
#define MIN_SPEED		0.1f


typedef struct{
    float4 pos;
    float4 vel;
    float mass;
    float dummy[3];		// need this to make sure the float2 vel is aligned to a 16 byte boundary
} Particle;



__kernel void updateParticle(__global Particle* particles, const float4 mousePos, const float4 dimensions){
    int id = get_global_id(0);
    __global Particle *p = &particles[id];
    
    float4 diff = mousePos - p->pos;
    float invDistSQ = 1.0f / dot(diff, diff);
    diff *= MOUSE_FORCE * invDistSQ;
    
    p->vel += (dimensions*0.5f - p->pos) * CENTER_FORCE - diff* p->mass;
    
    float speed2 = dot(p->vel, p->vel);
    if(speed2<MIN_SPEED) p->pos = mousePos + diff * (1.0f + p->mass);
    
    p->pos += p->vel;
    p->vel *= DAMP;
}

