

kernel void increment(global float3* counts)
{
    int i = get_local_id(0);

    counts[i].x += 1.0;
    counts[i].y += 1.0;
    counts[i].z += 1.0;
}