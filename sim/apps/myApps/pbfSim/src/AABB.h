/*******************************************************************************
 * AABB.h
 * - A Simple axis-aligned bounding box implementation
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef AABB_H
#define AABB_H

#include <iostream>
#include "MSAOpenCL.h"
#include "Definitions.h"

/******************************************************************************/

class AABB
{
    protected:
        float3 minV, maxV;

    public:
        AABB();
        AABB(float3 p, float3 q);
        AABB(const AABB& other);
        AABB& operator=(const AABB& other);
    
        const float3& getMinExtent() const;
        const float3& getMaxExtent() const;
    
        friend std::ostream& operator<<(std::ostream& stream, const AABB& aabb);
};

/******************************************************************************/

#endif
