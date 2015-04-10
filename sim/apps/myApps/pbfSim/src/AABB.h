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
#include "Definitions.h"

/******************************************************************************/

class AABB
{
    protected:
        EigenVector3 minV, maxV;

    public:
        AABB();
        AABB(EigenVector3 p, EigenVector3 q);
        AABB(const AABB& other);
        AABB& operator=(const AABB& other);
    
        const EigenVector3& getMinExtent() const;
        const EigenVector3& getMaxExtent() const;
    
        friend std::ostream& operator<<(std::ostream& stream, const AABB& aabb);
};

/******************************************************************************/

#endif
