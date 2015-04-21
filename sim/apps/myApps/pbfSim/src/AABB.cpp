/*******************************************************************************
 * AABB.cpp
 * - A Simple axis-aligned bounding box implementation
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include <algorithm>
#include "AABB.h"

/******************************************************************************/

using namespace std;

/******************************************************************************/

ostream& operator<<(ostream& os, const AABB& aabb)
{
    return os << "AABB { "
              << aabb.getMinExtent() << ", " << aabb.getMaxExtent()
              << " }";
}

/******************************************************************************/

AABB::AABB():
    minV(float3(0.0f, 0.0f, 0.0f)),
    maxV(float3(0.0f, 0.0f, 0.0f))
{
    
}

AABB::AABB(float3 _p, float3 _q) :
    minV(float3(std::min(_p[0], _q[0]), std::min(_p[1], _q[1]), std::min(_p[2], _q[2]))),
    maxV(float3(std::max(_p[0], _q[0]), std::max(_p[1], _q[1]), std::max(_p[2], _q[2])))
{
    
}

AABB::AABB(const AABB& other) :
    minV(other.minV),
    maxV(other.maxV)
{
    
}

AABB& AABB::operator=(const AABB& other)
{
    if (this == &other) {
        return *this;
    }
    
    this->minV = other.minV;
    this->maxV = other.maxV;
    
    return *this;
}

const float3& AABB::getMinExtent() const
{
    return this->minV;
}

const float3& AABB::getMaxExtent() const
{
    return this->maxV;
}

void AABB::addMin(ofVec3f amount)
{
    this->minV += amount;
}

void AABB::addMax(ofVec3f amount)
{
    this->maxV += amount;
}

void AABB::addMinX(float amount)
{
    this->minV.x += amount;
}

void AABB::addMaxX(float amount)
{
    this->maxV.x += amount;
}

void AABB::addMinY(float amount)
{
    this->minV.y += amount;
}

void AABB::addMaxY(float amount)
{
    this->maxV.y += amount;
}

void AABB::addMinZ(float amount)
{
    this->minV.z += amount;
}

void AABB::addMaxZ(float amount)
{
    this->maxV.z += amount;
}


/******************************************************************************/
