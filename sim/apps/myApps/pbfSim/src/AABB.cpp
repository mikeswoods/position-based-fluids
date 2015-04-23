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

ostream& operator<<(ostream& os, AABB& aabb)
{
    return os << "AABB { "
              << aabb.getMinExtent() << ", " << aabb.getMaxExtent()
              << " }";
}

/******************************************************************************/

AABB::AABB():
    minV(ofVec3f(0.0f, 0.0f, 0.0f)),
    maxV(ofVec3f(0.0f, 0.0f, 0.0f))
{
    
}

AABB::AABB(float3 _p, float3 _q) :
    minV(float3(std::min(_p.x, _q.x), std::min(_p.y, _q.y), std::min(_p.z, _q.z))),
    maxV(float3(std::max(_p.x, _q.x), std::max(_p.y, _q.y), std::max(_p.z, _q.z)))
{
    
}

AABB::AABB(const AABB& other) :
    minV(other.minV),
    maxV(other.maxV)
{
    
}

ofVec3f& AABB::getMinExtent()
{
    return this->minV;
}

ofVec3f& AABB::getMaxExtent()
{
    return this->maxV;
}

float AABB::width() const
{
    return this->maxV.x - this->minV.x;
}

float AABB::height() const
{
    return this->maxV.y - this->minV.y;
}

float AABB::depth() const
{
    return this->maxV.z - this->minV.z;
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
