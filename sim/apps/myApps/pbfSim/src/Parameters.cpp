/*******************************************************************************
 * Parameters.cpp
 * - The heart of the position-based fluids simulator. This class encapsulates
 *   the current state of the simulation
 *
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#include "Parameters.h"

/******************************************************************************/

using namespace std;

/******************************************************************************/

//Parameters(0.5f, 1.0f, 1.1f, 0.005f, 0.001f, 4.0f, 0.1f, 0.01f);
Parameters::Parameters() :
    particleRadius(0.5f),
    smoothingRadius(1.1f),
    relaxation(0.005f),
    artificialPressureK(0.001f),
    artificialPressureN(4.0f),
    vorticityEpsilon(0.1f),
    viscosityCoeff(0.01f)
{

}

Parameters::Parameters(float _particleRadius
                      ,float _smoothingRadius
                      ,float _relaxation
                      ,float _artificialPressureK
                      ,float _artificialPressureN
                      ,float _vorticityEpsilon
                      ,float _viscosityCoeff) :
    particleRadius(_particleRadius),
    smoothingRadius(_smoothingRadius),
    relaxation(_relaxation),
    artificialPressureK(_artificialPressureK),
    artificialPressureN(_artificialPressureN),
    vorticityEpsilon(_vorticityEpsilon),
    viscosityCoeff(_viscosityCoeff)
{
    
}

Parameters::Parameters(const Parameters& other) :
    particleRadius(other.particleRadius),
    smoothingRadius(other.smoothingRadius),
    relaxation(other.relaxation),
    artificialPressureK(other.artificialPressureK),
    artificialPressureN(other.artificialPressureN),
    vorticityEpsilon(other.vorticityEpsilon),
    viscosityCoeff(other.viscosityCoeff)
{
    
}

ostream& operator<<(ostream& os, Parameters p)
{
    os << "Parameters {" << endl
    << "  particleRadius:\t"      << p.particleRadius      << endl
    << "  smoothingRadius:\t"     << p.smoothingRadius     << endl
    << "  relaxation:\t"          << p.relaxation          << endl
    << "  artificialPressureK:\t" << p.artificialPressureK << endl
    << "  artificialPressureN:\t" << p.artificialPressureN << endl
    << "  vorticityEpsilon:\t"    << p.vorticityEpsilon    << endl
    << "  viscosityCoeff:\t"      << p.viscosityCoeff    << endl
    << "}";
    return os;
}

/******************************************************************************/

