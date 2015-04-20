/*******************************************************************************
 * Parameters.h
 * CIS563: Physically Based Animation final project
 * Created by Michael Woods & Michael O'Meara
 ******************************************************************************/

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <iostream>

/*******************************************************************************
 * Tuneable parameters for the simulation
 ******************************************************************************/

typedef struct Parameters
{
    float particleRadius;      // Particle radius
    
    float particleMass;        // Particle mass
    
    float smoothingRadius;     // Kernel smoothing radius
    
    float relaxation;          // Pressure relaxation coefficient (epsilon), as
                               // described in the section 3 "Enforcing Incompressibility"
                               // of the Position Based Fluids paper
    float artificialPressureK; // Artificial pressure coefficient K
    
    float epsilonVorticity;    // Vorticity coefficient
    
    float viscosityCoeff;      // Viscosity coefficient
    
    float __padding1[3];
    
    int artificialPressureN;   // Artificial pressure coefficient N
    
    int __padding2[3];

    Parameters();
    Parameters(float _particleRadius
              ,float _particleMass
              ,float _smoothingRadius
              ,float _relaxation
              ,float _artificialPressureK
              ,int _artificialPressureN
              ,float _epsilonVorticity
              ,float _viscosityCoeff);
    Parameters(const Parameters& other);
    
    friend std::ostream& operator<<(std::ostream& os, Parameters p);

} Parameters;

#endif
