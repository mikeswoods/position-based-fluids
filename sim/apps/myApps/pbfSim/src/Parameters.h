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

    float smoothingRadius;     // Kernel smoothing radius
    
    float relaxation;          // Pressure relaxation coefficient (epsilon), as
                               // described in the section 3 "Enforcing Incompressibility"
                               // of the Position Based Fluids paper
    float artificialPressureK; // Artificial pressure coefficient K
    
    float artificialPressureN;   // Artificial pressure coefficient N
    
    float vorticityEpsilon;    // Vorticity coefficient

    float viscosityCoeff;      // Viscosity coefficient

    Parameters();
    Parameters(float _particleRadius
              ,float _smoothingRadius
              ,float _relaxation
              ,float _artificialPressureK
              ,float _artificialPressureN
              ,float _vorticityEpsilon
              ,float _viscosityCoeff);
    Parameters(const Parameters& other);
    
    friend std::ostream& operator<<(std::ostream& os, Parameters p);

} Parameters;

#endif
