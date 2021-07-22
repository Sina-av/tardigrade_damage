#pragma once
#include "InternalForce.h"

//Forward declarations
class GradientEnhancedDamagedInternalForce;

template <>
InputParameters validParams<GradientEnhancedDamagedInternalForce>();

class GradientEnhancedDamagedInternalForce: public InternalForce{
    public:
      GradientEnhancedDamagedInternalForce(const InputParameters & parameters);

    protected:

      virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 

      /// The MOOSE variable number of the nonlocal damage variable
      unsigned int _nonlocal_damage_var;

      const MaterialProperty< Real > & _domega_dk;
};
