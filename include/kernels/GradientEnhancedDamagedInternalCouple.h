#pragma once
#include "InternalCouple.h"

//Forward declarations
class GradientEnhancedDamagedInternalCouple;

template <>
InputParameters validParams<GradientEnhancedDamagedInternalCouple>();

class GradientEnhancedDamagedInternalCouple: public InternalCouple{
    public:
      GradientEnhancedDamagedInternalCouple(const InputParameters & parameters);

    protected:

      virtual Real computeQpOffDiagJacobian(unsigned int jvar) override; 
      /// The MOOSE variable number of the nonlocal damage variable
      unsigned int _nonlocal_damage_var;

      const MaterialProperty< Real > & _domega_dk;
};
