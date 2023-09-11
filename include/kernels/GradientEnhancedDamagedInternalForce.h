#pragma once
#include "InternalForce.h"

class GradientEnhancedDamagedInternalForce : public InternalForce
{
public:
  GradientEnhancedDamagedInternalForce(const InputParameters & parameters);

  static InputParameters validParams();

protected:
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The MOOSE variable number of the nonlocal damage variable
  unsigned int _nonlocal_damage_var;

  const MaterialProperty<Real> & _mat_omega;
  const MaterialProperty<Real> & _domega_dk;
};
