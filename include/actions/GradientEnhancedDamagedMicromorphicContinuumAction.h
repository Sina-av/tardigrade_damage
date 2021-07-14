#pragma once

#include "MicromorphicContinuumAction.h"

class GradientEnhancedDamagedMicromorphicContinuumAction : public MicromorphicContinuumAction
{
public:
  static InputParameters validParams();

  GradientEnhancedDamagedMicromorphicContinuumAction( const InputParameters & params );

  void act();

protected:

  void addKernels();
  void addMaterial();

};
