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
  void addGradientEnhancedMicromorphicDamageKernel(const std::string& gradient_enhanced_damage_kernel);

};
