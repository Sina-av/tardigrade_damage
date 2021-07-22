#pragma once

#include "Action.h"

class MicromorphicContinuumAction : public Action
{
public:
  static InputParameters validParams();

  MicromorphicContinuumAction( const InputParameters & params );

  void act();

protected:
  void addKernels();
  void addMaterial();

  /// these parameters are not passed to the invoked kernels
  const static std::vector< std::string > excludedParameters;

  void addInternalForceKernels(const std::string& internal_force_kernel_name);
  void addInternalCoupleKernels(const std::string& internal_couple_kernel_name);

  const unsigned int _ndisp;
  const unsigned int _nmicro_disp_grad;
};
