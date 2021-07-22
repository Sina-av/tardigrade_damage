#include "GradientEnhancedDamagedMicromorphicContinuumAction.h"
#include <string>
#include <vector>
#include "FEProblem.h"
#include "Factory.h"

registerMooseAction( "tardigradeApp", GradientEnhancedDamagedMicromorphicContinuumAction, "add_kernel" );

registerMooseAction( "tardigradeApp", GradientEnhancedDamagedMicromorphicContinuumAction, "add_material" );

InputParameters
GradientEnhancedDamagedMicromorphicContinuumAction::validParams()
{
  InputParameters params = MicromorphicContinuumAction::validParams();
  params.addRequiredCoupledVar( "nonlocal_damage", "The nonlocal damage field" );

  params.addRequiredParam<std::vector<Real>>( "gradient_enhanced_damage_fparameters", "The vector of floating point material parameters for the gradient-enhanced damage");
  params.addRequiredParam<Real>( "nonlocal_radius", "the nonlocal radius");

  return params;
}

GradientEnhancedDamagedMicromorphicContinuumAction::GradientEnhancedDamagedMicromorphicContinuumAction(
    const InputParameters & parameters )
  : MicromorphicContinuumAction( parameters )
{
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::act()
{
  if ( _current_task == "add_kernel" )
    addKernels();
  else if ( _current_task == "add_material" )
    addMaterial();
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::addGradientEnhancedMicromorphicDamageKernel(const std::string& gradient_enhanced_damage_kernel)
{
  InputParameters gradient_enhanced_damage_kernel_params = _factory.getValidParams( gradient_enhanced_damage_kernel );

  gradient_enhanced_damage_kernel_params.applyParameters( parameters(), excludedParameters );

  const std::string kernel_name = name() + "_gradient_enhanced_damage";

  gradient_enhanced_damage_kernel_params.set< NonlinearVariableName >( "variable" ) =
      getParam< std::vector< VariableName > >( "nonlocal_damage" )[0];

  _problem->addKernel( gradient_enhanced_damage_kernel, kernel_name, gradient_enhanced_damage_kernel_params );
}

void
GradientEnhancedDamagedMicromorphicContinuumAction::addKernels()
{

  MicromorphicContinuumAction::addInternalForceKernels("GradientEnhancedDamagedInternalForce");

  MicromorphicContinuumAction::addInternalCoupleKernels("GradientEnhancedDamagedInternalCouple");

  addGradientEnhancedMicromorphicDamageKernel( "GradientEnhancedMicromorphicDamage" );

}

void
GradientEnhancedDamagedMicromorphicContinuumAction::addMaterial()
{
  std::string materialType = "GradientEnhancedDamagedMicromorphicMaterial";

  auto materialParameters = _factory.getValidParams( materialType );
  materialParameters.applyParameters( parameters() , {"nonlocal_radius"} );

  materialParameters.set< std::string >( "model_name" ) =
      getParam< std::string >( "model_name" );

  materialParameters.set< std::vector< Real > >( "material_fparameters" ) =
      getParam< std::vector< Real > >( "material_fparameters" );

  materialParameters.set< std::vector< Real > >( "gradient_enhanced_damage_fparameters" ) = getParam< std::vector< Real > >( "gradient_enhanced_damage_fparameters" );

  _problem->addMaterial( materialType, name() + "_material", materialParameters );
}
