#include<GradientEnhancedDamagedInternalForce.h>

//We define the valid parameters for this kernel and their default values
registerMooseObject("tardigradeApp", GradientEnhancedDamagedInternalForce);

template<>
InputParameters
validParams<GradientEnhancedDamagedInternalForce>(){
    InputParameters params = validParams<InternalForce>();
    params.addRequiredCoupledVar( "nonlocal_damage", "The nonlocal damage field" );
    return params;
}

GradientEnhancedDamagedInternalForce::GradientEnhancedDamagedInternalForce(const InputParameters & parameters)
    : InternalForce(parameters),
    _nonlocal_damage_var( coupled( "nonlocal_damage" ) ),
    _domega_dk( getMaterialPropertyByName< Real >( "domega_dk" ) )
{

}

Real GradientEnhancedDamagedInternalForce::computeQpOffDiagJacobian(unsigned int jvar){

    if ( jvar == _nonlocal_damage_var) 
        return - InternalForce::computeQpResidual() * _domega_dk[_qp] * _phi[_j][_qp];

    return InternalForce::computeQpOffDiagJacobian(jvar);
}
