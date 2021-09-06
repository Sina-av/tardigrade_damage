###############################################################################
#                                stretch_y.i                                  #
###############################################################################
# A ``sign of life'' test which makes sure that a simple stretching test      #
# problem runs to completion.                                                 #
###############################################################################
[Mesh]
  type = FileMesh
  file = mesh_coarse.inp
  construct_node_list_from_side_list = false # prevents from erronously adding side nodes to the alphabetically first nodeset
[]

[GlobalParams]
  volumetric_locking_correction = false
  displacements = 'disp_x disp_y disp_z'
  nonlocal_damage = nonlocal_damage
  order = SECOND
  family = LAGRANGE
[]

#[Mesh]
#  type = GeneratedMesh
#  displacements = 'disp_x disp_y disp_z'
#  dim = 3
#  nx = 1
#  ny = 1
#  nz = 1
##  file = unit_cube.e
#[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./phi_xx]
  [../]
  [./phi_yy]
  [../]
  [./phi_zz]
  [../]
  [./phi_yz]
  [../]
  [./phi_xz]
  [../]
  [./phi_xy]
  [../]
  [./phi_zy]
  [../]
  [./phi_zx]
  [../]
  [./phi_yx]
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable  = disp_x

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = disp_y

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = disp_z

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalCouple
    component_i = 0
    component_j = 0
    dof_num     = 3
    variable    = phi_xx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 4
    variable    = phi_xy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 5
    variable    = phi_xz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 6
    variable    = phi_yx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 7
    variable    = phi_yy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 8
    variable    = phi_yz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 9
    variable    = phi_zx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 10
    variable    = phi_zy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 11
    variable    = phi_zz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
[]

[AuxVariables]
  [./test]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pk2_22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pk2_33]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_33]
    order = FIRST
    family = MONOMIAL
  [../]
  [./macro_isv]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_isv]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_2]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_3]
    order = FIRST
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./test]
    type = MaterialStdVectorAux
    property = PK2
    index = 0
    variable = test
  [../]
[]

[AuxKernels]
  [./pk2_11]
    type = MaterialStdVectorAux
    property = PK2
    index = 0
    variable = pk2_11
  [../]
[]

[AuxKernels]
  [./pk2_22]
    type = MaterialStdVectorAux
    property = PK2
    index = 4
    variable = pk2_22
  [../]
[]

[AuxKernels]
  [./pk2_33]
    type = MaterialStdVectorAux
    property = PK2
    index = 8
    variable = pk2_33
  [../]
[]

[AuxKernels]
  [./sigma_11]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 0
    variable = sigma_11
  [../]
[]

[AuxKernels]
  [./sigma_22]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 4
    variable = sigma_22
  [../]
[]

[AuxKernels]
  [./sigma_33]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 8
    variable = sigma_33
  [../]
[]

[AuxKernels]
  [./macro_isv]
    type = MaterialStdVectorAux
    property = SDVS
    index = 0
    variable = macro_isv
  [../]
[]

[AuxKernels]
  [./micro_isv]
    type = MaterialStdVectorAux
    property = SDVS
    index = 1
    variable = micro_isv
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_1]
    type = MaterialStdVectorAux
    property = SDVS
    index = 2
    variable = micro_gradient_isv_1
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_2]
    type = MaterialStdVectorAux
    property = SDVS
    index = 3
    variable = micro_gradient_isv_2
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_3]
    type = MaterialStdVectorAux
    property = SDVS
    index = 4
    variable = micro_gradient_isv_3
  [../]
[]

# [BCs]
#   active = 'left_x back_z bottom_y top_y'
# #  active = 'left_x back_z bottom_y bottom_x top_y top_x'
#   [./left_x]
#     type = DirichletBC
#     #type = PresetBC
#     variable = disp_x
#     boundary = 'left'
#     #boundary = 'left right bottom top front back'
#     preset = true
#     value = 0
#   [../]
#   [./back_z]
#     type = DirichletBC
#     #type = PresetBC
#     variable = disp_z
#     boundary = 'back'
#     #boundary = 'left right bottom top front back'
#     preset = true
#     value = 0
#   [../]
#   [./bottom_x]
#     type = DirichletBC
#     #type = PresetBC
#     variable = disp_x
#     boundary = 'bottom'
#     #boundary = 'left right bottom top front back'
#     preset = true
#     value = 0
#   [../]
#   [./bottom_y]
#     type = DirichletBC
#     #type = PresetBC
#     variable = disp_y
#     boundary = 'bottom'
#     #boundary = 'left right bottom top front back'
#     preset = true
#     value = 0
#   [../]
#   [./top_x]
#     type     = DirichletBC
#     #type     = PresetBC
#     variable = disp_x
#     boundary = 'top'
#     preset = true
#     value    = 0
#   [../]
#   [./top_y]
#     #type = DirichletBC
#     #type = PresetBC
#     type = FunctionDirichletBC
#     variable = disp_y
#     boundary = 'top'
#     #boundary = 'left right bottom top front back'
#     preset = true
#     function = top_bc
#   [../]
# []

#[Functions]
#  [./top_bc]
#    type  = ParsedFunction
#    value = 0.1*t
#  [../]
#[]

[BCs]
  [bottom_x]
    type = DirichletBC
    variable = disp_x
    boundary = bottom
    value = 0
  []
  [bottom_y]
    type = DirichletBC
    variable = disp_y
    boundary = bottom
    value = 0
  []
  [bottom_z]
    type = DirichletBC
    variable = disp_z
    boundary = bottom 
    value = 0
  []
  [sym_z]
    type = DirichletBC
    variable = disp_z
    boundary = zsym
    value = 0 
  []
  [load]
    type = FunctionDirichletBC
    variable = disp_y
    boundary = load
    function = '1 * t'
  []
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '2 1e2        0 
                            2 2e2        0 
                            2 1e2        0 
                            
                            2 0.56 0.2
                            2 0. 0. 
                            2 0. 0. 

                            2 0.56 0.2
                            2 0. 0. 
                            2 0. 0. 

                            2 29480 25480 
                            5 1000 400 -1500 -1400 -3000 
                            11 0 0 0 0 0 0 1e+06 0 0 0 0 
                            2 400 -3000 
                            0.5 0.5 0.5 
                            1e-9 1e-9'

    number_SDVS = 55
    model_name = "LinearElasticityDruckerPragerPlasticity"

    #Coupled variables
    u1     = 'disp_x'
    u2     = 'disp_y'
    u3     = 'disp_z'
    phi_11 = 'phi_xx'
    phi_22 = 'phi_yy'
    phi_33 = 'phi_zz'
    phi_23 = 'phi_yz'
    phi_13 = 'phi_xz'
    phi_12 = 'phi_xy'
    phi_32 = 'phi_zy'
    phi_31 = 'phi_zx'
    phi_21 = 'phi_yx'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
    full = true
    petsc_options_iname = '     -pc_type
                                -pc_hypre_type
                                -ksp_type
                                -ksp_gmres_restart
                                -pc_hypre_boomeramg_relax_type_all
                                -pc_hypre_boomeramg_strong_threshold
                                -pc_hypre_boomeramg_agg_nl
                                -pc_hypre_boomeramg_agg_num_paths
                                -pc_hypre_boomeramg_max_levels
                                -pc_hypre_boomeramg_coarsen_type
                                -pc_hypre_boomeramg_interp_type
                                -pc_hypre_boomeramg_P_max
                                -pc_hypre_boomeramg_truncfactor' 

    petsc_options_value = '     hypre
                                boomeramg
                                gmres
                                301
                                symmetric-SOR/Jacobi
                                0.75
                                4 
                                2
                                25
                                Falgout
                                ext+i
                                0
                                0.1 '
  [../]
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  nl_rel_tol = 1e-12
  nl_abs_tol = 5e-8
  l_tol = 1e-5
  l_max_its = 150
  nl_max_its = 12
  nl_div_tol = 1e3

  automatic_scaling=true
  compute_scaling_once =true
  verbose=false

  line_search = none

  dtmin = 1e-7
  dtmax= 5e-3
  
  start_time = 0.0
  end_time = 1.0 

  num_steps = 2000

  [TimeStepper]
    type = IterationAdaptiveDT
    optimal_iterations = 8
    iteration_window = 3
    linear_iteration_ratio = 1000
    growth_factor=1.2
    cutback_factor=0.5
    dt = 2.5e-3
  []
  [Quadrature]
    order=SECOND
  []
  [Predictor]
    type = SimplePredictor
    scale = 1.0
    skip_after_failed_timestep = true
  []
[] 

[Outputs]
  interval = 1
  print_linear_residuals = false
  csv = true
  checkpoint = true
  [exodus]
   type=Exodus
    execute_on = 'timestep_end final'  
  []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'timestep_end final'  # Default is "final"
    level = 2             # Default is 1
  []
[]
