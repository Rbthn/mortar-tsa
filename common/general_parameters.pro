/** General parameters - independent of geometry or PDE */

// TSA
DefineConstant[ USE_REFERENCE = {0, Choices{0, 1}, Name "Input/99Custom/Use reference implementation"} ];
DefineConstant[ N_ele = {5, Min 0, Max 10, Name "Input/30TSA/No. of thin shell elements", Visible !USE_REFERENCE} ];
DefineConstant[ GaussOrder1D = {2, Min 1, Max 10, Name "Input/30TSA/Gauss interpolation order for thin shells", Visible !USE_REFERENCE} ];
DefineConstant[ USE_NONLIN_MATERIAL = {0, Choices{0, 1}, Name "Input/99Custom/Use nonlinear material model"}];

//Timestepping
DefineConstant[ STEADY_STATE = {0, Choices{0, 1}, Name "Input/20Timestepping/Compute steady-state solution"} ];
DefineConstant[ t_0 = {0., Name "Input/20Timestepping/Start time", Visible !STEADY_STATE }];
DefineConstant[ t_end = {1., Name "Input/20Timestepping/Stop time", Visible !STEADY_STATE } ];
DefineConstant[ t_step = {1., Name "Input/20Timestepping/Step size", Visible !STEADY_STATE } ];
DefineConstant[ t_theta = {1., Name "Input/20Timestepping/Theta", Visible !STEADY_STATE } ];

// NL iteration
NMaxIt      = 100; // max number of its
nl_relTol   = 1e-4; // rel stopping tolerance
nl_absTol   = 1e-2; // abs stopping tolerance
relaxFactor = 0.5; // relax factor