/** Geometry-specific parameters */
DefineConstant[ coil_xdim = {0.00134, Name "Geometry/Coil width"} ];
DefineConstant[ coil_ydim = {0.011, Name "Geometry/Coil height"} ];
DefineConstant[ coil_wIns = {2*8.4e-5, Name "Geometry/Coil insulation thickness"} ];
DefineConstant[ cl_coil = {0.001, Name "Geometry/Coil mesh size"} ];
DefineConstant[ yoke_xdim = {2.5*coil_xdim, Name "Geometry/Yoke width"} ];
DefineConstant[ yoke_ydim = {0.005, Name "Geometry/Yoke height"} ];
DefineConstant[ yoke_wIns = {8.4e-5, Name "Geometry/Yoke insulation thickness"} ];
DefineConstant[ yoke_dist = {0.001, Name "Geometry/Distance from yoke to insulation"} ];
DefineConstant[ cl_yoke = {0.001, Name "Geometry/Yoke mesh size"} ];
DefineConstant[ cl_air = {0.001, Name "Geometry/Air mesh size"} ];


/** Physical tags */

phys_tag_coil_left = 10;
phys_tag_coil_right = 11;
phys_tag_yoke = 12;
phys_tag_yoke_gap = 13;

// reference
phys_tag_coil_iso = 21;
phys_tag_yoke_iso = 22;

// mortar + TSA
phys_tag_coil_aux_left = 31;
phys_tag_coil_left_interface_right = 32;
phys_tag_coil_aux_right = 33;
phys_tag_coil_right_interface_left = 34;
phys_tag_coil_aux_center = 35;

phys_tag_yoke_left_aux_top = 41;
phys_tag_coil_left_interface_bottom = 42;
phys_tag_yoke_left_aux_bottom= 43;
phys_tag_yoke_left_interface_top = 44;
phys_tag_yoke_left_aux_center = 45;

phys_tag_yoke_right_aux_top = 51;
phys_tag_coil_right_interface_bottom = 52;
phys_tag_yoke_right_aux_bottom= 53;
phys_tag_yoke_right_interface_top = 54;
phys_tag_yoke_right_aux_center = 55;
