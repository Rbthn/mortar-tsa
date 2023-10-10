/** Geometry-specific parameters */
DefineConstant[L1 = {0.5, Min 0.01, Max 0.99, Step 0.01, Name "Dimensions/L1"}];
DefineConstant[L2 = {0.5, Min 0.01, Max 0.99, Step 0.01, Name "Dimensions/L2"}];
DefineConstant[cl1 = {0.15, Name "Dimensions/Size 1"}];
DefineConstant[cl2 = {0.25, Name "Dimensions/Size 2"}];
DefineConstant[w_Ins = {0.1, Name "Dimensions/Width of insulator"} ];

// size of domain
size_x = 1.;
size_y = 1.;
size_z = 1.;

w_y = size_y;
w_z = size_z;

// physical tags to assign
phys_tag_air = 3;
phys_tag_left = 1;
phys_tag_right = 2;
phys_tag_port_left = 10;
phys_tag_port_right = 11;
phys_tag_interface_right = 22;
phys_tag_interface_left = 21;
phys_tag_shell_aux = 31;
phys_tag_left_aux = 32;
phys_tag_right_aux = 33;