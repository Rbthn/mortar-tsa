DefineConstant[w_Ins = {0.001, Name "Dimensions/Width of insulator"} ];

// size of domain
size_x = 1.;
size_y = 1.;
size_z = 1.;

// size of conductor
w_x = 0.4;
w_y = 0.4;
w_z = size_z;

// physical tags to assign
phys_tag_coil_top = 10;
phys_tag_coil_bottom = 20;
phys_tag_air_top = 11;
phys_tag_air_bottom = 21;
phys_tag_coil_interface_top = 15;
phys_tag_coil_interface_bottom = 25;
phys_tag_air_interface_top = 16;
phys_tag_air_interface_bottom = 26;

// Cuts
phys_tag_cut_base = 100;


// Reference
phys_tag_coil = 10;
phys_tag_air = 11;
