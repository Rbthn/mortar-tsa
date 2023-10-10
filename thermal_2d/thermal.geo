Include "../common/general_parameters.pro";
Include "geometry_parameters.pro";

SetFactory("OpenCASCADE");
cl_min = Min(cl1, cl2);

If (!USE_REFERENCE)
Point(1) = {0, 0, 0, cl1};
Point(2) = {L1, 0, 0, cl1};
Point(3) = {L2, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1 : 4};
Plane Surface(1) = {1};

Point(5) = {L1+w_Ins, 0, 0, cl2};
Point(6) = {1, 0, 0, cl2};
Point(7) = {1, 1, 0, cl2};
Point(8) = {L2+w_Ins, 1, 0, cl2};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {5 : 8};
Plane Surface(2) = {2};
Else // meshed insulation
Point(1) = {0, 0, 0, cl_min};
Point(2) = {L1, 0, 0, cl_min};
Point(3) = {L2, 1, 0, cl_min};
Point(4) = {0, 1, 0, cl_min};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1 : 4};
Plane Surface(1) = {1};

Point(5) = {L1+w_Ins, 0, 0, cl_min};
Point(6) = {1, 0, 0, cl_min};
Point(7) = {1, 1, 0, cl_min};
Point(8) = {L2+w_Ins, 1, 0, cl_min};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Curve Loop(2) = {5 : 8};
Plane Surface(2) = {2};

Line(9) = {2, 5};
Line(10) = {8, 3};
Curve Loop(3) = {9, -8, 10, -2};
Plane Surface(3) = {3};
EndIf

Physical Surface("LEFT", 1) = 1;
Physical Surface("RIGHT", 2) = 2;
Physical Line("A", 11) = 4;
Physical Line("B", 12) = 2;
Physical Line("C", 13) = 8;
Physical Line("D", 14) = 6;

// Aux. Surface for TSA
If (!USE_REFERENCE)
// Left aux
Point(21) = {L1+w_Ins/2, 0, 0, cl1};
Point(22) = {L2+w_Ins/2, 1, 0, cl1};
Line(21) = {21, 22};
Physical Line("Left aux. surface for mortaring", 1001) = {21};

// Right aux
Point(31) = {L1+w_Ins/2, 0, 0, cl2};
Point(32) = {L2+w_Ins/2, 1, 0, cl2};
Line(31) = {31, 32};
Physical Line("Right aux. surface for mortaring", 1002) = {31};

// Middle aux
Point(41) = {L1+w_Ins/2, 0, 0, cl_min};
Point(42) = {L2+w_Ins/2, 1, 0, cl_min};
Line(41) = {41, 42};
Physical Line("Middle aux. surface for mortaring", 1003) = {41};

// Middle aux
Point(101) = {L1+w_Ins/2, 0, 0, cl_min};
Point(102) = {L2+w_Ins/2, 1, 0, cl_min};
Line(101) = {101, 102};
Physical Line("Aux. surface for TSA", 10000) = {101};
Else
Physical Surface("Insulation", 3) = 3;
EndIf
