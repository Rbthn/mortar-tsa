Include "../common/general_parameters.pro";
Include "geometry_parameters.pro";

Group {
    If (USE_REFERENCE)
        Coil = Region[phys_tag_coil];
        dom_H_1 = Region[{}];
        dom_H_2 = Region[{}];
        dom_H = Region[{Coil}];
        Air = Region[phys_tag_air];
        dom_Phi = Region[{Air}];
    Else
        Coil_Top = Region[phys_tag_coil_top];
        Coil_Bottom = Region[phys_tag_coil_bottom];
        Coil = Region[{Coil_Top, Coil_Bottom}];
        Air_Top = Region[phys_tag_air_top];
        Air_Bottom = Region[phys_tag_air_bottom];
        Air = Region[{Air_Top, Air_Bottom}];

        dom_H_1 = Region[{Coil_Top}];
        dom_H_2 = Region[{Coil_Bottom}];
        dom_H = Region[{dom_H_1, dom_H_2}];
        dom_Phi = Region[{Air}];
    EndIf

    n_cuts = 1;
    Cut = Region[{}];
    For i In {1 : n_cuts}
        idx = phys_tag_cut_base+i;
        Cut~{i} = Region[idx];
        Cut = Region[{Cut, Cut~{i}}];
    EndFor

    Coil_Interface_Bottom = Region[phys_tag_coil_interface_bottom];
    Coil_Interface_Top = Region[phys_tag_coil_interface_top];
    Air_Interface_Bottom = Region[phys_tag_air_interface_bottom];
    Air_Interface_Top = Region[phys_tag_air_interface_top];

    Shell_Top = Region[phys_tag_coil_interface_top];
    Shell_Bottom = Region[phys_tag_coil_interface_bottom];
    Shell_Aux = Region[{Shell_Top}];
}

Function {
    mu_0 = 4*Pi*1e-7;
    mu[] = 1. * mu_0;
    I_max = 50.;
    t_rise = 1.;
    t_hold = 5.;
    t_fall = 0.5;
    I[] = InterpolationLinear[$1]{
        List[{
            0,0,
            t_rise,I_max,
            t_rise+t_hold,I_max,
            t_rise+t_hold+t_fall,0,
            t_rise+t_hold+t_fall+1,0}]
    };
    sigma[] = 1e5;

    // ------------ 1D FE Matrices -----------------------------------------------

    // element width of 1D FE element
    delta = w_Ins / (N_ele == 0 ? 1 : N_ele);

    sigma_stiffness[Shell_Aux] = (($1 == $2)? 1 : -1 ) * 1/sigma[]/delta;

    sigma_mass[Shell_Aux] = (($1 == $2)? 2: 1) * delta * 1/sigma[]/6;

    mu_mass[Shell_Aux] = (($1 == $2)? 2: 1) * delta * mu[]/6;
}


Constraint {
    { Name LinkAir; Type Link;
        Case {
            { Region Air_Interface_Top; RegionRef Air_Interface_Bottom;
                Function Vector[X[], Y[], Z[]]; Coefficient 1.; }
        }
    }
    { Name Current;
        Case {
            If (USE_REFERENCE)
                { Region Cut; Type Assign; Value 1; TimeFunction I[$Time]; }
            Else
                { Region Cut; Type Assign; Value 1; TimeFunction -I[$Time]; }
            EndIf
        }
    }
    { Name Voltage;
        Case {
        }
    }
}

FunctionSpace {
    { Name Hcurl_reference; Type Form1;
        BasisFunction {
            { Name sn; NameOfCoef hn; Function BF_Edge;
            Support Region[{dom_H}]; Entity EdgesOf[ All, Not {dom_Phi} ]; }
            { Name vn; NameOfCoef phi_n; Function BF_GradNode;
            Support Region[{dom_H, dom_Phi}]; Entity NodesOf[ dom_Phi ]; }
            { Name ci; NameOfCoef Ii; Function BF_GroupOfEdges;
            Support Region[{dom_H, dom_Phi}]; Entity GroupsOfEdgesOf[ Cut ]; }
        }
        GlobalQuantity {
            { Name I; Type AliasOf; NameOfCoef Ii; }
            { Name V; Type AssociatedWith; NameOfCoef Ii; }
        }
        Constraint {
            { NameOfCoef I; EntityType GroupsOfEdgesOf; NameOfConstraint Current; }
            { NameOfCoef V; EntityType GroupsOfEdgesOf; NameOfConstraint Voltage; }
        }
    }
    { Name Hcurl; Type Form1;
        BasisFunction {
            { Name vn; NameOfCoef phi_n; Function BF_GradNode;
            Support Region[{dom_Phi, dom_H, Shell_Top, Shell_Bottom}]; Entity NodesOf[ {dom_Phi} ]; }
            { Name sn2; NameOfCoef hn2; Function BF_Edge;
            Support Region[{dom_H~{2}, Shell_Bottom}]; Entity EdgesOf[ All, Not {dom_Phi} ]; }
            { Name sn1; NameOfCoef hn1; Function BF_Edge;
            Support Region[{dom_H~{1}, Shell_Top}]; Entity EdgesOf[ All, Not {dom_Phi} ]; }
            { Name ci; NameOfCoef Ii; Function BF_GroupOfEdges;
            Support Region[{dom_H, dom_Phi}]; Entity GroupsOfEdgesOf[Cut]; }
        }
        GlobalQuantity {
            { Name I; Type AliasOf; NameOfCoef Ii; }
            { Name V; Type AssociatedWith; NameOfCoef Ii; }
        }
        Constraint {
            { NameOfCoef phi_n; EntityType NodesOf; NameOfConstraint LinkAir; }
            { NameOfCoef I; EntityType GroupsOfEdgesOf; NameOfConstraint Current; }
            { NameOfCoef V; EntityType GroupsOfEdgesOf; NameOfConstraint Voltage; }
        }

        SubSpace {
            { Name Top_Coil; NameOfBasisFunction {sn1, vn};}
            { Name Bottom_Coil; NameOfBasisFunction {sn2, vn};}
        }
    }
    // TSA
    For i In {0:N_ele}
        { Name HCurl_shell~{i} ; Type Form0 ;
        BasisFunction {
            { Name sn ; NameOfCoef vn~{i} ; Function BF_Node ;
            Support Shell_Aux ; Entity EdgesOf[ All , Not dom_Phi ] ; }
        }
        }
    EndFor
    // Mortar
    { Name Hcurl_lambda_Top; Type Form1;
        BasisFunction {
            { Name sn; NameOfCoef vn; Function BF_Edge;
            Support Shell_Top; Entity EdgesOf [ All, Not dom_Phi ]; }
            { Name sn_grad; NameOfCoef vn_grad; Function BF_GradNode;
            Support Shell_Top; Entity NodesOf [ dom_Phi ]; }
        }
    }
    { Name Hcurl_lambda_Bottom; Type Form1;
        BasisFunction {
            { Name sn; NameOfCoef vn; Function BF_Edge;
            Support Shell_Bottom; Entity EdgesOf [ All, Not dom_Phi ]; }
            { Name sn_grad; NameOfCoef vn_grad; Function BF_GradNode;
            Support Shell_Bottom; Entity NodesOf [ dom_Phi ]; }
        }
    }
}

Jacobian {
    { Name Vol;
    Case {
        { Region All; Jacobian Vol;}
    }
    }
    { Name Sur;
    Case {
        { Region All; Jacobian Sur;}
    }
    }
}

Integration {
    { Name Int;
    Case {
        { Type Gauss;
            Case {
              { GeoElement Point;        NumberOfPoints  1; }
              { GeoElement Line;         NumberOfPoints  4; }
              { GeoElement Triangle;     NumberOfPoints  3; }
              { GeoElement Quadrangle;   NumberOfPoints  4; }
              { GeoElement Tetrahedron;  NumberOfPoints  5; }
              { GeoElement Hexahedron;   NumberOfPoints  6; }
              { GeoElement Prism;        NumberOfPoints  5; }
              { GeoElement Pyramid;      NumberOfPoints  8; }
            }
        }
    }
    }
}

Formulation {
    { Name mag; Type FemEquation;
      Quantity {
        If (USE_REFERENCE)
            { Name h; Type Local; NameOfSpace Hcurl_reference; }
            { Name I; Type Global; NameOfSpace Hcurl_reference[I]; }
            { Name V; Type Global; NameOfSpace Hcurl_reference[V]; }
        Else
            { Name h; Type Local; NameOfSpace Hcurl; }
            { Name h_1; Type Local; NameOfSpace Hcurl[Top_Coil]; }
            { Name h_2; Type Local; NameOfSpace Hcurl[Bottom_Coil]; }
            { Name I; Type Global; NameOfSpace Hcurl[I]; }
            { Name V; Type Global; NameOfSpace Hcurl[V]; }
            // TSA
            For i In {0:N_ele}
                { Name hi~{i}; Type Local; NameOfSpace HCurl_shell~{i}; }
            EndFor
            // Mortar
            { Name lambda_top; Type Local; NameOfSpace Hcurl_lambda_Top; }
            { Name lambda_bottom; Type Local; NameOfSpace Hcurl_lambda_Bottom; }
        EndIf
      }
      Equation {
          // REFERENCE
          If (USE_REFERENCE)
            // PDE
            Integral{[1/sigma[] * Dof{d h}, {d h}];
            In dom_H; Jacobian Vol; Integration Int; }
            Integral{DtDof[ mu[] * Dof{h}, {h}];
            In Region[{dom_H, dom_Phi}]; Jacobian Vol; Integration Int; }
            // excitation
            GlobalTerm{[ Dof{V}, {I}]; In Cut; }
          Else // Mortar
            // PDE
            Integral{[1/sigma[] * Dof{d h}, {d h}];
            In dom_H; Jacobian Vol; Integration Int; }
            Integral{DtDof[ mu[] * Dof{h}, {h}];
            In Region[{dom_H, dom_Phi}]; Jacobian Vol; Integration Int; }
            // excitation
            GlobalTerm{[Dof{V}, {I}]; In Cut; }
            // TSA
            For i In {0:N_ele-1} // loop over 1D FE elements
                For k In {1:2} // row of 1D FE matrix
                    For l In {1:2} // column of 1D FE matrix
                        Integral { [ sigma_mass[k, l] * Dof{d hi~{i + k - 1}} , {d hi~{i + l - 1}} ];
                            In Shell_Aux; Integration Int; Jacobian Sur; }
                        Integral { DtDof[ mu_mass[k, l] * Dof{hi~{i + k - 1}} , {hi~{i + l - 1}} ];
                            In Shell_Aux; Integration Int; Jacobian Sur; }
                        Integral { [ sigma_stiffness[k, l] * Dof{hi~{i + k - 1}} , {hi~{i + l - 1}} ];
                            In Shell_Aux; Integration Int; Jacobian Sur; }
                    EndFor
                EndFor
            EndFor
            // mortar top
            Integral{[Trace[ Dof{lambda_top}, Shell_Top], {hi~{0}}];
            In Shell_Aux; Jacobian Sur; Integration Int; }
            Integral{[- Dof{lambda_top}, {h_1}];
            In Shell_Top; Jacobian Sur; Integration Int; }
            Integral{[ Dof{h_1}, {lambda_top}];
            In Shell_Top; Jacobian Sur; Integration Int; }
            Integral{[- Trace[ Dof{hi~{0}}, Shell_Aux], {lambda_top}];
            In Shell_Top; Jacobian Sur; Integration Int; }

            // mortar bottom
            Integral{[Trace[ Dof{lambda_bottom}, Shell_Bottom], {hi~{N_ele}}];
            In Shell_Aux; Jacobian Sur; Integration Int; }
            Integral{[- Dof{lambda_bottom}, {h_2}];
            In Shell_Bottom; Jacobian Sur; Integration Int; }
            Integral{[ Dof{h_2}, {lambda_bottom}];
            In Shell_Bottom; Jacobian Sur; Integration Int; }
            Integral{[- Trace[ Dof{hi~{N_ele}}, Shell_Aux], {lambda_bottom}];
            In Shell_Bottom; Jacobian Sur; Integration Int; }
        EndIf
      }
    }
}

Resolution {
    { Name res;
        System {
            If (USE_REFERENCE)
                { Name sys_mag; NameOfFormulation mag; NameOfMesh "mag_reference.msh"; }
            Else
                { Name sys_mag; NameOfFormulation mag; NameOfMesh "mag.msh"; }
            EndIf
        }
        Operation {
            InitSolution[sys_mag];

            If (STEADY_STATE)
                // iterate to solve nonlinear problem, Picard iteration in this case
                IterativeLoopN[NMaxIt, relaxFactor,
                System { { sys_mag, nl_relTol, nl_absTol, Solution LinfNorm } } ]
                {
                // Generate and solve system
                Generate[sys_mag]; Solve[sys_mag];
                }
                SaveSolution[sys_mag];
            Else
                TimeLoopTheta [t_0, t_end, t_step, t_theta] {
                IterativeLoopN[NMaxIt, relaxFactor,
                    System { { sys_mag, nl_relTol, nl_absTol, Solution LinfNorm } } ]
                {
                    // Generate and solve system
                    Generate[sys_mag]; Solve[sys_mag];
                }
                SaveSolution[sys_mag];
                }
            EndIf
        }
    }
}

PostProcessing {
    { Name mag; NameOfFormulation mag;
    Quantity {
        { Name h; Value {
            Local{[ {h}]; In Region[{dom_H, dom_Phi}]; Jacobian Vol; }
        }}
        { Name b; Value {
            Local{[ mu[] * {h} ]; In Region[{dom_H, dom_Phi}]; Jacobian Vol; }
        }}
        { Name j; Value {
            Local{[ {d h} ]; In Region[{dom_H, dom_Phi}]; Jacobian Vol; }
        }}
        { Name I; Value {
            Local {[ {I} ]; In Cut; Jacobian Vol; }
        }}
        { Name V; Value {
            Local {[ {V} ]; In Cut; Jacobian Vol; }
        }}
    }}
}

PostOperation {
    { Name map; NameOfPostProcessing mag;
        If (USE_REFERENCE)
        Operation {
            Print[ h, OnElementsOf Region[{dom_Phi, dom_H}], File "h_ref.pos"];
            Print[ b, OnElementsOf Region[{dom_Phi, dom_H}], File "b_ref.pos"];
            Print[ j, OnElementsOf Region[{dom_H}], File "j_ref.pos"];
            Print[ I, OnRegion Cut, Format Table, File "I_ref.csv"];
            Print[ V, OnRegion Cut, Format Table, File "V_ref.csv"];
        }
        Else
        Operation {
            Print[ h, OnElementsOf Region[{dom_Phi, dom_H}], File "h.pos"];
            Print[ b, OnElementsOf Region[{dom_Phi, dom_H}], File "b.pos"];
            Print[ j, OnElementsOf Region[{dom_H}], File "j.pos"];
            Print[ I, OnRegion Cut, Format Table, File "I.csv"];
            Print[ V, OnRegion Cut, Format Table, File "V.csv"];
        }
        EndIf
    }
}