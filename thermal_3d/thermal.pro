Include "../common/general_parameters.pro";
Include "geometry_parameters.pro";
If (USE_NONLIN_MATERIAL)
  Include "../common/material.pro";
Else
  Function {
    magnField = 0.;
    currDens = 0.;
    Jc = 0.;
    Ic = 0.;
  }
EndIf

Group{
    Left = Region[phys_tag_left];
    Right = Region[phys_tag_right];
    Bare = Region[{Left, Right}];
    If (USE_REFERENCE)
      Insulation = Region[phys_tag_air];
    Else
      Insulation = Region[{}];
    EndIf
    Vol = Region[{Bare, Insulation}];

    Port_1 = Region[phys_tag_port_left];
    Port_2 = Region[phys_tag_port_right];
    Interface_left = Region[phys_tag_interface_left];
    Interface_right = Region[phys_tag_interface_right];

    Bnd_Neu = Region[{}];
    If (USE_REFERENCE)
      Left_aux = Region[{}];
      Right_aux = Region[{}];
      Shell_aux = Region[{}];
    Else
      Left_aux = Region[phys_tag_left_aux];
      Right_aux = Region[phys_tag_right_aux];
      Shell_aux = Region[phys_tag_shell_aux];
    EndIf
}

Function {
  // --------------- Material --------------------------------------------------
  If (!USE_NONLIN_MATERIAL)
    kappa_Bare[] = 400.;          // W/m/K
    kappa_Insulation[] = .25;     // Ws/K/m^3
    heatCap_Bare[] = 3.45e+06;    // W/m/K
    heatCap_Insulation[] = 3.e+06;// Ws/K/m^3
  EndIf

  kappa[Left] = kappa_Bare[$1, $2];
  kappa[Right] = kappa_Bare[$1, $2];
  kappa[Insulation] = kappa_Insulation[$1, $2];
  heatCap[Left] = heatCap_Bare[$1, $2, $3];
  heatCap[Right] = heatCap_Bare[$1, $2, $3];
  heatCap[Insulation] = heatCap_Insulation[$1, $2, $3];

  // -------------- Excitation -------------------------------------------------
  vol_heat_src[] = 0.;
  shell_heat_src[] = 0.;

  // ----------- Boundary conditions -------------------------------------------
  f_neu[Bnd_Neu] = 1.;


  // ------------ 1D FE Matrices -----------------------------------------------

  // element width of 1D FE element
  delta = w_Ins / (N_ele == 0 ? 1 : N_ele);

  kappa_stiffness[Shell_aux] = (($1 == $2)? 1 : -1 ) * kappa_Insulation[$3, $4]/delta;

  kappa_mass[Shell_aux] = (($1 == $2)? 2: 1) * delta * kappa_Insulation[$3, $4]/6;

  heatCap_mass[Shell_aux] = (($1 == $2)? 2: 1) * delta * heatCap_Insulation[$3, $4, $5]/6;

  heat_src[Shell_aux] = 1/2 * delta * shell_heat_src[];

  // ------------ Analytical solution  -----------------------------------------
  w_Bare = L1 - w_Ins/2;

  T_R = 4.2;
  T_L = 1.9;

  r_c[] = kappa_Bare[] / w_Bare;
  r_nc[] = kappa_Insulation[] / w_Ins;

  T_ana[Left] = T_L + X[] / w_Bare * r_nc[] / (2*r_nc[] + r_c[]) * (T_R - T_L);
  T_ana[Insulation] = 1 / (2*r_nc[]+r_c[]) * (r_nc[] * T_R + (r_c[] + r_nc[]) * T_L + (X[] - w_Bare)/w_Ins * r_c[] * (T_R-T_L));
  T_ana[Right] = 1 / (2*r_nc[]+r_c[]) * (r_nc[] * T_L + (r_c[] + r_nc[])*T_R + (X[] - w_Bare - w_Ins)/w_Bare * r_nc[] * (T_R-T_L));
}

Constraint {
  { Name Dirichlet_heat;
    Case {
      { Region Port_1; Value 1.9; }
      { Region Port_2; Value 4.2; }
      { Region Vol; Value 4.2; Type Init; }
    }
  }
  { Name InitTemp; Type Init;
    Case {
      { Region Shell_aux; Value 4.2; }
    }
  }
  { Name Link_Minus; Type Link;
    Case {
      { Region Left_aux; RegionRef Interface_left;
        Function Vector[X[]-w_Ins/2, Y[], Z[]]; Coefficient 1; }
    }
  }
  { Name Link_Plus; Type Link;
    Case {
      { Region Right_aux; RegionRef Interface_right;
        Function Vector[X[]+w_Ins/2, Y[], Z[]]; Coefficient 1; }
    }
  }
}

FunctionSpace{
  { Name H_reference; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Vol}]; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_heat; }
    }
  }
  { Name H_v1; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Left, Interface_left, Bnd_Neu}]; Entity NodesOf[ All, Not Right ]; }
      { Name sn_aux; NameOfCoef vn_aux; Function BF_Node;
        Support Region[{Interface_left, Left_aux}]; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_heat; }
      { NameOfCoef vn_aux; EntityType NodesOf; NameOfConstraint Link_Minus; }
    }
  }
  { Name H_v2; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Right, Interface_right, Bnd_Neu}]; Entity NodesOf[ All, Not Left ]; }
      { Name sn_aux; NameOfCoef vn_aux; Function BF_Node;
        Support Region[{Interface_right, Right_aux}]; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_heat; }
      { NameOfCoef vn_aux; EntityType NodesOf; NameOfConstraint Link_Plus; }
    }
  }
  For i In {0:N_ele}
    { Name H1_thinShell~{i} ; Type Form0 ;
      BasisFunction {
        { Name sn ; NameOfCoef vn~{i} ; Function BF_Node ;
        Support Shell_aux ; Entity NodesOf[ All ] ; }
      }
      Constraint {
        { NameOfCoef vn~{i}; EntityType NodesOf; NameOfConstraint InitTemp; }
      }
    }
  EndFor
  { Name H1_lambda_1; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Left_aux; Entity NodesOf[All]; }
    }
  }
  { Name H1_lambda_2; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Right_aux; Entity NodesOf[All]; }
    }
  }
}

Jacobian {
  { Name Vol ;
    Case {
      { Region All ; Jacobian Vol ; }
    }
  }
  { Name Sur ;
    Case {
      { Region All ; Jacobian Sur ; }
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
  { Name heat_eqn; Type FemEquation;
    Quantity {
      If (USE_REFERENCE)
      { Name v; Type Local; NameOfSpace H_reference; }
      Else
      { Name v1; Type Local; NameOfSpace H_v1; }
      { Name v2; Type Local; NameOfSpace H_v2; }
      For i In {0:N_ele}
        { Name vi~{i}; Type Local; NameOfSpace H1_thinShell~{i}; }
      EndFor
      { Name lambda_1; Type Local; NameOfSpace H1_lambda_1; }
      { Name lambda_2; Type Local; NameOfSpace H1_lambda_2; }
      EndIf
    }
    Equation {
      If (USE_REFERENCE)
      Integral { [ kappa[{v}, magnField] * Dof{d v}, {d v} ];
        In Vol; Jacobian Vol; Integration Int; }
      Integral { DtDof[ heatCap[{v}, magnField, currDens/Jc * Ic] * Dof{v}, {v} ];
        In Vol; Jacobian Vol; Integration Int; }
      Integral { [ - vol_heat_src[], {v} ];
        In Vol; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v} ];
        In Bnd_Neu; Jacobian Sur; Integration Int; }
      Else
      /** Left subdomain */
      Integral { [ kappa[{v1}, magnField] * Dof{d v1}, {d v1} ];
        In Left; Jacobian Vol; Integration Int; }
      Integral { DtDof[ heatCap[{v1}, magnField, currDens/Jc * Ic] * Dof{v1}, {v1} ];
        In Left; Jacobian Vol; Integration Int; }
      Integral { [ -1 * vol_heat_src[], {v1} ];
        In Left; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v1} ];
        In Bnd_Neu; Jacobian Sur; Integration Int; }


      /** Right subdomain */
      Integral { [ kappa[{v2}, magnField] * Dof{d v2}, {d v2} ];
        In Right; Jacobian Vol; Integration Int; }
      Integral { DtDof[ heatCap[{v2}, magnField, currDens/Jc * Ic] * Dof{v2}, {v2} ];
        In Right; Jacobian Vol; Integration Int; }
      Integral { [ -1 * vol_heat_src[], {v2} ];
        In Right; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v2} ];
        In Bnd_Neu; Jacobian Sur; Integration Int; }

      /** 1D FE */
      For i In {0:N_ele-1} // loop over 1D FE elements
        For k In {1:2} // row of 1D FE matrix
          For l In {1:2} // column of 1D FE matrix
            Integral { [ kappa_mass[k, l, {vi~{i + k - 1}}, magnField] * Dof{d vi~{i + k - 1}} , {d vi~{i + l - 1}} ];
                In Shell_aux; Integration Int; Jacobian Sur; }
            Integral { DtDof[ heatCap_mass[k, l, {vi~{i + k - 1}}, magnField, currDens/Jc * Ic] * Dof{vi~{i + k - 1}} , {vi~{i + l - 1}} ];
                In Shell_aux; Integration Int; Jacobian Sur; }
            Integral { [ kappa_stiffness[k, l, {vi~{i + k - 1}}, magnField] * Dof{vi~{i + k - 1}} , {vi~{i + l - 1}} ];
                In Shell_aux; Integration Int; Jacobian Sur; }
          EndFor
          Integral { [ -1 * shell_heat_src[k], {vi~{i + k - 1}}];
            In Shell_aux; Integration Int; Jacobian Sur; }
        EndFor
      EndFor

      // Mortar

      Integral { [ Dof{lambda_1}, {v1}];
        In Left_aux; Jacobian Sur; Integration Int; }
      Integral { [ - Trace[ Dof{lambda_1}, Left_aux], {vi~{0}}];
        In Shell_aux; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v1}, {lambda_1}];
        In Left_aux; Jacobian Sur; Integration Int; }
      Integral { [ - Trace [Dof{vi~{0}}, Shell_aux], {lambda_1}];
        In Left_aux; Jacobian Sur; Integration Int; }

      Integral { [ - Dof{lambda_2}, {v2}];
        In Right_aux; Jacobian Sur; Integration Int; }
      Integral { [ Trace[ Dof{lambda_2}, Right_aux], {vi~{N_ele}}];
        In Shell_aux; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v2}, {lambda_2}];
        In Right_aux; Jacobian Sur; Integration Int; }
      Integral { [ - Trace [Dof{vi~{N_ele}}, Shell_aux], {lambda_2}];
        In Right_aux; Jacobian Sur; Integration Int; }
    EndIf
    }
  }
}

Resolution {
  { Name heat_v;
    System {
      If (USE_REFERENCE)
        { Name sys_heat; NameOfFormulation heat_eqn; NameOfMesh "thermal_reference.msh"; }
      Else
        { Name sys_heat; NameOfFormulation heat_eqn; NameOfMesh "thermal.msh"; }
      EndIf
    }
    Operation {
      // reset max. temp.
      PostOperation[ClearMaxTemp];
      // initial state
      InitSolution[sys_heat];
      // print init. temp.
      PostOperation[PrintMaxTemp];
      // timestepping
      If (STEADY_STATE)
        // iterate to solve nonlinear problem, Picard iteration in this case
        IterativeLoopN[NMaxIt, relaxFactor,
          System { { sys_heat, nl_relTol, nl_absTol, Solution LinfNorm } } ]
        {
          // Generate and solve system
          Generate[sys_heat]; Solve[sys_heat];
        }
        SaveSolution[sys_heat];
        PostOperation[PrintMaxTemp];
      Else
        TimeLoopTheta [t_0, t_end, t_step, t_theta] {
          IterativeLoopN[NMaxIt, relaxFactor,
            System { { sys_heat, nl_relTol, nl_absTol, Solution LinfNorm } } ]
          {
            // Generate and solve system
            Generate[sys_heat]; Solve[sys_heat];
          }
          SaveSolution[sys_heat];
          PostOperation[PrintMaxTemp];
        }
      EndIf
    }
  }
}

PostProcessing {
  { Name heat_v; NameOfFormulation heat_eqn;
    Quantity {
      { Name v; Value {
          If (!USE_REFERENCE)
            Local { [ {v1} ]; In Left; Jacobian Vol; }
            Local { [ {v2} ]; In Right; Jacobian Vol; }
          Else
            Local { [ {v} ]; In Vol; Jacobian Vol; }
          EndIf
        }
      }
      { Name v_ana; Value {
          Local { [ T_ana[] ]; In Bare; Jacobian Vol; }
        }
      }
      // Maximum temperature in bare part from register 1 (saved in
      // post-operation by StoreMaxInRegister)
      { Name Tmax; Value{
        Term{
          Type Global; [#1]; In Bare; Jacobian Vol;}
        }
      }
    }
  }
}

PostOperation {
  { Name Map; NameOfPostProcessing heat_v;
     Operation {
      If (USE_REFERENCE)
        Print [ v, OnElementsOf Vol, File "v_ref.pos" ];
      Else
        Print [ v, OnElementsOf Vol, File "v.pos" ];
        Print [v_ana, OnElementsOf Vol, File "v_ana.pos" ];
      EndIf
     }
  }
  // print max. Temp. to file
  { Name PrintMaxTemp; NameOfPostProcessing heat_v;
    Operation {
      // get max., store in Reg. 1
      Print[ v, OnElementsOf Bare, StoreMaxInRegister 1, Format Table,
        LastTimeStepOnly 1, File "dummy.txt"];
      //
      If (USE_REFERENCE)
        Print[ Tmax, OnRegion Bare, Format TimeTable,
          File "Tmax_ref.txt", LastTimeStepOnly 1, AppendToExistingFile 1];
      Else
        Print[ Tmax, OnRegion Bare, Format TimeTable,
          File "Tmax.txt", LastTimeStepOnly 1, AppendToExistingFile 1];
      EndIf
    }
  }
  // clear max. Temp. output file
  { Name ClearMaxTemp; NameOfPostProcessing heat_v;
    Operation {
      If (USE_REFERENCE)
        Echo["", Format Table, File "Tmax_ref.txt"];
      Else
        Echo["", Format Table, File "Tmax.txt"];
      EndIf
    }
  }
}