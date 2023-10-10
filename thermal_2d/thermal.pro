Include "../common/general_parameters.pro";
Include "geometry_parameters.pro";
Include "pde_parameters.pro";
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

Group {
  Left = Region[1];
  Right = Region[2];
  If (USE_REFERENCE)
    Insulation = Region[3];
  Else
    Insulation = Region[{}];
  EndIf

  A = Region[11]; // left
  B = Region[12]; // cut, left side
  C = Region[13]; // cut, right side. Ovelaps B if TSA is used
  D = Region[14]; // right

  Bare = Region[{ Left, Right }];
  Vol_heat = Region[ {Bare, Insulation} ];
  Bnd_Dir = Region[ {A, D} ];
  Bnd_Neu = Region[ {} ];

  // Mortar
  Left_aux = Region[1001];
  Right_aux = Region[1002];
  // used for TSA
  Shell_Aux = Region[10000];
}

Function {
  // --------------- Material --------------------------------------------------
  If (!USE_NONLIN_MATERIAL)
    kappa_Bare[] = 1.;
    kappa_Insulation[] = 1.;
    heatCap_Bare[] = 1.;
    heatCap_Insulation[] = 1.;
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

  kappa_stiffness[Shell_Aux] = (($1 == $2)? 1 : -1 ) * kappa_Insulation[$3, $4]/delta;

  kappa_mass[Shell_Aux] = (($1 == $2)? 2: 1) * delta * kappa_Insulation[$3, $4]/6;

  heatCap_mass[Shell_Aux] = (($1 == $2)? 2: 1) * delta * heatCap_Insulation[$3, $4, $5]/6;

  heat_src[Shell_Aux] = 1/2 * delta * shell_heat_src[];
}

Constraint {
  { Name Dirichlet_heat;
    Case {
      { Region A; Value 0.; }
      { Region D; Value InitTemp; }
      { Region Vol_heat; Value InitTemp; Type Init; }
    }
  }
  { Name InitTemp; Type Init;
    Case {
      { Region Shell_Aux; Value InitTemp; }
    }
  }
  { Name Link_Minus; Type Link;
    Case {
      { Region Left_aux; RegionRef B;
        Function Vector[X[]-w_Ins/2, Y[], Z[]]; Coefficient 1; }
    }
  }
  { Name Link_Plus; Type Link;
    Case {
      { Region Right_aux; RegionRef C;
        Function Vector[X[]+w_Ins/2, Y[], Z[]]; Coefficient 1; }
    }
  }
}

FunctionSpace {
  // TSA
  { Name H_v1; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Left, B, Bnd_Neu}]; Entity NodesOf[ All, Not Right ]; }
      { Name sn_aux; NameOfCoef vn_aux; Function BF_Node;
        Support Region[{B, Left_aux}]; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_heat; }
      { NameOfCoef vn_aux; EntityType NodesOf; NameOfConstraint Link_Minus; }
    }
  }
  { Name H_v2; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Right, C, Bnd_Neu}]; Entity NodesOf[ All, Not Left ]; }
      { Name sn_aux; NameOfCoef vn_aux; Function BF_Node;
        Support Region[{C, Right_aux}]; Entity NodesOf[ All ]; }
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
        Support Shell_Aux ; Entity NodesOf[ All ] ; }
      }
      Constraint {
        { NameOfCoef vn~{i}; EntityType NodesOf; NameOfConstraint InitTemp; }
      }
    }
  EndFor
  // lagrange multiplier
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
  // meshed insulation
  { Name H_v; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{Vol_heat, Bnd_Neu}]; Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint Dirichlet_heat; }
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
  { Name Int ;
    Case { { Type Gauss ;
             Case { { GeoElement Line        ; NumberOfPoints  4 ; }
                    { GeoElement Triangle    ; NumberOfPoints  4 ; }
                    { GeoElement Quadrangle  ; NumberOfPoints  4 ; } }
      }
    }
  }
}

Formulation {
  If (USE_REFERENCE)
  // ------------------------ Meshed Insulation --------------------------------
  { Name heat_eqn; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace H_v; }
    }
    Equation {
      Integral { [ kappa[{v}, magnField] * Dof{d v}, {d v} ];
      In Vol_heat; Jacobian Vol; Integration Int; }
      Integral { DtDof [ heatCap[{v}, magnField, currDens/Jc * Ic] * Dof{v}, {v} ];
      In Vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - vol_heat_src[], {v} ];
      In Vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v} ];
      In Bnd_Neu; Jacobian Sur; Integration Int; }
    }
  }
  Else
  // ------------------------ Mortar & TSA -------------------------------------
  { Name heat_eqn; Type FemEquation;
    Quantity {
      { Name v1; Type Local; NameOfSpace H_v1; }
      { Name v2; Type Local; NameOfSpace H_v2; }
      For i In {0:N_ele}
        { Name vi~{i}; Type Local; NameOfSpace H1_thinShell~{i}; }
      EndFor
      { Name lambda_1; Type Local; NameOfSpace H1_lambda_1; }
      { Name lambda_2; Type Local; NameOfSpace H1_lambda_2; }
    }
    Equation {
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
                In Shell_Aux; Integration Int; Jacobian Sur; }
            Integral { DtDof[ heatCap_mass[k, l, {vi~{i + k - 1}}, magnField, currDens/Jc * Ic] * Dof{vi~{i + k - 1}} , {vi~{i + l - 1}} ];
                In Shell_Aux; Integration Int; Jacobian Sur; }
            Integral { [ kappa_stiffness[k, l, {vi~{i + k - 1}}, magnField] * Dof{vi~{i + k - 1}} , {vi~{i + l - 1}} ];
                In Shell_Aux; Integration Int; Jacobian Sur; }
          EndFor
          Integral { [ -1 * shell_heat_src[k], {vi~{i + k - 1}}];
            In Shell_Aux; Integration Int; Jacobian Sur; }
        EndFor
      EndFor


      /** Mortar */
      Integral { [ Dof{lambda_1}, {v1} ];
        In Left_aux; Jacobian Sur; Integration Int; }
      Integral { [ Trace[ - Dof{lambda_1}, Left_aux], {vi~{0}} ];
        In Shell_Aux; Jacobian Sur; Integration Int; }
      Integral { [ - Dof{lambda_2}, {v2} ];
        In Right_aux; Jacobian Sur; Integration Int; }
      Integral { [ Trace[ Dof{lambda_2}, Right_aux], {vi~{N_ele}} ];
        In Shell_Aux; Jacobian Sur; Integration Int; }

      Integral { [ Dof{v1}, {lambda_1} ];
        In Left_aux; Jacobian Sur; Integration Int; }
      Integral { [ Trace[ - Dof{vi~{0}}, Shell_Aux], {lambda_1} ];
        In Left_aux; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v2}, {lambda_2} ];
        In Right_aux; Jacobian Sur; Integration Int; }
      Integral { [ Trace[ - Dof{vi~{N_ele}}, Shell_Aux], {lambda_2}];
        In Right_aux; Jacobian Sur; Integration Int; }
    }
  }
  EndIf
}

Resolution {
  { Name heat_v;
    System {
      { Name sys_heat; NameOfFormulation heat_eqn; }
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
            Local { [ {v} ]; In Vol_heat; Jacobian Vol; }
          EndIf
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
       Print [ v, OnElementsOf Vol_heat, File "v.pos" ];
     }
  }
  // print max. Temp. to file
  { Name PrintMaxTemp; NameOfPostProcessing heat_v;
    Operation {
      // get max., store in Reg. 1
      Print[ v, OnElementsOf Bare, StoreMaxInRegister 1, Format Table,
        LastTimeStepOnly 1, File "dummy.txt"];
      //
      Print[ Tmax, OnRegion Bare, Format TimeTable,
        File "Tmax.txt", LastTimeStepOnly 1, AppendToExistingFile 1];
    }
  }
  // clear max. Temp. output file
  { Name ClearMaxTemp; NameOfPostProcessing heat_v;
    Operation {
      Echo["", Format Table, File "Tmax.txt"];
    }
  }
}
