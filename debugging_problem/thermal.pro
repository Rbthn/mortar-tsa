Include "../common/general_parameters.pro";
Include "../common/material.pro";
Include "geometry_parameters.pro";
Include "pde_parameters.pro";

Group {
  coil_left = Region[phys_tag_coil_left];
  coil_right = Region[phys_tag_coil_right];
  yoke = Region[phys_tag_yoke];
  yoke_gap = Region[phys_tag_yoke_gap];

  coil_iso = Region[phys_tag_coil_iso];
  yoke_iso = Region[phys_tag_yoke_iso];

  // Aux. surfaces for TSA
  // definitions
  coil_aux_left = Region[phys_tag_coil_aux_left]; // Link constraint to left coil
  coil_aux_right = Region[phys_tag_coil_aux_right]; // Link constraint to right coil
  coil_aux = Region[phys_tag_coil_aux_center]; // TSA functions live here
  coil_left_interface_right = Region[phys_tag_coil_left_interface_right];
  coil_right_interface_left = Region[phys_tag_coil_right_interface_left];

  yoke_left_aux_top = Region[phys_tag_yoke_left_aux_top];
  yoke_left_aux_bottom = Region[phys_tag_yoke_left_aux_bottom];
  yoke_left_aux = Region[phys_tag_yoke_left_aux_center];
  coil_left_interface_bottom = Region[phys_tag_coil_left_interface_bottom];
  yoke_left_interface_top = Region[phys_tag_yoke_left_interface_top];

  yoke_right_aux_top = Region[phys_tag_yoke_right_aux_top];
  yoke_right_aux_bottom = Region[phys_tag_yoke_right_aux_bottom];
  yoke_right_aux = Region[phys_tag_yoke_right_aux_center];
  coil_right_interface_bottom = Region[phys_tag_coil_right_interface_bottom];
  yoke_right_interface_top = Region[phys_tag_yoke_right_interface_top];

  // assign indices for easier handling later
  shell_aux_0 = Region[{coil_aux}];
  shell_left_0 = Region[{coil_aux_left}];
  shell_right_0 = Region[{coil_aux_right}];
  shell_0 = Region[{shell_left_0, shell_aux_0, shell_right_0}];

  shell_aux_1 = Region[{yoke_left_aux}];
  shell_left_1 = Region[{yoke_left_aux_top}];
  shell_right_1 = Region[{yoke_left_aux_bottom}];
  shell_1 = Region[{shell_left_1, shell_aux_1, shell_right_1}];

  shell_aux_2 = Region[{yoke_right_aux}];
  shell_left_2 = Region[{yoke_right_aux_top}];
  shell_right_2 = Region[{yoke_right_aux_bottom}];
  shell_2 = Region[{shell_left_2, shell_aux_2, shell_right_2}];


  coils = Region[{ coil_left, coil_right }];
  vol_heat = Region[{coils, coil_iso, yoke, yoke_gap, yoke_iso}];

  Bnd_Neu_Left = Region[{}];
  Bnd_Neu_Right = Region[{}];
  Bnd_Neu_Yoke = Region[{}];
  Bnd_Neu = Region[{Bnd_Neu_Left, Bnd_Neu_Right, Bnd_Neu_Yoke}];

  Bnd_Rob_Left = Region[{}];
  Bnd_Rob_Right = Region[{102}];
  Bnd_Rob_Yoke = Region[{}];
  Bnd_Rob = Region[{Bnd_Rob_Left, Bnd_Rob_Right, Bnd_Rob_Yoke}];

  N_thinshells = 3;

}

Function {
  // Cooling condition
  alpha[] = CFUN_hHe[$1, $2];
  Tref = InitTemp;


  // --------------- Material --------------------------------------------------
  If (!USE_NONLIN_MATERIAL)
    kappa_B[] = 1.;
    kappa_I[] = 1.;
    heatCap_B[] = 1.;
    heatCap_I[] = 1.;

    kappa[coils] = kappa_B[];
    kappa[coil_iso] = kappa_I[];
    kappa[yoke_iso] = kappa_I[];
    kappa[yoke] = 1.;
    kappa[yoke_gap] = kappa_I[];
    heatCap[coils] = heatCap_B[];
    heatCap[coil_iso] = heatCap_I[];
    heatCap[yoke_iso] = heatCap_I[];
    heatCap[yoke] = 1.;
    heatCap[yoke_gap] = heatCap_I[];
  Else
    kappa[coils] = kappa_Bare[$1, $2];
    kappa[coil_iso] = kappa_Insulation[$1];
    kappa[yoke_iso] = kappa_Insulation[$1];
    kappa[yoke] = CFUN_kSteel[$1];
    kappa[yoke_gap] = kappa_Insulation[$1];
    heatCap[coils] = heatCap_Bare[$1, $2, $3];
    heatCap[coil_iso] = heatCap_Insulation[$1];
    heatCap[yoke_iso] = heatCap_Insulation[$1];
    heatCap[yoke] = CFUN_CvSteel[$1];
    heatCap[yoke_gap] = heatCap_Insulation[$1];
  EndIf

  // -------------- Excitation -------------------------------------------------
  // Quench heating
  vol_heat_src[All] = 0.;
  vol_heat_src[coil_left] = Q;

  shell_heat_src_0[] = 0.;
  shell_heat_src_1[] = 0.;
  shell_heat_src_2[] = 0.;

  // ----------- Boundary conditions -------------------------------------------
  f_neu[Bnd_Neu] = 1.;

  // ------------ 1D FE Matrices -----------------------------------------------

  /**
   * Indices of TSAs:
   * 0: Left coil - Right coil
   * 1: Left coil - yoke_gap
   * 2: Right coil - yoke_gap
   */
  w_Ins_0 = coil_wIns;
  w_Ins_1 = yoke_wIns;
  w_Ins_2 = yoke_wIns;
  N_ele_0 = 2;
  N_ele_1 = 2;
  N_ele_2 = 2;
  kappa_insulation_0[] = kappa_Insulation[$1, $2];
  cv_insulation_0[] = heatCap_Insulation[$1, $2, $3];
  kappa_insulation_1[] = kappa_Insulation[$1, $2];
  cv_insulation_1[] = heatCap_Insulation[$1, $2, $3];
  kappa_insulation_2[] = kappa_Insulation[$1, $2];
  cv_insulation_2[] = heatCap_Insulation[$1, $2, $3];

  For m In {0:N_thinshells-1}
  // element width of 1D FE element
    delta~{m} = w_Ins~{m} / (N_ele~{m} == 0 ? 1 : N_ele~{m});

    kappa_stiffness~{m}[shell_aux~{m}] = kKaptonStiffness[$1, $2, delta~{m}]{GaussOrder1D};

    kappa_mass~{m}[shell_aux~{m}] = kKaptonMass[$1, $2, $4, $5, delta~{m}]{GaussOrder1D};

    heatCap_mass~{m}[shell_aux~{m}] = CvKaptonMass[$1, $2, $4, $5, delta~{m}]{GaussOrder1D};

    heat_src~{m}[shell_aux~{m}] = 1/2 * delta~{m} * shell_heat_src~{m}[];
  EndFor
}

Constraint {
  { Name InitTemp; Type Init;
    Case {
      { Region shell_0; Value InitTemp; }
      { Region shell_1; Value InitTemp; }
      { Region shell_2; Value InitTemp; }
      { Region vol_heat; Value InitTemp; Type Init; }
    }
  }
  // coil TSA: Link constraints
  {
    Name Link_TSA; Type Link;
    Case{
      { Region coil_aux_left; RegionRef coil_left_interface_right;
        Function Vector[X[]-coil_wIns/2, Y[], Z[]]; Coefficient 1; }
      { Region coil_aux_right; RegionRef coil_right_interface_left;
        Function Vector[X[]+coil_wIns/2, Y[], Z[]]; Coefficient 1; }
      { Region yoke_left_aux_top; RegionRef coil_left_interface_bottom;
        Function Vector[X[], Y[]+yoke_wIns/2, Z[]]; Coefficient 1; }
      { Region yoke_left_aux_bottom; RegionRef yoke_left_interface_top;
        Function Vector[X[], Y[]-yoke_wIns/2, Z[]]; Coefficient 1; }
      { Region yoke_right_aux_top; RegionRef coil_right_interface_bottom;
        Function Vector[X[], Y[]+yoke_wIns/2, Z[]]; Coefficient 1; }
      { Region yoke_right_aux_bottom; RegionRef yoke_right_interface_top;
        Function Vector[X[], Y[]-yoke_wIns/2, Z[]]; Coefficient 1; }
    }
  }
}

FunctionSpace {
  // Reference
  {
    Name H1_reference; Type Form0;
    BasisFunction {
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Region[{vol_heat, Bnd_Neu, Bnd_Rob}]; Entity NodesOf[All]; }
    }
    Constraint {
      { NameOfCoef vn; EntityType NodesOf; NameOfConstraint InitTemp; }
    }
  }
  // TSA
  { Name H1; Type Form0;
    BasisFunction {
      { Name sn1; NameOfCoef vn1; Function BF_Node;
        Support Region[{coil_left, Bnd_Neu_Left, Bnd_Rob_Left}]; Entity NodesOf[ All ]; }
      { Name sn1_aux; NameOfCoef vn1_aux; Function BF_Node;
        Support Region[{coil_left_interface_right, coil_aux_left,
          coil_left_interface_bottom, yoke_left_aux_top}]; Entity NodesOf[ All ]; }
      { Name sn2; NameOfCoef vn2; Function BF_Node;
        Support Region[{coil_right, Bnd_Neu_Right, Bnd_Rob_Right}]; Entity NodesOf[ All ]; }
      { Name sn2_aux; NameOfCoef vn2_aux; Function BF_Node;
        Support Region[{coil_right_interface_left, coil_aux_right,
          coil_right_interface_bottom, yoke_right_aux_top}];
        Entity NodesOf[ All ]; }
      { Name sn3; NameOfCoef vn3; Function BF_Node;
        Support Region[{yoke, yoke_gap, Bnd_Neu_Yoke, Bnd_Rob_Yoke}]; Entity NodesOf[ All ]; }
      { Name sn3_aux; NameOfCoef vn3_aux; Function BF_Node;
        Support Region[{yoke_left_interface_top, yoke_left_aux_bottom,
          yoke_right_interface_top, yoke_right_aux_bottom}];
        Entity NodesOf[ All ]; }
    }
    Constraint {
      { NameOfCoef vn1; EntityType NodesOf; NameOfConstraint InitTemp; }
      { NameOfCoef vn2; EntityType NodesOf; NameOfConstraint InitTemp; }
      { NameOfCoef vn3; EntityType NodesOf; NameOfConstraint InitTemp; }
      { NameOfCoef vn1_aux; EntityType NodesOf; NameOfConstraint Link_TSA; }
      { NameOfCoef vn2_aux; EntityType NodesOf; NameOfConstraint Link_TSA; }
      { NameOfCoef vn3_aux; EntityType NodesOf; NameOfConstraint Link_TSA; }
    }
    SubSpace {
      { Name Left_Coil; NameOfBasisFunction {sn1, sn1_aux};}
      { Name Right_Coil; NameOfBasisFunction {sn2, sn2_aux};}
      { Name YokeAndGap; NameOfBasisFunction {sn3, sn3_aux}; }
    }
  }

  // TSA functions
  For m In {0:N_thinshells-1}
    For i In {0:N_ele~{m}}
      { Name H1_thinShell~{m}~{i} ; Type Form0 ;
        BasisFunction {
          { Name sn ; NameOfCoef vn~{i} ; Function BF_Node ;
          Support shell_aux~{m} ; Entity NodesOf[ All ] ; }
        }
        Constraint {
          { NameOfCoef vn~{i}; EntityType NodesOf; NameOfConstraint InitTemp; }
        }
      }
    EndFor
    // lagrange multiplier
    { Name H1_lambda~{m}~{1}; Type Form0;
      BasisFunction {
        { Name sn; NameOfCoef vn; Function BF_Node;
          Support shell_left~{m}; Entity NodesOf[All]; }
      }
    }
    { Name H1_lambda~{m}~{2}; Type Form0;
      BasisFunction {
        { Name sn; NameOfCoef vn; Function BF_Node;
          Support shell_right~{m}; Entity NodesOf[All]; }
      }
    }
  EndFor

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
      { Name v; Type Local; NameOfSpace H1_reference; }
    }
    Equation {
      Integral { [ kappa[{v}, magnField] * Dof{d v}, {d v} ];
      In vol_heat; Jacobian Vol; Integration Int; }
      Integral { DtDof [ heatCap[{v}, magnField, currDens/Jc * Ic] * Dof{v}, {v} ];
      In vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - vol_heat_src[], {v} ];
      In vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v} ];
      In Bnd_Neu; Jacobian Sur; Integration Int; }
      Integral { [ alpha[{v}, Tref] * Dof{v}, {v}];
      In Bnd_Rob; Jacobian Sur; Integration Int; }
      Integral { [ - alpha[{v}, Tref] * Tref, {v}];
      In Bnd_Rob; Jacobian Sur; Integration Int; }
    }
  }
  Else
  // ------------------------ Mortar & TSA -------------------------------------
  { Name heat_eqn; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace H1; }
      For m In {0:N_thinshells-1}
        For i In {0:N_ele_0}
          { Name vi~{m}~{i}; Type Local; NameOfSpace H1_thinShell~{m}~{i}; }
        EndFor
        { Name lambda~{m}~{1}; Type Local; NameOfSpace H1_lambda~{m}~{1}; }
        { Name lambda~{m}~{2}; Type Local; NameOfSpace H1_lambda~{m}~{2}; }
      EndFor
      { Name v1; Type Local; NameOfSpace H1[Left_Coil]; }
      { Name v2; Type Local; NameOfSpace H1[Right_Coil]; }
      { Name v3; Type Local; NameOfSpace H1[YokeAndGap]; }
    }
    Equation {
      /** PDE */
      Integral { [ kappa[{v}, magnField] * Dof{d v}, {d v} ];
        In vol_heat; Jacobian Vol; Integration Int; }
      Integral { DtDof[ heatCap[{v}, magnField, currDens/Jc * Ic] * Dof{v}, {v} ];
        In vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - vol_heat_src[], {v} ];
        In vol_heat; Jacobian Vol; Integration Int; }
      Integral { [ - f_neu[], {v} ];
        In Bnd_Neu; Jacobian Sur; Integration Int; }
      Integral { [ alpha[{v}, Tref] * Dof{v}, {v} ];
        In Bnd_Rob; Jacobian Sur; Integration Int; }
      Integral { [ - alpha[{v}, Tref] * Tref, {v} ];
        In Bnd_Rob; Jacobian Sur; Integration Int; }

      For m In {0:N_thinshells-1}
        /** 1D FE */
        For i In {0:N_ele~{m}-1} // loop over 1D FE elements
          For k In {1:2} // row of 1D FE matrix
            For l In {1:2} // column of 1D FE matrix
              Integral { [ kappa_mass~{m}[{vi~{m}~{i}}, {vi~{m}~{i+1}}, magnField, k, l] * Dof{d vi~{m}~{i + k - 1}} , {d vi~{m}~{i + l - 1}} ];
                  In shell_aux~{m}; Integration Int; Jacobian Sur; }
              Integral { DtDof[ heatCap_mass~{m}[{vi~{m}~{i}}, {vi~{m}~{i+1}}, magnField, k, l] * Dof{vi~{m}~{i + k - 1}} , {vi~{m}~{i + l - 1}} ];
                  In shell_aux~{m}; Integration Int; Jacobian Sur; }
              Integral { [ ((k==l) ? 1 : -1) * kappa_stiffness~{m}[{vi~{m}~{i}}, {vi~{m}~{i+1}}, magnField] * Dof{vi~{m}~{i + k - 1}} , {vi~{m}~{i + l - 1}} ];
                  In shell_aux~{m}; Integration Int; Jacobian Sur; }
            EndFor
            Integral { [ -1 * shell_heat_src~{m}[k], {vi~{m}~{i + k - 1}}];
              In shell_aux~{m}; Integration Int; Jacobian Sur; }
          EndFor
        EndFor
      EndFor

      /** Mortar for shell between coils */
      Integral { [ Dof{lambda~{0}~{1}}, {v1} ];
        In shell_left~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{lambda~{0}~{1}}, shell_left~{0}], {vi~{0}~{0}} ];
        In shell_aux~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v1}, {lambda~{0}~{1}} ];
        In shell_left~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{vi~{0}~{0}}, shell_aux~{0} ], {lambda~{0}~{1}} ];
        In shell_left~{0}; Jacobian Sur; Integration Int; }

      Integral { [ Dof{lambda~{0}~{2}}, {v2} ];
        In shell_right~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{lambda~{0}~{2}}, shell_right~{0} ], {vi~{0}~{N_ele~{0}}} ];
        In shell_aux~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v2}, {lambda~{0}~{2}} ];
        In shell_right~{0}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{vi~{0}~{N_ele~{0}}}, shell_aux~{0} ], {lambda~{0}~{2}}];
        In shell_right~{0}; Jacobian Sur; Integration Int; }

      /** Mortar for left yoke TSA */
      Integral { [ Dof{lambda~{1}~{1}}, {v1} ];
        In shell_left~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{lambda~{1}~{1}}, shell_left~{1}], {vi~{1}~{0}} ];
        In shell_aux~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v1}, {lambda~{1}~{1}} ];
        In shell_left~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{vi~{1}~{0}}, shell_aux~{1} ], {lambda~{1}~{1}} ];
        In shell_left~{1}; Jacobian Sur; Integration Int; }

      Integral { [ Dof{lambda~{1}~{2}}, {v3} ];
        In shell_right~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{lambda~{1}~{2}}, shell_right~{1} ], {vi~{1}~{N_ele~{1}}} ];
        In shell_aux~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Dof{v3}, {lambda~{1}~{2}} ];
        In shell_right~{1}; Jacobian Sur; Integration Int; }
      Integral { [ Trace [ - Dof{vi~{1}~{N_ele~{1}}}, shell_aux~{1} ], {lambda~{1}~{2}}];
        In shell_right~{1}; Jacobian Sur; Integration Int; }

      /** Mortar for right yoke TSA */
      Integral{ [ Dof{lambda~{2}~{1}}, {v2} ];
        In shell_left~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Trace [ - Dof{lambda~{2}~{1}}, shell_left~{2}], {vi~{2}~{0}} ];
        In shell_aux~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Dof{v2}, {lambda~{2}~{1}} ];
        In shell_left~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Trace [ - Dof{vi~{2}~{0}}, shell_aux~{2} ], {lambda~{2}~{1}} ];
        In shell_left~{2}; Jacobian Sur; Integration Int; }

      Integral{ [ Dof{lambda~{2}~{2}}, {v3} ];
        In shell_right~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Trace [ - Dof{lambda~{2}~{2}}, shell_right~{2}], {vi~{2}~{N_ele~{2}}} ];
        In shell_aux~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Dof{v3}, {lambda~{2}~{2}} ];
        In shell_right~{2}; Jacobian Sur; Integration Int; }
      Integral{ [ Trace [ - Dof{vi~{2}~{N_ele~{2}}}, shell_aux~{2} ], {lambda~{2}~{2}} ];
        In shell_right~{2}; Jacobian Sur; Integration Int; }
    }
  }
  EndIf
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
          Local { [ {v} ]; In vol_heat; Jacobian Vol; }
        }
      }
      // Maximum temperature in bare part from register 1 (saved in
      // post-operation by StoreMaxInRegister)
      { Name Tmax; Value{
        Term{
          Type Global; [#1]; In coils; Jacobian Vol;}
        }
      }
    }
  }
}

PostOperation {
  { Name Map; NameOfPostProcessing heat_v;
     Operation {
      If (USE_REFERENCE)
        Print [ v, OnElementsOf vol_heat, File "v_ref.pos" ];
      Else
        Print [ v, OnElementsOf vol_heat, File "v.pos" ];
      EndIf
     }
  }
  // print max. Temp. to file
  { Name PrintMaxTemp; NameOfPostProcessing heat_v;
    Operation {
      // get max., store in Reg. 1
      Print[ v, OnElementsOf coils, StoreMaxInRegister 1, Format Table,
        LastTimeStepOnly 1, File "dummy.txt"];
      //
      If (USE_REFERENCE)
        Print[ Tmax, OnRegion coils, Format TimeTable,
          File "Tmax_ref.txt", LastTimeStepOnly 1, AppendToExistingFile 1];
      Else
        Print[ Tmax, OnRegion coils, Format TimeTable,
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
