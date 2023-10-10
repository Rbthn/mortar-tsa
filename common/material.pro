// ------------ Material Functions -------------------------------------------
Function {
  //  material correcetions due to TSA. Not used.
  f_kappa[] = 1.;
  f_heatCap[] = 1.;

  // Parameters for the material fitting functions
  RRR_Cu     = 200; // Residual resistance ratio
  Tup_RRR_Cu = 295; // Residual resistance ratio

  C1 = 3336.363785; // for fitting IC of NbTi
  C2 = -248.917389; // for fitting IC of NbTi

  f_SC_strand_inGroup = 0.2; // percentage of SC mat inside strand

  partIc = 0.5; // normalized current (normalized by critical current)
  Jc     = 1; // critical current density (cancels out anyways)
  Ic     = 750; // critical current in A

  magnField = 5; // applied magnetic field in T
  currDens  = 1/2 * Jc; // current density in A/mm^2

  // thermal conductivity of the bare part w/o correction in W/mK
  kappa_SCSt[] = (1 - f_SC_strand_inGroup) *
    CFUN_kCu[$1, CFUN_rhoCu[$1, 0]{RRR_Cu, Tup_RRR_Cu},
    CFUN_rhoCu[$1, Norm[$2]]{RRR_Cu, Tup_RRR_Cu}]{RRR_Cu} +
    f_SC_strand_inGroup * 0;

  // thermal conductivity of the bare part after correction in W/mK
  kappa_Bare[] = f_kappa[] * kappa_SCSt[$1, $2];

  // thermal conductivity of the insulation in W/mK
  kappa_Insulation[] = CFUN_kKapton[$1];

  // heat capacity of bare part w/o correction in J/m^3 K
  heatCap_SCSt[] = (1 - f_SC_strand_inGroup) * CFUN_CvCu_T[$1] +
                      f_SC_strand_inGroup * CFUN_CvNbTi[$1, $2, $3, C1, C2];

  // heat capacity of the bare part after correction in J/m^3 K
  heatCap_Bare[] = f_heatCap[] * heatCap_SCSt[$1, $2, $3];

  // heat capacity of insulation in J/m^3 K
  heatCap_Insulation[] = CFUN_CvKapton[$1];
}