// Define the addresses of the variables which will be added to the TTree
float eventWeight_S, eventWeight_B;
float trijetPt_S, trijetPt_B;
float trijetEta_S, trijetEta_B;
float trijetPhi_S, trijetPhi_B;
float trijetMass_S, trijetMass_B;
float trijetPFCharge_S, trijetPFCharge_B;
float trijetCombinedCvsL_S, trijetCombinedCvsL_B;
float trijetDeepCvsL_S, trijetDeepCvsL_B;
float trijetPtD_S, trijetPtD_B;
float trijetPtDR_S, trijetPtDR_B;
float trijetAxis1_S, trijetAxis1_B;
float trijetAxis2_S, trijetAxis2_B;
float trijetMult_S, trijetMult_B;
float trijetQGLikelihood_S, trijetQGLikelihood_B;
float trijetQGLikelihood_avg_S, trijetQGLikelihood_avg_B;
float trijetChiSquared_S, trijetChiSquared_B;
  
// Dijet Variables
float dijetPt_S, dijetPt_B;
float dijetEta_S, dijetEta_B;
float dijetPhi_S, dijetPhi_B;
float dijetMass_S, dijetMass_B;
float dijetPtDR_S, dijetPtDR_B;
float dijetPFCharge_S, dijetPFCharge_B;
float dijetCombinedCvsL_S, dijetCombinedCvsL_B;
float dijetDeepCvsL_S, dijetDeepCvsL_B;
float dijetPtD_S, dijetPtD_B;
float dijetAxis1_S, dijetAxis1_B;
float dijetAxis2_S, dijetAxis2_B;
int dijetMult_S, dijetMult_B;
float dijetQGLikelihood_S, dijetQGLikelihood_B;
float dijetQGLikelihood_avg_S, dijetQGLikelihood_avg_B;
float dijetChiSquared_S, dijetChiSquared_B;
  
// Leading jet from dijet system
float LdgJetPt_S, LdgJetPt_B;
float LdgJetEta_S, LdgJetEta_B;
float LdgJetPhi_S, LdgJetPhi_B;
float LdgJetMass_S, LdgJetMass_B;
float LdgJetPFCharge_S, LdgJetPFCharge_B;
float LdgJetBdisc_S, LdgJetBdisc_B;
float LdgJetCombinedCvsL_S, LdgJetCombinedCvsL_B;
float LdgJetDeepCvsL_S, LdgJetDeepCvsL_B;
float LdgJetPtD_S, LdgJetPtD_B;
float LdgJetAxis2_S, LdgJetAxis2_B;
float LdgJetAxis1_S, LdgJetAxis1_B;
int LdgJetMult_S, LdgJetMult_B;
float LdgJetQGLikelihood_S, LdgJetQGLikelihood_B;
//float LdgJetPullMagnitude_S, LdgJetPullMagnitude_B;
  
// Subleading jet from dijet system
float SubldgJetPt_S, SubldgJetPt_B;
float SubldgJetEta_S, SubldgJetEta_B;
float SubldgJetPhi_S, SubldgJetPhi_B;
float SubldgJetMass_S, SubldgJetMass_B;
float SubldgJetPFCharge_S, SubldgJetPFCharge_B;
float SubldgJetBdisc_S, SubldgJetBdisc_B;
float SubldgJetCombinedCvsL_S, SubldgJetCombinedCvsL_B;
float SubldgJetDeepCvsL_S, SubldgJetDeepCvsL_B;
float SubldgJetPtD_S, SubldgJetPtD_B;
float SubldgJetAxis2_S, SubldgJetAxis2_B;
float SubldgJetAxis1_S, SubldgJetAxis1_B;
int SubldgJetMult_S, SubldgJetMult_B;
float SubldgJetQGLikelihood_S, SubldgJetQGLikelihood_B;
//float SubldgJetPullMagnitude_S, LdgJetPullMagnitude_B;
  
// b-jet 
float bjetPt_S, bjetPt_B;
float bjetEta_S, bjetEta_B;
float bjetPhi_S, bjetPhi_B;
float bjetBdisc_S, bjetBdisc_B;
float bjetMass_S, bjetMass_B;
float bjetQGLikelihood_S, bjetQGLikelihood_B;
float bjetCombinedCvsL_S, bjetCombinedCvsL_B;
float bjetDeepCvsL_S, bjetDeepCvsL_B;
float bjetPtD_S, bjetPtD_B;
float bjetAxis2_S, bjetAxis2_B;
float bjetAxis1_S, bjetAxis1_B;
int bjetMult_S, bjetMult_B;
float bjetPFCharge_S, bjetPFCharge_B;
float bjetLdgJetMass_S, bjetLdgJetMass_B;
float bjetSubldgJetMass_S, bjetSubldgJetMass_B;
  
// Others (Soft-drop, distances, pull variables, etc)
float SoftDrop_n2_S, SoftDrop_n2_B;
//float PullAngleJ1J2_S, PullAngleJ1J2_B;
//float PullAngleJ2J1_S, PullAngleJ2J1_B;
  
float DEtaJ1withJ2_S, DEtaJ1withJ2_B;
float DEtaJ1withBJet_S, DEtaJ1withBJet_B;
float DEtaJ2withBJet_S, DEtaJ2withBJet_B;
float DEtaDijetwithBJet_S, DEtaDijetwithBJet_B;
float DEtaJ1BJetwithJ2_S, DEtaJ1BJetwithJ2_B;
float DEtaJ2BJetwithJ1_S, DEtaJ2BJetwithJ1_B;
float DPhiJ1withJ2_S, DPhiJ1withJ2_B;
float DPhiJ1withBJet_S, DPhiJ1withBJet_B;
float DPhiJ2withBJet_S, DPhiJ2withBJet_B;
float DPhiDijetwithBJet_S, DPhiDijetwithBJet_B;
float DPhiJ1BJetwithJ2_S, DPhiJ1BJetwithJ2_B;
float DPhiJ2BJetwithJ1_S, DPhiJ2BJetwithJ1_B;
float DRJ1withJ2_S, DRJ1withJ2_B;
float DRJ1withBJet_S, DRJ1withBJet_B;
float DRJ2withBJet_S, DRJ2withBJet_B;
float DRDijetwithBJet_S, DRDijetwithBJet_B;
float DRJ1BJetwithJ2_S, DRJ1BJetwithJ2_B;
float DRJ2BJetwithJ1_S, DRJ2BJetwithJ1_B;
  
float dijetMassOverTrijetMass_S, dijetMassOverTrijetMass_B;
