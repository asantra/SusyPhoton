{

  gStyle->SetOptStat(0000000);
  //Mass points
  Double_t XCen[5];
  XCen[0] = 150;
  XCen[1] = 175;
  XCen[2] = 200;
  XCen[3] = 225;
  XCen[4] = 250;
  Double_t XCenErr[5];
  XCenErr[0]=0;
  XCenErr[1]=0;
  XCenErr[2]=0;
  XCenErr[3]=0;
  XCenErr[4]=0;

  //These are expected limit points
  Double_t YCen[5];
  YCen[0] = 0.0347;
  YCen[1] = 0.0174;
  YCen[2] = 0.0184;
  YCen[3] = 0.0150;
  YCen[4] = 0.0148;

  //Upper yellow band, 2sigma expected minus central value
  Double_t YCenErrYellowUp[5];
  YCenErrYellowUp[0] = 0.12-YCen[0]; 
  YCenErrYellowUp[1] = 0.0567-YCen[1]; 
  YCenErrYellowUp[2] = 0.0526-YCen[2];
  YCenErrYellowUp[3] = 0.0419-YCen[3];
  YCenErrYellowUp[4] = 0.0412-YCen[4];


  //Lower yellow band, central value - 2 sigma low.
  Double_t YCenErrYellowDown[5];
  YCenErrYellowDown[0] = YCen[0] - 0.0223; 
  YCenErrYellowDown[1] = YCen[1] - 0.0084; 
  YCenErrYellowDown[2] = YCen[2] - 0.0089;
  YCenErrYellowDown[3] = YCen[3] - 0.0128;
  YCenErrYellowDown[4] = YCen[4] - 0.0123;

  //Upper green band, 1sigma expected minus central value
  Double_t YCenErrGreenUp[5];
  YCenErrGreenUp[0] = 0.0596-YCen[0]; 
  YCenErrGreenUp[1] = 0.0317-YCen[1]; 
  YCenErrGreenUp[2] = 0.0300-YCen[2];
  YCenErrGreenUp[3] = 0.0222-YCen[3];
  YCenErrGreenUp[4] = 0.0261-YCen[4];

  //Lower green band, central value minus 1 sigma low
  Double_t YCenErrGreenDown[5];
  YCenErrGreenDown[0] = YCen[0] - 0.0254; 
  YCenErrGreenDown[1] = YCen[1] - 0.0146; 
  YCenErrGreenDown[2] = YCen[2] - 0.0156;
  YCenErrGreenDown[3] = YCen[3] - 0.0138;
  YCenErrGreenDown[4] = YCen[4] - 0.0133;


  //Observed limit
  Double_t ObservedY[5];
  ObservedY[0] = 0.0835;
  ObservedY[1] = 0.0431;
  ObservedY[2] = 0.04095;
  ObservedY[3] = 0.0368;
  ObservedY[4] = 0.0359;

  //Vectorlike Confinement sigma x br, NR
  Double_t VectorTheory[5];
  VectorTheory[0] = 0.200 * 0.551;
  VectorTheory[1] = 0.113 * 0.443;
  VectorTheory[2] = 0.0664* 0.372;
  VectorTheory[3] = 0.0417* 0.323;
  VectorTheory[4] = 0.0274* 0.291;

  //Vectorlike Confinement sigma x br, M=1.5 TeV
   Double_t VectorTheory15[5];
  VectorTheory15[0] = 0.146;
  VectorTheory15[1] = 0.0723;
  VectorTheory15[2] = 0.0408;
  VectorTheory15[3] = 0.0253;
  VectorTheory15[4] = 0.0173;

  TH2F *ax = new TH2F("ax","",1000,130,275,1000,0.005,0.15);
  ax->Draw();
  ax->GetXaxis()->SetTitle("M_{#tilde#pi^{0}} (GeV)");
  ax->GetYaxis()->SetTitle("#sigma x Br(#tilde#pi^{0}#rightarrow#gamma#gamma) (pb)");
  ax->GetYaxis()->SetTitleOffset(1.2);
  TGraphAsymmErrors *garpY = new TGraphAsymmErrors(5, XCen, YCen, XCenErr, XCenErr,YCenErrYellowDown, YCenErrYellowUp);
  garpY->SetFillColor(5);
  garpY->Draw("E3");

  TGraphAsymmErrors *garpG = new TGraphAsymmErrors(5, XCen, YCen, XCenErr, XCenErr,YCenErrGreenDown, YCenErrGreenUp);
  garpG->SetFillColor(8);
  garpG->Draw("SAMEE3");

  TGraphAsymmErrors *garpL = new TGraphAsymmErrors(5, XCen, YCen, XCenErr, XCenErr,XCenErr, XCenErr);
  garpL->SetLineWidth(2);
  garpL->SetLineStyle(2);
  garpL->Draw("SAMEC");
  
  TGraphAsymmErrors *garpO = new TGraphAsymmErrors(5, XCen, ObservedY, XCenErr, XCenErr, XCenErr, XCenErr);
  garpO->SetLineWidth(2);
  garpO->Draw("SAMEC");

  TGraphAsymmErrors *theo = new TGraphAsymmErrors(5, XCen, VectorTheory, XCenErr, XCenErr, XCenErr, XCenErr);
  theo->SetLineWidth(2);
  theo->SetLineColor(2);
  theo->Draw("SAMEC");

  TGraphAsymmErrors *theoR15 = new TGraphAsymmErrors(5, XCen, VectorTheory15, XCenErr, XCenErr, XCenErr, XCenErr);
  theoR15->SetLineWidth(2);
  theoR15->SetLineStyle(2);
  theoR15->SetLineColor(2);
  theoR15->Draw("SAMEC");

  TLegend *legend = new TLegend(0.4525862,0.5826271,0.8864943,0.8940678,NULL,"brNDC");
  legend->SetFillColor(0);
  legend->SetLineWidth(1);
  legend->SetTextFont(62);
  legend->AddEntry(garpO,"Observed Limit @95%CL","l");
  legend->AddEntry(garpL,"Expected Limit @95%CL","l");
  legend->AddEntry(garpG,"+/-1 #sigma expected","f");
  legend->AddEntry(garpY,"+/-2 #sigma expected","f");
  legend->AddEntry(theo,"Vectorlike Confinement, non-resonant","l");
  legend->AddEntry(theoR15,"Vectorlike Confinement, M_{#tilde#rho}=1.5 TeV","l");
  legend->SetBorderSize(0);
  legend->Draw();


}
