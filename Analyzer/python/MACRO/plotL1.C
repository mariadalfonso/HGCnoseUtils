void plotH3D(TString histoname1, TString histoname2,float rangeMin=0.f, float rangeMax=0.f) {

  TString dir = "/afs/cern.ch/work/d/dalfonso/HGCnose/CMSSW_11_0_0_pre8_L1/src/HGCnoseUtils/Analyzer/python/"; 
  //  TString dir = "/afs/cern.ch/work/d/dalfonso/HGCnose/CMSSW_11_0_X_2019-09-03-1100_L1/src/HGCnoseUtils/Analyzer/python/";

  gStyle->SetOptStat(00000000);  
  TFile *f_th_g = new TFile(Form("%s/histo_gamma_Pt20_3D_STC.root",dir.Data()));
  TFile *f_th_p = new TFile(Form("%s/test_pion_Pt20_3D_STC.root",dir.Data()));
  //  TFile *f_th_k = new TFile(Form("%s/test_K0L_Pt20.root",dir.Data()));

  TH1F* h_th_g_lead = (TH1F*) f_th_g->Get(Form("comparisonPlots/%s",histoname1.Data()));
  TH1F* h_th_p_lead = (TH1F*) f_th_p->Get(Form("comparisonPlots/%s",histoname1.Data())); 
  //  TH1F* h_th_k_lead = (TH1F*) f_th_k->Get(Form("comparisonPlots/%s",histoname1.Data())); 

  TH1F* h_th_g = (TH1F*) f_th_g->Get(Form("comparisonPlots/%s",histoname2.Data()));
  TH1F* h_th_p = (TH1F*) f_th_p->Get(Form("comparisonPlots/%s",histoname2.Data())); 
  //  TH1F* h_th_k = (TH1F*) f_th_k->Get(Form("comparisonPlots/%s",histoname2.Data())); 

  TCanvas *c1 = new TCanvas("c", "c", 700, 600);

  //
  h_th_g_lead->SetTitle("histoMaxVariableDR_C3d");
  h_th_g_lead->GetXaxis()->SetTitle("p_{T}^{L1} / p_{T}^{Gen} ");
  h_th_g_lead->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
  h_th_g_lead->SetLineWidth(3);
  h_th_g_lead->SetLineColor(kBlue);
  //  h_th_g_lead->SetLineColor(kCyan);
  h_th_g_lead->DrawNormalized("hist");
  //
  h_th_p_lead->SetLineWidth(3);
  h_th_p_lead->SetLineColor(kRed);
  //  h_th_p_lead->SetLineColor(kOrange);
  h_th_p_lead->DrawNormalized("hist sames");
  //
  //  h_th_k_lead->SetLineWidth(3);
  //  h_th_k_lead->SetLineColor(kGreen+1);
  //  h_th_k_lead->SetLineColor(kGreen);
  //  h_th_k_lead->DrawNormalized("hist sames");

  //
  if(h_th_g) h_th_g->SetLineWidth(3);
  if(h_th_g) h_th_g->SetLineColor(4);
  if(h_th_g) h_th_g->DrawNormalized("hist");
  //
  if(h_th_p) h_th_p->SetLineWidth(3);
  if(h_th_p) h_th_p->SetLineColor(kRed);
  if(h_th_p) h_th_p->DrawNormalized("hist sames");
  //
  //  if(h_th_k) h_th_k->SetLineWidth(3);
  //  if(h_th_k) h_th_k->SetLineColor(kGreen+1);
  //  if(h_th_k) h_th_k->DrawNormalized("hist sames");


  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);                                                                                                                               
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->DrawLatex(0.6,0.75,"#gamma, leading 3D clusters");
  //  tex->DrawLatex(0.6,0.87,"Threshould 1.5 GeV");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.6,0.70,"#pi^{#pm}, leading 3D clusters");
  //  tex->DrawLatex(0.6,0.75,"Super TC");
  /*
  tex->SetTextColor(kCyan);
  tex->DrawLatex(0.6,0.75,"#gamma, all cluster");
  //  tex->DrawLatex(0.6,0.87,"Threshould 1.5 GeV");
  tex->SetTextColor(kOrange);
  tex->DrawLatex(0.6,0.70,"#pi^{#pm}, all cluster");
  //  tex->DrawLatex(0.6,0.75,"Super TC");
  */

  c1->SaveAs(Form("~/www/HGCAL/L1Plot/%s.png",histoname1.Data()));

}


void plotH2D(TString histoname, float rangeMin=0.f, float rangeMax=0.f) {
  
  //  TString dir2 = "/afs/cern.ch/work/d/dalfonso/HGCnose/CMSSW_11_0_X_2019-09-03-1100_L1/src/HGCnoseUtils/Analyzer/python/";
  TString dir = "/afs/cern.ch/work/d/dalfonso/HGCnose/CMSSW_11_0_0_pre8_L1/src/HGCnoseUtils/Analyzer/python/"; 

  gStyle->SetOptStat(00000000);  
  TFile *f_th_g = new TFile(Form("%s/histo_gamma_Pt20.root",dir.Data()));
  TFile *f_th_pi = new TFile(Form("%s/histo_pion_Pt20.root",dir.Data()));
  TFile *f_th_k = new TFile(Form("%s/histo_K0L_Pt20.root",dir.Data()));
  //  TFile *f_th_pi = new TFile(Form("%s/test_L1_PionPt20_checkScale.root",dir2.Data()));
  //  TFile *f_tc_g = new TFile(Form("../test_L1_PhotonPt20_STC.root",dir.Data()));
  TFile *f_tc_g = new TFile(Form("%s/histo_gamma_Pt20_STC.root",dir.Data()));
  TFile *f_tc_pi = new TFile(Form("%s/histo_pion_Pt20_STC.root",dir.Data()));
  TFile *f_tc_k = new TFile(Form("%s/histo_K0L_Pt20_STC.root",dir.Data()));

  //  TFile *f_tc = new TFile(Form("../test_L1_PhotonPt20.root"));
  //  TFile *f_th = new TFile(Form("../test_Photon10_threshold.root"));
  //  TFile *f_tc = new TFile(Form("../test_Photon10_superTC.root"));

  TH1F* h_th_g = (TH1F*) f_th_g->Get(Form("comparisonPlots/%s",histoname.Data()));
  TH1F* h_tc_g = (TH1F*) f_tc_g->Get(Form("comparisonPlots/%s",histoname.Data()));
  TH1F* h_th_pi = (TH1F*) f_th_pi->Get(Form("comparisonPlots/%s",histoname.Data()));
  TH1F* h_tc_pi = (TH1F*) f_tc_pi->Get(Form("comparisonPlots/%s",histoname.Data()));

  f_th_g->ls();

  TCanvas *c1 = new TCanvas("c", "c", 700, 600);

  h_th_g->SetTitle("Concentrator");
  if(histoname.Contains("hClusterLayer")) {
    h_th_g->GetXaxis()->SetTitle("layer # ");
  } else  {
    h_th_g->GetXaxis()->SetTitle("p_{T}^{L1} / p_{T}^{Gen} ");
  }
  h_th_g->GetXaxis()->SetRangeUser(rangeMin,rangeMax);
  h_th_g->SetLineWidth(3);
  h_th_g->SetLineColor(4);
  h_th_g->DrawNormalized("hist");
  //
  h_th_pi->SetLineWidth(3);
  h_th_pi->SetLineColor(kRed);
  h_th_pi->DrawNormalized("hist sames");

  //  h_tc_g->SetLineStyle(3);
  if(h_tc_g) h_tc_g->SetLineWidth(3);
  if(h_tc_g) h_tc_g->SetLineColor(38);
  if(h_tc_g) h_tc_g->DrawNormalized("hist sames");
  //
  //  h_tc_pi->SetLineStyle(3);
  if(h_tc_pi) h_tc_pi->SetLineWidth(3);
  if(h_tc_pi) h_tc_pi->SetLineColor(kOrange);
  if(h_tc_pi) h_tc_pi->DrawNormalized("hist sames");

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);                                                                                                                               
  tex->SetTextSize(0.05);
  tex->SetTextColor(4);
  tex->DrawLatex(0.6,0.75,"#gamma - th 1.35 #MIPsT");
  //  tex->DrawLatex(0.6,0.87,"Threshould 1.5 GeV");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.6,0.70,"#pi^{#pm} - th 1.35 #MIPsT");
  //  tex->DrawLatex(0.6,0.75,"Super TC");
  tex->SetTextColor(38);
  tex->DrawLatex(0.6,0.65,"#gamma - STC");
  //  tex->DrawLatex(0.6,0.87,"Threshould 1.5 GeV");
  tex->SetTextColor(kOrange);
  tex->DrawLatex(0.6,0.60,"#pi^{#pm} - STC");
  //  tex->DrawLatex(0.6,0.75,"Super TC");

  
  c1->SaveAs(Form("~/www/HGCAL/L1Plot/%s.png",histoname.Data()));

}

void plotL1() {

  plotH3D("hMaxClResponse","hClusterResponse");
  plotH3D("hClusterConstituents","hClusterConstituents");


  //  plotH2D("hClusterLayer",0,10.);

  return;
  //  plotH2D("hClusterResponse",0,1.5);
  plotH2D("hResponse",0,1.5);

  /**/
  plotH2D("hClSigmaEtaEtaMax");
  plotH2D("hClSigmaEtaEtaTot");
  /**/
  plotH2D("hClSigmaPhiPhiMax");
  plotH2D("hClSigmaPhiPhiTot");
  /**/
  plotH2D("hClSigmaRRMax");
  plotH2D("hClSigmaRRTot");
  /**/
  plotH2D("hClSigmaZZ");

  return;


}
