void plotH(TString histoname) {

  gStyle->SetOptStat(00000000);  

  TFile *f_th = new TFile(Form("../test_Photon10_threshold.root"));
  TFile *f_tc = new TFile(Form("../test_Photon10_superTC.root"));

  TH1F* h_th = (TH1F*) f_th->Get(Form("comparisonPlots/%s",histoname.Data()));
  TH1F* h_tc = (TH1F*) f_tc->Get(Form("comparisonPlots/%s",histoname.Data()));

  TCanvas c;

  h_th->GetXaxis()->SetTitle("Reco/Gen  p_{T} ");
  h_th->SetLineWidth(3);
  h_th->SetMaximum(200);
  h_th->Draw("hist");
  //
  h_tc->SetLineWidth(3);
  h_tc->SetLineColor(2);
  //
  h_tc->Draw("hist sames");

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);                                                                                                                               tex->SetTextSize(0.08);
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.6,0.87,"Threshould 1.5 GeV");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.6,0.75,"Super TC");

  
  c.SaveAs(Form("~/www/HGCAL/L1Plot/%s.png",histoname.Data()));

}

void plotL1() {

  plotH("hClusterResponse");
  /**/
  plotH("hClSigmaEtaEtaMax");
  plotH("hClSigmaEtaEtaTot");
  /**/
  plotH("hClSigmaPhiPhiMax");
  plotH("hClSigmaPhiPhiTot");
  /**/
  plotH("hClSigmaRRMax");
  plotH("hClSigmaRRTot");
  /**/
  plotH("hClSigmaZZ");


}
