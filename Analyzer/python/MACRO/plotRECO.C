TString dir = "/afs/cern.ch/work/d/dalfonso/HGCnose/CMSSW_11_1_0_pre7_PRhfnticl/src/HGCnoseUtils/Analyzer/python/"; 

void plotCLsize() {

  gStyle->SetOptStat(00000000);  

  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_MuonPt20_kappa5.root",dir.Data()));

  TH1F* h_cl1 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L1");
  TH1F* h_cl2 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L2");
  TH1F* h_cl3 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L3");
  TH1F* h_cl4 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L4");
  TH1F* h_cl5 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L5");
  TH1F* h_cl6 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L6");
  TH1F* h_cl7 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L7");
  TH1F* h_cl8 = (TH1F*) f_reco->Get("comparisonPlots/ClusterSize_L8");

  TCanvas *c1 = new TCanvas("c", "c", 700, 600);

  //  h_cl1->SetTitle("#gamma p_{T} = 20 GeV |#eta|=3.5 ");
  h_cl1->SetTitle("L1");
  h_cl1->GetXaxis()->SetTitle("Cluster size");
  h_cl1->SetMaximum(800);
  h_cl1->SetLineWidth(3);
  h_cl1->SetLineColor(1);
  h_cl1->Draw("hist");

  h_cl2->SetTitle("L2");
  h_cl2->SetLineColor(2);
  h_cl2->SetLineWidth(3);
  h_cl2->Draw("hist same");

  h_cl3->SetTitle("L3");
  h_cl3->SetLineColor(3);
  h_cl3->SetLineWidth(3);
  h_cl3->Draw("hist same");

  h_cl4->SetTitle("L4");
  h_cl4->SetLineColor(4);
  h_cl4->SetLineWidth(3);
  h_cl4->Draw("hist same");

  h_cl5->SetTitle("L5");
  h_cl5->SetLineColor(5);
  h_cl5->SetLineWidth(3);
  h_cl5->Draw("hist same");

  h_cl6->SetTitle("L6");
  h_cl6->SetLineColor(6);
  h_cl6->SetLineWidth(3);
  h_cl6->Draw("hist same");

  h_cl7->SetTitle("L7");
  h_cl7->SetLineColor(7);
  h_cl7->SetLineWidth(3);
  h_cl7->Draw("hist same");

  h_cl8->SetTitle("L8");
  h_cl8->SetLineColor(8);
  h_cl8->SetLineWidth(3);
  h_cl8->Draw("hist same");

  c1->BuildLegend();
  c1->SaveAs(Form("~/www/HGCAL/TICLPlot/Response_Muon_ClusterSize.png"));

}

void plotResp() {

  gStyle->SetOptStat(00000000);  
  TFile *f_sim = new TFile(Form("%s/plots_simHits_D47_GammaPt4.root",dir.Data()));
  //  TFile *f_sim = new TFile(Form("%s/plots_simHits_D47_GammaPt20.root",dir.Data()));
  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_PionPt20.root",dir.Data()));
  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_GammaPt20.root",dir.Data()));
  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_GammaPt20_fixPR.root",dir.Data()));
  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_GammaPt4.root",dir.Data()));
  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_GammaPt20_minCL4_miss3.root",dir.Data()));

  TFile *f_study = new TFile(Form("%s/plots_simHits_D47_GammaPt20_th05MIP.root",dir.Data()));
  TFile *f_CLstudy = new TFile(Form("%s/plots_simHits_D47_GammaPt20_LCecut0.root",dir.Data()));
  TFile *f_CLstudy2 = new TFile(Form("%s/plots_recHits_D47_GammaPt20_kappa5.root",dir.Data()));

  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_GammaPt20_Plus2Layer.root",dir.Data()));

  //  TFile *f_sim = new TFile(Form("%s/plots_simHits_D47_PionPt20.root",dir.Data()));
  //  TFile *f_reco = new TFile(Form("%s/plots_recHits_D47_PionPt20.root",dir.Data()));

  TH1F* h_simH = (TH1F*) f_sim->Get("comparisonPlots/hRecHitResponse");
  TH1F* h_studySim = (TH1F*) f_study->Get("comparisonPlots/hRecHitResponse");
  TH1F* h_recH = (TH1F*) f_reco->Get("comparisonPlots/hRecHitResponse"); 
  TH1F* h_recC = (TH1F*) f_reco->Get("comparisonPlots/Cluster_Response"); 
  TH1F* h_recT = (TH1F*) f_reco->Get("comparisonPlots/Trackster_Response"); 
  TH1F* h_CLstudy = (TH1F*) f_CLstudy->Get("comparisonPlots/Cluster_Response");
  TH1F* h_CLstudy2 = (TH1F*) f_CLstudy2->Get("comparisonPlots/Cluster_Response");

  TCanvas *c1 = new TCanvas("c", "c", 700, 600);

  h_simH->SetTitle("#gamma p_{T}=4 GeV |#eta|=3.5");
  //  h_simH->SetTitle("#gamma p_{T}=20 GeV |#eta|=3.5");
  //  if(h_simH) h_simH->SetTitle("#pi p_{T}=20 GeV |#eta|=3.5");

  if(h_simH) h_simH->GetXaxis()->SetTitle("E^{reco} / E^{caloParticle} ");
  if(h_simH) h_simH->SetMaximum(500);
  //  if(h_simH) h_simH->SetMaximum(400);
  //  if(h_simH) h_simH->GetXaxis()->SetRangeUser(0.4,1.5);
  if(h_simH) h_simH->GetXaxis()->SetRangeUser(0.,1.5);
  if(h_simH) h_simH->SetLineColor(kRed);
  if(h_simH) h_simH->SetLineWidth(3);
  if(h_simH) h_simH->Draw("hist");
  /*
  if(h_studySim) h_studySim->SetLineColor(kOrange+1);
  if(h_studySim) h_studySim->SetLineWidth(3);
  if(h_studySim) h_studySim->Draw("hist sames");
  */
  if(h_recH) h_recH->SetLineWidth(3);
  if(h_recH) h_recH->SetLineColor(kBlue);
  if(h_recH) h_recH->Draw("hist sames");
  if(h_recC) h_recC->SetLineWidth(3);
  if(h_recC) h_recC->SetLineColor(kGreen+1);
  if(h_recC) h_recC->Draw("hist sames");
  /*
  if(h_CLstudy) h_CLstudy->SetLineColor(kMagenta+1);
  if(h_CLstudy) h_CLstudy->SetLineWidth(3);
  if(h_CLstudy) h_CLstudy->Draw("hist sames");
  if(h_CLstudy2) h_CLstudy2->SetLineColor(kMagenta-9);
  if(h_CLstudy2) h_CLstudy2->SetLineWidth(3);
  if(h_CLstudy2) h_CLstudy2->Draw("hist sames");
  //  if(h_recT) h_recT->SetLineWidth(3);
  //  if(h_recT) h_recT->SetLineColor(kMagenta);
  //  if(h_recT) h_recT->Draw("hist sames");
  */

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextSize(0.035);
  tex->SetTextColor(kRed);
  if(h_simH) tex->DrawLatex(0.65,0.9-0.05,Form("simHit #mu = %.2f", h_simH->GetMean()));
  //  tex->SetTextColor(kOrange+1);
  //  if(h_studySim) tex->DrawLatex(0.6,0.85-0.05,Form("simHit (>0.5MIP) #mu = %.2f", h_studySim->GetMean()));
  ////
  tex->SetTextColor(kBlue);
  if(h_recH) tex->DrawLatex(0.65,0.8-0.05,Form("recHit #mu = %.2f", h_recH->GetMean()));
  tex->SetTextColor(kGreen+1);
  if(h_recC) tex->DrawLatex(0.65,0.75-0.05,Form("layerCl #mu = %.2f", h_recC->GetMean()));
  tex->SetTextColor(kMagenta+1);
  /*
  if(h_CLstudy) tex->DrawLatex(0.6,0.7-0.05,Form("layerCl(ecut=0) #mu = %.2f", h_CLstudy->GetMean()));
  tex->SetTextColor(kMagenta-9);
  if(h_CLstudy2) tex->DrawLatex(0.6,0.65-0.05,Form("layerCl(kappa=5) #mu = %.2f", h_CLstudy2->GetMean()));
  */
  //  tex->SetTextColor(kMagenta);
  //  if(h_recT) tex->DrawLatex(0.65,0.7-0.05,Form("tracksterEM"));
  //                                                                                                                                                                                            
  //  c1->SaveAs(Form("~/www/HGCAL/TICLPlot/Response.png"));
  //  c1->SaveAs(Form("~/www/HGCAL/TICLPlot/Response_Pion.png"));
  c1->SaveAs(Form("~/www/HGCAL/TICLPlot/Response_Gamma_StudyResp_pt4.png"));

}


void plotRECO() {

  //  plotCLsize();
  plotResp();
  return;

}
