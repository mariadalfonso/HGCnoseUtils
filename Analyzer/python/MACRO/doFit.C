#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"



void curveTime() {

  TMultiGraph *mg = new TMultiGraph();

  TCanvas c;
  c.cd();

  double E[] = { 1, 2, 3, 5, 10, 20};
  double sigma[] = { 0.0488089, 0.0404168, 0.0332543, 0.0264095, 0.0195331, 0.015022 };

  double Egamma[] = { 1, 2, 3, 5, 10, 20};
  double sigmagamma[] = { 0.0291707, 0.0221313, 0.0192523, 0.0166872 , 0.0130426 , 0.0104512};

  /// this is the NOSE K0L
  TGraphErrors *gr = new TGraphErrors(6,E,sigma);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(20);

  /// this is the NOSE GAMMA
  TGraphErrors *gr2 = new TGraphErrors(6,Egamma,sigmagamma);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(21);

  mg->Add(gr2);
  mg->Add(gr);
  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("#sigma_{time} (ns)");
  mg->GetXaxis()->SetTitle("p_{T} (GeV)");
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetXaxis()->SetTitleOffset(0.9);
  mg->GetYaxis()->SetTitleOffset(1);
  mg->GetYaxis()->SetRangeUser(0., 0.1);

  c.Update();

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);
  tex->SetTextSize(0.08);
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.6,0.87,"K_{L}^{0}");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.6,0.75,"#gamma");
  //  TString resultHF = Form("#frac{#sigma(E)}{E} = %.2f%% + #frac{%.2f%%}{#sqrt{E}}", 100*func2->GetParameter(1), 100*func2->GetParameter(0));

  /*
    Energy=5 mean=0.00134558 sigma=0.0641475
    Energy=10 mean=0.00123703 sigma=0.0636641
    Energy=20 mean=0.000472522 sigma=0.0603595
    [dalfonso@lxplus752 MACRO]$ cat fitTime.txt | grep Energy
    Energy=5 mean=0.00225554 sigma=0.074664
    Energy=10 mean=0.00090718 sigma=0.0757135
    Energy=20 mean=0.000214311 sigma=0.0749297
  */

  c.SaveAs(Form("~/www/HGCAL/fitPlot/Curve_Time.png"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/Curve_Time.pdf"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/Curve_Time.C"));

  }

void curveH() {

  //  gStyle->SetOptFit();
  TMultiGraph *mg = new TMultiGraph();

  TCanvas c;
  c.cd();

  double E[] = {10, 20, 40, 80, 100, 150};
  double Eerr[] = {0, 0, 0, 0, 0, 0};
  double sigma[] = { 0.469791, 0.405783, 0.30167, 0.23561, 0.218876, 0.192147 };
  double mean[] = { 0.714755, 0.780795, 0.848058, 0.889448, 0.892589, 0.909765 };
  double sigmaNorm[] = { 0.469791/0.714755,
			 0.405783/0.780795,
			 0.30167/0.848058,
			 0.23561/0.889448, 
			 0.218876/0.892589,
			 0.192147/0.909765
  };

  /// this is the NOSE
  TGraphErrors *gr = new TGraphErrors(6,E,sigmaNorm);
  //  TGraphErrors *gr = new TGraphErrors(6,E,sigma);
  gr->SetMarkerColor(kBlue);
  gr->SetMarkerStyle(20);
  TF1 *func = new TF1("func","[1] + [0]*TMath::Power(x,-0.5)",0,200);
  func->SetParNames("stochasticHGCNose", "constantHGCNose");
  func->SetParameters(0.30, 0.01);
  func->SetLineColor(kBlue);

  mg->Add(gr);
  gr->Fit("func","r");

  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("#sigma(E)/E");
  mg->GetXaxis()->SetTitle("E (GeV)");
  mg->GetYaxis()->SetRangeUser(0.,1.);

  c.Update();

  double EHF[] = {40, 60, 80, 100, 150, 200};
  double EerrHF[] = {0, 0, 0, 0 , 0, 0};
  double sigmaHF[] = {0.403678, 0.335126, 0.310582, 0.313421, 0.285184, 0.242201 };
  double meanHF[] = { 0.852772, 0.907938, 0.929439, 0.942506, 0.978554, 1.01269 };
  double sigmaNormHF[] = { 0.403678/0.852772,
			   0.335126/0.907938,
			   0.310582/0.929439,
			   0.313421/0.942506,
			   0.285184/0.978554,
			   0.242201/1.01269
  };

  /// this is the HF
  TGraphErrors *gr2 = new TGraphErrors(6,EHF,sigmaNormHF);
  //  TGraphErrors *gr2 = new TGraphErrors(6,EHF,sigmaHF);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(21);
  TF1 *func2 = new TF1("func2","[1] + [0]*TMath::Power(x,-0.5)",0,200);
  func2->SetParNames("stochasticHF", "constantHF");
  func2->SetParameters(0.30, 0.01);
  func2->SetLineColor(kRed);

  mg->Add(gr2);
  gr2->Fit("func2","r");

  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("#sigma(E)/E" );
  mg->GetXaxis()->SetTitle("E (GeV)" );

  c.Update();


  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);
  tex->SetTextSize(0.0425);
  tex->DrawLatex(0.45,0.95,"#pi");
  tex->SetTextColor(kRed);
  TString resultHF = Form("#frac{#sigma(E)}{E} = %.2f%% + #frac{%.2f%%}{#sqrt{E}}", 100*func2->GetParameter(1), 100*func2->GetParameter(0));
  tex->DrawLatex(0.35,0.66+0.15,Form(" HF:    %s", resultHF.Data()));
  tex->SetTextColor(kBlue);
  TString resultHGC = Form("#frac{#sigma(E)}{E} = %.2f%% + #frac{%.2f%%}{#sqrt{E}}", 100*func->GetParameter(1), 100*func->GetParameter(0));
  tex->DrawLatex(0.35,0.56+0.15,Form(" HF+HGCnose:    %s", resultHGC.Data()));

  //  tex->DrawLatex(0.5,0.66+0.1,Form(" HF: constant=%.2f%%", 100* func2->GetParameter(1)));
  //  tex->DrawLatex(0.5,0.6+0.1,Form(" HF: stochastic %.2f%%", 100*func2->GetParameter(0)));
  //  tex->SetTextColor(kBlue);
  //  tex->DrawLatex(0.5,0.55+0.1,Form(" HGCnose: constant=%.2f%%", 100* func->GetParameter(1)));
  //  tex->DrawLatex(0.5,0.5+0.1,Form(" HGCnose: stochastic %.2f%%", 100*func->GetParameter(0)));


  c.SaveAs(Form("~/www/HGCAL/fitPlot/pion_Nose.png"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/pion_Nose.pdf"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/pion_Nose.C"));


}


void curveE() {

  TCanvas c;
  c.cd();

  //  gStyle->SetOptFit();

  TMultiGraph *mg = new TMultiGraph();

  double E[] = {10, 20, 40, 80, 100, 150 };
  //  double E[] = {TMath::Power(10,-0.5), TMath::Power(20,-0.5), TMath::Power(40,-0.5), TMath::Power(80,-0.5), TMath::Power(100,-0.5), TMath::Power(150,-0.5) };

  double Eerr[] = {0, 0, 0, 0, 0, 0 };
  double mean[] = {0.657502 , 0.666607, 0.678935, 0.688737, 0.691885, 0.697249};
  double sigma[] = { 0.150206, 0.107712, 0.0759965, 0.0559665, 0.0517583, 0.044394};
  double sigmaNorm[] = {0.150206/0.657502, 
			0.107712/0.666607, 
			0.0759965/0.678935,
			0.0559665/0.688737,
			0.0517583/0.69188,
			0.044394/0.697249
                        };

  //  double EHF[] = {10, 20, 40, 60, 80, 100, 150 };
  double EHF[] =  { 40, 60, 80, 100, 150 };
  //  double EHF[] = {TMath::Power(10,-0.5), TMath::Power(20,-0.5), TMath::Power(40,-0.5), TMath::Power(60,-0.5), TMath::Power(80,-0.5), TMath::Power(100,-0.5), TMath::Power(150,-0.5) };
  double meanHF[] = { 0.716012, 0.738176 , 0.825067, 0.852763, 0.865134, 0.866213, 0.890558 };
  double sigmaHF[] = { 0.686558, 0.507386, 0.323526, 0.277076, 0.254026, 0.243544, 0.224338};
  

  double sigmaNormHF[] = {//0.686558/0.716012,
			  //0.507386/0.738176,
			  0.323526/0.825067,
			  0.277076/0.852763,
			  0.254026/0.865134,
			  0.243544/0.866213,
			  0.224338/0.890558
  };

  /*
  //redefined
  double sigmaNormHF[] = {0.583468/0.730908,
			  0.539367/0.648432,
			  0.309463/0.787747,
			  0.297712/0.806747,
			  0.246464/0.84729,
			  0.236622/0.833874,
			  0.197621/0.863276
  };
  */

  /*
  double EHF[] = { 20, 40, 60, 80, 100, 150 };
  double meanHF[] = { 0.738176 , 0.825067, 0.852763, 0.865134, 0.866213, 0.890558 };
  double sigmaHF[] = { 0.507386, 0.323526, 0.277076, 0.254026, 0.243544, 0.224338};
  
  double sigmaNormHF[] = {//0.686558/0.716012,
			  0.507386/0.738176,
			  0.323526/0.825067,
			  0.277076/0.852763,
			  0.254026/0.865134,
			  0.243544/0.866213,
			  0.224338/0.890558
  };
  */

  /*
  Double_t * sigmaNorm;
  sigmaNorm[0] = sigma[0]/mean[0];
  sigmaNorm[1] = sigma[1]/mean[1];
  sigmaNorm[2] = sigma[2]/mean[2];
  sigmaNorm[3] = sigma[3]/mean[3];
  sigmaNorm[4] = sigma[4]/mean[4];
  sigmaNorm[5] = sigma[5]/mean[5];
  */

  /// this is the NOSE
  TGraphErrors *gr1 = new TGraphErrors(6,E,sigmaNorm);
  //  TGraphErrors *gr1 = new TGraphErrors(6,E,sigma);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  TF1 *func = new TF1("func","[1] + [0]*TMath::Power(x,-0.5)",0,200);
  //  TF1 *func = new TF1("func","[1] + [0]*x",0,200);
  func->SetParNames("stochasticHGCNose", "constantHGCNose");
  func->SetParameters(0.70, 0.03);
  func->SetParLimits(1,0.01, 0.10);
  func->SetLineColor(4);

  mg->Add(gr1);
  //  gr1->Fit("func","r");
  gr1->Fit("func");

  c.Update();

  /// this is the HF
  TGraphErrors *gr2 = new TGraphErrors(5,EHF,sigmaNormHF);
  //  TGraphErrors *gr2 = new TGraphErrors(7,EHF,sigmaHF);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);
  TF1 *func2 = new TF1("func2","[1] + [0]*TMath::Power(x,-0.5)",0,200);
  //  TF1 *func2 = new TF1("func2","[1] + [0]*x",0,200);
  func2->SetParNames("stochasticHF", "constantHF");
  func2->SetParameters(3, 0.03);
  func2->SetParLimits(1,0.05, 0.15);
  func2->SetLineColor(kRed);

  mg->Add(gr2);
  gr2->Fit("func2","r");
  //  mg->Fit("gauss");

  mg->Draw("ap");
  mg->GetYaxis()->SetTitle("#sigma(E)/E" );
  mg->GetXaxis()->SetTitle("E (GeV)" );
  mg->GetYaxis()->SetRangeUser(0.,1.);

  TString resultNose = Form("#frac{#sigma(E)}{E} = %.1f%% + #frac{%.1f%%}{#sqrt{E}}", 100*func->GetParameter(1), 100*func->GetParameter(0));
  TString resultHF = Form("#frac{#sigma(E)}{E} = %.1f%% + #frac{%.1f%%}{#sqrt{E}}", 100*func2->GetParameter(1), 100*func2->GetParameter(0));

  //  TString result = Form("%.2f%% + #frac{%.2f%%}{#sqrt{E}}", 100*func->GetParameter(1), 100*func->GetParameter(0));
  //  mg->SetTitle(result.Data());
  //  TString result = Form("%.2f%% + #sqrt E ", 100*func->GetParameter(1), 100*func->GetParameter(0));
  //  tex->SetLineWidth(2);
  //  tex->DrawLatex(0.6,0.6,Form("%.2f%% + #frac{%.2f%%}{#sqrt{E}}:", 100*func->GetParameter(1), 100*func->GetParameter(0)));


  c.Update();

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  //  tex->SetTextFont(61);
  tex->SetTextSize(0.0425);
  tex->DrawLatex(0.45,0.95,"#gamma");
  tex->SetTextColor(kRed);
  tex->DrawLatex(0.45,0.66+0.2,Form(" HF: %s", resultHF.Data()));
  tex->SetTextColor(kBlue);
  tex->DrawLatex(0.45,0.5+0.2,Form(" HGCnose: %s", resultNose.Data()));


  c.SaveAs(Form("~/www/HGCAL/fitPlot/gamma_Nose.png"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/gamma_Nose.pdf"));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/gamma_Nose.C"));

}

void doFitTime(int Energy) {

  gStyle->SetOptFit(1);

  TCanvas c;
  c.cd();
  c.SetLogy(1);
  /*
    -rw-r--r--. 1 dalfonso zh 54757 Apr 22 17:08 test_K0L_Pt5.root
    -rw-r--r--. 1 dalfonso zh 56199 Apr 22 17:12 test_K0L_Pt10.root
    -rw-r--r--. 1 dalfonso zh 57764 Apr 22 17:22 test_K0L_Pt20.root
  */

  //  -rw-r--r--. 1 dalfonso zh 30671 Apr 22 19:49 ../TIME/test_Timing_smallCone_Gamma_Pt100.root
  //  -rw-r--r--. 1 dalfonso zh 30671 Apr 22 20:55 ../TIME/test_Timing_Gamma_Pt100.root
  //  TFile *f = new TFile(Form("../TIME/test_Timing_smallCone_Gamma_Pt%d.root",Energy));

  //  TFile *f = new TFile(Form("../TIME/test_K0L_smallCone_Pt%d.root",Energy));

  TFile *f = new TFile(Form("../TIME2/test_K0L_Pt%d.root",Energy));
  //  TFile *f = new TFile(Form("../TIME2/test_gamma_Pt%d.root",Energy));
  TH1F* h1 = (TH1F*) f->Get("comparisonPlots/hParticleTime");

  h1->Draw();
  h1->SetXTitle("hitTime - clusterTime");
  TF1 func("func","gaus",h1->GetMean()-3*h1->GetRMS(),h1->GetMean()+3*h1->GetRMS());
  h1->Fit("func","r");

  std::cout << "Energy=" << Energy << " mean=" << h1->GetFunction("func")->GetParameter(1) << " sigma=" << h1->GetFunction("func")->GetParameter(2)<< std::endl;

  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitTimeGamma_smallCone_Pt_%d.png",Energy));
  //  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitTimeK0L_smallCone_Pt_%d.png",Energy));
  //  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitTimeK0L_smallCone_Pt_%d.png",Energy));


}

void doFitG(int Energy) {

  TCanvas c;
  c.cd();

  /*
  TString filename= Form("../RESOLUTION/test_Pion_Nose_E%d.root",Energy);

  std::cout << "filename" << filename.Data() << std::endl;
  TFile *f = new TFile(filename.Data());
  TH1F* h1 = (TH1F*) f->Get("comparisonPlots/hJetResponse");
  */
  /*
  TString filename= Form("../RESOLUTION/test_Pion_HF_E%d.root",Energy);

  std::cout << "filename" << filename.Data() << std::endl;
  TFile *f = new TFile(filename.Data());
  TH1F* h1 = (TH1F*) f->Get("comparisonPlots/hJetResponseHF");
  */

  TFile *f = new TFile(Form("../RESOLUTION/test_Gamma_HF_E%d.root",Energy));
  TH1F* h1 = (TH1F*) f->Get("comparisonPlots/hJetResponseHF");

  //  TFile *f = new TFile(Form("../RESOLUTION/test_Photons_Nose_E%d.root",Energy));
  //  TH1F* h1 = (TH1F*) f->Get("comparisonPlots/hJetResponseNose");

  h1->Draw();
  //  h1->Fit("gaus");
  TF1 func("func","gaus",h1->GetMean()-h1->GetRMS(),h1->GetMean()+h1->GetRMS());
  h1->Fit("func","r");

  std::cout << "Energy=" << Energy << " mean=" << h1->GetFunction("func")->GetParameter(1) << " sigma=" << h1->GetFunction("func")->GetParameter(2)<< std::endl;

  //  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitPionNose_%d.png",Energy));
  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitGammHF_%d.png",Energy));
  

}

void deriveWeight() {


  TCanvas c;
  c.cd();

  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  TFile *f1 = new TFile("SCAN/test_Pion_Nose_E10.root");
  TH1F* h1 = (TH1F*) f1->Get("comparisonPlots/hJetResponse");
  h1->Rebin(2);
  h1->Draw();
  h1->GetXaxis()->SetTitle("reco Pion Energy/ gen Pion Energy");
  h1->SetTitle("#pi 10 GeV");
  h1->SetMaximum(500.);
  //  h1->SetTitle("#pi 100 GeV");
  //  h1->SetMaximum(2*500.);
  h1->SetLineWidth(3);
  TF1 func1("func1","gaus",h1->GetMean()-h1->GetRMS(),h1->GetMean()+h1->GetRMS());
  //  TF1 func1("func1","gaus",h1->GetMean()-0.1, h1->GetMean()+0.1);
  func1.SetLineColor(1);
  h1->Fit("func1","r");

  TFile *f2 = new TFile("SCAN/test_Pion_Nose_E10_W08.root");
  TH1F* h2 = (TH1F*) f2->Get("comparisonPlots/hJetResponse");
  h2->Rebin(2);
  h2->SetLineColor(kRed);
  h2->SetLineWidth(3);
  h2->Draw("sames");
  TF1 func2("func2","gaus",h2->GetMean()-h2->GetRMS(),h2->GetMean()+h2->GetRMS());
  func2.SetLineColor(2);
  h2->Fit("func2","r");


  TFile *f3 = new TFile("SCAN/test_Pion_Nose_E10_W06.root");
  TH1F* h3 = (TH1F*) f3->Get("comparisonPlots/hJetResponse");
  h3->Rebin(2);
  h3->SetLineColor(kGreen+1);
  h3->SetLineWidth(3);
  h3->Draw("sames");
  TF1 func3("func3","gaus",h3->GetMean()-h3->GetRMS(),h3->GetMean()+h3->GetRMS());
  func3.SetLineColor(kGreen+1);
  h3->Fit("func3","r");


  TFile *f4 = new TFile("SCAN/test_Pion_Nose_E10_W12.root");
  TH1F* h4 = (TH1F*) f4->Get("comparisonPlots/hJetResponse");
  h4->Rebin(2);
  h4->SetLineColor(kBlue);
  h4->SetLineWidth(3);
  h4->Draw("sames");
  TF1 func4("func4","gaus",h4->GetMean()-h4->GetRMS(),h4->GetMean()+h4->GetRMS());
  func4.SetLineColor(kBlue);
  h4->Fit("func4","r");

  /*
  TFile *f5 = new TFile("test_Pion_Nose_E10_W14.root");
  TH1F* h5 = (TH1F*) f5->Get("comparisonPlots/hJetResponse");
  h5->SetLineColor(kViolet);
  h5->SetLineWidth(3);
  h5->Rebin(2);
  h5->Draw("sames");
  TF1 func5("func5","gaus",h5->GetMean()-h5->GetRMS(),h5->GetMean()+h5->GetRMS());
  func5.SetLineColor(kViolet);
  h5->Fit("func5","r");
  */

  TLatex *  tex = new TLatex();
  tex->SetNDC();
  tex->SetTextAlign(13);
  tex->SetTextSize(0.035);
  //  tex->DrawLatex(0.45,0.95,"#pi");
  //
  tex->SetTextColor(kBlack);
  double NormSigm=func1.GetParameter(2)/func1.GetParameter(1);
  TString result = Form("CF-E + CF-H + HF #sigma^{eff}= %.2f",NormSigm);
  tex->DrawLatex(0.5,0.85,Form("%s", result.Data()));
  //
  tex->SetTextColor(kRed);
  NormSigm=func2.GetParameter(2)/func2.GetParameter(1);
  result = Form("CF-E + 0.8*CF-H + 0.8*HF #sigma^{eff}= %.2f",NormSigm);
  tex->DrawLatex(0.5,0.8,Form("%s", result.Data()));
  //
  tex->SetTextColor(kGreen+1);
  NormSigm=func3.GetParameter(2)/func3.GetParameter(1);
  result = Form("CF-E + 0.6*CF-H + 0.6*HF #sigma^{eff}= %.2f",NormSigm);
  tex->DrawLatex(0.5,0.75,Form("%s", result.Data()));
  //
  tex->SetTextColor(kBlue+1);
  NormSigm=func4.GetParameter(2)/func4.GetParameter(1);
  result = Form("CF-E + 1.2*CF-H + 1.2*HF #sigma^{eff}= %.2f",NormSigm);
  tex->DrawLatex(0.5,0.7,Form("%s", result.Data()));

  c.SaveAs("~/www/HGCAL/fitPlot/fitPionNose_Weight_E10.png");
  //  c.SaveAs("~/www/HGCAL/fitPlot/fitPionNose_Weight_E100.png");
  //  c.SaveAs(Form("~/www/HGCAL/fitPlot/fitPionNose_Weight.png",Energy));

}


void doFit() {

  curveTime();
  return;

  //  doFitTime(1);
  //  doFitTime(2);
  doFitTime(3);
  //  doFitTime(5);
  //  doFitTime(10);
  //  doFitTime(20);

  return;


  /*
  doFitG(10);
  doFitG(20);
  doFitG(40);
  doFitG(60);
  doFitG(80);
  doFitG(100);
  doFitG(150);
  //  doFitG(200);
  */

  //  curveE();
 
  curveH();
  //  deriveWeight();

}
