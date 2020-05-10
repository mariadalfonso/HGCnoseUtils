// -*- C++ -*-
//
// Package:    ComparisonPlots/HCALGPUAnalyzer
// Class:      HCALGPUAnalyzer
//
/**\class HCALGPUAnalyzer HCALGPUAnalyzer.cc ComparisonPlots/HCALGPUAnalyzer/plugins/HCALGPUAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mariarosaria D'Alfonso
//         Created:  Mon, 17 Dec 2018 16:22:58 GMT
//
//


// system include files                                                                                                                                                                                                                                               
#include <memory>
#include <string>
#include <map>
#include <iostream>
using namespace std;


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
// #include "DataFormats/TrackReco/interface/Track.h"
// #include "DataFormats/TrackReco/interface/TrackFwd.h"

//#include "FWCore/Framework/interface/Frameworkfwd.h"
//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/Framework/interface/Event.h"

//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include "DataFormats/Math/interface/deltaR.h"

// more HGC classes
#include "DataFormats/ForwardDetId/interface/ForwardSubdetector.h"
#include "DataFormats/ForwardDetId/interface/HFNoseDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCSiliconDetId.h"
#include "DataFormats/ForwardDetId/interface/HGCScintillatorDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"
#include "DataFormats/HGCDigi/interface/HGCDigiCollections.h"
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/CaloTest/interface/HGCalTestNumbering.h"
#include "SimDataFormats/CaloHit/interface/PCaloHit.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

//#include "CLHEP/Units/PhysicalConstants.h"
//#include "CLHEP/Units/SystemOfUnits.h"

#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

//#include "RecoParticleFlow/PFClusterProducer/plugins/SimMappers/ComputeClusterTime.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class GenAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit GenAnalyzer(const edm::ParameterSet&);
  ~GenAnalyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void analyzeDigi(edm::Handle<HGCalDigiCollection> digiNose, const math::XYZTLorentzVector & initialP4, const HGCalGeometry* geom);

  double analyzeHits(const std::vector<PCaloHit>&, const HGCRecHitCollection & recHitNose, const HFRecHitCollection& , const math::XYZTLorentzVector & initialP4 , const double & jetE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal, bool );

  void analyzeClusters(std::vector<reco::CaloCluster> recoNoseClusters,const math::XYZTLorentzVector & initialP4, const double & jetE, const HGCalGeometry* geom);

  std::vector<const reco::GenJet*> doVBSselection(edm::Handle<reco::GenParticleCollection> genParticles, edm::Handle<reco::GenJetCollection> genJets);

  void getSingle(edm::Handle<reco::GenParticleCollection> genParticles, std::vector<CaloParticle>, std::vector<reco::CaloCluster> recoNoseClusters, edm::Handle<HGCalDigiCollection> digiNose, const HGCRecHitCollection & recHitNose, const std::vector<PCaloHit>&, const HFRecHitCollection&, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  // ----------member data ---------------------------
  //  void ClearVariables();

  bool verbose = false;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTag_;
  edm::EDGetTokenT<reco::GenJetCollection> genJetsTag_;
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitTag_;
  edm::EDGetTokenT<edm::PCaloHitContainer> simHcalTag_;
  edm::EDGetTokenT<HGCalDigiCollection> digiNose_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitNose_;
  
  edm::EDGetTokenT<std::vector<SimCluster>> simClusterTag_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticleTag_;

  edm::EDGetTokenT<vector<reco::CaloCluster>> recoClusterTag_;
  edm::EDGetTokenT<vector<reco::CaloCluster>> recoNoseClusterTag_;

  edm::EDGetTokenT<HFRecHitCollection> tok_hf_;

  TH1F *hDRpi0;

  TH1F *hLepEta;
  TH1F *hLepPt;
  TH1F *hLepCharge;

  TH1F *hJetEta;
  TH1F *hJetEtaLarge;

  TH1F *hGenParticleEta;
  TH1F *hGenParticlePt;
  TH1F *hGenParticleE;

  TH1F *hJetDeltaRap;
  TH1F *hJetDiMass;

  TH1F *hJetPtMin;
  TH1F *hJetPtMax;

  TH1F *hEnergyE;
  TH1F *hEnergyH;
  TH1F *hEnergy;
  TH1F *hEnergyFrac;

  TH1F *hParticleTime;
  TH1F *hDigiTime;
  TH1F *hDigiTimeL1;
  TH1F *hDigiTimeL2;
  TH1F *hDigiTimeL3;
  TH1F *hDigiTimeL4;
  TH1F *hDigiTimeL5;
  TH1F *hDigiTimeL6;
  TH1F *hDigiTimeL7;
  TH1F *hDigiTimeL8;
  TH1F *hDigiCharge;
  TH2F *hDigiTimeCharge;

  TH1F *hEnergy1;
  TH1F *hEnergy2;
  TH1F *hEnergy3;
  TH1F *hEnergy4;
  TH1F *hEnergy5;
  TH1F *hEnergy6;
  TH1F *hEnergy7;
  TH1F *hEnergy8;
  TH1F *hEnergyHF;

  TH2F *hRecHitPosition;

  TH2F *hEnergyPosition;
  TH2F *hEnergyEtaPhi;
  TH2F *hEnergyEtaPhiHF;

  TH1F *hNmips;
  TH1F *hEnergyMIPS;

  TH1F *hNHitsEE;
  TH1F *hNHitsHE;

  TH1F *hJetResponse;
  TH1F *hJetResponseHF;
  TH1F *hJetResponseNose;

  TH1F *hJetResponseNoseEE;
  TH1F *hJetResponseNoseEH;

  TH2F *hJetResponseE;
  TH2F *hJetResponseEta;

  TH2F *hJetResponseHFEta;
  TH2F *hJetResponseNoseEta;

  TH1F *hSimClusterEta;
  TH1F *hRecoClusterEta;

  TH1F *hRecoHFNoseClusterEta;
  TH1F *hRecoHFNoseClusterE;
  TH1F *hRecoHFNoseClusterEee;
  TH1F *hRecoHFNoseClusterEhe;

  TH1F *hMpi0;

  //$$$$// MONEY PLOTS
  TH2F *hClusterResponseE;
  TH2F *hRecHitResponseE;
  TH1F *hClusterResponse;
  TH1F *hRecHitResponse;
  //$$$$//

  //  double coneSize=0.1;
  double coneSize=0.4; // good for jets
  double hadWeight=1.;
  bool doSingle=false;

  edm::Service<TFileService> FileService;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

void GenAnalyzer::analyzeDigi(edm::Handle<HGCalDigiCollection> digiNose, const math::XYZTLorentzVector & initialP4 , const HGCalGeometry* geom ) {

  float timeCl = -99.;
  std::vector<float> time;

  for(HGCalDigiCollection::const_iterator digi =  digiNose->begin(); digi != digiNose->end() ; ++digi ) {

    constexpr int iSample=2;
    const auto& sample = digi->sample(iSample);

    HFNoseDetId detId = digi->id();    
    int layer   = detId.layer();
    int zside   = detId.zside();

    //    if(layer!=2) continue;
    //    if(zside<0) continue;
    //    cout << " layer=" << layer << " zside=" << zside << endl;

    GlobalPoint global = geom->getPosition(detId);

    double thisDeltaR1 = ::deltaR(global.eta(), global.phi() , initialP4.eta(), initialP4.phi());
    double Dx = (global.x()- initialP4.x());
    double Dy = (global.y()- initialP4.y());
    double radius=sqrt( Dx*Dx + Dy*Dy ); 

    // 0.1 MeV                                                                                                                                                                       
    if ( thisDeltaR1 > coneSize ) continue;

    ///++++++++++++
    //  https://github.com/cms-sw/cmssw/blob/254c8cd20c45dedf06dff74efbd42773c9b89a6a/SimCalorimetry/HGCalSimProducers/python/hgcalDigitizer_cfi.py#L248
    //  https://github.com/amartelli/cmssw/blob/866397f183f7da0abecd69dcece939de52991164/RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h#L72
    ///++++++++++++

    double toaLSBToNS_ = 0.0244;
    double jitter_ = double(sample.toa()) * toaLSBToNS_;

    auto nBits=10;
    double adcSaturation_fC=100;
    float adcLSB_ = adcSaturation_fC/pow(2.,nBits); 
    double amplitude_ = double(sample.data()) * adcLSB_;

    if(sample.getToAValid()) {
      //      cout << "jitter=" << jitter_ << endl;
      //      cout << " sample=" << sample.data()  << endl;
      hDigiTime->Fill(jitter_-5);
      if(layer==1) hDigiTimeL1->Fill(jitter_-5);
      if(layer==2) hDigiTimeL2->Fill(jitter_-5);
      if(layer==3) hDigiTimeL3->Fill(jitter_-5);
      if(layer==4) hDigiTimeL4->Fill(jitter_-5);
      if(layer==5) hDigiTimeL5->Fill(jitter_-5);
      if(layer==6) hDigiTimeL6->Fill(jitter_-5);
      if(layer==7) hDigiTimeL7->Fill(jitter_-5);
      if(layer==8) hDigiTimeL8->Fill(jitter_-5);

      hDigiCharge->Fill(amplitude_);
      hDigiTimeCharge->Fill(amplitude_,jitter_);

      if(jitter_<0) continue;
      time.push_back(jitter_-5);

    }

  }

  //  if(time.size() >= 3) timeCl = hgcalsimclustertime::fixSizeHighestDensity(time);

  hParticleTime->Fill(timeCl);

}

void GenAnalyzer::analyzeClusters(std::vector<reco::CaloCluster> recoNoseClusters,const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom) {

  TLorentzVector EsumVectorNose(0.,0.,0.,0.);

    for (auto const& cl : recoNoseClusters) {
      hRecoHFNoseClusterEta->Fill(cl.eta());
      hRecoHFNoseClusterE->Fill(cl.energy());
      hRecoHFNoseClusterEee->Fill(cl.energy());
      hRecoHFNoseClusterEhe->Fill(cl.energy());

      double energy      = cl.energy();
      double   thisDeltaR1 = ::deltaR(cl.eta(), cl.phi() , initialP4.eta(), initialP4.phi());

      // 0.1 MeV
      if ( thisDeltaR1 < coneSize ) {

        double theta = 2 * atan(exp(-cl.eta()));
        double phi = cl.phi();
        double px = sin(theta)*cos(phi);
        double py = sin(theta)*sin(phi);
        double pz = cos(theta);
        TLorentzVector constituentbase(px,py,pz,1);

	EsumVectorNose+= energy*constituentbase;

      }
    }

    //    hClusterResponseE->Fill(trueE, (EsumVectorNose).Energy()/trueE);
    hClusterResponse->Fill((EsumVectorNose).Energy()/trueE);

}


double GenAnalyzer::analyzeHits(std::vector<PCaloHit> const& simHits,  const HGCRecHitCollection & recHits, const HFRecHitCollection & hfhits, const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal, bool isNose) {

  // HF geometry  http://cds.cern.ch/record/896897/files/note05_016.pdf
  // HF in the eta 2.85 - 4.2 should have 8 cells in r

  TLorentzVector EsumVectorNose(0.,0.,0.,0.);
  TLorentzVector EsumVectorHF(0.,0.,0.,0.);
  TLorentzVector EsumVectorE(0.,0.,0.,0.);
  TLorentzVector EsumVectorH(0.,0.,0.,0.);

  double EnergySumHF=0.;

  if(true) {   

    for (auto const& hfhit : hfhits) {

      uint32_t id      = hfhit.id();
      HcalDetId detId  = HcalDetId(id);
      int subdet       = detId.subdetId();
      auto depth       = detId.depth();

      hRecHitPosition->Fill(detId.ieta(),detId.iphi());

      //      cout << " detId.det() = " << detId.det() << " detId.subdetId() = " << detId.subdetId()  << " ieta = " << detId.ieta() << " iphi = " << detId.iphi()  << " depth =" << depth  << endl;

      if(detId.det() != DetId::Hcal) continue;

      if ( detId.subdet() == 4 ) {

	std::shared_ptr<const CaloCellGeometry> thisCell= geomHcal->getGeometry(detId);

	if(thisCell) {
	  const GlobalPoint & global = thisCell->getPosition();

	  
	  if(verbose) cout << " position ("  << global.x() << ", " 
	       << global.y() << ", " << global.z() << " )   ";
	  
	  double   thisDeltaR1 = ::deltaR(global.eta(), global.phi() , initialP4.eta(), initialP4.phi());
	  
	  // 0.1 MeV
	  if ( thisDeltaR1 < coneSize ) {
	    

	    EnergySumHF += hfhit.energy() ;
	    double radius=sqrt(global.x()*global.x() + global.y()*global.y());
	    
	    hEnergyPosition->Fill(abs(global.z()),radius);
	    hEnergyEtaPhiHF->Fill(global.eta(),global.phi());

	    double theta = 2 * atan(exp(-global.eta()));
	    double phi = global.phi();
	    double px = sin(theta)*cos(phi);
	    double py = sin(theta)*sin(phi);
	    double pz = cos(theta);
	    TLorentzVector constituentbase(px,py,pz,1);

	    EsumVectorHF += hfhit.energy()*constituentbase;

	  }
	
	}
      }

    }

    if(doSingle && EnergySumHF>0) hJetResponseHFEta->Fill(initialP4.eta(), EnergySumHF/initialP4.energy());
    if(doSingle && EnergySumHF>0) hJetResponseHF->Fill(EnergySumHF/trueE);

    if(verbose) cout << "   --------   end of HF recHits   " << endl;

  }

  double thicknessCorrectionEM = 0.759;
  double thicknessCorrectionHad = 0.84;
  // old v8 thicknessCorrection = 1.132
  // https://github.com/cms-sw/cmssw/blob/71bcdf3775b9bbc8911257bfc67fce073b44a0b0/RecoLocalCalo/HGCalRecProducers/python/HGCalRecHit_cfi.py#L233;

  double EnergySumNOSE = 0.;

  //  double keV2MIP = 0.044259;
  //  double fCPerMIP = 1.25;

  //      hit.energy() in KeV
  //      double numberOfMIPin1GeV = (45./1000000);
  double numberOfMIPin1GeV = (40./1000000);
  //      double numberOfMIPin1GeV = (32./1000000);

  double dedx1 = 39.500245, dedx2 = 39.756638, dedx3 = 39.756638, dedx4 = 39.756638, dedx5 = 39.756638, dedx6 = 66.020266, dedx7 = 92.283895, dedx8 = 92.283895;

  double  Esum=0.;
  double  EsumE=0.;
  double  EsumH=0.;
  //
  double  Esum1=0., Esum2=0., Esum3=0., Esum4=0., Esum5=0., Esum6=0., Esum7=0., Esum8=0.;

  int NhitsEE = 0 ;
  int NhitsHE = 0 ;

  if(true)  { 

    //    for (auto const& hit : simHits) {
    for (auto const& hit : recHits) {

      double energy      = hit.energy();
      double nMIPs = (hit.energy()/numberOfMIPin1GeV);
      //      double nMIPs = hit.energy()/numberOfMIPin1GeV;
      double time        = hit.time();
      uint32_t id        = hit.id();
      //    HepGeom::Point3D<float> gcoord;
      
      HFNoseDetId detId = HFNoseDetId(id);
      int subdet  = detId.subdetId();
      int cell    = detId.cellU();
      //    int cell2   = detId.cellV();
      int sector  = detId.waferU();
      //    int sector2 = detId.waferV();
      //    int type    = detId.type();
      int layer   = detId.layer();
      int zside   = detId.zside();
      
      GlobalPoint global = geom->getPosition(id);
      
      if(verbose) cout << " position ("  << global.x() << ", " 
		       << global.y() << ", " << global.z()  << " )";

      double   thisDeltaR1 = ::deltaR(global.eta(), global.phi() , initialP4.eta(), initialP4.phi());
      
      // 0.1 MeV
      if ( thisDeltaR1 < coneSize ) {
	
	if(false) {	
	  std::cout << "detId = " << detId ;
	  std::cout << "  subdet = " << subdet ;
	  std::cout << "  cell = " << cell ;
	  std::cout << "  sector = " << sector ;
	  std::cout << "  layer = " << layer ;
	  std::cout << "  zside = " << zside;
	  if(detId.isEE()) std::cout << " EE ";
	  if(detId.isHE()) std::cout << " HE ";
	  //       std::cout << " ETA = " << etaphi.first;
	  //       std::cout << " PHI = " << etaphi.second;
	  std::cout << "  energy = " << energy ;
	  std::cout << "  nMIPs = " << nMIPs ;
	  std::cout << "  time = " << time << std::endl;
	}
      
	double theta = 2 * atan(exp(-global.eta()));
	double phi = global.phi();
	double px = sin(theta)*cos(phi);
	double py = sin(theta)*sin(phi);
	double pz = cos(theta);
	TLorentzVector constituentbase(px,py,pz,1);

	if(detId.layer()==1) { double e=hit.energy();  Esum1 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx1)>1) NhitsEE++; };
	if(detId.layer()==2) { double e=hit.energy();  Esum2 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx2)>1) NhitsEE++; };
	if(detId.layer()==3) { double e=hit.energy();  Esum3 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx3)>1) NhitsEE++; };
	if(detId.layer()==4) { double e=hit.energy();  Esum4 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx4)>1) NhitsEE++; };
	if(detId.layer()==5) { double e=hit.energy();  Esum5 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx5)>1) NhitsEE++; };
	if(detId.layer()==6) { double e=hit.energy();  Esum6 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx6)>1) NhitsEE++; };
	if(detId.layer()==7) { double e=hit.energy();  Esum7 += e ; Esum  += hadWeight*e ; EsumH += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx7)>1) NhitsHE++; };
	if(detId.layer()==8) { double e=hit.energy();  Esum8 += e ; Esum  += hadWeight*e ; EsumH += e ; EsumVectorNose+= e*constituentbase; if((nMIPs/dedx8)>1) NhitsHE++; };

	/* below for simHits
	if(detId.layer()==1) { double e=nMIPs*dedx1*0.001/thicknessCorrectionEM;  Esum1 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==2) { double e=nMIPs*dedx2*0.001/thicknessCorrectionEM;  Esum2 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==3) { double e=nMIPs*dedx3*0.001/thicknessCorrectionEM;  Esum3 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==4) { double e=nMIPs*dedx4*0.001/thicknessCorrectionEM;  Esum4 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==5) { double e=nMIPs*dedx5*0.001/thicknessCorrectionEM;  Esum5 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==6) { double e=nMIPs*dedx6*0.001/thicknessCorrectionEM;  Esum6 += e ; Esum  += e ; EsumE += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==7) { double e=nMIPs*dedx7*0.001/thicknessCorrectionEM;  Esum7 += e ; Esum  += hadWeight*e ; EsumH += e ; EsumVectorNose+= e*constituentbase; };
	if(detId.layer()==8) { double e=nMIPs*dedx8*0.001/thicknessCorrectionEM;  Esum8 += e ; Esum  += hadWeight*e ; EsumH += e ; EsumVectorNose+= e*constituentbase; };
	*/
	double radius=sqrt(global.x()*global.x() + global.y()*global.y());
	
	hEnergyPosition->Fill(abs(global.z()),radius);
	hEnergyEtaPhi->Fill(global.eta(),global.phi());
	hNmips->Fill(nMIPs);
	hEnergyMIPS->Fill(energy);
	
      }

    }

  }

  hNHitsEE->Fill(NhitsEE);
  hNHitsHE->Fill(NhitsHE);


  EnergySumNOSE = Esum;

  if(Esum>0) {

    hEnergyE->Fill(EsumE);
    hEnergyH->Fill(EsumH);
    hEnergy->Fill(Esum);
    hEnergyFrac->Fill(EsumE/Esum);
    //
    hEnergy1->Fill(Esum1);
    hEnergy2->Fill(Esum2);
    hEnergy3->Fill(Esum3);
    hEnergy4->Fill(Esum4);
    hEnergy5->Fill(Esum5);
    hEnergy6->Fill(Esum6);
    hEnergy7->Fill(Esum7);
    hEnergy8->Fill(Esum8);

  }

  if(Esum>0) hJetResponseNose->Fill((EsumE+EsumH)/trueE);
  if(EsumE>0) hJetResponseNoseEE->Fill(EsumE/trueE);
  if(EsumH>0) hJetResponseNoseEH->Fill(EsumH/trueE);
  if(EnergySumNOSE>0) hJetResponseNoseEta->Fill(initialP4.eta(), EnergySumNOSE/trueE);

  // SCALAR SUM
  //  hJetResponseNoseE->Fill(trueE, EnergySumNOSE/trueE);
  //  if(doSingle) hJetResponseEta->Fill(initialP4.eta(),(EnergySumNOSE+hadWeight*EnergySumHF)/trueE);
  //  if(doSingle)  hJetResponse->Fill((EnergySumNOSE+hadWeight*EnergySumHF)/trueE);


  // VECTOR SUM
  hRecHitResponse->Fill((EsumVectorNose).Energy()/trueE);
  //  hRecHitResponseE->Fill(trueE, (EsumVectorNose).Energy()/trueE);
  hJetResponseE->Fill(trueE, (EsumVectorNose+EsumVectorHF).Energy()/trueE);
  hJetResponseEta->Fill(initialP4.eta(), (EsumVectorNose+EsumVectorHF).Energy()/trueE);
  hJetResponse->Fill((EsumVectorNose+EsumVectorHF).Energy()/trueE);

  if (isNose) { return EnergySumNOSE;
  } else if (!isNose) { return EnergySumHF;
  } else { return 0.; }

}

void GenAnalyzer::getSingle(edm::Handle<reco::GenParticleCollection> genParticles, std::vector<CaloParticle> caloParticles, std::vector<reco::CaloCluster> recoNoseClusters, edm::Handle<HGCalDigiCollection> digiNose, const HGCRecHitCollection& recHitNose, const std::vector<PCaloHit>& hits, const HFRecHitCollection& hfhits, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {
  
  doSingle=true;
  //  int pdgId=22;
  //  int pdgId=211;

  if(false) {
    for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

      // 22 photon, 130 k0L , 211 pi
      //   bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==130 );
      //    bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==22 );
      bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==211 );
      if (!isPhoton) continue;

      /*
	bool isGluon = ( genpart->status()==23 and abs(genpart->pdgId())==21 and abs(genpart->eta())<4.2 and abs(genpart->eta())>3  );
	if (!isGluon) continue;
      */
      /*
	bool isQuark = ( genpart->status()==23 and (abs(genpart->pdgId())==1 or abs(genpart->pdgId())==2 or abs(genpart->pdgId())==3) and abs(genpart->eta())<4.2 and abs(genpart->eta())>3  );
	if (!isQuark) continue;
      */
      //   cout << ">>>>>>> pid,status,px,py,pz,e,eta,phi= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << " , ( " << genpart->eta() << " , " << genpart->phi() << ") " << endl;

      hGenParticleEta->Fill(abs(genpart->eta()));
      hGenParticlePt->Fill(genpart->pt());
      hGenParticleE->Fill(genpart->energy());

      // this is DIGI
      double E = analyzeHits(hits, recHitNose, hfhits, genpart->p4(), genpart->energy(), geom, geomHcal, true);
      //   analyzeDigi(digiNose, genpart->p4(), geom);
      analyzeClusters(recoNoseClusters, genpart->p4(), genpart->energy(), geom);

    }
  }

  if(true) {
    for (auto const& cl : caloParticles) {
      //     cout << " cl.pdgId() = " << cl.pdgId() << " cl.energy() = " << cl.energy() << " cl.eta() = " << cl.eta() << endl;
      bool isPhoton = ( abs(cl.pdgId())==22 );
      //      bool isPhoton = ( abs(cl.pdgId())==211 );
      if (!isPhoton) continue;

      const math::XYZTLorentzVector p4=(math::XYZTLorentzVector) cl.p4();

      // this is DIGI
      double E = analyzeHits(hits, recHitNose, hfhits, p4, cl.energy(), geom, geomHcal, true);
      //   analyzeDigi(digiNose, genpart->p4(), geom);
      analyzeClusters(recoNoseClusters, p4, cl.energy(), geom);

    }
  }
}


std::vector<const reco::GenJet*> GenAnalyzer::doVBSselection(edm::Handle<reco::GenParticleCollection> genParticles, edm::Handle<reco::GenJetCollection> genJets) {

  std::vector<const reco::GenJet*> selectedGenJets;
  const reco::GenParticle * genP1 = NULL;
  const reco::GenParticle * genP2 = NULL;
  
  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){
    
    bool isEmu = ( genpart->isPromptFinalState() and genpart->fromHardProcessFinalState() and (abs(genpart->pdgId())==11 or abs(genpart->pdgId())==13 ) and genpart->pt()>10 );
    bool isNu = ( abs(genpart->pdgId())==12 or abs(genpart->pdgId())==14 or abs(genpart->pdgId())==16);
    bool isTau = ( abs(genpart->pdgId()==15) );
    bool isWZ  = ( abs(genpart->pdgId()==24) or abs(genpart->pdgId()==23 ) );
    
    if(verbose and (isEmu or isTau or isWZ)) cout << ">>>>>>> pid,status,px,py,pz,e,eta,phi= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << " , ( " << genpart->eta() << " , " << genpart->phi() << ") " << endl;
    
    if(isEmu and genP1==NULL) genP1 = &*genpart;
    if(isEmu and genP2==NULL and &*genpart!=genP1) genP2 = &*genpart;
    
    if(isEmu) hLepEta->Fill(genpart->eta());
    if(isEmu) hLepPt->Fill(genpart->pt());
    
  }
  
  if(verbose && genP1) cout << ">>> genP1 >>>> pid,status,px,py,pz,e,eta,phi= "  << genP1->pdgId() << " , " << genP1->status() << " , " << genP1->px() << " , " << genP1->py() << " , " << genP1->pz() << " , " << genP1->energy() << " , ( " << genP1->eta() << " , " << genP1->phi() << ") " << endl;
  
  if((genP1 == NULL) or (genP2 == NULL)) return selectedGenJets;

  hLepCharge->Fill(genP1->charge()*genP2->charge());  
  
  if(verbose) cout << "-----------------------------------------------------------------------------------------" << endl;
  

  for(reco::GenJetCollection::const_iterator genJ = genJets->begin(); genJ != genJets->end(); ++genJ){
    
    if( verbose && genJ->energy()>100) cout << " JETS >>>>>>> status,  px,py,pz,e,eta= "  << " , " << genJ->status() << " , " << genJ->px() << " , " << genJ->py() << " , " << genJ->pz() << " , " << genJ->energy() << " , (" << genJ->eta() << " , " << genJ->phi() << ") " << endl;
    if( verbose && abs(genJ->eta())>3 ) cout << " this is candidate " << endl; 
    
    double thisDeltaR1 = 0.;
    double thisDeltaR2 = 0.;
    
    //     if(genJ->DeltaR(genP1)<0.4 or genJ->DeltaR(genP2)<0.4)  cout << " this is lepton " << endl;
    if(genP1)  {
      thisDeltaR1 = ::deltaR(genP1->eta(), genP1->phi(), genJ->eta(), genJ->phi());
      //       if( thisDeltaR1<0.4 ) cout << "this is a duplicate " << endl;
    }     
    
    if(genP2)  {
      thisDeltaR2 = ::deltaR(genP2->eta(), genP2->phi(), genJ->eta(), genJ->phi());
      //       if( thisDeltaR2<0.4 )  cout << "this is a duplicate " << endl;
    }
    
    //    if( genJ->pt()>30 and abs(genJ->eta())<4.7 
    if( genJ->pt()>15 and abs(genJ->eta())<4.7 
	and (genJ->eta() < std::min(genP1->eta(),genP2->eta()) or genJ->eta() > std::max(genP1->eta(),genP2->eta())) 
	and thisDeltaR1 > 0.4 and thisDeltaR2 > 0.4 ) {
      selectedGenJets.push_back( &(*genJ) ); 
    }
  }
  
   /*
   // require the leptonEta between the two jets
   if( selectedGenJets.at(indexRapMax)->rapidity() < std::max(genP1->eta(),genP2->eta()) ) return;
   if( selectedGenJets.at(indexRapMin)->rapidity() > std::min(genP1->eta(),genP2->eta()) ) return;
   */

  return selectedGenJets;

}

//
// constructors and destructor
//
GenAnalyzer::GenAnalyzer(const edm::ParameterSet& iConfig)
// :
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{

  usesResource("TFileService");  
  genParticlesTag_ = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("GenParticleTag", edm::InputTag("genParticles")));
  genJetsTag_ = consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("GenJetTag", edm::InputTag("ak4GenJets")));

  simHitTag_ = mayConsume<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HFNoseHits"));
  //  simHcalTag_ = mayConsume<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  simHcalTag_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  tok_hf_ = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<string>("HFRecHits","hfreco"));

  //  digiNose_ = consumes<HGCalDigiCollection>(iConfig.getParameter<edm::InputTag>("simHFNoseUnsuppressedDigis","HFNose"));
  digiNose_ = consumes<HGCalDigiCollection>(iConfig.getUntrackedParameter<edm::InputTag>("DIGITAG",edm::InputTag("simHFNoseUnsuppressedDigis:HFNose")));

  recHitNose_ = consumes<HGCRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("REChitTAG",edm::InputTag("HGCalRecHit:HGCHFNoseRecHits")));

  simClusterTag_ = mayConsume<std::vector<SimCluster>>(edm::InputTag("mix","MergedCaloTruth"));
  caloParticleTag_ = mayConsume<std::vector<CaloParticle>>(edm::InputTag("mix","MergedCaloTruth"));

  recoClusterTag_ = consumes<std::vector<reco::CaloCluster>>(edm::InputTag("hgcalLayerClusters"));
  recoNoseClusterTag_ = consumes<std::vector<reco::CaloCluster>>(edm::InputTag("hgcalLayerClustersHFNose"));

  /*
    cout << "---------------------------------------------------" << endl; 
    cout << "-------------------  Token done  ------------------" << endl;
    cout << "-----------------  book-histo done  ---------------" << endl;
    cout << "---------------------------------------------------" << endl; 
  */   
  
  hSimClusterEta = FileService->make<TH1F>("hSimClusterEta","hSimClusterEta", 100, -5. , 5.);
  hMpi0 = FileService->make<TH1F>("hMpi0","hMpi0", 100, 0. , 1.);
  hDRpi0 = FileService->make<TH1F>("hDRpi0","hDRpi0", 100, 0. , 0.5);

  hRecoClusterEta = FileService->make<TH1F>("hRecoClusterEta","hRecoClusterEta", 100, -5. , 5.);
  hRecoHFNoseClusterEta = FileService->make<TH1F>("hRecoHFNoseClusterEta","hRecoHFNoseClusterEta", 100, -5. , 5.);

  hRecoHFNoseClusterE = FileService->make<TH1F>("hRecoHFNoseClusterE","hRecoHFNoseClusterE", 100, 0. , 1000.);
  hRecoHFNoseClusterEee = FileService->make<TH1F>("hRecoHFNoseClusterEee","hRecoHFNoseClusterEee", 100, 0. , 1000.);
  hRecoHFNoseClusterEhe = FileService->make<TH1F>("hRecoHFNoseClusterEhe","hRecoHFNoseClusterEhe", 100, 0. , 1000.);

  hLepPt = FileService->make<TH1F>("hLepPt","hLepPt", 100, 0. , 100.);
  hLepEta = FileService->make<TH1F>("hLepEta","hLepEta", 100, -5. , 5.);
  hLepCharge = FileService->make<TH1F>("hLepCharge","hLepCharge", 4, -2. , 2.);

  hGenParticleEta = FileService->make<TH1F>("hGenParticleEta","hGenParticleEta", 100, 0. , 5.);
  hGenParticlePt = FileService->make<TH1F>("hGenParticlePt","hGenParticlePt", 100, 0. , 100.);
  hGenParticleE = FileService->make<TH1F>("hGenParticleE","hGenParticleE", 100, 0. , 1000.);

  hJetEta = FileService->make<TH1F>("hJetEta","hJetEta", 100, 0. , 5.);
  hJetEtaLarge = FileService->make<TH1F>("hJetEtaLarge","hJetEtaLarge", 100, 0. , 5.);
  hJetDeltaRap = FileService->make<TH1F>("hJetDeltaRap","hJetDeltaRap", 100, 0. , 20.);

  hJetPtMin = FileService->make<TH1F>("hJetPtMin","hJetPtMin", 100, 0. , 200.);
  hJetPtMax = FileService->make<TH1F>("hJetPtMax","hJetPtMax", 100, 0. , 200.);

  hJetDiMass = FileService->make<TH1F>("hJetDiMass","hJetDiMass", 100, 0. , 1000.);

  hParticleTime = FileService->make<TH1F>("hParticleTime","hParticleTime", 100, -1. , 1.);
  hDigiTime = FileService->make<TH1F>("hDigiTime","hDigiTime", 500, -10. , 10.);
  hDigiTimeL1 = FileService->make<TH1F>("hDigiTimeL1","hDigiTimeL1", 500, -10. , 10.);
  hDigiTimeL2 = FileService->make<TH1F>("hDigiTimeL2","hDigiTimeL2", 500, -10. , 10.);
  hDigiTimeL3 = FileService->make<TH1F>("hDigiTimeL3","hDigiTimeL3", 500, -10. , 10.);
  hDigiTimeL4 = FileService->make<TH1F>("hDigiTimeL4","hDigiTimeL4", 500, -10. , 10.);
  hDigiTimeL5 = FileService->make<TH1F>("hDigiTimeL5","hDigiTimeL5", 500, -10. , 10.);
  hDigiTimeL6 = FileService->make<TH1F>("hDigiTimeL6","hDigiTimeL6", 500, -10. , 10.);
  hDigiTimeL7 = FileService->make<TH1F>("hDigiTimeL7","hDigiTimeL7", 500, -10. , 10.);
  hDigiTimeL8 = FileService->make<TH1F>("hDigiTimeL8","hDigiTimeL8", 500, -10. , 10.);
  hDigiCharge = FileService->make<TH1F>("hDigiCharge","hDigiCharge", 100, 0. , 100.);
  hDigiTimeCharge = FileService->make<TH2F>("hDigiTimeCharge","hDigiTimeCharge", 100, 0. , 100., 500, -10. , 10.);

  hEnergyE = FileService->make<TH1F>("hEnergyE","hEnergyE", 100, 0. , 1.);
  hEnergyH = FileService->make<TH1F>("hEnergyH","hEnergyH", 100, 0. , 1.);
  hEnergy = FileService->make<TH1F>("hEnergy","hEnergy", 100, 0. , 1.);
  hEnergyFrac = FileService->make<TH1F>("hEnergyFrac","hEnergyFrac", 100, 0. , 1.);

  hRecHitPosition = FileService->make<TH2F>("hRecHitPosition","hRecHitPosition", 82, -41., 41., 73 , 0., 73.);

  hEnergyPosition = FileService->make<TH2F>("hEnergyPosition","hEnergyPosition", 600, 1040., 1200., 500, 0., 150.);
  hEnergyEtaPhi = FileService->make<TH2F>("hEnergyEtaPhi","hEnergyEtaPhi", 100, -5., 5., 18*2, -3.14, 3.14);
  // modularity of the HF is 18 tile in phi
  hEnergyEtaPhiHF = FileService->make<TH2F>("hEnergyEtaPhiHF","hEnergyEtaPhiHF", 100, -5., 5., 628, -3.14, 3.14);
  hNmips = FileService->make<TH1F>("hNmips","hNmips", 500, 0. , 500.);
  hEnergyMIPS = FileService->make<TH1F>("hEnergyMIPS","hEnergyMIPS", 100, 0. , 0.0001);
  hNHitsEE = FileService->make<TH1F>("hNHitsEE","hNHitsEE", 500, 0. , 1000.);
  hNHitsHE = FileService->make<TH1F>("hNHitsHE","hNHitsHE", 500, 0. , 1000.);

  double myInt=10.;

  hEnergy1 = FileService->make<TH1F>("E Layer1","E Layer1", 100, 0. , 10*myInt);
  hEnergy2 = FileService->make<TH1F>("E Layer2","E Layer2", 100, 0. , 10*myInt);
  hEnergy3 = FileService->make<TH1F>("E Layer3","E Layer3", 100, 0. , 10*myInt);
  hEnergy4 = FileService->make<TH1F>("E Layer4","E Layer4", 100, 0. , 10*myInt);
  hEnergy5 = FileService->make<TH1F>("E Layer5","E Layer5", 100, 0. , 10*myInt);
  hEnergy6 = FileService->make<TH1F>("E Layer6","E Layer6", 100, 0. , 10*myInt);
  hEnergy7 = FileService->make<TH1F>("E Layer7","E Layer7", 100, 0. , 10*myInt);
  hEnergy8 = FileService->make<TH1F>("E Layer8","E Layer8", 100, 0. , 10*myInt);
  hEnergyHF = FileService->make<TH1F>("E LayerHF","E LayerHF", 100, 0. , 10*myInt);

  hJetResponse = FileService->make<TH1F>("hJetResponse","hJetResponse", 100, 0. , 2.);
  hJetResponseHF = FileService->make<TH1F>("hJetResponseHF","hJetResponseHF", 100, 0. , 2.);
  hJetResponseNose = FileService->make<TH1F>("hJetResponseNose","hJetResponseNose", 100, 0. , 2.);
  hJetResponseNoseEE = FileService->make<TH1F>("hJetResponseNoseEE","hJetResponseNoseEE", 100, 0. , 2.);
  hJetResponseNoseEH = FileService->make<TH1F>("hJetResponseNoseEH","hJetResponseNoseEH", 100, 0. , 2.);

  hJetResponseEta = FileService->make<TH2F>("hJetResponseEta","hJetResponseEta", 100., -5., 5., 100, 0. , 2.);
  hJetResponseHFEta = FileService->make<TH2F>("hJetResponseHFEta","hJetResponseHFEta", 100., -5., 5., 100, 0. , 2.);
  hJetResponseNoseEta = FileService->make<TH2F>("hJetResponseNoseEta","hJetResponseNoseEta", 100., -5., 5., 100, 0. , 2.);

  hJetResponseE = FileService->make<TH2F>("hJetResponseE","hJetResponseE", 1000., 0, 1000., 100, 0. , 2.);

  hClusterResponseE = FileService->make<TH2F>("hClusterResponseE","hClusterResponseE", 1000., 0, 1000., 100, 0. , 2.);
  hRecHitResponseE = FileService->make<TH2F>("hRecHitResponseE","hRecHitResponseE", 1000., 0, 1000., 100, 0. , 2.);

  hClusterResponse = FileService->make<TH1F>("hClusterResponse","hClusterResponse", 100, 0. , 2.);
  hRecHitResponse = FileService->make<TH1F>("hRecHitResponse","hRecHitResponse", 100, 0. , 2.);

}


GenAnalyzer::~GenAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  bool doVBS=false;
  bool runSingle=true;
  
  using namespace edm;
  
  edm::ESHandle<HGCalGeometry> geomHGcalHandle;
  iSetup.get<IdealGeometryRecord>().get("HGCalHFNoseSensitive",geomHGcalHandle);
  const HGCalGeometry* geom = (geomHGcalHandle.product());
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry *geoHcal = geoHandle->getSubdetectorGeometry(DetId::Hcal, HcalForward);

  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByToken(genParticlesTag_, genParticles);
  
  edm::Handle<HGCalDigiCollection> digiNose;
  iEvent.getByToken(digiNose_, digiNose);

  edm::Handle<HGCRecHitCollection> recHitNose;
  iEvent.getByToken(recHitNose_, recHitNose);
  const HGCRecHitCollection noseRecHits = *(recHitNose.product());

  edm::Handle<reco::GenJetCollection> genJets;
  iEvent.getByToken(genJetsTag_, genJets);

  edm::Handle<edm::PCaloHitContainer> theCaloHitContainers;
  iEvent.getByToken(simHitTag_, theCaloHitContainers);
  std::vector<PCaloHit>               caloHits;
  
  ///###

  edm::Handle<std::vector<SimCluster>> simClustersH;
  iEvent.getByToken(simClusterTag_, simClustersH);
  std::vector<SimCluster>  simClusters = *(simClustersH.product());;

  edm::Handle<std::vector<CaloParticle>> caloParticlesH;
  iEvent.getByToken(caloParticleTag_, caloParticlesH);
  std::vector<CaloParticle>  caloParticles = *(caloParticlesH.product());;

  ///###

  edm::Handle<std::vector<reco::CaloCluster>> recoClustersH;
  iEvent.getByToken(recoClusterTag_, recoClustersH);
  std::vector<reco::CaloCluster>  recoClusters = *(recoClustersH.product());;

  edm::Handle<std::vector<reco::CaloCluster>> recoNoseClustersH;
  iEvent.getByToken(recoNoseClusterTag_, recoNoseClustersH);
  std::vector<reco::CaloCluster>  recoNoseClusters = *(recoNoseClustersH.product());;

  if(theCaloHitContainers.isValid()) {
   caloHits.insert(caloHits.end(), theCaloHitContainers->begin(), 
		   theCaloHitContainers->end());
  }   

   edm::Handle<HFRecHitCollection> hf_hits;
   iEvent.getByToken(tok_hf_,hf_hits);
   const HFRecHitCollection Hithf = *(hf_hits.product());

   /*
   cout << "------------------------------------------------" << endl; 
   cout << "-----------------  Handles done ----------------" << endl;
   cout << "------------------------------------------------" << endl; 
   */   

   //   std::vector<CaloParticle> clPos;
   //   std::vector<CaloParticle> clNeg;
   //   for (auto const& cl : caloParticles) {

   std::vector<SimCluster> clPos;
   std::vector<SimCluster> clNeg;

   for (auto const& cl : simClusters) {

     hSimClusterEta->Fill(cl.eta());
     //     cout << " cl.pdgId() = " << cl.pdgId() << " cl.energy() = " << cl.energy() << " cl.eta() = " << cl.eta() << endl;

     if(cl.eta()>0) clPos.push_back(cl);
     if(cl.eta()<0) clNeg.push_back(cl);

   }

   /*
   double mpi0Pos = (clPos.at(0).p4() + clPos.at(1).p4()).M();
   double mpi0Neg = (clNeg.at(0).p4() + clNeg.at(1).p4()).M();

   double   thisDeltaPos = ::deltaR(clPos.at(0).eta(), clPos.at(0).phi() , clPos.at(1).eta(), clPos.at(1).phi());
   double   thisDeltaNeg = ::deltaR(clNeg.at(0).eta(), clNeg.at(0).phi() , clNeg.at(1).eta(), clNeg.at(1).phi());

   hMpi0->Fill(mpi0Pos);
   hMpi0->Fill(mpi0Neg);

   hDRpi0->Fill(thisDeltaPos);
   hDRpi0->Fill(thisDeltaNeg);

   */

   for (auto const& cl : recoClusters) {
     hRecoClusterEta->Fill(cl.eta());
   }

   if(runSingle) getSingle(genParticles, caloParticles, recoNoseClusters, digiNose, noseRecHits, caloHits, Hithf, geom, geoHcal );

   return;
   
   if(not doVBS) return;
   
   std::vector<const reco::GenJet*> selectedGenJets = doVBSselection(genParticles, genJets);
   
   if(doVBS and selectedGenJets.size()<2) return;

  //////////
  
  int indexRapMin=0;
  int indexRapMax=0;
   
   double rapMax = selectedGenJets.at(indexRapMax)->rapidity();
   double rapMin = selectedGenJets.at(indexRapMin)->rapidity();

   //   int njets = 0;
   for (unsigned iGenJet = 0; iGenJet < selectedGenJets.size(); ++iGenJet) {
     const reco::GenJet* genJet = selectedGenJets.at(iGenJet);

     if( genJet->rapidity() > rapMax ) { rapMax= genJet->rapidity(); indexRapMax = iGenJet; }
     if( genJet->rapidity() < rapMin ) { rapMin= genJet->rapidity(); indexRapMin = iGenJet; }

     //     double genJetE   = genJet->energy();
     double genJetPt  = genJet->pt();
     //     double genJetEta = genJet->eta();
     //     double genJetPhi = genJet->phi();

     if (verbose) cout << "pt(j)=" << genJetPt << endl;

   }

   //////////

   double maxAbsJetEta=std::max(
				abs(selectedGenJets.at(indexRapMin)->rapidity()), 
				abs(selectedGenJets.at(indexRapMax)->rapidity())
				);

   hJetEta->Fill(maxAbsJetEta);

   double rapGap = selectedGenJets.at(indexRapMax)->rapidity()-selectedGenJets.at(indexRapMin)->rapidity();
   if(rapGap>4) hJetEtaLarge->Fill( maxAbsJetEta ) ;
   hJetDeltaRap->Fill(rapGap);

   double invMass = (selectedGenJets.at(indexRapMax)->p4()+selectedGenJets.at(indexRapMin)->p4()).M();
   hJetDiMass->Fill(invMass);

   hJetPtMin->Fill(selectedGenJets.at(indexRapMin)->pt() );
   hJetPtMax->Fill(selectedGenJets.at(indexRapMax)->pt() );


   // need to add the VBS selection 
   // mjj> 500
   // Delta rap > 2.5

   //   cout << "----------------------------------------------  starting HFnose caloHits !! HGC !! -------------------------------------------" << endl;

   bool isNose=true;

   double EMaxNose= analyzeHits(caloHits, noseRecHits, Hithf, selectedGenJets.at(indexRapMax)->p4(), selectedGenJets.at(indexRapMax)->energy(), geom, geoHcal, isNose);
   double EMinNose= analyzeHits(caloHits, noseRecHits, Hithf, selectedGenJets.at(indexRapMin)->p4(), selectedGenJets.at(indexRapMin)->energy(), geom, geoHcal, isNose);


   //   cout << "----------------------------------------------  starting HF caloHits -------------------------------------------" << endl;

   isNose=false;

   double EMax= analyzeHits(caloHits, noseRecHits, Hithf, selectedGenJets.at(indexRapMax)->p4(), selectedGenJets.at(indexRapMax)->energy(), geom, geoHcal, isNose);
   if(EMax>0 and fabs(selectedGenJets.at(indexRapMax)->eta())>3.2 ) hJetResponseHF->Fill(EMax/selectedGenJets.at(indexRapMax)->energy());
   if(EMax>0) hJetResponseHFEta->Fill(selectedGenJets.at(indexRapMax)->eta(), EMax/selectedGenJets.at(indexRapMax)->energy());
   
   double EMin= analyzeHits(caloHits, noseRecHits, Hithf, selectedGenJets.at(indexRapMin)->p4(), selectedGenJets.at(indexRapMin)->energy(), geom, geoHcal, isNose);
   if(EMin>0 and fabs(selectedGenJets.at(indexRapMin)->eta())>3.2 ) hJetResponseHF->Fill(EMin/selectedGenJets.at(indexRapMin)->energy());
   if(EMin>0) hJetResponseHFEta->Fill(selectedGenJets.at(indexRapMin)->eta(),EMin/selectedGenJets.at(indexRapMin)->energy());

   if((EMinNose+EMin)>0) {
     hJetResponseEta->Fill(selectedGenJets.at(indexRapMin)->eta(),(EMinNose+EMin)/selectedGenJets.at(indexRapMin)->energy());
     if(fabs(selectedGenJets.at(indexRapMin)->eta())>3.2) hJetResponse->Fill((EMinNose+EMin)/selectedGenJets.at(indexRapMin)->energy());
     hEnergyHF->Fill(EMin);	    

   }
   if((EMaxNose+EMax)>0) {
     hJetResponseEta->Fill(selectedGenJets.at(indexRapMax)->eta(),(EMaxNose+EMax)/selectedGenJets.at(indexRapMax)->energy());
     if(fabs(selectedGenJets.at(indexRapMax)->eta())>3.2) hJetResponse->Fill((EMaxNose+EMax)/selectedGenJets.at(indexRapMax)->energy());
     hEnergyHF->Fill(EMax);	    
   }

}


// ------------ method called once each job just before starting event loop  ------------
void
GenAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
GenAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}


//define this as a plug-in
DEFINE_FWK_MODULE(GenAnalyzer);
