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
//#include "FWCore/Framework/interface/ConsumesCollector.h"

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
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

// L1
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"
#include "DataFormats/L1THGCal/interface/HGCalTower.h"

//
#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

//
// class declaration
//

using namespace std;
using namespace edm;
using namespace reco;
using namespace l1t;


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


//using reco::TrackCollection;

class L1Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit L1Analyzer(const edm::ParameterSet&);
  ~L1Analyzer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  double analyze2Dcl(const l1t::HGCalClusterBxCollection& L12dNose,
		     const math::XYZTLorentzVector & initialP4 , const double & jetE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  double analyze3Dcl(const l1t::HGCalMulticlusterBxCollection& L13dNose,
		     const math::XYZTLorentzVector & initialP4 , const double & jetE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  double analyzeTowers(const l1t::HGCalTowerBxCollection & towers,
		       const math::XYZTLorentzVector & initialP4 , const double & jetE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  void getSingle(edm::Handle<reco::GenParticleCollection> genParticles, /*std::vector<CaloParticle> caloParticles,*/ const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const l1t::HGCalTowerBxCollection & towers, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  // ----------member data ---------------------------
  //  void ClearVariables();

  TH1F *hGenParticleEta;
  TH1F *hGenParticlePt;
  TH1F *hGenParticleE;

  bool verbose = false;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTag_;
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticleTag_;

  edm::EDGetToken clusters_token_, multiclusters_token_, towers_token_;
  HGCalTriggerTools triggerTools_;
  int pdgId_=22; // default do photons

  TH1F *hResponse;
  TH1F *hResponse2D;
  TH1F *hResponse3D;
  TH1F *hResponseTow;
  TH2F *hResponseE;
  TH2F *hResponseEta;

  TH2F *hClusterEtaPhi;
  TH1F *hClusterEnergy;
  TH1F *hClusterLayer;
  TH1F *hClusterSize;

  TH2F *h3DClusterEtaPhi;
  TH1F *h3DClusterEnergy;
  TH1F *h3DClusterLayer;
  TH1F *h3DClusterSize;

  TH1F *hTowClusterEnergy;

  TH2F *hClusterLayerRho;
  TH2F *hRecHitPositionNull;

  TH1F *hMax2DClResponse;
  TH1F *hMax3DClResponse;

  TH1F *hEemTow;
  TH1F *hEhadTow;

  TH1F * hClSigmaEtaEtaTot;
  TH1F * hClSigmaEtaEtaMax;
  TH1F * hClSigmaPhiPhiTot;
  TH1F * hClSigmaPhiPhiMax;
  TH1F * hClSigmaRRTot;
  TH1F * hClSigmaRRMax;
  TH1F * hClZZ;
  TH1F * hHoverE;

  double coneSize=0.1;
  //  double coneSize=0.4; // good for jets
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

double L1Analyzer::analyze2Dcl(const l1t::HGCalClusterBxCollection& L12dNose,
			     const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  TLorentzVector EsumVectorL1Nose(0.,0.,0.,0.);
  TLorentzVector EsumVectorL1HF(0.,0.,0.,0.);
  double maxECLuster = 0;

  for (auto const& cl : L12dNose) {

    double   thisDeltaR1 = ::deltaR(cl.eta(), cl.phi() , initialP4.eta(), initialP4.phi());

    if ( thisDeltaR1 < coneSize ) {

      if(cl.energy()>maxECLuster) maxECLuster=cl.energy();

      hClusterEnergy->Fill(cl.energy());
      hClusterEtaPhi->Fill(cl.eta(),cl.phi());

      double theta = 2 * atan(exp(cl.eta()));
      double phi = cl.phi();
      double px = sin(theta)*cos(phi);
      double py = sin(theta)*sin(phi);
      double pz = cos(theta);
      TLorentzVector constituentbase(px,py,pz,1);
      EsumVectorL1Nose += cl.energy()*constituentbase;

    }
  }

  hMax2DClResponse->Fill(maxECLuster/trueE);
  hResponse2D->Fill((EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);

  //  if(maxECLuster/trueE<0.01) hRecHitPositionNull->Fill(initialP4.eta(),initialP4.phi());

  return EsumVectorL1Nose.Energy();

}


double L1Analyzer::analyze3Dcl(const l1t::HGCalMulticlusterBxCollection& L13dNose,
			     const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  TLorentzVector EsumVectorL1Nose(0.,0.,0.,0.);
  TLorentzVector EsumVectorL1HF(0.,0.,0.,0.);
  
  double maxECLuster = 0;

  for (auto const& cl : L13dNose) {

    double   thisDeltaR1 = ::deltaR(cl.eta(), cl.phi() , initialP4.eta(), initialP4.phi());

    //    cout << " thisDeltaR = " <<  thisDeltaR1 <<    "( " << cl.eta() << " " << cl.phi() << "   ) "<< endl;

    if ( thisDeltaR1 < coneSize ) {

      //      cout << " layerWithOffset " << triggerTools_.layerWithOffset(cl.detId()) << endl;

      if(cl.energy()>maxECLuster) maxECLuster=cl.energy();

      h3DClusterEnergy->Fill(cl.energy());
      h3DClusterEtaPhi->Fill(cl.eta(),cl.phi());

      if(cl.energy() > 5 ) {
	//      unsigned layer = 9999.;
      unsigned layer = triggerTools_.layerWithOffset(cl.detId());
      h3DClusterLayer->Fill(layer);
      h3DClusterSize->Fill(cl.constituents().size());

      //	cout << " constituents " << cl.constituents().size() << " layer = " << layer << " eta = " << cl.eta() << " phi = " << cl.phi() << " cl.energy() = " << cl.energy() << " trueEnergy = " << initialP4.energy() << endl;
      }

      hClSigmaEtaEtaTot->Fill(cl.sigmaEtaEtaTot());
      hClSigmaEtaEtaMax->Fill(cl.sigmaEtaEtaMax());
      hClSigmaPhiPhiTot->Fill(cl.sigmaPhiPhiTot());
      hClSigmaPhiPhiMax->Fill(cl.sigmaEtaEtaMax());
      hClSigmaRRTot->Fill(cl.sigmaRRTot());
      hClSigmaRRMax->Fill(cl.sigmaRRMax());
      hClZZ->Fill(cl.sigmaZZ());
      hHoverE->Fill(cl.hOverE());

      //      GlobalPoint global= triggerTools_.getTCPosition(DetId(cl.detId()));
      //      cout << " global.x()=" << global.x() << " global.y()=" << global.y()  << endl;
      //      double rho = sqrt(global.x()*global.x()+global.y()*global.y());
      //      hClusterLayerRho->Fill(layer,rho);
      
      double theta = 2 * atan(exp(cl.eta()));
      double phi = cl.phi();
      double px = sin(theta)*cos(phi);
      double py = sin(theta)*sin(phi);
      double pz = cos(theta);
      TLorentzVector constituentbase(px,py,pz,1);
      EsumVectorL1Nose += cl.energy()*constituentbase;
      
    }

  }
  
  /*
    hResponseE->Fill(trueE, (EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);
    hResponseEta->Fill(initialP4.eta(), (EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);
  */


  hMax3DClResponse->Fill(maxECLuster/trueE);
  hResponse3D->Fill((EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);

  return EsumVectorL1Nose.Energy();

}

double L1Analyzer::analyzeTowers(const l1t::HGCalTowerBxCollection & towers,
				 const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  TLorentzVector EsumVectorL1Nose(0.,0.,0.,0.);
  TLorentzVector EsumVectorL1HF(0.,0.,0.,0.);

  TLorentzVector EsumVectorL1NoseEM(0.,0.,0.,0.);
  TLorentzVector EsumVectorL1NoseHAD(0.,0.,0.,0.);

  for (auto const& tow : towers) {

    double   thisDeltaR1 = ::deltaR(tow.eta(), tow.phi() , initialP4.eta(), initialP4.phi());

    if ( thisDeltaR1 < coneSize ) {

      hTowClusterEnergy->Fill(tow.energy());

      double theta = 2 * atan(exp(tow.eta()));
      double phi = tow.phi();
      double px = sin(theta)*cos(phi);
      double py = sin(theta)*sin(phi);
      double pz = cos(theta);
      TLorentzVector constituentbase(px,py,pz,1);
      EsumVectorL1Nose += tow.energy()*constituentbase;
      EsumVectorL1NoseEM += tow.etEm()*constituentbase;
      EsumVectorL1NoseHAD += tow.etHad()*constituentbase;

    } 
  }

  hResponseTow->Fill((EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);
  hEemTow->Fill(EsumVectorL1NoseEM.Energy());
  hEhadTow->Fill(EsumVectorL1NoseHAD.Energy());

  //  if(maxECLuster/trueE<0.01) hRecHitPositionNull->Fill(initialP4.eta(),initialP4.phi());

  return EsumVectorL1Nose.Energy();

}


void L1Analyzer::getSingle(edm::Handle<reco::GenParticleCollection> genParticles, /*std::vector<CaloParticle> caloParticles,*/ const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const l1t::HGCalTowerBxCollection & towers, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  doSingle=true;

  if(true) {
  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){
    
    // 22 photon, 130 k0L , 211 pi
    //    bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==130 );
    bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==pdgId_);
    //    bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==211 );
    if (!isPhoton) continue;
    
    /*
      bool isGluon = ( genpart->status()==23 and abs(genpart->pdgId())==21 and abs(genpart->eta())<4.2 and abs(genpart->eta())>3  );
      if (!isGluon) continue;
      bool isQuark = ( genpart->status()==23 and (abs(genpart->pdgId())==1 or abs(genpart->pdgId())==2 or abs(genpart->pdgId())==3) and abs(genpart->eta())<4.2 and abs(genpart->eta())>3  );
      if (!isQuark) continue;
    */
    
    //   cout << ">>>>>>> pid,status,px,py,pz,e,eta,phi= "  << genpart->pdgId() << " , " << genpart->status() << " , " << genpart->px() << " , " << genpart->py() << " , " << genpart->pz() << " , " << genpart->energy() << " , ( " << genpart->eta() << " , " << genpart->phi() << ") " << endl;

   hGenParticleEta->Fill(abs(genpart->eta()));
   hGenParticlePt->Fill(genpart->pt());
   hGenParticleE->Fill(genpart->energy());
   
   //   analyzeL1(L12dNose,L13dNose,towers, genpart->p4(), genpart->energy(), geom, geomHcal);
   double E2cl = analyze2Dcl(L12dNose, genpart->p4(), genpart->energy(), geom, geomHcal);
   double E3cl = analyze3Dcl(L13dNose, genpart->p4(), genpart->energy(), geom, geomHcal);
   double Etow = analyzeTowers(towers, genpart->p4(), genpart->energy(), geom, geomHcal);
   
  }
  }
}

//
// constructors and destructor
//
L1Analyzer::L1Analyzer(const edm::ParameterSet& iConfig)
  : pdgId_(iConfig.getUntrackedParameter<int>("pdgId"))
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))
{

  usesResource("TFileService");  
  genParticlesTag_ = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("GenParticleTag", edm::InputTag("genParticles")));
  //  caloParticleTag_ = mayConsume<std::vector<CaloParticle>>(edm::InputTag("mix","MergedCaloTruth"));

  clusters_token_ = consumes<l1t::HGCalClusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("Clusters", edm::InputTag("hgcalBackEndLayer1ProducerHFNose:HGCalBackendLayer1Processor2DClustering")));
  multiclusters_token_ = consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("Multiclusters", edm::InputTag("hgcalBackEndLayer2ProducerHFNose:HGCalBackendLayer2Processor3DClustering")));

  towers_token_ = consumes<l1t::HGCalTowerBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("Towers", edm::InputTag("hgcalTowerProducerHFNose:HGCalTowerProcessor")));

  /*
    cout << "---------------------------------------------------" << endl; 
    cout << "-------------------  Token done  ------------------" << endl;
    cout << "-----------------  book-histo done  ---------------" << endl;
    cout << "---------------------------------------------------" << endl; 
  */   

  hGenParticleEta = FileService->make<TH1F>("hGenParticleEta","hGenParticleEta", 100, 0. , 5.);
  hGenParticlePt = FileService->make<TH1F>("hGenParticlePt","hGenParticlePt", 100, 0. , 100.);
  hGenParticleE = FileService->make<TH1F>("hGenParticleE","hGenParticleE", 100, 0. , 1000.);

  //  hClusterLayerRho = FileService->make<TH2F>("hClusterLayerRho","hClusterLayerRho", 10, 0. , 10., 150, 0, 150);
  hClusterEnergy = FileService->make<TH1F>("hClusterEnergy","hClusterEnergy", 100, 0. , 100.);
  hClusterEtaPhi = FileService->make<TH2F>("hClusterEtaPhi","hClusterEtaPhi", 100, -5., 5., 628, -3.14, 3.14);
  hClusterLayer = FileService->make<TH1F>("hClusterLayer","hClusterLayer", 10, 0. , 10.);
  hClusterSize = FileService->make<TH1F>("hClusterConstituents","hClusterConstituents", 10, 0. , 10.);

  h3DClusterEnergy = FileService->make<TH1F>("h3DClusterEnergy","h3DClusterEnergy", 100, 0. , 100.);
  h3DClusterEtaPhi = FileService->make<TH2F>("h3DClusterEtaPhi","h3DClusterEtaPhi", 100, -5., 5., 628, -3.14, 3.14);
  h3DClusterLayer = FileService->make<TH1F>("h3DClusterLayer","h3DClusterLayer", 10, 0. , 10.);
  h3DClusterSize = FileService->make<TH1F>("h3DClusterConstituents","h3DClusterConstituents", 10, 0. , 10.);

  hTowClusterEnergy = FileService->make<TH1F>("hTowClusterEnergy","hTowClusterEnergy", 100, 0. , 100.);
  hEemTow = FileService->make<TH1F>("hTowEmEnergy","hTowEmEnergy", 100, 0. , 100.);
  hEhadTow = FileService->make<TH1F>("hTowHadEnergy","hTowHadEnergy", 100, 0. , 100.);

  hRecHitPositionNull = FileService->make<TH2F>("hRecHitPositionNull","hRecHitPositionNull", 100, -5., 5., 628, -3.14, 3.14);
  hResponse = FileService->make<TH1F>("hResponse","hResponse", 100, 0. , 2.);
  hResponse2D = FileService->make<TH1F>("hResponse2D","hResponse2D", 100, 0. , 2.);
  hResponse3D = FileService->make<TH1F>("hResponse3D","hResponse3D", 100, 0. , 2.);
  hResponseTow = FileService->make<TH1F>("hResponseTow","hResponseTow", 100, 0. , 2.);
  hMax2DClResponse = FileService->make<TH1F>("hMax2DClResponse","hMax2DClResponse", 100, 0. , 2.);
  hMax3DClResponse = FileService->make<TH1F>("hMax3DClResponse","hMax3DClResponse", 100, 0. , 2.);

  hClSigmaEtaEtaTot = FileService->make<TH1F>("hClSigmaEtaEtaTot","hClSigmaEtaEtaTot", 100, 0. , 0.05);
  hClSigmaEtaEtaMax = FileService->make<TH1F>("hClSigmaEtaEtaMax","hClSigmaiEtaiEtaMax", 100, 0. , 0.05);
  hClSigmaPhiPhiTot = FileService->make<TH1F>("hClSigmaPhiPhiTot","hClSigmaEtaEtaTot", 100, 0. , 0.05);
  hClSigmaPhiPhiMax = FileService->make<TH1F>("hClSigmaPhiPhiMax","hClSigmaPhiPhiMax", 100, 0. , 0.05);
  hClSigmaRRTot = FileService->make<TH1F>("hClSigmaRRTot","hClSigmaRRTot", 100, 0. , 0.005);
  hClSigmaRRMax = FileService->make<TH1F>("hClSigmaRRMax","hClSigmaRRMax", 100, 0. , 0.005);
  hClZZ = FileService->make<TH1F>("hClSigmaZZ","hClSigmaZZ", 100, 0. , 12.5);
  hHoverE = FileService->make<TH1F>("hHoverE","hHoverE", 100, -2. , 4.);

}


L1Analyzer::~L1Analyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

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

  /*
  edm::Handle<std::vector<CaloParticle>> caloParticlesH;
  iEvent.getByToken(caloParticleTag_, caloParticlesH);
  std::vector<CaloParticle>  caloParticles = *(caloParticlesH.product()); 
  */

  edm::Handle<l1t::HGCalClusterBxCollection> L12dNose;
  iEvent.getByToken(clusters_token_, L12dNose);
  const l1t::HGCalClusterBxCollection noseL12d = *(L12dNose.product());

  edm::Handle<l1t::HGCalMulticlusterBxCollection> L13dNose;
  iEvent.getByToken(multiclusters_token_, L13dNose);
  const l1t::HGCalMulticlusterBxCollection noseL13d = *(L13dNose.product());

  /*  

  edm::ESHandle<HGCalTriggerGeometryBase> L1geometryHandle;
  iSetup.get<CaloGeometryRecord>().get(L1geometryHandle);
  const HGCalTriggerGeometryBase* geomTrigger = (L1geometryHandle.product());

  */

  edm::Handle<l1t::HGCalTowerBxCollection> towers_h;
  iEvent.getByToken(towers_token_, towers_h);
  const l1t::HGCalTowerBxCollection& towers = *towers_h;

  /*
  for (auto tower_itr = towers.begin(0); tower_itr != towers.end(0); tower_itr++) {
    if(tower_itr->pt() > 0.5) std::cout << " pt = " << tower_itr->pt() << " eta " << tower_itr->eta() << " energy = " << tower_itr->energy() << " etEm = " << tower_itr->etEm() << " etHad = " << tower_itr->etHad()  << std::endl;
  }
  */



   /*
   cout << "------------------------------------------------" << endl; 
   cout << "-----------------  Handles done ----------------" << endl;
   cout << "------------------------------------------------" << endl; 
   */   

  if(runSingle) getSingle(genParticles, /*caloParticles,*/ noseL12d, noseL13d, towers, geom, geoHcal );

   return;

}


// ------------ method called once each job just before starting event loop  ------------
void
L1Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
L1Analyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(L1Analyzer);
