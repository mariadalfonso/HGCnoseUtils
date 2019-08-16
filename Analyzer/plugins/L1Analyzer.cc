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
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"

// L1
#include "DataFormats/L1THGCal/interface/HGCalCluster.h"
#include "DataFormats/L1THGCal/interface/HGCalMulticluster.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerTools.h"
#include "L1Trigger/L1THGCal/interface/HGCalTriggerGeometryBase.h"


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

  void analyzeL1(const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const math::XYZTLorentzVector & initialP4 , const double & jetE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);

  void getSingle(edm::Handle<reco::GenParticleCollection> genParticles, const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal);


  // ----------member data ---------------------------
  //  void ClearVariables();

  TH1F *hGenParticleEta;
  TH1F *hGenParticlePt;
  TH1F *hGenParticleE;

  bool verbose = false;

  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesTag_;

  edm::EDGetToken clusters_token_, multiclusters_token_;
  HGCalTriggerTools triggerTools_;

  //https://github.com/cms-sw/cmssw/blob/master/L1Trigger/L1THGCalUtilities/plugins/ntuples/HGCalTriggerNtupleHGCClusters.cc

  /*
  edm::EDGetTokenT<reco::GenJetCollection> genJetsTag_;
  edm::EDGetTokenT<edm::PCaloHitContainer> simHitTag_;
  edm::EDGetTokenT<edm::PCaloHitContainer> simHcalTag_;
  edm::EDGetTokenT<HGCalDigiCollection> digiNose_;
  edm::EDGetTokenT<HGCRecHitCollection> recHitNose_;
  
  edm::EDGetTokenT<HFRecHitCollection> tok_hf_;

  TH1F *hLepEta;
  TH1F *hLepPt;
  TH1F *hLepCharge;

  TH1F *hJetEta;
  TH1F *hJetEtaLarge;

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

  TH1F *hNmips;
  TH1F *hEnergyMIPS;
  TH1F *hJetResponseHF;
  TH1F *hJetResponseNose;

  TH1F *hJetResponseNoseEE;
  TH1F *hJetResponseNoseEH;

  TH2F *hJetResponseHFEta;
  TH2F *hJetResponseNoseEta;

  */

  TH2F *hClusterEtaPhi;
  TH1F *hClusterEnergy;
  TH1F *hClusterLayer;

  TH1F *hJetResponse;
  TH2F *hJetResponseE;
  TH2F *hJetResponseEta;

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


void L1Analyzer::analyzeL1(const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const math::XYZTLorentzVector & initialP4, const double & trueE, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  TLorentzVector EsumVectorL1Nose(0.,0.,0.,0.);
  TLorentzVector EsumVectorL1HF(0.,0.,0.,0.);

  for (auto const& cl : L12dNose) {

    double   thisDeltaR1 = ::deltaR(cl.eta(), cl.phi() , initialP4.eta(), initialP4.phi());

    //    cout << " thisDeltaR = " <<  thisDeltaR1 <<    "( " << cl.eta() << " " << cl.phi() << "   ) "<< endl;

    if ( thisDeltaR1 < coneSize ) {

      //      cout << " layerWithOffset " << triggerTools_.layerWithOffset(cl.detId()) << endl;

    hClusterEnergy->Fill(cl.energy());
    hClusterEtaPhi->Fill(cl.eta(),cl.phi());
    hClusterLayer->Fill(triggerTools_.layerWithOffset(cl.detId()));

    double theta = 2 * atan(exp(cl.eta()));
    double phi = cl.phi();
    double px = sin(theta)*cos(phi);
    double py = sin(theta)*sin(phi);
    double pz = cos(theta);
    TLorentzVector constituentbase(px,py,pz,1);
    EsumVectorL1Nose += cl.energy()*constituentbase;

    }
    //    cout << cluster.energy() << endl;

  }

  /*
  hJetResponseE->Fill(trueE, (EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);
  hJetResponseEta->Fill(initialP4.eta(), (EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);
  */
  hJetResponse->Fill((EsumVectorL1Nose+EsumVectorL1HF).Energy()/trueE);

}

void L1Analyzer::getSingle(edm::Handle<reco::GenParticleCollection> genParticles, const l1t::HGCalClusterBxCollection& L12dNose, const l1t::HGCalMulticlusterBxCollection& L13dNose, const HGCalGeometry* geom, const CaloSubdetectorGeometry *geomHcal) {

  doSingle=true;

  for(reco::GenParticleCollection::const_iterator genpart = genParticles->begin(); genpart != genParticles->end(); ++genpart){

   // 22 photon, 130 k0L , 211 pi
   //   bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==130 );
    bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==22 );
   //   bool isPhoton = ( genpart->isPromptFinalState() and abs(genpart->pdgId())==211 );
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
   
   analyzeL1(L12dNose,L13dNose,genpart->p4(), genpart->energy(), geom, geomHcal);
   
 }

}

//
// constructors and destructor
//
L1Analyzer::L1Analyzer(const edm::ParameterSet& iConfig)
// :
  //  tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

{

  usesResource("TFileService");  
  genParticlesTag_ = consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("GenParticleTag", edm::InputTag("genParticles")));

  clusters_token_ = consumes<l1t::HGCalClusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("Clusters", edm::InputTag("hgcalBackEndLayer1Producer:HGCalBackendLayer1Processor2DClustering")));
  multiclusters_token_ = consumes<l1t::HGCalMulticlusterBxCollection>(iConfig.getUntrackedParameter<edm::InputTag>("Multiclusters", edm::InputTag("hgcalBackEndLayer2Producer:HGCalBackendLayer2Processor3DClustering")));



  /*                                                                                                                                                                                 
Type                                  Module                      Label             Process
----------------------------------------------------------------------------------------------                                                                                       
  BXVector<l1t::HGCalCluster>           "hgcalBackEndLayer1Producer"   "HGCalBackendLayer1Processor2DClustering"   "DIGI"
  BXVector<l1t::HGCalMulticluster>      "hgcalBackEndLayer2Producer"   "HGCalBackendLayer2Processor3DClustering"   "DIGI"
  BXVector<l1t::HGCalTower>             "hgcalTowerProducer"        "HGCalTowerProcessor"   "DIGI"


  edm::SortedCollection<HGCRecHit,edm::StrictWeakOrdering<HGCRecHit> >    "HGCalRecHit"               "HGCHFNoseRecHits"   "RECO"

  */


  /*
  genJetsTag_ = consumes<reco::GenJetCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("GenJetTag", edm::InputTag("ak4GenJets")));

  simHitTag_ = mayConsume<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HFNoseHits"));
  //  simHcalTag_ = mayConsume<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  simHcalTag_ = consumes<edm::PCaloHitContainer>(edm::InputTag("g4SimHits","HcalHits"));
  tok_hf_ = consumes<HFRecHitCollection>(iConfig.getUntrackedParameter<string>("HFRecHits","hfreco"));

  //  digiNose_ = consumes<HGCalDigiCollection>(iConfig.getParameter<edm::InputTag>("simHFNoseUnsuppressedDigis","HFNose"));
  digiNose_ = consumes<HGCalDigiCollection>(iConfig.getUntrackedParameter<edm::InputTag>("DIGITAG",edm::InputTag("simHFNoseUnsuppressedDigis:HFNose")));

  recHitNose_ = consumes<HGCRecHitCollection>(iConfig.getUntrackedParameter<edm::InputTag>("REChitTAG",edm::InputTag("HGCalRecHit:HGCHFNoseRecHits")));
  */

  /*
    cout << "---------------------------------------------------" << endl; 
    cout << "-------------------  Token done  ------------------" << endl;
    cout << "-----------------  book-histo done  ---------------" << endl;
    cout << "---------------------------------------------------" << endl; 
  */   

  hGenParticleEta = FileService->make<TH1F>("hGenParticleEta","hGenParticleEta", 100, 0. , 5.);
  hGenParticlePt = FileService->make<TH1F>("hGenParticlePt","hGenParticlePt", 100, 0. , 100.);
  hGenParticleE = FileService->make<TH1F>("hGenParticleE","hGenParticleE", 100, 0. , 1000.);

  hClusterEtaPhi = FileService->make<TH2F>("hClusterEtaPhi","hClusterEtaPhi", 100, -5., 5., 628, -3.14, 3.14);
  hClusterEnergy = FileService->make<TH1F>("hClusterEnergy","hClusterEnergy", 100, 0. , 100.);
  hClusterLayer = FileService->make<TH1F>("hClusterLayer","hClusterLayer", 10, 0. , 10.);

  hJetResponse = FileService->make<TH1F>("hJetResponse","hJetResponse", 100, 0. , 2.);

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

   /*
   cout << "------------------------------------------------" << endl; 
   cout << "-----------------  Handles done ----------------" << endl;
   cout << "------------------------------------------------" << endl; 
   */   

  if(runSingle) getSingle(genParticles, noseL12d, noseL13d, geom, geoHcal );

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
