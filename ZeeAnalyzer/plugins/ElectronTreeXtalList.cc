// system include files
#include <memory>
#include <vector>
#include <fstream>
#include <map>
#include <iostream>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"


#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"


#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#define NELE 3
#define initFloat     {-999.,-999.,-999.}
#define initInt       {0,0,0}
#define initIntCharge {-100,-100,-100}

class ElectronTreeXtalList : public edm::EDAnalyzer {
public:
  explicit ElectronTreeXtalList(const edm::ParameterSet&);
  ~ElectronTreeXtalList();
  
  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON,// The last does not include tau parents
			  TRUE_ELECTRON_FROM_Z,
			  TRUE_ELECTRON_FROM_GAMMA, 
                          };
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  std::pair<int,double> matchToTruth(const reco::GsfElectron* el, 
		   edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
  
  void findFirstNonElectronMother(const reco::Candidate *particle,
				  int &ancestorPID, int &ancestorStatus);
  const float getEffectiveArea(float eta) const;
  
  edm::EDGetTokenT<double> theRhoToken;
  edm::Handle< double > _rhoHandle;

  // ----------member data ---------------------------
  // edm::LumiReWeighting LumiWeights_;  

  std::map<std::string, std::vector<long int> > xtalMap;

  bool isRecovered2;
  bool isDead2;
  bool isRecovered;
  bool isDead;

  bool isMC=false; 
  // Data members
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventInfoProduct_;
  edm::EDGetToken electronsToken_;
  edm::EDGetToken patElectronsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> rechits_EB_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_ ;

  TFile *tree_file;
  TTree * tree; 

  UInt_t     	runNumber;		///< run number
  UShort_t      lumiBlock;		///< lumi section
  Long64_t    eventNumber;	///< event number
  UInt_t	eventTime;		///< unix time of the event
  UShort_t	      nBX;			///< bunch crossing
  int    	nTruePU;		///< pu
  
  float  MyWeight;
  Float_t     mcGenWeight; ///< weight in generator for MC
    
  int e1Charge; 
  int e2Charge; 
  
  Float_t   e1Eta;
  Float_t   e1Phi;
  Float_t   e1R9;
  Float_t   e2Eta;
  Float_t   e2Phi;
  Float_t   e2R9;
  float   e1SeedIEta;
  float   e1SeedIPhi;
  float   e2SeedIEta;
  float   e2SeedIPhi;
  float e1SigmaIetaIeta;
  float e2SigmaIetaIeta;
  Float_t   e1EtaSC;
  Float_t   e2EtaSC;
  Float_t   e1PhiSC;
  Float_t   e2PhiSC;
  
  Float_t e1Energy;
  Float_t e1RawEnergy;
  Float_t e1Energy3x3SC;
  Float_t e1Energy5x5SC;
  Float_t e2Energy;
  Float_t e2RawEnergy;
  Float_t e2Energy3x3SC;
  Float_t e2Energy5x5SC;
  
  Float_t e1Pt; 
  Float_t e2Pt; 

  Float_t invMass;
  Float_t invMass_rawSC;
  
  Float_t  e1GenEnergy, e2GenEnergy;
  Float_t  invMass_MC;

  Int_t e1IsDead, e1IsRecovered;
  Int_t e2IsDead, e2IsRecovered;

  std::vector<bool> vkGood,    vkPoorReco,    vkOutOfTimE,    vkFaultyHardware,    vkNoisy,    vkPoorCalib,    vkSaturated,    vkLeadingEdgeRecovered,    vkNeighboursRecovered,    vkTowerRecovered,    vkDead,    vkKilled,    vkTPSaturated,    vkL1SpikeFlag,    vkWeird,    vkDiWeird,    vkHasSwitchToGain6,    vkHasSwitchToGain1,    vkUnknown;
  std::vector<float> vXtalEn, vSum8, vAve, vNeigh0, vNeigh1, vNeigh2, vNeigh3, vNeigh5, vNeigh6, vNeigh7, vNeigh8,vEnFr;
  std::vector<int> vIc, vIeta, vIphi, vIsm;
  std::vector<unsigned long int> vRawId;

  // Histos
  TH1F *h_eta, *h_phi, *h_isTrue, *h_eID;
  TH1F *h_EB_pt, *h_EE_pt;
  TH1F *h_EB_rawEne, *h_EE_rawEne;
  TH1F *h_EB_sigmaIeIe, *h_EE_sigmaIeIe;
  TH1F *h_EB_r9, *h_EE_r9;
  TH1F *h_EB_r9uz, *h_EE_r9uz;
  TH1F *h_EB_hoe, *h_EE_hoe;
  TH1F *h_EB_dz, *h_EE_dz;
  TH1F *h_EB_conv, *h_EE_conv;
  TH1F *hAll_eta, *hAll_phi, *hAll_isTrue, *hAll_eID;
  TH1F *hAll_EB_pt, *hAll_EE_pt;
  TH1F *hAll_EB_rawEne, *hAll_EE_rawEne;
  TH1F *hAll_EB_sigmaIeIe, *hAll_EE_sigmaIeIe;
  TH1F *hAll_EB_r9, *hAll_EE_r9;
  TH1F *hAll_EB_r9uz, *hAll_EE_r9uz;
  TH1F *hAll_EB_hoe, *hAll_EE_hoe;
  TH1F *hAll_EB_dz, *hAll_EE_dz;
  TH1F *hAll_EB_conv, *hAll_EE_conv;
  TH1F *hRecovered_eta, *hRecovered_phi, *hRecovered_isTrue, *hRecovered_eID;
  TH1F *hRecovered_EB_pt, *hRecovered_EE_pt;
  TH1F *hRecovered_EB_rawEne, *hRecovered_EE_rawEne;
  TH1F *hRecovered_EB_sigmaIeIe, *hRecovered_EE_sigmaIeIe;
  TH1F *hRecovered_EB_r9, *hRecovered_EE_r9;
  TH1F *hRecovered_EB_r9uz, *hRecovered_EE_r9uz;
  TH1F *hRecovered_EB_hoe, *hRecovered_EE_hoe;
  TH1F *hRecovered_EB_dz, *hRecovered_EE_dz;
  TH1F *hRecovered_EB_conv, *hRecovered_EE_conv;

  TH1F * hRecovered_EB_rawGenEnRatio;
  TH1F * hRecovered_EB_eGenEn;
  TH1F * hRecovered_EB_eEnGenEnRatio;
  TH1F * h_EB_rawGenEnRatio;
  TH1F * h_EB_eGenEn;
  TH1F * h_EB_eEnGenEnRatio;
  TH1F * hAll_EB_rawGenEnRatio;
  TH1F * hAll_EB_eGenEn;
  TH1F * hAll_EB_eEnGenEnRatio;

  TH1F * h_Flags, *h_FlagsGood; 

  TH1F *h_ZeeMass, *h_ZMassRawEn,*hAll_ZeeMass, *hAll_ZMassRawEn, *hRecovered_ZeeMass, *hRecoverd_ZMassRawEn, *hEventCount; 


};

ElectronTreeXtalList::ElectronTreeXtalList(const edm::ParameterSet& iConfig) {

  isMC=iConfig.getParameter<bool>("isMC");

  beamSpotToken_    = consumes<reco::BeamSpot> 
    (iConfig.getParameter <edm::InputTag>
     ("beamSpot"));

  genEventInfoProduct_ = consumes<GenEventInfoProduct> 
    (iConfig.getParameter <edm::InputTag>
     ("genEventInfoProduct"));

  electronsToken_    = mayConsume<reco::GsfElectronCollection> 
    (iConfig.getParameter<edm::InputTag>
     ("electrons"));
  // patElectronsToken_    = mayConsume<pat::ElectronCollection >
  //   (iConfig.getParameter<edm::InputTag>
  //    ("patelectrons"));

  genParticlesToken_ = mayConsume<edm::View<reco::GenParticle> >
    (iConfig.getParameter<edm::InputTag>
     ("genParticles"));

  vtxToken_          = mayConsume<reco::VertexCollection>
    (iConfig.getParameter<edm::InputTag>
     ("vertices"));

  conversionsToken_ = mayConsume< reco::ConversionCollection >
    (iConfig.getParameter<edm::InputTag>
     ("conversions"));

  rechits_EB_=consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("inputRecHitsEB"));

  pileupSummaryToken_ = consumes<std::vector<PileupSummaryInfo> >(edm::InputTag("addPileupInfo"));

  theRhoToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralCalo"));

  edm::Service<TFileService> fs;


 
	//now do what ever initialization is needed
  tree = fs->make<TTree>("selected", "selected");
	//tree = fs->make<TTree>("selected","selected"); //no otherwise you have the extraCalib in the same file


  tree->Branch("runNumber",     &runNumber,   "runNumber/i");
  tree->Branch("lumiBlock",     &lumiBlock,   "lumiBlock/s");
  tree->Branch("eventNumber",   &eventNumber, "eventNumber/l");
  tree->Branch("eventTime",     &eventTime,   "eventTime/i");
  tree->Branch("nBX",           &nBX,         "nBX/s");
  tree->Branch("nTruePU",       &nTruePU,         "nTruePU/I");
  tree->Branch("puweight",       &MyWeight,         "puWeight/F");
  

  // ele
  tree->Branch("e1Charge",   &e1Charge,    "e1Charge/I");
  tree->Branch("e1Eta",      &e1Eta,       "e1Eta/F");
  tree->Branch("e1Phi",      &e1Phi,       "e1Phi/F");
  tree->Branch("e1R9",       &e1R9,       "e1R9/F");
  tree->Branch("e2Charge",   &e2Charge,    "e2Charge/I");
  tree->Branch("e2Eta",      &e2Eta,       "e2Eta/F");
  tree->Branch("e2Phi",      &e2Phi,       "e2Phi/F");
  tree->Branch("e2R9",       &e2R9,       "e2R9/F");
  tree->Branch("e1SigmaIetaIeta", &e1SigmaIetaIeta, "e1SigmaIetaIeta/F");
  tree->Branch("e2SigmaIetaIeta", &e2SigmaIetaIeta, "e2SigmaIetaIeta/F");

  tree->Branch("e1GenEnergy", &e1GenEnergy, "e1GenEnergy/F");
  tree->Branch("e2GenEnergy", &e2GenEnergy, "e2GenEnergy/F");
	// SC
  tree->Branch("e1EtaSC",      &e1EtaSC,       "e1EtaSC/F");
  tree->Branch("e1PhiSC",      &e1PhiSC,       "e1PhiSC/F");
  tree->Branch("e1Energy",     &e1Energy,      "e1EnergySC/F");
  tree->Branch("e1RawEnergy", &e1RawEnergy, "e1RawEnergySC/F");
  // tree->Branch("e1energy5x5SC", e1energy5x5SC, "e1energy5x5SC/F");
  // tree->Branch("e1energy3x3SC", e1energy3x3SC, "e1energy3x3SC/F");
  tree->Branch("e2EtaSC",      &e2EtaSC,       "e2EtaSC/F");
  tree->Branch("e2PhiSC",      &e2PhiSC,       "e2PhiSC/F");
  tree->Branch("e2Energy",     &e2Energy,      "e2EnergySC/F");
  tree->Branch("e2RawEnergy", &e2RawEnergy, "e2RawEnergySC/F");
  // tree->Branch("e2energy5x5SC", e2energy5x5SC, "e2energy5x5SC/F");
  // tree->Branch("e2energy3x3SC", e2energy3x3SC, "e2energy3x3SC/F");

  tree->Branch("invMass", &invMass, "invMass/F");
  tree->Branch("invMass_rawSC", &invMass_rawSC,   "invMass_rawSC/F");

  tree->Branch("e1IsDead", &e1IsDead, "e1IsDead/I");
  tree->Branch("e1IsRecovered", &e1IsRecovered, "e1IsRecovered/I");
  tree->Branch("e2IsDead", &e2IsDead, "e2IsDead/I");
  tree->Branch("e2IsRecovered", &e2IsRecovered, "e2IsRecovered/I");

  tree->Branch("e1SeedIEta", &e1SeedIEta, "e1SeedIEta/F");
  tree->Branch("e2SeedIEta", &e2SeedIEta, "e2SeedIEta/F");
  tree->Branch("e1SeedIPhi", &e1SeedIPhi, "e1SeedIPhi/F");
  tree->Branch("e2SeedIPhi", &e2SeedIPhi, "e2SeedIPhi/F");

  tree->Branch("kGood"                ,&vkGood                );
  tree->Branch("kPoorReco"            ,&vkPoorReco            );
  tree->Branch("kOutOfTimE"           ,&vkOutOfTimE           );
  tree->Branch("kFaultyHardware"      ,&vkFaultyHardware      );
  tree->Branch("kNoisy"               ,&vkNoisy               );
  tree->Branch("kPoorCalib"           ,&vkPoorCalib           );
  tree->Branch("kSaturated"           ,&vkSaturated           );
  tree->Branch("kLeadingEdgeRecovered",&vkLeadingEdgeRecovered);
  tree->Branch("kNeighboursRecovered" ,&vkNeighboursRecovered );
  tree->Branch("kTowerRecovered"      ,&vkTowerRecovered      );
  tree->Branch("kDead"                ,&vkDead                );
  tree->Branch("kKilled"              ,&vkKilled              );
  tree->Branch("kTPSaturated"         ,&vkTPSaturated         );
  tree->Branch("kL1SpikeFlag"         ,&vkL1SpikeFlag         );
  tree->Branch("kWeird"               ,&vkWeird               );
  tree->Branch("kDiWeird"             ,&vkDiWeird             );
  tree->Branch("kHasSwitchToGain6"    ,&vkHasSwitchToGain6    );
  tree->Branch("kHasSwitchToGain1"    ,&vkHasSwitchToGain1    );
  tree->Branch("kUnknown"             ,&vkUnknown             );  

  tree->Branch("xtalIeta",&vIeta );  
  tree->Branch("xtalIphi",&vIphi );  
  tree->Branch("xtalIsm",&vIsm );  
  tree->Branch("xtalIc",&vIc );  
  tree->Branch("xtalRawId",&vRawId );  
  tree->Branch("xtalEn",&vXtalEn);
  tree->Branch("vEnFr", &vEnFr); 


  tree->Branch("vSum8"     , &vSum8);	 
  tree->Branch("vAve"	    , &vAve);	 
  tree->Branch("vNeigh0"   , &vNeigh0);
  tree->Branch("vNeigh1"   , &vNeigh1);
  tree->Branch("vNeigh2"   , &vNeigh2);
  tree->Branch("vNeigh3"   , &vNeigh3);
  tree->Branch("vNeigh5"   , &vNeigh5);
  tree->Branch("vNeigh6"   , &vNeigh6);
  tree->Branch("vNeigh7"   , &vNeigh7);
  tree->Branch("vNeigh8"   , &vNeigh8);

  h_Flags= fs->make<TH1F>("hFlags", "hFlags", 20, 0, 20);
  h_FlagsGood= fs->make<TH1F>("hFlagsGood", "hFlagsGood", 20, 0, 20);

  hEventCount = fs->make<TH1F>("hEventCount", "hEventCount", 3, -0.5, 2.5);
  h_eta = fs->make<TH1F>("eta", "eta", 50,-2.5,2.5);
  h_phi = fs->make<TH1F>("phi", "phi", 50,-2.5,2.5);
  h_isTrue = fs->make<TH1F>("isTrue", "isTrue", 2,-0.5,1.5);
  h_eID = fs->make<TH1F>("eIDSummer16VetoWP", "eIDSummer16VetoWP", 2, -0.5, 1.5);
  h_EB_pt = fs->make<TH1F>("EB_pt", "EB_pt", 50,0.,200.);
  h_EE_pt = fs->make<TH1F>("EE_pt", "EE_pt", 50,0.,200.);
  h_EB_rawEne = fs->make<TH1F>("EB_rawEne", "EB_rawEne", 500,0.,500.);
  h_EE_rawEne = fs->make<TH1F>("EE_rawEne", "EE_rawEne", 100,0.,400.);
  h_EB_sigmaIeIe = fs->make<TH1F>("EB_sigmaIeIe", "EB_sigmaIeIe", 50,0.,0.05);
  h_EE_sigmaIeIe = fs->make<TH1F>("EE_sigmaIeIe", "EE_sigmaIeIe", 50,0.,0.1);
  h_EB_r9 = fs->make<TH1F>("EB_r9", "EB_r9", 50,0.5,1.);
  h_EE_r9 = fs->make<TH1F>("EE_r9", "EE_r9", 50,0.5,1.);
  h_EB_r9uz = fs->make<TH1F>("EB_r9uz", "EB_r9uz", 100,0.,1.);
  h_EE_r9uz = fs->make<TH1F>("EE_r9uz", "EE_r9uz", 100,0.,1.);
  h_EB_hoe = fs->make<TH1F>("EB_hoe", "EB_hoe", 50,0.,0.1);
  h_EE_hoe = fs->make<TH1F>("EE_hoe", "EE_hoe", 50,0.,0.1);
  h_EB_dz = fs->make<TH1F>("EB_dz", "EB_dz", 50,0.,0.1);
  h_EE_dz = fs->make<TH1F>("EE_dz", "EE_dz", 50,0.,0.1);
  h_EB_conv = fs->make<TH1F>("EB_conv", "EB_conv", 2,-0.5,1.5);
  h_EE_conv = fs->make<TH1F>("EE_conv", "EE_conv", 2,-0.5,1.5); 

  hAll_eta = fs->make<TH1F>("all_eta", "eta all", 50,-2.5,2.5);
  hAll_phi = fs->make<TH1F>("all_phi", "phi all", 50,-2.5,2.5);
  hAll_isTrue = fs->make<TH1F>("all_isTrue", "isTrue all", 2,-0.5,1.5);
  hAll_eID = fs->make<TH1F>("all_eIDSummer16VetoWP", "eIDSummer16VetoWP all", 2, -0.5, 1.5);
  hAll_EB_pt = fs->make<TH1F>("all_EB_pt", "EB_pt all", 50,0.,200.);
  hAll_EE_pt = fs->make<TH1F>("all_EE_pt", "EE_pt all", 50,0.,200.);
  hAll_EB_rawEne = fs->make<TH1F>("all_EB_rawEne", "EB_rawEne all", 500,0.,500.);
  hAll_EE_rawEne = fs->make<TH1F>("all_EE_rawEne", "EE_rawEne all", 100,0.,400.);
  hAll_EB_sigmaIeIe = fs->make<TH1F>("all_EB_sigmaIeIe", "EB_sigmaIeIe all", 50,0.,0.05);
  hAll_EE_sigmaIeIe = fs->make<TH1F>("all_EE_sigmaIeIe", "EE_sigmaIeIe all", 50,0.,0.1);
  hAll_EB_r9 = fs->make<TH1F>("all_EB_r9", "EB_r9 all", 50,0.5,1.);
  hAll_EE_r9 = fs->make<TH1F>("all_EE_r9", "EE_r9 all", 50,0.5,1.);
  hAll_EB_r9uz = fs->make<TH1F>("all_EB_r9uz", "EB_r9uz all", 100,0.,1.);
  hAll_EE_r9uz = fs->make<TH1F>("all_EE_r9uz", "EE_r9uz all", 100,0.,1.);
  hAll_EB_hoe = fs->make<TH1F>("all_EB_hoe", "EB_hoe all", 50,0.,0.1);
  hAll_EE_hoe = fs->make<TH1F>("all_EE_hoe", "EE_hoe all", 50,0.,0.1);
  hAll_EB_dz = fs->make<TH1F>("all_EB_dz", "EB_dz all", 50,0.,0.1);
  hAll_EE_dz = fs->make<TH1F>("all_EE_dz", "EE_dz all", 50,0.,0.1);
  hAll_EB_conv = fs->make<TH1F>("all_EB_conv", "EB_conv all", 2,-0.5,1.5);
  hAll_EE_conv = fs->make<TH1F>("all_EE_conv", "EE_conv all", 2,-0.5,1.5); 

  hRecovered_eta = fs->make<TH1F>("recovered_eta", "eta", 50,-2.5,2.5);
  hRecovered_phi = fs->make<TH1F>("recovered_phi", "phi", 50,-2.5,2.5);
  hRecovered_isTrue = fs->make<TH1F>("recovered_isTrue", "isTrue", 2,-0.5,1.5);
  hRecovered_eID = fs->make<TH1F>("recovered_eIDSummer16VetoWP", "eIDSummer16VetoWP", 2, -0.5, 1.5);
  hRecovered_EB_pt = fs->make<TH1F>("recovered_EB_pt", "EB_pt", 50,0.,200.);
  hRecovered_EE_pt = fs->make<TH1F>("recovered_EE_pt", "EE_pt", 50,0.,200.);
  hRecovered_EB_rawEne = fs->make<TH1F>("recovered_EB_rawEne", "EB_rawEne", 500,0.,500.);
  hRecovered_EE_rawEne = fs->make<TH1F>("recovered_EE_rawEne", "EE_rawEne", 100,0.,400.);
  hRecovered_EB_sigmaIeIe = fs->make<TH1F>("recovered_EB_sigmaIeIe", "EB_sigmaIeIe", 50,0.,0.05);
  hRecovered_EE_sigmaIeIe = fs->make<TH1F>("recovered_EE_sigmaIeIe", "EE_sigmaIeIe", 50,0.,0.1);
  hRecovered_EB_r9 = fs->make<TH1F>("recovered_EB_r9", "EB_r9", 50,0.5,1.);
  hRecovered_EE_r9 = fs->make<TH1F>("recovered_EE_r9", "EE_r9", 50,0.5,1.);
  hRecovered_EB_r9uz = fs->make<TH1F>("recovered_EB_r9uz", "EB_r9uz", 100,0.,1.);
  hRecovered_EE_r9uz = fs->make<TH1F>("recovered_EE_r9uz", "EE_r9uz", 100,0.,1.);
  hRecovered_EB_hoe = fs->make<TH1F>("recovered_EB_hoe", "EB_hoe", 50,0.,0.1);
  hRecovered_EE_hoe = fs->make<TH1F>("recovered_EE_hoe", "EE_hoe", 50,0.,0.1);
  hRecovered_EB_dz = fs->make<TH1F>("recovered_EB_dz", "EB_dz", 50,0.,0.1);
  hRecovered_EE_dz = fs->make<TH1F>("recovered_EE_dz", "EE_dz", 50,0.,0.1);
  hRecovered_EB_conv = fs->make<TH1F>("recovered_EB_conv", "EB_conv", 2,-0.5,1.5);
  hRecovered_EE_conv = fs->make<TH1F>("recovered_EE_conv", "EE_conv", 2,-0.5,1.5); 

  h_ZeeMass=fs->make<TH1F>("h_ZeeMass", "Zee Inv Mass", 100, 0, 200);
  h_ZMassRawEn=fs->make<TH1F>("h_ZMassRawEn", "Zee Inv Mass, rawEn", 100, 0, 200);
  hAll_ZeeMass=fs->make<TH1F>("hAll_ZeeMass", "Zee Inv Mass, all", 100, 0, 200);
  hAll_ZMassRawEn=fs->make<TH1F>("hAll_ZMassRawEn", "Zee Inv Mass, all, rawEn", 100, 0, 200);
  hRecovered_ZeeMass=fs->make<TH1F>("hRecovered_ZeeMass", "Zee Inv Mass, recov", 100, 0, 200);
  hRecoverd_ZMassRawEn=fs->make<TH1F>("hRecovered_ZMassRawEn", "Zee Inv Mass, recov, rawEn", 100, 0, 200);

  hRecovered_EB_rawGenEnRatio=fs->make<TH1F>("hRecovered_EB_rawGenEnRatio", "e raw En / e gen En", 100, -5, 5);
  hRecovered_EB_eGenEn=fs->make<TH1F>("hRecovered_EB_eGenEn", "e Enenergy", 100, 0, 100);
  hRecovered_EB_eEnGenEnRatio=fs->make<TH1F>("hRecovered_EB_EnGenEnRatio", "e En / e gen En", 100, -5, 5);
  h_EB_rawGenEnRatio=fs->make<TH1F>("h_EB_rawGenEnRatio", "e raw En / e gen En", 100, -5, 5);
  h_EB_eGenEn=fs->make<TH1F>("h_EB_eGenEn", "e Enenergy", 100, 0, 100);
  h_EB_eEnGenEnRatio=fs->make<TH1F>("h_EB_EnGenEnRatio", "e raw En / e gen En", 100, -5, 5);
  hAll_EB_rawGenEnRatio=fs->make<TH1F>("hAll_EB_rawGenEnRatio", "e raw En / e gen En", 100, -5, 5);
  hAll_EB_eGenEn=fs->make<TH1F>("hAll_EB_eGenEn", "e Enenergy", 100, 0, 100);
  hAll_EB_eEnGenEnRatio=fs->make<TH1F>("hAll_EB_EnGenEnRatio", "e En / e gen En", 100, -5, 5);


  typedef reco::Candidate::LorentzVector LorentzVector;
  
  // LumiWeights_ = edm::LumiReWeighting("/afs/cern.ch/work/t/taroni/private/newDeadCh102X/src/ZeeAnalyzer/ZeeAnalyzer/test/puMC_JuneProjectionFull18_PoissonOOTPU.root",
  //                                     "/afs/cern.ch/work/t/taroni/private/newDeadCh102X/src/ZeeAnalyzer/ZeeAnalyzer/test/MyDataPileupHistogram.root",
  //                                     "puMC",
  //                                     "pileup");


}


ElectronTreeXtalList::~ElectronTreeXtalList() { }
const float ElectronTreeXtalList::getEffectiveArea(float eta) const{
  std::vector<double> absEtaMin_={ 0.0000, 1.0000, 1.4790, 2.0000, 2.2000, 2.3000, 2.4000};
  std::vector<double> absEtaMax_={ 1.0000,  1.4790, 2.0000,  2.2000, 2.3000, 2.4000, 5.0000};
  std::vector<double> effectiveAreaValues_={ 0.1703, 0.1715, 0.1213, 0.1230, 0.1635, 0.1937, 0.2393};


  float effArea = 0;
  uint nEtaBins = absEtaMin_.size();
  for(uint iEta = 0; iEta<nEtaBins; iEta++){
    if( std::abs(eta) >= absEtaMin_[iEta]
    && std::abs(eta) < absEtaMax_[iEta] ){
      effArea = effectiveAreaValues_[iEta];
      break;
    }
  }

  return effArea;
}

void ElectronTreeXtalList::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace cms;
  runNumber=iEvent.id().run();
  lumiBlock=(UShort_t)iEvent.id().luminosityBlock();
  eventNumber=(Long64_t) iEvent.id().event();
  edm::Timestamp time = iEvent.eventAuxiliary().time();
  eventTime= (UInt_t) time.unixTime();
  nBX=(UShort_t) iEvent.bunchCrossing();	

  vector<bool> vkGood_,    vkPoorReco_,    vkOutOfTimE_,    vkFaultyHardware_,    vkNoisy_,    vkPoorCalib_,    vkSaturated_,    vkLeadingEdgeRecovered_,    vkNeighboursRecovered_,    vkTowerRecovered_,    vkDead_,    vkKilled_,    vkTPSaturated_,    vkL1SpikeFlag_,    vkWeird_,    vkDiWeird_,    vkHasSwitchToGain6_,    vkHasSwitchToGain1_,    vkUnknown_;
  vector<float> vXtalEn_,vNeigh0_,vNeigh1_,vNeigh2_,vNeigh3_,vNeigh5_,vNeigh6_,vNeigh7_,vNeigh8_,vSum8_, vAve_,vEnFr_;
  vector<int >  vIeta_, vIphi_, vIsm_, vIc_;
  vector<unsigned long int> vRawId_;
      
  e1Charge=-2.; 
  e2Charge=-2.; 
   
  e1Eta=-99.;
  e1Phi=-99.;
  e1Phi=-99.;
  e1R9=-99.;
  e1SigmaIetaIeta=-99.;
  e2SigmaIetaIeta=-99.;
  
  e1EtaSC=-99.;
  e2EtaSC=-99.;
  e1PhiSC=-99.;
  e2PhiSC=-99.;
  
  e1Energy=-99.;
  e1RawEnergy=-99.;
  e1Energy3x3SC=-99.;
  e1Energy5x5SC=-99.;
  e2Energy=-99.;
  e2RawEnergy=-99.;
  e2Energy3x3SC=-99.;
  e2Energy5x5SC=-99.;
  

  e1SeedIEta=-99.;
  e2SeedIEta=-99.;
  e1SeedIPhi=-599.;
  e2SeedIPhi=-599.;

  e1Pt=-99.; 
  e2Pt=-99.; 

  invMass=-99.;
  invMass_rawSC=-99.;
  
  e1GenEnergy=-99.;
  e2GenEnergy=-99.;
  invMass_MC=-99.;

  e1IsDead=-99.;
  e1IsRecovered=-99.;
  e2IsDead=-99.;
  e2IsRecovered=-99.;


  edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
  


  std::vector<PileupSummaryInfo>::const_iterator PVI;
  nTruePU=1;
  if (isMC==true){
    iEvent.getByToken(pileupSummaryToken_, PupInfo);

    //float Tnpv = -1;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      
      int BX = PVI->getBunchCrossing();
      
      if(BX == 0) { 
	nTruePU = PVI->getTrueNumInteractions();
	continue;
      }
      
    }
  }
  MyWeight=1;
  // if (isMC==true) MyWeight = LumiWeights_.weight( nTruePU );
  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ <<  " " << MyWeight << std::endl;

  // edm::ESHandle<EcalChannelStatus> chStatus;
  // iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  // const EcalChannelStatus* ecs = nullptr;
  // if( chStatus.isValid() ) ecs = chStatus.product();
  // Get the beam spot
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(beamSpotToken_,theBeamSpot);  
  
  // Get electrons
  edm::Handle<reco::GsfElectronCollection > electrons;
  iEvent.getByToken(electronsToken_, electrons);
  // Get electrons
   // edm::Handle<pat::ElectronCollection > patelectrons;
   // iEvent.getByToken(patElectronsToken_, patelectrons);
   
  // Get MC collection
  Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_,genParticles);

  // Get PV
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(vtxToken_, vertices);

  edm::Handle<EcalRecHitCollection> rechit_EB_col;
  iEvent.getByToken(rechits_EB_,rechit_EB_col);
  // // get the geometry
  // edm::ESHandle<CaloSubdetectorGeometry> theBarrelGeometry_handle;
  // iSetup.get<EcalBarrelGeometryRecord>().get("EcalBarrel",theBarrelGeometry_handle);
  // const CaloSubdetectorGeometry *theBarrelGeometry;
  // theBarrelGeometry = &(*theBarrelGeometry_handle);
  
  // get the topology
  edm::ESHandle<CaloTopology> theCaloTopo;
  iSetup.get<CaloTopologyRecord>().get(theCaloTopo);
  const CaloTopology *topology = theCaloTopo.product();


  hEventCount->Fill(1);

  if (vertices->empty()) return; // skip the event if no PV found
  

  // Find the first vertex in the collection that passes
  VertexCollection::const_iterator firstGoodVertex = vertices->end();
  int firstGoodVertexIdx = 0;
  for (VertexCollection::const_iterator vtx = vertices->begin(); 
       vtx != vertices->end(); ++vtx, ++firstGoodVertexIdx) {
    bool isFake = vtx->isFake();
    // Check the goodness
    if ( !isFake
	 &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0
	 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==vertices->end() ) return; // skip event if there are no good PVs


  iEvent.getByToken(theRhoToken,_rhoHandle);
  

  // // Get the conversions collection
  // edm::Handle<reco::ConversionCollection> conversions;
  // iEvent.getByToken(conversionsToken_, conversions);

  int i=0; 
  vIeta_.clear();

  for (reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++){
    vIphi_.clear();
    vIsm_.clear();
    vIc_.clear();
    vRawId_.clear();
    vXtalEn_.clear();
    vEnFr_.clear();
    vkGood_.clear();
    vkPoorReco_.clear();
    vkOutOfTimE_.clear();
    vkFaultyHardware_.clear();
    vkNoisy_.clear();
    vkPoorCalib_.clear();
    vkSaturated_.clear();
    vkLeadingEdgeRecovered_.clear();
    vkNeighboursRecovered_.clear();
    vkTowerRecovered_.clear();
    vkDead_.clear();
    vkKilled_.clear();
    vkTPSaturated_.clear();
    vkL1SpikeFlag_.clear();
    vkWeird_.clear();
    vkDiWeird_.clear();
    vkHasSwitchToGain6_.clear();
    vkHasSwitchToGain1_.clear();
    vkUnknown_.clear();

    vNeigh0_.clear();
    vNeigh1_.clear();
    vNeigh2_.clear();
    vNeigh3_.clear();
    vNeigh5_.clear();
    vNeigh6_.clear();
    vNeigh7_.clear();
    vNeigh8_.clear();

    vSum8_.clear();
    vAve_.clear();
    std::vector<DetId> matrix;

    std::cout << __LINE__<< " vXtalEn_, first ele iterator " << vXtalEn_.size() << std::endl;
    isRecovered=false;
    isDead = false;

    //std::cout << "cut based id " << it->userInt("cutbasedID_veto")<< std::endl;
    //h_eID ->Fill (it->userInt("cutbasedID_veto"));  
    //  }

    //  // Loop over electrons
    //  for (size_t i = 0; i < electrons->size(); ++i){
    //    const auto el = electrons->ptrAt(i);
    const auto el=it; 
    // acceptance
    if( el->pt() < 1 ) continue;
    if( fabs(el->eta()) > 2.5 ) continue;
    double pt_e=el->pt();
    auto etrack  = el->gsfTrack(); 
    float abseta = fabs((el->superCluster().get())->position().eta());
    if (abseta<1.479){
      if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>2) continue;
      if (el->full5x5_sigmaIetaIeta() >  0.0128) continue;
      if (fabs(el->deltaPhiSuperClusterTrackAtVtx())> 0.159) continue;
      if (fabs(el->deltaEtaSeedClusterTrackAtVtx())> 0.00523) continue;
      if (el->hadronicOverEm()>0.247) continue;
      const float  eA = getEffectiveArea( abseta );
      const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle.product()) : 0; 
      if (( el->pfIsolationVariables().sumChargedHadronPt +  std::max(float(0.0), el->pfIsolationVariables().sumNeutralHadronEt +  el->pfIsolationVariables().sumPhotonEt -  eA*rho )) > 0.168*pt_e) continue;
      const float ecal_energy_inverse = 1.0/el->ecalEnergy();
      const float eSCoverP = el->eSuperClusterOverP();
      if ( std::abs(1.0 - eSCoverP)*ecal_energy_inverse > 0.193 )continue;
    }else{
      if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>3) continue;
      if (el->full5x5_sigmaIetaIeta() >  0.0445) continue;
      if (fabs(el->deltaPhiSuperClusterTrackAtVtx())>0.157 ) continue;
      if (fabs(el->deltaEtaSeedClusterTrackAtVtx())> 0.00984) continue;
      if (el->hadronicOverEm()>0.0982) continue;
      const float  eA = getEffectiveArea( abseta );
      const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle.product()) : 0; 
      if (( el->pfIsolationVariables().sumChargedHadronPt + std::max(float(0.0), el->pfIsolationVariables().sumNeutralHadronEt + el->pfIsolationVariables().sumPhotonEt -  eA*rho )) > 0.185*pt_e) continue;
      const float ecal_energy_inverse = 1.0/el->ecalEnergy();
      const float eSCoverP = el->eSuperClusterOverP();
      if ( std::abs(1.0 - eSCoverP)*ecal_energy_inverse > 0.0962 )continue;
    }
    std::cout << "=============" << std::endl;
    std::cout << "event: "<<    runNumber<< " "<< lumiBlock << " " << eventNumber << std::endl;

    const std::vector< std::pair< DetId, float > > hitAndFr_v = el->superCluster()->hitsAndFractions() ;
    vector<EBDetId> vIc1;
    
    // std::cout << __LINE__ << std::endl;
    for (unsigned int ii=0; ii<hitAndFr_v.size(); ii++){
      EBDetId idCurrent= hitAndFr_v[ii].first ;
      float enFr = hitAndFr_v[ii].second ;
      edm::SortedCollection<EcalRecHit>::const_iterator hit = rechit_EB_col->find( idCurrent );
      if (hit ==  rechit_EB_col->end()) continue;
      
      //std::cout << __LINE__ << endl ; 
      //if ( bool(hit->checkFlag(EcalRecHit::kNeighboursRecovered))==true) {
      std::string myXtalString = std::to_string(runNumber)+":"+std::to_string(lumiBlock)+":"+std::to_string(eventNumber);
      map<std::string, std::vector<long int> >::iterator xtalIt = xtalMap.find(myXtalString);
      if (xtalIt==xtalMap.end()) continue; 
      std::vector<long int>::iterator itRawId;
      itRawId = find (xtalIt->second.begin(), xtalIt->second.end(),idCurrent.rawId());
      if (itRawId!=xtalIt->second.end()){
	//std::cout << __LINE__ <<  "xtal found "<< std::endl;


	float neigh0=0.;
	float neigh1=0.;
	float neigh2=0.;
	float neigh3=0.;
	float neigh5=0.;
	float neigh6=0.;
	float neigh7=0.;
	float neigh8=0.;
	float sum8=0.;
	
	matrix = EcalClusterTools::matrixDetId( topology, idCurrent, -1, 1, -1, 1 );
	
	for ( size_t iM = 0; iM < matrix.size(); ++iM ) {
	  edm::SortedCollection<EcalRecHit>::const_iterator hitNeigh = rechit_EB_col->find(matrix[iM] );
	  if (hitNeigh ==  rechit_EB_col->end()) continue;
	  if (hitNeigh == hit) continue;
	  
	  if (hitNeigh->checkFlag(EcalRecHit::kL1SpikeFlag)==true) std::cout<< "SPIKE!" << std::endl;
	  if (hitNeigh->checkFlag(EcalRecHit::kSaturated)==true)continue;
	  //std::cout << "neighbour energy :"<<  i << " " <<  hitNeigh->energy() << ", sum8 " << sum8 << ", ave "<< sum8/8. << std::endl;
	  sum8+=hitNeigh->energy();
	  if (iM==0)  neigh0=hitNeigh->energy();
	  if (iM==1)  neigh1=hitNeigh->energy();
	  if (iM==2)  neigh2=hitNeigh->energy();
	  if (iM==3)  neigh3=hitNeigh->energy();
	  if (iM==5)  neigh5=hitNeigh->energy();
	  if (iM==6)  neigh6=hitNeigh->energy();
	  if (iM==7)  neigh7=hitNeigh->energy();
	  if (iM==8)  neigh8=hitNeigh->energy();
	}
	
	vSum8_.push_back(sum8);
	vAve_.push_back(sum8/8.);
	vNeigh0_.push_back(neigh0);
	vNeigh1_.push_back(neigh1);
	vNeigh2_.push_back(neigh2);
	vNeigh3_.push_back(neigh3);
	vNeigh5_.push_back(neigh5);
	vNeigh6_.push_back(neigh6);
	vNeigh7_.push_back(neigh7);
	vNeigh8_.push_back(neigh8);
	
	vkGood_.push_back(bool(hit->checkFlag(EcalRecHit::kGood)));
	vkPoorReco_.push_back(hit->checkFlag(EcalRecHit::kPoorReco));
	vkOutOfTimE_.push_back(hit->checkFlag(EcalRecHit::kOutOfTime));   
	vkFaultyHardware_.push_back(hit->checkFlag(EcalRecHit::kFaultyHardware));
	vkNoisy_.push_back(hit->checkFlag(EcalRecHit::kNoisy));
	vkPoorCalib_.push_back(hit->checkFlag(EcalRecHit::kPoorCalib));
	vkSaturated_.push_back(hit->checkFlag(EcalRecHit::kSaturated));
	vkLeadingEdgeRecovered_.push_back(hit->checkFlag(EcalRecHit::kLeadingEdgeRecovered));
	vkNeighboursRecovered_.push_back(hit->checkFlag(EcalRecHit::kNeighboursRecovered));
	vkTowerRecovered_.push_back(hit->checkFlag(EcalRecHit::kTowerRecovered));
	vkDead_.push_back(hit->checkFlag(EcalRecHit::kDead));    
	vkKilled_.push_back(hit->checkFlag(EcalRecHit::kKilled));
	vkTPSaturated_.push_back(hit->checkFlag(EcalRecHit::kTPSaturated));
	vkL1SpikeFlag_.push_back(hit->checkFlag(EcalRecHit::kL1SpikeFlag));
	vkWeird_.push_back(hit->checkFlag(EcalRecHit::kWeird));    
	vkDiWeird_.push_back(hit->checkFlag(EcalRecHit::kDiWeird));
	vkHasSwitchToGain6_.push_back(hit->checkFlag(EcalRecHit::kHasSwitchToGain6));    
	vkHasSwitchToGain1_.push_back(hit->checkFlag(EcalRecHit::kHasSwitchToGain1));
	vkUnknown_.push_back(hit->checkFlag(EcalRecHit::kUnknown));

	vIeta_.push_back(idCurrent.ieta());
	vIphi_.push_back(idCurrent.iphi());
	vIsm_.push_back(idCurrent.ism());
	vIc_.push_back(idCurrent.ic());
	vRawId_.push_back(idCurrent.rawId());
	vXtalEn_.push_back(hit->energy());
	vEnFr_.push_back(enFr); 
	isRecovered=true;
	//break;
      }//checking the hit is recovered
    }//running on the hit in the supercluster of the electron
    for (reco::GsfElectronCollection::const_iterator jt = it+1; jt != electrons->end(); jt++){
      isRecovered2=false;
      isDead2=false;

      if (jt->pt()< 1 ) continue;
      if( fabs(jt->eta()) > 2.5 ) continue;
      double pt_e=jt->pt();
      auto etrack  = jt->gsfTrack(); 
      float abseta = fabs((jt->superCluster().get())->position().eta());
      if (abseta<1.479){
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>2) continue;
	if (jt->full5x5_sigmaIetaIeta() >  0.0128) continue;
	if (fabs(jt->deltaPhiSuperClusterTrackAtVtx())> 0.159) continue;
	if (fabs(jt->deltaEtaSeedClusterTrackAtVtx())> 0.00523) continue;
	if (jt->hadronicOverEm()>0.247) continue;
	const float  eA = getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle.product()) : 0; 
	if (( jt->pfIsolationVariables().sumChargedHadronPt + 
	      std::max(float(0.0), jt->pfIsolationVariables().sumNeutralHadronEt +  
		       jt->pfIsolationVariables().sumPhotonEt -  eA*rho )
	      ) > 0.168*pt_e) continue;
	const float ecal_energy_inverse = 1.0/jt->ecalEnergy();
	const float eSCoverP = jt->eSuperClusterOverP();
	if ( std::abs(1.0 - eSCoverP)*ecal_energy_inverse > 0.193 )continue;
      }else{
	if (etrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)>3) continue;
	if (jt->full5x5_sigmaIetaIeta() >  0.0445) continue;
	if (fabs(jt->deltaPhiSuperClusterTrackAtVtx())>0.157 ) continue;
	if (fabs(jt->deltaEtaSeedClusterTrackAtVtx())> 0.00984) continue;
	if (jt->hadronicOverEm()>0.0982) continue;
	const float  eA = getEffectiveArea( abseta );
	const float rho = _rhoHandle.isValid() ? (float)(*_rhoHandle.product()) : 0; 
	if (( jt->pfIsolationVariables().sumChargedHadronPt + 
	      std::max(float(0.0), jt->pfIsolationVariables().sumNeutralHadronEt +  
		       jt->pfIsolationVariables().sumPhotonEt -  eA*rho )
	      ) > 0.185*pt_e) continue;
	const float ecal_energy_inverse = 1.0/jt->ecalEnergy();
	const float eSCoverP = jt->eSuperClusterOverP();
	if ( std::abs(1.0 - eSCoverP)*ecal_energy_inverse > 0.0962 )continue;

      }
      const std::vector< std::pair< DetId, float > > hitAndFr_v = jt->superCluster()->hitsAndFractions() ;
      vector<EBDetId> vIc2;
      // std::cout << __LINE__ << std::endl;
      for (unsigned int ii=0; ii<hitAndFr_v.size(); ii++){
	EBDetId idCurrent= hitAndFr_v[ii].first ;
	float enFr = hitAndFr_v[ii].second;
	edm::SortedCollection<EcalRecHit>::const_iterator hit = rechit_EB_col->find( idCurrent );
	if (hit ==  rechit_EB_col->end()) continue;

	std::string myXtalString = std::to_string(runNumber)+":"+std::to_string(lumiBlock)+":"+std::to_string(eventNumber);
	map<std::string, std::vector<long int> >::iterator xtalIt = xtalMap.find(myXtalString);
	//std::cout << __LINE__ << " map size " << xtalMap.size() << std::endl ; 
	if (xtalIt==xtalMap.end()) continue; 
	std::vector<long int>::iterator itRawId;
	itRawId = find (xtalIt->second.begin(), xtalIt->second.end(),idCurrent.rawId());
	if (itRawId!=xtalIt->second.end()){
	  //std::cout << __LINE__ <<  " xtal found "<< std::endl;

	  float neigh0=0.;
	  float neigh1=0.;
	  float neigh2=0.;
	  float neigh3=0.;
	  float neigh5=0.;
	  float neigh6=0.;
	  float neigh7=0.;
	  float neigh8=0.;
	  float sum8=0.;
	  
	  matrix = EcalClusterTools::matrixDetId( topology, idCurrent, -1, 1, -1, 1 );
	  
	  for ( size_t iM = 0; iM < matrix.size(); ++iM ) {
	    edm::SortedCollection<EcalRecHit>::const_iterator hitNeigh = rechit_EB_col->find(matrix[iM] );
	    if (hitNeigh ==  rechit_EB_col->end()) continue;
	    if (hitNeigh == hit) continue;
	  
	    if (hitNeigh->checkFlag(EcalRecHit::kL1SpikeFlag)==true) std::cout<< "SPIKE!" << std::endl;
	    if (hitNeigh->checkFlag(EcalRecHit::kSaturated)==true)continue;
	    sum8+=hitNeigh->energy();
	    //std::cout << "neighbour energy :"<<  i << " " <<  hitNeigh->energy() << ", sum8 " << sum8 << ", ave "<< sum8/8. << std::endl;
	   
	    if (iM==0)  neigh0=hitNeigh->energy();
	    if (iM==1)  neigh1=hitNeigh->energy();
	    if (iM==2)  neigh2=hitNeigh->energy();
	    if (iM==3)  neigh3=hitNeigh->energy();
	    if (iM==5)  neigh5=hitNeigh->energy();
	    if (iM==6)  neigh6=hitNeigh->energy();
	    if (iM==7)  neigh7=hitNeigh->energy();
	    if (iM==8)  neigh8=hitNeigh->energy();
	  }
	  
	  vSum8_.push_back(sum8);
	  vAve_.push_back(sum8/8.);
	  vNeigh0_.push_back(neigh0);
	  vNeigh1_.push_back(neigh1);
	  vNeigh2_.push_back(neigh2);
	  vNeigh3_.push_back(neigh3);
	  vNeigh5_.push_back(neigh5);
	  vNeigh6_.push_back(neigh6);
	  vNeigh7_.push_back(neigh7);
	  vNeigh8_.push_back(neigh8);

	  vkGood_.push_back(bool(hit->checkFlag(EcalRecHit::kGood)));
	  vkPoorReco_.push_back(hit->checkFlag(EcalRecHit::kPoorReco));
	  vkOutOfTimE_.push_back(hit->checkFlag(EcalRecHit::kOutOfTime));   
	  vkFaultyHardware_.push_back(hit->checkFlag(EcalRecHit::kFaultyHardware));
	  vkNoisy_.push_back(hit->checkFlag(EcalRecHit::kNoisy));
	  vkPoorCalib_.push_back(hit->checkFlag(EcalRecHit::kPoorCalib));
	  vkSaturated_.push_back(hit->checkFlag(EcalRecHit::kSaturated));
	  vkLeadingEdgeRecovered_.push_back(hit->checkFlag(EcalRecHit::kLeadingEdgeRecovered));
	  vkNeighboursRecovered_.push_back(hit->checkFlag(EcalRecHit::kNeighboursRecovered));
	  vkTowerRecovered_.push_back(hit->checkFlag(EcalRecHit::kTowerRecovered));
	  vkDead_.push_back(hit->checkFlag(EcalRecHit::kDead));    
	  vkKilled_.push_back(hit->checkFlag(EcalRecHit::kKilled));
	  vkTPSaturated_.push_back(hit->checkFlag(EcalRecHit::kTPSaturated));
	  vkL1SpikeFlag_.push_back(hit->checkFlag(EcalRecHit::kL1SpikeFlag));
	  vkWeird_.push_back(hit->checkFlag(EcalRecHit::kWeird));    
	  vkDiWeird_.push_back(hit->checkFlag(EcalRecHit::kDiWeird));
	  vkHasSwitchToGain6_.push_back(hit->checkFlag(EcalRecHit::kHasSwitchToGain6));    
	  vkHasSwitchToGain1_.push_back(hit->checkFlag(EcalRecHit::kHasSwitchToGain1));
	  vkUnknown_.push_back(hit->checkFlag(EcalRecHit::kUnknown));
	  
	  vIeta_.push_back(idCurrent.ieta());
	  vIphi_.push_back(idCurrent.iphi());
	  vIsm_.push_back(idCurrent.ism());
	  vIc_.push_back(idCurrent.ic());
	  vRawId_.push_back(idCurrent.rawId());
	  vXtalEn_.push_back(hit->energy());
	  vEnFr_.push_back(enFr); 

	  isRecovered2=true;
	  //break;
	}//checking it is recovered
      }//running on the rechit from the electron supercluster
      if (isRecovered==false && isRecovered2==false) continue;
      double t1 = TMath::Exp(-it->eta());
      double t2 = TMath::Exp(-jt->eta());
      double t1q = t1 * t1;
      double t2q = t2 * t2;
      double angle = 1 - ( (1 - t1q) * (1 - t2q) + 4 * t1 * t2 * cos(it->phi() - jt->phi())) / ((1 + t1q) * (1 + t2q) );
      //raw Energy invariant mass 
      invMass_rawSC = sqrt(2 * it->superCluster()->rawEnergy() * jt->superCluster()->rawEnergy() * angle); 
      reco::Candidate::LorentzVector e1p4(0,0,0,0);
      reco::Candidate::LorentzVector e2p4(0,0,0,0);
      e1p4.SetPxPyPzE(it->gsfTrack()->momentum().x(),it->gsfTrack()->momentum().y(),it->gsfTrack()->momentum().z(), it->superCluster()->rawEnergy() );
      e2p4.SetPxPyPzE(jt->gsfTrack()->momentum().x(),jt->gsfTrack()->momentum().y(),jt->gsfTrack()->momentum().z(), jt->superCluster()->rawEnergy() );
      reco::Candidate::LorentzVector zp4 = it->p4() + jt->p4();
      //reco::Candidate::LorentzVector zRawp4=e1p4+e2p4;
      if (abs( zp4.M() - 91)<40){
       	hAll_ZeeMass->Fill( zp4.M(),MyWeight);
       	hAll_ZMassRawEn->Fill(invMass_rawSC,MyWeight);
      	if (isRecovered==true ||  isRecovered2==true) {
      	  hRecovered_ZeeMass->Fill( zp4.M(),MyWeight);
       	  hRecoverd_ZMassRawEn->Fill(invMass_rawSC,MyWeight);
      	}else {
      	  h_ZeeMass->Fill( zp4.M(),MyWeight);
       	  h_ZMassRawEn->Fill(invMass_rawSC,MyWeight);
	}
	std::cout << "raw Mass " << invMass_rawSC << " e1 raw En " <<  it->superCluster()->rawEnergy() << " " << jt->superCluster()->rawEnergy() << std::endl;
 	mcGenWeight=-1.;
       	e1Charge=(int)it->charge();
	e2Charge=(int)jt->charge();
	std::cout << "e1 charge "<< e1Charge << " "<< it->charge() << ", e2 charge "<< e2Charge << " "<< jt->charge() << std::endl;
  	e1Eta=it->eta();
	e2Eta=jt->eta();
  	e1Phi=it->phi();
	e2Phi=jt->phi();
	e1R9=it->r9();
	e2R9=jt->r9();
	e1SigmaIetaIeta=it->full5x5_sigmaIetaIeta();
	e2SigmaIetaIeta=jt->full5x5_sigmaIetaIeta();
	  
	e1EtaSC=it->superCluster()->eta();
	e2EtaSC=jt->superCluster()->eta();
	e1PhiSC=it->superCluster()->phi();
	e2PhiSC=jt->superCluster()->phi();

	//if (it->superCluster()->caloID().detector(reco::CaloID::DET_ECAL_BARREL)==true){
	if(it->superCluster()->seed()->seed().subdetId()==EcalBarrel){
	  EBDetId e1SeedDetId(it->superCluster()->seed()->seed());
	  std::cout << "Seed " <<  e1SeedDetId << " "<< e1SeedDetId.ieta()<< " " << e1SeedDetId.iphi()<< std::endl;
	  e1SeedIEta=(float)e1SeedDetId.ieta();
	  e1SeedIPhi=(float)e1SeedDetId.iphi();
	}
	//if (jt->superCluster()->caloID().detector(reco::CaloID::DET_ECAL_BARREL)==true){
	if(jt->superCluster()->seed()->seed().subdetId()==EcalBarrel){
	  EBDetId e2SeedDetId(jt->superCluster()->seed()->seed());
	  std::cout << "Seed " <<  e2SeedDetId << " "<< e2SeedDetId.ieta()<< " " << e2SeedDetId.iphi()<< std::endl;
	  e2SeedIEta=(float)e2SeedDetId.ieta();
	  e2SeedIPhi=(float)e2SeedDetId.iphi();
	}
	std::cout << runNumber << " " << lumiBlock << " " << eventNumber << " " << eventTime << " " << nBX << std::endl;

	std::cout <<e1SeedIEta << " " <<  e1SeedIPhi << " " << e2SeedIEta <<  " "<< e2SeedIPhi << std::endl;
  	e1Energy=it->energy();
	e1RawEnergy=(float)it->superCluster()->rawEnergy();
	//e1Energy3x3SC;
	//e1Energy5x5SC;
	e2Energy=jt->energy();
	e2RawEnergy=(float)jt->superCluster()->rawEnergy();
	//e2Energy3x3SC;
	//e2Energy5x5SC;
  	e1Pt=it->pt();
	e2Pt=jt->pt();
	invMass= zp4.M();
	
  	if (isMC==true){
	  std::pair<int,double> e1GenMatchEn=std::make_pair(0,0);
	  std::pair<int,double> e2GenMatchEn=std::make_pair(0,0);
	  e1GenMatchEn=matchToTruth( &*it, genParticles);
	  e2GenMatchEn=matchToTruth( &*jt, genParticles);
	  e1GenEnergy=e1GenMatchEn.second;
	  e2GenEnergy=e2GenMatchEn.second;
	}
	// e1IsDead=(int)isDead;
	e1IsRecovered=(int)isRecovered;
	// e2IsDead=(int) isDead2;
	e2IsRecovered=(int)isRecovered2;
	std::cout << "kGood size " << vkGood_.size() << std::endl;
	vkGood=		 vkGood_;
	vkPoorReco=		 vkPoorReco_;
	vkOutOfTimE=		 vkOutOfTimE_;
	vkFaultyHardware=	 vkFaultyHardware_;
	vkNoisy=		 vkNoisy_;
	vkPoorCalib=		 vkPoorCalib_;
	vkSaturated=		 vkSaturated_;
	vkLeadingEdgeRecovered= vkLeadingEdgeRecovered_;
	vkNeighboursRecovered= vkNeighboursRecovered_;
	vkTowerRecovered=	 vkTowerRecovered_;
	vkDead=		 vkDead_;
	vkKilled=		 vkKilled_;
	vkTPSaturated=		 vkTPSaturated_;
	vkL1SpikeFlag=		 vkL1SpikeFlag_;
	vkWeird=		 vkWeird_;
	vkDiWeird=		 vkDiWeird_;
	vkHasSwitchToGain6=	 vkHasSwitchToGain6_;
	vkHasSwitchToGain1=	 vkHasSwitchToGain1_;
	vkUnknown=              vkUnknown_;              

	vSum8=vSum8_;
	vAve=vAve_;

	vIeta=vIeta_;
	vIphi=vIphi_;
	vIsm=vIsm_;
	vIc=vIc_;
	vRawId=vRawId_;
	vXtalEn=vXtalEn_;
	vEnFr=vEnFr_;

	
	vNeigh0=vNeigh0_;
	vNeigh1=vNeigh1_;
	vNeigh2=vNeigh2_;
	vNeigh3=vNeigh3_;
	vNeigh5=vNeigh5_;
	vNeigh6=vNeigh6_;
	vNeigh7=vNeigh7_;
	vNeigh8=vNeigh8_;



	std::cout << __LINE__ << " " << vXtalEn_.size() << std::endl;
	//filling the tree
	tree->Fill();
      }      
    }

    // MC truth match
    //    if (isMC==true) h_isTrue -> Fill(matchToTruth( el, genParticles));
    int eMatch=0;
    double eGenEn=0;
    std::pair<int,double> eGenMatchEn=std::make_pair(0,0);

    if (isMC==true){
      eGenMatchEn=matchToTruth( &*it, genParticles);
      eMatch=eGenMatchEn.first;
      eGenEn=eGenMatchEn.second;
      //if (isRecovered==true )std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " gen match  " << eMatch<< " " << eGenEn << std::endl;
    }
      
    // Kinematics and shower shapes
    float scEta = el->superCluster()->eta();
    float scPhi = el->superCluster()->phi();
    hAll_eta -> Fill(scEta,MyWeight);
    hAll_phi -> Fill(scPhi,MyWeight);
    if (fabs(scEta)<1.5) {
      hAll_EB_pt -> Fill(el->pt(),MyWeight);
      hAll_EB_rawEne -> Fill( el->superCluster()->rawEnergy() ,MyWeight);
      if (isMC==true && eMatch>0 && eGenEn!=0.) {
	hAll_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn,MyWeight );
	hAll_EB_eGenEn->Fill(eGenEn,MyWeight);
	hAll_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn ,MyWeight);
      }
      hAll_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight); 
      hAll_EB_r9 -> Fill( el->r9() ,MyWeight);  
      hAll_EB_r9uz -> Fill( el->r9() ,MyWeight);  
      hAll_EB_hoe -> Fill( el->full5x5_hcalOverEcal() ,MyWeight);  
    } else {
      hAll_EE_pt -> Fill(el->pt(),MyWeight);
      hAll_EE_rawEne -> Fill( el->superCluster()->rawEnergy() ,MyWeight);
      hAll_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight); 
      hAll_EE_r9 -> Fill( el->r9(),MyWeight );  
      hAll_EE_r9uz -> Fill( el->r9() ,MyWeight);  
      hAll_EE_hoe -> Fill( el->full5x5_hcalOverEcal() ,MyWeight);  
    }

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    if (fabs(scEta)<1.5)
      hAll_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) ,MyWeight);
    else
      hAll_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) ,MyWeight);
    

    if (isRecovered==true){
      hRecovered_eta -> Fill(scEta,MyWeight);
      hRecovered_phi -> Fill(scPhi,MyWeight);
      if (fabs(scEta)<1.5) {
	hRecovered_EB_pt -> Fill(el->pt(),MyWeight);
	hRecovered_EB_rawEne -> Fill( el->superCluster()->rawEnergy() ,MyWeight);
	hRecovered_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight); 
	hRecovered_EB_r9 -> Fill( el->r9() ,MyWeight);  
	hRecovered_EB_r9uz -> Fill( el->r9() ,MyWeight);  
	hRecovered_EB_hoe -> Fill( el->full5x5_hcalOverEcal(),MyWeight );  
	if (isMC==true && eMatch>0 && eGenEn!=0.) {
	  hRecovered_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn ,MyWeight);
	  hRecovered_EB_eGenEn->Fill(eGenEn,MyWeight);
	  hRecovered_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn ,MyWeight );
	}

      } else {
	hRecovered_EE_pt -> Fill(el->pt(),MyWeight);
	hRecovered_EE_rawEne -> Fill( el->superCluster()->rawEnergy() ,MyWeight);
	hRecovered_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight); 
	hRecovered_EE_r9 -> Fill( el->r9(),MyWeight );  
	hRecovered_EE_r9uz -> Fill( el->r9(),MyWeight );  
	hRecovered_EE_hoe -> Fill( el->full5x5_hcalOverEcal(),MyWeight );  
      }
      
      if (fabs(scEta)<1.5)
	hRecovered_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ),MyWeight );
      else
	hRecovered_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ),MyWeight );
    }else{
      h_eta -> Fill(scEta);
      h_phi -> Fill(scPhi);
      if (fabs(scEta)<1.5) {
	h_EB_pt -> Fill(el->pt(),MyWeight);
	h_EB_rawEne -> Fill( el->superCluster()->rawEnergy(),MyWeight );
	h_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight ); 
	h_EB_r9 -> Fill( el->r9() ,MyWeight);  
	h_EB_r9uz -> Fill( el->r9(),MyWeight );  
	h_EB_hoe -> Fill( el->full5x5_hcalOverEcal() ,MyWeight);  
	if (isMC==true && eMatch>0 && eGenEn!=0.) {
	  h_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn ,MyWeight);
	  h_EB_eGenEn->Fill(eGenEn,MyWeight);
	  h_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn ,MyWeight);
	}

      } else {
	h_EE_pt -> Fill(el->pt(),MyWeight);
	h_EE_rawEne -> Fill( el->superCluster()->rawEnergy() ,MyWeight);
	h_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ,MyWeight); 
	h_EE_r9 -> Fill( el->r9() ,MyWeight);  
	h_EE_r9uz -> Fill( el->r9(),MyWeight );  
	h_EE_hoe -> Fill( el->full5x5_hcalOverEcal(),MyWeight );  
      }

      // Impact parameter
      reco::GsfTrackRef theTrack = el->gsfTrack();
      if (fabs(scEta)<1.5)
	h_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ),MyWeight );
      else
	h_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ),MyWeight );
      
    }

    // // Conversion rejection
    // bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, theBeamSpot->position());
    // if (fabs(scEta)<1.5)   
    //   h_EB_conv -> Fill( (int) passConvVeto ); 
    // else
    //   h_EE_conv -> Fill( (int) passConvVeto ); 
    i++;
  } // Loop over electrons
   

}


void ElectronTreeXtalList::beginJob() { 
  std::string intXtal; 
  std::ifstream infile("xtals_ZeePierreTag_sum8gt20.txt");
  std::cout << __LINE__ << " file open? " <<infile.is_open()<< std::endl; 
  if (infile.is_open()){
    int lineCount=0;
    while(getline(infile, intXtal)){
      lineCount++;
      if (lineCount<11)std::cout << intXtal << std::endl;
      std::string mystring=intXtal.substr(0, intXtal.find_last_of(":"));
      long int myRawId=std::stol(intXtal.substr( intXtal.find_last_of(":")+1));
      if(xtalMap.find(mystring)!=xtalMap.end()){
	xtalMap[mystring].push_back(myRawId);
      }else{
	std::vector<long int> xtals;
	xtals.push_back(myRawId); 
	xtalMap.insert(std::pair<std::string, std::vector<long int> >(mystring, xtals));
      }
    }
    infile.close();
  }

}

void ElectronTreeXtalList::endJob() { }

std::pair<int,double> ElectronTreeXtalList::matchToTruth(const reco::GsfElectron* el, edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
  double genEn=0; 
  const reco::Candidate *closestElectron = 0;
  for(size_t i=0; i<prunedGenParticles->size();i++){
    const reco::Candidate *particle = &(*prunedGenParticles)[i];
    // Drop everything that is not electron or not status 1
    if( abs(particle->pdgId()) != 11 || particle->status() != 1 )
      continue;
    //
    double dRtmp = ROOT::Math::VectorUtil::DeltaR( el->p4(), particle->p4() );
    if( dRtmp < dR ){
      dR = dRtmp;
      closestElectron = particle;
      genEn=particle->energy();
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return std::make_pair(UNMATCHED,genEn);
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);
  //std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " << ancestorPID << " " << ancestorStatus << " " << genEn << std::endl;

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return std::make_pair(UNMATCHED,genEn);
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 62 )
    return std::make_pair(TRUE_NON_PROMPT_ELECTRON, genEn);

  if( abs(ancestorPID) == 23 && ancestorStatus == 62 )
    return std::make_pair(TRUE_ELECTRON_FROM_Z, genEn);

  if( abs(ancestorPID) == 22 && ancestorStatus == 62 )
    return std::make_pair(TRUE_ELECTRON_FROM_GAMMA, genEn);

  if( abs(ancestorPID) == 15 && ancestorStatus == 62 )
    return std::make_pair(TRUE_ELECTRON_FROM_TAU,genEn);

  // What remains is true prompt electrons
  return std::make_pair(TRUE_PROMPT_ELECTRON, genEn);
}

void ElectronTreeXtalList::findFirstNonElectronMother(const reco::Candidate *particle,
						 int &ancestorPID, int &ancestorStatus){

  if( particle == 0 ){
    printf("ElectronNtupler: ERROR! null candidate pointer, this should never happen\n");
    return;
  }

  // Is this the first non-electron parent? If yes, return, otherwise
  // go deeper into recursion
  if( abs(particle->pdgId()) == 11 ){
    findFirstNonElectronMother(particle->mother(0), ancestorPID, ancestorStatus);
  }else{
    ancestorPID = particle->pdgId();
    ancestorStatus = particle->status();
  }

  return;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronTreeXtalList);
