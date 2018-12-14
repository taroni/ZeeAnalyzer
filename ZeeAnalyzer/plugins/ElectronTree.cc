// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
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

#include "TTree.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"


#define NELE 3
#define initFloat     {-999.,-999.,-999.}
#define initInt       {0,0,0}
#define initIntCharge {-100,-100,-100}

class ElectronTree : public edm::EDAnalyzer {
public:
  explicit ElectronTree(const edm::ParameterSet&);
  ~ElectronTree();
  
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
  
  // ----------member data ---------------------------
  
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

  TFile *tree_file;
  TTree * tree; 

  UInt_t     	runNumber;		///< run number
  UShort_t      lumiBlock;		///< lumi section
  Long64_t    eventNumber;	///< event number
  UInt_t	eventTime;		///< unix time of the event
  UShort_t	      nBX;			///< bunch crossing
  
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


  TH1F *h_ZeeMass, *h_ZMassRawEn,*hAll_ZeeMass, *hAll_ZMassRawEn, *hRecovered_ZeeMass, *hRecoverd_ZMassRawEn, *hEventCount; 
};

ElectronTree::ElectronTree(const edm::ParameterSet& iConfig) {

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

  edm::Service<TFileService> fs;
 
	//now do what ever initialization is needed
  tree = fs->make<TTree>("selected", "selected");
	//tree = fs->make<TTree>("selected","selected"); //no otherwise you have the extraCalib in the same file


  tree->Branch("runNumber",     &runNumber,   "runNumber/i");
  tree->Branch("lumiBlock",     &lumiBlock,   "lumiBlock/s");
  tree->Branch("eventNumber",   &eventNumber, "eventNumber/l");
  tree->Branch("eventTime",     &eventTime,   "eventTime/i");
  tree->Branch("nBX",           &nBX,         "nBX/s");
  

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
}


ElectronTree::~ElectronTree() { }

void ElectronTree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;
  using namespace edm;
  using namespace reco;

  runNumber=iEvent.id().run();
  lumiBlock=(UShort_t)iEvent.id().luminosityBlock();
  eventNumber=(Long64_t) iEvent.id().event();
  edm::Timestamp time = iEvent.eventAuxiliary().time();
  eventTime= (UInt_t) time.unixTime();
  nBX=(UShort_t) iEvent.bunchCrossing();	

 
      
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

  


  // // Get the conversions collection
  // edm::Handle<reco::ConversionCollection> conversions;
  // iEvent.getByToken(conversionsToken_, conversions);

  int i=0; 
  for (reco::GsfElectronCollection::const_iterator it = electrons->begin(); it != electrons->end(); it++){
    bool isRecovered=false;
    bool isDead = false;
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

    const std::vector< std::pair< DetId, float > > hitAndFr_v = el->superCluster()->hitsAndFractions() ;
    for (unsigned int ii=0; ii<hitAndFr_v.size(); ii++){
      EBDetId idCurrent= hitAndFr_v[ii].first ;
      edm::SortedCollection<EcalRecHit>::const_iterator hit = rechit_EB_col->find( idCurrent );
      if ( hit->checkFlag(EcalRecHit::kNeighboursRecovered)) {
	  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<  idCurrent << " " << hit->checkFlag(EcalRecHit::kNeighboursRecovered) <<  " " << hit->energy()<< std::endl;

	  
	isRecovered=true;
	//break;
      }
      if( hit->checkFlag(EcalRecHit::kDead)) {
	isDead=true;
      }//kDead
    }
    bool isRecovered2=false;
    bool isDead2=false;
    for (reco::GsfElectronCollection::const_iterator jt = it+1; jt != electrons->end(); jt++){
      if (jt->pt()< 1 ) continue;
      if( fabs(jt->eta()) > 2.5 ) continue;
      const std::vector< std::pair< DetId, float > > hitAndFr_v = jt->superCluster()->hitsAndFractions() ;
      for (unsigned int ii=0; ii<hitAndFr_v.size(); ii++){
	EBDetId idCurrent= hitAndFr_v[ii].first ;
	edm::SortedCollection<EcalRecHit>::const_iterator hit = rechit_EB_col->find( idCurrent );
	if ( hit->checkFlag(EcalRecHit::kNeighboursRecovered)) {
	  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<  idCurrent << " " << hit->checkFlag(EcalRecHit::kNeighboursRecovered) <<  " " << hit->energy()<< std::endl;
	  isRecovered2=true;
	  //break;
	}
	if( hit->checkFlag(EcalRecHit::kDead)) {
	  isDead2=true;
        }//kDead
      }
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
       	hAll_ZeeMass->Fill( zp4.M());
       	hAll_ZMassRawEn->Fill(invMass_rawSC);
      	if (isRecovered==true ||  isRecovered2==true) {
      	  hRecovered_ZeeMass->Fill( zp4.M());
       	  hRecoverd_ZMassRawEn->Fill(invMass_rawSC);
      	}else {
      	  h_ZeeMass->Fill( zp4.M());
       	  h_ZMassRawEn->Fill(invMass_rawSC);
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
	e1IsDead=(int)isDead;
	e1IsRecovered=(int)isRecovered;
	e2IsDead=(int) isDead2;
	e2IsRecovered=(int)isRecovered2;
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
    hAll_eta -> Fill(scEta);
    hAll_phi -> Fill(scPhi);
    if (fabs(scEta)<1.5) {
      hAll_EB_pt -> Fill(el->pt());
      hAll_EB_rawEne -> Fill( el->superCluster()->rawEnergy() );
      if (isMC==true && eMatch>0 && eGenEn!=0.) {
	hAll_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn );
	hAll_EB_eGenEn->Fill(eGenEn);
	hAll_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn );
      }
      hAll_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
      hAll_EB_r9 -> Fill( el->r9() );  
      hAll_EB_r9uz -> Fill( el->r9() );  
      hAll_EB_hoe -> Fill( el->full5x5_hcalOverEcal() );  
    } else {
      hAll_EE_pt -> Fill(el->pt());
      hAll_EE_rawEne -> Fill( el->superCluster()->rawEnergy() );
      hAll_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
      hAll_EE_r9 -> Fill( el->r9() );  
      hAll_EE_r9uz -> Fill( el->r9() );  
      hAll_EE_hoe -> Fill( el->full5x5_hcalOverEcal() );  
    }

    // Impact parameter
    reco::GsfTrackRef theTrack = el->gsfTrack();
    if (fabs(scEta)<1.5)
      hAll_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
    else
      hAll_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
    

    if (isRecovered==true){
      hRecovered_eta -> Fill(scEta);
      hRecovered_phi -> Fill(scPhi);
      if (fabs(scEta)<1.5) {
	hRecovered_EB_pt -> Fill(el->pt());
	hRecovered_EB_rawEne -> Fill( el->superCluster()->rawEnergy() );
	hRecovered_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
	hRecovered_EB_r9 -> Fill( el->r9() );  
	hRecovered_EB_r9uz -> Fill( el->r9() );  
	hRecovered_EB_hoe -> Fill( el->full5x5_hcalOverEcal() );  
	if (isMC==true && eMatch>0 && eGenEn!=0.) {
	  hRecovered_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn );
	  hRecovered_EB_eGenEn->Fill(eGenEn);
	  hRecovered_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn );
	}

      } else {
	hRecovered_EE_pt -> Fill(el->pt());
	hRecovered_EE_rawEne -> Fill( el->superCluster()->rawEnergy() );
	hRecovered_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
	hRecovered_EE_r9 -> Fill( el->r9() );  
	hRecovered_EE_r9uz -> Fill( el->r9() );  
	hRecovered_EE_hoe -> Fill( el->full5x5_hcalOverEcal() );  
      }
      
      if (fabs(scEta)<1.5)
	hRecovered_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
      else
	hRecovered_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
    }else{
      h_eta -> Fill(scEta);
      h_phi -> Fill(scPhi);
      if (fabs(scEta)<1.5) {
	h_EB_pt -> Fill(el->pt());
	h_EB_rawEne -> Fill( el->superCluster()->rawEnergy() );
	h_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
	h_EB_r9 -> Fill( el->r9() );  
	h_EB_r9uz -> Fill( el->r9() );  
	h_EB_hoe -> Fill( el->full5x5_hcalOverEcal() );  
	if (isMC==true && eMatch>0 && eGenEn!=0.) {
	  h_EB_rawGenEnRatio -> Fill( el->superCluster()->rawEnergy()/eGenEn );
	  h_EB_eGenEn->Fill(eGenEn);
	  h_EB_eEnGenEnRatio-> Fill( el->energy()/eGenEn );
	}

      } else {
	h_EE_pt -> Fill(el->pt());
	h_EE_rawEne -> Fill( el->superCluster()->rawEnergy() );
	h_EE_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
	h_EE_r9 -> Fill( el->r9() );  
	h_EE_r9uz -> Fill( el->r9() );  
	h_EE_hoe -> Fill( el->full5x5_hcalOverEcal() );  
      }

      // Impact parameter
      reco::GsfTrackRef theTrack = el->gsfTrack();
      if (fabs(scEta)<1.5)
	h_EB_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
      else
	h_EE_dz -> Fill( theTrack->dz( firstGoodVertex->position() ) );
      
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


void ElectronTree::beginJob() { }

void ElectronTree::endJob() { }

std::pair<int,double> ElectronTree::matchToTruth(const reco::GsfElectron* el, edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

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

void ElectronTree::findFirstNonElectronMother(const reco::Candidate *particle,
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
DEFINE_FWK_MODULE(ElectronTree);
