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

class ElectronPlots : public edm::EDAnalyzer {
public:
  explicit ElectronPlots(const edm::ParameterSet&);
  ~ElectronPlots();
  
  enum ElectronMatchType {UNMATCHED = 0, 
			  TRUE_PROMPT_ELECTRON, 
			  TRUE_ELECTRON_FROM_TAU,
			  TRUE_NON_PROMPT_ELECTRON}; // The last does not include tau parents
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  int matchToTruth(const edm::Ptr<reco::GsfElectron> el, 
		   const edm::Handle<edm::View<reco::GenParticle>>  &genParticles);
  
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
  TH1F *hRecovered_eta, *hRecovered_phi, *hRecovered_isTrue, *hRecovered_eID;
  TH1F *hRecovered_EB_pt, *hRecovered_EE_pt;
  TH1F *hRecovered_EB_rawEne, *hRecovered_EE_rawEne;
  TH1F *hRecovered_EB_sigmaIeIe, *hRecovered_EE_sigmaIeIe;
  TH1F *hRecovered_EB_r9, *hRecovered_EE_r9;
  TH1F *hRecovered_EB_r9uz, *hRecovered_EE_r9uz;
  TH1F *hRecovered_EB_hoe, *hRecovered_EE_hoe;
  TH1F *hRecovered_EB_dz, *hRecovered_EE_dz;
  TH1F *hRecovered_EB_conv, *hRecovered_EE_conv;

  TH1F *hZeeMass, *hZMassRawEn, *hRecovered_ZeeMass, *hRecoverd_ZMassRawEn, *hEventCount; 
};

ElectronPlots::ElectronPlots(const edm::ParameterSet& iConfig) {

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

  hZeeMass=fs->make<TH1F>("hZeeMass", "Zee Inv Mass", 100, 0, 200);
  hZMassRawEn=fs->make<TH1F>("hZMassRawEn", "Zee Inv Mass, rawEn", 100, 0, 200);
  hRecovered_ZeeMass=fs->make<TH1F>("hRecovered_ZeeMass", "Zee Inv Mass, recov", 100, 0, 200);
  hRecoverd_ZMassRawEn=fs->make<TH1F>("hRecoverd_ZMassRawEn", "Zee Inv Mass, recov, rawEn", 100, 0, 200);

  typedef reco::Candidate::LorentzVector LorentzVector;
}


ElectronPlots::~ElectronPlots() { }

void ElectronPlots::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace std;
  using namespace edm;
  using namespace reco;

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
	std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<  idCurrent << " " << hit->checkFlag(EcalRecHit::kNeighboursRecovered) << std::endl;
	isRecovered=true;
	break;
      }

    }
     

    bool isRecovered2=false;
    for (reco::GsfElectronCollection::const_iterator jt = it+1; jt != electrons->end(); jt++){
      if (jt->pt()< 1 ) continue;
      if( fabs(jt->eta()) > 2.5 ) continue;
      const std::vector< std::pair< DetId, float > > hitAndFr_v = jt->superCluster()->hitsAndFractions() ;
      for (unsigned int ii=0; ii<hitAndFr_v.size(); ii++){
	EBDetId idCurrent= hitAndFr_v[ii].first ;
	edm::SortedCollection<EcalRecHit>::const_iterator hit = rechit_EB_col->find( idCurrent );
	if ( hit->checkFlag(EcalRecHit::kNeighboursRecovered)) {
	  std::cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " " <<  idCurrent << " " << hit->checkFlag(EcalRecHit::kNeighboursRecovered) << std::endl;
	  isRecovered2=true;
	  break;
	}
	  
      }
      
      //pat::CompositeCandidate diele;
      //diele.addDaughter(*it);
      //diele.addDaughter(*jt);
      reco::Candidate::LorentzVector e1p4(0,0,0,0);
      reco::Candidate::LorentzVector e2p4(0,0,0,0);
      
      // if (isRecovered==true){
	e1p4.SetPxPyPzE(it->gsfTrack()->momentum().x(),it->gsfTrack()->momentum().y(),it->gsfTrack()->momentum().z(), it->superCluster()->rawEnergy() );
      // } else {
      //   e1p4= it->p4();
      // }
      // if (isRecovered2==true){
	e2p4.SetPxPyPzE(jt->gsfTrack()->momentum().x(),jt->gsfTrack()->momentum().y(),jt->gsfTrack()->momentum().z(), jt->superCluster()->rawEnergy() );
      // } else {
      // 	e2p4= jt->p4();
      // }
      reco::Candidate::LorentzVector zp4 = it->p4() + jt->p4();
      reco::Candidate::LorentzVector zRawp4=e1p4+e2p4;
      //diele.setP4(zp4);
      if (abs( zp4.M() - 91)<20){
	//std::cout << iEvent.id().event() << " electronsize "<< electrons->size()<<  ", mass  " << zp4.M() << ", electron "<< i<< ", pt "<< it->pt() << ", electron "<< i+1 << jt->pt() << std::endl;   
	hZeeMass->Fill( zp4.M());
	hZMassRawEn->Fill( zRawp4.M());
	if (isRecovered==true ||  isRecovered2==true) {
	  hRecovered_ZeeMass->Fill( zp4.M());
	  hRecoverd_ZMassRawEn->Fill( zRawp4.M());
	}

      }
      
    }

    // MC truth match
    //    if (isMC==true) h_isTrue -> Fill(matchToTruth( el, genParticles));

    // Kinematics and shower shapes
    float scEta = el->superCluster()->eta();
    float scPhi = el->superCluster()->phi();
    h_eta -> Fill(scEta);
    h_phi -> Fill(scPhi);
    if (fabs(scEta)<1.5) {
      h_EB_pt -> Fill(el->pt());
      h_EB_rawEne -> Fill( el->superCluster()->rawEnergy() );
      h_EB_sigmaIeIe -> Fill( el->full5x5_sigmaIetaIeta() ); 
      h_EB_r9 -> Fill( el->r9() );  
      h_EB_r9uz -> Fill( el->r9() );  
      h_EB_hoe -> Fill( el->full5x5_hcalOverEcal() );  
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


void ElectronPlots::beginJob() { }

void ElectronPlots::endJob() { }

int ElectronPlots::matchToTruth(const edm::Ptr<reco::GsfElectron> el, const edm::Handle<edm::View<reco::GenParticle>> &prunedGenParticles){

  // Find the closest status 1 gen electron to the reco electron
  double dR = 999;
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
    }
  }
  // See if the closest electron (if it exists) is close enough.
  // If not, no match found.
  if( !(closestElectron != 0 && dR < 0.1) ) {
    return UNMATCHED;
  }

  // 
  int ancestorPID = -999; 
  int ancestorStatus = -999;
  findFirstNonElectronMother(closestElectron, ancestorPID, ancestorStatus);

  if( ancestorPID == -999 && ancestorStatus == -999 ){
    // No non-electron parent??? This should never happen.
    // Complain.
    printf("ElectronNtupler: ERROR! Electron does not apper to have a non-electron parent\n");
    return UNMATCHED;
  }
  
  if( abs(ancestorPID) > 50 && ancestorStatus == 2 )
    return TRUE_NON_PROMPT_ELECTRON;

  if( abs(ancestorPID) == 15 && ancestorStatus == 2 )
    return TRUE_ELECTRON_FROM_TAU;

  // What remains is true prompt electrons
  return TRUE_PROMPT_ELECTRON;
}

void ElectronPlots::findFirstNonElectronMother(const reco::Candidate *particle,
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
DEFINE_FWK_MODULE(ElectronPlots);
