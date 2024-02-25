// -*- C++ -*-
//
// Package:    EventAnalyzer/Analysis
// Class:      Analysis
// 
/**\class Analysis Analysis.cc EventAnalyzer/Analysis/plugins/Analysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 15 Feb 2024 08:14:27 GMT
//
//


// system include files
#include <memory>
#include <string>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//
// class declaration
//
#include "EventAnalyzer/Analysis/interface/TreeWriter.h"
#include "EventAnalyzer/Analysis/src/TreeWriter.cc"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/BTauReco/interface/TaggingVariable.h"


// simulated vertices,..., add <use name=SimDataFormats/Vertex> and <../Track>
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/TrackIPTagInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TH1F.h"
#include "TTree.h"
#include "TString.h"
#include "TFile.h"


// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class Analysis : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit Analysis(const edm::ParameterSet&);
		~Analysis();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

// --------------member data ---------------------------

//		std::vector<std::string> triggerList_;
		edm::Service<TFileService> fs;

		edm::Handle<unsigned int> lumiBlock;
		edm::Handle<unsigned int> runNumber;
		edm::Handle<unsigned int> eventNumber;
				
// Gen Partile Info
//		edm::EDGetTokenT<std::vector<reco::GenParticle>> genPartToken_;
//		edm::Handle<std::vector<reco:;GenParticle>> genParticles;		

// PU info
		edm::EDGetTokenT<std::vector<PileupSummaryInfo>> PUsrc_;
		edm::EDGetTokenT<std::vector<reco::Vertex>> PVsrc_;

		edm::Handle<std::vector<PileupSummaryInfo>> PUInfo;
		edm::Handle<std::vector<reco::Vertex>> PVInfo; 
		

// gen jets
		edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsToken_;
		edm::EDGetTokenT<std::vector<reco::GenJet> > genJetsAK8Token_;

		edm::Handle<std::vector<reco::GenJet> > genJets;
		edm::Handle<std::vector<reco::GenJet> > genJetsAK8;

// reco jets token
		edm::EDGetTokenT<std::vector<pat::Jet> > pfJetsToken_;
		edm::EDGetTokenT<std::vector<pat::Jet> > puppiJetsToken_;
		edm::EDGetTokenT<std::vector<pat::Jet> > recoJetsAK8Token_;		

//		edm::EDGetTokenT<edm::OwnVector<reco::BaseTagInfo> > tagInfoToken_;

		edm::Handle<std::vector<pat::Jet> > pfJets;
		edm::Handle<std::vector<pat::Jet> > puppiJets;
		edm::Handle<std::vector<pat::Jet> > recoJetsAK8;
		
//		edm::Handle<edm::OwnVector<reco::BaseTagInfo> > tagInfo;

// PF Cands
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> pfCandsToken_;
		edm::Handle<std::vector<pat::PackedCandidate>> pfCands;

// MET
		edm::EDGetTokenT<std::vector<pat::MET>> pfmetToken_;
		edm::EDGetTokenT<std::vector<pat::MET>> puppimetToken_;
		edm::Handle<std::vector<pat::MET>> pfMET;
		edm::Handle<std::vector<pat::MET>> puppiMET;
		                      
		
		float dRthreshold = 0.2;
		unsigned int nPV, nTruePU;
		 
		bool isMC;
		double pf_met, puppi_met;

		double pt_cand, pt_jet;
		
		TTree* tree;
		TreeWriter t;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analysis::Analysis(const edm::ParameterSet& iConfig)

{
			
	isMC =iConfig.getUntrackedParameter<bool>("isMC",false);
	pt_cand =iConfig.getParameter<double>("Cand_Pt_THR");
	pt_jet =iConfig.getParameter<double>("Jet_Pt_THR");

// get PU info token
   PVsrc_ = consumes<std::vector<reco::Vertex>>(edm::InputTag("offlineSlimmedPrimaryVertices"));
   
// get PF cand token
	pfCandsToken_ = consumes<std::vector<pat::PackedCandidate>>(edm::InputTag("packedPFCandidates"));

	if(isMC){
		PUsrc_ = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
//  gen particle token		
//		genPartToken_ = consumes<std::vector<reco::GenParticle>>((edm::InputTag("prunedGenParticles"));
		
// gen jet token
		genJetsToken_ = consumes<std::vector<reco::GenJet>>(edm::InputTag("slimmedGenJets"));
		genJetsAK8Token_ = consumes<std::vector<reco::GenJet>>(edm::InputTag("slimmedGenJetsAK8"));
	}	

// jet token
	pfJetsToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJets"));
	puppiJetsToken_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsPuppi"));
	recoJetsAK8Token_ = consumes<std::vector<pat::Jet>>(edm::InputTag("slimmedJetsAK8")); 

//  MET token
	pfmetToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETs"));
	puppimetToken_ = consumes<std::vector<pat::MET>>(edm::InputTag("slimmedMETsPuppi"));
	
//now do what ever initialization is needed
	usesResource("TFileService");
	tree = (TTree*) fs->make<TTree>("Events", "events");

// event ids
    t.InitIntBranch("runNumber");
    t.InitIntBranch("eventNumber");
    t.InitIntBranch("lumiBlock");

// PU info
	t.InitIntBranch("nPV");

// PF Cands	
	t.InitIntBranch("npfCands");
	t.InitFloatArray("pfCand_Pt", "npfCands");
	t.InitFloatArray("pfCand_Eta", "npfCands");
	t.InitFloatArray("pfCand_Phi", "npfCands");
	t.InitFloatArray("pfCand_E", "npfCands");
	t.InitFloatArray("pfCand_M", "npfCands");

	
// GenJets
	if(isMC){
		t.InitIntBranch("nTruePU");

		t.InitIntBranch("nGenJets");
		t.InitFloatArray("GenJet_Pt",  "nGenJets");
		t.InitFloatArray("GenJet_Eta", "nGenJets");
		t.InitFloatArray("GenJet_Phi", "nGenJets");
		t.InitFloatArray("GenJet_E", "nGenJets");    
		t.InitIntArray("GenJet_ID",  "nGenJets");

		t.InitIntBranch("nGenJetsAK8");
		t.InitFloatArray("GenJetAK8_Pt",  "nGenJetsAK8");
		t.InitFloatArray("GenJetAK8_Eta", "nGenJetsAK8");
		t.InitFloatArray("GenJetAK8_Phi", "nGenJetsAK8");
		t.InitIntArray(  "GenJetAK8_E",  "nGenJetsAK8");
	}
	
    // PFjet values
	t.InitIntBranch("nPFJets");
	t.InitFloatArray("PFJet_Pt",  "nPFJets");
	t.InitFloatArray("PFJet_Eta", "nPFJets");
	t.InitFloatArray("PFJet_Phi", "nPFJets");
	t.InitFloatArray("PFJet_E", "nPFJets");    
	t.InitIntArray("PFJet_partonFlav", "nPFJets");
	t.InitIntArray("PFJet_hadronFlav", "nPFJets");
	t.InitIntArray("PFJet_hasGenJet","nPFJets");
	t.InitIntArray("PFJet_GenJet_idx","nPFJets");
	t.InitIntArray("PFJet_GenJet_dR","nPFJets");	
	
	// Puppi jet values

	t.InitIntBranch("nPuppiJets");
	t.InitFloatArray("PuppiJet_Pt",  "nPuppiJets");
	t.InitFloatArray("PuppiJet_Eta", "nPuppiJets");
	t.InitFloatArray("PuppiJet_Phi", "nPuppiJets");
	t.InitFloatArray("PuppiJet_E", "nPuppiJets");    
//  t.InitIntArray("PuppiJet_partonFlav", "nPuppiJets");
//  t.InitIntArray("PuppiJet_hadronFlav", "nPuppiJets");


	// AK8 jet values
		
    t.InitIntBranch("nrecoJetsAK8");
    t.InitFloatArray("recoJetAK8_Pt",  "nrecoJetsAK8");
    t.InitFloatArray("recoJetAK8_Eta", "nrecoJetsAK8");
    t.InitFloatArray("recoJetAK8_Phi", "nrecoJetsAK8");
    t.InitFloatArray("recoJetAK8_E", "nrecoJetsAK8");    

    t.InitIntArray("recoJetAK8_hasGenJet","nrecoJetsAK8");
    t.InitIntArray("recoJetAK8_GenJetAK8_idx","nrecoJetsAK8");
    t.InitIntArray("recoJetAK8_GenJetAK8_dR","nrecoJetsAK8");	
    

/*
	// DR jet values
    t.InitIntBranch("nrecoJetsDR");
    t.InitFloatArray("recoJetDR_Pt",  "nrecoJetsDR");
    t.InitFloatArray("recoJetDR_Eta", "nrecoJetsDR");
    t.InitFloatArray("recoJetDR_Phi", "nrecoJetsDR");
    t.InitFloatArray("recoJetDR_E", "nrecoJetsDR");
    t.InitFloatArray("recoJetDR_Radius", "nrecoJetsDR");
*/
	// MET values
	
	t.InitFloatBranch("pfMET");
	t.InitFloatBranch("puppiMET");

}


Analysis::~Analysis()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	using namespace edm;
	t.SetDeafultValues();
	
	t.FillIntBranch("eventNumber",iEvent.id().event());
	t.FillIntBranch("runNumber",iEvent.id().run());
	t.FillIntBranch("lumiBlock",iEvent.luminosityBlock());   

// get gen jets
	if(isMC){
//		iEvent.getByToken(genPartToken_, genParticles);
		iEvent.getByToken(genJetsToken_,    genJets);
		iEvent.getByToken(genJetsAK8Token_,    genJetsAK8);
	}
	
// get reco jets
	iEvent.getByToken(pfJetsToken_, pfJets);
	iEvent.getByToken(puppiJetsToken_, puppiJets);
	iEvent.getByToken(recoJetsAK8Token_, recoJetsAK8);
	
// get pf cand
	iEvent.getByToken(pfCandsToken_, pfCands);
	
// get MET info
	iEvent.getByToken(pfmetToken_, pfMET);
	iEvent.getByToken(puppimetToken_, puppiMET);
	
// get PU info
	iEvent.getByToken(PUsrc_, PUInfo);
	iEvent.getByToken(PVsrc_, PVInfo);
	
	
//	iEvent.getByToken(tagInfoToken_,       tagInfo);

    int nPFJets = 0;
    for(auto j = pfJets->begin(); j != pfJets->end(); ++j){
        // jet selection
        if(j->pt() < 10. || abs(j->eta()) > 4.0) continue;
        t.FillFloatArray("PFJet_Pt",   nPFJets, j->pt());
        t.FillFloatArray("PFJet_Eta",  nPFJets, j->eta());
        t.FillFloatArray("PFJet_Phi",  nPFJets, j->phi());
        t.FillFloatArray("PFJet_E",  nPFJets, j->energy());
        t.FillFloatArray("PFJet_M",  nPFJets, j->mass());
           
        t.FillFloatArray("PFJet_partonFlav",  nPFJets, j->partonFlavour());   
        t.FillFloatArray("PFJet_hadronFlav",  nPFJets, j->hadronFlavour()); 
             
        int foundGenJet = 0;
        for (auto gj = genJets->begin(); gj != genJets->end(); ++gj)
        {
            float dR = ROOT::Math::VectorUtil::DeltaR(j->p4(), gj->p4());
            if(dR < dRthreshold)
            {
                foundGenJet = 1;
                // fill gen jet info
                t.FillFloatArray("PFJet_GenJet_DeltaR",  nPFJets, dR);
                t.FillIntArray(  "PFJet_GenJet_idx", nPFJets, gj - genJets->begin());
                break;
            }
        }
        t.FillIntArray("PFJet_hasGenJet", nPFJets, foundGenJet);
        nPFJets++;
	}
	t.FillIntBranch("nPFJets", nPFJets);
	
    int nPuppiJets = 0;
    for(auto j = puppiJets->begin(); j != puppiJets->end(); ++j){
        // jet selection
        if(j->pt() < pt_jet || abs(j->eta()) > 4.0) continue;
        t.FillFloatArray("PuppiJet_Pt",   nPuppiJets, j->pt());
        t.FillFloatArray("PuppiJet_Eta",  nPuppiJets, j->eta());
        t.FillFloatArray("PuppiJet_Phi",  nPuppiJets, j->phi());
        t.FillFloatArray("PuppiJet_E",  nPuppiJets, j->energy());
        t.FillFloatArray("PuppiJet_M",  nPuppiJets, j->mass());        
        nPuppiJets++;
    }
    t.FillIntBranch("nPuppiJets", nPuppiJets);

    int nrecoJetsAK8 = 0;
    for(auto j = recoJetsAK8->begin(); j != recoJetsAK8->end(); ++j){
        // jet selection
        if(j->pt() < 100. || abs(j->eta()) > 4.0) continue;
        t.FillFloatArray("recoJetAK8_Pt",   nrecoJetsAK8, j->pt());
        t.FillFloatArray("recoJetAK8_Eta",  nrecoJetsAK8, j->eta());
        t.FillFloatArray("recoJetAK8_Phi",  nrecoJetsAK8, j->phi());
        t.FillFloatArray("recoJetAK8_E",  nrecoJetsAK8, j->energy());
        t.FillFloatArray("recoJetAK8_M",  nrecoJetsAK8, j->mass());        
        int foundGenJet = 0;
        for (auto gj = genJetsAK8->begin(); gj != genJetsAK8->end(); ++gj)
        {
            float dR = ROOT::Math::VectorUtil::DeltaR(j->p4(), gj->p4());
            if(dR < dRthreshold)
            {
                foundGenJet = 1;
                // fill gen jet info
                t.FillFloatArray("recoJetAK8_GenJetAK8_dR",  nrecoJetsAK8, dR);
                t.FillIntArray(  "recoJetAK8_GenJetAK8_idx", nrecoJetsAK8, gj - genJets->begin());
                break;
            }
        }
        t.FillIntArray("recoJetAK8_hasGenJet", nrecoJetsAK8, foundGenJet);
        
        nrecoJetsAK8++;
    }
    t.FillIntBranch("nrecoJetsAK8", nrecoJetsAK8); 
    
    int npfcands=0;
    for(auto cand=pfCands->begin(); cand != pfCands->end(); ++cand){
	if(cand->pt() < pt_cand || abs(cand->eta())> 4.0) continue;
        t.FillFloatArray("pfCand_Pt",   npfcands, cand->pt());
        t.FillFloatArray("pfCand_Eta",  npfcands, cand->eta());
        t.FillFloatArray("pfCand_Phi",  npfcands, cand->phi());
        t.FillFloatArray("pfCand_E",  npfcands, cand->energy());
        t.FillFloatArray("pfCand_M",  npfcands, cand->mass());
	npfcands++;
   }
   t.FillIntBranch("npfCands",npfcands);
		    
/*    int nrecoJetsDR = 0;
    for(auto j = recoJetsDR->begin(); j != recoJetsDR->end(); ++j)
    {
        // jet selection
        if(j->pt() < 10. || abs(j->eta()) > 4.0) continue;
        double radius = j->radius();
        t.FillFloatArray("recoJetDR_Pt",   nrecoJetsDR, j->pt());
        t.FillFloatArray("recoJetDR_Eta",  nrecoJetsDR, j->eta());
        t.FillFloatArray("recoJetDR_Phi",  nrecoJetsDR, j->phi());
        t.FillFloatArray("recoJetDR_E",  nrecoJetsDR, j->energy());
        t.FillFloatArray("recoJetDR_M",  nrecoJetsDR, j->mass());        
        t.FillFloatArray("recoJetDR_Radius", nrecoJetsDR, radius);
        nrecoJetsDR++;
    }
    t.FillIntBranch("nrecoJetsDR", nrecoJetsDR);
*/      
    // loop over gen jets
	if(isMC){
		int nGenJets = 0;
		for(auto j = genJets->begin(); j != genJets->end(); ++j){
			int i = j - genJets->begin();
			t.FillFloatArray("GenJet_Pt",  i, j->pt());
			t.FillFloatArray("GenJet_Eta", i, j->eta());
			t.FillFloatArray("GenJet_Phi", i, j->phi());
			t.FillFloatArray("GenJet_E", i, j->energy());        
			t.FillIntArray(  "GenJet_ID",  i, j->pdgId());
			nGenJets++;
		}
		t.FillIntBranch("nGenJets", nGenJets);
		
		int nGenJetsAK8 = 0;
		for(auto j = genJetsAK8->begin(); j != genJetsAK8->end(); ++j){
			int i = j - genJetsAK8->begin();
			t.FillFloatArray("GenJetAK8_Pt",  i, j->pt());
			t.FillFloatArray("GenJetAK8_Eta", i, j->eta());
			t.FillFloatArray("GenJetAK8_Phi", i, j->phi());
			t.FillFloatArray("GenJetAK8_E", i, j->energy());        
			nGenJetsAK8++;
		}
		t.FillIntBranch("nGenJetsAK8", nGenJetsAK8);
	}
	
	pf_met = pfMET->front().pt();
	puppi_met = puppiMET->front().pt();

	t.FillFloatBranch("pfMET", pf_met);
	t.FillFloatBranch("puppiMET",puppi_met);
	
	if(isMC){
		for(auto PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) {
			nTruePU = PVI->getTrueNumInteractions(); 
		}
		t.FillIntBranch("nTruePU", nTruePU);
	}

	nPV = PVInfo->size();	
	t.FillIntBranch("nPV", nPV);
	tree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
Analysis::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analysis::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analysis);
