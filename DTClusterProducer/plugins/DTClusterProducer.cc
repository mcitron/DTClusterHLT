// -*- C++ -*-
//
// Package:    DTClusterHLT/DTClusterProducer
// Class:      DTClusterProducer
//
/**\class DTClusterProducer DTClusterProducer.cc DTClusterHLT/DTClusterProducer/plugins/DTClusterProducer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Matthew Citron <mcitron@ucsb.edu> 10/19/2017
//         Created:  Thu, 12 Aug 2021 21:27:41 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "RecoJets/JetProducers/interface/JetSpecific.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "RecoJets/JetProducers/plugins/FastjetJetProducer.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "DataFormats/JetReco/interface/Jet.h"
//
// class declaration
//

class DTClusterProducer : public edm::stream::EDProducer<> {
    public:
	explicit DTClusterProducer(const edm::ParameterSet&);
	~DTClusterProducer();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
	void beginStream(edm::StreamID) override;
	void produce(edm::Event&, const edm::EventSetup&) override;
	void endStream() override;

	edm::InputTag dtRechitLabel_;
	edm::EDGetTokenT<DTRecHitCollection> dtRechitInputToken_;
	double _rParam;
	bool _skipMB1;

	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

	// ----------member data ---------------------------
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
DTClusterProducer::DTClusterProducer(const edm::ParameterSet& iConfig)
{
    dtRechitLabel_ = iConfig.getParameter<edm::InputTag>("dtRechits");
    dtRechitInputToken_ = consumes<DTRecHitCollection>(edm::InputTag(dtRechitLabel_));
    _rParam= iConfig.getParameter<double>("rParam");
    _skipMB1 = iConfig.getParameter<bool>("skipMB1");
    produces<reco::BasicJetCollection>();

    //register your products
    /* Examples
       produces<ExampleData2>();

    //if do put with a label
    produces<ExampleData2>("label");

    //if you want to put into the Run
    produces<ExampleData2,InRun>();
    */
    //now do what ever other initialization is needed
}

DTClusterProducer::~DTClusterProducer() {
    // do anything here that needs to be done at destruction time
    // (e.g. close files, deallocate resources etc.)
    //
    // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void DTClusterProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    using namespace edm;
    edm::Handle<DTRecHitCollection> dtRechits;
    iEvent.getByToken(dtRechitInputToken_,dtRechits);
    std::vector<fastjet::PseudoJet> fjInputs;
    edm::ESHandle<DTGeometry> dtG;
    iSetup.get<MuonGeometryRecord>().get(dtG); 
    for(DTRecHit1DPair dtRechit: *dtRechits){
	LocalPoint  localPosition       = dtRechit.localPosition();
	DetId geoid = dtRechit.geographicalId();
	DTChamberId dtdetid = DTChamberId(geoid);
	const DTChamber * dtchamber = dtG->chamber(dtdetid);
	if (_skipMB1 && dtdetid.station() == 1) continue;
	if (dtchamber) {
	    GlobalPoint globalPosition = dtchamber->toGlobal(localPosition);
	    double norm = sqrt(globalPosition.x()*globalPosition.x()+globalPosition.y()*globalPosition.y()+globalPosition.z()*globalPosition.z());
	    fastjet::PseudoJet j(globalPosition.x()/norm,globalPosition.y()/norm,globalPosition.z()/norm,1);
	    fjInputs.push_back(j);
	}
    }
    std::shared_ptr<fastjet::JetDefinition> fjJetDefinition = std::make_shared<fastjet::JetDefinition>(fastjet::cambridge_algorithm, _rParam);
    std::shared_ptr<fastjet::ClusterSequence> fjClusterSeq = std::make_shared<fastjet::ClusterSequence>(fjInputs, *fjJetDefinition);
    std::vector<fastjet::PseudoJet> fjJets = fastjet::sorted_by_E(fjClusterSeq->inclusive_jets(0));
    auto jetCollection = std::make_unique<reco::BasicJetCollection>();
    for (fastjet::PseudoJet fjJet: fjJets){
	reco::Particle::Point point(0, 0, 0);
	double norm = fjJet.e()/sqrt(fjJet.px()*fjJet.px()+fjJet.py()*fjJet.py()+fjJet.pz()*fjJet.pz());
	reco::BasicJet toput(math::XYZTLorentzVector(fjJet.px()*norm, fjJet.py()*norm, fjJet.pz()*norm, fjJet.e()), point);
	// std::cout << fjJet.e() << " " << toput.energy() << std::endl;
	jetCollection->push_back(toput);
    }
    iEvent.put(std::move(jetCollection));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void DTClusterProducer::beginStream(edm::StreamID) {
    // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void DTClusterProducer::endStream() {
    // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
   void
   DTClusterProducer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void
   DTClusterProducer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void
   DTClusterProducer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void
   DTClusterProducer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DTClusterProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("dtRechits", edm::InputTag("hltDt1DRecHits"));
    desc.add<double>("rParam", 0.2);
    desc.add<bool>("skipMB1", true);
    descriptions.add("dtClusterProducer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DTClusterProducer);
