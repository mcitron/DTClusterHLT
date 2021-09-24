// -*- C++ -*-
//
// Package:    DTCluster/DTClusterFilter
// Class:      DTClusterFilter
//
/**\class DTClusterFilter DTClusterFilter.cc DTCluster/DTClusterFilter/plugins/DTClusterFilter.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Matthew Daniel Citron
//         Created:  Thu, 15 Jul 2021 03:19:13 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/JetReco/interface/BasicJet.h"

#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
//
// class declaration
//
namespace edm {
  class ConfigurationDescriptions;
}

class DTClusterFilter : public HLTFilter {
    public:
	explicit DTClusterFilter(const edm::ParameterSet& iConfig);
	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
	bool hltFilter(edm::Event&, const edm::EventSetup&, trigger::TriggerFilterObjectWithRefs& filterproduct) const;
    private:


	//virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
	//virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	//virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
	edm::InputTag dtClusterLabel_;
        unsigned int minClusterSize_;
        unsigned int minClusterNum_;
	edm::EDGetTokenT<reco::BasicJetCollection> dtClusterInputToken;

	// ----------member data ---------------------------
#ifdef THIS_IS_AN_EVENT_EXAMPLE
	edm::EDGetTokenT<ExampleData> exampleToken_;
#endif
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
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
DTClusterFilter::DTClusterFilter(const edm::ParameterSet& iConfig) : HLTFilter(iConfig){
    dtClusterLabel_= iConfig.getParameter<edm::InputTag>("dtClusters");
    minClusterSize_= iConfig.getParameter<unsigned int>("minClusterSize");
    minClusterNum_= iConfig.getParameter<unsigned int>("minClusterNum");
    dtClusterInputToken = consumes<std::vector<reco::BasicJet>>(dtClusterLabel_);
    //now do what ever initialization is needed
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool DTClusterFilter::hltFilter(edm::Event& iEvent, const edm::EventSetup& iSetup,trigger::TriggerFilterObjectWithRefs& filterproduct) const {
    bool accept = false;
    edm::Handle<reco::BasicJetCollection> dtClusters;
    iEvent.getByToken(dtClusterInputToken, dtClusters);
    unsigned int ndtClusters = 0;
    for (auto const& c : *dtClusters) {
	if (c.energy() >= minClusterSize_) ndtClusters++;
    }
    accept = ndtClusters >= minClusterNum_;
    return accept;
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DTClusterFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    makeHLTFilterDescription(desc);
    desc.add<edm::InputTag>("dtClusters", edm::InputTag(""));
    desc.add<unsigned int>("minClusterSize", 1);
    desc.add<unsigned int>("minClusterNum", 1);
    descriptions.add("dtClusterFilter", desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(DTClusterFilter);
