#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalCommonData/interface/HcalHitRelabeller.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonClient.h"
#include "HeterogeneousCore/SonicCore/interface/SonicEDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"

//#include "RecoLocalCalo/HcalRecProducers/src/HcalPhase1Reconstructor_FACILE2.h"
class HcalPhase1Reconstructor_FACILE : public SonicEDProducer<TritonClient>
{
public:
    explicit HcalPhase1Reconstructor_FACILE(edm::ParameterSet const& cfg) : 
        SonicEDProducer<TritonClient>(cfg),
        //nInputs_(cfg.getParameter<int>("nInputs")),
        //batchSize_(cfg.getParameter<int>("batchSize")),
        fChannelInfoName(cfg.getParameter<edm::InputTag>("ChannelInfoName")),
        fTokChannelInfo(this->template consumes<HBHEChannelInfoCollection>(fChannelInfoName))
     
    {
        this->template produces<HBHERecHitCollection>();
        this->setDebugName("HcalPhase1Reconstructor_FACILE");
    }

    void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, Input& iInput) override {

	const HBHEChannelInfoCollection *channelInfo   = 0;
	edm::Handle<HBHEChannelInfoCollection> hChannelInfo;
	iEvent.getByToken(fTokChannelInfo, hChannelInfo); 
	channelInfo = hChannelInfo.product();

	//iInput = Input(1000, 0.f);
	auto& input1 = iInput.begin()->second;
	auto data1 = std::make_shared<TritonInput<float>>();
	data1->reserve(input1.batchSize());
	
        unsigned int iBatch = 0;
	for(HBHEChannelInfoCollection::const_iterator itC = channelInfo->begin(); itC != channelInfo->end(); itC++){

	    //keeptrack of batching for inputs
	    iBatch = std::distance(channelInfo->begin(),itC);

	    const HBHEChannelInfo& pChannel(*itC);
  	    const HcalDetId        pDetId = pChannel.id();
 	    hcalIds.push_back(pDetId);    

	    //FACILE uses iphi as a continuous variable
	    //iInput[iBatch*nInputs_ + 0] = (float)pDetId.iphi();
	    //iInput[iBatch*nInputs_ + 1] = (float)pChannel.tsGain(0);
	    data1->emplace_back(iBatch, (float)pDetId.iphi());
	    
	    for (unsigned int itTS=0; itTS < pChannel.nSamples(); ++itTS) {
		//iInput[iBatch*nInputs_ + itTS + 2] = (float)pChannel.tsRawCharge(itTS);
		data1->emplace_back(iBatch, (float)pChannel.tsRawCharge(itTS));
	    }
  
	    //FACILE considers 7 Hcal depths as binary variables
            for (int itDepth=1; itDepth < 8; itDepth++){
		if (pDetId.depth() == itDepth)  data1->emplace_back(iBatch, 1.);//iInput[iBatch*nInputs_ + itDepth + 9] = 1.;
		else				data1->emplace_back(iBatch, 0.);//iInput[iBatch*nInputs_ + itDepth + 9] = 0.;
	    }

	    //ieta is also encoded as a binary variable
	    for (int itIeta = 0; itIeta < 30; itIeta++){
		if (std::abs(pDetId.ieta()) == itIeta)  data1->emplace_back(iBatch, 1.);//iInput[iBatch*nInputs_ + itIeta + 17] = 1.;
		else					data1->emplace_back(iBatch, 0.);//iInput[iBatch*nInputs_ + itIeta + 17] = 0.;
	    }
	}

	input1.toServer(data1);
    }

    void produce(edm::Event& iEvent, edm::EventSetup const& iSetup, Output const& iOutput) override {
	std::unique_ptr<HBHERecHitCollection> out;
	out = std::make_unique<HBHERecHitCollection>();
	out->reserve(hcalIds.size());

	unsigned int iBatch = 0;
	const auto& output1 = iOutput.begin()->second;
	const auto& tmp = output1.fromServer<float>();

	for(HBHERecHitCollection::const_iterator itRH = out->begin(); itRH != out->end(); itRH++){

	    //float rhE = iOutput[iBatch];
	    float rhE = tmp[0][iBatch];
	    //FACILE uses rectified linear activation function =>should be positive definite
	    if(rhE < 0.) rhE = 0.;
	    //throw cms exception?
	    if(std::isnan(rhE)) rhE = 0; 
	    if(std::isinf(rhE)) rhE = 0; 

	    //FACILE does no time reco 
	    HBHERecHit rh = HBHERecHit(hcalIds[iBatch],rhE,0.f,0.f);
	    out->push_back(rh); //hcalIds[iBatch],rhE,0.f,0.f);

	    iBatch++;
	}
	iEvent.put(std::move(out));	
    }

    /*static void fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
        edm::ParameterSetDescription desc;
        TritonClient::fillPSetDescription(desc);
        //add producer-specific parameters
 	desc.add<edm::InputTag>("ChannelInfoName","hbheprereco");
        descriptions.add("HcalPhase1Reconstructor_FACILE",desc);
    }*/


private:
    //int nInputs_, batchSize_;
    edm::InputTag fChannelInfoName;   
    edm::EDGetTokenT<HBHEChannelInfoCollection> fTokChannelInfo;
    std::vector<HcalDetId> hcalIds;
};

DEFINE_FWK_MODULE(HcalPhase1Reconstructor_FACILE);
