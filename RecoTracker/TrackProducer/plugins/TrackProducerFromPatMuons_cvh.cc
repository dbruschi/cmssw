#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

class TrackProducerFromPatMuons_cvh : public edm::global::EDProducer<>
{
public:
  explicit TrackProducerFromPatMuons_cvh(const edm::ParameterSet &);
  ~TrackProducerFromPatMuons_cvh() override = default;

//   static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  
  void produce(edm::StreamID, edm::Event &, const edm::EventSetup &) const override;
  
  edm::EDGetTokenT<std::vector<pat::Muon>> inputMuons_;
  edm::EDPutTokenT<reco::TrackCollection> outputTrack_;
  edm::EDPutTokenT<edm::Association<std::vector<pat::Muon>>> outputAssoc_;
  bool innerTrackOnly_;
  double ptMin_;

};


TrackProducerFromPatMuons_cvh::TrackProducerFromPatMuons_cvh(const edm::ParameterSet &iConfig)

{
  inputMuons_ = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("src"));
  outputTrack_ = produces<reco::TrackCollection>();
  outputAssoc_ = produces<edm::Association<std::vector<pat::Muon>>>();
  innerTrackOnly_ = iConfig.getParameter<bool>("innerTrackOnly");
  ptMin_ = iConfig.getParameter<double>("ptMin");
}

// ------------ method called for each event  ------------
void TrackProducerFromPatMuons_cvh::produce(edm::StreamID, edm::Event &iEvent, const edm::EventSetup &iSetup) const {

  using namespace edm;
  
  Handle<std::vector<pat::Muon>> muons;
  iEvent.getByToken(inputMuons_, muons);

  reco::TrackCollection tracksOut;

  edm::Association<std::vector<pat::Muon>> association;
  edm::Association<std::vector<pat::Muon>>::Filler assocfiller(association);

  edm::RefProd<std::vector<pat::Muon>> muonRefProd(muons);
  association.setRef(muonRefProd);
  
  std::vector<int> associdxs;

  for (unsigned int iMuon = 0; iMuon < muons->size(); ++iMuon) {
    auto const &muon = (*muons)[iMuon];

    if (muon.pt() < ptMin_) continue;

    const reco::TrackRef trackRef = innerTrackOnly_ ? muon.innerTrack() : muon.muonBestTrack();
    if (trackRef.isNonnull() && trackRef->extra().isAvailable()) {
      tracksOut.emplace_back(*trackRef);
      associdxs.push_back(iMuon);
    }
  }
  
  auto trackouth = iEvent.emplace(outputTrack_, std::move(tracksOut));
  assocfiller.insert(trackouth, associdxs.begin(), associdxs.end());
  assocfiller.fill();
  iEvent.emplace(outputAssoc_, std::move(association));


}

DEFINE_FWK_MODULE(TrackProducerFromPatMuons_cvh);
