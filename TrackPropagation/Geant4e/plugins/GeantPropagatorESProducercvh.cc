#include "GeantPropagatorESProducercvh.h"

#include "TrackPropagation/Geant4e/interface/Geant4ePropagator.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include <memory>
#include <string>

using namespace edm;

GeantPropagatorESProducercvh::GeantPropagatorESProducercvh(const edm::ParameterSet &p)
    : fieldlabel_(p.getParameter<std::string>("MagneticFieldLabel")),
      magFieldToken_(setWhatProduced(this, p.getParameter<std::string>("ComponentName"))
                         .consumesFrom<MagneticField, IdealMagneticFieldRecord>(edm::ESInputTag("", fieldlabel_))) {
  pset_ = p;
  plimit_ = pset_.getParameter<double>("PropagationPtotLimit");
}

std::unique_ptr<Propagator> GeantPropagatorESProducercvh::produce(const TrackingComponentsRecord &iRecord) {
  std::string pdir = pset_.getParameter<std::string>("PropagationDirection");
  std::string particleName = pset_.getParameter<std::string>("ParticleName");

  PropagationDirection dir = alongMomentum;

  if (pdir == "oppositeToMomentum") {
    dir = oppositeToMomentum;
  } else if (pdir == "alongMomentum") {
    dir = alongMomentum;
  } else if (pdir == "anyDirection") {
    dir = anyDirection;
  }

  return std::make_unique<Geant4ePropagator>(&(iRecord.get(magFieldToken_)), particleName, dir, plimit_);
}
