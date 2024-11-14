#ifndef TrackPropagators_ESProducers_GeantPropagatorESProducercvh_h
#define TrackPropagators_ESProducers_GeantPropagatorESProducercvh_h

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include <memory>

/*
 * GeantPropagatorESProducercvh
 *
 * Produces an Geant4ePropagator for track propagation
 *
 */

class GeantPropagatorESProducercvh : public edm::ESProducer {
public:
  GeantPropagatorESProducercvh(const edm::ParameterSet &p);
  ~GeantPropagatorESProducercvh() override = default;

  std::unique_ptr<Propagator> produce(const TrackingComponentsRecord &);

private:
  edm::ParameterSet pset_;
  const std::string fieldlabel_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
  double plimit_;
};

#endif
