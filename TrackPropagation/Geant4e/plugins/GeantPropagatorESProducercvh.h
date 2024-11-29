#ifndef TrackPropagators_ESProducers_GeantPropagatorESProducercvh_h
#define TrackPropagators_ESProducers_GeantPropagatorESProducercvh_h

#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include <memory>

/*
 * GeantPropagatorESProducer
 *
 * Produces an Geant4ePropagator for track propagation
 *
 */

class GeantPropagatorESProducercvh : public edm::ESProducer {
public:
  GeantPropagatorESProducercvh(const edm::ParameterSet &p);
  ~GeantPropagatorESProducercvh() override;

  std::unique_ptr<Propagator> produce(const TrackingComponentsRecord &);

private:
  edm::ParameterSet pset_;
  double plimit_;
  std::string fieldlabel_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magFieldToken_;
};

#endif

