#include "MagneticField/Engine/interface/MagneticField.h"
#include "SimG4Core/MagneticField/interface/Field.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "G4Mag_UsualEqRhs.hh"

#include <CLHEP/Units/SystemOfUnits.h>

using namespace sim;

Field::Field(const MagneticField *f, double d) : G4MagneticField(), theCMSMagneticField(f), theDelta(d), dxi(0.) {
  for (int i = 0; i < 3; ++i) {
    oldx[i] = 1.0e12;
    oldb[i] = 0.0;
    offset[i] = 0.0;
  }
}

Field::~Field() {}

void Field::SetOffset(double x, double y, double z) {
  offset[0] = x * CLHEP::tesla;
  offset[1] = y * CLHEP::tesla;
  offset[2] = z * CLHEP::tesla;
}
