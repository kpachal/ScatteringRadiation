#ifndef __ACTION__
#define __ACTION__

#include "G4VUserActionInitialization.hh"
#include "ScatteringRadiation/Beam.h"

class Action : public G4VUserActionInitialization {

public :

   Action() {};
   virtual ~Action() {};

   virtual void Build() const;


};

#endif