#ifndef __ACTION__
#define __ACTION__

#include "G4VUserActionInitialization.hh"

#include "ScatteringRadiation/Beam.h"
#include "ScatteringRadiation/Event.h"
#include "ScatteringRadiation/Stepping.h"
#include "ScatteringRadiation/Run.h"

class Action : public G4VUserActionInitialization {

public :

   Action() {};
   virtual ~Action() {};

   virtual void Build() const;


};

#endif