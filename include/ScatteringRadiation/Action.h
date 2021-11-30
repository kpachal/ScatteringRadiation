#ifndef __ACTION__
#define __ACTION__

#include "G4VUserActionInitialization.hh"

#include "ScatteringRadiation/Beam.h"
#include "ScatteringRadiation/Event.h"
#include "ScatteringRadiation/Stepping.h"
#include "ScatteringRadiation/Run.h"
#include "ScatteringRadiation/Track.h"

class Action : public G4VUserActionInitialization {

public :

   Action(std::string outputFilename) {m_outputFilename = outputFilename; };
   virtual ~Action() {};

   virtual void Build() const;

private : 

   std::string m_outputFilename;

};

#endif