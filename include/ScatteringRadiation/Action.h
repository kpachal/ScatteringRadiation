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

   Action(std::string outputFilename, bool doSaveOnly, int PDGID_only, bool saveHighAngleOnly) {
      m_outputFilename = outputFilename; 
      m_do_save_only = doSaveOnly;
      m_PDGID_only = PDGID_only;
      m_save_highangle_only = saveHighAngleOnly;
   };
   
   virtual ~Action() {};

   virtual void Build() const;

private : 

   std::string m_outputFilename;
   bool m_do_save_only;
   int m_PDGID_only;
   bool m_save_highangle_only;

};

#endif