#ifndef __EVENT__
#define __EVENT__

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "g4root.hh"

#include "ScatteringRadiation/Run.h"

class EventAction : public G4UserEventAction {

  public :
    EventAction();
    ~EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

  private :

};

#endif

