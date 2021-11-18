#ifndef __TRACK__
#define __TRACK__

#include "G4UserTrackingAction.hh"
#include "ScatteringRadiation/Event.h"

class TrackingAction : public G4UserTrackingAction {

  public :
    TrackingAction(EventAction* eventAction);
    ~TrackingAction();

    virtual void PreUserTrackingAction(const G4Track* track);
    virtual void PostUserTrackingAction(const G4Track* track);

  private :

    EventAction * m_EventAction;

};

#endif

