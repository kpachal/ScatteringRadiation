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

    void setTrackHighAngleOnly(bool saveAngleOnly, float minimumAngle);

  private :

    EventAction * m_EventAction;

    bool m_save_highangle_only; // Keep only final state particles at high angles? This is useful for estimating physics backgrounds.
    float m_minimum_angle; // Minimum angle from z axis below which we will ignore particles

};

#endif

