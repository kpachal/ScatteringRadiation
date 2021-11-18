#include "ScatteringRadiation/Track.h"
#include "TFile.h"

TrackingAction::TrackingAction(EventAction * eventAction) {
    m_EventAction = eventAction;
}

TrackingAction::~TrackingAction() {

}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {

}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {

  // Want to fill info at end of track
  auto position = track->GetPosition();
  auto momentum = track->GetMomentumDirection();
  auto pdgID = track->GetParticleDefinition()->GetPDGEncoding();
  m_EventAction->FillTrackInfo(pdgID, position[0], position[1], position[2],
                                 momentum[0], momentum[1], momentum[2]);
}
