#include "TFile.h"
#include "ScatteringRadiation/Track.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

TrackingAction::TrackingAction(EventAction * eventAction) {
    m_EventAction = eventAction;
    m_save_highangle_only = false;
    m_minimum_angle = 0;
}

TrackingAction::~TrackingAction() {

}

void TrackingAction::setTrackHighAngleOnly(bool saveAngleOnly, float minimumAngle) {
    m_save_highangle_only = saveAngleOnly;
    m_minimum_angle = minimumAngle*CLHEP::degree;
}

void TrackingAction::PreUserTrackingAction(const G4Track* track) {

}

void TrackingAction::PostUserTrackingAction(const G4Track* track) {

  // If we only want high angle tracks, kill track if it isn't.
  if (m_save_highangle_only) {
    CLHEP::Hep3Vector zdir(0, 0, 1);
    float angle = zdir.angle(track->GetMomentum());
    if (angle < m_minimum_angle) return;
  } 


  // Want to fill info at end of track
  auto position = track->GetPosition();
  auto momentum = track->GetMomentum();
  auto pdgID = track->GetParticleDefinition()->GetPDGEncoding();
  auto mass = track->GetParticleDefinition()->GetPDGMass();
  m_EventAction->FillTrackInfo(pdgID, position[0], position[1], position[2],
                                 momentum[0], momentum[1], momentum[2], mass);
}
